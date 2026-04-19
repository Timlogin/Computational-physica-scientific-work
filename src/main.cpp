#include "gmsh_box.hpp"
#include "simulation.hpp"
#include "vtk_writer.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace
{
// Формирует имя файла кадра
std::string frame_name(const std::string& output_dir,
                       const std::string& prefix,
                       const std::string& extension,
                       int frame_id)
{
  std::ostringstream out;
  out << output_dir << "/" << prefix << "/" << prefix << "_" << std::setw(4) << std::setfill('0') << frame_id
      << extension;
  return out.str();
}

// Формирует относительное имя файла для записи в PVD
std::string relative_frame_name(const std::string& prefix, const std::string& extension, int frame_id)
{
  std::ostringstream out;
  out << prefix << "/" << prefix << "_" << std::setw(4) << std::setfill('0') << frame_id << extension;
  return out.str();
}
}

int main()
{
  try
  {
    // Результаты каждого запуска кладём в одну и ту же папку,
    // полностью пересоздавая её
    const std::string output_dir = "results/latest_run";
    std::filesystem::remove_all(output_dir);
    std::filesystem::create_directories(output_dir + "/water");
    std::filesystem::create_directories(output_dir + "/box");
    std::filesystem::create_directories(output_dir + "/ground");
    std::filesystem::create_directories(output_dir + "/mesh");

    // Создаём конфигурацию и отдельно строим сетку бокса.
    SimulationConfig config;
    const TetMesh box_mesh
        = build_box_shell_mesh(config.target, config.box_mesh_size, output_dir + "/mesh/box_shell.msh");

    write_box_vtu(
        output_dir + "/mesh/box_shell.vtu",
        box_mesh,
        {0.0, 0.0, 0.0},
        config.target,
        {0.0, 0.0, 0.0},
        0.0,
        {});

    // После подготовки геометрии запускаем физическую симуляцию
    SphSimulation simulation(config);
    simulation.initialize_fluid_block();

    std::vector<FrameInfo> water_frames;
    std::ofstream box_loads(output_dir + "/box_loads.csv");
    box_loads << "frame_id,time,tx,ty,tz,roof_load,roof_center_x,roof_center_y,side_x_load,side_y_load\n";

    // Между сохранениями кадра накапливаем максимум нагрузки и средний центр удара
    double interval_roof_peak_load = 0.0;
    double interval_roof_weighted_x = 0.0;
    double interval_roof_weighted_y = 0.0;
    double interval_roof_total_weight = 0.0;
    double interval_side_x_peak_load = 0.0;
    double interval_side_y_peak_load = 0.0;

    // Земля не меняется по времени, поэтому записываем её один раз.
    write_ground_vtp(output_dir + "/ground/ground.vtp", simulation.ground());

    // Нулевой кадр
    int frame_id = 0;
    const std::string first_water_file = frame_name(output_dir, "water", ".vtp", frame_id);
    const Vec3 box_translation = simulation.box().center - config.target.center;
    write_particles_vtp(first_water_file, simulation.particles(), config.visual_particle_radius);
    water_frames.push_back({0.0, relative_frame_name("water", ".vtp", frame_id)});
    box_loads << frame_id << ",0.0,"
              << box_translation.x << "," << box_translation.y << "," << box_translation.z << ","
              << simulation.roof_impact_load() << ","
              << simulation.roof_impact_center_x() << ","
              << simulation.roof_impact_center_y() << ","
              << simulation.side_x_impact_load() << ","
              << simulation.side_y_impact_load() << "\n";

    // Цикл по времени
    for (int step = 0; step < config.total_steps; ++step)
    {
      simulation.step();

      // Удар воды короткий, поэтому сохраняем максимум по интервалу,
      // а не только значение в одном конкретном шаге
      const double roof_load = simulation.roof_impact_load();
      const double side_x_load = simulation.side_x_impact_load();
      const double side_y_load = simulation.side_y_impact_load();

      interval_roof_peak_load = std::max(interval_roof_peak_load, roof_load);
      interval_side_x_peak_load = std::max(interval_side_x_peak_load, side_x_load);
      interval_side_y_peak_load = std::max(interval_side_y_peak_load, side_y_load);

      if (roof_load > 1e-9)
      {
        interval_roof_total_weight += roof_load;
        interval_roof_weighted_x += roof_load * simulation.roof_impact_center_x();
        interval_roof_weighted_y += roof_load * simulation.roof_impact_center_y();
      }

      // Сохраняем кадры не на каждом шаге, а через заданный интервал
      if ((step + 1) % config.output_every == 0)
      {
        ++frame_id;
        const std::string water_file = frame_name(output_dir, "water", ".vtp", frame_id);
        const Vec3 box_step_translation = simulation.box().center - config.target.center;
        const double roof_center_x = interval_roof_total_weight > 1e-9
            ? interval_roof_weighted_x / interval_roof_total_weight
            : 0.0;
        const double roof_center_y = interval_roof_total_weight > 1e-9
            ? interval_roof_weighted_y / interval_roof_total_weight
            : 0.0;

        write_particles_vtp(water_file, simulation.particles(), config.visual_particle_radius);
        water_frames.push_back({simulation.current_time(), relative_frame_name("water", ".vtp", frame_id)});
        box_loads << frame_id << "," << simulation.current_time() << ","
                  << box_step_translation.x << "," << box_step_translation.y << "," << box_step_translation.z << ","
                  << interval_roof_peak_load << ","
                  << roof_center_x << ","
                  << roof_center_y << ","
                  << interval_side_x_peak_load << ","
                  << interval_side_y_peak_load << "\n";

        std::cout << "Saved frame " << frame_id << ", time = " << simulation.current_time() << " s\n";

        interval_roof_peak_load = 0.0;
        interval_roof_weighted_x = 0.0;
        interval_roof_weighted_y = 0.0;
        interval_roof_total_weight = 0.0;
        interval_side_x_peak_load = 0.0;
        interval_side_y_peak_load = 0.0;
      }
    }

    box_loads.close();

    write_pvd_collection(output_dir + "/water_series.pvd", water_frames);

    // После расчёта воды запускаем FEniCSx-постобработку для деформации бокса
    const char* fenics_python_env = std::getenv("FENICSX_PYTHON");
    const std::string fenics_python = fenics_python_env != nullptr
        ? fenics_python_env
        : "/opt/anaconda3/envs/fenicsx-clean/bin/python";
    const std::string fenics_script = std::filesystem::absolute("src/fenics_box_solver.py").string();
    const std::string fenics_command =
        "env FI_PROVIDER=tcp FI_TCP_IFACE=en0 " + fenics_python + " " + fenics_script + " " + output_dir;

    std::cout << "Running FEniCSx box deformation postprocess...\n";
    const int fenics_exit_code = std::system(fenics_command.c_str());
    if (fenics_exit_code != 0)
      throw std::runtime_error("FEniCSx postprocess failed.");
  }
  catch (const std::exception& error)
  {
    std::cerr << "Simulation error: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
