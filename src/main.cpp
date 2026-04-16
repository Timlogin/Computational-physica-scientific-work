#include "simulation.hpp"
#include "vtk_writer.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace
{
std::string frame_name(const std::string& output_dir, const std::string& prefix, int frame_id)
{
  std::ostringstream out;
  out << output_dir << "/" << prefix << "/" << prefix << "_" << std::setw(4) << std::setfill('0') << frame_id
      << ".vtp";
  return out.str();
}

std::string relative_frame_name(const std::string& prefix, int frame_id)
{
  std::ostringstream out;
  out << prefix << "/" << prefix << "_" << std::setw(4) << std::setfill('0') << frame_id << ".vtp";
  return out.str();
}
}

int main()
{
  try
  {
    const std::string output_dir = "results/latest_run";
    std::filesystem::remove_all(output_dir);
    std::filesystem::create_directories(output_dir + "/water");
    std::filesystem::create_directories(output_dir + "/box");

    SimulationConfig config;
    SphSimulation simulation(config);
    simulation.initialize_fluid_block();

    std::cout << "SPH demo started\n";
    std::cout << "Particles: " << simulation.particles().size() << "\n";
    std::cout << "Time step: " << config.dt << " s\n";
    std::cout << "Total steps: " << config.total_steps << "\n";
    std::cout << "Effective box density: " << simulation.box().density << " kg/m^3\n";
    std::cout << "Effective box mass: " << simulation.box_mass() << " kg\n";

    std::vector<FrameInfo> water_frames;
    std::vector<FrameInfo> box_frames;

    int frame_id = 0;
    const std::string first_box_file = frame_name(output_dir, "box", frame_id);
    const std::string first_water_file = frame_name(output_dir, "water", frame_id);
    write_box_vtp(first_box_file, simulation.box(), simulation.box_velocity(), simulation.box_mass());
    write_particles_vtp(first_water_file, simulation.particles(), config.visual_particle_radius);
    box_frames.push_back({0.0, relative_frame_name("box", frame_id)});
    water_frames.push_back({0.0, relative_frame_name("water", frame_id)});

    for (int step = 0; step < config.total_steps; ++step)
    {
      simulation.step();

      if ((step + 1) % config.output_every == 0)
      {
        ++frame_id;
        const std::string box_file = frame_name(output_dir, "box", frame_id);
        const std::string water_file = frame_name(output_dir, "water", frame_id);
        write_box_vtp(box_file, simulation.box(), simulation.box_velocity(), simulation.box_mass());
        write_particles_vtp(water_file, simulation.particles(), config.visual_particle_radius);
        box_frames.push_back({simulation.current_time(), relative_frame_name("box", frame_id)});
        water_frames.push_back({simulation.current_time(), relative_frame_name("water", frame_id)});

        std::cout << "Saved frame " << frame_id << ", time = " << simulation.current_time() << " s\n";
      }
    }

    write_pvd_collection(output_dir + "/box_series.pvd", box_frames);
    write_pvd_collection(output_dir + "/water_series.pvd", water_frames);

    std::cout << "Simulation finished. Open files from " << output_dir << " in ParaView.\n";
    std::cout << "Hint for ParaView: for water use Glyph or Point Gaussian and scale by the 'radius' array.\n";
  }
  catch (const std::exception& error)
  {
    std::cerr << "Simulation error: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
