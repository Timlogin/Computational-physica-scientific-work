#include "simulation.hpp"

#include <gmsh.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolygon.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// В этом файле специально оставлены все простые функции проекта:
// построение сетки Gmsh, запись VTK и основной цикл расчёта.
// Так код легче читать как учебный пример, без прыжков по множеству файлов.

struct TetMesh
{
  std::vector<Vec3> nodes;                    // Координаты узлов сетки бокса.
  std::vector<std::array<int, 4>> tetrahedra; // Индексы четырёх вершин каждого тетраэдра.
};

struct FrameInfo
{
  double time = 0.0;      // Время кадра.
  std::string file_name;  // Относительный путь к файлу кадра для ParaView.
};

namespace
{
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

std::string relative_frame_name(const std::string& prefix, const std::string& extension, int frame_id)
{
  std::ostringstream out;
  out << prefix << "/" << prefix << "_" << std::setw(4) << std::setfill('0') << frame_id << extension;
  return out.str();
}

void ensure_parent_directory(const std::string& file_path)
{
  const std::filesystem::path path(file_path);
  if (path.has_parent_path())
    std::filesystem::create_directories(path.parent_path());
}

vtkSmartPointer<vtkDoubleArray> make_scalar_array(const char* name, const std::vector<double>& values)
{
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetName(name);
  for (double value : values)
    array->InsertNextValue(value);
  return array;
}

vtkSmartPointer<vtkDoubleArray> make_vector_array(const char* name, const std::vector<Vec3>& values)
{
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetName(name);
  array->SetNumberOfComponents(3);

  for (const Vec3& value : values)
  {
    const double tuple[3] = {value.x, value.y, value.z};
    array->InsertNextTuple(tuple);
  }

  return array;
}

vtkSmartPointer<vtkDoubleArray> make_cell_scalar_array(const char* name, std::size_t count, double value)
{
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetName(name);
  for (std::size_t i = 0; i < count; ++i)
    array->InsertNextValue(value);
  return array;
}

vtkSmartPointer<vtkUnsignedCharArray> make_cell_rgb_array(const char* name, std::size_t count, int r, int g, int b)
{
  auto array = vtkSmartPointer<vtkUnsignedCharArray>::New();
  array->SetName(name);
  array->SetNumberOfComponents(3);

  const unsigned char color[3] = {
      static_cast<unsigned char>(r),
      static_cast<unsigned char>(g),
      static_cast<unsigned char>(b)};

  for (std::size_t i = 0; i < count; ++i)
    array->InsertNextTypedTuple(color);

  return array;
}

TetMesh build_box_shell_mesh(const BoxShell& box, double mesh_size, const std::string& msh_file_path)
{
  const Vec3 inner_size = {
      box.outer_size.x - 2.0 * box.wall_thickness,
      box.outer_size.y - 2.0 * box.wall_thickness,
      box.outer_size.z - 2.0 * box.wall_thickness};

  if (inner_size.x <= 0.0 || inner_size.y <= 0.0 || inner_size.z <= 0.0)
    throw std::runtime_error("Box wall thickness is too large.");

  TetMesh mesh;
  gmsh::initialize();

  try
  {
    gmsh::option::setNumber("General.Terminal", 0.0);
    gmsh::model::add("box_shell");

    // Строим внешний параллелепипед и вырезаем из него внутренний.
    const Vec3 outer_min = box.center - box.outer_size * 0.5;
    const Vec3 inner_min = outer_min + Vec3{box.wall_thickness, box.wall_thickness, box.wall_thickness};

    const int outer_tag = gmsh::model::occ::addBox(
        outer_min.x, outer_min.y, outer_min.z, box.outer_size.x, box.outer_size.y, box.outer_size.z);
    const int inner_tag = gmsh::model::occ::addBox(
        inner_min.x, inner_min.y, inner_min.z, inner_size.x, inner_size.y, inner_size.z);

    std::vector<std::pair<int, int>> out;
    std::vector<std::vector<std::pair<int, int>>> out_map;
    gmsh::model::occ::cut({{3, outer_tag}}, {{3, inner_tag}}, out, out_map);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMin", mesh_size);
    gmsh::option::setNumber("Mesh.MeshSizeMax", mesh_size);
    gmsh::option::setNumber("Mesh.Algorithm3D", 1.0);
    gmsh::option::setNumber("Mesh.Optimize", 1.0);
    gmsh::model::mesh::generate(3);

    if (!msh_file_path.empty())
    {
      ensure_parent_directory(msh_file_path);
      gmsh::write(msh_file_path);
    }

    std::vector<std::size_t> node_tags;
    std::vector<double> coordinates;
    std::vector<double> parametric_coordinates;
    gmsh::model::mesh::getNodes(node_tags, coordinates, parametric_coordinates);

    std::map<std::size_t, int> node_index;
    mesh.nodes.resize(node_tags.size());
    for (std::size_t i = 0; i < node_tags.size(); ++i)
    {
      mesh.nodes[i] = {coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]};
      node_index[node_tags[i]] = static_cast<int>(i);
    }

    std::vector<int> element_types;
    std::vector<std::vector<std::size_t>> element_tags;
    std::vector<std::vector<std::size_t>> element_nodes;
    gmsh::model::mesh::getElements(element_types, element_tags, element_nodes, 3);

    for (std::size_t type_id = 0; type_id < element_types.size(); ++type_id)
    {
      std::string name;
      int dim = 0;
      int order = 0;
      int num_nodes = 0;
      int num_primary_nodes = 0;
      std::vector<double> local_coordinates;
      gmsh::model::mesh::getElementProperties(
          element_types[type_id], name, dim, order, num_nodes, local_coordinates, num_primary_nodes);

      if (dim != 3 || num_primary_nodes != 4)
        continue;

      const auto& connectivity = element_nodes[type_id];
      for (std::size_t i = 0; i + 3 < connectivity.size(); i += 4)
      {
        mesh.tetrahedra.push_back({node_index.at(connectivity[i]),
                                   node_index.at(connectivity[i + 1]),
                                   node_index.at(connectivity[i + 2]),
                                   node_index.at(connectivity[i + 3])});
      }
    }

    gmsh::finalize();
  }
  catch (...)
  {
    gmsh::finalize();
    throw;
  }

  if (mesh.nodes.empty() || mesh.tetrahedra.empty())
    throw std::runtime_error("Gmsh did not generate a tetrahedral mesh.");

  return mesh;
}

void write_particles_vtp(const std::string& file_path,
                         const std::vector<Particle>& particles,
                         double particle_radius)
{
  ensure_parent_directory(file_path);

  auto data = vtkSmartPointer<vtkPolyData>::New();
  auto points = vtkSmartPointer<vtkPoints>::New();
  auto vertices = vtkSmartPointer<vtkCellArray>::New();

  std::vector<double> density;
  std::vector<double> pressure;
  std::vector<double> speed;
  std::vector<double> radius;
  std::vector<Vec3> velocity;

  for (const Particle& particle : particles)
  {
    const vtkIdType id = points->InsertNextPoint(
        particle.position.x, particle.position.y, particle.position.z);
    vertices->InsertNextCell(1, &id);

    density.push_back(particle.density);
    pressure.push_back(particle.pressure);
    velocity.push_back(particle.velocity);
    speed.push_back(norm(particle.velocity));
    radius.push_back(particle_radius);
  }

  data->SetPoints(points);
  data->SetVerts(vertices);
  data->GetPointData()->AddArray(make_scalar_array("density", density));
  data->GetPointData()->AddArray(make_scalar_array("pressure", pressure));
  data->GetPointData()->AddArray(make_vector_array("velocity", velocity));
  data->GetPointData()->AddArray(make_scalar_array("speed", speed));
  data->GetPointData()->AddArray(make_scalar_array("radius", radius));

  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(file_path.c_str());
  writer->SetInputData(data);
  writer->SetDataModeToAscii();
  writer->Write();
}

void write_box_mesh_vtu(const std::string& file_path, const TetMesh& box_mesh, const BoxShell& box)
{
  ensure_parent_directory(file_path);

  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  auto points = vtkSmartPointer<vtkPoints>::New();

  for (const Vec3& point : box_mesh.nodes)
    points->InsertNextPoint(point.x, point.y, point.z);

  grid->SetPoints(points);

  for (const auto& tet : box_mesh.tetrahedra)
  {
    auto vtk_tet = vtkSmartPointer<vtkTetra>::New();
    vtk_tet->GetPointIds()->SetId(0, tet[0]);
    vtk_tet->GetPointIds()->SetId(1, tet[1]);
    vtk_tet->GetPointIds()->SetId(2, tet[2]);
    vtk_tet->GetPointIds()->SetId(3, tet[3]);
    grid->InsertNextCell(vtk_tet->GetCellType(), vtk_tet->GetPointIds());
  }

  const std::size_t count = box_mesh.tetrahedra.size();
  grid->GetCellData()->AddArray(make_cell_scalar_array("wall_thickness", count, box.wall_thickness));
  grid->GetCellData()->AddArray(make_cell_scalar_array("density", count, box.density));
  grid->GetCellData()->AddArray(make_cell_rgb_array("box_rgb", count, 160, 170, 182));

  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(file_path.c_str());
  writer->SetInputData(grid);
  writer->SetDataModeToAscii();
  writer->Write();
}

void write_ground_vtp(const std::string& file_path, const GroundPlane& ground)
{
  ensure_parent_directory(file_path);

  auto data = vtkSmartPointer<vtkPolyData>::New();
  auto points = vtkSmartPointer<vtkPoints>::New();
  auto polygon = vtkSmartPointer<vtkPolygon>::New();
  auto polys = vtkSmartPointer<vtkCellArray>::New();

  const std::vector<Vec3> corners = {
      {-ground.half_size_x, -ground.half_size_y, ground.z},
      {ground.half_size_x, -ground.half_size_y, ground.z},
      {ground.half_size_x, ground.half_size_y, ground.z},
      {-ground.half_size_x, ground.half_size_y, ground.z}};

  polygon->GetPointIds()->SetNumberOfIds(4);
  for (vtkIdType i = 0; i < 4; ++i)
  {
    points->InsertNextPoint(corners[i].x, corners[i].y, corners[i].z);
    polygon->GetPointIds()->SetId(i, i);
  }

  polys->InsertNextCell(polygon);
  data->SetPoints(points);
  data->SetPolys(polys);

  data->GetCellData()->AddArray(make_cell_scalar_array("normal_stiffness", 1, ground.normal_stiffness));
  data->GetCellData()->AddArray(make_cell_scalar_array("max_sink", 1, ground.max_sink));
  data->GetCellData()->AddArray(
      make_cell_rgb_array("ground_rgb", 1, ground.color_rgb[0], ground.color_rgb[1], ground.color_rgb[2]));

  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(file_path.c_str());
  writer->SetInputData(data);
  writer->SetDataModeToAscii();
  writer->Write();
}

void write_pvd_collection(const std::string& file_path, const std::vector<FrameInfo>& frames)
{
  ensure_parent_directory(file_path);

  std::ofstream out(file_path);
  out << std::fixed << std::setprecision(6);
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <Collection>\n";
  for (const FrameInfo& frame : frames)
  {
    out << "    <DataSet timestep=\"" << frame.time << "\" group=\"\" part=\"0\" file=\"" << frame.file_name
        << "\"/>\n";
  }
  out << "  </Collection>\n";
  out << "</VTKFile>\n";
}

void run_mesh_preview()
{
  const std::string output_dir = "results/mesh_preview";
  std::filesystem::remove_all(output_dir);
  std::filesystem::create_directories(output_dir);

  SimulationConfig config;
  const TetMesh mesh = build_box_shell_mesh(config.target, config.box_mesh_size, output_dir + "/box_shell.msh");
  write_box_mesh_vtu(output_dir + "/box_shell.vtu", mesh, config.target);

  std::cout << "Mesh preview saved to " << output_dir << "\n";
}

void run_full_simulation()
{
  const std::string output_dir = "results/latest_run";
  std::filesystem::remove_all(output_dir);
  std::filesystem::create_directories(output_dir + "/water");
  std::filesystem::create_directories(output_dir + "/box");
  std::filesystem::create_directories(output_dir + "/ground");
  std::filesystem::create_directories(output_dir + "/mesh");

  SimulationConfig config;
  const TetMesh box_mesh = build_box_shell_mesh(config.target, config.box_mesh_size, output_dir + "/mesh/box_shell.msh");
  write_box_mesh_vtu(output_dir + "/mesh/box_shell.vtu", box_mesh, config.target);

  SphSimulation simulation(config);
  simulation.initialize_fluid_block();

  std::vector<FrameInfo> water_frames;
  std::ofstream box_loads(output_dir + "/box_loads.csv");

  // FEniCSx потом прочитает этот CSV и превратит удар воды в деформацию сетки бокса.
  box_loads << "frame_id,time,tx,ty,tz,roof_load,roof_impulse,roof_center_x,roof_center_y,"
               "side_x_load,side_y_load,side_x_impulse,side_y_impulse\n";

  double roof_peak_load = 0.0;
  double roof_impulse = 0.0;
  double roof_weighted_x = 0.0;
  double roof_weighted_y = 0.0;
  double roof_weight = 0.0;
  double side_x_peak_load = 0.0;
  double side_y_peak_load = 0.0;
  double side_x_impulse = 0.0;
  double side_y_impulse = 0.0;

  write_ground_vtp(output_dir + "/ground/ground.vtp", simulation.ground());

  int frame_id = 0;
  write_particles_vtp(frame_name(output_dir, "water", ".vtp", frame_id),
                      simulation.particles(),
                      config.visual_particle_radius);
  water_frames.push_back({0.0, relative_frame_name("water", ".vtp", frame_id)});
  box_loads << "0,0,0,0,0,0,0,0,0,0,0,0,0\n";

  for (int step = 0; step < config.total_steps; ++step)
  {
    simulation.step();

    roof_peak_load = std::max(roof_peak_load, simulation.roof_impact_load());
    side_x_peak_load = std::max(side_x_peak_load, simulation.side_x_impact_load());
    side_y_peak_load = std::max(side_y_peak_load, simulation.side_y_impact_load());

    roof_impulse += simulation.roof_impact_impulse();
    side_x_impulse += simulation.side_x_impact_impulse();
    side_y_impulse += simulation.side_y_impact_impulse();

    if (simulation.roof_impact_impulse() > 1e-9)
    {
      roof_weight += simulation.roof_impact_impulse();
      roof_weighted_x += simulation.roof_impact_impulse() * simulation.roof_impact_center_x();
      roof_weighted_y += simulation.roof_impact_impulse() * simulation.roof_impact_center_y();
    }

    if ((step + 1) % config.output_every == 0)
    {
      ++frame_id;
      write_particles_vtp(frame_name(output_dir, "water", ".vtp", frame_id),
                          simulation.particles(),
                          config.visual_particle_radius);
      water_frames.push_back({simulation.current_time(), relative_frame_name("water", ".vtp", frame_id)});

      const Vec3 translation = simulation.box().center - config.target.center;
      const double center_x = roof_weight > 1e-9 ? roof_weighted_x / roof_weight : 0.0;
      const double center_y = roof_weight > 1e-9 ? roof_weighted_y / roof_weight : 0.0;

      box_loads << frame_id << "," << simulation.current_time() << ","
                << translation.x << "," << translation.y << "," << translation.z << ","
                << roof_peak_load << "," << roof_impulse << ","
                << center_x << "," << center_y << ","
                << side_x_peak_load << "," << side_y_peak_load << ","
                << side_x_impulse << "," << side_y_impulse << "\n";

      std::cout << "Saved frame " << frame_id << ", time = " << simulation.current_time() << " s\n";

      roof_peak_load = 0.0;
      roof_impulse = 0.0;
      roof_weighted_x = 0.0;
      roof_weighted_y = 0.0;
      roof_weight = 0.0;
      side_x_peak_load = 0.0;
      side_y_peak_load = 0.0;
      side_x_impulse = 0.0;
      side_y_impulse = 0.0;
    }
  }

  box_loads.close();
  write_pvd_collection(output_dir + "/water_series.pvd", water_frames);

  const char* fenics_python_env = std::getenv("FENICSX_PYTHON");
  const std::string fenics_python = fenics_python_env != nullptr
      ? fenics_python_env
      : "/opt/anaconda3/envs/fenicsx-clean/bin/python";
  const std::string fenics_script = std::filesystem::absolute("src/fenics_box_solver.py").string();
  const std::string fenics_command =
      "env FI_PROVIDER=tcp FI_TCP_IFACE=en0 " + fenics_python + " " + fenics_script + " " + output_dir;

  std::cout << "Running FEniCSx box deformation postprocess...\n";
  if (std::system(fenics_command.c_str()) != 0)
    throw std::runtime_error("FEniCSx postprocess failed.");
}
} // namespace

int main(int argc, char** argv)
{
  try
  {
    if (argc > 1 && std::string(argv[1]) == "--mesh-only")
      run_mesh_preview();
    else
      run_full_simulation();
  }
  catch (const std::exception& error)
  {
    std::cerr << "Simulation error: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
