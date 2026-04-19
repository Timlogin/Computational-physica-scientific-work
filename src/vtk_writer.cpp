#include "vtk_writer.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
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

namespace
{
// Создаёт родительскую директорию перед записью файла.
void ensure_parent_directory(const std::string& file_path)
{
  const std::filesystem::path path(file_path);
  if (path.has_parent_path())
    std::filesystem::create_directories(path.parent_path());
}

// Центр сетки нужен, чтобы прикладывать простую деформацию в локальной системе координат бокса.
Vec3 compute_mesh_center(const TetMesh& box_mesh)
{
  if (box_mesh.nodes.empty())
    return {};

  Vec3 min_point = box_mesh.nodes.front();
  Vec3 max_point = box_mesh.nodes.front();

  for (const Vec3& point : box_mesh.nodes)
  {
    min_point.x = std::min(min_point.x, point.x);
    min_point.y = std::min(min_point.y, point.y);
    min_point.z = std::min(min_point.z, point.z);
    max_point.x = std::max(max_point.x, point.x);
    max_point.y = std::max(max_point.y, point.y);
    max_point.z = std::max(max_point.z, point.z);
  }

  return {(min_point.x + max_point.x) * 0.5,
          (min_point.y + max_point.y) * 0.5,
          (min_point.z + max_point.z) * 0.5};
}

// Простая упругая деформация для режима preview:
// стенки чуть втягиваются внутрь, крыша мнётся вниз, днище немного выгибается вверх.
Vec3 deform_box_point(const Vec3& point,
                      const Vec3& reference_center,
                      const Vec3& outer_size,
                      const BoxElasticDeformation& deformation)
{
  Vec3 local = point - reference_center;
  const Vec3 half_outer = outer_size * 0.5;
  const double roof_patch_x = std::max(0.35, 0.22 * outer_size.x);
  const double roof_patch_y = std::max(0.25, 0.28 * outer_size.y);

  const double wx = half_outer.x > 1e-9 ? std::abs(local.x) / half_outer.x : 0.0;
  const double wy = half_outer.y > 1e-9 ? std::abs(local.y) / half_outer.y : 0.0;

  if (local.x > 0.0)
    local.x -= deformation.side_x * wx * wx;
  else
    local.x += deformation.side_x * wx * wx;

  if (local.y > 0.0)
    local.y -= deformation.side_y * wy * wy;
  else
    local.y += deformation.side_y * wy * wy;

  if (local.z > 0.0)
  {
    const double roof_weight = half_outer.z > 1e-9 ? local.z / half_outer.z : 0.0;
    const double dx = local.x - deformation.roof_center_x;
    const double dy = local.y - deformation.roof_center_y;
    const double roof_shape = std::exp(
        -(dx * dx) / (2.0 * roof_patch_x * roof_patch_x)
        -(dy * dy) / (2.0 * roof_patch_y * roof_patch_y));
    local.z -= deformation.roof * roof_weight * roof_weight * roof_shape;
  }
  else
  {
    const double floor_weight = half_outer.z > 1e-9 ? std::abs(local.z) / half_outer.z : 0.0;
    local.z += deformation.floor * floor_weight * floor_weight;
  }

  return reference_center + local;
}

// Добавляет скалярный массив в PointData.
vtkSmartPointer<vtkDoubleArray> make_point_scalar_array(
    const char* name, const std::vector<double>& values)
{
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetName(name);
  for (double value : values)
    array->InsertNextValue(value);
  return array;
}

// Добавляет векторный массив в PointData.
vtkSmartPointer<vtkDoubleArray> make_point_vector_array(
    const char* name, const std::vector<Vec3>& values)
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

// Создаёт скалярный CellData-массив, одинаковый для всех ячеек.
vtkSmartPointer<vtkDoubleArray> make_cell_scalar_array(
    const char* name, std::size_t count, double value)
{
  auto array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetName(name);
  for (std::size_t i = 0; i < count; ++i)
    array->InsertNextValue(value);
  return array;
}

// Создаёт RGB-массив для ячеек.
vtkSmartPointer<vtkUnsignedCharArray> make_cell_rgb_array(
    const char* name, std::size_t count, int r, int g, int b)
{
  auto array = vtkSmartPointer<vtkUnsignedCharArray>::New();
  array->SetName(name);
  array->SetNumberOfComponents(3);

  const unsigned char tuple[3] = {
      static_cast<unsigned char>(r),
      static_cast<unsigned char>(g),
      static_cast<unsigned char>(b)};

  for (std::size_t i = 0; i < count; ++i)
    array->InsertNextTypedTuple(tuple);

  return array;
}
} // namespace

void write_particles_vtp(const std::string& file_path,
                         const std::vector<Particle>& particles,
                         double particle_radius)
{
  ensure_parent_directory(file_path);

  auto poly_data = vtkSmartPointer<vtkPolyData>::New();
  auto points = vtkSmartPointer<vtkPoints>::New();
  auto verts = vtkSmartPointer<vtkCellArray>::New();

  std::vector<double> density;
  std::vector<double> pressure;
  std::vector<double> speed;
  std::vector<double> radius;
  std::vector<Vec3> velocity;

  density.reserve(particles.size());
  pressure.reserve(particles.size());
  speed.reserve(particles.size());
  radius.reserve(particles.size());
  velocity.reserve(particles.size());

  for (std::size_t i = 0; i < particles.size(); ++i)
  {
    const Particle& particle = particles[i];
    points->InsertNextPoint(particle.position.x, particle.position.y, particle.position.z);
    verts->InsertNextCell(1);
    verts->InsertCellPoint(static_cast<vtkIdType>(i));

    density.push_back(particle.density);
    pressure.push_back(particle.pressure);
    velocity.push_back(particle.velocity);
    speed.push_back(norm(particle.velocity));
    radius.push_back(particle_radius);
  }

  poly_data->SetPoints(points);
  poly_data->SetVerts(verts);
  poly_data->GetPointData()->AddArray(make_point_scalar_array("density", density));
  poly_data->GetPointData()->AddArray(make_point_scalar_array("pressure", pressure));
  poly_data->GetPointData()->AddArray(make_point_vector_array("velocity", velocity));
  poly_data->GetPointData()->AddArray(make_point_scalar_array("speed", speed));
  poly_data->GetPointData()->AddArray(make_point_scalar_array("radius", radius));

  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(file_path.c_str());
  writer->SetInputData(poly_data);
  writer->SetDataModeToAscii();
  if (writer->Write() == 0)
    throw std::runtime_error("Could not write VTP file for particles: " + file_path);
}

void write_box_vtu(const std::string& file_path,
                   const TetMesh& box_mesh,
                   const Vec3& translation,
                   const BoxShell& box,
                   const Vec3& box_velocity,
                   double box_mass,
                   const BoxElasticDeformation& deformation)
{
  ensure_parent_directory(file_path);

  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  auto points = vtkSmartPointer<vtkPoints>::New();

  const Vec3 reference_center = compute_mesh_center(box_mesh);

  for (const Vec3& point : box_mesh.nodes)
  {
    const Vec3 deformed = deform_box_point(point, reference_center, box.outer_size, deformation);
    const Vec3 shifted = deformed + translation;
    points->InsertNextPoint(shifted.x, shifted.y, shifted.z);
  }

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

  const std::size_t num_cells = box_mesh.tetrahedra.size();
  const double side_deflection = 0.5 * (deformation.side_x + deformation.side_y);
  const double cell_mass = num_cells == 0 ? 0.0 : box_mass / static_cast<double>(num_cells);

  grid->GetCellData()->AddArray(make_cell_scalar_array("wall_thickness", num_cells, box.wall_thickness));
  grid->GetCellData()->AddArray(make_cell_scalar_array("density", num_cells, box.density));
  grid->GetCellData()->AddArray(make_cell_scalar_array("mass", num_cells, cell_mass));
  grid->GetCellData()->AddArray(make_cell_scalar_array("speed", num_cells, norm(box_velocity)));
  grid->GetCellData()->AddArray(make_cell_scalar_array("roof_deflection", num_cells, deformation.roof));
  grid->GetCellData()->AddArray(make_cell_scalar_array("roof_center_x", num_cells, deformation.roof_center_x));
  grid->GetCellData()->AddArray(make_cell_scalar_array("roof_center_y", num_cells, deformation.roof_center_y));
  grid->GetCellData()->AddArray(make_cell_scalar_array("floor_deflection", num_cells, deformation.floor));
  grid->GetCellData()->AddArray(make_cell_scalar_array("side_deflection", num_cells, side_deflection));
  grid->GetCellData()->AddArray(make_cell_rgb_array("box_rgb", num_cells, 160, 170, 182));

  auto total_mass = vtkSmartPointer<vtkDoubleArray>::New();
  total_mass->SetName("box_total_mass");
  total_mass->InsertNextValue(box_mass);
  grid->GetFieldData()->AddArray(total_mass);

  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(file_path.c_str());
  writer->SetInputData(grid);
  writer->SetDataModeToAscii();
  if (writer->Write() == 0)
    throw std::runtime_error("Could not write VTU file for box: " + file_path);
}

void write_ground_vtp(const std::string& file_path, const GroundPlane& ground)
{
  ensure_parent_directory(file_path);

  auto poly_data = vtkSmartPointer<vtkPolyData>::New();
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
  poly_data->SetPoints(points);
  poly_data->SetPolys(polys);

  poly_data->GetCellData()->AddArray(make_cell_scalar_array("normal_stiffness", 1, ground.normal_stiffness));
  poly_data->GetCellData()->AddArray(make_cell_scalar_array("lateral_stiffness", 1, ground.lateral_stiffness));
  poly_data->GetCellData()->AddArray(make_cell_scalar_array("max_sink", 1, ground.max_sink));
  poly_data->GetCellData()->AddArray(
      make_cell_rgb_array("ground_rgb", 1, ground.color_rgb[0], ground.color_rgb[1], ground.color_rgb[2]));

  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(file_path.c_str());
  writer->SetInputData(poly_data);
  writer->SetDataModeToAscii();
  if (writer->Write() == 0)
    throw std::runtime_error("Could not write VTP file for ground: " + file_path);
}

void write_pvd_collection(const std::string& file_path, const std::vector<FrameInfo>& frames)
{
  ensure_parent_directory(file_path);

  std::ofstream out(file_path);
  if (!out)
    throw std::runtime_error("Could not open PVD file: " + file_path);

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
