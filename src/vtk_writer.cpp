#include "vtk_writer.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace
{
// Создаёт родительскую директорию перед записью файла
void ensure_parent_directory(const std::string& file_path)
{
  const std::filesystem::path path(file_path);
  if (path.has_parent_path())
    std::filesystem::create_directories(path.parent_path());
}

// Центр сетки нужен, чтобы деформацию прикладывать в локальной системе координат бокса
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

// Простая упругая деформация: сдвигаем узлы внутрь около стенок,
// вниз около крыши и вверх около днища.
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
}

void write_particles_vtp(const std::string& file_path,
                         const std::vector<Particle>& particles,
                         double particle_radius)
{
  ensure_parent_directory(file_path);

  std::ofstream out(file_path);
  if (!out)
    throw std::runtime_error("Could not open VTP file for particles: " + file_path);

  // Основной XML-заголовок VTK
  out << std::fixed << std::setprecision(6);
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <PolyData>\n";
  out << "    <Piece NumberOfPoints=\"" << particles.size() << "\" NumberOfVerts=\"" << particles.size()
      << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const Particle& particle : particles)
    out << "          " << particle.position.x << " " << particle.position.y << " " << particle.position.z
        << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";

  // Сохраняем основные физические поля каждой частицы
  out << "      <PointData Scalars=\"density\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"density\" format=\"ascii\">\n";
  for (const Particle& particle : particles)
    out << "          " << particle.density << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\">\n";
  for (const Particle& particle : particles)
    out << "          " << particle.pressure << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const Particle& particle : particles)
    out << "          " << particle.velocity.x << " " << particle.velocity.y << " " << particle.velocity.z
        << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float64\" Name=\"speed\" format=\"ascii\">\n";
  for (const Particle& particle : particles)
  {
    const double speed = norm(particle.velocity);
    out << "          " << speed << "\n";
  }
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float64\" Name=\"radius\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < particles.size(); ++i)
    out << "          " << particle_radius << "\n";
  out << "        </DataArray>\n";
  out << "      </PointData>\n";

  out << "      <Verts>\n";
  // Каждая частица считается отдельной вершиной
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < particles.size(); ++i)
    out << "          " << i << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < particles.size(); ++i)
    out << "          " << i + 1 << "\n";
  out << "        </DataArray>\n";
  out << "      </Verts>\n";

  out << "    </Piece>\n";
  out << "  </PolyData>\n";
  out << "</VTKFile>\n";
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

  std::ofstream out(file_path);
  if (!out)
    throw std::runtime_error("Could not open VTU file for box: " + file_path);

  const std::size_t num_cells = box_mesh.tetrahedra.size();
  const Vec3 reference_center = compute_mesh_center(box_mesh);

  out << std::fixed << std::setprecision(6);
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <UnstructuredGrid>\n";
  out << "    <Piece NumberOfPoints=\"" << box_mesh.nodes.size() << "\" NumberOfCells=\"" << num_cells
      << "\">\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const Vec3& point : box_mesh.nodes)
  {
    // Сначала деформируем узел в локальной системе бокса,
    // затем переносим всю геометрию в текущее положение.
    const Vec3 deformed = deform_box_point(point, reference_center, box.outer_size, deformation);
    const Vec3 shifted = deformed + translation;
    out << "          " << shifted.x << " " << shifted.y << " " << shifted.z << "\n";
  }
  out << "        </DataArray>\n";
  out << "      </Points>\n";

  out << "      <CellData Scalars=\"wall_thickness\">\n";
  // Параметры материала одинаковы для всех ячеек бокса.
  out << "        <DataArray type=\"Float64\" Name=\"wall_thickness\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < num_cells; ++i)
    out << "          " << box.wall_thickness << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float64\" Name=\"density\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < num_cells; ++i)
    out << "          " << box.density << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float64\" Name=\"mass\" format=\"ascii\">\n";
  const double cell_mass = num_cells == 0 ? 0.0 : box_mass / static_cast<double>(num_cells);
  for (std::size_t i = 0; i < num_cells; ++i)
    out << "          " << cell_mass << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float64\" Name=\"speed\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < num_cells; ++i)
    out << "          " << norm(box_velocity) << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Float64\" Name=\"roof_deflection\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < num_cells; ++i)
    out << "          " << deformation.roof << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Float64\" Name=\"roof_center_x\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < num_cells; ++i)
    out << "          " << deformation.roof_center_x << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Float64\" Name=\"roof_center_y\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < num_cells; ++i)
    out << "          " << deformation.roof_center_y << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Float64\" Name=\"floor_deflection\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < num_cells; ++i)
    out << "          " << deformation.floor << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Float64\" Name=\"side_deflection\" format=\"ascii\">\n";
  const double side_deflection = 0.5 * (deformation.side_x + deformation.side_y);
  for (std::size_t i = 0; i < num_cells; ++i)
    out << "          " << side_deflection << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"UInt8\" Name=\"box_rgb\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < num_cells; ++i)
    out << "          160 170 182\n";
  out << "        </DataArray>\n";
  out << "      </CellData>\n";

  // Ниже задаём сами тетраэдры через вершины.
  out << "      <Cells>\n";
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (const auto& tet : box_mesh.tetrahedra)
    out << "          " << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (std::size_t i = 1; i <= num_cells; ++i)
    out << "          " << i * 4 << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < num_cells; ++i)
    out << "          10\n";
  out << "        </DataArray>\n";
  out << "      </Cells>\n";

  out << "      <FieldData>\n";
  out << "        <DataArray type=\"Float64\" Name=\"box_total_mass\" NumberOfTuples=\"1\" format=\"ascii\">\n";
  out << "          " << box_mass << "\n";
  out << "        </DataArray>\n";
  out << "      </FieldData>\n";

  out << "    </Piece>\n";
  out << "  </UnstructuredGrid>\n";
  out << "</VTKFile>\n";
}

void write_ground_vtp(const std::string& file_path, const GroundPlane& ground)
{
  ensure_parent_directory(file_path);

  std::ofstream out(file_path);
  if (!out)
    throw std::runtime_error("Could not open VTP file for ground: " + file_path);

  const std::vector<Vec3> points = {
      {-ground.half_size_x, -ground.half_size_y, ground.z},
      {ground.half_size_x, -ground.half_size_y, ground.z},
      {ground.half_size_x, ground.half_size_y, ground.z},
      {-ground.half_size_x, ground.half_size_y, ground.z}};

  out << std::fixed << std::setprecision(6);
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <PolyData>\n";
  out << "    <Piece NumberOfPoints=\"4\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"1\">\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const Vec3& point : points)
    out << "          " << point.x << " " << point.y << " " << point.z << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";

  // Храним параметры контакта и цвет земли
  out << "      <CellData Scalars=\"normal_stiffness\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"normal_stiffness\" format=\"ascii\">\n";
  out << "          " << ground.normal_stiffness << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Float64\" Name=\"lateral_stiffness\" format=\"ascii\">\n";
  out << "          " << ground.lateral_stiffness << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Float64\" Name=\"max_sink\" format=\"ascii\">\n";
  out << "          " << ground.max_sink << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"UInt8\" Name=\"ground_rgb\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  out << "          " << ground.color_rgb[0] << " " << ground.color_rgb[1] << " " << ground.color_rgb[2]
      << "\n";
  out << "        </DataArray>\n";
  out << "      </CellData>\n";

  out << "      <Polys>\n";
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  out << "          0 1 2 3\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  out << "          4\n";
  out << "        </DataArray>\n";
  out << "      </Polys>\n";
  out << "    </Piece>\n";
  out << "  </PolyData>\n";
  out << "</VTKFile>\n";
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
