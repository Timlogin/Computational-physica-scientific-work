#include "vtk_writer.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace
{
void ensure_parent_directory(const std::string& file_path)
{
  const std::filesystem::path path(file_path);
  if (path.has_parent_path())
    std::filesystem::create_directories(path.parent_path());
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

void write_box_vtp(const std::string& file_path,
                   const BoxShell& box,
                   const Vec3& box_velocity,
                   double box_mass)
{
  ensure_parent_directory(file_path);

  std::ofstream out(file_path);
  if (!out)
    throw std::runtime_error("Could not open VTP file for box: " + file_path);

  const Vec3 h = box.outer_size * 0.5;

  const std::vector<Vec3> points = {
      {box.center.x - h.x, box.center.y - h.y, box.center.z - h.z},
      {box.center.x + h.x, box.center.y - h.y, box.center.z - h.z},
      {box.center.x + h.x, box.center.y + h.y, box.center.z - h.z},
      {box.center.x - h.x, box.center.y + h.y, box.center.z - h.z},
      {box.center.x - h.x, box.center.y - h.y, box.center.z + h.z},
      {box.center.x + h.x, box.center.y - h.y, box.center.z + h.z},
      {box.center.x + h.x, box.center.y + h.y, box.center.z + h.z},
      {box.center.x - h.x, box.center.y + h.y, box.center.z + h.z}};

  const int faces[6][4] = {
      {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
      {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};

  out << std::fixed << std::setprecision(6);
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <PolyData>\n";
  out << "    <Piece NumberOfPoints=\"8\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" "
         "NumberOfPolys=\"6\">\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const Vec3& point : points)
    out << "          " << point.x << " " << point.y << " " << point.z << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";

  out << "      <CellData Scalars=\"wall_thickness\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"wall_thickness\" format=\"ascii\">\n";
  for (int i = 0; i < 6; ++i)
    out << "          " << box.wall_thickness << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float64\" Name=\"density\" format=\"ascii\">\n";
  for (int i = 0; i < 6; ++i)
    out << "          " << box.density << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float64\" Name=\"mass\" format=\"ascii\">\n";
  for (int i = 0; i < 6; ++i)
    out << "          " << box_mass << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float64\" Name=\"speed\" format=\"ascii\">\n";
  for (int i = 0; i < 6; ++i)
    out << "          " << norm(box_velocity) << "\n";
  out << "        </DataArray>\n";
  out << "      </CellData>\n";

  out << "      <Polys>\n";
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (const auto& face : faces)
    out << "          " << face[0] << " " << face[1] << " " << face[2] << " " << face[3] << "\n";
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (int i = 1; i <= 6; ++i)
    out << "          " << i * 4 << "\n";
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
