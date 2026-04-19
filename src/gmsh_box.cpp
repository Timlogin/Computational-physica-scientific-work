#include "gmsh_box.hpp"

#include <gmsh.h>

#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
// Векторное произведение нужно для вычисления объёма тетраэдра.
Vec3 cross(const Vec3& a, const Vec3& b)
{
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
}

TetMesh build_box_shell_mesh(const BoxShell& box, double mesh_size, const std::string& msh_file_path)
{
  // Внутренний бокс нужен, чтобы вырезать полость из внешнего объёма.
  const Vec3 inner_size = {box.outer_size.x - 2.0 * box.wall_thickness,
                           box.outer_size.y - 2.0 * box.wall_thickness,
                           box.outer_size.z - 2.0 * box.wall_thickness};

  // Если толщина стенки слишком велика, полость построить уже нельзя.
  if (inner_size.x <= 0.0 || inner_size.y <= 0.0 || inner_size.z <= 0.0)
    throw std::runtime_error("Box wall thickness is too large for the selected outer size.");

  TetMesh mesh;

  // Gmsh инициализируем только на время построения сетки.
  gmsh::initialize();

  try
  {
    gmsh::option::setNumber("General.Terminal", 0.0);
    gmsh::model::add("box_shell");

    // Нижние углы внешнего и внутреннего боксов.
    const Vec3 outer_min = box.center - box.outer_size * 0.5;
    const Vec3 inner_min = outer_min + Vec3{box.wall_thickness, box.wall_thickness, box.wall_thickness};

    // Строим два параллелепипеда через OpenCASCADE
    const int outer_box_tag = gmsh::model::occ::addBox(
        outer_min.x, outer_min.y, outer_min.z, box.outer_size.x, box.outer_size.y, box.outer_size.z);
    const int inner_box_tag = gmsh::model::occ::addBox(
        inner_min.x, inner_min.y, inner_min.z, inner_size.x, inner_size.y, inner_size.z);

    std::vector<std::pair<int, int>> outer_box = {{3, outer_box_tag}};
    std::vector<std::pair<int, int>> inner_box = {{3, inner_box_tag}};
    std::vector<std::pair<int, int>> shell_box;

    std::vector<std::vector<std::pair<int, int>>> shell_parts;
    // Логика такая же, как в torus:
    gmsh::model::occ::cut(outer_box, inner_box, shell_box, shell_parts);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMin", mesh_size);
    gmsh::option::setNumber("Mesh.MeshSizeMax", mesh_size);
    gmsh::option::setNumber("Mesh.Algorithm3D", 1.0);
    gmsh::option::setNumber("Mesh.Optimize", 1.0);
    gmsh::model::mesh::generate(3);

    if (!msh_file_path.empty())
    {
      const std::filesystem::path path(msh_file_path);
      if (path.has_parent_path())
        std::filesystem::create_directories(path.parent_path());
      gmsh::write(msh_file_path);
    }

    std::vector<std::size_t> node_tags;
    std::vector<double> coordinates;
    std::vector<double> parametric_coordinates;
    gmsh::model::mesh::getNodes(node_tags, coordinates, parametric_coordinates);
    mesh.nodes.resize(node_tags.size());

    std::map<std::size_t, int> node_index;
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
      int dimension = 0;
      int order = 0;
      int num_nodes = 0;
      std::vector<double> local_coordinates;
      int num_primary_nodes = 0;
      gmsh::model::mesh::getElementProperties(
          element_types[type_id], name, dimension, order, num_nodes, local_coordinates, num_primary_nodes);

      // Нас интересуют только линейные тетраэдры с 4 вершинами
      if (dimension != 3 || num_primary_nodes != 4)
        continue;

      const std::vector<std::size_t>& connectivity = element_nodes[type_id];
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
    throw std::runtime_error("Gmsh did not generate a tetrahedral box mesh.");

  return mesh;
}

double compute_tet_mesh_volume(const TetMesh& mesh)
{
  double volume = 0.0;

  for (const auto& tet : mesh.tetrahedra)
  {
    const Vec3& a = mesh.nodes[tet[0]];
    const Vec3& b = mesh.nodes[tet[1]];
    const Vec3& c = mesh.nodes[tet[2]];
    const Vec3& d = mesh.nodes[tet[3]];

    volume += std::abs(dot(a - d, cross(b - d, c - d))) / 6.0;
  }

  return volume;
}
