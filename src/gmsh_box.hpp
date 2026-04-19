#pragma once

#include "simulation.hpp"

#include <array>
#include <string>
#include <vector>

// Cтруктура для хранения узлов и тетраэдров сетки в памяти
struct TetMesh
{
  std::vector<Vec3> nodes;                 // Координаты всех узлов сетки.
  std::vector<std::array<int, 4>> tetrahedra; // Связность тетраэдров по индексам узлов.
};

// Строит тетраэдральную сетку тонкостенного бокса.
// Если указан путь msh_file_path, Gmsh дополнительно сохраняет исходную сетку в .msh.
TetMesh build_box_shell_mesh(const BoxShell& box, double mesh_size, const std::string& msh_file_path = "");

// Суммарный объём сетки как сумма объёмов всех тетраэдров.
double compute_tet_mesh_volume(const TetMesh& mesh);
