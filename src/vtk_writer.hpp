#pragma once

#include "gmsh_box.hpp"
#include "simulation.hpp"

#include <string>
#include <vector>

struct FrameInfo
{
  double time = 0.0;      // Время, которому соответствует кадр.
  std::string file_name;  // Относительный путь до файла этого кадра внутри серии
};

// Записывает облако частиц воды в VTP
void write_particles_vtp(const std::string& file_path,
                         const std::vector<Particle>& particles,
                         double particle_radius);

// Записывает тетраэдральную сетку бокса в VTU
void write_box_vtu(const std::string& file_path,
                   const TetMesh& box_mesh,
                   const Vec3& translation,
                   const BoxShell& box,
                   const Vec3& box_velocity,
                   double box_mass,
                   const BoxElasticDeformation& deformation);

// Записывает плоскость земли в VTP
void write_ground_vtp(const std::string& file_path, const GroundPlane& ground);

// Записывает PVD-коллекцию кадров для ParaView
void write_pvd_collection(const std::string& file_path, const std::vector<FrameInfo>& frames);
