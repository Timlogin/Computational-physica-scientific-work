#pragma once

#include "simulation.hpp"

#include <string>
#include <vector>

struct FrameInfo
{
  double time = 0.0;
  std::string file_name;
};

void write_particles_vtp(const std::string& file_path,
                         const std::vector<Particle>& particles,
                         double particle_radius);
void write_box_vtp(const std::string& file_path,
                   const BoxShell& box,
                   const Vec3& box_velocity,
                   double box_mass);
void write_pvd_collection(const std::string& file_path, const std::vector<FrameInfo>& frames);
