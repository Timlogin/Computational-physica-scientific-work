#pragma once

#include <string>
#include <vector>

struct Vec3
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

Vec3 operator+(const Vec3& a, const Vec3& b);
Vec3 operator-(const Vec3& a, const Vec3& b);
Vec3 operator*(const Vec3& v, double s);
Vec3 operator/(const Vec3& v, double s);
Vec3& operator+=(Vec3& a, const Vec3& b);
Vec3& operator-=(Vec3& a, const Vec3& b);
Vec3& operator*=(Vec3& a, double s);

double dot(const Vec3& a, const Vec3& b);
double norm(const Vec3& v);

struct Particle
{
  Vec3 position;
  Vec3 velocity;
  Vec3 acceleration;
  double density = 0.0;
  double pressure = 0.0;
};

struct BoxShell
{
  Vec3 center;
  Vec3 outer_size;
  double wall_thickness = 0.02;
  double density = 7800.0;
};

struct SimulationConfig
{
  double dt = 0.008;
  int total_steps = 1000;
  int output_every = 5;

  double particle_spacing = 0.28;
  double smoothing_length = 0.44;
  double total_water_mass = 2000.0;
  double particle_mass = 1.0;
  double rest_density = 1000.0;
  double pressure_stiffness = 30.0;
  double viscosity = 18.0;
  double air_drag = 0.18;
  double visual_particle_radius = 0.28;

  Vec3 gravity{0.0, 0.0, -9.81};

  Vec3 fluid_block_min{-1.4, -0.7, 53.5};
  Vec3 fluid_block_max{1.4, 0.7, 55.0};

  BoxShell target{{0.0, 0.0, 1.0}, {4.0, 2.0, 1.8}, 0.06, 550.0};
  double box_support_stiffness = 90000.0;
  double box_support_damping = 24000.0;

  Vec3 world_min{-6.0, -4.0, -0.5};
  Vec3 world_max{6.0, 4.0, 60.0};
  double collision_damping = 0.35;
};

class SphSimulation
{
public:
  explicit SphSimulation(SimulationConfig config);

  void initialize_fluid_block();
  void step();

  const SimulationConfig& config() const { return _config; }
  const std::vector<Particle>& particles() const { return _particles; }
  const BoxShell& box() const { return _box; }
  const Vec3& box_velocity() const { return _box_velocity; }
  double box_mass() const { return _box_mass; }
  int current_step() const { return _current_step; }
  double current_time() const { return _current_time; }

private:
  void compute_density_and_pressure();
  void compute_forces();
  void integrate_particles();
  void integrate_box();
  void resolve_world_collision(Particle& particle);
  void resolve_box_collision(Particle& particle, const Vec3& previous_position);
  double compute_box_mass() const;

  double poly6_kernel(double distance) const;
  double spiky_gradient(double distance) const;
  double viscosity_laplacian(double distance) const;

  SimulationConfig _config;
  std::vector<Particle> _particles;
  BoxShell _box;
  Vec3 _box_velocity{};
  Vec3 _box_acceleration{};
  Vec3 _box_force{};
  Vec3 _box_reference_center{};
  double _box_mass = 1.0;
  int _current_step = 0;
  double _current_time = 0.0;
};
