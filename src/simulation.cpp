#include "simulation.hpp"

#include <algorithm>
#include <cmath>

namespace
{
constexpr double kPi = 3.14159265358979323846;

double clamp(double value, double low, double high)
{
  return std::max(low, std::min(value, high));
}
}

Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 operator*(const Vec3& v, double s) { return {v.x * s, v.y * s, v.z * s}; }
Vec3 operator/(const Vec3& v, double s) { return {v.x / s, v.y / s, v.z / s}; }
Vec3& operator+=(Vec3& a, const Vec3& b)
{
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  return a;
}
Vec3& operator-=(Vec3& a, const Vec3& b)
{
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
  return a;
}
Vec3& operator*=(Vec3& a, double s)
{
  a.x *= s;
  a.y *= s;
  a.z *= s;
  return a;
}

double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }

SphSimulation::SphSimulation(SimulationConfig config)
    : _config(config), _box(config.target), _box_reference_center(config.target.center)
{
  _box_mass = compute_box_mass();
}

void SphSimulation::initialize_fluid_block()
{
  _particles.clear();

  const double h = _config.particle_spacing;
  for (double x = _config.fluid_block_min.x; x <= _config.fluid_block_max.x + 1e-9; x += h)
  {
    for (double y = _config.fluid_block_min.y; y <= _config.fluid_block_max.y + 1e-9; y += h)
    {
      for (double z = _config.fluid_block_min.z; z <= _config.fluid_block_max.z + 1e-9; z += h)
      {
        Particle particle;
        particle.position = {x, y, z};
        _particles.push_back(particle);
      }
    }
  }

  if (!_particles.empty())
    _config.particle_mass = _config.total_water_mass / static_cast<double>(_particles.size());

  // Сразу считаем начальные плотности и давления, чтобы нулевой кадр
  // не содержал пустых значений в визуализации.
  compute_density_and_pressure();
}

void SphSimulation::step()
{
  _box_force = {};
  compute_density_and_pressure();
  compute_forces();
  integrate_particles();
  integrate_box();
  ++_current_step;
  _current_time += _config.dt;
}

void SphSimulation::compute_density_and_pressure()
{
  const double self_kernel = poly6_kernel(0.0);

  for (Particle& pi : _particles)
  {
    double density = _config.particle_mass * self_kernel;
    for (const Particle& pj : _particles)
    {
      const Vec3 r = pi.position - pj.position;
      const double distance = norm(r);
      if (distance <= 1e-9)
        continue;
      density += _config.particle_mass * poly6_kernel(distance);
    }

    pi.density = std::max(density, _config.rest_density * 0.5);
    pi.pressure = _config.pressure_stiffness * (pi.density - _config.rest_density);
  }
}

void SphSimulation::compute_forces()
{
  for (Particle& pi : _particles)
  {
    Vec3 pressure_force{};
    Vec3 viscosity_force{};

    for (const Particle& pj : _particles)
    {
      const Vec3 r = pi.position - pj.position;
      const double distance = norm(r);
      if (distance <= 1e-9 || distance >= _config.smoothing_length)
        continue;

      const Vec3 direction = r / distance;
      const double avg_pressure = (pi.pressure + pj.pressure) / 2.0;
      const double pressure_term
          = -_config.particle_mass * avg_pressure * spiky_gradient(distance) / pj.density;
      pressure_force += direction * pressure_term;

      const Vec3 velocity_diff = pj.velocity - pi.velocity;
      const double viscosity_term
          = _config.viscosity * _config.particle_mass * viscosity_laplacian(distance) / pj.density;
      viscosity_force += velocity_diff * viscosity_term;
    }

    const Vec3 drag_force = pi.velocity * (-_config.air_drag);
    pi.acceleration = _config.gravity + (pressure_force + viscosity_force + drag_force) / pi.density;
  }
}

void SphSimulation::integrate_particles()
{
  for (Particle& particle : _particles)
  {
    const Vec3 previous_position = particle.position;
    particle.velocity += particle.acceleration * _config.dt;
    particle.position += particle.velocity * _config.dt;

    resolve_box_collision(particle, previous_position);
    resolve_world_collision(particle);
  }
}

void SphSimulation::integrate_box()
{
  const Vec3 spring_force = (_box_reference_center - _box.center) * _config.box_support_stiffness;
  const Vec3 damping_force = _box_velocity * (-_config.box_support_damping);
  const Vec3 total_force = _box_force + spring_force + damping_force;

  _box_acceleration = total_force / _box_mass;
  _box_velocity += _box_acceleration * _config.dt;
  _box.center += _box_velocity * _config.dt;

  // Пока бокс моделируется как одно твёрдое тело, а не как деформируемая оболочка.
  // Поэтому ограничим его общий сдвиг, чтобы визуализация оставалась разумной.
  const double max_shift_xy = 0.35;
  const double max_shift_down = 0.35;
  const double max_shift_up = 0.10;

  const double min_x = _box_reference_center.x - max_shift_xy;
  const double max_x = _box_reference_center.x + max_shift_xy;
  const double min_y = _box_reference_center.y - max_shift_xy;
  const double max_y = _box_reference_center.y + max_shift_xy;
  const double min_z = _box_reference_center.z - max_shift_down;
  const double max_z = _box_reference_center.z + max_shift_up;

  const double clamped_x = clamp(_box.center.x, min_x, max_x);
  const double clamped_y = clamp(_box.center.y, min_y, max_y);
  const double clamped_z = clamp(_box.center.z, min_z, max_z);

  if (clamped_x != _box.center.x)
    _box_velocity.x = 0.0;
  if (clamped_y != _box.center.y)
    _box_velocity.y = 0.0;
  if (clamped_z != _box.center.z)
    _box_velocity.z = 0.0;

  _box.center = {clamped_x, clamped_y, clamped_z};
}

void SphSimulation::resolve_world_collision(Particle& particle)
{
  auto apply_plane = [&](double& coord, double& velocity, double low, double high)
  {
    if (coord < low)
    {
      coord = low;
      if (velocity < 0.0)
        velocity *= -_config.collision_damping;
    }
    else if (coord > high)
    {
      coord = high;
      if (velocity > 0.0)
        velocity *= -_config.collision_damping;
    }
  };

  apply_plane(particle.position.x, particle.velocity.x, _config.world_min.x, _config.world_max.x);
  apply_plane(particle.position.y, particle.velocity.y, _config.world_min.y, _config.world_max.y);
  apply_plane(particle.position.z, particle.velocity.z, _config.world_min.z, _config.world_max.z);
}

void SphSimulation::resolve_box_collision(Particle& particle, const Vec3& previous_position)
{
  const Vec3 half_outer = _box.outer_size * 0.5;
  const double surface_offset = 0.02;
  const Vec3 local = particle.position - _box.center;
  const Vec3 previous_local = previous_position - _box.center;
  const bool inside_outer = std::abs(local.x) <= half_outer.x && std::abs(local.y) <= half_outer.y
                            && std::abs(local.z) <= half_outer.z;

  // Пока бокс не деформируется как оболочка, считаем его замкнутым внешним объёмом.
  // Это заметно лучше удерживает воду снаружи и убирает "пролетание" сквозь стенки.
  if (!inside_outer)
    return;

  const double dx = half_outer.x - std::abs(local.x);
  const double dy = half_outer.y - std::abs(local.y);
  const double dz = half_outer.z - std::abs(local.z);

  int axis = -1;
  double sign = 1.0;

  if (std::abs(previous_local.x) > half_outer.x)
  {
    axis = 0;
    sign = previous_local.x > 0.0 ? 1.0 : -1.0;
  }
  else if (std::abs(previous_local.y) > half_outer.y)
  {
    axis = 1;
    sign = previous_local.y > 0.0 ? 1.0 : -1.0;
  }
  else if (std::abs(previous_local.z) > half_outer.z)
  {
    axis = 2;
    sign = previous_local.z > 0.0 ? 1.0 : -1.0;
  }
  else
  {
    axis = 0;
    double min_distance = dx;
    if (dy < min_distance)
    {
      axis = 1;
      min_distance = dy;
    }
    if (dz < min_distance)
      axis = 2;

    if (axis == 0)
      sign = local.x >= 0.0 ? 1.0 : -1.0;
    else if (axis == 1)
      sign = local.y >= 0.0 ? 1.0 : -1.0;
    else
      sign = local.z >= 0.0 ? 1.0 : -1.0;
  }

  if (axis == 0)
  {
    const double old_velocity = particle.velocity.x;
    const double relative_velocity = particle.velocity.x - _box_velocity.x;
    particle.position.x = _box.center.x + sign * (half_outer.x + surface_offset);
    if (relative_velocity * sign < 0.0)
      particle.velocity.x = _box_velocity.x - relative_velocity * _config.collision_damping;
    _box_force.x += _config.particle_mass * (old_velocity - particle.velocity.x) / _config.dt;
  }
  else if (axis == 1)
  {
    const double old_velocity = particle.velocity.y;
    const double relative_velocity = particle.velocity.y - _box_velocity.y;
    particle.position.y = _box.center.y + sign * (half_outer.y + surface_offset);
    if (relative_velocity * sign < 0.0)
      particle.velocity.y = _box_velocity.y - relative_velocity * _config.collision_damping;
    _box_force.y += _config.particle_mass * (old_velocity - particle.velocity.y) / _config.dt;
  }
  else
  {
    const double old_velocity = particle.velocity.z;
    const double relative_velocity = particle.velocity.z - _box_velocity.z;
    particle.position.z = _box.center.z + sign * (half_outer.z + surface_offset);
    if (relative_velocity * sign < 0.0)
      particle.velocity.z = _box_velocity.z - relative_velocity * _config.collision_damping;
    _box_force.z += _config.particle_mass * (old_velocity - particle.velocity.z) / _config.dt;
  }
}

double SphSimulation::compute_box_mass() const
{
  const Vec3 inner_size
      = _box.outer_size - Vec3{2.0 * _box.wall_thickness, 2.0 * _box.wall_thickness, 2.0 * _box.wall_thickness};

  const double outer_volume = _box.outer_size.x * _box.outer_size.y * _box.outer_size.z;
  const double inner_volume = std::max(0.0, inner_size.x) * std::max(0.0, inner_size.y)
                              * std::max(0.0, inner_size.z);
  const double shell_volume = std::max(0.0, outer_volume - inner_volume);

  return std::max(1.0, _box.density * shell_volume);
}

double SphSimulation::poly6_kernel(double distance) const
{
  const double h = _config.smoothing_length;
  if (distance >= h)
    return 0.0;

  const double term = (h * h - distance * distance);
  return 315.0 / (64.0 * kPi * std::pow(h, 9)) * term * term * term;
}

double SphSimulation::spiky_gradient(double distance) const
{
  const double h = _config.smoothing_length;
  if (distance <= 1e-9 || distance >= h)
    return 0.0;

  const double term = (h - distance);
  return -45.0 / (kPi * std::pow(h, 6)) * term * term;
}

double SphSimulation::viscosity_laplacian(double distance) const
{
  const double h = _config.smoothing_length;
  if (distance >= h)
    return 0.0;

  return 45.0 / (kPi * std::pow(h, 6)) * (h - distance);
}
