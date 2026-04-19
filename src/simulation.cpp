#include "simulation.hpp"

#include <algorithm>
#include <cmath>

namespace
{
// Число pi
constexpr double kPi = 3.14159265358979323846;

// Удобная функция для ограничения значения заданным интервалом
double clamp(double value, double low, double high)
{
  return std::max(low, std::min(value, high));
}
}

// 3D-векторы
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

// При создании симуляции сразу вычисляем эффективную массу бокса
SphSimulation::SphSimulation(SimulationConfig config)
    : _config(config), _box(config.target), _box_reference_center(config.target.center)
{
  _box_mass = compute_box_mass();
}

void SphSimulation::initialize_fluid_block()
{
  _particles.clear();

  // Равномерно заполняем стартовый блок воды частицами
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

  // Делим суммарную массу воды поровну между всеми частицами
  if (!_particles.empty())
    _config.particle_mass = _config.total_water_mass / static_cast<double>(_particles.size());

  // Сразу считаем начальные плотности и давления
  compute_density_and_pressure();
}

void SphSimulation::step()
{
  // Импульс от воды по боксу пересчитывается заново каждый шаг
  _box_force = {};
  _roof_impact_load = 0.0;
  _roof_impact_x_moment = 0.0;
  _roof_impact_y_moment = 0.0;
  _side_x_impact_load = 0.0;
  _side_y_impact_load = 0.0;
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
    // Полная плотность складывается из вкладов соседних частиц.
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

    // Давление считаем по линейной связи с плотностью
    pi.pressure = _config.pressure_stiffness * (pi.density - _config.rest_density);
  }
}

void SphSimulation::compute_forces()
{
  for (Particle& pi : _particles)
  {
    // Силы давления и вязкости будем накапливать от всех соседей
    Vec3 pressure_force{};
    Vec3 viscosity_force{};

    for (const Particle& pj : _particles)
    {
      const Vec3 r = pi.position - pj.position;
      const double distance = norm(r);
      if (distance <= 1e-9 || distance >= _config.smoothing_length)
        continue;

      // Направление от соседней частицы к текущей
      const Vec3 direction = r / distance;

      // Сила давления стремится раздвигать частицы
      const double avg_pressure = (pi.pressure + pj.pressure) / 2.0;
      const double pressure_term
          = -_config.particle_mass * avg_pressure * spiky_gradient(distance) / pj.density;
      pressure_force += direction * pressure_term;

      // Вязкость сглаживает разницу скоростей
      const Vec3 velocity_diff = pj.velocity - pi.velocity;
      const double viscosity_term
          = _config.viscosity * _config.particle_mass * viscosity_laplacian(distance) / pj.density;
      viscosity_force += velocity_diff * viscosity_term;
    }

    // А почему бы не добавить сопротивление воздуха
    const Vec3 drag_force = pi.velocity * (-_config.air_drag);
    pi.acceleration = _config.gravity + (pressure_force + viscosity_force + drag_force) / pi.density;
  }
}

void SphSimulation::integrate_particles()
{
  for (Particle& particle : _particles)
  {
    const Vec3 previous_position = particle.position;

    // Выполняем шаг интегрирования
    particle.velocity += particle.acceleration * _config.dt;
    particle.position += particle.velocity * _config.dt;

    // После перемещения проверяем столкновения
    resolve_box_collision(particle, previous_position);
    resolve_world_collision(particle);
  }
}

void SphSimulation::integrate_box()
{
  // Бокс получает силу от удара воды и реакцию со стороны земли
  Vec3 total_force = _box_force;
  total_force += _box_velocity * (-_config.box_air_drag);

  const double half_height = 0.5 * _box.outer_size.z;
  const double bottom_z = _box.center.z - half_height;
  const double penetration = _config.ground.z - bottom_z;
  const bool in_ground_contact = penetration > 0.0;
  double ground_normal_force = 0.0;

  if (in_ground_contact)
  {
    // Вертикальная реакция земли против вдавливания бокса
    const double normal_force
        = std::max(0.0, _config.ground.normal_stiffness * penetration - _config.ground.normal_damping * _box_velocity.z);
    ground_normal_force = normal_force;
    total_force.z += normal_force;

    // Горизонтальная часть нужна, чтобы бокс не уезжал по плоскости
    const Vec3 lateral_shift = {_box_reference_center.x - _box.center.x, _box_reference_center.y - _box.center.y, 0.0};
    total_force.x += _config.ground.lateral_stiffness * lateral_shift.x
                     - _config.ground.lateral_damping * _box_velocity.x;
    total_force.y += _config.ground.lateral_stiffness * lateral_shift.y
                     - _config.ground.lateral_damping * _box_velocity.y;
  }

  // Обновляем движение бокса как одного жёсткого тела
  _box_acceleration = total_force / _box_mass;
  _box_velocity += _box_acceleration * _config.dt;
  _box.center += _box_velocity * _config.dt;

  // Ограничиваем максимальное вдавливание в землю
  const double min_bottom = _config.ground.z - _config.ground.max_sink;
  const double min_center_z = min_bottom + half_height;
  if (_box.center.z < min_center_z)
  {
    _box.center.z = min_center_z;
    if (_box_velocity.z < 0.0)
      _box_velocity.z = 0.0;
  }

  // Бокс не должен отскакивать от земли вверх как мячик
  if (_box.center.z > _box_reference_center.z)
  {
    _box.center.z = _box_reference_center.z;
    if (_box_velocity.z > 0.0)
      _box_velocity.z = 0.0;
  }

  const double max_shift_xy = 0.03;
  _box.center.x = clamp(_box.center.x, _box_reference_center.x - max_shift_xy, _box_reference_center.x + max_shift_xy);
  _box.center.y = clamp(_box.center.y, _box_reference_center.y - max_shift_xy, _box_reference_center.y + max_shift_xy);

  // После движения обновляем простую упругую форму бокса
  update_box_elasticity(ground_normal_force);
}

void SphSimulation::resolve_world_collision(Particle& particle)
{
  // Контакт воды с землёй делаем почти неупругим
  // частица не подпрыгивает, а оседает
  if (particle.position.z < _config.ground.z)
  {
    particle.position.z = _config.ground.z;
    if (particle.velocity.z < 0.0)
      particle.velocity.z *= (1.0 - _config.water_ground_absorption);

    const double slide_factor = std::max(0.0, 1.0 - _config.water_ground_friction);
    particle.velocity.x *= slide_factor;
    particle.velocity.y *= slide_factor;
  }

  // Универсальная проверка отражения от одной пары плоскостей
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
  apply_plane(particle.position.z, particle.velocity.z, _config.ground.z, _config.world_max.z);
}

void SphSimulation::resolve_box_collision(Particle& particle, const Vec3& previous_position)
{
  // Работаем в локальной системе координат бокса
  const Vec3 half_outer = _box.outer_size * 0.5;
  const double surface_offset = _config.box_collision_margin;
  const Vec3 local = particle.position - _box.center;
  const Vec3 previous_local = previous_position - _box.center;
  const bool inside_outer = std::abs(local.x) <= half_outer.x && std::abs(local.y) <= half_outer.y
                            && std::abs(local.z) <= half_outer.z;

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
    // Столкновение со стенкой, нормаль которой направлена вдоль оси X
    const double old_velocity = particle.velocity.x;
    const double relative_velocity = particle.velocity.x - _box_velocity.x;
    particle.position.x = _box.center.x + sign * (half_outer.x + surface_offset);
    if (relative_velocity * sign < 0.0)
      particle.velocity.x = _box_velocity.x - relative_velocity * _config.collision_damping;
    const double force_delta = _config.particle_mass * (old_velocity - particle.velocity.x) / _config.dt;
    _side_x_impact_load += std::abs(force_delta);
    _box_force.x += force_delta * _config.rigid_body_impact_fraction;
  }
  else if (axis == 1)
  {
    // Столкновение со стенкой, нормаль которой направлена вдоль оси Y
    const double old_velocity = particle.velocity.y;
    const double relative_velocity = particle.velocity.y - _box_velocity.y;
    particle.position.y = _box.center.y + sign * (half_outer.y + surface_offset);
    if (relative_velocity * sign < 0.0)
      particle.velocity.y = _box_velocity.y - relative_velocity * _config.collision_damping;
    const double force_delta = _config.particle_mass * (old_velocity - particle.velocity.y) / _config.dt;
    _side_y_impact_load += std::abs(force_delta);
    _box_force.y += force_delta * _config.rigid_body_impact_fraction;
  }
  else
  {
    // Столкновение со стенкой, нормаль которой направлена вдоль оси Z
    const double old_velocity = particle.velocity.z;
    const double relative_velocity = particle.velocity.z - _box_velocity.z;
    particle.position.z = _box.center.z + sign * (half_outer.z + surface_offset);
    if (relative_velocity * sign < 0.0)
      particle.velocity.z = _box_velocity.z - relative_velocity * _config.collision_damping;
    const double force_delta = _config.particle_mass * (old_velocity - particle.velocity.z) / _config.dt;

    if (sign > 0.0)
    {
      const double roof_load = std::max(0.0, -force_delta);
      _roof_impact_load += roof_load;
      _roof_impact_x_moment += roof_load * local.x;
      _roof_impact_y_moment += roof_load * local.y;
    }
    else
    {
      _box_force.z += force_delta * _config.rigid_body_impact_fraction;
    }
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

void SphSimulation::update_box_elasticity(double ground_normal_force)
{
  // Целевая деформация пропорциональна действующим силам,
  // но ограничивается сверху
  const double roof_target = clamp(
      _roof_impact_load / _config.roof_elastic_stiffness,
      0.0,
      _config.max_roof_deflection);
  const double floor_target = clamp(
      std::max(0.0, ground_normal_force) / _config.floor_elastic_stiffness,
      0.0,
      _config.max_floor_deflection);
  const double side_x_target = clamp(
      _side_x_impact_load / _config.wall_elastic_stiffness,
      0.0,
      _config.max_wall_deflection);
  const double side_y_target = clamp(
      _side_y_impact_load / _config.wall_elastic_stiffness,
      0.0,
      _config.max_wall_deflection);

  // Деформация не прыгает мгновенно, а плавно тянется к ней
  const double alpha = clamp(_config.elastic_relaxation_rate * _config.dt, 0.0, 1.0);
  _box_deformation.roof += alpha * (roof_target - _box_deformation.roof);
  _box_deformation.floor += alpha * (floor_target - _box_deformation.floor);
  _box_deformation.side_x += alpha * (side_x_target - _box_deformation.side_x);
  _box_deformation.side_y += alpha * (side_y_target - _box_deformation.side_y);

  const Vec3 half_outer = _box.outer_size * 0.5;
  const double center_alpha = clamp(_config.roof_center_relaxation_rate * _config.dt, 0.0, 1.0);
  double roof_center_x_target = 0.0;
  double roof_center_y_target = 0.0;

  if (_roof_impact_load > 1e-9)
  {
    roof_center_x_target = clamp(
        _roof_impact_x_moment / _roof_impact_load,
        -0.45 * half_outer.x,
        0.45 * half_outer.x);
    roof_center_y_target = clamp(
        _roof_impact_y_moment / _roof_impact_load,
        -0.45 * half_outer.y,
        0.45 * half_outer.y);
  }

  _box_deformation.roof_center_x += center_alpha * (roof_center_x_target - _box_deformation.roof_center_x);
  _box_deformation.roof_center_y += center_alpha * (roof_center_y_target - _box_deformation.roof_center_y);
}

double SphSimulation::poly6_kernel(double distance) const
{
  // Ядро плотности
  const double h = _config.smoothing_length;
  if (distance >= h)
    return 0.0;

  const double term = (h * h - distance * distance);
  return 315.0 / (64.0 * kPi * std::pow(h, 9)) * term * term * term;
}

double SphSimulation::spiky_gradient(double distance) const
{
  // Производная ядра давления
  const double h = _config.smoothing_length;
  if (distance <= 1e-9 || distance >= h)
    return 0.0;

  const double term = (h - distance);
  return -45.0 / (kPi * std::pow(h, 6)) * term * term;
}

double SphSimulation::viscosity_laplacian(double distance) const
{
  // Лапласиан ядра вязкости
  const double h = _config.smoothing_length;
  if (distance >= h)
    return 0.0;

  return 45.0 / (kPi * std::pow(h, 6)) * (h - distance);
}
