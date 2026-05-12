#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Эта реализация специально сделана компактной:
// вода - SPH-частицы, здание - частицы с пружинами между соседями.
// Цель не в точном инженерном расчёте, а в наглядной particle-based модели.

struct Vec3
{
  double x = 0.0; // Координата или компонента вектора по оси X.
  double y = 0.0; // Координата или компонента вектора по оси Y.
  double z = 0.0; // Координата или компонента вектора по оси Z.
};

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

double dot(const Vec3& a, const Vec3& b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

double norm(const Vec3& v)
{
  return std::sqrt(dot(v, v));
}

Vec3 normalized(const Vec3& v)
{
  const double length = norm(v);
  if (length < 1e-12)
    return {};
  return v / length;
}

struct WaterParticle
{
  Vec3 position;      // Текущее положение частицы воды.
  Vec3 velocity;      // Текущая скорость частицы воды.
  Vec3 acceleration;  // Ускорение частицы на текущем шаге.
  double density = 1000.0;  // Локальная SPH-плотность около частицы.
  double pressure = 0.0;    // Давление, вычисленное по отклонению плотности от 1000 кг/м^3.
};

struct SolidParticle
{
  Vec3 position;          // Текущее положение частицы здания.
  Vec3 initial_position;  // Начальное положение, нужно для поля displacement в ParaView.
  Vec3 velocity;          // Скорость частицы здания.
  Vec3 force;             // Суммарная сила от пружин и удара воды.
  bool fixed = false;     // true для нижнего слоя: это условный фундамент здания.
};

struct Spring
{
  int a = 0;                 // Номер первой частицы здания.
  int b = 0;                 // Номер второй частицы здания.
  double rest_length = 0.0;  // Текущая "естественная" длина пружины.
  double initial_length = 0.0; // Исходная длина, нужна для ограничения остаточного смятия.
};

struct Config
{
  double dt = 0.009;       // Шаг по времени в секундах.
  int steps = 1000;        // Общее число шагов.
  int output_every = 2;    // Сохраняем каждый 2-й шаг: получится 500 кадров.

  double water_mass = 2000.0; // Полная масса воды, кг.
  double water_spacing = 0.1; // Расстояние между начальными частицами воды.
  double water_particle_mass = 1.0; // Масса одной частицы, пересчитывается после генерации воды.
  double smoothing_length = 0.45;   // Радиус влияния соседей в SPH.
  double rest_density = 1000.0;     // Опорная плотность воды, кг/м^3.
  double pressure_stiffness = 500.0; // Насколько быстро растёт давление при увеличении плотности.
  double viscosity = 0.5;          // Численная вязкость, сглаживает разницу скоростей.
  double air_drag = 0.04;           // Простое сопротивление воздуха для воды.

  double solid_spacing = 0.3;       // Расстояние между частицами здания.
  double solid_particle_mass = 12.0; // Масса одной частицы здания.
  double spring_stiffness = 8500.0;  // Жёсткость пружин между частицами здания.
  double spring_damping = 55.0;      // Демпфирование пружин, чтобы здание не дрожало бесконечно.
  double plastic_compression_ratio = 0.82; // Порог сжатия: ниже него пружина получает остаточную деформацию.
  double plastic_rate = 0.04;        // Скорость накопления остаточного смятия пружины.
  double solid_velocity_damping = 0.998; // Небольшое общее демпфирование скоростей здания.

  double contact_radius = 0.18;      // Радиус контакта воды с частицами здания.
  double contact_stiffness = 18000.0; // Жёсткость контактной силы вода-здание.
  double contact_damping = 160.0;    // Демпфирование контакта вода-здание.
  double ground_friction = 0.35;     // Потеря горизонтальной скорости воды на земле.

  Vec3 gravity{0.0, 0.0, -9.81};     // Гравитация действует на воду.
};

class Simulation
{
public:
  explicit Simulation(Config config) : _config(config)
  {
    create_water();
    create_building();
    compute_water_density_pressure();
  }

  void step()
  {
    compute_water_density_pressure();
    compute_water_forces();
    compute_spring_forces();
    compute_contact_forces();
    integrate();
    ++_step;
    _time += _config.dt;
  }

  int step_id() const { return _step; }
  double time() const { return _time; }
  const std::vector<WaterParticle>& water() const { return _water; }
  const std::vector<SolidParticle>& solid() const { return _solid; }
  const std::vector<Spring>& springs() const { return _springs; }

private:
  void create_water()
  {
    // Компактный блок воды объёмом 2 м^3 на высоте около 50 м.
    for (double x = -1.0; x <= 1.0 + 1e-9; x += _config.water_spacing)
    {
      for (double y = -0.5; y <= 0.5 + 1e-9; y += _config.water_spacing)
      {
        for (double z = 50.0; z <= 51.0 + 1e-9; z += _config.water_spacing)
        {
          WaterParticle p;
          p.position = {x, y, z};
          _water.push_back(p);
        }
      }
    }

    if (!_water.empty())
      _config.water_particle_mass = _config.water_mass / static_cast<double>(_water.size());
  }

  void create_building()
  {
    // Здание - не тонкая оболочка, а массивный блок частиц.
    // Так проще увидеть деформацию без сложной механики тонких стенок.
    const int nx = 11;
    const int ny = 7;
    const int nz = 10;
    const double h = _config.solid_spacing;

    auto index = [=](int i, int j, int k)
    {
      return i + nx * (j + ny * k);
    };

    _solid.resize(nx * ny * nz);
    for (int k = 0; k < nz; ++k)
    {
      for (int j = 0; j < ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          SolidParticle p;
          p.position = {
              (i - 0.5 * (nx - 1)) * h,
              (j - 0.5 * (ny - 1)) * h,
              k * h};
          p.initial_position = p.position;
          p.fixed = (k == 0);
          _solid[index(i, j, k)] = p;
        }
      }
    }

    // Связываем соседние частицы пружинами.
    // Кроме прямых соседей добавляем диагонали в ячейках, чтобы блок меньше
    // разваливался на "слои" и лучше сопротивлялся сдвигу.
    for (int k = 0; k < nz; ++k)
    {
      for (int j = 0; j < ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          add_spring_if_inside(index(i, j, k), i + 1, j, k, nx, ny, nz, index);
          add_spring_if_inside(index(i, j, k), i, j + 1, k, nx, ny, nz, index);
          add_spring_if_inside(index(i, j, k), i, j, k + 1, nx, ny, nz, index);
          add_spring_if_inside(index(i, j, k), i + 1, j + 1, k, nx, ny, nz, index);
          add_spring_if_inside(index(i, j, k), i + 1, j, k + 1, nx, ny, nz, index);
          add_spring_if_inside(index(i, j, k), i, j + 1, k + 1, nx, ny, nz, index);
        }
      }
    }
  }

  template <class IndexFunction>
  void add_spring_if_inside(int a, int i, int j, int k, int nx, int ny, int nz, IndexFunction index)
  {
    if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz)
      return;

    const int b = index(i, j, k);
    const double length = norm(_solid[b].position - _solid[a].position);
    _springs.push_back({a, b, length, length});
  }

  double poly6(double r) const
  {
    const double h = _config.smoothing_length;
    if (r >= h)
      return 0.0;
    const double q = h * h - r * r;
    return 315.0 / (64.0 * pi * std::pow(h, 9)) * q * q * q;
  }

  double spiky_grad(double r) const
  {
    const double h = _config.smoothing_length;
    if (r <= 1e-12 || r >= h)
      return 0.0;
    const double q = h - r;
    return -45.0 / (pi * std::pow(h, 6)) * q * q;
  }

  double viscosity_laplace(double r) const
  {
    const double h = _config.smoothing_length;
    if (r >= h)
      return 0.0;
    return 45.0 / (pi * std::pow(h, 6)) * (h - r);
  }

  void compute_water_density_pressure()
  {
    for (WaterParticle& pi : _water)
    {
      double density = _config.water_particle_mass * poly6(0.0);
      for (const WaterParticle& pj : _water)
      {
        const double r = norm(pi.position - pj.position);
        if (r > 1e-12)
          density += _config.water_particle_mass * poly6(r);
      }

      pi.density = std::max(density, 0.5 * _config.rest_density);
      pi.pressure = std::max(0.0, _config.pressure_stiffness * (pi.density - _config.rest_density));
    }
  }

  void compute_water_forces()
  {
    for (WaterParticle& pi : _water)
    {
      Vec3 pressure_force;
      Vec3 viscosity_force;

      for (const WaterParticle& pj : _water)
      {
        const Vec3 rij = pi.position - pj.position;
        const double r = norm(rij);
        if (r <= 1e-12 || r >= _config.smoothing_length)
          continue;

        const Vec3 dir = rij / r;
        const double p = 0.5 * (pi.pressure + pj.pressure);
        pressure_force += dir * (-_config.water_particle_mass * p * spiky_grad(r) / pj.density);

        const Vec3 dv = pj.velocity - pi.velocity;
        viscosity_force += dv * (_config.viscosity * _config.water_particle_mass * viscosity_laplace(r) / pj.density);
      }

      const Vec3 drag = pi.velocity * (-_config.air_drag);
      pi.acceleration = _config.gravity + (pressure_force + viscosity_force + drag) / pi.density;
    }
  }

  void compute_spring_forces()
  {
    // В этой демонстрации собственный вес здания не учитываем.
    // Иначе оно начинает деформироваться под гравитацией ещё до удара воды,
    // и становится непонятно, какая деформация вызвана именно падением воды.
    for (SolidParticle& p : _solid)
      p.force = {};

    for (Spring& spring : _springs)
    {
      SolidParticle& a = _solid[spring.a];
      SolidParticle& b = _solid[spring.b];

      const Vec3 dx = b.position - a.position;
      const double length = norm(dx);
      if (length < 1e-12)
        continue;

      const Vec3 dir = dx / length;
      const double stretch = length - spring.rest_length;
      const double relative_speed = dot(b.velocity - a.velocity, dir);
      const double force_value = _config.spring_stiffness * stretch + _config.spring_damping * relative_speed;
      const Vec3 force = dir * force_value;

      if (!a.fixed)
        a.force += force;
      if (!b.fixed)
        b.force -= force;

      // Простая остаточная деформация: если пружина сильно сжалась,
      // её "естественная" длина немного уменьшается.
      if (length < _config.plastic_compression_ratio * spring.rest_length)
      {
        const double new_rest = (1.0 - _config.plastic_rate) * spring.rest_length + _config.plastic_rate * length;
        spring.rest_length = std::max(0.55 * spring.initial_length, new_rest);
      }
    }
  }

  void compute_contact_forces()
  {
    for (WaterParticle& water : _water)
    {
      for (SolidParticle& solid : _solid)
      {
        const Vec3 delta = water.position - solid.position;
        const double distance = norm(delta);
        if (distance >= _config.contact_radius || distance < 1e-12)
          continue;

        const Vec3 normal = delta / distance;
        const double overlap = _config.contact_radius - distance;
        const double relative_speed = dot(water.velocity - solid.velocity, normal);
        const double force_value = _config.contact_stiffness * overlap - _config.contact_damping * relative_speed;
        if (force_value <= 0.0)
          continue;

        const Vec3 force = normal * force_value;
        water.acceleration += force / _config.water_particle_mass;
        if (!solid.fixed)
          solid.force -= force;
      }
    }
  }

  void integrate()
  {
    for (WaterParticle& p : _water)
    {
      p.velocity += p.acceleration * _config.dt;
      p.position += p.velocity * _config.dt;
      collide_water_with_ground(p);
    }

    for (SolidParticle& p : _solid)
    {
      if (p.fixed)
      {
        p.position = p.initial_position;
        p.velocity = {};
        continue;
      }

      const Vec3 acceleration = p.force / _config.solid_particle_mass;
      p.velocity += acceleration * _config.dt;
      p.velocity *= _config.solid_velocity_damping; // Небольшое численное демпфирование здания.
      p.position += p.velocity * _config.dt;
      collide_solid_with_ground(p);
    }
  }

  void collide_water_with_ground(WaterParticle& p) const
  {
    if (p.position.z >= 0.0)
      return;

    p.position.z = 0.0;
    if (p.velocity.z < 0.0)
      p.velocity.z = 0.0;
    p.velocity.x *= (1.0 - _config.ground_friction);
    p.velocity.y *= (1.0 - _config.ground_friction);
  }

  void collide_solid_with_ground(SolidParticle& p) const
  {
    if (p.position.z >= 0.0)
      return;

    p.position.z = 0.0;
    if (p.velocity.z < 0.0)
      p.velocity.z = 0.0;
    p.velocity.x *= 0.8;
    p.velocity.y *= 0.8;
  }

  static constexpr double pi = 3.14159265358979323846;

  Config _config;
  std::vector<WaterParticle> _water;
  std::vector<SolidParticle> _solid;
  std::vector<Spring> _springs;
  int _step = 0;
  double _time = 0.0;
};

std::string frame_path(const std::string& dir, const std::string& name, int frame)
{
  std::ostringstream out;
  out << dir << "/" << name << "_" << std::setw(4) << std::setfill('0') << frame << ".vtk";
  return out.str();
}

void write_water_vtk(const std::string& path, const std::vector<WaterParticle>& particles)
{
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Can not write " + path);

  out << "# vtk DataFile Version 3.0\n";
  out << "water particles\n";
  out << "ASCII\n";
  out << "DATASET POLYDATA\n";
  out << "POINTS " << particles.size() << " float\n";
  for (const WaterParticle& p : particles)
    out << p.position.x << " " << p.position.y << " " << p.position.z << "\n";

  out << "VERTICES " << particles.size() << " " << 2 * particles.size() << "\n";
  for (std::size_t i = 0; i < particles.size(); ++i)
    out << "1 " << i << "\n";

  out << "POINT_DATA " << particles.size() << "\n";
  out << "SCALARS density float 1\nLOOKUP_TABLE default\n";
  for (const WaterParticle& p : particles)
    out << p.density << "\n";

  out << "SCALARS pressure float 1\nLOOKUP_TABLE default\n";
  for (const WaterParticle& p : particles)
    out << p.pressure << "\n";

  out << "SCALARS speed float 1\nLOOKUP_TABLE default\n";
  for (const WaterParticle& p : particles)
    out << norm(p.velocity) << "\n";

  out << "VECTORS velocity float\n";
  for (const WaterParticle& p : particles)
    out << p.velocity.x << " " << p.velocity.y << " " << p.velocity.z << "\n";
}

void write_building_vtk(const std::string& path,
                        const std::vector<SolidParticle>& particles,
                        const std::vector<Spring>& springs)
{
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Can not write " + path);

  out << "# vtk DataFile Version 3.0\n";
  out << "particle building\n";
  out << "ASCII\n";
  out << "DATASET POLYDATA\n";
  out << "POINTS " << particles.size() << " float\n";
  for (const SolidParticle& p : particles)
    out << p.position.x << " " << p.position.y << " " << p.position.z << "\n";

  out << "VERTICES " << particles.size() << " " << 2 * particles.size() << "\n";
  for (std::size_t i = 0; i < particles.size(); ++i)
    out << "1 " << i << "\n";

  out << "LINES " << springs.size() << " " << 3 * springs.size() << "\n";
  for (const Spring& spring : springs)
    out << "2 " << spring.a << " " << spring.b << "\n";

  out << "POINT_DATA " << particles.size() << "\n";
  out << "SCALARS fixed float 1\nLOOKUP_TABLE default\n";
  for (const SolidParticle& p : particles)
    out << (p.fixed ? 1.0 : 0.0) << "\n";

  out << "SCALARS displacement float 1\nLOOKUP_TABLE default\n";
  for (const SolidParticle& p : particles)
    out << norm(p.position - p.initial_position) << "\n";

  out << "SCALARS speed float 1\nLOOKUP_TABLE default\n";
  for (const SolidParticle& p : particles)
    out << norm(p.velocity) << "\n";

  out << "VECTORS velocity float\n";
  for (const SolidParticle& p : particles)
    out << p.velocity.x << " " << p.velocity.y << " " << p.velocity.z << "\n";
}

int main()
{
  try
  {
    const std::string output_dir = "sph_building_demo/results/latest_run";
    std::filesystem::remove_all(output_dir);
    std::filesystem::create_directories(output_dir);

    Config config;
    Simulation simulation(config);

    int frame = 0;
    write_water_vtk(frame_path(output_dir, "water", frame), simulation.water());
    write_building_vtk(frame_path(output_dir, "building", frame), simulation.solid(), simulation.springs());

    for (int step = 0; step < config.steps; ++step)
    {
      simulation.step();

      if ((step + 1) % config.output_every == 0)
      {
        ++frame;
        write_water_vtk(frame_path(output_dir, "water", frame), simulation.water());
        write_building_vtk(frame_path(output_dir, "building", frame), simulation.solid(), simulation.springs());
        std::cout << "Saved frame " << frame << ", time = " << simulation.time() << " s\n";
      }
    }

    std::cout << "Done. Open numbered .vtk files from " << output_dir << " in ParaView.\n";
    std::cout << "Water particles: " << simulation.water().size() << "\n";
    std::cout << "Building particles: " << simulation.solid().size() << "\n";
    std::cout << "Springs: " << simulation.springs().size() << "\n";
  }
  catch (const std::exception& error)
  {
    std::cerr << "Error: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
