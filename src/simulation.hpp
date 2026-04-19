#pragma once

#include <array>
#include <string>
#include <vector>

struct Vec3
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

// Базовые операции над векторами
Vec3 operator+(const Vec3& a, const Vec3& b);
Vec3 operator-(const Vec3& a, const Vec3& b);
Vec3 operator*(const Vec3& v, double s);
Vec3 operator/(const Vec3& v, double s);
Vec3& operator+=(Vec3& a, const Vec3& b);
Vec3& operator-=(Vec3& a, const Vec3& b);
Vec3& operator*=(Vec3& a, double s);

double dot(const Vec3& a, const Vec3& b);
double norm(const Vec3& v);

// Одна частица воды
struct Particle
{
  Vec3 position;      // Текущее положение частицы воды.
  Vec3 velocity;      // Текущая скорость частицы воды.
  Vec3 acceleration;  // Ускорение, посчитанное на текущем шаге.
  double density = 0.0;   // Локальная плотность, полученная через SPH-ядро.
  double pressure = 0.0;  // Давление частицы, связанное с плотностью.
};

// Параметры тонкостенного бокса.
struct BoxShell
{
  Vec3 center;            // Геометрический центр бокса.
  Vec3 outer_size;        // Внешние размеры бокса: длина, ширина и высота.
  double wall_thickness = 0.02; // Толщина стенки.
  double density = 7800.0;      // Плотность материала стенок.
};

// Параметры плоскости земли.
struct GroundPlane
{
  double z = 0.0;                     // Высота плоскости земли.
  double half_size_x = 7.0;          // Полуразмер земли по X для визуализации.
  double half_size_y = 5.0;          // Полуразмер земли по Y для визуализации.
  double normal_stiffness = 280000.0; // Жёсткость реакции земли по нормали.
  double normal_damping = 46000.0;    // Демпфирование по нормали.
  double lateral_stiffness = 120000.0; // Сопротивление сдвигу по плоскости.
  double lateral_damping = 52000.0;    // Демпфирование бокового сдвига.
  double max_sink = 0.03;              // Максимально допустимое продавливание в землю.
  std::array<int, 3> color_rgb{82, 122, 74}; // Цвет земли в ParaView.
};

// Простая модель упругой деформации бокса
struct BoxElasticDeformation
{
  double roof = 0.0;           // Прогиб крыши вниз.
  double floor = 0.0;          // Прогиб днища вверх.
  double side_x = 0.0;         // Смятие боковых стенок, нормальных к оси X.
  double side_y = 0.0;         // Смятие боковых стенок, нормальных к оси Y.
  double roof_center_x = 0.0;  // Центр основного удара по крыше по X.
  double roof_center_y = 0.0;  // Центр основного удара по крыше по Y.
};

// Все основные параметры симуляции
struct SimulationConfig
{
  double dt = 0.004;     // Шаг по времени.
  int total_steps = 2000; // Общее число шагов расчёта.
  int output_every = 20;  // Как часто сохранять кадр в результаты.

  double particle_spacing = 0.15;    // Расстояние между стартовыми частицами воды.
  double smoothing_length = 0.28;    // Радиус влияния частиц в SPH-ядре.
  double total_water_mass = 2000.0;  // Полная масса воды.
  double particle_mass = 1.0;        // Масса одной частицы, позже пересчитывается.
  double rest_density = 1000.0;      // Опорная плотность воды.
  double pressure_stiffness = 120.0; // Коэффициент связи давления и плотности.
  double viscosity = 24.0;           // Коэффициент вязкости.
  double air_drag = 0.05;            // Сопротивление воздуха для частиц воды.
  double visual_particle_radius = 0.22; // Радиус частицы только для визуализации.
  double water_ground_friction = 0.18;  // Трение воды о землю.
  double water_ground_absorption = 1.0; // Потеря вертикальной скорости на земле.

  Vec3 gravity{0.0, 0.0, -9.81}; // Вектор ускорения свободного падения.

  Vec3 fluid_block_min{-1.0, -0.5, 50.0}; // Нижний угол стартового объёма воды.
  Vec3 fluid_block_max{1.0, 0.5, 51.0};   // Верхний угол стартового объёма воды.

  BoxShell target{{0.0, 0.0, 0.9}, {4.0, 2.0, 1.8}, 0.0008, 7800.0}; // Параметры кузова.
  double box_collision_margin = 0.05;  // Небольшой запас при выталкивании воды из стенки.
  double box_mesh_size = 0.16;         // Типичный размер тетраэдра в сетке бокса.
  double box_air_drag = 30.0;          // Сопротивление воздуха для движения бокса как тела.
  double rigid_body_impact_fraction = 0.10; // Какая доля удара идёт в перенос бокса как целого тела.
  double roof_elastic_stiffness = 180000.0; // Жёсткость крыши в упрощённой модели.
  double floor_elastic_stiffness = 950000.0; // Жёсткость днища в упрощённой модели.
  double wall_elastic_stiffness = 420000.0;  // Жёсткость боковых стенок в упрощённой модели.
  double elastic_relaxation_rate = 5.0;      // Скорость подстройки упрощённой деформации.
  double roof_center_relaxation_rate = 6.0;  // Скорость смещения центра удара по крыше.
  double max_roof_deflection = 0.30;         // Ограничение на прогиб крыши.
  double max_floor_deflection = 0.05;        // Ограничение на прогиб днища.
  double max_wall_deflection = 0.14;         // Ограничение на прогиб стенок.
  GroundPlane ground{};                      // Параметры земли.

  Vec3 world_min{-6.0, -4.0, 0.0}; // Нижняя граница нашего мира
  Vec3 world_max{6.0, 4.0, 60.0};  // Верхняя граница нашего мира.
  double collision_damping = 0.35; // Коэффициент потери скорости при столкновении.
};

// Основной класс симуляции воды и реакции бокса
class SphSimulation
{
public:
  // Создаёт симуляцию с заданной конфигурацией
  explicit SphSimulation(SimulationConfig config);

  // Формирует стартовый объём воды в виде набора частиц
  void initialize_fluid_block();

  // Выполняет один шаг по времени
  void step();

  const SimulationConfig& config() const { return _config; }
  const std::vector<Particle>& particles() const { return _particles; }
  const BoxShell& box() const { return _box; }
  const GroundPlane& ground() const { return _config.ground; }
  const Vec3& box_velocity() const { return _box_velocity; }
  const BoxElasticDeformation& box_deformation() const { return _box_deformation; }
  double roof_impact_load() const { return _roof_impact_load; }
  double roof_impact_center_x() const
  {
    return _roof_impact_load > 1e-9 ? _roof_impact_x_moment / _roof_impact_load : 0.0;
  }
  double roof_impact_center_y() const
  {
    return _roof_impact_load > 1e-9 ? _roof_impact_y_moment / _roof_impact_load : 0.0;
  }
  double side_x_impact_load() const { return _side_x_impact_load; }
  double side_y_impact_load() const { return _side_y_impact_load; }
  double box_mass() const { return _box_mass; }
  int current_step() const { return _current_step; }
  double current_time() const { return _current_time; }

private:
  // Основные этапы одного шага SPH-модели
  void compute_density_and_pressure();
  void compute_forces();
  void integrate_particles();
  void integrate_box();
  void update_box_elasticity(double ground_normal_force);

  // Обработка столкновений частиц с внешними границами и с боксом
  void resolve_world_collision(Particle& particle);
  void resolve_box_collision(Particle& particle, const Vec3& previous_position);

  // Эффективная масса бокса по его геометрии и плотности материала
  double compute_box_mass() const;

  // SPH-ядра для плотности, давления и вязкости
  double poly6_kernel(double distance) const;
  double spiky_gradient(double distance) const;
  double viscosity_laplacian(double distance) const;

  // Текущее состояние симуляции.
  SimulationConfig _config;
  std::vector<Particle> _particles; // Все частицы воды.
  BoxShell _box;                    // Текущее состояние бокса.
  Vec3 _box_velocity{};             // Скорость бокса как жёсткого тела.
  Vec3 _box_acceleration{};         // Ускорение бокса как жёсткого тела.
  Vec3 _box_force{};                // Суммарная сила от ударов воды на текущем шаге.
  Vec3 _box_reference_center{};     // Исходное положение центра бокса.
  double _roof_impact_load = 0.0;   // Нагрузка на крышу на текущем шаге.
  double _roof_impact_x_moment = 0.0; // Момент нагрузки по X для поиска центра удара.
  double _roof_impact_y_moment = 0.0; // Момент нагрузки по Y для поиска центра удара.
  double _side_x_impact_load = 0.0; // Боковая нагрузка на стенки, нормальные к X.
  double _side_y_impact_load = 0.0; // Боковая нагрузка на стенки, нормальные к Y.
  BoxElasticDeformation _box_deformation{}; // Упрощённая C++-деформация для отладки и экспорта.
  double _box_mass = 1.0;           // Эффективная масса стенок бокса.
  int _current_step = 0;            // Номер текущего шага по времени.
  double _current_time = 0.0;       // Текущее физическое время расчёта.
};
