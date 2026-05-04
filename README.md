## О проекте
Теория и физ модель описана в папке `real_physics`

Это микропроект по моделированию удара 2 тонн воды по кузову автомобиля с 50 метров.

Постановка звучит так:

- на автомобиль падает `2 тонны` воды;
- высота падения составляет `50 метров`;
- вместо реального минивэна на текущем этапе используется упрощённая геометрия;
- кузов моделируется как полый тонкостенный параллелепипед;
- основная цель текущей версии проекта — получить правдоподобную визуализацию процесса удара и смятия кузова в `ParaView`.


## Текущая постановка задачи

Рассматривается следующая упрощённая задача:

- вода - компактный объём `2 м^3`, что соответствует массе `2000 кг`;
- этот объём расположен на высоте `50 м`;
- кузов заменён на тонкостенный прямоугольный бокс;
- толщина стенок бокса на текущем этапе равна `0.8 мм` (число из интернета, то есть толщина кузова);
- материал кузова задаётся как высокопрочная сталь с плотностью `7800 кг/м^3`;
- бокс расположен на плоскости земли;
- земля не является полноценной моделью грунта, но работает как опора и не даёт кузову бесконечно уходить вниз;
- после удара часть деформации сохраняется, то есть бокс не возвращается полностью к исходной форме.

## Этап

### `C++`

Основная логика симуляции воды написана на `C++`.

На стороне C++ реализованы:

- структура проекта;
- запуск симуляции;
- модель воды на частицах;
- контакт воды с боксом;
- контакт воды с землёй;
- накопление эквивалентных нагрузок на кузов;
- экспорт результатов воды и служебных данных.

### `SPH`

Вода моделируется упрощённым вариантом метода частиц, близким к `SPH` (`Smoothed Particle Hydrodynamics`).

Сейчас для воды учитываются:

- гравитация;
- плотность;
- давление;
- вязкость;
- сопротивление воздуха;
- контакт с кузовом;
- контакт с землёй.

### `Gmsh`

`Gmsh` нужен для построения геометрии и сетки кузова.

С его помощью:

- строится внешний параллелепипед;
- из него вычитается внутренняя полость;
- получается тонкостенная объёмная геометрия;
- геометрия разбивается на тетраэдры.

Именно эта сетка потом используется в расчёте деформации.

### `FEniCSx`

`FEniCSx` используется для расчёта деформации кузова по тетраэдральной сетке.

На текущем этапе:

- вода сама не считается в `FEniCSx`;
- `FEniCSx` отвечает именно за деформацию бокса;
- из C++ в него передаются эквивалентные нагрузки от удара воды;
- затем решается задача линейной упругости;

### `VTK` и `ParaView`

Используются форматы:

- `.vtp` для воды и земли;
- `.vtu` для объёмной сетки кузова;
- `.pvd` для временных серий.

## Структура проекта

- [src](/Users/timofejloginov/Documents/Computational-physica-scientific-work/src) — исходный код проекта
- [build](/Users/timofejloginov/Documents/Computational-physica-scientific-work/build) — директория сборки
- [results](/Users/timofejloginov/Documents/Computational-physica-scientific-work/results) — выходные данные симуляции

Основные файлы:

- [src/main.cpp](/Users/timofejloginov/Documents/Computational-physica-scientific-work/src/main.cpp) — главный C++-файл: построение сетки `Gmsh`, запись `VTK`, запуск симуляции и режим просмотра сетки
- [src/simulation.hpp](/Users/timofejloginov/Documents/Computational-physica-scientific-work/src/simulation.hpp) — основные структуры данных и параметры модели
- [src/simulation.cpp](/Users/timofejloginov/Documents/Computational-physica-scientific-work/src/simulation.cpp) — реализация движения воды, столкновений и сбора нагрузок
- [src/fenics_box_solver.py](/Users/timofejloginov/Documents/Computational-physica-scientific-work/src/fenics_box_solver.py) — расчёт деформации кузова через `FEniCSx`
- [real_physics/physics_model.tex](/Users/timofejloginov/Documents/Computational-physica-scientific-work/real_physics/physics_model.tex) — описание физической модели с формулами

## Модели

### Вода

Вода представлена набором частиц.

Для каждой частицы считаются:

- положение;
- скорость;
- ускорение;
- плотность;
- давление.

Частицы взаимодействуют между собой и с внешними поверхностями.

### Кузов

Кузов представлен как:

- полый тонкостенный объём;
- тетраэдральная сетка;
- стальной материал;
- деформируемое тело с остаточным смятием.

### Земля

Земля моделируется как плоскость:

- не является отдельной сложной средой;
- ограничивает движение кузова вниз;
- не даёт воде улетать под модель;
- служит визуальной и механической опорой.

## Ограничения

Сейчас не учитываются:

- настоящая пластичность стали;
- разрушение кузова;
- трещины и отрыв элементов;
- нелинейная динамика реального грунта;
- точная геометрия минивэна;

Поэтому результаты чисто ознакомительные.

## Сборка проекта

Сборка выполняется стандартно через `CMake`.

```bash
cmake -S . -B build
cmake --build build -j2
```

После сборки появляется основной исполняемый файл:

- `build/water_drop_demo`

## Запуск проекта

### Основной расчёт

Полный запуск симуляции:

```bash
./build/water_drop_demo
```

### Просмотр только сетки бокса

Если нужно отдельно посмотреть сетку кузова без воды:

```bash
./build/water_drop_demo --mesh-only
```

Результаты сохраняются в:

- [results/mesh_preview](/Users/timofejloginov/Documents/Computational-physica-scientific-work/results/mesh_preview)

## Выходные файлы

После основного запуска используются:

- [results/latest_run/water_series.pvd](/Users/timofejloginov/Documents/Computational-physica-scientific-work/results/latest_run/water_series.pvd) — временная серия воды
- [results/latest_run/box_series.pvd](/Users/timofejloginov/Documents/Computational-physica-scientific-work/results/latest_run/box_series.pvd) — временная серия деформированного кузова
- [results/latest_run/ground/ground.vtp](/Users/timofejloginov/Documents/Computational-physica-scientific-work/results/latest_run/ground/ground.vtp) — плоскость земли
- [results/latest_run/mesh/box_shell.vtu](/Users/timofejloginov/Documents/Computational-physica-scientific-work/results/latest_run/mesh/box_shell.vtu) — исходная сетка кузова
- [results/latest_run/mesh/box_shell.msh](/Users/timofejloginov/Documents/Computational-physica-scientific-work/results/latest_run/mesh/box_shell.msh) — сетка в формате `gmsh`
- [results/latest_run/box_loads.csv](/Users/timofejloginov/Documents/Computational-physica-scientific-work/results/latest_run/box_loads.csv) — эквивалентные нагрузки от воды по кадрам

Папка `results/latest_run` пересоздаётся при каждом новом запуске.

## Улучшения в будующем (если надо)

- переход к более реалистичной модели пластичности металла;
- улучшение связи между SPH-водой и деформацией кузова;
- переход от Python-описания `FEniCSx`-части к `C++`-реализации;
- использование более реалистичной геометрии кузова;
- добавление численных графиков: прогиб, скорость, энергия, давление;
- более серьёзная модель контакта с землёй. (хотя неуверен, что это сильно как-то повлияет на результат)
