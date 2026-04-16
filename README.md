# Микропроект: Падение воды на минивэн

## Первая итерация

На текущем этапе:

- минивэн упрощен до полого прямоугольного параллелепипеда;
- вода моделируется базовым методом частиц, похожим на SPH;
- параллелепипед пока не является полноценной деформируемой оболочкой (FEM-оболочкой);
- получил первую визуализацию в ParaView.

## Текущая физическая модель

Текущая модель сильно упрощена.

- Вода представлена частицами с учётом гравитации, давления, вязкости и сопротивления воздуха.
- параллелепипед представлен как жёсткое полое тело с плотностью и эффективной массой.
- параллелепипед реагирует на удар как жёсткое тело на упругой опоре.
- Обработка столкновений настроена так, чтобы вода визуально не проходила сквозь ящик.

## Структура проекта

- [src](/Users/timofejloginov/Documents/Computational-physica-scientific-work/src) — исходный код
- [build](/Users/timofejloginov/Documents/Computational-physica-scientific-work/build) — директория сборки
- [results](/Users/timofejloginov/Documents/Computational-physica-scientific-work/results) — директория с результатами симуляции

Текущие исходные файлы:

- [src/main.cpp](/Users/timofejloginov/Documents/Computational-physica-scientific-work/src/main.cpp) — цикл запуска и сохранения кадров
- [src/simulation.hpp](/Users/timofejloginov/Documents/Computational-physica-scientific-work/src/simulation.hpp) — основные структуры и класс симуляции
- [src/simulation.cpp](/Users/timofejloginov/Documents/Computational-physica-scientific-work/src/simulation.cpp) — SPH-логика и реакция жёсткого параллелепипеда
- [src/vtk_writer.hpp](/Users/timofejloginov/Documents/Computational-physica-scientific-work/src/vtk_writer.hpp) — объявления функций для экспорта
- [src/vtk_writer.cpp](/Users/timofejloginov/Documents/Computational-physica-scientific-work/src/vtk_writer.cpp) — реализация экспорта в форматы VTK/VTP

## Сборка и запуск

```bash
cmake -S . -B build
cmake --build build -j2
./build/water_drop_demo
```

## Выходные данные

Файлы симуляции записываются в директорию:

- [results/latest_run](/Users/timofejloginov/Documents/Computational-physica-scientific-work/results/latest_run)

Эта папка пересоздаётся при каждом запуске.

Для лучшей видимости частиц воды в ParaView:

- можно использовать фильтры `Glyph` или `Point Gaussian`
- настроить масштабирование по массиву `radius`
