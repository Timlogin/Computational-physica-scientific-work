cmake_minimum_required(VERSION 3.16)
project(water_drop_microproject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(water_drop_demo
  src/main.cpp
  src/simulation.cpp
  src/vtk_writer.cpp
)

target_include_directories(water_drop_demo PRIVATE src)

if(MSVC)
  target_compile_options(water_drop_demo PRIVATE /W4)
else()
  target_compile_options(water_drop_demo PRIVATE -Wall -Wextra -Wpedantic)
endif()
