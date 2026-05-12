cmake_minimum_required(VERSION 3.16)

project(sph_building_demo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(sph_building_demo main.cpp)

if(MSVC)
  target_compile_options(sph_building_demo PRIVATE /W4)
else()
  target_compile_options(sph_building_demo PRIVATE -Wall -Wextra -Wpedantic)
endif()
