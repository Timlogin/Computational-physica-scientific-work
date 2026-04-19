cmake_minimum_required(VERSION 3.16)

project(water_drop_microproject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_path(GMSH_INCLUDE_DIR gmsh.h
  PATHS
    /opt/homebrew/include
    /usr/local/include
    /opt/local/include
)

find_library(GMSH_LIBRARY gmsh
  PATHS
    /opt/homebrew/lib
    /usr/local/lib
    /opt/local/lib
)

add_executable(water_drop_demo
  src/main.cpp
  src/gmsh_box.cpp
  src/simulation.cpp
  src/vtk_writer.cpp
)

add_executable(box_mesh_preview
  src/box_mesh_preview.cpp
  src/gmsh_box.cpp
  src/simulation.cpp
  src/vtk_writer.cpp
)

target_include_directories(water_drop_demo PRIVATE src)
target_include_directories(box_mesh_preview PRIVATE src)

target_include_directories(water_drop_demo PRIVATE ${GMSH_INCLUDE_DIR})
target_include_directories(box_mesh_preview PRIVATE ${GMSH_INCLUDE_DIR})

target_link_libraries(water_drop_demo PRIVATE ${GMSH_LIBRARY})
target_link_libraries(box_mesh_preview PRIVATE ${GMSH_LIBRARY})

if(MSVC)
  target_compile_options(water_drop_demo PRIVATE /W4)
  target_compile_options(box_mesh_preview PRIVATE /W4)
else()
  target_compile_options(water_drop_demo PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(box_mesh_preview PRIVATE -Wall -Wextra -Wpedantic)
endif()
