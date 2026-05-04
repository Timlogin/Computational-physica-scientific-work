cmake_minimum_required(VERSION 3.16)

project(water_drop_microproject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(VTK REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  IOXML
)

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
  src/simulation.cpp
)

target_include_directories(water_drop_demo PRIVATE src)

target_include_directories(water_drop_demo PRIVATE ${GMSH_INCLUDE_DIR})

target_link_libraries(water_drop_demo PRIVATE ${GMSH_LIBRARY})
target_link_libraries(water_drop_demo PRIVATE VTK::CommonCore VTK::CommonDataModel VTK::IOXML)

vtk_module_autoinit(
  TARGETS water_drop_demo
  MODULES VTK::CommonCore VTK::CommonDataModel VTK::IOXML
)

if(MSVC)
  target_compile_options(water_drop_demo PRIVATE /W4)
else()
  target_compile_options(water_drop_demo PRIVATE -Wall -Wextra -Wpedantic)
endif()
