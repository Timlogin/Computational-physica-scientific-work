cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11")

project(vtkdemo)
add_executable(sample_cubic2d sample_cubic2d.cpp)
add_executable(sample_cubic3d sample_cubic3d.cpp)
add_executable(sample_tetr3d sample_tetr3d.cpp)

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/opt/homebrew/lib/libgmsh.dylib"
    INTERFACE_INCLUDE_DIRECTORIES "/opt/homebrew/include"
)
set(GMSH_LIBRARIES gmsh)

find_package(VTK REQUIRED)

include(${VTK_USE_FILE})

target_link_libraries(sample_cubic2d ${VTK_LIBRARIES})
target_link_libraries(sample_cubic3d ${VTK_LIBRARIES})

target_link_libraries(sample_tetr3d ${VTK_LIBRARIES})
target_link_libraries(sample_tetr3d ${GMSH_LIBRARIES})