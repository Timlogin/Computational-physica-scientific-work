#include "gmsh_box.hpp"
#include "simulation.hpp"
#include "vtk_writer.hpp"

#include <filesystem>
#include <iostream>

int main()
{
  try
  {
    // Этот режим нужен только для просмотра самой сетки бокса
    // без запуска полной симуляции воды.
    const std::string output_dir = "results/mesh_preview";
    std::filesystem::remove_all(output_dir);
    std::filesystem::create_directories(output_dir);

    // Берём ту же конфигурацию, что и у основного проекта.
    SimulationConfig config;

    // Сохраняем родной файл Gmsh, чтобы было удобно смотреть именно сетку.
    const TetMesh mesh
        = build_box_shell_mesh(config.target, config.box_mesh_size, output_dir + "/box_shell.msh");

    // И дополнительно сохраняем эту же сетку в VTU для ParaView.
    write_box_vtu(output_dir + "/box_shell.vtu", mesh, {0.0, 0.0, 0.0}, config.target, {0.0, 0.0, 0.0}, 0.0, {});
  }
  catch (const std::exception& error)
  {
    std::cerr << "Mesh preview error: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
