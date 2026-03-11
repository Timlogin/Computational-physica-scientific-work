#include <set>
#include <gmsh.h>
#include <vector>
#include <utility>

// было принято стратегическое решение пользоваться ядром
// OpenCASCADE так как в geo задавать тор в моём понимании
// это 8 окружностей и много много точек, то есть код очень большой

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("t8");

  double x = 0;
  double y = 0;
  double z = 0;

  double R = 1;
  double r = 0.5;

  double h = 0.1;

  int first_torus = gmsh::model::occ::addTorus(x, y, z, R, r);
  int second_torus = gmsh::model::occ::addTorus(x, y, z, R, r - h);

  std::vector<std::pair<int, int>> big_torus = {{3, first_torus}};
  std::vector<std::pair<int, int>> small_torus = {{3, second_torus}};
  std::vector<std::pair<int, int>> result_torus;

  std::vector<std::vector<std::pair<int, int>>> results_figurs;
  // я долго не врубался зачем эта штука нужна, но без неё у меня не запускалось
  // results_figurs нужен если у нас например big_torus состоит из 5 фигур
  // и после вычитания у нас будет минимум 5 кусочков
  // так вот results_figurs будет все эти кусочки хранить по отдельности

  gmsh::model::occ::cut(big_torus, small_torus, result_torus, results_figurs);

  gmsh::model::occ::synchronize();

  double lc = h / 4;
  // тут я таким образом пытаюсь удовлетворить условие в 3 - 4 тетрайдера на толщину тора
  // неуверен, что это правильная статегия
  gmsh::option::setNumber("Mesh.MeshSizeMin", lc * 1);
  gmsh::option::setNumber("Mesh.MeshSizeMax", lc * 1.3);

  gmsh::model::mesh::generate(3);

  gmsh::write("t8.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}
