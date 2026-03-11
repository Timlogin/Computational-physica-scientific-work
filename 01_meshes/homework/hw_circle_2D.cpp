#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("t6");

  double lc = 1e-2;
  gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
  gmsh::model::geo::addPoint(0.1, 0, 0, lc, 2);
  gmsh::model::geo::addPoint(-0.1, 0, 0, lc, 3);

  gmsh::model::geo::addCircleArc(3, 1, 2, 1);
  gmsh::model::geo::addCircleArc(2, 1, 3, 2);

  gmsh::model::geo::addCurveLoop({1, 2}, 1);

  gmsh::model::geo::addPlaneSurface({1}, 1);

  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(2);

  gmsh::write("t6.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}
