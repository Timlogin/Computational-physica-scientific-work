#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("t7");

  double lc = 1e-2;
  gmsh::model::geo::addPoint(0, 0, 0.2, lc, 1);
  gmsh::model::geo::addPoint(0, -0.1, 0.2, lc, 2);
  gmsh::model::geo::addPoint(0, 0.1, 0.2, lc, 3);
  gmsh::model::geo::addPoint(0, 0, 0, lc, 4);
  gmsh::model::geo::addPoint(0, -0.1, 0, lc, 5);
  gmsh::model::geo::addPoint(0, 0.1, 0, lc, 6);

  gmsh::model::geo::addLine(2, 5, 1);
  gmsh::model::geo::addLine(6, 3, 2);
  gmsh::model::geo::addCircleArc(2, 1, 3, 3);
  gmsh::model::geo::addCircleArc(5, 4, 6, 4);
  gmsh::model::geo::addCircleArc(3, 1, 2, 5);
  gmsh::model::geo::addCircleArc(6, 4, 5, 6);

  // я нашёл такую функцию как addSurfaceFilling
  // она натягивает поверхность на заданный контур
  // в отличие от addPlaneSurface, которая работает для плоских контуров

  gmsh::model::geo::addCurveLoop({1, 4, 2, -3}, 1);
  gmsh::model::geo::addSurfaceFilling({1}, 1);

  gmsh::model::geo::addCurveLoop({2, 5, 1, -6}, 2);
  gmsh::model::geo::addSurfaceFilling({2}, 2);

  gmsh::model::geo::addCurveLoop({3, 5}, 3);
  gmsh::model::geo::addPlaneSurface({3}, 3);

  gmsh::model::geo::addCurveLoop({4, 6}, 4);
  gmsh::model::geo::addPlaneSurface({4}, 4);

  gmsh::model::geo::addSurfaceLoop({1, 2, 3, 4}, 1);
  gmsh::model::geo::addVolume({1}, 1);

  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(3);

  gmsh::write("t7.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}
