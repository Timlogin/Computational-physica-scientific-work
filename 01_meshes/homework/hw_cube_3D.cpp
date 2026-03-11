#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("t5");

  double lc = 1e-2;
  gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
  gmsh::model::geo::addPoint(0.1, 0, 0, lc, 2);
  gmsh::model::geo::addPoint(0, 0.1, 0, lc, 3);
  gmsh::model::geo::addPoint(0.1, 0.1, 0, lc, 4);
  gmsh::model::geo::addPoint(0, 0, 0.1, lc, 5);
  gmsh::model::geo::addPoint(0.1, 0, 0.1, lc, 6);
  gmsh::model::geo::addPoint(0, 0.1, 0.1, lc, 7);
  gmsh::model::geo::addPoint(0.1, 0.1, 0.1, lc, 8);

  gmsh::model::geo::addLine(1, 2, 1);
  gmsh::model::geo::addLine(2, 4, 2);
  gmsh::model::geo::addLine(4, 3, 3);
  gmsh::model::geo::addLine(3, 1, 4);
  gmsh::model::geo::addLine(1, 5, 5);
  gmsh::model::geo::addLine(5, 6, 6);
  gmsh::model::geo::addLine(6, 2, 7);
  gmsh::model::geo::addLine(6, 8, 8);
  gmsh::model::geo::addLine(8, 7, 9);
  gmsh::model::geo::addLine(7, 5, 10);
  gmsh::model::geo::addLine(7, 3, 11);
  gmsh::model::geo::addLine(8, 4, 12);


  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);

  gmsh::model::geo::addCurveLoop({12, 3, -11, -9}, 2);
  gmsh::model::geo::addPlaneSurface({2}, 2);

  gmsh::model::geo::addCurveLoop({10, 6, 8, 9}, 3);
  gmsh::model::geo::addPlaneSurface({3}, 3);

  gmsh::model::geo::addCurveLoop({-10, 5, 4, 11}, 4);
  gmsh::model::geo::addPlaneSurface({4}, 4);

  gmsh::model::geo::addCurveLoop({7, 2, -12, -8}, 5);
  gmsh::model::geo::addPlaneSurface({5}, 5);

  gmsh::model::geo::addCurveLoop({5, 6, 7, -1}, 6);
  gmsh::model::geo::addPlaneSurface({6}, 6);


  gmsh::model::geo::addSurfaceLoop({1, 2, 3, 4, 5, 6}, 1);
  gmsh::model::geo::addVolume({1}, 1);

  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(3);

  gmsh::write("t5.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}
