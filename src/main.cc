#include "grid.h"
#include "vtk_writer.h"
#include <fstream>

using namespace std;

int main()
{
   Grid grid;
   grid.read_gmsh("cylinder.msh");
   VTKWriter writer("cylinder.vtk", grid);
   writer.close();

   grid.construct_esup();
   grid.construct_psup(Grid::psup_edge);
   grid.construct_iface();
   grid.construct_esuf();
   grid.construct_esue(Grid::esue_moore);
   grid.compute_cell_area();
   grid.compute_face_normal();
   grid.compute_cell_centroid();
   grid.compute_face_centroid();
   grid.write_msh("test_mesh.msh");
}
