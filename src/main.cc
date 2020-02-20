#include "grid.h"
#include "vtk_writer.h"

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
   grid.construct_esue(Grid::esue_neumann);
   grid.compute_cell_area();
   grid.compute_face_normal();
   grid.compute_cell_centroid();
   grid.compute_face_centroid();
}
