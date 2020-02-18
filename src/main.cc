#include "grid.h"

int main()
{
   Grid grid;
   grid.read_gmsh("cylinder.msh");
   grid.write_vtk("cylinder.vtk");
   grid.construct_esup();
   grid.construct_psup(true);
   grid.construct_iface();
   grid.construct_esuf();
   grid.construct_esue(Grid::esue_neumann);
   grid.compute_cell_area();
   grid.compute_face_normal();
   grid.compute_cell_centroid();
   grid.compute_face_centroid();
}
