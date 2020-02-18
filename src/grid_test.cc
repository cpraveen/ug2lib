#include "grid.h"

// Function to dump all the constructed data from reference mesh
void grid_test1()
{
   //Read the reference mesh file ("test1.msh")
   Grid grid;
   grid.read_gmsh("test1.msh");

   std::cout << "Constructing dump data ... \n";
   grid.construct_esuf();
   grid.construct_esue(grid.esue_neumann);

   // Write the reference data to test1.dat
   std::cout << "Writing to test1.dat ... ";
   std::ofstream file("test1.dat");
   unsigned int n_bface = grid.get_n_bface()+1; // to sync cell numbering with gmsh numbering
   // Interior Faces
   file << "// Interior Faces\n";
   for(unsigned int iface=0; iface<grid.get_n_iface(); ++iface)
   {
      auto face_nodes = grid.get_iface_vertices(iface);
      auto face_nbr = grid.get_iface_cell(iface);
      file << "Face no.: " << iface+1 << "; Nodes: " << face_nodes[0]+1 << ", " << face_nodes[1]+1
           << "; left cell: " << face_nbr[0]+n_bface << ", right cell: " << face_nbr[1]+n_bface << std::endl;
   }

   // Boundary Faces
   file << "//Boundary Faces\n";
   for(unsigned int bface=0; bface<grid.get_n_bface(); ++bface)
   {
      auto face_nodes = grid.get_bface_vertices(bface);
      auto face_nbr = grid.get_bface_cell(bface);
      file << "Face no.: " << bface+1 << "; Nodes: " << face_nodes[0]+1 << ", " << face_nodes[1]+1
           << "; left cell: " << face_nbr+n_bface << std::endl;
   }

   // Elements surrounding point
   file << "//Elements surrounding point\n";
   for(unsigned int node=0; node<grid.get_n_vertex(); ++node)
   {
      auto esup = grid.get_esup(node);
      file << "Node: " << node+1 << "; Surrounding Element(s): ";
      for(unsigned int elem=0; elem<esup.first; ++elem)
         file << esup.second[elem]+n_bface << " ";
      file << std::endl;
   }

   // Points surrounding point (all_point = true)
   file << "//Points surrounding point\n";
   for(unsigned int node=0; node<grid.get_n_vertex(); ++node)
   {
      auto psup = grid.get_psup(node);
      file << "Node: " << node+1 << "; Surrounding Point(s): ";
      for(unsigned int jnode=0; jnode<psup.first; ++jnode)
         file << psup.second[jnode]+1 << " ";
      file << std::endl;
   }

   // Elements surrounding element (Von-Neumann)
   file << "//Elements surrounding element\n";
   for(unsigned int ielem=0; ielem<grid.get_n_cell(); ++ielem)
   {
      auto esue = grid.get_esue(ielem);
      file << "Cell no.: " << ielem+n_bface << "; Surrounding cell(s): ";
      for(unsigned int jelem=0; jelem<esue.first; ++jelem)
         file << esue.second[jelem]+n_bface << " ";
      file << std::endl;
   }

   file.close();
   std::cout << "Done\n";
}
