#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <utility> // std::pair, std::make_pair
#include <algorithm>
#include <cmath>

#include "grid.h"

using namespace std;

// constructor
Grid::Grid()
{
   has_esup = false;
   has_psup_all = false;
   has_psup_edge = false;
   has_iface = false;
   has_esuf = false;
   has_esue_moore = false;
   has_esue_neumann = false;
}

// destructor
Grid::~Grid()
{
   if(!coord) delete [] coord;
   if(!cell1) delete [] cell1;
   if(!cell2) delete [] cell2;
   if(!bface) delete [] bface;
   if(!bface_type) delete [] bface_type;
   if(!bface_cell) delete [] bface_cell;
   if(!esup1) delete [] esup1;
   if(!esup2) delete [] esup2;
   if(!psup1) delete [] psup1;
   if(!psup2) delete [] psup2;
   if(!iface) delete [] iface;
   if(!iface_cell)  delete [] iface_cell;
   if(!esue1) delete [] esue1;
   if(!esue2) delete [] esue2;
}

// Read mesh from file in gmsh format.
// Currently it can only read msh2 format.
void Grid::read_gmsh(const string grid_file)
{
   // Temporary arrays
   unsigned int n_elem;
   vector<unsigned int> elem1, elem2, elem_type, elem_phy;

   cout << "Reading gmsh grid file " << grid_file << endl;

   ifstream file;
   file.open(grid_file.c_str());
   assert(file.is_open());

   string line;
   int tag;
   unsigned int count, ntags;
   double rdummy, version;

   // Skip some lines
   file >> line;
   file >> version >> tag >> tag;
   assert(version == 2.2);
   file >> line;
   file >> line;

   // Read vertices
   file >> n_vertex;
   assert(n_vertex > 0);
   cout << "Number of vertices = " << n_vertex << endl;
   coord = new double[dim*n_vertex];

   for(unsigned int i=0; i<n_vertex; ++i)
      file >> count
           >> coord[i*dim]
           >> coord[i*dim+1]
           >> rdummy; // 2d, dont need z coordinate

   // Skip two lines
   file >> line;
   file >> line;

   // Read cells
   file >> n_elem;
   assert(n_elem > 0);
   cout << "Numbers of elements = " << n_elem << endl;

   // for triangle and quad only
   elem1.resize(n_elem+1);
   elem_type.resize(n_elem);
   elem_phy.resize(n_elem);
   elem1[0] = 0;

   n_tri = 0;
   n_quad = 0;
   n_bface = 0;
   for(unsigned int i=0; i<n_elem; ++i)
   {
      file >> count
           >> elem_type[i]
           >> ntags;

      // Some gmsh files have 2 and some have 3 tags
      assert(ntags==2 || ntags==3);

      if(elem_type[i] == Gmsh::Line2) // Line face
      {
         file >> elem_phy[i]; // First tag is physical type
         // Dummy tags
         for(unsigned int p=1; p<ntags; ++p)
            file >> tag;

         elem1[i+1] = elem1[i] + 2;
         elem2.resize(elem1[i+1]);
         file >> elem2[elem1[i]]
              >> elem2[elem1[i]+1];
         ++n_bface;
      }
      else if(elem_type[i] == Gmsh::Triangle3) // Triangular cell
      {
         file >> elem_phy[i]; // First tag is physical type
         for(unsigned int p = 1; p<ntags; ++p)
            file >> tag ; // Dummy tags

         elem1[i+1] = elem1[i] + 3;
         elem2.resize(elem1[i+1]);
         file >> elem2[elem1[i]]
              >> elem2[elem1[i]+1]
              >> elem2[elem1[i]+2];
         ++n_tri;
      }
      else if(elem_type[i] == Gmsh::Quadrilateral4) // Quadrilateral cell
      {
         file >> elem_phy[i]; // First tag is physical type
         for(unsigned int p = 1; p<ntags; ++p)
            file >> tag ; // Dummy tags

         elem1[i+1] = elem1[i] + 4;
         elem2.resize(elem1[i+1]);
         file >> elem2[elem1[i]]
              >> elem2[elem1[i]+1]
              >> elem2[elem1[i]+2]
              >> elem2[elem1[i]+3];
         ++n_quad;
      }
      else
      {
         cout << "Unknown element type !!!" << endl;
         cout << "   Element type =" << elem_type[i] << endl;
         exit(0);
      }
   }
   file.close();

   n_cell = n_tri + n_quad;

   cout << "Numbers of boundary faces = " << n_bface << endl;
   cout << "Number of triangles       = " << n_tri << endl;
   cout << "Number of quadrilaterals  = " << n_quad << endl;
   cout << "Total number of cells     = " << n_cell << endl;

   assert(n_cell > 0);
   assert(n_bface > 0);
   assert(elem1[n_elem] == 2*n_bface + 3*n_tri + 4*n_quad);

   // gmsh vertex numbering starts at 1, decrease by 1
   for(unsigned int i=0; i<elem1[n_elem]; ++i)
      --elem2[i];

   cell1 = new unsigned int[3*n_tri + 4*n_quad];
   cell2 = new unsigned int[n_cell+1];
   bface = new unsigned int[2*n_bface];
   bface_type = new int[n_bface];
   ctype = new int[n_cell];

   cell2[0] = 0;
   unsigned int bface_count = 0;
   unsigned int cell_count = 0;

   for(unsigned int i=0; i<n_elem; ++i)
   {
      if(elem_type[i] == 1) // Line face
      {
         bface_type[bface_count] = elem_phy[i];
         bface[2*bface_count] = elem2[elem1[i]];
         bface[2*bface_count+1] = elem2[elem1[i]+1];
         ++bface_count;
      }
      else if(elem_type[i] == 2 || elem_type[i] == 3) // tri and quad
      {
         cell2[cell_count+1] = cell2[cell_count] + elem1[i+1] - elem1[i];
         for(unsigned int j=0; j<elem1[i+1]-elem1[i]; ++j)
            cell1[cell2[cell_count]+j] = elem2[elem1[i]+j];
         ctype[cell_count] = elem_type[i];
         ++cell_count;
      }
      else
      {
         cout << "Unknown element type = " << elem_type[i] << endl;
         exit(0);
      }
   }

   assert(bface_count == n_bface);
   assert(cell_count == n_cell);

   // Delete local arrays
   elem1.resize(0);
   elem2.resize(0);
   elem_type.resize(0);
   elem_phy.resize(0);
}

// Write grid in vtk format
// ofstream must be opened before calling this function.
void Grid::write_vtk(std::ofstream& vtk)
{
   assert(vtk.is_open());

   // Save vtk file
   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << "Test file" << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET UNSTRUCTURED_GRID" << endl;
   vtk << "POINTS  " << n_vertex << "  float" << endl;

   for(unsigned int i=0; i<n_vertex; ++i)
   {
      auto x = get_coord(i);
      vtk << x[0] << " " << x[1] << " " << 0.0 << endl;
   }

   vtk << "CELLS  " << n_cell << " " << 4 * n_tri + 5 * n_quad << endl;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      auto cell = get_cell_vertices(i);
      vtk << cell.first;
      for(unsigned int j=0; j<cell.first; ++j)
         vtk << "  " << cell.second[j];
      vtk << endl;
   }

   vtk << "CELL_TYPES  " << n_cell << endl;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      if(ctype[i] == 2) vtk << 5 << endl;
      else if(ctype[i] == 3) vtk << 9 << endl;
      else
      {
         cout << "Unknown ctype = " << ctype[i] << endl;
         exit(0);
      }
   }

   vtk.close();
}

// Write grid into file specifed by filename argument
void Grid::write_vtk(const string filename)
{
   std::ofstream vtk;
   vtk.open(filename);
   write_vtk(vtk);
   vtk.close();
   cout << "Wrote " << filename << endl;
}

// Find cells surrounding a point
// See Lohner: section xxx
void Grid::construct_esup()
{
   if(has_esup) return;
   cout << "Constructing elements surrounding point ... ";

   esup2 = new unsigned int[n_vertex+1];
   for(unsigned int i=0; i<n_vertex+1; ++i) // Initialize esup2
   {
      esup2[i] = 0;
   }

   for(unsigned int icell=0; icell<n_cell; ++icell) // Loop over all the cells
   {
      auto cell = get_cell_vertices(icell);
      for(unsigned int ipoint=0; ipoint<cell.first; ++ipoint)  //Loop over nodes of the cell
      {
         ++esup2[cell.second[ipoint]+1]; // Count the cells surrounding the node
      }
   }

   for(unsigned int i=0; i<n_vertex; ++i) // Add the previous element to create the esup2
   {
      esup2[i+1] += esup2[i];
   }

   esup1 = new unsigned int[esup2[n_vertex]]; // Create and initialize esup1
   for(unsigned int i=0; i<esup2[n_vertex]; ++i)
   {
      esup1[i] = 0; // not really required
   }

   for(unsigned int icell=0; icell<n_cell; ++icell)
   {
      auto cell = get_cell_vertices(icell);
      for(unsigned int ipoint=0; ipoint<cell.first; ++ipoint)
      {
         unsigned int gnode = cell.second[ipoint]; // Get the global node number
         unsigned int istor = esup2[gnode]; // location in esup1 where the element will be stored
         esup1[istor] = icell;
         ++esup2[gnode];
      }
   }

   for(unsigned int i=n_vertex; i>0; --i)
   {
      esup2[i] = esup2[i-1]; // reshuffle the esup2
   }
   esup2[0] = 0;
   has_esup = true;
   cout << "Done\n";
}

// Find points surrounding a point
// If all_points==false, find only points connected by an edge.
void Grid::construct_psup(const psup_type type)
{
   if(type == psup_all && has_psup_all) return;
   if(type == psup_edge && has_psup_edge) return;
   construct_esup(); // psup needs esup data

   // Different type of psup requested
   // Delete existing memory
   if(!psup1) delete[] psup1;
   if(!psup2) delete[] psup2;

   if(type == psup_all)
      cout << "Constructing points surrounding point: all ... ";
   else
      cout << "Constructing points surrounding point: edge ... ";

   psup2 = new unsigned int[n_vertex+1];
   unsigned int* lpoin = new unsigned int[n_vertex]; // Help array to avoid duplication from neighbouring cells

   // Initialize
   for(unsigned int i=0; i<n_vertex; ++i)
   {
      psup2[i] = 0;
      lpoin[i] = 0;
   }
   psup2[n_vertex] = 0;

   std::vector<unsigned int> psup1_temp(0);  // temporary vector to allow resize function
   unsigned int istor = 0;
   for(unsigned int ipoint=0; ipoint<n_vertex; ++ipoint) // Loop over all the nodes
   {
      auto esup = get_esup(ipoint); // get cells surrounding the node
      for(unsigned int icell=0; icell<esup.first; ++icell) // Loop over cells surrounding node
      {
         auto cell = get_cell_vertices(esup.second[icell]);
         for(unsigned int jpoint=0; jpoint<cell.first; ++jpoint) // Loop over nodes of cell
         {
            unsigned int gnode = cell.second[jpoint];  // global node number
            if(gnode != ipoint && lpoin[gnode] != ipoint+1) // check for duplication
            {
               if(type == psup_all) // get all the points surrounding the node
               {
                  psup1_temp.resize(psup1_temp.size()+1);
                  psup1_temp[istor] = gnode;
                  lpoin[gnode] = ipoint+1; // using ipoint+1 as 0 is the initialized value
                  ++istor;
               }
               else
               {
                  int prev_node = jpoint-1, next_node = (jpoint+1) % cell.first;
                  if(prev_node < 0) prev_node = cell.first-1;
                  // Check if connected by edge
                  if(cell.second[prev_node] == ipoint || cell.second[next_node] == ipoint)
                  {
                     psup1_temp.resize(psup1_temp.size()+1);
                     psup1_temp[istor] = gnode;
                     lpoin[gnode] = ipoint+1;
                     ++istor;
                  }
               }
            }
         }
      }
      psup2[ipoint+1] = istor;
   }
   // Copy data from psup1_temp to psup1
   psup1 = new unsigned int[psup2[n_vertex]];
   for(unsigned int i=0; i<psup1_temp.size(); ++i)
      psup1[i] = psup1_temp[i];

   psup1_temp.resize(0);
   delete [] lpoin;
   has_psup_all = (type == psup_all);
   has_psup_edge = (type == psup_edge);
   cout << "Done\n";
}

// Find all the interior faces with corresponding node numbers
void Grid::construct_iface()
{
   if(has_iface) return;
   construct_psup(psup_edge); // construct psup data
   cout << "Constructing faces ... ";
   n_iface = 0;
   std::vector<unsigned int> face_temp(0); // to enable resize function

   for(unsigned int ipoint=0; ipoint<n_vertex; ++ipoint) // Loop over all the nodes
   {
      auto psup = get_psup(ipoint);    // get the nodes connected to ipoint
      for(unsigned int jpoint=0; jpoint<psup.first; ++jpoint)  // loop over the surrounding points
      {
         // check if the face is not already defined with a previous ipoint
         if(ipoint<psup.second[jpoint])
         {
            // Get esup for both nodes and take intersection
            // 2 cells => directly connected and interior face,
            // 1 cell => either not connected by an edge or a boundary face
            auto esup_node1 = get_esup(ipoint);  // get elements surrounding node1
            auto esup_node2 = get_esup(psup.second[jpoint]);  // get elements surrounding node2
            std::vector<unsigned int> cell_node1(esup_node1.first);
            std::vector<unsigned int> cell_node2(esup_node2.first);
            for(unsigned int i=0; i<esup_node1.first; ++i)
               cell_node1[i] = esup_node1.second[i];
            for(unsigned int i=0; i<esup_node2.first; ++i)
               cell_node2[i] = esup_node2.second[i];

            // Sort the arrays, needed for intersection
            std::sort(cell_node1.begin(), cell_node1.end());
            std::sort(cell_node2.begin(), cell_node2.end());

            // find the intersection
            std::vector<unsigned int> inter(std::min(esup_node1.first, esup_node2.first));
            auto ls = std::set_intersection(cell_node1.begin(), cell_node1.end(),
                                            cell_node2.begin(), cell_node2.end(), inter.begin());
            unsigned int n_esuf = ls - inter.begin(); // number of faces surrounding the element

            if(n_esuf == 2)   // interior face
            {
               unsigned int istor = face_temp.size();
               face_temp.resize(face_temp.size()+2);
               face_temp[istor] = ipoint;
               face_temp[istor+1] = psup.second[jpoint];
               ++n_iface;
            }
         }
      }
   }
   iface = new unsigned int[face_temp.size()];
   for(unsigned int i=0; i<face_temp.size(); ++i)  // copy the data from temporary vector
      iface[i] = face_temp[i];

   face_temp.resize(0);
   has_iface = true;
   cout << "Done\n";
   cout << "Total number of interior faces = " << n_iface << endl;
}

// Find neighbouring cells for a face
void Grid::construct_esuf()
{
   if(has_esuf) return;
   construct_iface();
   cout << "Constructing elements surrounding face ... ";

   //Interior Faces
   iface_cell = new unsigned int[2*n_iface];
   for(unsigned int jface=0; jface<n_iface; ++jface) // loop over all the faces
   {
      auto nface = get_iface_vertices(jface); // get the nodes for the face
      const unsigned int node1 = nface[0];
      const unsigned int node2 = nface[1];
      auto esup_node1 = get_esup(node1);  // get elements surrounding node1
      auto esup_node2 = get_esup(node2);  // get elements surrounding node2
      std::vector<unsigned int> cell_node1(esup_node1.first);
      std::vector<unsigned int> cell_node2(esup_node2.first);
      for(unsigned int i=0; i<esup_node1.first; ++i)
         cell_node1[i] = esup_node1.second[i];
      for(unsigned int i=0; i<esup_node2.first; ++i)
         cell_node2[i] = esup_node2.second[i];

      // Sort the arrays, needed for intersection
      std::sort(cell_node1.begin(), cell_node1.end());
      std::sort(cell_node2.begin(), cell_node2.end());

      // Find the intersection
      std::vector<unsigned int> inter(std::min(esup_node1.first, esup_node2.first));
      std::set_intersection(cell_node1.begin(), cell_node1.end(),
                            cell_node2.begin(), cell_node2.end(), inter.begin());

      // get the nodes of the first cell (to check the orientation)
      auto cell = get_cell_vertices(inter[0]);
      auto a = get_coord(node1);
      auto b = get_coord(node2);
      const double* c = 0;
      // loop over the nodes of the element and get node c (other than a and b)
      for(unsigned int inode=0; inode<cell.first; ++inode)
      {
         if(cell.second[inode] != node1 && cell.second[inode] != node2)
         {
            c = get_coord(cell.second[inode]);
            break;
         }
      }
      if(orient(a, b, c) == 1)  // c is to the left of ab
      {
         iface_cell[2*jface] = inter[0];
         iface_cell[2*jface+1] = inter[1];
      }
      else  // c is to the right of ab
      {
         iface_cell[2*jface] = inter[1];
         iface_cell[2*jface+1] = inter[0];
      }
   }

   //Boundary Faces
   bface_cell = new unsigned int[n_bface];   // only 1 neighbouring cell
   for(unsigned int jface=0; jface<n_bface; ++jface)  // loop over all the faces
   {
      auto bf = get_bface_vertices(jface);   // get the nodes for jface
      const unsigned int node1 = bf[0];
      const unsigned int node2 = bf[1];

      auto esup_node1 = get_esup(node1);  // get elements surrounding node1
      auto esup_node2 = get_esup(node2);  // get elements surrounding node2
      std::vector<unsigned int> cell_node1(esup_node1.first);
      std::vector<unsigned int> cell_node2(esup_node2.first);
      for(unsigned int i=0; i<esup_node1.first; ++i)
         cell_node1[i] = esup_node1.second[i];
      for(unsigned int i=0; i<esup_node2.first; ++i)
         cell_node2[i] = esup_node2.second[i];

      // Sort the arrays, needed for intersection
      std::sort(cell_node1.begin(), cell_node1.end());
      std::sort(cell_node2.begin(), cell_node2.end());

      // find the intersection
      std::vector<unsigned int> inter(std::min(esup_node1.first, esup_node2.first));
      std::set_intersection(cell_node1.begin(), cell_node1.end(),
                            cell_node2.begin(), cell_node2.end(), inter.begin());

      bface_cell[jface] = inter[0];
      // get the nodes of the first cell (to check the orientation)
      auto cell = get_cell_vertices(inter[0]);
      auto a = get_coord(node1);
      auto b = get_coord(node2);
      const double* c = 0;
      // loop over the nodes of the element and get node c (other than a and b)
      for(unsigned int inode=0; inode<cell.first; ++inode)
      {
         if(cell.second[inode] != node1 && cell.second[inode] != node2)
         {
            c = get_coord(cell.second[inode]);
            break;
         }
      }
      if(orient(a, b, c) == 0)  // c is to the right of ab
      {
         // change the face direction so that the boundary cell is on the left
         std::swap(bface[2*jface], bface[2*jface+1]);
      }
   }

   has_esuf = true;
   cout << "Done\n";
}

// Find cells surrounding a cell
// type = 0 => face connected neighbours only, 1 => all neighbours including the node connected ones
void Grid::construct_esue(const esue_type type)
{
   if(type == esue_neumann && has_esue_neumann) return;
   if(type == esue_moore && has_esue_moore) return;

   if(type == esue_moore)  // Moore neighbours
   {
      construct_esup();
      std::cout << "Constructing cells surrounding cell: moore ... ";
   }
   else if(type == esue_neumann)  // Von-Neumann neighbours
   {
      construct_esuf();
      std::cout << "Constructing cells surrounding cell: neumann ... ";
   }

   // Delete if memory is already allocated
   // This is needed if e.g., we already have moore neighbours but now we
   // want neumann neighbours.
   if(!esue1) delete[] esue1;
   if(!esue2) delete[] esue2;

   esue2 = new unsigned int[n_cell+1];
   // Initialize esue2
   for(unsigned int i=0; i<=n_cell; ++i)
      esue2[i] = 0;

   if(type == esue_moore)  // Moore neighbours (we construct this similar to psup)
   {
      unsigned int* lelem = new unsigned int[n_cell]; // help array to avoid duplication
      for(unsigned int i=0; i<n_cell; ++i)
         lelem[i] = 0;  // Initialize lelem

      unsigned int istor = 0;
      std::vector<unsigned int> esue1_temp(0);  // temporary vector to enable resize function
      for(unsigned int icell=0; icell<n_cell; ++icell)   // loop over all the cells
      {
         auto cell = get_cell_vertices(icell);
         for(unsigned int inode=0; inode<cell.first; ++inode)  //loop over the nodes of icell
         {
            auto esup = get_esup(cell.second[inode]);  // get the cells surrounding inode
            for(unsigned int jcell=0; jcell<esup.first; ++jcell)  //loop over the cell surrounding inode
            {
               // check if not already stored (similar to what we have done for psup)
               if(esup.second[jcell] != icell && lelem[esup.second[jcell]] != icell+1)
               {
                  esue1_temp.resize(esue1_temp.size()+1);
                  esue1_temp[istor] = esup.second[jcell];
                  lelem[esup.second[jcell]] = icell+1;
                  ++istor;
               }
            }
         }
         esue2[icell+1] = istor;
      }
      // Copy the data from esue1_temp to esue1
      esue1 = new unsigned int[esue2[n_cell]];
      for(unsigned int i=0; i<esue2[n_cell]; ++i)
         esue1[i] = esue1_temp[i];

      esue1_temp.resize(0);
      delete [] lelem;
      has_esue_moore = true;
   }
   else  // Von-Neumann neighbours (we construct this similar to esup)
   {
      // construct esue2
      for(unsigned int iface=0; iface<n_iface; ++iface)  // loop over all the faces
      {
         auto face_cell = get_iface_cell(iface); // get the face neighbours
         ++esue2[face_cell[0]+1];  // add one neighbour to left cell
         ++esue2[face_cell[1]+1];  // add one neighbour to right cell
      }

      for(unsigned int i=1; i<=n_cell; ++i)
      {
         esue2[i] += esue2[i-1];
      }
      //construct esue1
      esue1 = new unsigned int[esue2[n_cell]];
      for(unsigned int iface=0; iface<n_iface; ++iface)
      {
         auto face_cell = get_iface_cell(iface);      // get the face neighbours
         unsigned int istor1 = esue2[face_cell[0]];   // left cell
         unsigned int istor2 = esue2[face_cell[1]];   // right cell
         esue1[istor1] = face_cell[1]; // store right cell as a neighbour of left cell
         esue1[istor2] = face_cell[0]; // store left cell as a neighbour of right cell
         ++esue2[face_cell[0]];
         ++esue2[face_cell[1]];
      }
      // reshuffle esue2
      for(unsigned int i=n_cell; i>0; --i)
         esue2[i] = esue2[i-1];
      esue2[0] = 0;
      has_esue_neumann = true;
   }

   std::cout << "Done\n";
}

// Compute the area of all cells
void Grid::compute_cell_area()
{
   carea = new double[n_cell];
   for(unsigned int icell=0; icell<n_cell; ++icell)
   {
      double area = 0;
      auto cell = get_cell_vertices(icell);
      auto p_n = get_coord(cell.second[cell.first-1]); // get coordinates of last point
      for(unsigned int ipoint=0; ipoint<cell.first-1; ++ipoint)
      {
         auto p_i = get_coord(cell.second[ipoint]);   // get coordinates of ipoint
         auto p_i_1 = get_coord(cell.second[ipoint+1]);  // get coordinates of next ipoint
         area += 0.5*((p_i[0]-p_n[0])*(p_i_1[1]-p_n[1])-(p_i[1]-p_n[1])*(p_i_1[0]-p_n[0]));
      }
      carea[icell] = std::abs(area);
   }
}

// Compute the unit normal and the length for all faces
void Grid::compute_face_normal()
{
   //Interior faces
   iface_norm = new double[2*n_iface];
   iface_len = new double[n_iface];
   for(unsigned int iface=0; iface<get_n_iface(); ++iface)
   {
      auto face = get_iface_vertices(iface); // get the nodes of iface
      auto node1 = get_coord(face[0]); // coordinates of first node
      auto node2 = get_coord(face[1]); // coordinates of second node
      //calculate length of the face
      iface_len[iface] = std::sqrt((node2[1]-node1[1])*(node2[1]-node1[1])
                                   + (node2[0]-node1[0])*(node2[0]-node1[0]));
      iface_norm[2*iface] = (node2[1]-node1[1])/iface_len[iface];
      iface_norm[2*iface+1] = (node1[0]-node2[0])/iface_len[iface];
   }

   // Boundary faces
   bface_norm = new double[2*n_bface];
   bface_len =new double[n_bface];
   for(unsigned int bface=0; bface<get_n_bface(); ++bface)
   {
      auto face = get_bface_vertices(bface); // get the nodes of bface
      auto node1 = get_coord(face[0]); // coordinates of first node
      auto node2 = get_coord(face[1]); // coordinates of second node
      // calculate length of the face
      bface_len[bface] = std::sqrt((node2[1]-node1[1])*(node2[1]-node1[1])
                                   + (node2[0]-node1[0])*(node2[0]-node1[0]));
      bface_norm[2*bface] = (node2[1]-node1[1])/bface_len[bface];
      bface_norm[2*bface+1] = (node1[0]-node2[0])/bface_len[bface];
   }
}

// Compute cell centroid
void Grid::compute_cell_centroid()
{
   cell_centroid = new double[2*n_cell];
   for(unsigned int icell=0; icell<n_cell; ++icell)
   {
      double x_centroid = 0;
      double y_centroid = 0;
      auto cell = get_cell_vertices(icell);  // get the nodes of the cell
      for(unsigned int inode=0; inode<cell.first; ++inode) //loop over the nodes of the cell
      {
         auto node = get_coord(cell.second[inode]);   // get the node coordinates
         // calculate the centroid
         x_centroid += node[0];
         y_centroid += node[1];
      }
      cell_centroid[2*icell] = x_centroid/cell.first;
      cell_centroid[2*icell+1] = y_centroid/cell.first;
   }
}

// Compute face centroid
void Grid::compute_face_centroid()
{
   // Interior faces
   iface_centroid = new double[2*n_iface];
   for(unsigned int iface=0; iface<n_iface; ++iface)
   {
      auto face = get_iface_vertices(iface); // get the nodes of iface
      auto node1 = get_coord(face[0]); // coordinates of node1
      auto node2 = get_coord(face[1]); // coordinates of node2
      iface_centroid[2*iface] = 0.5*(node1[0]+node2[0]);
      iface_centroid[2*iface+1] = 0.5*(node1[1]+node2[1]);
   }

   // Boundary faces
   bface_centroid = new double[2*n_bface];
   for(unsigned int bface=0; bface<n_bface; ++bface)
   {
      auto face = get_bface_vertices(bface); // get the nodes of bface
      auto node1 = get_coord(face[0]); // coordinates of node1
      auto node2 = get_coord(face[1]); // coordinates of node2
      bface_centroid[2*bface] = 0.5*(node1[0]+node2[0]);
      bface_centroid[2*bface+1] = 0.5*(node1[1]+node2[1]);
   }
}
