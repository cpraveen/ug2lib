#include <iostream>
#include <fstream>
#include <vector>
#include <utility> // std::pair, std::make_pair

using namespace std;

class Grid
{
public:
   Grid();
   ~Grid();
   void read_gmsh(const string grid_file);
   void write_vtk(const string grid_file);
   void construct_esup();
   void construct_psup(const bool all_points=true);
   void construct_iface();
   void construct_esuf();
   void construct_esue(const bool type=0);
   void compute_carea();

   inline unsigned int get_n_vertex()
   {
      return n_vertex;
   }

   inline unsigned int get_n_cell()
   {
      return n_cell;
   }

   inline unsigned int get_n_tri()
   {
      return n_tri;
   }

   inline unsigned int get_n_quad()
   {
      return n_quad;
   }

   inline unsigned int get_n_bface()
   {
      return n_bface;
   }

   inline unsigned int get_n_iface()
   {
      return n_iface;
   }

   inline const double* get_coord(unsigned int i)
   {
      return &coord[i*dim];
   }

   inline const double& get_cell_area(unsigned int i)
   {
      return carea[i];
   }

   // We are assuming linear edge
   inline const unsigned int* get_iface_vertices(unsigned int i)
   {
      return &iface[2*i];
   }

   inline const double& get_iface_length(unsigned int i)
   {
      return iface_len[i];
   }

   inline const double* get_iface_normal(unsigned int i)
   {
      return &iface_norm[i*dim];
   }

   inline const double* get_bface_normal(unsigned int i)
   {
      return &bface_norm[i*dim];
   }

   inline const unsigned int* get_bface_vertices(unsigned int i)
   {
      return &bface[2*i];
   }

   inline std::pair<unsigned int,const unsigned int*> get_cell_vertices(unsigned int i)
   {
      unsigned int start = cell2[i];
      unsigned int end = cell2[i+1];
      return std::make_pair(end-start,&cell1[start]);
   }

   inline std::pair<unsigned int,const unsigned int*> get_esup(unsigned int i)
   {
      unsigned int start = esup2[i];
      unsigned int end = esup2[i+1];
      return std::make_pair(end-start,&esup1[start]);
   }

   inline std::pair<unsigned int,const unsigned int*> get_psup(unsigned int i)
   {
      unsigned int start = psup2[i];
      unsigned int end = psup2[i+1];
      return std::make_pair(end-start,&psup1[start]);
   }

   inline std::pair<unsigned int,const unsigned int*> get_esue(unsigned int i)
   {
      unsigned int start = esue2[i];
      unsigned int end = esue2[i+1];
      return std::make_pair(end-start,&esue1[start]);
   }

   inline const unsigned int* get_iface_cell(unsigned int i)
   {
      return &iface_cell[2*i];
   }

   inline unsigned int get_bface_cell(unsigned int i)
   {
      return bface_cell[i];
   }

private:
   const int    dim = 2;
   unsigned int n_vertex, n_cell, n_tri, n_quad, n_bface, n_iface;
   double*       coord;

   // cell data
   unsigned int* cell1, *cell2;
   int*          ctype;
   double*       carea;

   // connectivity information
   unsigned int* esup1, *esup2;
   unsigned int* psup1, *psup2;
   unsigned int* esue1, *esue2;

   // boundary face data
   unsigned int* bface;      // vertices forming the face
   unsigned int* bface_cell; // cell adjacent to boundary face
   int*          bface_type; // type read from grid file, used for bc
   double*       bface_norm; // unit outward normal
   double*       bface_len;  // length of face

   // Interior faces
   double*       iface_len;  // length of face
   unsigned int* iface;      // vertex numbers for each face
   double*       iface_norm; // unit normal to face
   unsigned int* iface_cell; // cells adjacent to face

   bool         has_esup;
   bool         has_psup;
   bool         has_iface;
   bool         has_esuf;
   bool         has_esue;
};

// Calculates the cross product (ca x cb)
inline bool orient(const double* a, const double* b, const double* c)
{
   double ar = (a[0]-c[0])*(b[1]-c[1]) - (b[0]-c[0])*(a[1]-c[1]);
   if(ar < 0)
      return 0;  // c is to the right of ab vector
   else
      return 1;  // c is to the left of ab vector
}

