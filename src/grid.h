#ifndef __GRID_H__
#define __GRID_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <utility> // std::pair, std::make_pair

using namespace std;

namespace Gmsh
{

// Gmsh element types
// Linear elements
   const int Line2 = 1;
   const int Triangle3 = 2;
   const int Quadrilateral4 = 3;

// Quadratic elements
   const int Line3 = 8;
   const int Triangle6 = 9;
   const int Quadrilateral9 = 10;
}

class Grid
{
public:
   enum esue_type {esue_neumann, esue_moore, esue_none};
   enum psup_type {psup_all, psup_edge, psup_none};

   Grid();
   ~Grid();
   void read_gmsh(const string grid_file);
   void write_vtk(std::ofstream& vtk);
   void write_vtk(const string filename);
   void construct_esup();
   void construct_psup(const psup_type type);
   void construct_iface();
   void construct_esuf();
   void construct_esue(const esue_type type);
   void compute_cell_area();
   void compute_face_normal();
   void compute_cell_centroid();
   void compute_face_centroid();

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

   inline const double* get_cell_centroid(unsigned int i)
   {
      return &cell_centroid[2*i];
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

   inline const double* get_iface_centroid(unsigned int i)
   {
      return &iface_centroid[2*i];
   }

   inline const double* get_bface_centroid(unsigned int i)
   {
      return &bface_centroid[2*i];
   }

   inline const double* get_bface_normal(unsigned int i)
   {
      return &bface_norm[i*dim];
   }

   inline const unsigned int* get_bface_vertices(unsigned int i)
   {
      return &bface[2*i];
   }

   inline const double& get_bface_length(unsigned int i)
   {
      return bface_len[i];
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

   inline const unsigned int& get_bface_cell(unsigned int i)
   {
      return bface_cell[i];
   }

   inline esue_type get_esue_type()
   {
      if(has_esue_moore) return esue_moore;
      if(has_esue_neumann) return esue_neumann;
      return esue_none;
   }

   inline psup_type get_psup_type()
   {
      if(has_psup_all) return psup_all;
      if(has_psup_edge) return psup_edge;
      return psup_none;
   }

private:
   const int     dim = 2;
   unsigned int  n_vertex, n_cell, n_tri, n_quad, n_bface, n_iface;
   double*       coord;

   // cell data
   unsigned int* cell1, *cell2;
   int*          ctype;
   double*       carea;
   double*       cell_centroid;  //centroid of the cells

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
   double*       bface_centroid; //centroid of bface

   // Interior faces
   double*       iface_len;  // length of face
   unsigned int* iface;      // vertex numbers for each face
   double*       iface_norm; // unit normal to face
   unsigned int* iface_cell; // cells adjacent to face
   double*       iface_centroid; //centroid of iface

   bool         has_esup;
   bool         has_psup_all;
   bool         has_psup_edge;
   bool         has_iface;
   bool         has_esuf;
   bool         has_esue_moore;
   bool         has_esue_neumann;
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

void grid_test1(); // Function to dump the constructed data (for reference mesh)

#endif
