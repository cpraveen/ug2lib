#ifndef __GRID_H__
#define __GRID_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <utility> // std::pair, std::make_pair

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
   enum esue_type {esue_neumann, esue_moore, esue_symmetric, esue_none};
   enum psup_type {psup_all, psup_edge, psup_none};

   Grid();
   ~Grid();
   void read_gmsh(const std::string grid_file);
   void write_vtk(std::ofstream& vtk);
   void write_vtk(const std::string filename);
   void write_msh(const std::string filename);
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

   inline unsigned int get_n_total_vertex()
   {
      return n_total_vertex;
   }

   inline unsigned int get_n_cell()
   {
      return n_cell;
   }

   inline unsigned int get_n_total_cell()
   {
      return n_total_cell;
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

   inline int get_cell_type(unsigned int i)
   {
      return ctype[i];
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

   inline std::pair<unsigned int,unsigned int> get_esue_range(unsigned int i)
   {
      return std::make_pair(esue2[i],esue2[i+1]);
   }

   inline const unsigned int* get_iface_cell(unsigned int i)
   {
      return &iface_cell[2*i];
   }

   inline const unsigned int get_bface_cell(unsigned int i)
   {
      return bface_cell[i];
   }

   inline unsigned int map_ghost_cell(unsigned int i)
   {
      return ghost_cell[i-n_cell];
   }

   inline unsigned int map_ghost_vertex(unsigned int i)
   {
      return ghost_vertex[i-n_vertex];
   }

   inline bool periodic()
   {
      return has_periodic;
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

   inline bool cell_is_real(unsigned int i)
   {
      return (i < n_cell) ? true : false;
   }

   inline bool vertex_is_real(unsigned int i)
   {
      return (i < n_vertex) ? true : false;
   }

private:
   const int     dim = 2;
   unsigned int  n_vertex, n_cell, n_tri, n_quad, n_bface, n_iface;
   unsigned int  n_ghost_cell, n_ghost_vertex, n_total_cell, n_total_vertex;
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

   // Periodic boundary data
   unsigned int* periodicity;  // periodicity of nodes (0, 1, 2)
   unsigned int* prd_node;     // map of periodic nodes
   unsigned int* ghost_cell;   // map of ghost and real cells
   unsigned int* ghost_vertex; // map of ghost and real vertices

   bool         has_esup;
   bool         has_psup_all;
   bool         has_psup_edge;
   bool         has_iface;
   bool         has_esuf;
   bool         has_esue_moore;
   bool         has_esue_neumann;
   bool         has_esue_symmetric;
   bool         has_periodic;

   static double distance(const double* a, const double* b);
   static double dot(const double* a, const double* b);
   static bool orient(const double* a, const double* b, const double* c);
};

void grid_test1(); // Function to dump the constructed data (for reference mesh)

#endif
