#ifndef __VTK_WRITER_H__
#define __VTK_WRITER_H__

#include "grid.h"

class VTKWriter
{
public:
   VTKWriter(std::string filename, Grid& grid,
             double time=0.0, int cycle=0);
   ~VTKWriter();

   void write_cell_scalar(const double *data, std::string name);
   void close();

private:
   void write_grid();

   std::string   filename;
   Grid*         grid;
   double        time;
   int           cycle;
   std::ofstream vtk;
};
#endif
