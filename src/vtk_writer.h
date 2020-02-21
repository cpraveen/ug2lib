#ifndef __VTK_WRITER_H__
#define __VTK_WRITER_H__

#include "grid.h"

class VTKWriter
{
public:
   VTKWriter(std::string filename, Grid& grid);
   ~VTKWriter();

   void write_cell_scalar(const double *data, std::string name);
   void close();

private:
   void write_grid();

   std::string   filename;
   Grid*         grid;
   std::ofstream vtk;
};
#endif
