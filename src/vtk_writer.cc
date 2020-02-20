#include "vtk_writer.h"

// constructor
VTKWriter::VTKWriter(string filename, Grid& grid)
:
filename (filename),
grid (&grid)
{
   vtk.open(filename);
   assert(vtk.is_open());
   write_grid();
}

// destructor
VTKWriter::~VTKWriter()
{
}

// Write grid to file
void VTKWriter::write_grid()
{
   auto n_vertex = grid->get_n_vertex();
   auto n_cell = grid->get_n_cell();
   auto n_tri = grid->get_n_tri();
   auto n_quad = grid->get_n_quad();

   // Save vtk file
   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << "Test file" << endl;
   vtk << "ASCII" << endl;
   vtk << "DATASET UNSTRUCTURED_GRID" << endl;
   vtk << "POINTS  " << n_vertex << "  float" << endl;

   for(unsigned int i=0; i<n_vertex; ++i)
   {
      auto x = grid->get_coord(i);
      vtk << x[0] << " " << x[1] << " " << 0.0 << endl;
   }

   vtk << "CELLS  " << n_cell << " " << 4 * n_tri + 5 * n_quad << endl;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      auto cell = grid->get_cell_vertices(i);
      vtk << cell.first;
      for(unsigned int j=0; j<cell.first; ++j)
         vtk << "  " << cell.second[j];
      vtk << endl;
   }

   vtk << "CELL_TYPES  " << n_cell << endl;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      auto ctype = grid->get_cell_type(i);
      if(ctype == 2) vtk << 5 << endl;
      else if(ctype == 3) vtk << 9 << endl;
      else
      {
         cout << "Unknown ctype = " << ctype << endl;
         exit(0);
      }
   }
}

// Write cell centered scalar data
void VTKWriter::write_cell_scalar(const double *data, string name)
{
   auto n_cell = grid->get_n_cell();

   vtk << "CELL_DATA " << n_cell << std::endl;
   vtk << "SCALARS " << name << " double 1\n";
   vtk << "LOOKUP_TABLE default\n";
   for (unsigned int i=0; i<n_cell; ++i)
      vtk << data[i] << std::endl;
}

void VTKWriter::close()
{
   vtk.close();
   std::cout << "Wrote file " << filename << std::endl;
}
