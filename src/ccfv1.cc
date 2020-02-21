// Cell centered first order FV for linear advection in 2d
#include <cmath>
#include <algorithm>

#include "grid.h"
#include "vtk_writer.h"

using namespace std;

double norm(const double *a)
{
   return sqrt(pow(a[0],2) + pow(a[1],2));
}

double distance(const double *a, const double *b)
{
   return sqrt(pow(a[0]-b[0],2) + pow(a[1]-b[1],2));
}

double dot(const double*a, const double* b)
{
   return a[0]*b[0] + a[1]*b[1];
}

void advection_velocity(const double *x, double *v)
{
   v[0] =  x[1];
   v[1] = -x[0];
}

double initial_condition(const double *x)
{
   double r2 = pow(x[0]-0.5,2) + pow(x[1],2);
   return exp(-50.0*r2);
}

double compute_time_step(Grid& grid)
{
   auto n_cell = grid.get_n_cell();
   vector<double> hmin(n_cell);
   for(unsigned int i=0; i<n_cell; ++i)
      hmin[i] = 1.0e20;

   auto n_iface = grid.get_n_iface();
   for(unsigned int i=0; i<n_iface; ++i)
   {
      auto cell = grid.get_iface_cell(i);
      auto xl = grid.get_coord(cell[0]);
      auto xr = grid.get_coord(cell[1]);
      auto xf = grid.get_iface_centroid(i);
      auto hl = distance(xl, xf);
      auto hr = distance(xr, xf);
      hmin[cell[0]] = fmin(hmin[cell[0]], hl);
      hmin[cell[1]] = fmin(hmin[cell[1]], hr);
   }

   double dt = 1.0e20;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      auto x = grid.get_cell_centroid(i);
      double v[2];
      advection_velocity(x, v);
      dt = fmin(dt, hmin[i]/(norm(v)+1.0e-14));
   }
   cout << "Time step = " << dt << endl;
   return dt;
}

double num_flux(const double ul, const double ur,
                const double* v, const double* normal)
{
   double vn = dot(v, normal);
   return (vn > 0.0) ? (vn*ul) : (vn*ur);
}

void compute_residual(Grid& grid, const double* u, double* R)
{
   auto n_cell = grid.get_n_cell();
   for(unsigned int i=0; i<n_cell; ++i)
      R[i] = 0.0;

   // interior faces
   auto n_iface = grid.get_n_iface();
   for(unsigned int i=0; i<n_iface; ++i)
   {
      auto cell = grid.get_iface_cell(i);
      auto ul = u[cell[0]];
      auto ur = u[cell[1]];
      auto x = grid.get_iface_centroid(i);
      double v[2];
      advection_velocity(x, v);
      auto normal = grid.get_iface_normal(i);
      double flux = num_flux(ul, ur, v, normal);
      auto ds = grid.get_iface_length(i);
      R[cell[0]] += flux * ds;
      R[cell[1]] -= flux * ds;
   }

   //boundary faces
   auto n_bface = grid.get_n_bface();
   for(unsigned int i=0; i<n_bface; ++i)
   {
      auto cell = grid.get_bface_cell(i);
      auto ul = u[cell];
      auto ur = 0.0;
      auto x = grid.get_bface_centroid(i);
      double v[2];
      advection_velocity(x, v);
      auto normal = grid.get_bface_normal(i);
      double flux = num_flux(ul, ur, v, normal);
      auto ds = grid.get_bface_length(i);
      R[cell] += flux * ds;
   }
}

// update solution using forward euler
void update_solution(Grid& grid, double* u, const double* R, const double dt)
{
   auto n_cell = grid.get_n_cell();
   for(unsigned int i=0; i<n_cell; ++i)
      u[i] = u[i] - (dt/grid.get_cell_area(i)) * R[i];
}

// set initial condition into solution array
void set_initial_condition(Grid& grid, double* u)
{
   auto n_cell = grid.get_n_cell();
   for(unsigned int i=0; i<n_cell; ++i)
   {
      auto x = grid.get_cell_centroid(i);
      u[i] = initial_condition(x);
   }
}

// save solution to file with different name
void save_solution(Grid& grid, const double* u)
{
   static int counter = 0;
   string filename = "sol";
   if(counter < 10)
      filename += "00" + std::to_string(counter);
   else if(counter < 100)
      filename += "0" + std::to_string(counter);
   else if(counter < 1000)
      filename += std::to_string(counter);
   else
   {
      cout << "counter is too large !!!\n";
      exit(0);
   }

   filename += ".vtk";
   VTKWriter writer(filename, grid);
   writer.write_cell_scalar(u,"u");
   ++counter;
}

int main()
{
   Grid grid;
   grid.read_gmsh("ccfv1.msh");

   grid.construct_esuf();
   grid.compute_cell_area();
   grid.compute_cell_centroid();
   grid.compute_face_normal();
   grid.compute_face_centroid();

   auto n_cell = grid.get_n_cell();
   double* u = new double[n_cell];

   // Set initial condition
   set_initial_condition(grid, u);
   save_solution(grid, u);

   double dt = 0.40 * compute_time_step(grid);
   double* R = new double[n_cell];
   double t = 0.0;
   unsigned int iter = 0;
   while(t < 2.0*M_PI)
   {
      compute_residual(grid, u, R);
      update_solution(grid, u, R, dt);
      t += dt; ++iter;
      cout << "iter, t = " << iter << " " << t << endl;
      if(iter%100 == 0) save_solution(grid, u);
   }
}
