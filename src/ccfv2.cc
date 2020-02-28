#include <iostream>
#include <cmath>
#include <cassert>
#include <cstring>
#include <fstream>

#include "grid.h"
#include "vtk_writer.h"

using namespace std;

// return euclidean norm of 2d vector
double norm(const double* a)
{
   return sqrt(pow(a[0],2) + pow(a[1],2));
}

// return euclidean distance b/w 2d vectors
double distance(const double* a, const double* b)
{
   return sqrt(pow(a[0]-b[0],2) + pow(a[1]-b[1],2));
}

// return dot product of 2d vectors
double dot(const double* a, const double* b)
{
   return a[0]*b[0] + a[1]*b[1];
}

// advection velocity in linear advection problem
void advection_velocity(const double* x, double* v)
{
   v[0] =  x[1];
   v[1] = -x[0];
}

// initial condition as a function of x,y
double initial_condition(const double* x)
{
   double r2 = pow(x[0]-0.5,2) + pow(x[1],2);
   return exp(-50.0*r2);
}

// upwind flux
double num_flux(const double ul, const double ur,
                const double* v, const double* normal)
{
   double vn = dot(v, normal);
   return (vn > 0.0) ? (vn*ul) : (vn*ur);
}

struct Parameters
{
   string grid_file;
   double tfinal;
   int save_freq;
   string esue_type;
   int order;
};

// class for finite volume method
class FVM
{
public:
   FVM(Parameters& param);
   ~FVM();
   void run();

private:
   void read_grid_and_preprocess();
   void allocate_memory();
   void set_initial_condition();
   void compute_time_step();
   void compute_gradient();
   void compute_residual();
   void update_solution_rk1();
   void update_solution_rk2();
   void save_solution();
   void compute_lscoef();
   void compute_error();

   Parameters*  param;
   Grid         grid;
   double*      uold; // solution at time level n
   double*      u;    // solution at time level n+1
   double*      res;  // residual
   double*      du;   // gradient of u
   double*      lscoef; // coefficients for least squares

   double       t, dt;
   unsigned int iter;
};

// constructor
FVM::FVM(Parameters& param)
   :
   param(&param)
{
   t = 0.0;
   iter = 0;
}

// destructor
FVM::~FVM()
{
   if(!u)      delete[] u;
   if(!uold)   delete[] uold;
   if(!res)    delete[] res;
   if(!du)     delete[] du;
   if(!lscoef) delete[] lscoef;
}

void FVM::read_grid_and_preprocess()
{
   grid.read_gmsh(param->grid_file);
   grid.construct_esuf();
   if(param->esue_type == "neumann")
      grid.construct_esue(Grid::esue_neumann);
   else if(param->esue_type == "moore")
      grid.construct_esue(Grid::esue_moore);
   else
   {
      cout << "Invalid esue_type, options: moore/neumann\n";
      exit(0);
   }
   grid.compute_cell_area();
   grid.compute_cell_centroid();
   grid.compute_face_normal();
   grid.compute_face_centroid();
}

void FVM::allocate_memory()
{
   auto n_cell = grid.get_n_cell();

   u = new double[n_cell];
   uold = new double[n_cell];
   res = new double[n_cell];
   du = new double[2*n_cell];

   // memory needed for storing least squares coefficients
   auto range = grid.get_esue_range(n_cell-1);
   lscoef = new double[2*range.second];
}

// set initial condition into solution array
void FVM::set_initial_condition()
{
   auto n_cell = grid.get_n_cell();
   for(unsigned int i=0; i<n_cell; ++i)
   {
      auto x = grid.get_cell_centroid(i);
      u[i] = initial_condition(x);
   }
}

// compute dt using cfl condition based on maximum principle
void FVM::compute_time_step()
{
   auto n_cell = grid.get_n_cell();
   vector<double> tmp(n_cell);
   for(unsigned int i=0; i<n_cell; ++i)
      tmp[i] = 0.0;

   auto n_iface = grid.get_n_iface();
   for(unsigned int i=0; i<n_iface; ++i)
   {
      auto xf = grid.get_iface_centroid(i);
      double v[2];
      advection_velocity(xf, v);
      auto normal = grid.get_iface_normal(i);
      auto vn = dot(v, normal);
      auto ds = grid.get_iface_length(i);
      auto cell = grid.get_iface_cell(i);
      if(vn > 0)
         tmp[cell[0]] += vn * ds;
      else
         tmp[cell[1]] -= vn * ds;
   }

   auto n_bface = grid.get_n_bface();
   for(unsigned int i=0; i<n_bface; ++i)
   {
      auto xf = grid.get_bface_centroid(i);
      double v[2];
      advection_velocity(xf, v);
      auto normal = grid.get_bface_normal(i);
      auto vn = dot(v, normal);
      auto ds = grid.get_bface_length(i);
      auto cell = grid.get_bface_cell(i);
      if(vn > 0) tmp[cell] += vn * ds;
   }

   dt = 1.0e20;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      auto area = grid.get_cell_area(i);
      dt = fmin(dt, area/(tmp[i]+1.0e-14));
   }
   cout << "Time step = " << dt << endl;
}

// compute least squares coefficients
void FVM::compute_lscoef()
{
   auto n_cell = grid.get_n_cell();
   for(unsigned int i=0; i<n_cell; ++i)
   {
      auto cent_i = grid.get_cell_centroid(i);
      auto esue = grid.get_esue(i); // get surrounding cells
      double a11 = 0, a12 = 0, a22 = 0; // elements of A matrix
      double* w = new double[esue.first]; // weights matrix W (1/r_ij)
      double* dx = new double[esue.first];
      double* dy = new double[esue.first]; // D matrix

      // compute the elements of A, W and D matrix
      for(unsigned int j=0; j<esue.first; ++j)
      {
         auto cent_j = grid.get_cell_centroid(esue.second[j]);
         w[j] = 1.0/distance(cent_i, cent_j);
         dx[j] = cent_j[0] - cent_i[0];
         dy[j] = cent_j[1] - cent_i[1];
         a11 += pow(w[j]*dx[j],2);
         a22 += pow(w[j]*dy[j],2);
         a12 += pow(w[j],2)*dx[j]*dy[j];
      }

      // compute lscoef
      auto range = grid.get_esue_range(i);
      double det = a11*a22 - a12*a12;
      assert(det > 0);
      for(unsigned int j=0; j<esue.first; ++j)
      {
         unsigned int k = 2*(range.first + j);
         lscoef[k]   = pow(w[j],2)*(a22*dx[j] - a12*dy[j])/det; // alpha
         lscoef[k+1] = pow(w[j],2)*(a11*dy[j] - a12*dx[j])/det; // beta
      }

      delete[] w;
      delete[] dx;
      delete[] dy;
   }
}

// compute gradient based on least squares approach
void FVM::compute_gradient()
{
   auto n_cell = grid.get_n_cell();
   if(param->order == 1)
   {
      for(unsigned int i=0; i<n_cell; ++i)
      {
         du[2*i]   = 0;
         du[2*i+1] = 0;
      }
      return;
   }

   for(unsigned int i=0; i<n_cell; ++i)
   {
      du[2*i] = du[2*i+1] = 0; // Initialize
      auto esue = grid.get_esue(i); // get surrounding cells
      auto range = grid.get_esue_range(i);
      for(unsigned int j=0; j<esue.first; ++j)
      {
         unsigned int k = 2*(range.first + j); // location in lscoef
         double dU = u[esue.second[j]] - u[i];
         du[2*i]   += lscoef[k]*dU;   // u_x
         du[2*i+1] += lscoef[k+1]*dU; // u_y
      }
   }
}

// compute finite volume residual R in the semi-discrete equation
// A*du/dt + R = 0, A = cell area
void FVM::compute_residual()
{
   compute_gradient();

   auto n_cell = grid.get_n_cell();
   for(unsigned int i=0; i<n_cell; ++i)
      res[i] = 0.0;

   // interior faces
   auto n_iface = grid.get_n_iface();
   for(unsigned int i=0; i<n_iface; ++i)
   {
      auto cell = grid.get_iface_cell(i);

      // linear reconstruction of solution at face
      auto xf = grid.get_iface_centroid(i);
      auto xl = grid.get_cell_centroid(cell[0]); // centroid of left cell
      auto xr = grid.get_cell_centroid(cell[1]); // right cell
      double rl[2], rr[2]; // vector from cell to face centroid
      rl[0] = xf[0] - xl[0];
      rl[1] = xf[1] - xl[1];
      rr[0] = xf[0] - xr[0];
      rr[1] = xf[1] - xr[1];
      auto ul = u[cell[0]] + dot(&du[2*cell[0]], rl);
      auto ur = u[cell[1]] + dot(&du[2*cell[1]], rr);

      double v[2];
      advection_velocity(xf, v);
      auto normal = grid.get_iface_normal(i);
      double flux = num_flux(ul, ur, v, normal);
      auto ds = grid.get_iface_length(i);
      res[cell[0]] += flux * ds;
      res[cell[1]] -= flux * ds;
   }

   // boundary faces
   auto n_bface = grid.get_n_bface();
   for(unsigned int i=0; i<n_bface; ++i)
   {
      auto cell = grid.get_bface_cell(i);

      // linear reconstruction
      auto xf = grid.get_bface_centroid(i);
      auto xl = grid.get_cell_centroid(cell);
      double rl[2];
      rl[0] = xf[0] - xl[0];
      rl[1] = xf[1] - xl[1];
      auto ul = u[cell] + dot(&du[2*cell], rl);
      auto ur = 0.0;

      double v[2];
      advection_velocity(xf, v);
      auto normal = grid.get_bface_normal(i);
      double flux = num_flux(ul, ur, v, normal);
      auto ds = grid.get_bface_length(i);
      res[cell] += flux * ds;
   }
}

// save solution to file with different name
void FVM::save_solution()
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
   VTKWriter writer(filename, grid, t, iter);
   writer.write_cell_scalar(u,"u");
   writer.close();
   ++counter;
}

// updates solution with 1st order forward Euler method
void FVM::update_solution_rk1()
{
   auto n_cell = grid.get_n_cell();
   compute_residual();
   for(unsigned int i=0; i<n_cell; ++i)
      u[i] = u[i] - (dt/grid.get_cell_area(i)) * res[i];
}

// updates solution with 2nd order Runge-Kutta method (Heun's method)
void FVM::update_solution_rk2()
{
   auto n_cell = grid.get_n_cell();

   // 1st stage
   compute_residual();
   for(unsigned int i=0; i<n_cell; ++i)
      u[i] = uold[i] - (dt/grid.get_cell_area(i)) * res[i];

   // 2nd stage
   compute_residual();
   for(unsigned int i=0; i<n_cell; ++i)
      u[i] = 0.5*uold[i] + 0.5*(u[i] - (dt/grid.get_cell_area(i)) * res[i]);
}

// compute l1 and l2 norm of error
void FVM::compute_error()
{
   double l1error = 0.0;
   double l2error = 0.0;
   double total_area = 0.0;
   auto n_cell = grid.get_n_cell();
   for(unsigned int i=0; i<n_cell; ++i)
   {
      auto x = grid.get_cell_centroid(i);
      auto uexact = initial_condition(x);
      auto area = grid.get_cell_area(i);
      l1error += abs(uexact - u[i]) * area;
      l2error += pow(uexact-u[i],2) * area;
      total_area += area;
   }
   l1error /= total_area;
   l2error = sqrt(l2error/total_area);
   cout << "l1 norm, l2 norm = " << l1error << " " << l2error << endl;
}

void FVM::run()
{
   read_grid_and_preprocess();
   allocate_memory();
   if(param->order == 2) compute_lscoef();
   set_initial_condition();
   save_solution();
   compute_time_step();

   while(t < param->tfinal)
   {
      if(t+dt > param->tfinal) dt = param->tfinal - t;
      for(unsigned int i=0; i<grid.get_n_cell(); ++i)
         uold[i] = u[i];
      if(param->order == 1)
         update_solution_rk1();
      else
         update_solution_rk2();
      t += dt;
      ++iter;
      cout << "iter, t = " << iter << " " << t << endl;
      if(iter % param->save_freq == 0 || abs(t-param->tfinal) < 1.0e-14)
         save_solution();
   }

   compute_error();
}

void read_parameters(Parameters& param, char* param_file)
{
   ifstream file(param_file);
   assert(file.is_open());
   file >> param.grid_file
        >> param.tfinal
        >> param.save_freq
        >> param.esue_type
        >> param.order;
   file.close();
}

int main(int argc, char* argv[])
{
   assert(argc == 2);
   Parameters param;
   read_parameters(param, argv[1]);
   FVM fvm(param);
   fvm.run();
}
