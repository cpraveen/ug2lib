// return euclidean norm of 2d vector
double norm(const double *a)
{
   return sqrt(pow(a[0],2) + pow(a[1],2));
}

// return euclidean distance b/w 2d vectors
double distance(const double *a, const double *b)
{
   return sqrt(pow(a[0]-b[0],2) + pow(a[1]-b[1],2));
}

// return dot product of 2d vectors
double dot(const double*a, const double* b)
{
   return a[0]*b[0] + a[1]*b[1];
}

// advection velocity in linear advection problem
void advection_velocity(const double *x, double *v)
{
   v[0] =  x[1];
   v[1] = -x[0];
}

// initial condition as a function of x,y
double initial_condition(const double *x)
{
   double r2 = pow(x[0]-0.5,2) + pow(x[1],2);
   return exp(-50.0*r2);
}

// class for finite volume method
class FVM
{
public:
   FVM();
   ~FVM();
   void run();

private:
   void read_grid_and_preprocess();
   void allocate_memory();
   void set_initial_condition();
   void compute_time_step();
   void compute_gradient();
   void compute_residual();
   void update_solution();
   void save_solution();

   Grid         grid;
   double*      uold; // solution at time level n
   double*      u;    // solution at time level n+1
   double*      res;  // residual
   double*      du;   // gradient of u
   double*      lscoef; // coefficients for least squares

   double       t, dt;
   unsigned int iter;
};

void FVM::allocate_memory()
{
   auto n_cell = grid.get_n_cell();

   u = new double[n_cell];
   uold = new double[n_cell];
   res = new double[n_cell];
   du = new double[2*n_cell];

   // TODO: count memory needed for storing least squares coefficients
   unsigned int c = 0;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      auto cell = grid.get_esue(i);
      c += cell.first;
   }
   lscoef = new double[2*c];
}

void FVM::run()
{
}

int main()
{
   FVM fvm;
   fvm.run();
}
