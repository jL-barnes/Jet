enum{RHO,PPP,UU1,UU2,UU3}; // density, pressure, velocity vector
enum{DEN,TAU,SS1,SS2,SS3}; // mass, energy-rest mass energy, momentum vector
enum{C_FIXED,C_WRIEMANN,C_WCELL};

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define MOVE_CELLS C_WCELL

#define NUM_C 4 // in two dimensions, you only use 2 velocity/momentum coordinates
#define NUM_N 4 // number of extra things you're keeping track of
// For the jet-drive sne: We have three passive scalars: max pressure, initial mass coordinate,
// and X56 (calculated in situ based on pressure thresholds)
// For collapsar mixing: we have only one (whether or not it's a wind) and I'm giving that it's own index,
// rather than reusing the collapsar model
#define NUM_Q (NUM_C+NUM_N)
#define NUM_G 2

struct param_list{
   double t_min, t_max;
   int Num_R, Num_T, Num_P;
   int NumRepts, NumSnaps, NumChecks;
   int Out_LogTime;

   double rmin, rmax;
   double thmin, thmax;
   double phimax;

   double LogZoning, MaxShort, MaxLong, Target_X, Target_W, ShockPos;
   int Absorb_BC,Move_BCs,Initial_Regrid,Initial_Cons;
   int Reset_Entropy,Add_Cooling,Make_Nickel;
   double P56_critical;

   double CFL, PLM;
   double Density_Floor, Pressure_Floor;

   double Adiabatic_Index;
   int Gravity_Switch, Output_Mass;
   double PointMass;

  int Nozzle_Switch, Nozzle_is_Wind;
   double Nozzle_Power, Nozzle_Gamma, Nozzle_Eta, Nozzle_r0, Nozzle_th0, Nozzle_Time;
  double gamH_rel_Th;

  double Wind_Nozzle_Beta, Wind_Mass, Wind_t0, Wind_dt;
  int Start_Wind_tmin;

   double Explosion_Energy,Gam_0,Gam_Boost;

   int restart_flag;
};

struct domain{
   struct cell ** theCells;
   int * Nr;
   int Nt,Np,Ng;
   double * t_jph;
   double * p_kph;

   double g_point_mass;

   time_t Wallt_init;
   int count_steps;

   int rank,size;
   int dim_rank[2];
   int dim_size[2];
   int left_rank[2];
   int right_rank[2];
   MPI_Comm theComm;

   struct param_list theParList;

   double t;
   double t_init, t_fin;
   int nrpt;
   int N_rpt;
   int nsnp;
   int N_snp;
   int nchk;
   int N_chk;

   int final_step;
};

struct cell{
   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double grad[NUM_Q];
   double gradr[NUM_Q];
   double riph;
   double RKriph;
   double dr;
   double wiph;
};

struct face{
   struct cell * L;
   struct cell * R;
   double dxL;
   double dxR;
   double cm[3];
   double dr;
   double dA;
};

