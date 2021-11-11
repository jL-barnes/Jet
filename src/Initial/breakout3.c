
#include "../paul.h"

static double gam;
static double exp_energy;
static double rhow;

void setICparams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
   exp_energy = theDomain->theParList.Explosion_Energy;
   rhow = theDomain->theParList.Gam_0;
}

void initial( double * prim , double * x ){

   double M0 = 1.0;
   double R0 = 1.0;
   double r  = x[0];

   double rho0 = 0.0615*M0/R0/R0/R0;
   double rho_wind = rhow*M0/pow(R0,3.);
   double n = 3.5; 

   double R_cav = 1.5e-3/1.6;
   double rho_cav_boundary = rho0 * pow(R0/R_cav,2.65) * pow(1-R_cav/R0,3.5 );
   double R_atm = 0.02*R0;
   double Rcut = 0.98*R0;

   double rho = rho0*pow(R0/r,2.65);
   rho *= pow(1.-r/R0,n);
   double rhocut = rho0*pow(R0/Rcut,2.65);
   rhocut *= pow(1.-Rcut/R0,n);
   if( r>Rcut ) rho = rho_wind*pow(R0/r,2.) + rhocut*exp(-(r-Rcut)/R_atm );
   else if (r < R_cav){ 
     rho *= 1.0e-3; }// * rho_cav_boundary; }

   double Pp = 1e-6*rho;
   if( Pp > 1.*M0/R0/R0/R0 ) Pp = 1.*M0/R0/R0/R0;

   double Rexp = 1e-2*R0;
   double E = exp_energy;
   double Pexp = (gam-1.)*E*exp(-r*r/Rexp/Rexp)/pow(sqrt(M_PI)*Rexp,3.);

   Pp += Pexp;
  
   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}
