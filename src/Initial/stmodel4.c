
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

double mass_appx( double r , double rin ){
   double f1 = 2.*(pow(r,.35)-pow(rin,.35));
   double n = 6.;
   double f = pow( pow(f1,-n) + 1. , -1./n );
   if( r<rin ) f = 0.0;
   if( r>1. ) f = 1.;
   return(f);
}

void initial( double * prim , double * x ){

   double M0 = 1.0;
   double R0 = 1.0;
   double r  = x[0];

   double rho0 = 0.0615*M0/R0/R0/R0;
   double rho_wind = 1e-9*M0/pow(R0,3.);
   double n = 3.5; 

   double R_cav = 1.5e-3/1.6;
   //double rho_cav_boundary = rho0 * pow(R0/R_cav,2.65) * pow(1-R_cav/R0,3.5 );

   double rho;
   if( r <= R0 ){ // within the boundary of the star
     rho = rho0*pow(R0/r,2.65);
     rho *= pow(1.0-r/R0,n);
     rho += rho_wind;
     if( r < R_cav ){ // within the cavity
       //   rho = rho_cav_boundary * 1e-3;
       rho *= 1e-3;
     }
   }
   else rho = rho_wind*pow(R0/r,2.);

   double Pp = 1e-6*rho;
   if( Pp > 1.*M0/R0/R0/R0 ) Pp = 1.*M0/R0/R0/R0;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

   // first passive scalar is the maximum Pressure
   if( NUM_N > 0 ) prim[NUM_C] = 0.0;
   // second passive scalar is the initial mass coordinate
   if( NUM_N > 1 ){
     double rin = 0;
      double X = mass_appx( r , rin );
      prim[NUM_C+1] = X;
   }
   // third passive scalar is the Ni-56 mass fraction 
   // (calculated according to the input pressure threshold)
   if( NUM_N > 2 ){
      prim[NUM_C+2] = 0.0;
   }

}
