#include "paul.h"

void report( struct domain * theDomain , double t ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Ng = theDomain->Ng;
   int * Nr = theDomain->Nr;
   int rank = theDomain->rank;
   int size = theDomain->size;
   MPI_Comm grid_comm = theDomain->theComm;
   double * t_jph = theDomain->t_jph;
   double Pth = theDomain->theParList.P56_critical;
   double ghth = theDomain->theParList.gamH_rel_Th;


   int jmin = Ng;
   int jmax = Nt-Ng;
   if( rank==0 ) jmin = 0;
   if( rank==size-1 ) jmax = Nt;

   double L1 = 0.0;
   double uSum = 0.0;
   double uXSum = 0.0;
   double ESum = 0.0;
   double ETherm = 0.0;
   double EXSum = 0.0;
   double EJet = 0.0;
   double XSum = 0.0;
   double MSum = 0.0;
   double Ni_direct   = 0.0;
   double Ni_viaPPP = 0.0;
   double uMax = 0.0;
   double gh_max = 0.0;
   double rMax = 0.0;
   double rMin = HUGE_VAL;
   double rMin_rel_N = HUGE_VAL;
   double rMax_rel_N = 0.0;
   double rMin_rel_S = HUGE_VAL;
   double rMax_rel_S = 0.0;
   double I_th = 0.0;
   int i,j;
   for( j=jmin ; j<jmax ; ++j ){
      int jk = j;
      double dOmega = 2.*M_PI*(cos(t_jph[j-1])-cos(t_jph[j]));
      double dE = 0.0;
      for( i=0 ; i<Nr[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         double ur = c->prim[UU1];
         double up = c->prim[UU2];
         double X = 0.0; // X56 calculated with method 1
	 double XP = 0.0; // max pressure for calculation
	 // currently, the third passive scalar is the one
	 // that tracks nickel mass fraction
	 if( NUM_Q > NUM_C ) XP = c->prim[NUM_C];  
	 if( NUM_Q > NUM_C+2 ) X = c->prim[NUM_C+2];
         double u = sqrt(ur*ur+up*up);
         double E = c->cons[TAU];
         double M = c->cons[DEN];
         double h = 1.+4.*c->prim[PPP]/c->prim[RHO];
         double gam = sqrt( 1. + u*u );
	 double gh = gam*h;
	 if (gh >= ghth){ //relativistic flow
	   double r = c->riph;
	   if (t_jph[jk] <= M_PI/2.0){ // Northern hemisphere
	     if (r < rMin_rel_N) rMin_rel_N = r;
	     if (r > rMax_rel_N) rMax_rel_N = r;
	   }
	   else{ // Southern hemisphere
	     printf("in south (th = %.2f):\n", t_jph[jk] );
	     printf( "r, rmin, rmax = %.3e %.3e %.3e\n", r, rMin_rel_S, rMax_rel_S );
	     if (r < rMin_rel_S) rMin_rel_S = r;
	     if (r > rMax_rel_S) rMax_rel_S = r;
	   }
	 }
         uSum += u*E;
         uXSum += u*X*E;
         ESum += E;
         ETherm += (3.*c->prim[PPP])*c->cons[DEN]/c->prim[RHO];
         EXSum += E*X;
         dE += E;
         double u_eff = sqrt( gam*h*gam*h - 1.0 );
         if( u_eff > 1.0 ) EJet += E;
         XSum += X*E;
         MSum += M;
	 //if( X > 485.0 ) Ni += M;
         //Ni += X*M;
	 if( XP > Pth ) Ni_viaPPP += M;
	 Ni_direct += M*X;
         if( uMax < u ) uMax = u;
         if( gh_max < gam*h ) gh_max = gam*h;
      }
      I_th += dE*dE/dOmega;
   }

   MPI_Allreduce( MPI_IN_PLACE , &uMax   , 1 , MPI_DOUBLE , MPI_MAX , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &gh_max , 1 , MPI_DOUBLE , MPI_MAX , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &uSum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &uXSum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &ESum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &ETherm , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &EXSum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &EJet , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &XSum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &MSum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Ni_direct   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Ni_viaPPP   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &I_th   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );

   double uAv = uSum/ESum;
   double sintj2 = ESum/sqrt( 4.*M_PI*I_th );

   for( j=jmin ; j<jmax ; ++j ){
      int jk = j;
      for( i=0 ; i<Nr[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         double ur = c->prim[UU1];
         double up = c->prim[UU2];
         double u = sqrt(ur*ur+up*up);
         double r = c->riph;
         if( rMax < r && u>.5*uAv ) rMax = r;
         if( rMin > r && u>.5*uAv ) rMin = r;
         L1 += fabs(c->prim[PPP]/pow(c->prim[RHO],5./3.)-1.)*c->dr;
      }
   }

   MPI_Allreduce( MPI_IN_PLACE , &rMax , 1 , MPI_DOUBLE , MPI_MAX , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &rMin , 1 , MPI_DOUBLE , MPI_MIN , grid_comm );

   double v_phot = 0.0;
   if( rank==0 ){
      int jk = 0;
      double tau = 0.0;
      for( i=Nr[jk]-1 ; i>=0 && tau<1. ; --i ){
         struct cell * c = &(theCells[jk][i]);
         double rho = c->prim[RHO];
         double kappa = 3.*t*t;
         double dr = c->dr;
         //double r = c->riph;
         tau += rho*kappa*dr;
         double ur = c->prim[UU1];
         double up = c->prim[UU2];
         double gam = sqrt(1.+ur*ur+up*up);
         v_phot = ur/gam;
      }    
   }

   if( rank==0 ){
      FILE * rFile = fopen("report.dat","a");
      if ( theDomain->count_steps == 0 ){
	fprintf( rFile, "# time R_forward_shock Rmin uAverage uMax E_total M_total E_jet_rel gamma_h_max Rmin_ghRel_N Rmax_ghRel_N Rmin_ghRel_S Rmax_ghRel_S v_photosphere M_Ni56_direct M_Ni56_PPP sint_jet_2 E_therm_tot\n" ); }
      fprintf(rFile,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",t,rMax,rMin,uAv,uMax,ESum,MSum,EJet,gh_max,rMin_rel_N,rMax_rel_N,rMin_rel_S,rMax_rel_S,v_phot,Ni_direct, Ni_viaPPP,sintj2,ETherm);
      //fprintf(rFile,"%e %e %e %e %e\n",t,uSum,ESum,uXSum,EXSum);
      fclose(rFile);
   }

}
