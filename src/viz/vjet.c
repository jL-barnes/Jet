
//
// This code was created by Jeff Molofee '99 (ported to Linux/GLUT by Richard Campbell '99)
// If you've found this code useful, please let me know.
//
// Visit me at www.demonews.com/hosted/nehe 
// (email Richard Campbell at ulmont@bellsouth.net)
//

//#include <list>
#include <stdlib.h>
#include <math.h>
//#include <fstream>
#include <string.h>

#include <hdf5.h>

#include <GLUT/glut.h>    // Header File For The GLUT Library 
#include <OpenGL/gl.h>	// Header File For The OpenGL32 Library
#include <OpenGL/glu.h>	// Header File For The GLu32 Library
#include <unistd.h>     // needed to sleep

/* ASCII code for the escape key. */
#define ESCAPE 27

#define VAL_FLOOR 0 //(-HUGE_VAL) //0// -4 //-11 //-10 //(-HUGE_VAL)//-10// 0.0 //-20.0 //0.0 //-8.0 //(-HUGE_VAL)
#define VAL_CEIL  0.35 //HUGE_VAL //5e-11 //0 //-7 //0 //(HUGE_VAL)//1// 2.0 //2.0 //5e-4 //3.0 //(HUGE_VAL)
#define FIXMAXMIN true
#define COLORMAX 8
#define CAM_BACKUP  1.5

static int WindowWidth  = 970;
static int WindowHeight = 600;

int CommandMode;
int FullScreenMode=0;

int cmap = 4;
int cmap2 = 4;
int draw_bar = false;
int draw_t   = true;
int draw_lines = false;
int draw_jet = false;
int draw_scale = true;
int draw_detail = false;
int reflect  = true;
int invert = 0;
int valq=0;
int val2=-1;
int draw_border = 0;
int logscale = 0;
int floors=1;
int help_screen=0;
int print_vals=0;
int contours=0;

double rotate_angle = M_PI/2.;

double offx, offy, rescale, maxval, minval, maxval2, minval2;

double t;
int Nt,Np,Nq;
int * Nr;
double * t_jph;
double ** r_iph;
double *** theZones;

void get_rgb( double , float * , float * , float * , int , int );

double getval( double * thisZone , int q , double r , double theta ){
   if( q!=-1 ) return( thisZone[q] );
   double rho = thisZone[0];
   double P   = thisZone[1];
   double ur  = thisZone[2];
   double up  = thisZone[3];
   double s   = thisZone[4];
   double gam = sqrt(1.+ur*ur+up*up);
   double h = 1. + 4.*P/rho;
   double e = (rho+4.*P)*gam*gam-P - rho*gam;
   double v = sqrt( 1. - 1./gam/gam );

   double p = 2.5;
   double emm = rho*pow(P/rho,p-1.)*pow(P,(p+1.)/4.);
   double u = sqrt( ur*ur + up*up );

   double f = exp(-pow(r/cos(2.*theta)/1e3,20.));
   if( r>1e3 || cos(2.*theta) < 0.0 ) f=0.0;
   f *= pow( r , -4.5 )*3.;

   double dP = P - exp(s)*pow(rho,4./3.);

   if( dP < 1e-6 ) dP = 1e-6;

   double X = thisZone[5]*.5;
   if( s>0.5 ) X = 1.;

   double kappa = 3.*t*t;

   double shocked = 0.0;
   if( P/rho > 0.01 ) shocked = 1.0;
 
   double rs = r*sin(theta);
   //return( P*shocked + 1e-10 );
   //return( P/rho );
   //return( P/pow(rho,5./3.)-1. );
   //return( exp(-10.*pow(0.05*rho*kappa*r-1.,2.)) );// log(dP)/log(10.) );//sqrt(P)*rho );//e/gam/rho );// fabs(thisZone[1]/pow(thisZone[0],5./3.)-1.) );
   double T0 = 6.5e8;
   //double T = T0*pow(s,0.25);
   double T = T0*pow(P,0.25);
   double day = 86400.;
   T /= .5*day/1e3;
   double M0 = 0.04*2e33;
   double c = 3e10;
   double R = .5*day*c;
//   rho *= M0/R/R/R;
   //if( T>6e8 ) T=6e8;
   //return( T );
   //return(log(rho)/log(10.));
   double mix = s*(1.-s);
   return( mix*rho );

   double check = 0.0;
//   if( P/pow(rho,4./3.) > 1e-3 && rho*t*t*t > 1e-2 ) check = 1.0;
   if( ur > 10. ) check = 1.0;
   double tau = (4*M_PI*r*r)*rho*ur*gam*h;//rho*gam*(gam*h-1.);
//   if( ur < 0.001 ) tau = log(-1);
   double taumin = rho/1e8;
   if( tau < taumin ) tau = taumin;

   if( rho > 1e10 ) rho = 1e10;
   if( rho < 1e-14 ) rho = 1e-14;
//   return( log(tau/rho)/log(10.) );
//   return(log(rho)/log(10.));
//   return( tau );
//   return( ur/gam );
//   return( log(tau)/log(10.) );
//   return(log(r*r*rho)/log(10.));
   //return( log(rho*r*r*r) );

//   return( check );
//   return( fabs(P/pow(rho,4./3.)-1.)+1e-10 );
//   if( rho > .04 ) rho = .04;
//   return( rho );//log(rho)/log(10.) );

//   return( log(gam*h)/log(10.) );
 
//   return(log(P/pow(rho,4./3.))/log(10.));
   //return(u/sqrt(1.+u*u));
 
}

void getMaxMin(void){
   
   maxval = -HUGE_VAL;
   minval = HUGE_VAL;
   maxval2 = maxval;
   minval2 = minval;
   int q  = valq;
   int q2 = val2;
   double val,val2;
   for( int j=0 ; j<Nt ; ++j ){
      double tp = t_jph[j];
      double tm = t_jph[j-1];
      double th = .5*(tm+tp);
      for( int i=0 ; i<Nr[j]-1 ; ++i ){
         double rp = r_iph[j][i];
         double rm = r_iph[j][i-1];
         val  = getval( theZones[j][i], q , .5*(rp+rm) , th );
         val2 = getval( theZones[j][i], q2, .5*(rp+rm) , th );
         if( logscale ) val = log(val)/log(10.);
         if( maxval < val ) maxval = val;
         if( maxval2 < val2 ) maxval2 = val2;
         if( minval > val ) minval = val;
         if( minval2 > val2 ) minval2 = val2;
      }
   }
//   double midval = .5*(maxval+minval);
//   minval = midval-1.5;
//   maxval = midval+1.5;
//   minval = getval( theZones[Nt-1][Nr[Nt-1]-1] , q , r_iph[0][Nr[0]-1] );
//   if( logscale ) minval = log(minval)/log(10.);
//   maxval = minval+3.;
   if( floors ){
      if( maxval > VAL_CEIL  ) maxval = VAL_CEIL;
      if( minval < VAL_FLOOR ) minval = VAL_FLOOR;
   }
   if( FIXMAXMIN && floors ){
      maxval = VAL_CEIL;
      minval = VAL_FLOOR;
   }
//   maxval2 = 40;

//   minval2 = 0.0;
//   maxval2 = 3.0;
//   maxval = log(3e4/t/t/t)/log(10.);

//   minval = maxval-30;
//   minval2 = maxval2-30;

//   maxval2 = maxval;
//   minval2 = minval;
   //if( maxval2 > 2.0 ) maxval2 = 2.0;
   //maxval2 = 2.3;
   //minval2 = -5.0;
   //if( floors ) minval = maxval - 15.0;
   //minval = maxval - 12.0;
   //minval = maxval - 10.0;
//   maxval2 = 1.0;
//   minval2 = 0.0;
   //if( floors ) minval = maxval-5.0;
   //if( floors ) minval = log( getval( theZones[0][Nr[0]-1] , q ) )/log(10.);
/*
if( floors ){
  maxval  = 5;
  minval  = -15;

  maxval2 = 5;
  minval2 = -15;
}*/


}

void getH5dims( char * file , char * group , char * dset , hsize_t * dims ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );
   hid_t h5spc = H5Dget_space( h5dst );

   H5Sget_simple_extent_dims( h5spc , dims , NULL);

   H5Sclose( h5spc );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readSimple( char * file , char * group , char * dset , void * data , hid_t type ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readPatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   hsize_t mdims[dim];
   hsize_t fdims[dim];

   hsize_t fstart[dim];
   hsize_t fstride[dim];
   hsize_t fcount[dim];
   hsize_t fblock[dim];

   int d;
   for( d=0 ; d<dim ; ++d ){
      mdims[d] = loc_size[d];
      fdims[d] = glo_size[d];

      fstart[d]  = start[d];
      fstride[d] = 1;
      fcount[d]  = loc_size[d];
      fblock[d]  = 1;
   }
   hid_t mspace = H5Screate_simple(dim,mdims,NULL);
   hid_t fspace = H5Screate_simple(dim,fdims,NULL);

   H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

   H5Dread( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

/* The number of our GLUT window */
int window; 

// Here are the fonts: 
void ** glutFonts[7] = { 
    GLUT_BITMAP_9_BY_15, 
    GLUT_BITMAP_8_BY_13, 
    GLUT_BITMAP_TIMES_ROMAN_10, 
    GLUT_BITMAP_TIMES_ROMAN_24, 
    GLUT_BITMAP_HELVETICA_10, 
    GLUT_BITMAP_HELVETICA_12, 
    GLUT_BITMAP_HELVETICA_18 
}; 

// This function prints some text wherever you want it. 
void glutPrint(float x, float y, float z , void ** font, char* text, float r, float g, float b, float a) 
{ 
    if(!text || !strlen(text)) return; 
    bool blending = false; 
    if(glIsEnabled(GL_BLEND)) blending = true; 
    glEnable(GL_BLEND); 
    glColor4f(r,g,b,a); 
    glRasterPos3f(x,y,z); 
    while (*text) { 
        glutBitmapCharacter(font, *text); 
        text++; 
    } 
    if(!blending) glDisable(GL_BLEND); 
}  
 
// A general OpenGL initialization function.  Sets all of the initial parameters. 
void InitGL(int Width, int Height)	        // We call this right after our OpenGL window is created.
{
  glClearColor(0.8f, 0.8f, 0.8f, 0.0f);		// This Will Clear The Background Color To Black
  glClearDepth(1.0);				// Enables Clearing Of The Depth Buffer
  glDepthFunc(GL_LESS);			        // The Type Of Depth Test To Do
  glEnable(GL_DEPTH_TEST);		        // Enables Depth Testing
  glShadeModel(GL_SMOOTH);			// Enables Smooth Color Shading
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();				// Reset The Projection Matrix
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);	// Calculate The Aspect Ratio Of The Window
  glMatrixMode(GL_MODELVIEW);
}

/* The function called when our window is resized (which shouldn't happen, because we're fullscreen) */
void ReSizeGLScene(int Width, int Height)
{
  if (Height==0)				// Prevent A Divide By Zero If The Window Is Too Small
    Height=1;

  glViewport(0, 0, Width, Height);		// Reset The Current Viewport And Perspective Transformation

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
  glMatrixMode(GL_MODELVIEW);
}

void TakeScreenshot(const char *useless){
   int dimx = WindowWidth;
   int dimy = WindowHeight;
   printf("Writing test.ppm... ");
   float pixels[3*dimx*dimy];

   glReadBuffer(GL_BACK);
   glPixelStorei(GL_PACK_ALIGNMENT,1); // byte aligned output 
   glReadPixels(720-dimx/2, 450-dimy/2, dimx, dimy, GL_RGB, GL_FLOAT, pixels);

   FILE * pFile = fopen("test.ppm","w");
   fprintf(pFile,"P3\n%d %d\n255\n",dimx,dimy);

   int i,j;
   for( i=dimy-1 ; i>=0 ; i--){
      for( j=0 ; j<dimx ; j++ ){
         int pixelPos = (i*dimx+j)*3;
         int r = (int)(pixels[pixelPos+0]*255.0);
         int g = (int)(pixels[pixelPos+1]*255.0);
         int b = (int)(pixels[pixelPos+2]*255.0);
         fprintf(pFile,"%d %d %d\n",r,g,b);
      }
   }
   fclose(pFile);
   printf("done!\n");
}

/* The main drawing function. */
void DrawGLScene()
{
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
  glLoadIdentity();

   int q  = valq;
   int q2 = val2;

   double uMax = 0.0;

   getMaxMin();
   if( print_vals==0 ){
      printf("Max = %e Min = %e\n",maxval,minval);
      print_vals=1;
   }

   double camdist = -CAM_BACKUP;
   double xoff    = offx;
   double yoff    = offy;
   double zoff    = 0.0;

   int Num_Draws = 1+reflect+draw_border;
   for( int count = 0 ; count < Num_Draws ; ++count ){
      int draw_border_now = 0;
      int draw_reflected_now = 0;
      if( count==reflect+1 ) draw_border_now = 1;
      if( reflect && count==1 ) draw_reflected_now = 1;
      if( draw_border_now ) zoff -= .0001;
   for( int j=0 ; j<Nt ; ++j ){
      double tm = t_jph[j-1];
      double tp = t_jph[j];
      double th = .5*(tp+tm);
      double dt = tp-tm;
      
      if( !draw_border ){
         tm = th-.55*dt;
         tp = th+.55*dt;
         dt *= 1.1;
      }
      if( draw_reflected_now ){
         tm = -tm;//M_PI-tm;
         tp = -tp;//M_PI-tp;
         th = -th;//M_PI-th;
         dt = -dt;
      }
      tm += rotate_angle;
      tp += rotate_angle;
      th += rotate_angle;
      for( int i=1 ; i<Nr[j] ; ++i ){

         double rp = r_iph[j][i]/rescale;
         double rm = r_iph[j][i-1]/rescale;
         if( !draw_border ){
            double r = .5*(rp+rm);
            double dr = rp-rm;
            rp = r+.51*dr;
            rm = r-.51*dr;
         }

         double v0 = getval(theZones[j][i],q,.5*(rp+rm)*rescale, th-rotate_angle );
         double val = (v0-minval)/(maxval-minval);
         if(logscale) val = (log(v0)/log(10.)-minval)/(maxval-minval);
         if(draw_reflected_now) val = (getval(theZones[j][i],q2,.5*(rp+rm)*rescale, th-rotate_angle )-minval2)/(maxval2-minval2);
     //    if( draw_reflected_now ) val = (log(v0)/log(10.)-minval)/(maxval-minval);
         if( val > 1.0 ) val = 1.0;
         if( val < 0.0 ) val = 0.0;
         if( contours ){
            val = (double)( ((int)(16.*val)) )/16.;
         }
         double u = getval( theZones[j][i] , 2 , .5*(rp+rm)*rescale , th-rotate_angle );
         if( uMax < u ) uMax = u;

         float rrr,ggg,bbb;
         get_rgb( val , &rrr , &ggg , &bbb , cmap , invert );
         if(draw_reflected_now) get_rgb( val , &rrr , &ggg , &bbb , cmap2 , invert );

         if(draw_detail){
            double saturate = pow(sin(25.*val),4.);
            rrr *= 1.-saturate;
            ggg *= 1.-saturate;
            bbb *= 1.-saturate;
         }

         if( !draw_border_now ){ 
            glColor3f( rrr , ggg , bbb );
            glBegin(GL_POLYGON);
         }else{
            //glLineWidth(2.0f);
            glLineWidth(3.0f);
            glColor3f(0.0,0.0,0.0);
            glBegin(GL_LINE_LOOP);
         }
    
         double c0 = rm*cos(tm);
         double c1 = rm*sin(tm);
         glVertex3f( c0-xoff, c1-yoff, camdist-zoff );

         c0 = rp*cos(tm);
         c1 = rp*sin(tm);
         glVertex3f( c0-xoff, c1-yoff, camdist-zoff );

         c0 = rp*cos(th)/cos(.5*dt);
         c1 = rp*sin(th)/cos(.5*dt);
         glVertex3f( c0-xoff, c1-yoff, camdist-zoff );

         c0 = rp*cos(tp);
         c1 = rp*sin(tp);
         glVertex3f( c0-xoff, c1-yoff, camdist-zoff );

         c0 = rm*cos(tp);
         c1 = rm*sin(tp);
         glVertex3f( c0-xoff, c1-yoff, camdist-zoff );

         glEnd();
      }
      }

   }


   if( draw_bar ){
      double xb = 0.6;
      double hb = 1.0;
      double wb = 0.02;
      int Nb = 1000;
      double dy = hb/(double)Nb;
      for( int k=0 ; k<Nb ; ++k ){
         double y = (double)k*dy - .5*hb;
         double val = (double)k/(double)Nb;
         float rrr,ggg,bbb;
         get_rgb( val , &rrr , &ggg , &bbb , cmap , invert );
         glLineWidth(0.0f);
         glColor3f( rrr , ggg , bbb );
         glBegin(GL_POLYGON);
         glVertex3f( xb    , y    , camdist + .001 ); 
         glVertex3f( xb+wb , y    , camdist + .001 ); 
         glVertex3f( xb+wb , y+dy , camdist + .001 ); 
         glVertex3f( xb    , y+dy , camdist + .001 ); 
         glEnd();
      }
      int Nv = 8;
      for( int k=0 ; k<Nv ; ++k ){
         double y = (double)k*hb/(double)(Nv-1) - .5*hb;
         double val = (double)k/(double)(Nv-1)*(maxval-minval) + minval;
         char valname[256];
         sprintf(valname,"%+.2e",val);
         glLineWidth(1.0f);
         glColor3f(0.0,0.0,0.0);
         glBegin(GL_LINE_LOOP);
         glVertex3f( xb    , y , camdist + .0011 );
         glVertex3f( xb+wb , y , camdist + .0011 );
         glEnd();
         glutPrint( xb+1.5*wb , y-.007 , camdist + .001 ,glutFonts[6] , valname , 0.0f, 0.0f , 0.0f , 0.5f );
      }
   } 
   
   if( draw_scale ){
      double thmax = t_jph[Nt-1] + rotate_angle;
      double thmin = rotate_angle;
      double r_min = r_iph[0][0]/rescale;
      double r_max = r_iph[0][Nr[0]-1]/rescale;

      glLineWidth(2.0f);
      glColor3f( 0.0 , 0.0 , 0.0 );
      glBegin(GL_LINES);
      glVertex3f( r_min*cos(thmax) - xoff , r_min*sin(thmax) - yoff , camdist + .001 );
      glVertex3f( r_max*cos(thmax) - xoff , r_max*sin(thmax) - yoff , camdist + .001 );
      glEnd();

      r_min = 0.0;
      r_max = r_iph[0][Nr[0]-1];
      int logmax =  (int) ( log(r_max)/log(10.) );
      if( r_max > pow(10.,logmax) ) ++logmax;
      r_max = pow( 10. , logmax );
      r_max /= rescale;

      int Nv = 10;
      int kmin = (int) ( (r_iph[0][0]/rescale-r_min)/(r_max-r_min)*(double)Nv );
      int kmax = (int) ( (r_iph[0][Nr[0]-1]/rescale-r_min)/(r_max-r_min)*(double)Nv );
      if( kmax > 10 ) kmax = 10;
      for( int k=kmin ; k<kmax ; ++k ){
         double dr = r_max-r_min;
         double r = (double)(k+1)*dr/(double)Nv + r_min;
         double dt = 0.02*1.0/r;
         char valname[256];
         sprintf(valname,"%.0e",r*rescale);
         glLineWidth(2.0f);
         glColor3f(0.0,0.0,0.0);
         glBegin(GL_LINES);
         glVertex3f( r*cos(thmax)    - xoff , r*sin(thmax)    - yoff , camdist + .001 );
         glVertex3f( r*cos(thmax+dt) - xoff , r*sin(thmax+dt) - yoff , camdist + .001 );
         glEnd();
         glutPrint( r*cos(thmax+4.*dt) - xoff , r*sin(thmax+4.*dt) - yoff , camdist + .001 , glutFonts[6] , valname , 0.0f, 0.0f , 0.0f , 0.5f );
            int Nt = 20;
            for( int l=0 ; l<Nt ; ++l ){
               double t0 = ((double)l    )*(thmax-thmin)/(double)Nt + thmin;
               double t1 = ((double)l+0.5)*(thmax-thmin)/(double)Nt + thmin;
      	       glLineWidth(2.0f);
               glBegin(GL_LINES);
               glVertex3f( r*cos(t0) - xoff , r*sin(t0) - yoff , camdist + .001 );
               glVertex3f( r*cos(t1) - xoff , r*sin(t1) - yoff , camdist + .001 );
               glEnd();
            }
      }  
      if( kmin == 0 ){
         //Draw Minors
         r_max = (r_max-r_min)/(double)Nv;
         r_min = 0.0;
         kmin = (int) ( (r_iph[0][0]/rescale-r_min)/(r_max-r_min)*(double)Nv );
         for( int k=kmin ; k<Nv ; ++k ){
            double dr = r_max-r_min;
            double r = (double)(k+1)*dr/(double)Nv + r_min;
            double dt = 0.01*1.0/r;
            glLineWidth(1.0f);
            glColor3f(0.0,0.0,0.0);
            glBegin(GL_LINES);
            glVertex3f( r*cos(thmax)    - xoff , r*sin(thmax)    - yoff , camdist + .001 );
            glVertex3f( r*cos(thmax+dt) - xoff , r*sin(thmax+dt) - yoff , camdist + .001 );
            glEnd();
         } 
      } 
   }
 
   if( draw_t ){
      char tprint[256];
      sprintf(tprint,"t = %.2e",t);
      glutPrint( -.6 , .35 , camdist + .001 , glutFonts[6] , tprint , 0.0f, 0.0f , 0.0f , 0.5f );
      sprintf(tprint,"gMax = %.1f",sqrt(1.+uMax*uMax));
      glutPrint( -.6 , .25 , camdist + .001 , glutFonts[6] , tprint , 0.0f, 0.0f , 0.0f , 0.5f );
   }
   if( draw_lines ){
      glLineWidth( 2.0f );
      glColor3f( 1.0 , 1.0 , 1.0 );
      double opening_angle = 0.14/sqrt(2.);
      double th0 = opening_angle+rotate_angle;
      double r_min = r_iph[0][0];
      double r_max = r_iph[0][Nr[0]-1];
      int Nt = 100;
      glBegin(GL_LINE_LOOP);
      for( int k=0 ; k<Nt ; ++k ){
         double r = ( (double)k*(r_max-r_min)/(double)Nt + r_min )/rescale;
         double theta_new = ( (double)k*(M_PI)/(double)Nt ) + rotate_angle;
         r = (3.*r_min-r_min*theta_new*theta_new )/rescale;
         glVertex3f( r*cos(theta_new)-xoff , r*sin(theta_new)-yoff , camdist      + .001 );
         if( k%2 == 1 ){
            glEnd();
            glBegin(GL_LINE_LOOP);
         }
      }
      glEnd();
      th0 = -opening_angle+rotate_angle;
      glBegin(GL_LINE_LOOP);
      for( int k=0 ; k<Nt ; ++k ){
         double r = ( (double)k*(r_max-r_min)/(double)Nt + r_min )/rescale;
         glVertex3f( r*cos(th0)-xoff , r*sin(th0)-yoff , camdist      + .001 );
         if( k%2 == 1 ){
            glEnd();
            glBegin(GL_LINE_LOOP);
         }
      }
      glEnd();
   }
   if( draw_jet ){
      glLineWidth( 4.0f );
      glColor3f( 0.0 , 0.0 , 0.0 );
      double th0 = 0.15; 
      double r0  = 2.*r_iph[0][0];//2.*1e-3*t;
      double K = 1.0;

      int NR = 30; 
      glBegin(GL_LINE_LOOP);

      for( int k=0 ; k<NR ; ++k ){
/*
         double w = 1.3;
         double th = -(double)k*M_PI/(double)NR/2.;
         double r = (0.9+0.09*t)*sqrt(pow(w,-2./3.)*cos(th)*cos(th) + pow(w,4./3.)*sin(th)*sin(th));
         r /= rescale;
         glVertex3f( -r*sin(th)-xoff , r*cos(th)-yoff , camdist      + .001 );
         if( k%2 == 1 ){ 
            glEnd();
            glBegin(GL_LINE_LOOP);
         }    
*/

         double r = 2.*(double)k*r0/(double)NR;
         double fr = log(r/r0)-.5*(r*r/r0/r0);
         double cost = 1.-th0*th0*(K+fr);
         if( cost >  1. ) cost =  1.;
         if( cost < -1. ) cost = -1.;
         double th = acos( cost );
         r /= rescale;
         glVertex3f( -r*sin(th)-xoff , r*cos(th)-yoff , camdist      + .001 );
         if( k%2 == 1 ){ 
            glEnd();
            glBegin(GL_LINE_LOOP);
         }    

      }    
      glEnd();
   }
   if( help_screen ){
      int NLines = 9;
      double dy = -.04;
      char help[NLines][256];
      sprintf(help[0],"Help Display:");
      sprintf(help[1],"b - Toggle Colorbar");
      sprintf(help[2],"c - Change Colormap");
      sprintf(help[3],"f - Toggle Max/Min Floors");
      sprintf(help[4],"p - Toggle Planet Data");
      sprintf(help[5],"1-9 - Choose Primitive Variable to Display");
      sprintf(help[6],"wasd - Move Camera");
      sprintf(help[7],"z/x - Zoom in/out");
      sprintf(help[8],"h - Toggle Help Screen");
      for( int i=0 ; i<NLines ; ++i ){
         glutPrint( -.8 , i*dy , camdist + .001 , glutFonts[6] , help[i] , 0.0f, 0.0f , 0.0f , 0.5f );
      }
   } 

   if( CommandMode ){
      TakeScreenshot("test.ppm");
      exit(1);
   }

   glutSwapBuffers();

}

/* The function called whenever a key is pressed. */
void keyPressed(unsigned char key, int x, int y) 
{
   usleep(100);
   if (key == ESCAPE){
      glutDestroyWindow(window); 
      exit(0);                   
   }
   if (key == '[' ){
       --cmap;
       cmap = (cmap+COLORMAX) % COLORMAX;
   }  
   if (key == ']' ){
       ++cmap;
       cmap = cmap % COLORMAX;
   }
   if( key >= (int)'0' && key < (int)'1'+Nq ) valq = (int)key-(int)'1';
   if( key == 'b' ) draw_bar = !draw_bar;
   if( key == 'f' ) floors = !floors;
   if( key == 'g' ) draw_border = !draw_border;
   if( key == 'h' ) help_screen = !help_screen;
   if( key == 'i' ) invert = !invert;
   if( key == 'l' ) logscale = !logscale;
   if( key == 'r' ) reflect = !reflect;
   if( key == 's' ) draw_scale = !draw_scale;
   if( key == 't' ) draw_t = !draw_t;
   if( key == '/' ) draw_detail = !draw_detail;
   if( key == 'x' ){
      rescale *= 1.3;
      offx /= 1.3;
      offy /= 1.3;
   }
   if( key == 'z' ){
      rescale /= 1.3;
      offx *= 1.3;
      offy *= 1.3;
   }
   if (key == 'F') {
      if (FullScreenMode) {
         glutReshapeWindow(WindowWidth, WindowHeight);
         glutPositionWindow(0,0);
         FullScreenMode = 0; 
      }else{
         glutFullScreen();
         FullScreenMode = 1; 
      }
   }
   if( key == 'L' ) draw_lines = !draw_lines;
   if( key == 'c' ) contours = !contours;
   glutPostRedisplay();
}

void specialKeyPressed(int key, int x, int y){
   usleep(100);
   if (key == GLUT_KEY_LEFT ) offx -= .1;
   if (key == GLUT_KEY_RIGHT) offx += .1;

   if (key == GLUT_KEY_UP   ) offy += .1;
   if (key == GLUT_KEY_DOWN ) offy -= .1;
}

int main(int argc, char **argv) 
{
   if( argc < 2 ){
      printf("Please specify the input file.\n");
      exit(1);
   }
   char filename[256];
   if( argv[1] ){
      strcpy( filename , argv[1] );
   }
   CommandMode=0;
   if( argc>2 ){ CommandMode=1; FullScreenMode=1; }
   char group1[256];
   char group2[256];
   strcpy( group1 , "Grid" );
   strcpy( group2 , "Data" );

   hsize_t dims[3];
   
   readSimple( filename , group1 , (char *)"T" , &t , H5T_NATIVE_DOUBLE );
   getH5dims( filename , group1 , (char *)"t_jph" , dims );
   Nt = dims[0]-1;
   getH5dims( filename , group1 , (char *)"p_kph" , dims );
   Np = dims[0]-1;

   Nr = (int *) malloc( Nt*sizeof(int) );
   t_jph = (double *) malloc( (Nt+1)*sizeof(double) );
   int Tindex[Nt];
   
   printf("t = %.2f, Nt = %d\n",t,Nt);

   readSimple( filename , group1 , (char *)"t_jph" , t_jph , H5T_NATIVE_DOUBLE );

   int start[2]    = {Nt*(Np/2),0};
   int loc_size[2] = {Nt       ,1};
   int glo_size[2] = {Nt*Np    ,1};

   readPatch( filename , group1 , (char *)"Nr"    , Nr     , H5T_NATIVE_INT , 2 , start , loc_size , glo_size);
   readPatch( filename , group1 , (char *)"Index" , Tindex , H5T_NATIVE_INT , 2 , start , loc_size , glo_size);

   getH5dims( filename , group2 , (char *)"Cells" , dims );
   int Nc = dims[0];
   Nq = dims[1]-1;
   printf("Nc = %d Nt = %d Nq=%d\n",Nc,Nt,Nq);
   printf("tm = %e, tp = %e\n",t_jph[0],t_jph[Nt]);

   theZones = (double ***) malloc( Nt*sizeof(double **) );
   for( int j=0 ; j<Nt ; ++j ){
      theZones[j] = (double **) malloc( Nr[j]*sizeof( double * ) );
      for( int i=0 ; i<Nr[j] ; ++i ){
         theZones[j][i] = (double *) malloc( Nq*sizeof( double ) );
      }
   }
   r_iph = (double **) malloc( Nt*sizeof( double * ) );
   for( int j=0 ; j<Nt ; ++j ){
      r_iph[j] = (double *) malloc( Nr[j]*sizeof( double ) );
   }

   printf("Zones Allocated\n");
   loc_size[1] = Nq+1;
   glo_size[0] = Nc;
   glo_size[1] = Nq+1;

   for( int j=0 ; j<Nt ; ++j ){
      loc_size[0] = Nr[j];
      start[0] = Tindex[j];
      double TrackData[Nr[j]*(Nq+1)];
      readPatch( filename , group2 , (char *)"Cells" , TrackData , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size);
      for( int i=0 ; i<Nr[j] ; ++i ){
         r_iph[j][i] = TrackData[i*(Nq+1) + Nq];
         for( int q=0 ; q<Nq ; ++q ){
            theZones[j][i][q] = TrackData[i*(Nq+1) + q];
         }
      }
   }

   t_jph++;

   double RR = r_iph[0][Nr[0]-1];
   double r_max = .98*RR;
   double r_min = r_iph[0][0];
   //double r_min = .82*RR;

   double t_bo = 1.05;
   double r_shock = (asinh(cosh(t/t_bo)) - asinh(1.)) + .025;
   r_shock = 0.2*t + 1.5;
   if( r_shock < 1.2*t-10. ) r_shock = 1.2*t-10.;
   //if( t<.25 ) draw_border = 1;
   //r_shock = 1.0;

   printf("Rmin = %.2e Rmax = %.2e\n",r_min,r_max);
   //double thalf = .5*t_jph[Nt-1];
//   rescale = 0.25*(r_max-r_min); //.3
//   rescale = 0.95*(r_max-r_min);
   rescale = 2*(r_max-r_min);

   offy = 0.4;
//   rescale *= 0.5;
//   rescale = 0.2*t;

//   rescale = 0.95*t;
//   offx = -1.0;

/*
   rescale = 0.4*r_max;
   offx = 0.0;//.5*(r_min+r_max)/rescale;
   offy = 0.15*r_max/rescale;//.5*r_max/rescale;//.5*(r_min+r_max)*sin(thalf)/rescale;
*/

   //rescale = r_max;
//   rescale = 0.1*r_max;
//   offx = -0.5;
//   offx = -0.8;
   //offx = -0.1/rescale;
//   offy = 0.5; //0.4*r_max/rescale;
   //offy = 1.0/rescale;
//   offy = 8.0;

//   offx = .9*(.25*r_min+.75*r_max)/rescale;
//   offy = -.2*offx;
   offx = 0.0;
   offy = 0.0;
   printf("rescale = %e\n",rescale);

//////////////////////////////
  glutInit(&argc, argv);  
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);  
  glutInitWindowSize(WindowWidth, WindowHeight);
  glutInitWindowPosition(0, 0);
  window = glutCreateWindow("Flying Grid");
  glutDisplayFunc(&DrawGLScene);  
  if(FullScreenMode) glutFullScreen();
  glutIdleFunc(&DrawGLScene);
  glutReshapeFunc(&ReSizeGLScene);
  glutKeyboardFunc(&keyPressed);
  glutSpecialFunc(&specialKeyPressed);
  InitGL(WindowWidth, WindowHeight);
  glutMainLoop();  

   for( int j=0 ; j<Nt ; ++j ){
      for( int i=0 ; i<Nr[j] ; ++i ){
         free( theZones[j][i] );
      }
      free( theZones[j] );
   }
   free( theZones );
   for( int j=0 ; j<Nt ; ++j ){
      free( r_iph[j] );
   }
   free( r_iph );
   free( Nr );
   t_jph--;
   free( t_jph );

  return (0);
}
