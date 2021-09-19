#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <chrono>
#include <sys/time.h>

#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TVector3.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>

using namespace std;

const double pi=4.0*atan(1.0); /* Gives back value of Pi */
const double spedc=299792458.0; /* Speed of Light in m/s */

/* Set the value of the asymptotic parameter of the refractive index model */
const double A_ice=1.78;
const double TransitionBoundary=0;
// const double A_ice=1.775;
// const double TransitionBoundary=14.9;

/* Get the value of the B parameter for the refractive index model */
double GetB(double z){
  z=fabs(z);
  double B=0;

  B=-0.43;
  
  // if(z<=TransitionBoundary){
  //   B=-0.5019;
  // }else{
  //   B=-0.448023;
  // }
  
  return B;
}

/* Get the value of the C parameter for the refractive index model */
double GetC(double z){
  z=fabs(z);
  double C=0;

  C=0.0132;
 
  // if(z<=TransitionBoundary){
  //   C=0.03247;
  // }else{
  //   C=0.02469;
  // }
 
  return C;
}

/* Get the value of refractive index model for a given depth  */
double Getnz(double z){
  z=fabs(z);
  return A_ice+GetB(z)*exp(-GetC(z)*z);
}

/* E-feild Power Fresnel coefficient for S-polarised wave which is perpendicular to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R */
double Refl_S(double thetai){

  double Nair=1;
  double Nice=Getnz(0); 
  double n1=Nice;
  double n2=Nair;
  
  double sqterm=sqrt(1-pow((n1/n2)*(sin(thetai)),2));
  double num=n1*cos(thetai)-n2*sqterm;
  double den=n1*cos(thetai)+n2*sqterm;
  double RS=(num*num)/(den*den);

  if(std::isnan(RS)){
    RS=1;
  }
  return (RS);
}

/* E-feild Power Fresnel coefficient for P-polarised wave which is parallel to the plane of propogation/incidence. This function gives you back the reflectance. The transmittance is T=1-R */
double Refl_P(double thetai){
   
  double Nair=1;
  double Nice=Getnz(0); 
  double n1=Nice;
  double n2=Nair;

  double sqterm=sqrt(1-pow((n1/n2)*(sin(thetai)),2));
  double num=n1*sqterm-n2*cos(thetai);
  double den=n1*sqterm+n2*cos(thetai);
  double RP=(num*num)/(den*den);
  if(std::isnan(RP)){
    RP=1;
  }
  return (RP);
}

/* The temperature and attenuation model has been taken from AraSim which also took it from here http://icecube.wisc.edu/~araproject/radio/ . This is basically Matt Newcomb's icecube directory which has alot of information, plots and codes about South Pole Ice activities. Please read it if you find it interesting. */

/* Temperature model:The model takes in value of depth z in m and returns the value of temperature in Celsius.*/
double GetIceTemperature(double z){
  double depth=fabs(z);
  double t = 1.83415e-09*pow(depth,3) + (-1.59061e-08*pow(depth,2)) + 0.00267687*depth + (-51.0696 );
  return t;
}

/* Ice Attenuation Length model: Takes in value of frequency in Ghz and depth z and returns you the value of attenuation length in m */
double GetIceAttenuationLength(double z, double frequency){

  double t =GetIceTemperature(z);
  const double f0=0.0001, f2=3.16;
  const double w0=log(f0), w1=0.0, w2=log(f2), w=log(frequency);
  const double b0=-6.74890+t*(0.026709-t*0.000884);
  const double b1=-6.22121-t*(0.070927+t*0.001773);
  const double b2=-4.09468-t*(0.002213+t*0.000332);
  double a,bb;
  if(frequency<1.){
    a=(b1*w0-b0*w1)/(w0-w1);
    bb=(b1-b0)/(w1-w0);
  }
  else{
    a=(b2*w1-b1*w2)/(w1-w2);
    bb=(b2-b1)/(w2-w1);
  }
  double Lval=1./exp(a+bb*w);
  return Lval;
}

/* Setup the integrand to calculate the attenuation */
double AttenuationIntegrand (double x, void * params) {

  double *p=(double*)params;

  double A0=p[0];
  double Frequency=p[1];
  double L=p[2];

  double Integrand=(A0/GetIceAttenuationLength(x,Frequency))*sqrt(1+pow(tan(asin(L/Getnz(x))) ,2));
  return Integrand;
}

/* Integrate over the integrand to calculate the attenuation */
double IntegrateOverLAttn (double A0, double Frequency, double z0, double z1, double Lvalue) {
  gsl_integration_workspace * w= gsl_integration_workspace_alloc (1000);

  double result, error;
  double zlimit[2] = {z0,z1};
  double param[3] = {A0,Frequency,Lvalue};
  
  gsl_function F;
  F.function = &AttenuationIntegrand;
  F.params = &param;

  gsl_integration_qags (&F, zlimit[0], zlimit[1], 0, 1e-7, 1000,
                        w, &result, &error);

  // printf ("result          = % .18f\n", result);
  // printf ("estimated error = % .18f\n", error);
  // printf ("intervals       = %zu\n", w->size);

  gsl_integration_workspace_free (w);  

  return fabs(result);
}

/* Calculate the total attenuation for each type of ray */
double GetTotalAttenuationDirect (double A0, double frequency, double z0, double z1, double Lvalue) {
  z0=fabs(z0);
  z1=fabs(z1);
  return IntegrateOverLAttn(A0,frequency,z0,z1,Lvalue);
}

double GetTotalAttenuationReflected (double A0, double frequency, double z0, double z1, double Lvalue) {
  z0=fabs(z0);
  z1=fabs(z1);
  return IntegrateOverLAttn(A0,frequency,z0,0.000001,Lvalue) + IntegrateOverLAttn(A0,frequency,z1,0.000001,Lvalue);
}

double GetTotalAttenuationRefracted (double A0, double frequency, double z0, double z1, double zmax, double Lvalue) {
  z0=fabs(z0);
  z1=fabs(z1);
  return IntegrateOverLAttn(A0,frequency,z0,zmax,Lvalue) + IntegrateOverLAttn(A0,frequency,z1,zmax,Lvalue);
}

/* Use GSL minimiser which relies on calculating function deriavtives. This function uses GSL's Newton's algorithm to find root for a given function. */
double FindFunctionRootFDF(gsl_function_fdf FDF,double x_lo, double x_hi){
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  double x0 ,x = (x_lo+x_hi)/2;
  
  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, x);

  // printf ("using %s method\n",
  //         gsl_root_fdfsolver_name (s));

  // printf ("%-5s %10s %10s %10s\n",
  //         "iter", "root", "err", "err(est)");
  do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, 1e-6);
      
      // if (status == GSL_SUCCESS)
      //   printf ("Converged:\n");
      
      // printf ("%5d %10.7f %10.7f\n",
      //         iter, x, x - x0);
      
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fdfsolver_free (s);
  
  return x;
}

/* Use GSL minimiser which uses GSL's false position algorithm to find root for a given function. */
double FindFunctionRoot(gsl_function F,double x_lo, double x_hi)
{
  int status;
  int iter = 0, max_iter = 200;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
 
  T = gsl_root_fsolver_falsepos;
  s = gsl_root_fsolver_alloc (T);
  gsl_set_error_handler_off();
  status = gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      
      if(x_lo>x_hi){
      	swap(x_lo,x_hi);
      }
      double checkval=(*((F).function))(r,(F).params);
      status = gsl_root_test_residual (checkval,1e-6);
      //status = gsl_root_test_interval (x_lo, x_hi,0, 0.000001);
      
      if (status == GSL_SUCCESS){
	// printf ("Converged:");
	// printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi,r,x_hi - x_lo);
      }
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  //printf ("status = %s\n", gsl_strerror (status));
  gsl_root_fsolver_free (s);

  return r;
}

/* Use GSL minimiser which uses GSL's false position algorithm to find root for a given function. */
double FindFunctionRootZmax(gsl_function F,double x_lo, double x_hi)
{
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
 
  T = gsl_root_fsolver_falsepos;
  s = gsl_root_fsolver_alloc (T);
  gsl_set_error_handler_off();
  status = gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,1e-6, 1e-6);
	
      if (status == GSL_SUCCESS){
	// printf ("Converged:");
	// printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi,r,x_hi - x_lo);
      }
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  gsl_root_fsolver_free (s);

  return r;
}


/* Define the function that will be minimized to get the value of the depth of the turning point for a given refracted ray. This function basically requires the value of the L parameter to find the depth of the turning point.  This comes from the part of the fDnfR function where sqrt( n(z) - L ). This imposes the constraint then n(z)=L at the turning point of the ray from which we can find zmax. */
struct Minnz_params { double a,l; };
double GetMinnz(double x,void *params){
  struct Minnz_params *p= (struct Minnz_params *) params;
  double A = p->a;
  double L = p->l;
  return A+GetB(x)*exp(-GetC(x)*x)-L;
}

/* Get the value of the depth of the turning point for the refracted ray */
double GetZmax(double A, double L){
  gsl_function F1;
  struct Minnz_params params1= {A,L};
  F1.function = &GetMinnz;
  F1.params = &params1;
  double zmax=FindFunctionRootZmax(F1,0.0,5000);
  return zmax;
}

/* Analytical solution describing ray paths in ice as function of depth */
struct fDnfR_params { double a, b, c, l; };
double fDnfR(double x,void *params){
  
  struct fDnfR_params *p= (struct fDnfR_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double L = p->l;
  
  return (L/C)*(1.0/sqrt(A*A-L*L))*(C*x-log(A*Getnz(x)-L*L+sqrt(A*A-L*L)*sqrt(pow(Getnz(x),2)-L*L)));;
}

/* Analytical solution describing the ray path in ice as a function of the L parameter */
struct fDnfR_L_params { double a, b, c, z; };
double fDnfR_L(double x,void *params){
  
  struct fDnfR_L_params *p= (struct fDnfR_L_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Z = p->z;
  
  double result=(x/C)*(1.0/sqrt(A*A-x*x))*(C*Z-log(A*Getnz(Z)-x*x+sqrt(A*A-x*x)*sqrt(pow(Getnz(Z),2)-x*x)));
  
  return result;
}

/* The function used to calculate ray propogation time in ice */
struct ftimeD_params { double a, b, c, speedc,l; };
double ftimeD(double x,void *params){

  struct ftimeD_params *p= (struct ftimeD_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Speedc = p->speedc;
  double L = p->l;
  
  return (1.0/(Speedc*C*sqrt(pow(Getnz(x),2)-L*L)))*(pow(Getnz(x),2)-L*L+(C*x-log(A*Getnz(x)-L*L+sqrt(A*A-L*L)*sqrt(pow(Getnz(x),2)-L*L)))*(A*A*sqrt(pow(Getnz(x),2)-L*L))/sqrt(A*A-L*L) +A*sqrt(pow(Getnz(x),2)-L*L)*log(Getnz(x)+sqrt(pow(Getnz(x),2)-L*L)) );
}

/* The set of functions starting with the name "fDa" are used in the minimisation procedure to find the launch angle (or the L parameter) for the direct ray */
struct fDanfRa_params { double a, z0, x1, z1; };
double fDa(double x,void *params){
  struct fDanfRa_params *p= (struct fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  struct fDnfR_L_params params1a = {A, GetB(z1), GetC(z1), z1};
  struct fDnfR_L_params params1b = {A, GetB(z0), GetC(z0), z0};
  struct fDnfR_L_params params1c = {A, GetB(TransitionBoundary), GetC(TransitionBoundary), -TransitionBoundary};
  struct fDnfR_L_params params1d = {A, GetB(-(TransitionBoundary+0.000001)), GetC(-(TransitionBoundary+0.000001)), -(TransitionBoundary+0.000001)};

  double distancez0z1=0;
  if(TransitionBoundary!=0){
    if (fabs(z0)>TransitionBoundary && fabs(z1)>TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
    }
    if (fabs(z0)>TransitionBoundary && fabs(z1)<TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1c) + fDnfR_L(x,&params1d) - fDnfR_L(x,&params1b);
    }
    if (fabs(z0)<TransitionBoundary && fabs(z1)<TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
    }
    if (fabs(z0)==TransitionBoundary && fabs(z1)<TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
    }
    if (fabs(z0)==TransitionBoundary && fabs(z1)==TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
    }
    if (fabs(z0)>TransitionBoundary && fabs(z1)==TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1c) + fDnfR_L(x,&params1d) - fDnfR_L(x,&params1b);
    }
  }else{
    distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
  }
  // if(isnan(distancez0z1)){
  //   distancez0z1=1e9;
  // }
  double output=distancez0z1-x1;

  return output;
}

double fDa_df(double x,void *params){
  gsl_function F;
  F.function = &fDa;
  F.params = params;
 
  double result,abserr;
  gsl_deriv_central (&F, x, 1e-8, &result, &abserr);

  return result;
}

void fDa_fdf (double x, void *params,double *y, double *dy){ 
  *y = fDa(x,params);
  *dy = fDa_df(x,params);
}

/* The set of functions starting with the name "fRa" are used in the minimisation procedure to find the launch angle (or the L parameter) for the reflected ray */
double fRa(double x,void *params){
  struct fDanfRa_params *p= (struct fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  struct fDnfR_L_params params1a = {A, GetB(z1), -GetC(z1), -z1};
  struct fDnfR_L_params params1b = {A, GetB(z0), -GetC(z0), -z0};
  struct fDnfR_L_params params1c = {A, GetB(1e-7), -GetC(1e-7), 1e-7};
  struct fDnfR_L_params params1d = {A, GetB(TransitionBoundary), -GetC(TransitionBoundary), TransitionBoundary};
  struct fDnfR_L_params params1f = {A, GetB(TransitionBoundary+0.000001), -GetC(TransitionBoundary+0.000001), TransitionBoundary+0.000001};

  double distancez0z1=0;
  double distancez0surface=0;
  if(TransitionBoundary!=0){
    if (fabs(z0)>TransitionBoundary && fabs(z1)>TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
      distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1d) + fDnfR_L(x,&params1f) - fDnfR_L(x,&params1b);
    }
    if (fabs(z0)>TransitionBoundary && fabs(z1)<TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1d) + fDnfR_L(x,&params1f) - fDnfR_L(x,&params1b);
      distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1d) + fDnfR_L(x,&params1f) - fDnfR_L(x,&params1b);
    }
    if (fabs(z0)<TransitionBoundary && fabs(z1)<TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
      distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b); 
    }
    if (fabs(z0)==TransitionBoundary && fabs(z1)<TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
      distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b); 
    }
    if (fabs(z0)==TransitionBoundary && fabs(z1)==TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
      distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b); 
    }
    if (fabs(z0)>TransitionBoundary && fabs(z1)==TransitionBoundary){
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1d) + fDnfR_L(x,&params1f) - fDnfR_L(x,&params1b);
      distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1d) + fDnfR_L(x,&params1f) - fDnfR_L(x,&params1b); 
    }
  }else{
    distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
    distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b);
  }
  // if(isnan(distancez0z1)){
  //   distancez0z1=1e9;
  // }
  // if(isnan(distancez0surface)){
  //   distancez0surface=1e9;
  // }
  double output= distancez0z1 - 2*(distancez0surface) - x1;
  //cout<<distancez0z1<<" "<<2*(distancez0surface)<<" "<<x1<<" "<<output<<endl;
  
  return output;
}

double fRa_df(double x,void *params){
  gsl_function F;
  F.function = &fRa;
  F.params = params;
 
  double result,abserr;
  gsl_deriv_central (&F, x, 1e-8, &result, &abserr);

  return result;
}

void fRa_fdf (double x, void *params,double *y, double *dy){ 
  *y = fRa(x,params);
  *dy = fRa_df(x,params);
}

/* This function is minimised to find the launch angle (or the L parameter) for the refracted ray */
double fRaa(double x,void *params){
  struct fDanfRa_params *p= (struct fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;
  
  double zmax= GetZmax(A,x)+1e-7;
  double output=0;
  if(zmax>0){  
  
    struct fDnfR_L_params params1a = {A, GetB(z1), -GetC(z1), -z1};
    struct fDnfR_L_params params1b = {A, GetB(z0), -GetC(z0), -z0};
    struct fDnfR_L_params params1c = {A, GetB(zmax), -GetC(zmax), zmax};
    struct fDnfR_L_params params1d = {A, GetB(TransitionBoundary), -GetC(TransitionBoundary), TransitionBoundary};
    struct fDnfR_L_params params1f = {A, GetB(TransitionBoundary+1e-7), -GetC(TransitionBoundary+1e-7), TransitionBoundary+1e-7};

    double distancez0z1=0;
    double distancez0surface=0;
    if(TransitionBoundary!=0){
      if (fabs(z0)>TransitionBoundary && fabs(z1)>TransitionBoundary){
	if(zmax<=TransitionBoundary){
	  distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
	  distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1d) + fDnfR_L(x,&params1f) - fDnfR_L(x,&params1b);
	}else{
	  distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
	  distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b);
	}
      }
      if (fabs(z0)>TransitionBoundary && fabs(z1)<TransitionBoundary){
	distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1d) + fDnfR_L(x,&params1f) - fDnfR_L(x,&params1b);
	distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1d) + fDnfR_L(x,&params1f) - fDnfR_L(x,&params1b);
      }
      if (fabs(z0)<TransitionBoundary && fabs(z1)<TransitionBoundary){
	distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
	distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b); 
      }
      if (fabs(z0)==TransitionBoundary && fabs(z1)<TransitionBoundary){
	distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
	distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b); 
      }
      if (fabs(z0)==TransitionBoundary && fabs(z1)==TransitionBoundary){
	distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
	distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b); 
      }
      if (fabs(z0)>TransitionBoundary && fabs(z1)==TransitionBoundary){
	distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1d) + fDnfR_L(x,&params1f) - fDnfR_L(x,&params1b);
	distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1d) + fDnfR_L(x,&params1f) - fDnfR_L(x,&params1b); 
      }
    }else{
      distancez0z1=fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b);
      distancez0surface=fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b);
    }

    if(isnan(distancez0z1)){
      distancez0z1=1e9;
    }
    if(isnan(distancez0surface)){
      distancez0surface=1e9;
    }
    output= (distancez0z1 - 2*(distancez0surface) - x1);
  }else{
    output=1e9;
  }
  
  if(output<-1e8){
    output=1e9;
  }
  
  return output;
}

double fRaa_df(double x,void *params){
  gsl_function F;
  F.function = &fRaa;
  F.params = params;
 
  double result,abserr;
  gsl_deriv_central (&F, x, 1e-8, &result, &abserr);

  return result;
}

void fRaa_fdf (double x, void *params,double *y, double *dy){ 
  *y = fRaa(x,params);
  *dy = fRaa_df(x,params);
}

/* This functions works for the Direct ray and gives you back the launch angle, receive angle and propagation time of the ray together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
double* GetDirectRayPar(double z0, double x1, double z1){

  double *output=new double[5];
  
  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* First we setup the fDa function that will be minimised to get the launch angle (or the L parameter) for the direct ray. */
  gsl_function F1;
  struct fDanfRa_params params1= {A_ice, z0, x1, z1};
  F1.function = &fDa;
  F1.params = &params1;

  // gsl_function_fdf F1;
  // struct fDanfRa_params params1= {A_ice, z0, x1, z1};
  // F1.f = &fDa;
  // F1.df = &fDa_df;
  // F1.fdf = &fDa_fdf;
  // F1.params = &params1;
  
  /* In my raytracing solution given in the function fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0 and z1 and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */ 
  double UpLimnz[]={Getnz(z1),Getnz(z0)};
  double* UpperLimitL=min_element(UpLimnz,UpLimnz+2);

  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function. */
  double lvalueD=FindFunctionRoot(F1,1e-7,UpperLimitL[0]);
  double LangD=asin(lvalueD/Getnz(z0))*(180.0/pi);
  double checkzeroD=fDa(lvalueD,&params1);

  /* Get the propagation time for the direct ray using the ftimeD function after we have gotten the value of the L parameter. */
  struct ftimeD_params params2a = {A_ice, GetB(z0), -GetC(z0), spedc,lvalueD};
  struct ftimeD_params params2b = {A_ice, GetB(z1), -GetC(z1), spedc,lvalueD};
  struct ftimeD_params params2c = {A_ice, GetB(TransitionBoundary), -GetC(TransitionBoundary), spedc, lvalueD};
  struct ftimeD_params params2d = {A_ice, GetB(TransitionBoundary+1e-7), -GetC(TransitionBoundary+1e-7), spedc, lvalueD};

  /* we do the subtraction because we are measuring the time taken between the Tx and Rx positions */
  double timeD=0;
  if(TransitionBoundary!=0){
    if (fabs(z0)>TransitionBoundary && fabs(z1)>TransitionBoundary){
      timeD=ftimeD(-z0,&params2a) - ftimeD(-z1,&params2b);
    }
    if (fabs(z0)>TransitionBoundary && fabs(z1)<TransitionBoundary){
      timeD=ftimeD(-z0,&params2a) - ftimeD(TransitionBoundary+1e-7,&params2d) + ftimeD(TransitionBoundary,&params2c) - ftimeD(-z1,&params2b);
    }
    if (fabs(z0)<TransitionBoundary && fabs(z1)<TransitionBoundary){
      timeD=ftimeD(-z0,&params2a) - ftimeD(-z1,&params2b);
    }
    if (fabs(z0)==TransitionBoundary && fabs(z1)<TransitionBoundary){
      timeD=ftimeD(-z0,&params2a) - ftimeD(-z1,&params2b);
    }
    if (fabs(z0)==TransitionBoundary && fabs(z1)==TransitionBoundary){
      timeD=ftimeD(-z0,&params2a) - ftimeD(-z1,&params2b);
    }
    if (fabs(z0)>TransitionBoundary && fabs(z1)==TransitionBoundary){
      timeD=ftimeD(-z0,&params2a) - ftimeD(TransitionBoundary+1e-7,&params2d) + ftimeD(TransitionBoundary,&params2c) - ftimeD(-z1,&params2b);
    }
  }else{
    timeD=ftimeD(-z0,&params2a) - ftimeD(-z1,&params2b);
  }

  /* Setup the function that will be used to calculate the angle of reception for all the rays */
  gsl_function F5;
  struct fDnfR_params params5a = {A_ice, GetB(z1), -GetC(z1), lvalueD};
  double result, abserr;
  F5.function = &fDnfR;

  /* Calculate the recieve angle for direc rays by calculating the derivative of the function at the Rx position */
  F5.params = &params5a;
  gsl_deriv_central (&F5, -z1, 1e-8, &result, &abserr);
  double RangD=atan(result)*(180.0/pi);

  /* When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line. */
  if(z1==z0 && isnan(RangD)==true){
    RangD=180-LangD;
  }
  
  /* This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation. */
  if(z1!=z0 && isnan(RangD)==true){
    RangD=90;
  }
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangD;
  output[1]=LangD;
  output[2]=timeD;
  output[3]=lvalueD;
  output[4]=checkzeroD;

  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangD;
    output[1]=180-RangD;
  }
  
  return output;
}

/* This functions works for the Reflected ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. */
double *GetReflectedRayPar(double z0, double x1 ,double z1){

  double *output=new double[8];

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* First we setup the fRa function that will be minimised to get the launch angle (or the L parameter) for the reflected ray. */
  gsl_function F3;
  struct fDanfRa_params params3= {A_ice, z0, x1, z1};
  F3.function = &fRa;
  F3.params = &params3;

  // gsl_function_fdf F3;
  // struct fDanfRa_params params3= {A_ice, z0, x1, z1};
  // F3.f = &fRa;
  // F3.df = &fRa_df;
  // F3.fdf = &fRa_fdf;
  // F3.params = &params3;
  
  /* In my raytracing solution given in the function fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0, z1 and also 1e-7 m and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */
  double UpLimnz[]={Getnz(z1),Getnz(z0),Getnz(1e-7)};
  double *UpperLimitL=min_element(UpLimnz,UpLimnz+3);
  
  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fRa function. */
  double lvalueR=FindFunctionRoot(F3,1e-7,UpperLimitL[0]);
  double LangR=asin(lvalueR/Getnz(z0))*(180.0/pi);
  double checkzeroR=fRa(lvalueR,&params3); 

  /* Get the propagation time for the reflected ray using the ftimeD function after we have gotten the value of the L parameter. */
  struct ftimeD_params params3a = {A_ice, GetB(z0), GetC(z0), spedc,lvalueR};
  struct ftimeD_params params3b = {A_ice, GetB(z1), GetC(z1), spedc,lvalueR};
  struct ftimeD_params params3c = {A_ice, GetB(1e-7), GetC(1e-7), spedc,lvalueR};
  struct ftimeD_params params3d = {A_ice, GetB(TransitionBoundary), GetC(TransitionBoundary), spedc, lvalueR};
  struct ftimeD_params params3f = {A_ice, GetB(TransitionBoundary+1e-7), GetC(TransitionBoundary+1e-7), spedc, lvalueR};
  /* We do the subtraction because we are measuring the time taken between the Tx and Rx positions. In the reflected case we basically have two direct rays 1) from Tx to surface 2) from surface to Rx. Also get the time for the two individual direct rays separately */
  double timeR1=0;
  double timeR2=0;
  if(TransitionBoundary!=0){
    if (fabs(z0)>TransitionBoundary && fabs(z1)>TransitionBoundary){
      timeR1=ftimeD(-1e-7,&params3c) - ftimeD(-TransitionBoundary,&params3d) + ftimeD(-(TransitionBoundary+1e-7),&params3f) - ftimeD(z0,&params3a);
      timeR2=ftimeD(-1e-7,&params3c) - ftimeD(-TransitionBoundary,&params3d) + ftimeD(-(TransitionBoundary+1e-7),&params3f) - ftimeD(z1,&params3b);
    }
    if (fabs(z0)>TransitionBoundary && fabs(z1)<TransitionBoundary){
      timeR1=ftimeD(-1e-7,&params3c) - ftimeD(-TransitionBoundary,&params3d) + ftimeD(-(TransitionBoundary+1e-7),&params3f) - ftimeD(z0,&params3a);
      timeR2=ftimeD(-1e-7,&params3c) - ftimeD(z1,&params3b);
    }
    if (fabs(z0)<TransitionBoundary && fabs(z1)<TransitionBoundary){
      timeR1=ftimeD(-1e-7,&params3c) - ftimeD(z0,&params3a);
      timeR2=ftimeD(-1e-7,&params3c) - ftimeD(z1,&params3b); 
    }
    if (fabs(z0)==TransitionBoundary && fabs(z1)<TransitionBoundary){
      timeR1=ftimeD(-1e-7,&params3c) - ftimeD(z0,&params3a);
      timeR2=ftimeD(-1e-7,&params3c) - ftimeD(z1,&params3b);
    }
    if (fabs(z0)==TransitionBoundary && fabs(z1)==TransitionBoundary){
      timeR1=ftimeD(-1e-7,&params3c) - ftimeD(z0,&params3a);
      timeR2=ftimeD(-1e-7,&params3c) - ftimeD(z1,&params3b); 
    }
    if (fabs(z0)>TransitionBoundary && fabs(z1)==TransitionBoundary){
      timeR1=ftimeD(-1e-7,&params3c) - ftimeD(-TransitionBoundary,&params3d) + ftimeD(-(TransitionBoundary+1e-7),&params3f) - ftimeD(z0,&params3a);
      timeR2=ftimeD(-1e-7,&params3c) - ftimeD(z1,&params3b); 
    } 
  }else{
    timeR1=ftimeD(-1e-7,&params3c) - ftimeD(z0,&params3a);
    timeR2=ftimeD(-1e-7,&params3c) - ftimeD(z1,&params3b);
  }
  double timeR=timeR1 + timeR2;
  
  /* flip the times back if the original positions were flipped */
  if(Flip==true){
    double dumR=timeR2;
    timeR2=timeR1;
    timeR1=dumR;
  }
  timeR1=timeR1;
  timeR2=timeR2;

  /* Setup the function that will be used to calculate the angle of reception for all the rays */
  gsl_function F5;
  struct fDnfR_params params5b = {A_ice, GetB(z1), GetC(z1), lvalueR};
  double result, abserr;
  F5.function = &fDnfR;
  
  /* Calculate the recieve angle for reflected ray by calculating the derivative of the function at the Rx position */
  F5.params = &params5b;
  gsl_deriv_central (&F5, z1, 1e-8, &result, &abserr);
  double RangR=180-atan(result)*(180.0/pi);

  /* When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line. */
  if(z1==z0 && isnan(RangR)==true){
    RangR=180-LangR;
  }

  /* This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation. */
  if(z1!=z0 && isnan(RangR)==true){
    RangR=90;
  }

  if(std::isnan(checkzeroR)){
    checkzeroR=-1000;
  }
  
  /* Calculate the angle of incidence of the reflected ray at the surface ice. This will be used to calculate the Fresnel Coefficients. The angle is calculated by calculating the derivative of the ray path fucnction at the surface*/
  struct fDnfR_params paramsIAngB = {A_ice, GetB(1e-7), GetC(1e-7), lvalueR};
  F5.function = &fDnfR; 
  F5.params = &paramsIAngB;
  gsl_deriv_central (&F5, -1e-7, 1e-8, &result, &abserr);
  double IncidenceAngleInIce=atan(result)*(180.0/pi);
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangR;
  output[1]=LangR;
  output[2]=timeR;
  output[3]=lvalueR;
  output[4]=checkzeroR;
  output[5]=timeR1;
  output[6]=timeR2;
  output[7]=IncidenceAngleInIce;
  
  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangR;
    output[1]=180-RangR;
  } 
  
  return output;
}

/* This functions works for the Refracted ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter and checkzero variable. checkzero variable checks how close the minimiser came to 0. 0 is perfect and less than 0.5 is pretty good. more than that should not be acceptable. It requires the launch angle of the reflected ray as an input. */
double *GetRefractedRayPar(double z0, double x1 ,double z1, double LangR, double RangR){
  
  double *output=new double[8*2];

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }

  dsw=0;
  if(Flip==true){
    dsw=180-LangR;
    LangR=180-RangR;
    RangR=dsw;
  }
  
  /* Set up all the variables that will be used to get the parameters for the refracted ray */
  //double lvalueR=sin(LangR*(pi/180))*Getnz(z0);
  double lvalueRa[2]={0,0};
  double LangRa[2]={0,0};
  double checkzeroRa[2]={-1000,-1000};

  double timeRa[2]={0,0};
  double timeRa1[2]={0,0};
  double timeRa2[2]={0,0};
  double raytime[2]={0,0};
  double RangRa[2]={0,0};
  double zmax[2]={10,10};
  
  /* In my raytracing solution given in the function fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0, z1 and also 1e-7 m and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */
  double UpLimnz[]={Getnz(z0),Getnz(z1)};
  double *UpperLimitL=min_element(UpLimnz,UpLimnz+2);
  
  /* First we setup the fRa function that will be minimised to get the launch angle (or the L parameter) for the refracted ray. */
  gsl_function F4;
  struct fDanfRa_params params4= {A_ice, z0, x1, z1};
  F4.function = &fRaa;
  F4.params = &params4;
  //checkzeroRa[1]=fRaa(lvalueRa[1],&params4);
 
  gsl_function_fdf F4b;
  F4b.f = &fRaa;
  F4b.df = &fRaa_df;
  F4b.fdf = &fRaa_fdf;
  F4b.params = &params4;
  
  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fRaa function. The thing to note here is the lower limit of the minimisation function is set to the L value corresponding to the reflected ray launch angle. Since we know the refracted ray always has bigger launch angle the reflected ray this reduces our range and makes the function more efficient at finding the refracted ray launch angle. */
 double LowerLimit=0;
  LowerLimit=Getnz(z0)*sin((64.0*(pi/180.0)));
  if(LowerLimit>UpperLimitL[0]){
    LowerLimit=Getnz(z0)*sin((LangR*(pi/180.0)));
  }

  lvalueRa[0]=FindFunctionRoot(F4,LowerLimit,UpperLimitL[0]);
  LangRa[0]=asin(lvalueRa[0]/Getnz(z0))*(180.0/pi);
  checkzeroRa[0]=(fRaa(lvalueRa[0],&params4));
  zmax[0]=GetZmax(A_ice,lvalueRa[0])+1e-7;

  if(fabs(checkzeroRa[0])>0.5){
    lvalueRa[0]=FindFunctionRootFDF(F4b,LowerLimit,UpperLimitL[0]);
    LangRa[0]=asin(lvalueRa[0]/Getnz(z0))*(180.0/pi);
    checkzeroRa[0]=fRaa(lvalueRa[0],&params4);
    zmax[0]=GetZmax(A_ice,lvalueRa[0])+1e-7;
  }
  
  if(fabs(checkzeroRa[0])<0.5){

    lvalueRa[1]=FindFunctionRoot(F4,lvalueRa[0]-0.23,lvalueRa[0]-0.023);
    LangRa[1]=asin(lvalueRa[1]/Getnz(z0))*(180.0/pi);
    checkzeroRa[1]=fRaa(lvalueRa[1],&params4);
    zmax[1]=GetZmax(A_ice,lvalueRa[1])+1e-7;
  
    if(fabs(checkzeroRa[1])>0.5 || std::isnan(checkzeroRa[1])==true || fabs(lvalueRa[1]-lvalueRa[0])<1e-4 ){
      lvalueRa[1]=FindFunctionRoot(F4,lvalueRa[0]-0.15,lvalueRa[0]-0.023);
      LangRa[1]=asin(lvalueRa[1]/Getnz(z0))*(180.0/pi);
      checkzeroRa[1]=fRaa(lvalueRa[1],&params4);
      zmax[1]=GetZmax(A_ice,lvalueRa[1])+1e-7;
    }

    if(fabs(checkzeroRa[1])>0.5 || std::isnan(checkzeroRa[1])==true || fabs(lvalueRa[1]-lvalueRa[0])<1e-4 ){
      lvalueRa[1]=FindFunctionRoot(F4,lvalueRa[0]+0.005,UpperLimitL[0]);
      LangRa[1]=asin(lvalueRa[1]/Getnz(z0))*(180.0/pi);
      checkzeroRa[1]=fRaa(lvalueRa[1],&params4);
      zmax[1]=GetZmax(A_ice,lvalueRa[1])+1e-7;
    }
    
    if(fabs(checkzeroRa[1])>0.5 || std::isnan(checkzeroRa[1])==true || fabs(lvalueRa[1]-lvalueRa[0])<1e-4 ){  
      lvalueRa[1]=FindFunctionRootFDF(F4b,lvalueRa[0]-0.23,lvalueRa[0]-0.023);
      LangRa[1]=asin(lvalueRa[1]/Getnz(z0))*(180.0/pi);
      checkzeroRa[1]=fRaa(lvalueRa[1],&params4);
      zmax[1]=GetZmax(A_ice,lvalueRa[1])+1e-7;
    }

    if(fabs(checkzeroRa[1])>0.5 || std::isnan(checkzeroRa[1])==true || fabs(lvalueRa[1]-lvalueRa[0])<1e-4 ){  
      lvalueRa[1]=FindFunctionRootFDF(F4b,lvalueRa[0]-0.1,lvalueRa[0]-0.023);
      LangRa[1]=asin(lvalueRa[1]/Getnz(z0))*(180.0/pi);
      checkzeroRa[1]=fRaa(lvalueRa[1],&params4);
      zmax[1]=GetZmax(A_ice,lvalueRa[1])+1e-7;
    }
    
    if(LangRa[1]<LangRa[0] && fabs(checkzeroRa[0])<0.5 && fabs(checkzeroRa[1])<0.5){
      swap(lvalueRa[1],lvalueRa[0]);
      swap(LangRa[1],LangRa[0]);
      swap(checkzeroRa[1],checkzeroRa[0]);
      swap(zmax[1],zmax[0]);
    }
    
  }else{
    lvalueRa[1]=0;
    LangRa[1]=0;
    checkzeroRa[1]=-1000;
    zmax[1]=-1000;
  }
 
  
  for(int i=0;i<2;i++){////loop over the two refracted solutions
  
    /* If we still did not find a refracted ray then set the check zero parameter to 1000 to make sure my code does not output this as a possible solution */
    if(isnan(checkzeroRa[i])==true){
      checkzeroRa[i]=-1000;
    }

    /* If the turning point depth also came out to be zero then now we are sure that there is no refracted ray */
    if(zmax[i]==1e-7 || zmax[i]<=0){
      checkzeroRa[i]=-1000;
    }

    /* Set parameters for ftimeD function to get the propagation time for the refracted ray */
    struct ftimeD_params params4a = {A_ice, GetB(z0), GetC(z0), spedc,lvalueRa[i]};
    struct ftimeD_params params4b = {A_ice, GetB(z1), GetC(z1), spedc,lvalueRa[i]};
    struct ftimeD_params params4c = {A_ice, GetB(zmax[i]), GetC(zmax[i]), spedc,lvalueRa[i]};
    struct ftimeD_params params4d = {A_ice, GetB(TransitionBoundary), GetC(TransitionBoundary), spedc, lvalueRa[i]};
    struct ftimeD_params params4f = {A_ice, GetB(TransitionBoundary+1e-7), GetC(TransitionBoundary+1e-7), spedc, lvalueRa[i]};
    
    /* This if condition checks if the function has not gone crazy and given us a turning point of the ray which is lower than both Tx and Rx and is shallower in depth than both */
    if((z0<-zmax[i] || zmax[i]<-z1)){
      /* We do the subtraction because we are measuring the time taken between the Tx and Rx positions. In the refracted case we basically have two direct rays 1) from Tx to turning point 2) from turning point to Rx. Also get the time for the two individual direct rays separately */
     
      if(TransitionBoundary!=0){
	if (fabs(z0)>TransitionBoundary && fabs(z1)>TransitionBoundary){
	  if(zmax[i]<=TransitionBoundary){
	    timeRa1[i]=ftimeD(-zmax[i],&params4c) - ftimeD(-TransitionBoundary,&params4d) + ftimeD(-(TransitionBoundary+1e-7),&params4f) - ftimeD(z0,&params4a);
	    timeRa2[i]=ftimeD(-zmax[i],&params4c) - ftimeD(-TransitionBoundary,&params4d) + ftimeD(-(TransitionBoundary+1e-7),&params4f) - ftimeD(z1,&params4b);	
	  }else{
	    timeRa1[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z0,&params4a);
	    timeRa2[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z1,&params4b);
	  }  
	}
	if (fabs(z0)>TransitionBoundary && fabs(z1)<TransitionBoundary){
	  timeRa1[i]=ftimeD(-zmax[i],&params4c) - ftimeD(-TransitionBoundary,&params4d) + ftimeD(-(TransitionBoundary+1e-7),&params4f) - ftimeD(z0,&params4a);
	  timeRa2[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z1,&params4b);
	}
	if (fabs(z0)<TransitionBoundary && fabs(z1)<TransitionBoundary){
	  timeRa1[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z0,&params4a);
	  timeRa2[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z1,&params4b);
	} 
	if (fabs(z0)==TransitionBoundary && fabs(z1)<TransitionBoundary){
	  timeRa1[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z0,&params4a);
	  timeRa2[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z1,&params4b); 
	}
	if (fabs(z0)==TransitionBoundary && fabs(z1)==TransitionBoundary){
	  timeRa1[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z0,&params4a);
	  timeRa2[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z1,&params4b); 
	}
	if (fabs(z0)>TransitionBoundary && fabs(z1)==TransitionBoundary){
	  timeRa1[i]=ftimeD(-zmax[i],&params4c) - ftimeD(-TransitionBoundary,&params4d) + ftimeD(-(TransitionBoundary+1e-7),&params4f) - ftimeD(z0,&params4a);
	  timeRa2[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z1,&params4b);
	}
      }else{
	timeRa1[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z0,&params4a);
	timeRa2[i]=ftimeD(-zmax[i],&params4c) - ftimeD(z1,&params4b);
      }
      
      raytime[i]= timeRa1[i] + timeRa2[i];
      
      if(Flip==true){
	double dumRa=timeRa2[i];
	timeRa2[i]=timeRa1[i];
	timeRa1[i]=dumRa;
      }
    }
    timeRa[i]=raytime[i];

    /* Setup the function that will be used to calculate the angle of reception for all the rays */
    gsl_function F5;
    struct fDnfR_params params5c = {A_ice, GetB(z1), GetC(z1), lvalueRa[i]};
    double result, abserr;
    F5.function = &fDnfR;
    
    /* Calculate the recieve angle for refacted ray by calculating the derivative of the function at the Rx position */
    F5.params = &params5c;
    gsl_deriv_central (&F5, z1, 1e-8, &result, &abserr);
    RangRa[i]=180-atan(result)*(180.0/pi);

    /* When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line. */
    if(z1==z0 && isnan(RangRa[i])==true){
      RangRa[i]=180-LangRa[i];
    }

    /* This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation. */
    if(z1!=z0 && isnan(RangRa[i])==true){
      RangRa[i]=90;
    }

  }//// i loop
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangRa[0];
  output[1]=LangRa[0];
  output[2]=timeRa[0];
  output[3]=lvalueRa[0];
  output[4]=checkzeroRa[0];
  output[5]=timeRa1[0];
  output[6]=timeRa2[0];
  output[7]=zmax[0];
  
  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangRa[0];
    output[1]=180-RangRa[0];
  }

  output[8]=RangRa[1];
  output[9]=LangRa[1];
  output[10]=timeRa[1];
  output[11]=lvalueRa[1];
  output[12]=checkzeroRa[1];
  output[13]=timeRa1[1];
  output[14]=timeRa2[1];
  output[15]=zmax[1];
  
  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[8]=180-LangRa[1];
    output[9]=180-RangRa[1];
  }
  
  return output;
  }

/* This function returns the x and z values for the full Direct ray path in a TGraph and also prints out the ray path in a text file */
TGraph* GetFullDirectRayPath(double z0, double x1, double z1,double lvalueD){

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
   
  /* Set the name of the text files */
  ofstream aoutD("DirectRay.txt");
  /* Set the step size for plotting */
  double h=0.1;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;

  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;
  struct fDnfR_params params6a;
  struct fDnfR_params params6b;
  struct fDnfR_params params6c;
  struct fDnfR_params params6d;
  
  TGraph *gr1=new TGraph();
  for(int i=0;i<dmax;i++){
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueD};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueD};
    params6c = {A_ice, GetB(TransitionBoundary), GetC(TransitionBoundary), lvalueD};
    params6d = {A_ice, GetB(TransitionBoundary+1e-7), GetC(TransitionBoundary+1e-7), lvalueD};

    if(TransitionBoundary!=0){
      if (fabs(z0)<TransitionBoundary && fabs(zn)<TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)<TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(-TransitionBoundary,&params6c) + fDnfR(-(TransitionBoundary+1e-7),&params6d) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)>TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)==TransitionBoundary && fabs(zn)<TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)==TransitionBoundary && fabs(zn)==TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)==TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(-TransitionBoundary,&params6c) + fDnfR(-(TransitionBoundary+1e-7),&params6d) - fDnfR(z0,&params6b);
      }
    }else{
      xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
    }
    
    checknan=fDnfR(zn,&params6a);
    if(isnan(checknan)==false && Flip==false){
      gr1->SetPoint(npnt,xn,zn);
      aoutD<<npnt<<" "<<xn<<" "<<zn<<endl;;
      npnt++;
    }

    if(isnan(checknan)==false && Flip==true){
      gr1->SetPoint(npnt,x1-xn,zn);
      aoutD<<npnt<<" "<<x1-xn<<" "<<zn<<endl;;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;      
    }  
  }

  params6a = {A_ice, GetB(zn), GetC(zn), lvalueD};
  params6b = {A_ice, GetB(z0), GetC(z0), lvalueD};
  xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);  
  if(Flip==true){
    gr1->SetPoint(npnt,x1-xn,z0);
    aoutD<<npnt<<" "<<x1-xn<<" "<<zn<<endl;;
  }else{
    gr1->SetPoint(npnt,xn,z0);
    aoutD<<npnt<<" "<<xn<<" "<<zn<<endl;;
  }
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  return gr1;

}

/* This function returns the x and z values for the full Reflected ray path in a TGraph and also prints out the ray path in a text file */
TGraph* GetFullReflectedRayPath(double z0, double x1, double z1,double lvalueR){

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* Set the name of the text files */
  ofstream aoutR("ReflectedRay.txt");
  /* Set the step size for plotting. */
  double h=0.1;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;
  
  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;  
  struct fDnfR_params params6a;
  struct fDnfR_params params6b;
  struct fDnfR_params params6c;
  struct fDnfR_params params6f;
  struct fDnfR_params params6d;

  /* Map out the 1st part of the reflected ray */
  TGraph *gr2=new TGraph();
  for(int i=0;i<dmax;i++){
    params6a = {A_ice, GetB(zn), -GetC(zn), lvalueR};
    params6b = {A_ice, GetB(z0), -GetC(z0), lvalueR};
    params6c = {A_ice, GetB(1e-7), -GetC(1e-7), lvalueR};
    params6d = {A_ice, GetB(TransitionBoundary), -GetC(TransitionBoundary), lvalueR};
    params6f = {A_ice, GetB(TransitionBoundary+1e-7), -GetC(TransitionBoundary+1e-7), lvalueR};
    double distancez0z1=0;
    double distancez0surface=0;  

    if(TransitionBoundary!=0){
      if (fabs(z0)>TransitionBoundary && fabs(zn)>TransitionBoundary){
	distancez0z1=fDnfR(-zn,&params6a) - fDnfR(-z0,&params6b);
	distancez0surface=fDnfR(1e-7,&params6c) - fDnfR(TransitionBoundary,&params6d) + fDnfR(TransitionBoundary+1e-7,&params6f) - fDnfR(-z0,&params6b);
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)<TransitionBoundary){
	distancez0z1=fDnfR(-zn,&params6a) - fDnfR(TransitionBoundary,&params6d) + fDnfR(TransitionBoundary+1e-7,&params6f) - fDnfR(-z0,&params6b);
	distancez0surface=fDnfR(1e-7,&params6c) - fDnfR(TransitionBoundary,&params6d) + fDnfR(TransitionBoundary+1e-7,&params6f) - fDnfR(-z0,&params6b);
      }
      if (fabs(z0)<TransitionBoundary && fabs(zn)<TransitionBoundary){
	distancez0z1=fDnfR(-zn,&params6a) - fDnfR(-z0,&params6b);
	distancez0surface=fDnfR(1e-7,&params6c) - fDnfR(-z0,&params6b); 
      }
      if (fabs(z0)==TransitionBoundary && fabs(zn)<TransitionBoundary){
	distancez0z1=fDnfR(-zn,&params6a) - fDnfR(-z0,&params6b);
	distancez0surface=fDnfR(1e-7,&params6c) - fDnfR(-z0,&params6b); 
      }
      if (fabs(z0)==TransitionBoundary && fabs(zn)==TransitionBoundary){
	distancez0z1=fDnfR(-zn,&params6a) - fDnfR(-z0,&params6b);
	distancez0surface=fDnfR(1e-7,&params6c) - fDnfR(-z0,&params6b);
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)==TransitionBoundary){
	distancez0z1=fDnfR(-zn,&params6a) - fDnfR(TransitionBoundary,&params6d) + fDnfR(TransitionBoundary+1e-7,&params6f) - fDnfR(-z0,&params6b);
	distancez0surface=fDnfR(1e-7,&params6c) - fDnfR(TransitionBoundary,&params6d) + fDnfR(TransitionBoundary+1e-7,&params6f) - fDnfR(-z0,&params6b);
      }
    }else{
      distancez0z1=fDnfR(-zn,&params6a) - fDnfR(-z0,&params6b);
      distancez0surface=fDnfR(1e-7,&params6c) - fDnfR(-z0,&params6b);
    }

    xn= distancez0z1 - 2*(distancez0surface);
 
    checknan=fDnfR(-zn,&params6a);
    if(isnan(checknan)==false && zn<=0 && Flip==false){
      gr2->SetPoint(npnt,xn,zn);
      aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(isnan(checknan)==false && zn<=0 && Flip==true){
      gr2->SetPoint(npnt,x1-xn,zn);
      aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
      
    zn=zn+h;
    if(zn>0){
      i=dmax+2;      
    }
  }
  
  /* Map out the 2nd part of the reflected ray */
  zn=-1e-7;
  for(int i=0;i<dmax;i++){  
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueR};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueR};
    params6c = {A_ice, GetB(TransitionBoundary), GetC(TransitionBoundary), lvalueR};
    params6d = {A_ice, GetB(TransitionBoundary+1e-7), GetC(TransitionBoundary+1e-7), lvalueR};

    if(TransitionBoundary!=0){
      if (fabs(z0)<TransitionBoundary && fabs(zn)<TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)<TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(-TransitionBoundary,&params6c) + fDnfR(-(TransitionBoundary+1e-7),&params6d) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)>TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)==TransitionBoundary && fabs(zn)<TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)==TransitionBoundary && fabs(zn)==TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)==TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(-TransitionBoundary,&params6c) + fDnfR(-(TransitionBoundary+1e-6),&params6d) - fDnfR(z0,&params6b);
      }
    }else{
      xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
    }
    
    checknan=fDnfR(zn,&params6a);
    if(isnan(checknan)==false && Flip==false){
      gr2->SetPoint(npnt,xn,zn);
      aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }
      
    if(isnan(checknan)==false && Flip==true){
      gr2->SetPoint(npnt,x1-xn,zn);
      aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;
    }
  }

  params6a = {A_ice, GetB(zn), GetC(zn), lvalueR};
  params6b = {A_ice, GetB(z0), GetC(z0), lvalueR};
  xn=fDnfR(zn,&params6a) -fDnfR(z0,&params6b);
  if(Flip==true){
    gr2->SetPoint(npnt,x1-xn,zn);
    aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
  }else{
    gr2->SetPoint(npnt,xn,zn);
    aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
  }
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  return gr2;

}

/* This function returns the x and z values for the full Refracted ray path in a TGraph and also prints out the ray path in a text file */
TGraph* GetFullRefractedRayPath(double z0, double x1, double z1, double zmax, double lvalueRa, int raynumber){

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* Set the name of the text files */
  ofstream aoutRa;
  /* Set the name of the text files */
  if(raynumber==1){
    aoutRa.open("RefractedRay1.txt");
  }
  if(raynumber==2){
    aoutRa.open("RefractedRay2.txt");
  }
  /* Set the step size for plotting. */
  double h=0.1;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;

  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;
  struct fDnfR_params params6a;
  struct fDnfR_params params6b;
  struct fDnfR_params params6c;  
  struct fDnfR_params params6f;
  struct fDnfR_params params6d;

  /* Map out the 1st part of the refracted ray */
  TGraph *gr3=new TGraph();    
  for(int i=0;i<dmax;i++){
    params6a = {A_ice, GetB(zn), -GetC(zn), lvalueRa};
    params6b = {A_ice, GetB(z0), -GetC(z0), lvalueRa};
    params6c = {A_ice, GetB(zmax), -GetC(zmax), lvalueRa};
    params6d = {A_ice, GetB(TransitionBoundary), -GetC(TransitionBoundary), lvalueRa};
    params6f = {A_ice, GetB(TransitionBoundary+1e-6), -GetC(TransitionBoundary+1e-6), lvalueRa};

    double distancez0z1=0;
    double distancez0surface=0;  
    if(TransitionBoundary!=0){
      if (fabs(z0)>TransitionBoundary && fabs(zn)>TransitionBoundary){
	if(zmax<=TransitionBoundary){
	  distancez0z1=fDnfR(-zn,&params6a) - fDnfR(-z0,&params6b);
	  distancez0surface=fDnfR(zmax,&params6c) - fDnfR(TransitionBoundary,&params6d) + fDnfR(TransitionBoundary+1e-6,&params6f) - fDnfR(-z0,&params6b);
	}else{
	  distancez0z1=fDnfR(-zn,&params6a) - fDnfR(-z0,&params6b);
	  distancez0surface=fDnfR(zmax,&params6c) - fDnfR(-z0,&params6b);
	}
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)<TransitionBoundary){
	distancez0z1=fDnfR(-zn,&params6a) - fDnfR(TransitionBoundary,&params6d) + fDnfR(TransitionBoundary+1e-6,&params6f) - fDnfR(-z0,&params6b);
	distancez0surface=fDnfR(zmax,&params6c) - fDnfR(TransitionBoundary,&params6d) + fDnfR(TransitionBoundary+1e-6,&params6f) - fDnfR(-z0,&params6b);
      }
      if (fabs(z0)<TransitionBoundary && fabs(zn)<TransitionBoundary){
	distancez0z1=fDnfR(-zn,&params6a) - fDnfR(-z0,&params6b);
	distancez0surface=fDnfR(zmax,&params6c) - fDnfR(-z0,&params6b);
      }
      if (fabs(z0)==TransitionBoundary && fabs(zn)<TransitionBoundary){
	distancez0z1=fDnfR(-zn,&params6a) - fDnfR(-z0,&params6b);
	distancez0surface=fDnfR(zmax,&params6c) - fDnfR(-z0,&params6b); 
      }
      if (fabs(z0)==TransitionBoundary && fabs(zn)==TransitionBoundary){
	distancez0z1=fDnfR(-zn,&params6a) - fDnfR(-z0,&params6b);
	distancez0surface=fDnfR(zmax,&params6c) - fDnfR(-z0,&params6b); 
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)==TransitionBoundary){
	distancez0z1=fDnfR(-zn,&params6a) - fDnfR(TransitionBoundary,&params6d) + fDnfR(TransitionBoundary+1e-6,&params6f) - fDnfR(-z0,&params6b);
	distancez0surface=fDnfR(zmax,&params6c) - fDnfR(TransitionBoundary,&params6d) + fDnfR(TransitionBoundary+1e-6,&params6f) - fDnfR(-z0,&params6b); 
      }
    }else{
      distancez0z1=fDnfR(-zn,&params6a) - fDnfR(-z0,&params6b);
      distancez0surface=fDnfR(zmax,&params6c) - fDnfR(-z0,&params6b);
    }

    xn= distancez0z1 - 2*(distancez0surface);
    
    checknan=fDnfR(-zn,&params6a);
    if(isnan(checknan)==false && zn<=0 && Flip==false){
      gr3->SetPoint(npnt,xn,zn);
      aoutRa<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(isnan(checknan)==false && zn<=0 && Flip==true){
      gr3->SetPoint(npnt,x1-xn,zn);
      aoutRa<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
    
    zn=zn+h;
    if(zn>-zmax){
      i=dmax+2;      
    }
  }

  /* Map out the 2nd part of the refracted ray */
  zn=-zmax;
  for(int i=0;i<dmax;i++){  
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueRa};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueRa};
    params6c = {A_ice, GetB(TransitionBoundary), GetC(TransitionBoundary), lvalueRa};
    params6d = {A_ice, GetB(TransitionBoundary+1e-6), GetC(TransitionBoundary+1e-6), lvalueRa}; 
   
    if(TransitionBoundary!=0){
      if (fabs(z0)<TransitionBoundary && fabs(zn)<TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)<TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(-TransitionBoundary,&params6c) + fDnfR(-(TransitionBoundary+1e-6),&params6d) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)>TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)==TransitionBoundary && fabs(zn)<TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)==TransitionBoundary && fabs(zn)==TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
      }
      if (fabs(z0)>TransitionBoundary && fabs(zn)==TransitionBoundary){
	xn=fDnfR(zn,&params6a) - fDnfR(-TransitionBoundary,&params6c) + fDnfR(-(TransitionBoundary+1e-6),&params6d) - fDnfR(z0,&params6b);
      }
    }else{
      xn=fDnfR(zn,&params6a) - fDnfR(z0,&params6b);
    }
    
    checknan=fDnfR(zn,&params6a);
    if(isnan(checknan)==false && Flip==false){
      gr3->SetPoint(npnt,xn,zn);
      aoutRa<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(isnan(checknan)==false && Flip==true){
      gr3->SetPoint(npnt,x1-xn,zn);
      aoutRa<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
    
    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;
    }
  }

  params6a = {A_ice, GetB(zn), GetC(zn), lvalueRa};
  params6b = {A_ice, GetB(z0), GetC(z0), lvalueRa};  
  xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
  if(Flip==true){
    gr3->SetPoint(npnt,x1-xn,z0);
    aoutRa<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
  }else{
    gr3->SetPoint(npnt,xn,z0);
    aoutRa<<npnt<<" "<<xn<<" "<<zn<<endl;
  }  
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  return gr3;
  
}

/* function for plotting and storing all the rays */
void PlotAndStoreRays(double x0,double z0, double z1, double x1, double zmax[2], double lvalues[4], double checkzeroes[4]){
  
  double lvalueD=lvalues[0];
  double lvalueR=lvalues[1];
  double lvalueRa[2]={lvalues[2],lvalues[3]};

  double checkzeroD=checkzeroes[0];
  double checkzeroR=checkzeroes[1];
  double checkzeroRa[2]={checkzeroes[2],checkzeroes[3]}; 

  TMultiGraph *mg=new TMultiGraph();
  
  TGraph *gr1=GetFullDirectRayPath(z0,x1,z1,lvalueD);
  TGraph *gr2=GetFullReflectedRayPath(z0,x1,z1,lvalueR);
  TGraph *gr3A=new TGraph();
  TGraph *gr3B=new TGraph();
  
  if((fabs(checkzeroR)>0.5 || fabs(checkzeroD)>0.5) && fabs(checkzeroRa[0])<0.5){  
    gr3A=GetFullRefractedRayPath(z0,x1,z1,zmax[0],lvalueRa[0],1);
    if(fabs(checkzeroRa[1])<0.5){ 
      gr3B=GetFullRefractedRayPath(z0,x1,z1,zmax[1],lvalueRa[1],2);
    }
  }

  gr1->SetMarkerColor(kBlue);
  gr2->SetMarkerColor(kBlue);
  gr3A->SetMarkerColor(kBlue);
  gr3B->SetMarkerColor(kBlue);
  
  /* Plot the all the possible ray paths on the canvas */
  TGraph *gr4=new TGraph();
  gr4->SetPoint(0,x1,z1);
  gr4->SetMarkerStyle(20);
  gr4->SetMarkerColor(kRed);

  TGraph *gr4b=new TGraph();
  gr4b->SetPoint(0,0,z0);
  gr4b->SetMarkerStyle(20);
  gr4b->SetMarkerColor(kGreen);
  
  double zlower=z0;
  double zupper=z1;
  if(fabs(z0)<fabs(z1)){
    zlower=z1;
    zupper=z0;
  }
  if(fabs(z0)>fabs(z1)){
    zlower=z0;
    zupper=z1;
  }
  TGraph *gr5=new TGraph();
  gr5->SetPoint(0,0,zupper+100);
  gr5->SetPoint(1,0,zlower-100);
  gr5->SetPoint(2,x1+50,0);

  if(fabs(checkzeroD)<0.5){
    mg->Add(gr1);
  }
  if(fabs(checkzeroR)<0.5){
    mg->Add(gr2);
  }
  if(fabs(checkzeroRa[0])<0.5){
    mg->Add(gr3A);
  }
  if(fabs(checkzeroRa[1])<0.5){
    mg->Add(gr3B);
  }
  
  mg->Add(gr4);
  mg->Add(gr4b);
  //mg->Add(gr5);
  
  TString title="Depth vs Distance, Tx at x=";
  title+=x0;
  title+=" m,z=";
  title+=(int)z0;
  title+=" m, Rx at x=";
  title+=x1;
  title+=" m,z=";
  title+=(int)z1;
  title+=" m; Distance (m);Depth (m)";
  mg->SetTitle(title);
  
  TCanvas *c2=new TCanvas("c2","c2");
  c2->cd();
  mg->Draw("ALP");
  mg->GetXaxis()->SetNdivisions(20);
  c2->SetGridx();
  c2->SetGridy();

}

double *IceRayTracing(double x0, double z0, double x1, double z1){

  /* define a pointer to give back the output of raytracing */ 
  double *output=new double[25];

  /* Store the ray paths in text files */
  bool PlotRayPaths=true;
  /* calculate the attenuation (not included yet!) */
  bool attcal=false;
  
  double Txcor[2]={x0,z0};/* Tx positions */
  double Rxcor[2]={x1,z1};/* Rx Positions */
  
  /*  ********This part of the code will try to get the Direct ray between Rx and Tx.********** */
  double* GetDirectRay=GetDirectRayPar(z0,x1,z1);
  double RangD=GetDirectRay[0];
  double LangD=GetDirectRay[1];
  double timeD=GetDirectRay[2];
  double lvalueD=GetDirectRay[3];
  double checkzeroD=GetDirectRay[4];
  delete []GetDirectRay;
  
  /* ********This part of the code will try to get the Reflected ray between Rx and Tx.********** */
  double* GetReflectedRay=GetReflectedRayPar(z0,x1,z1);
  double RangR=GetReflectedRay[0];
  double LangR=GetReflectedRay[1];
  double timeR=GetReflectedRay[2];
  double lvalueR=GetReflectedRay[3];
  double checkzeroR=GetReflectedRay[4];
  double timeR1=GetReflectedRay[5];
  double timeR2=GetReflectedRay[6];
  double AngleOfIncidenceInIce=GetReflectedRay[7];
  delete []GetReflectedRay;

  /* ********This part of the code will try to get the Refracted ray between Rx and Tx.********** */
  double RangRa[2]={0,0};
  double LangRa[2]={0,0};
  double timeRa[2]={0,0};
  double lvalueRa[2]={0,0}; 
  double checkzeroRa[2]={-1000,-1000};
  double timeRa1[2]={0,0};
  double timeRa2[2]={0,0};
  double zmax[2]={0,0};
  
  /* This if condition makes sure that we only try to find a refracted ray if we don't get two possible ray paths from the direct and reflected case. This saves us alot of time since we know that between each Tx and Rx position we only expect 2 rays. */
  if(fabs(checkzeroR)>0.5 || fabs(checkzeroD)>0.5){
    double* GetRefractedRay=GetRefractedRayPar(z0,x1,z1,LangR,RangR);
    RangRa[0]=GetRefractedRay[0];
    LangRa[0]=GetRefractedRay[1];
    timeRa[0]=GetRefractedRay[2];
    lvalueRa[0]=GetRefractedRay[3]; 
    checkzeroRa[0]=GetRefractedRay[4];
    timeRa1[0]=GetRefractedRay[5];
    timeRa2[0]=GetRefractedRay[6];
    zmax[0]=GetRefractedRay[7];

    if(fabs(checkzeroR)>0.5 && fabs(checkzeroD)>0.5){
      RangRa[1]=GetRefractedRay[8];
      LangRa[1]=GetRefractedRay[9];
      timeRa[1]=GetRefractedRay[10];
      lvalueRa[1]=GetRefractedRay[11]; 
      checkzeroRa[1]=GetRefractedRay[12];
      timeRa1[1]=GetRefractedRay[13];
      timeRa2[1]=GetRefractedRay[14];
      zmax[1]=GetRefractedRay[15];
    }

    delete []GetRefractedRay;
  }
  
  /* This part of the code can be used if the user wants to plot the individual ray paths. This part of the code prints out the individual ray paths in text files and also plots them on a canvas */
  if(PlotRayPaths==true){
    double lvalues[4];
    lvalues[0]=lvalueD;
    lvalues[1]=lvalueR;
    lvalues[2]=lvalueRa[0];
    lvalues[3]=lvalueRa[1];

    double checkzeroes[4];
    checkzeroes[0]=checkzeroD;
    checkzeroes[1]=checkzeroR;
    checkzeroes[2]=checkzeroRa[0];
    checkzeroes[3]=checkzeroRa[1];
    
    PlotAndStoreRays(x0,z0,z1,x1,zmax,lvalues,checkzeroes);
  }
  
  /* print out all the output from the code */
  cout<<0<<" ,x0= "<<x0<<" ,z0= "<<z0<<" ,x1= "<<x1<<" ,z1= "<<z1<<" ,langRa= "<<LangRa[0]<<" ,langR= "<<LangR<<" ,langD= "<<LangD<<" ,langD-langR= "<<LangD-LangR<<" ,langD-langRa= "<<LangD-LangRa[0]<<" ,RangRa= "<<RangRa[0]<<" ,RangR= "<<RangR<<" ,RangD= "<<RangD<<" ,RangR-RangD= "<<RangR-RangD<<" ,RangRa-RangD= "<<RangRa[0]-RangD<<" ,timeRa= "<<timeRa[0]<<" ,timeR= "<<timeRa[0]<<" ,timeD= "<<timeD<<" ,timeR-timeD= "<<timeR-timeD<<" ,timeRa-timeD= "<<timeRa[0]-timeD<<" ,lvalueRa "<<lvalueRa[0]<<" ,lvalueR "<<lvalueR<<" "<<" ,lvalueD "<<lvalueD<<" ,checkzeroRa "<<checkzeroRa[0]<<" ,checkzeroR "<<checkzeroR<<" ,checkzeroD "<<checkzeroD<<endl;

  cout<<0<<" ,x0= "<<x0<<" ,z0= "<<z0<<" ,x1= "<<x1<<" ,z1= "<<z1<<" ,langRa= "<<LangRa[1]<<" ,langR= "<<LangR<<" ,langD= "<<LangD<<" ,langD-langR= "<<LangD-LangR<<" ,langD-langRa= "<<LangD-LangRa[1]<<" ,RangRa= "<<RangRa[1]<<" ,RangR= "<<RangR<<" ,RangD= "<<RangD<<" ,RangR-RangD= "<<RangR-RangD<<" ,RangRa-RangD= "<<RangRa[1]-RangD<<" ,timeRa= "<<timeRa[1]<<" ,timeR= "<<timeRa[1]<<" ,timeD= "<<timeD<<" ,timeR-timeD= "<<timeR-timeD<<" ,timeRa-timeD= "<<timeRa[1]-timeD<<" ,lvalueRa "<<lvalueRa[1]<<" ,lvalueR "<<lvalueR<<" "<<" ,lvalueD "<<lvalueD<<" ,checkzeroRa "<<checkzeroRa[1]<<" ,checkzeroR "<<checkzeroR<<" ,checkzeroD "<<checkzeroD<<endl;

  cout<<zmax[0]<<" "<<zmax[1]<<endl;
  /* Fill in the output pointer after calculating all the results */
  output[0]=LangD;
  output[1]=LangR;
  output[2]=LangRa[0];
  output[3]=LangRa[1];
  output[4]=timeD;
  output[5]=timeR;
  output[6]=timeRa[0];
  output[7]=timeRa[1];
  output[8]=RangD;
  output[9]=RangR;
  output[10]=RangRa[0];
  output[11]=RangRa[1];
  
  /* fill in the output array part where you fill in the times for the two parts of the reflected or refracted rays */
  if(fabs(checkzeroR)<0.5){
    output[12]=timeR1;
    output[13]=timeR2;
  }
  
  if(fabs(checkzeroRa[0])<0.5){
    output[14]=timeRa1[0];
    output[15]=timeRa2[0];
  }

  if(fabs(checkzeroRa[1])<0.5){
    output[16]=timeRa1[1];
    output[17]=timeRa2[1];
  }
  
  output[18]=AngleOfIncidenceInIce;
  output[19]=lvalueD;
  output[20]=lvalueR;
  output[21]=lvalueRa[0];
  output[22]=lvalueRa[1];
  output[23]=zmax[0];  
  output[24]=zmax[1];
  
  /* Set the recieve angle to be zero for a ray which did not give us a possible path between Tx and Rx. I use this as a flag to determine which two rays gave me possible ray paths. */
  if(fabs(checkzeroD)>0.5){
    output[8]=-1000;
  }
  if(fabs(checkzeroR)>0.5){
    output[9]=-1000;
  }
  if(fabs(checkzeroRa[0])>0.5){
    output[10]=-1000;
  }
  if(fabs(checkzeroRa[1])>0.5){
    output[11]=-1000;
  } 
  
  return output;
}

/* Analytical solution describing ray paths in ice as function of depth for constant refractive index*/
double fDnfR_Cnz(double x,void *params){
  
  struct fDnfR_params *p= (struct fDnfR_params *) params;
  double A = p->a;
  double L = p->l;
  
  return (L/sqrt(A*A-L*L))*x;
}

/* Analytical solution describing the ray path in ice as a function of the L parameter for constant refractive index*/
double fDnfR_L_Cnz(double x,void *params){
  
  struct fDnfR_L_params *p= (struct fDnfR_L_params *) params;
  double A = p->a;
  double Z = p->z;
  
  double out=0;
  if(A>x){
    out=(x/sqrt(A*A-x*x))*Z;
  }else{
    out=tan(asin(x/A))*Z;
  }
  return out;
}

/* This function is minimised to find the launch angle (or the L parameter) for the reflected ray for constant refractive index*/
double fRa_Cnz(double x,void *params){
  struct fDanfRa_params *p= (struct fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  struct fDnfR_L_params params1a = {A, 0, 0, -z1};
  struct fDnfR_L_params params1b = {A, 0, 0, -z0};
  struct fDnfR_L_params params1c = {A, 0, 0, 0.0};

  return fDnfR_L_Cnz(x,&params1a) - fDnfR_L_Cnz(x,&params1b) - 2*( fDnfR_L_Cnz(x,&params1c) - fDnfR_L_Cnz(x,&params1b) ) - x1;
}

double fRa_Cnz_df(double x,void *params){
  gsl_function F;
  F.function = &fRa_Cnz;
  F.params = params;
 
  double result,abserr;
  gsl_deriv_central (&F, x, 1e-8, &result, &abserr);

  return result;
}

void fRa_Cnz_fdf (double x, void *params,double *y, double *dy){ 
  *y = fRa_Cnz(x,params);
  *dy = fRa_Cnz_df(x,params);
}

/* This functions works for the Direct ray and gives you back the launch angle, receive angle and propagation time of the ray together with values of the L parameter. This for constant refractive index*/
double* GetDirectRayPar_Cnz(double z0, double x1, double z1, double A_ice_Cnz){

  double *output=new double[4];
  
  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }

  
  /* Calculate the launch angle and the value of the L parameter */
  double LangD=(pi*0.5-atan(fabs(z1-z0)/x1))*(180.0/pi);
  double lvalueD=A_ice_Cnz*sin(LangD*(pi/180.0));
  double timeD=(sqrt( pow(x1,2) + pow(z1-z0,2) )/spedc)*A_ice_Cnz;
  
  /* Calculate the recieve angle for direct rays by which is the same as the launch angle */
  double RangD=LangD;
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangD;
  output[1]=LangD;
  output[2]=timeD;
  output[3]=lvalueD;

  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangD;
    output[1]=180-RangD;
  }
  
  return output;
}

/* This functions works for the Reflected ray and gives you back the launch angle, receive angle and propagation times (of the whole ray and the two direct rays that make it up) together with values of the L parameter. This is for constant refractive index*/
double *GetReflectedRayPar_Cnz(double z0, double x1 , double z1, double A_ice_Cnz){

  double *output=new double[8];

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* First we setup the fRa function that will be minimised to get the launch angle (or the L parameter) for the reflected ray. */
  gsl_function F3;
  struct fDanfRa_params params3= {A_ice_Cnz, z0, x1, z1};
  F3.function = &fRa_Cnz;
  F3.params = &params3;

  // gsl_function_fdf F3;
  // struct fDanfRa_params params3= {A_ice_Cnz, z0, x1, z1};
  // F3.f = &fRa_Cnz;
  // F3.df = &fRa_Cnz_df;
  // F3.fdf = &fRa_Cnz_fdf;
  // F3.params = &params3;

  /* In my raytracing solution given in the function fDnfR_Cnz the launch angle (or the L parameter) has limit placed on it by this part in the solution that L<A . This sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit to be the angle of the direct ray as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */
  double UpperLimitL=A_ice_Cnz*sin(pi*0.5-atan(fabs(z1-z0)/x1));

  /* Do the minimisation and get the value of the L parameter and the launch angle */
  double lvalueR=FindFunctionRoot(F3,0.0,UpperLimitL);
  double LangR=asin(lvalueR/A_ice_Cnz)*(180.0/pi);
  
  /* In the reflected case we basically have two direct rays 1) from Tx to surface 2) from surface to Rx. . Also get the time for the two individual direct rays separately */
  double z2=0,x2=fabs(z0)*tan(LangR*(pi/180));///coordinates of point of incidence in ice at the surface
  double timeR1=(sqrt( pow(x2,2) + pow(z2-z0,2) )/spedc)*A_ice_Cnz;
  double timeR2=(sqrt( pow(x2-x1,2) + pow(z2-z1,2) )/spedc)*A_ice_Cnz;
  double timeR= timeR1 + timeR2;
  
  /* flip the times back if the original positions were flipped */
  if(Flip==true){
    double dumR=timeR2;
    timeR2=timeR1;
    timeR1=dumR;
  }
  timeR1=timeR1;
  timeR2=timeR2;
  
  /* Calculate the recieve angle for reflected ray using simple geometry*/
  double RangR=180-LangR;

   /* Calculate the angle of incidence of the reflected ray at the surface ice. This will be used to calculate the Fresnel Coefficients.*/
  double IncidenceAngleInIce=LangR;
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangR;
  output[1]=LangR;
  output[2]=timeR;
  output[3]=lvalueR;
  output[4]=0;
  output[5]=timeR1;
  output[6]=timeR2;
  output[7]=IncidenceAngleInIce;
  
  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangR;
    output[1]=180-RangR;
  } 
  
  return output;
}

/* This function returns the x and z values for the full Direct ray path in a TGraph and also prints out the ray path in a text file. This is for a constant refractive index. */
TGraph* GetFullDirectRayPath_Cnz(double z0, double x1, double z1, double lvalueD, double A_ice_Cnz){

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
   
  /* Set the name of the text files */
  ofstream aoutD("DirectRay_Cnz.txt");
  /* Set the step size for plotting */
  double h=0.1;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;

  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;
  struct fDnfR_params params6a;
  struct fDnfR_params params6b;
  
  TGraph *gr1=new TGraph();
  for(int i=0;i<dmax;i++){
    params6a = {A_ice_Cnz, GetB(zn), GetC(zn), lvalueD};
    params6b = {A_ice_Cnz, GetB(z0), GetC(z0), lvalueD};
    xn=fDnfR_Cnz(zn,&params6a)-fDnfR_Cnz(z0,&params6b);
    checknan=fDnfR(zn,&params6a);
    if(isnan(checknan)==false && Flip==false){
      gr1->SetPoint(npnt,xn,zn);
      aoutD<<npnt<<" "<<xn<<" "<<zn<<endl;;
      npnt++;
    }

    if(isnan(checknan)==false && Flip==true){
      gr1->SetPoint(npnt,x1-xn,zn);
      aoutD<<npnt<<" "<<x1-xn<<" "<<zn<<endl;;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;      
    }  
  }

  params6a = {A_ice_Cnz, 0, 0, lvalueD};
  params6b = {A_ice_Cnz, 0, 0, lvalueD};
  xn=fDnfR_Cnz(zn,&params6a)-fDnfR_Cnz(z0,&params6b); 
  if(Flip==true){
    gr1->SetPoint(npnt,x1-xn,z0);
    aoutD<<npnt<<" "<<x1-xn<<" "<<zn<<endl;;
  }else{
    gr1->SetPoint(npnt,xn,z0);
    aoutD<<npnt<<" "<<xn<<" "<<zn<<endl;;
  }
  npnt++;
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  return gr1;
}

/* This function returns the x and z values for the full Reflected ray path in a TGraph and also prints out the ray path in a text file. This is for a constant refractive index. */
TGraph* GetFullReflectedRayPath_Cnz(double z0, double x1, double z1, double lvalueR, double A_ice_Cnz){

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  // /* Set the name of the text files */
  ofstream aoutR("ReflectedRay_Cnz.txt");
  /* Set the step size for plotting. */
  double h=0.1;
  /* Set the total steps required for looping over the whole ray path */
  int dmax=100000;
  /* Set the values to start the rays from */
  double zn=z1;
  double xn=0;
  
  /* Map out the direct ray path */
  int npnt=0;
  double checknan=0;  
  struct fDnfR_params params6a;
  struct fDnfR_params params6b;
  struct fDnfR_params params6c;

  /* Map out the 1st part of the reflected ray */
  TGraph *gr2=new TGraph();
  for(int i=0;i<dmax;i++){
    params6a = {A_ice_Cnz, 0, 0, lvalueR};
    params6b = {A_ice_Cnz, 0, 0, lvalueR};
    params6c = {A_ice_Cnz, 0, 0, lvalueR};
    xn=(fDnfR_Cnz(-zn,&params6a)-fDnfR_Cnz(-z0,&params6b)+2*fabs(fDnfR_Cnz(0.0,&params6c)-fDnfR_Cnz(-z0,&params6b)));
    checknan=fDnfR_Cnz(-zn,&params6a);
    if(isnan(checknan)==false && zn<=0 && Flip==false){
      gr2->SetPoint(npnt,xn,zn);
      aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(isnan(checknan)==false && zn<=0 && Flip==true){
      gr2->SetPoint(npnt,x1-xn,zn);
      aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
      
    zn=zn+h;
    if(zn>0){
      i=dmax+2;      
    }
  }
  
  /* Map out the 2nd part of the reflected ray */
  zn=0.0;
  for(int i=0;i<dmax;i++){  
    params6a = {A_ice_Cnz, 0, 0, lvalueR};
    params6b = {A_ice_Cnz, 0, 0, lvalueR};
    xn=fDnfR_Cnz(zn,&params6a)-fDnfR_Cnz(z0,&params6b);
    checknan=fDnfR_Cnz(zn,&params6a);
    if(isnan(checknan)==false && Flip==false){
      gr2->SetPoint(npnt,xn,zn);
      aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }
      
    if(isnan(checknan)==false && Flip==true){
      gr2->SetPoint(npnt,x1-xn,zn);
      aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;
    }
  }

  params6a = {A_ice_Cnz, 0, 0, lvalueR};
  params6b = {A_ice_Cnz, 0, 0, lvalueR};
  xn=fDnfR_Cnz(zn,&params6a)-fDnfR_Cnz(z0,&params6b);
  if(Flip==true){
    gr2->SetPoint(npnt,x1-xn,zn);
    aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
  }else{
    gr2->SetPoint(npnt,xn,zn);
    aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
  }
  
  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }

  return gr2;

}

/* function for plotting and storing all the rays. This is for constant refractive index. */
void PlotAndStoreRays_Cnz(double x0,double z0, double z1, double x1, double lvalues[2], double A_ice_Cnz){
  
  double lvalueD=lvalues[0];
  double lvalueR=lvalues[1];

  TMultiGraph *mg=new TMultiGraph();
  
  TGraph *gr1=GetFullDirectRayPath_Cnz(z0,x1,z1,lvalueD,A_ice_Cnz);
  TGraph *gr2=GetFullReflectedRayPath_Cnz(z0,x1,z1,lvalueR,A_ice_Cnz);
 
  gr1->SetMarkerColor(kBlue);
  gr2->SetMarkerColor(kBlue); 
  
  /* Plot the all the possible ray paths on the canvas */
  TGraph *gr4=new TGraph();
  gr4->SetPoint(0,x1,z1);
  gr4->SetMarkerStyle(20);
  gr4->SetMarkerColor(kRed);

  TGraph *gr4b=new TGraph();
  gr4b->SetPoint(0,0,z0);
  gr4b->SetMarkerStyle(20);
  gr4b->SetMarkerColor(kGreen);
    
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(2);

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerColor(2);

  double zlower=z0;
  if(fabs(z0)<fabs(z1)){
    zlower=z1;
  }
  if(fabs(z0)>fabs(z1)){
    zlower=z0;
  }
  TGraph *gr5=new TGraph();
  gr5->SetPoint(0,0,zlower-50);
  gr5->SetPoint(1,x1+50,0);

  mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr4);
  mg->Add(gr4b);
  //mg->Add(gr5);

  TString title="Depth vs Distance, Tx at x=";
  title+=x0;
  title+=" m,z=";
  title+=(int)z0;
  title+=" m, Rx at x=";
  title+=x1;
  title+=" m,z=";
  title+=(int)z1;
  title+=" m; Distance (m);Depth (m)";
  mg->SetTitle(title);
  
  TCanvas *cRay2=new TCanvas("cRay2","cRay2");
  cRay2->cd();
  mg->Draw("AP");
  mg->GetXaxis()->SetNdivisions(20);
  cRay2->SetGridx();
  cRay2->SetGridy();
}

/* This is the main raytracing function. x0 always has to be zero. z0 is the Tx depth in m and z1 is the depth of the Rx in m. Both depths are negative. x1 is the distance between them. This functions works for a constant refractive index */
double *IceRayTracing_Cnz(double x0, double z0, double x1, double z1, double A_ice_Cnz){

  /* define a pointer to give back the output of raytracing */ 
  double *output=new double[9];

  /* Store the ray paths in text files */
  bool PlotRayPaths=false;
  
  /*  ********This part of the code will try to get the Direct ray between Rx and Tx.********** */
  double* GetDirectRay=GetDirectRayPar_Cnz(z0,x1,z1,A_ice_Cnz);
  double RangD=GetDirectRay[0];
  double LangD=GetDirectRay[1];
  double timeD=GetDirectRay[2];
  double lvalueD=GetDirectRay[3];
  double checkzeroD=GetDirectRay[4];
  delete []GetDirectRay;
  
  /* ********This part of the code will try to get the Reflected ray between Rx and Tx.********** */
  double* GetReflectedRay=GetReflectedRayPar_Cnz(z0,x1,z1,A_ice_Cnz);
  double RangR=GetReflectedRay[0];
  double LangR=GetReflectedRay[1];
  double timeR=GetReflectedRay[2];
  double lvalueR=GetReflectedRay[3];
  double checkzeroR=GetReflectedRay[4];
  double timeR1=GetReflectedRay[5];
  double timeR2=GetReflectedRay[6];
  double AngleOfIncidenceInIce=GetReflectedRay[7];
  delete []GetReflectedRay; 

  /* This part of the code can be used if the user wants to plot the individual ray paths. This part of the code prints out the individual ray paths in text files and also plots them on a canvas */
  if(PlotRayPaths==true){
    double lvalues[2];
    lvalues[0]=lvalueD;
    lvalues[1]=lvalueR;
    
    PlotAndStoreRays_Cnz(x0,z0,z1,x1,lvalues,A_ice_Cnz);
  }  

  /* Fill in the output pointer after calculating all the results */
  output[0]=LangD;
  output[1]=LangR;
  output[2]=timeD;
  output[3]=timeR;
  output[4]=RangD;
  output[5]=RangR;
  
  /* fill in the output array part where you fill in the times for the two parts of the reflected or refracted rays */  
  output[6]=timeR1;
  output[7]=timeR2;
  output[8]=AngleOfIncidenceInIce;

  
  return output;
}


void IceRayTracing_wROOTplot(double xTx, double yTx, double zTx, double xRx, double yRx, double zRx ){

  ////For recording how much time the process took
  auto t1b = std::chrono::high_resolution_clock::now();
  
  double TxCor[3]={xTx,yTx,zTx};
  double RxCor[3]={xRx,yRx,zRx};
  
  ////An example for Direct and Reflected Rays
  // double TxCor[3]={0,0,-100};
  // double RxCor[3]={0,100,-5};

  ////An example for Direct and Refracted Rays
  // double TxCor[3]={0,0,-67.5489};
  // double RxCor[3]={0,1338.3,-800};
  
  double x0=0;/////always has to be zero
  double z0=TxCor[2];
  double x1=sqrt(pow(TxCor[0]-RxCor[0],2)+pow(TxCor[1]-RxCor[1],2));
  double z1=RxCor[2];

  double * getresults=IceRayTracing(x0,z0,x1,z1);

  cout<<"*******For the Direct Ray********"<<endl;
  cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
  cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
  cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;
  cout<<"*******For the Refracted[1] Ray********"<<endl;
  cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
  cout<<"Recieve Angle: "<<getresults[10]<<" deg"<<endl;
  cout<<"Propogation Time: "<<getresults[6]*pow(10,9)<<" ns"<<endl;
  cout<<"*******For the Refracted[2] Ray********"<<endl;
  cout<<"Launch Angle: "<<getresults[3]<<" deg"<<endl;
  cout<<"Recieve Angle: "<<getresults[11]<<" deg"<<endl;
  cout<<"Propogation Time: "<<getresults[7]*pow(10,9)<<" ns"<<endl;
  cout<<"*******For the Reflected Ray********"<<endl;
  cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
  cout<<"Recieve Angle: "<<getresults[9]<<" deg"<<endl;
  cout<<"Propogation Time: "<<getresults[5]*pow(10,9)<<" ns"<<endl;   
  cout<<"Incident Angle in Ice on the Surface: "<<getresults[18]<<" deg"<<endl;
  
  cout<<" "<<endl;
  if(getresults[8]!=-1000 && getresults[10]!=-1000){ 
    cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;
    cout<<"Direct and Refracted[1]: dt(D,R)="<<(getresults[6]-getresults[4])*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Direct Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Refracted Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[10]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[6]*pow(10,9)<<" ns"<<endl;
  }

  if(getresults[8]!=-1000 && getresults[11]!=-1000){ 
    cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;
    cout<<"Direct and Refracted[2]: dt(D,R)="<<(getresults[7]-getresults[4])*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Direct Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Refracted Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[3]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[11]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[7]*pow(10,9)<<" ns"<<endl;
  }
      
  if(getresults[8]!=-1000 && getresults[9]!=-1000){
    cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;
    cout<<"Direct and Reflected: dt(D,R)="<<(getresults[5]-getresults[4])*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Direct Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Reflected Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[9]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[5]*pow(10,9)<<" ns"<<endl;   
    cout<<"Incident Angle in Ice on the Surface: "<<getresults[18]<<" deg"<<endl;
  }
      
  if(getresults[10]!=-1000 && getresults[9]!=-1000){
    cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;
    cout<<"Refracted[1] and Reflected: dt(D,R)="<<(getresults[5]-getresults[6])*pow(10,9)<<" ns "<<endl;
    cout<<"*******For the Reflected Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[9]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[5]*pow(10,9)<<" ns"<<endl;
    cout<<"Incident Angle on Ice Surface: "<<getresults[18]<<" deg"<<endl;
    cout<<"*******For the Refracted Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[10]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[6]*pow(10,9)<<" ns"<<endl;
  }

  if(getresults[11]!=-1000 && getresults[9]!=-1000){
    cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;
    cout<<"Refracted[2] and Reflected: dt(D,R)="<<(getresults[5]-getresults[7])*pow(10,9)<<" ns "<<endl;
    cout<<"*******For the Reflected Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[9]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[5]*pow(10,9)<<" ns"<<endl;
    cout<<"Incident Angle on Ice Surface: "<<getresults[18]<<" deg"<<endl;
    cout<<"*******For the Refracted Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[3]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[11]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[7]*pow(10,9)<<" ns"<<endl;
  }
 
  if(getresults[10]!=-1000 && getresults[11]!=-1000){
    cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;
    cout<<"Refracted[1] and Refracted[2]: dt(D,R)="<<(getresults[6]-getresults[7])*pow(10,9)<<" ns "<<endl;
    cout<<"*******For the Refracted[1] Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[10]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[6]*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Refracted[2] Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[3]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[11]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[7]*pow(10,9)<<" ns"<<endl;
  }

  // double * getresults2=IceRayTracing_Cnz(x0,z0,x1,z1,1.4);

  // cout<<" "<<endl;
  // cout<<" For constant refractive index"<<endl;
  // cout<<"Direct and Reflected: dt(D,R)="<<(getresults2[3]-getresults2[2])*pow(10,9)<<" ns"<<endl;
  // cout<<"*******For the Direct Ray********"<<endl;
  // cout<<"Launch Angle: "<<getresults2[0]<<" deg"<<endl;
  // cout<<"Recieve Angle: "<<getresults2[4]<<" deg"<<endl;
  // cout<<"Propogation Time: "<<getresults2[2]*pow(10,9)<<" ns"<<endl;
  // cout<<"*******For the Reflected Ray********"<<endl;
  // cout<<"Launch Angle: "<<getresults2[1]<<" deg"<<endl;
  // cout<<"Recieve Angle: "<<getresults2[5]<<" deg"<<endl;
  // cout<<"Propogation Time: "<<getresults2[3]*pow(10,9)<<" ns"<<endl;
  // cout<<"Incident Angle in Ice on the Surface: "<<getresults2[8]<<" deg"<<endl;

  // delete []getresults2; 
  
  auto t2b = std::chrono::high_resolution_clock::now();
  double durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();

  double Duration=durationb/1000;
  cout<<"total time taken by the script: "<<Duration<<" ms"<<endl; 

  // TGraph * gr=new TGraph();
  // struct fDanfRa_params params4= {A_ice, z0, x1, z1};
  // for(double i=0;i<1000;i++){
  //   double lvalue=1.3+0.001*i;
  //   double checkzeroRa=fRaa(lvalue,&params4);
  //   gr->SetPoint(i,lvalue,fabs(checkzeroRa));
  //   //cout<<i<<" "<<lvalue<<" "<<fabs(checkzeroRa)<<endl;
  //   if(lvalue>1.78){
  //     break;
  //   }
  // }
  // gr->Draw("ALP");
  
  delete []getresults;
  
}
