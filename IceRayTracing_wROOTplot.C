#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <chrono>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <sys/time.h>

#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TVector3.h"

using namespace std;

const double pi=4.0*atan(1.0); /* Gives back value of Pi */
const double spedc=299792458.0; /* Speed of Light in m/s */

/* Set the value of the asymptotic parameter of the refractive index model */
const double A_ice=1.78;

/* Get the value of the B parameter for the refractive index model */
double GetB(double z){
  z=fabs(z);
  double B=0;

  //B=-0.43;
  B=-0.48;
  return B;
}

/* Get the value of the C parameter for the refractive index model */
double GetC(double z){
  z=fabs(z);
  double C=0;

  //C=0.0132;
  C=0.02;
  return C;
}

/* Get the value of refractive index model for a given depth  */
double Getnz(double z){
  z=fabs(z);
  return A_ice+GetB(z)*exp(-GetC(z)*z);
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

/* Use GSL minimiser which uses Brent's Method to find root for a given function. This will be used to find roots wherever it is needed in my code. */
double FindFunctionRoot(gsl_function F,double x_lo, double x_hi)
{
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;

  T = gsl_root_fsolver_brent;
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
      status = gsl_root_test_interval (x_lo, x_hi,0, 0.000001);
	
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
  double zmax=FindFunctionRoot(F1,0.0,5000);
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

/* This function is minimised to find the launch angle (or the L parameter) for the direct ray */
struct fDanfRa_params { double a, z0, x1, z1; };
double fDa(double x,void *params){
  struct fDanfRa_params *p= (struct fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  struct fDnfR_L_params params1a = {A, GetB(z1), GetC(z1), z1};
  struct fDnfR_L_params params1b = {A, GetB(z0), GetC(z0), z0};
  
  return fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b) - x1;
}

/* This function is minimised to find the launch angle (or the L parameter) for the reflected ray */
double fRa(double x,void *params){
  struct fDanfRa_params *p= (struct fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  struct fDnfR_L_params params1a = {A, GetB(z1), -GetC(z1), -z1};
  struct fDnfR_L_params params1b = {A, GetB(z0), -GetC(z0), -z0};
  struct fDnfR_L_params params1c = {A, GetB(0.0000001), -GetC(0.0000001), 0.0000001};

  return fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b) - 2*( fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b) ) - x1;
}

/* This function is minimised to find the launch angle (or the L parameter) for the refracted ray */
double fRaa(double x,void *params){
  struct fDanfRa_params *p= (struct fDanfRa_params *) params;
  double A = p->a;
  double z0 = p->z0;
  double x1 = p->x1;
  double z1 = p->z1;

  double zmax= GetZmax(A,x)+0.0000001;

  struct fDnfR_L_params params1a = {A, GetB(z1), -GetC(z1), -z1};
  struct fDnfR_L_params params1b = {A, GetB(z0), -GetC(z0), -z0};
  struct fDnfR_L_params params1c = {A, GetB(zmax), -GetC(zmax), zmax};

  return fDnfR_L(x,&params1a) - fDnfR_L(x,&params1b) - 2*( fDnfR_L(x,&params1c) - fDnfR_L(x,&params1b) ) - x1;
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

  /* In my raytracing solution given in the function fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0 and z1 and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */ 
  double UpLimnz[]={Getnz(z1),Getnz(z0)};
  double* UpperLimitL=min_element(UpLimnz,UpLimnz+2);

  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function. */
  double lvalueD=FindFunctionRoot(F1,0.0000001,UpperLimitL[0]);
  double LangD=asin(lvalueD/Getnz(z0))*(180.0/pi);
  double checkzeroD=fDa(lvalueD,&params1);

  /* Get the propagation time for the direct ray using the ftimeD function after we have gotten the value of the L parameter. */
  struct ftimeD_params params2a = {A_ice, GetB(z0), -GetC(z0), spedc,lvalueD};
  struct ftimeD_params params2b = {A_ice, GetB(z1), -GetC(z1), spedc,lvalueD};

  /* we do the subtraction because we are measuring the time taken between the Tx and Rx positions */
  double timeD=+ftimeD(-z0,&params2a) - ftimeD(-z1,&params2b);

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

  /* In my raytracing solution given in the function fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0, z1 and also 0.0000001 m and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */
  double UpLimnz[]={Getnz(z1),Getnz(z0),Getnz(0.0000001)};
  double *UpperLimitL=min_element(UpLimnz,UpLimnz+3);
  
  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fRa function. */
  double lvalueR=FindFunctionRoot(F3,0.0000001,UpperLimitL[0]);
  double LangR=asin(lvalueR/Getnz(z0))*(180.0/pi);
  double checkzeroR=fRa(lvalueR,&params3); 

  /* Get the propagation time for the reflected ray using the ftimeD function after we have gotten the value of the L parameter. */
  struct ftimeD_params params3a = {A_ice, GetB(z0), GetC(z0), spedc,lvalueR};
  struct ftimeD_params params3b = {A_ice, GetB(z1), GetC(z1), spedc,lvalueR};
  struct ftimeD_params params3c = {A_ice, GetB(0.0000001), GetC(0.0000001), spedc,lvalueR};

  /* We do the subtraction because we are measuring the time taken between the Tx and Rx positions. In the reflected case we basically have two direct rays 1) from Tx to surface 2) from surface to Rx. Also get the time for the two individual direct rays separately */
  double timeR1=ftimeD(-0.0000001,&params3c) - ftimeD(z0,&params3a);
  double timeR2=ftimeD(-0.0000001,&params3c) - ftimeD(z1,&params3b);
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

  /* Calculate the angle of incidence of the reflected ray at the surface ice. This will be used to calculate the Fresnel Coefficients. The angle is calculated by calculating the derivative of the ray path fucnction at the surface*/
  struct fDnfR_params paramsIAngB = {A_ice, GetB(0.0000001), GetC(0.0000001), lvalueR};
  F5.function = &fDnfR; 
  F5.params = &paramsIAngB;
  gsl_deriv_central (&F5, -0.0000001, 1e-8, &result, &abserr);
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

  double *output=new double[8];

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
  double lvalueR=sin(LangR*(pi/180))*Getnz(z0);
  double lvalueRa=0;
  double LangRa=0;
  double checkzeroRa=-1000;

  double timeRa=0;
  double timeRa1=0;
  double timeRa2=0;
  double raytime=0;
  double RangRa=0;
  double zmax=10;

  /* In my raytracing solution given in the function fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0, z1 and also 0.0000001 m and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. */
  double UpLimnz[]={Getnz(z1),Getnz(z0)};
  double *UpperLimitL=min_element(UpLimnz,UpLimnz+2);
  
  /* First we setup the fRa function that will be minimised to get the launch angle (or the L parameter) for the refracted ray. */
  gsl_function F4;
  struct fDanfRa_params params4= {A_ice, z0, x1, z1};
  F4.function = &fRaa;
  F4.params = &params4;

  /* Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fRaa function. The thing to note here is the lower limit of the minimisation function is set to the L value corresponding to the reflected ray launch angle. Since we know the refracted ray always has bigger launch angle the reflected ray this reduces our range and makes the function more efficient at finding the refracted ray launch angle. */
  lvalueRa=FindFunctionRoot(F4,Getnz(z0)*sin((LangR*(pi/180.0))),UpperLimitL[0]);
  LangRa=asin(lvalueRa/Getnz(z0))*(180.0/pi);
  checkzeroRa=fRaa(lvalueRa,&params4);

  /* If the above strategy did not work then we start decreasing the reflected ray launch angle in steps of 5 degree and increase our range for minimisation to find the launch angle (or the L parameter). Sometimes the refracted and reflected rays come out to be the same in that case also I forced my solution to try harder by changing the minimisation range. */
  double iangstep=5;
  while( (isnan(checkzeroRa)==true || fabs(checkzeroRa)>0.5 || fabs(lvalueRa-lvalueR)<pow(10,-5) || fabs(LangRa-LangR)<pow(10,-1)) && LangR>iangstep && iangstep<90){
    //cout<<"2nd try to get Refracted ray "<<isnan(checkzeroRa)<<" "<<fabs(checkzeroRa)<<endl;
    lvalueRa=FindFunctionRoot(F4,Getnz(z0)*sin(((LangR-iangstep)*(pi/180.0))),UpperLimitL[0]);
    LangRa=asin(lvalueRa/Getnz(z0))*(180.0/pi);
    checkzeroRa=fRaa(lvalueRa,&params4);
    iangstep=iangstep+5;
  }///end the second attempt    

  if(fabs(LangRa-LangR)<pow(10,-1)){
    iangstep=5;
    while( (isnan(checkzeroRa)==true || fabs(checkzeroRa)>0.5 || fabs(lvalueRa-lvalueR)<pow(10,-5) || fabs(LangRa-LangR)<pow(10,-1)) && LangR>iangstep && iangstep<90){
      //cout<<"2nd try to get Refracted ray "<<isnan(checkzeroRa)<<" "<<fabs(checkzeroRa)<<endl;
      lvalueRa=FindFunctionRoot(F4,Getnz(z0)*sin(((LangR-iangstep)*(pi/180.0))),UpperLimitL[0]+0.01);
      LangRa=asin(lvalueRa/Getnz(z0))*(180.0/pi);
      checkzeroRa=fRaa(lvalueRa,&params4);
      iangstep=iangstep+5;
    }///end the second attempt    
  }
  
  /* If we still did not find a refracted ray then set the check zero parameter to 1000 to make sure my code does not output this as a possible solution */
  if(isnan(checkzeroRa)==true || LangRa==0){
    checkzeroRa=-1000;
  }

  /* If we did find a possible refracted ray then now we need to find the depth at which the ray turns back down without hitting the surface. */
  zmax=GetZmax(A_ice,lvalueRa)+0.0000001;

  if(fabs(LangR-LangRa)<0.1){
    cout<<LangR<<" "<<LangRa<<" "<<LangR-LangRa<<" "<<zmax<<" "<<Getnz(z1)<<" "<<Getnz(z0)<<" "<<UpperLimitL[0]<<endl;
  }

  
  /* If the turning point depth also came out to be zero then now we are sure that there is no refracted ray */
  if(zmax==0.0000001){
    checkzeroRa=-1000;
  }

  /* Set parameters for ftimeD function to get the propagation time for the refracted ray */
  struct ftimeD_params params4a = {A_ice, GetB(z0), GetC(z0), spedc,lvalueRa};
  struct ftimeD_params params4b = {A_ice, GetB(z1), GetC(z1), spedc,lvalueRa};
  struct ftimeD_params params4c = {A_ice, GetB(zmax), GetC(zmax), spedc,lvalueRa};

  /* This if condition checks if the function has not gone crazy and given us a turning point of the ray which is lower than both Tx and Rx and is shallower in depth than both */
  if((z0<-zmax || zmax<-z1)){
    /* we do the subtraction because we are measuring the time taken between the Tx and Rx positions. In the refracted case we basically have two direct rays 1) from Tx to turning point 2) from turning point to Rx. Also get the time for the two individual direct rays separately */
    timeRa1=ftimeD(-zmax,&params4c) - ftimeD(z0,&params4a);
    timeRa2=ftimeD(-zmax,&params4c) - ftimeD(z1,&params4b);
    raytime=timeRa1 + timeRa2;
    if(Flip==true){
      double dumRa=timeRa2;
      timeRa2=timeRa1;
      timeRa1=dumRa;
    }
    timeRa1=timeRa1;
    timeRa2=timeRa2;
  }
  timeRa=raytime;

  /* Setup the function that will be used to calculate the angle of reception for all the rays */
  gsl_function F5;
  struct fDnfR_params params5c = {A_ice, GetB(z1), GetC(z1), lvalueRa};
  double result, abserr;
  F5.function = &fDnfR;
    
  /* Calculate the recieve angle for refacted ray by calculating the derivative of the function at the Rx position */
  F5.params = &params5c;
  gsl_deriv_central (&F5, z1, 1e-8, &result, &abserr);
  RangRa=180-atan(result)*(180.0/pi);

  /* When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line. */
  if(z1==z0 && isnan(RangRa)==true){
    RangRa=180-LangRa;
  }

  /* This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation. */
  if(z1!=z0 && isnan(RangRa)==true){
    RangRa=90;
  }

  dsw=0;
  /* If the Tx and Rx depth were switched then put them back to their original position */
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }
  
  output[0]=RangRa;
  output[1]=LangRa;
  output[2]=timeRa;
  output[3]=lvalueRa;
  output[4]=checkzeroRa;
  output[5]=timeRa1;
  output[6]=timeRa2;
  output[7]=zmax;
  
  /* If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user. */
  if(Flip==true){
    output[0]=180-LangRa;
    output[1]=180-RangRa;
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
  
  TGraph *gr1=new TGraph();
  for(int i=0;i<dmax;i++){
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueD};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueD};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
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

  if(Flip==true){
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueD};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueD};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b); 
    gr1->SetPoint(npnt,x1-xn,z0);
  }else{
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueD};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueD};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
    gr1->SetPoint(npnt,xn,z0);
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

  /* Map out the 1st part of the reflected ray */
  TGraph *gr2=new TGraph();
  for(int i=0;i<dmax;i++){
    params6a = {A_ice, GetB(zn), -GetC(zn), lvalueR};
    params6b = {A_ice, GetB(z0), -GetC(z0), lvalueR};
    params6c = {A_ice, GetB(0.0000001), -GetC(0.0000001), lvalueR};
    xn=(fDnfR(-zn,&params6a)-fDnfR(-z0,&params6b)+2*fabs(fDnfR(0.0000001,&params6c)-fDnfR(-z0,&params6b)));
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
  zn=-0.0000001;
  for(int i=0;i<dmax;i++){  
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueR};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueR};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
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
  
  if(Flip==true){
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueR};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueR};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
    gr2->SetPoint(npnt,x1-xn,zn);
  }else{
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueR};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueR};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
    gr2->SetPoint(npnt,xn,zn);
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
TGraph* GetFullRefractedRayPath(double z0, double x1, double z1, double zmax, double lvalueRa){

  /* My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end */
  bool Flip=false;
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }
  
  /* Set the name of the text files */
  ofstream aoutRa("RefractedRay.txt");
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

  /* Map out the 1st part of the refracted ray */
  TGraph *gr3=new TGraph();    
  for(int i=0;i<dmax;i++){
    params6a = {A_ice, GetB(zn), -GetC(zn), lvalueRa};
    params6b = {A_ice, GetB(z0), -GetC(z0), lvalueRa};
    params6c = {A_ice, GetB(zmax), -GetC(zmax), lvalueRa};
    xn=(fDnfR(-zn,&params6a)-fDnfR(-z0,&params6b)+2*fabs(fDnfR(zmax,&params6c)-fDnfR(-z0,&params6b)));
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
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
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
  
  if(Flip==true){
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueRa};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueRa};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b); 
    gr3->SetPoint(npnt,x1-xn,z0);
  }else{
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueRa};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueRa};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
    gr3->SetPoint(npnt,xn,z0);
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
void PlotAndStoreRays(double x0,double z0, double z1, double x1, double zmax, double lvalues[3], double checkzeroes[3]){
  
  double lvalueD=lvalues[0];
  double lvalueR=lvalues[1];
  double lvalueRa=lvalues[2];

  double checkzeroD=checkzeroes[0];
  double checkzeroR=checkzeroes[1];
  double checkzeroRa=checkzeroes[2]; 

  TMultiGraph *mg=new TMultiGraph();
  
  TGraph *gr1=GetFullDirectRayPath(z0,x1,z1,lvalueD);
  TGraph *gr2=GetFullReflectedRayPath(z0,x1,z1,lvalueR);
  TGraph *gr3=new TGraph();
  if((fabs(checkzeroR)>0.5 || fabs(checkzeroD)>0.5) && fabs(checkzeroRa)<0.5){
    gr3=GetFullRefractedRayPath(z0,x1,z1,zmax,lvalueRa);
  }

  gr1->SetMarkerColor(kBlue);
  gr2->SetMarkerColor(kBlue);
  gr3->SetMarkerColor(kBlue);
  
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
  if(fabs(z0)<fabs(z1)){
    zlower=z1;
  }
  if(fabs(z0)>fabs(z1)){
    zlower=z0;
  }
  TGraph *gr5=new TGraph();
  gr5->SetPoint(0,0,zlower-50);
  gr5->SetPoint(1,x1+50,0);

  if(fabs(checkzeroD)<0.5){
    mg->Add(gr1);
  }
  if(fabs(checkzeroR)<0.5){
    mg->Add(gr2);
  }
  if(fabs(checkzeroRa)<0.5){
    mg->Add(gr3);
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
  double *output=new double[12+4];

  /* Store the ray paths in text files */
  bool PlotRayPaths=false;
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
  double RangRa=0;
  double LangRa=0;
  double timeRa=0;
  double lvalueRa=0; 
  double checkzeroRa=-1000;
  double timeRa1=0;
  double timeRa2=0;
  double zmax=0;
  
  /* This if condition makes sure that we only try to find a refracted ray if we don't get two possible ray paths from the direct and reflected case. This saves us alot of time since we know that between each Tx and Rx position we only expect 2 rays. */
  if(fabs(checkzeroR)>0.5 || fabs(checkzeroD)>0.5){
    double* GetRefractedRay=GetRefractedRayPar(z0,x1,z1,LangR,RangR);
    RangRa=GetRefractedRay[0];
    LangRa=GetRefractedRay[1];
    timeRa=GetRefractedRay[2];
    lvalueRa=GetRefractedRay[3]; 
    checkzeroRa=GetRefractedRay[4];
    timeRa1=GetRefractedRay[5];
    timeRa2=GetRefractedRay[6];
    zmax=GetRefractedRay[7];
    delete []GetRefractedRay;
  }

  /* Fill in the output pointer after calculating all the results */
  output[0]=LangD;
  output[1]=LangR;
  output[2]=LangRa;
  output[3]=timeD;
  output[4]=timeR;
  output[5]=timeRa;
  output[6]=RangD;
  output[7]=RangR;
  output[8]=RangRa;

  /* This part of the code can be used if the user wants to plot the individual ray paths. This part of the code prints out the individual ray paths in text files and also plots them on a canvas */
  if(PlotRayPaths==true){
    double lvalues[3];
    lvalues[0]=lvalueD;
    lvalues[1]=lvalueR;
    lvalues[2]=lvalueRa;

    double checkzeroes[3];
    checkzeroes[0]=checkzeroD;
    checkzeroes[1]=checkzeroR;
    checkzeroes[2]=checkzeroRa;
    
    PlotAndStoreRays(x0,z0,z1,x1,zmax,lvalues,checkzeroes);
  }

  if(fabs(LangR-LangRa)<0.1){
    cout<<LangR<<" "<<LangRa<<" "<<LangR-LangRa<<endl;
  }
  
  /* print out all the output from the code */
  //cout<<0<<" ,x0= "<<x0<<" ,z0= "<<z0<<" ,x1= "<<x1<<" ,z1= "<<z1<<" ,langRa= "<<output[2]<<" ,langR= "<<output[1]<<" ,langD= "<<output[0]<<" ,langD-langR= "<<output[0]-output[1]<<" ,langD-langRa= "<<output[0]-output[2]<<" ,RangRa= "<<output[8]<<" ,RangR= "<<output[7]<<" ,RangD= "<<output[6]<<" ,RangR-RangD= "<<output[7]-output[6]<<" ,RangRa-RangD= "<<output[8]-output[6]<<" ,timeRa= "<<output[5]<<" ,timeR= "<<output[4]<<" ,timeD= "<<output[3]<<" ,timeR-timeD= "<<output[4]-output[3]<<" ,timeRa-timeD= "<<output[5]-output[3]<<" ,lvalueRa "<<lvalueRa<<" ,lvalueR "<<lvalueR<<" "<<" ,lvalueD "<<lvalueD<<" ,checkzeroRa "<<checkzeroRa<<" ,checkzeroR "<<checkzeroR<<" ,checkzeroD "<<checkzeroD<<endl;

  /* fill in the output array part where you fill in the times for the two parts of the reflected or refracted rays */
  if(fabs(checkzeroR)<0.5){
    output[9]=timeR1;
    output[10]=timeR2;
  }
  
  if(fabs(checkzeroRa)<0.5){
    output[9]=timeRa1;
    output[10]=timeRa2;
  }

  output[11]=AngleOfIncidenceInIce;
  output[12]=lvalueD;
  output[13]=lvalueR;
  output[14]=lvalueRa;
  output[15]=zmax;
  
  /* Set the recieve angle to be zero for a ray which did not give us a possible path between Tx and Rx. I use this as a flag to determine which two rays gave me possible ray paths. */
  if(fabs(checkzeroD)>0.5){
    output[6]=-1000;
  }
  if(fabs(checkzeroR)>0.5){
    output[7]=-1000;
  }
  if(fabs(checkzeroRa)>0.5 || (fabs(checkzeroD)<0.5 && fabs(checkzeroR)<0.5)){
    output[8]=-1000;
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
  //ofstream aoutD("DirectRay.txt");
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
      //aoutD<<npnt<<" "<<xn<<" "<<zn<<endl;;
      npnt++;
    }

    if(isnan(checknan)==false && Flip==true){
      gr1->SetPoint(npnt,x1-xn,zn);
      //aoutD<<npnt<<" "<<x1-xn<<" "<<zn<<endl;;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;      
    }  
  }

  if(Flip==true){
    params6a = {A_ice_Cnz, 0, 0, lvalueD};
    params6b = {A_ice_Cnz, 0, 0, lvalueD};
    xn=fDnfR_Cnz(zn,&params6a)-fDnfR_Cnz(z0,&params6b); 
    gr1->SetPoint(npnt,x1-xn,z0);
  }else{
    params6a = {A_ice_Cnz, 0, 0, lvalueD};
    params6b = {A_ice_Cnz, 0, 0, lvalueD};
    xn=fDnfR_Cnz(zn,&params6a)-fDnfR_Cnz(z0,&params6b);
    gr1->SetPoint(npnt,xn,z0);
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
  // ofstream aoutR("ReflectedRay.txt");
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
      //aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(isnan(checknan)==false && zn<=0 && Flip==true){
      gr2->SetPoint(npnt,x1-xn,zn);
      //aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
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
      //aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }
      
    if(isnan(checknan)==false && Flip==true){
      gr2->SetPoint(npnt,x1-xn,zn);
      //aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;
    }
  }
  
  if(Flip==true){
    params6a = {A_ice_Cnz, 0, 0, lvalueR};
    params6b = {A_ice_Cnz, 0, 0, lvalueR};
    xn=fDnfR_Cnz(zn,&params6a)-fDnfR_Cnz(z0,&params6b);
    gr2->SetPoint(npnt,x1-xn,zn);
  }else{
    params6a = {A_ice_Cnz, 0, 0, lvalueR};
    params6b = {A_ice_Cnz, 0, 0, lvalueR};
    xn=fDnfR_Cnz(zn,&params6a)-fDnfR_Cnz(z0,&params6b);
    gr2->SetPoint(npnt,xn,zn);
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
  bool PlotRayPaths=true;
  
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

  cout<<" "<<endl;
  if(getresults[6]!=-1000 && getresults[8]!=-1000){ 
    cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;
    cout<<"Direct and Refracted: dt(D,R)="<<(getresults[5]-getresults[3])*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Direct Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[6]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[3]*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Refracted Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[5]*pow(10,9)<<" ns"<<endl;
  }
      
  if(getresults[6]!=-1000 && getresults[7]!=-1000){
    cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;
    cout<<"Direct and Reflected: dt(D,R)="<<(getresults[4]-getresults[3])*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Direct Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[6]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[3]*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Reflected Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[7]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;   
    cout<<"Incident Angle in Ice on the Surface: "<<getresults[11]<<" deg"<<endl;
  }
      
  if(getresults[8]!=-1000 && getresults[7]!=-1000){
    cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;
    cout<<"Refracted and Reflected: dt(D,R)="<<(getresults[4]-getresults[5])*pow(10,9)<<" ns "<<endl;
    cout<<"*******For the Reflected Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[7]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Refracted Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[5]*pow(10,9)<<" ns"<<endl;
    cout<<"Incident Angle on Ice Surface: "<<getresults[11]<<" deg"<<endl;
  }
  
  delete []getresults;  

  double * getresults2=IceRayTracing_Cnz(x0,z0,x1,z1,1.4);

  cout<<" "<<endl;
  cout<<" For constant refractive index"<<endl;
  cout<<"Direct and Reflected: dt(D,R)="<<(getresults2[3]-getresults2[2])*pow(10,9)<<" ns"<<endl;
  cout<<"*******For the Direct Ray********"<<endl;
  cout<<"Launch Angle: "<<getresults2[0]<<" deg"<<endl;
  cout<<"Recieve Angle: "<<getresults2[4]<<" deg"<<endl;
  cout<<"Propogation Time: "<<getresults2[2]*pow(10,9)<<" ns"<<endl;
  cout<<"*******For the Reflected Ray********"<<endl;
  cout<<"Launch Angle: "<<getresults2[1]<<" deg"<<endl;
  cout<<"Recieve Angle: "<<getresults2[5]<<" deg"<<endl;
  cout<<"Propogation Time: "<<getresults2[3]*pow(10,9)<<" ns"<<endl;
  cout<<"Incident Angle in Ice on the Surface: "<<getresults2[8]<<" deg"<<endl;

  delete []getresults2; 

  auto t2b = std::chrono::high_resolution_clock::now();
  double durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();

  double Duration=durationb/1000;
  cout<<"total time taken by the script: "<<Duration<<" ms"<<endl; 
  
}
