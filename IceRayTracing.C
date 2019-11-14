#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <sys/time.h>

using namespace std;

const double pi=4.0*atan(1.0); /**< Gives back value of Pi */
const double spedc=299792458.0; /**< Speed of Light in m/s */

////Set the value of the asymptotic parameter of the refractive index model
const double A_ice=1.78;

////Get the value of the B parameter for the refractive index model
double GetB(double z){
  double zabs=fabs(z);
  double B=0;

  //B=-0.375026*(1-(1.0/(1+exp(-(zabs-14.58)/8))))-0.60946*(1.0/(1+exp(-(zabs-14.58)/8)));
  //B=-0.36375*(1-(1.0/(1+exp(-(zabs-14.58)/8))))-0.60946*(1.0/(1+exp(-(zabs-14.58)/8)));
  B=-0.43;
  return B;
}

////Get the value of the C parameter for the refractive index model
double GetC(double z){
  double zabs=fabs(z);
  double C=0;
  
  //C=+0.0196219*(1-(1.0/(1+exp(-(zabs-14.58)/8))))+0.017177*(1.0/(1+exp(-(zabs-14.58)/8)));
  //C=+0.0215883*(1-(1.0/(1+exp(-(zabs-14.58)/8))))+0.017177*(1.0/(1+exp(-(zabs-14.58)/8)));
  C=0.0132;
  return C;
}

////Get the value of refractive index model for a given depth
double Getnz(double z){
  z=fabs(z);
  return A_ice+GetB(z)*exp(-GetC(z)*z);
}

////Use GSL minimiser which uses Brent's Method to find root for a given function. This will be used to find roots wherever it is needed in my code.
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
  // cout<<x_lo<<" "<<x_hi<<" "<<endl;
  // printf ("error: %s\n", gsl_strerror (status));
  
  // printf ("using %s method\n", gsl_root_fsolver_name (s));
  // printf ("%5s [%9s, %9s] %9s %9s\n","iter", "lower", "upper", "root", "err(est)");

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,0, 0.00001);
	
      if (status == GSL_SUCCESS){
	// printf ("Converged:");
	// printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi,r,x_hi - x_lo);
      }
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return r;
}

////Define the function that will be minimized to get the value of the depth of the turning point for a given refracted ray. This function basically requires the value of the L parameter to find the depth of the turning point.  This comes from the part of the fDnfR function where sqrt( n(z) - L ). This imposes the constraint then n(z)=L at the turning point of the ray from which we can find zmax.
struct Minnz_params { double a,l; };
double GetMinnz(double x,void *params){
  struct Minnz_params *p= (struct Minnz_params *) params;
  double A = p->a;
  double L = p->l;
  return A+GetB(x)*exp(-GetC(x)*x)-L;
}

////Get the value of the depth of the turning point for the refracted ray
double GetZmax(double A, double L){
  gsl_function F1;
  struct Minnz_params params1= {A,L};
  F1.function = &GetMinnz;
  F1.params = &params1;
  double zmax=FindFunctionRoot(F1,0.0,5000);
  return zmax;
}


////Analytical solution describing ray paths in ice as function of depth
struct fDnfR_params { double a, b, c, l; };
double fDnfR(double x,void *params){
  
  struct fDnfR_params *p= (struct fDnfR_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double L = p->l;
  
  return (L/C)*(1.0/sqrt(A*A-L*L))*(C*x-log(A*Getnz(x)-L*L+sqrt(A*A-L*L)*sqrt(pow(Getnz(x),2)-L*L)));;
}

////Analytical solution describing the ray path in ice as a function of the L parameter
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

////The function used to calculate ray propogation time in ice
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

////This function is minimised to find the launch angle (or the L parameter) for the direct ray
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

////This function is minimised to find the launch angle (or the L parameter) for the reflected ray
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

////This function is minimised to find the launch angle (or the L parameter) for the refracted ray
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

////(Not being used right now) This function can be minimised to calculate the angle of reciept (from the vertical) on the Rx antenna and the hit point of the ray on the ice surface given a ray incident angle. This is used to calculate at what distance a given a ray should hit the ice if we know the incoming zenith angle and the depth of the Rx antenna.
// struct fdxdz_params { double a, b, c, lang, z0,z1; };
// double fdxdz(double x,void *params){
  
//   struct fdxdz_params *p= (struct fdxdz_params *) params;
//   double A = p->a;
//   double B = p->b;
//   double C = p->c;
//   double Lang = p->lang;
//   double Z0 = p->z0;
//   double Z1 = p->z1;

//   return tan(asin( (Getnz(Z0)*sin(x))/Getnz(Z1) ) ) - tan(Lang);
// }

////function for plotting and storing the rays
void PlotAndStoreRays(double x0,double z0, double z1, double x1, double zmax, double lvalues[3], double checkzeroes[3], bool Flip){

  double lvalueD=lvalues[0];
  double lvalueR=lvalues[1];
  double lvalueRa=lvalues[2];

  double checkzeroD=checkzeroes[0];
  double checkzeroR=checkzeroes[1];
  double checkzeroRa=checkzeroes[2]; 
  
  ////Set the name of the text files
  ofstream aoutD("DirectRay.txt");
  ofstream aoutR("ReflectedRay.txt");
  ofstream aoutRa("RefractedRay.txt");

  ////Set the step size for plotting.
  double h=0.1;
  ////Set the total steps required for looping over the whole ray path
  int dmax=100000;

  ////Set the values to start the rays from
  double zn=z1;
  double xn=0;

  ////Map out the direct ray path
  int npnt=0;
  struct fDnfR_params params6a;
  struct fDnfR_params params6b;
  struct fDnfR_params params6c;
  double checknan=0;
    
  for(int i=0;i<dmax;i++){
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueD};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueD};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
    checknan=fDnfR(zn,&params6a);
    if(isnan(checknan)==false && Flip==false){
      aoutD<<npnt<<" "<<xn<<" "<<zn<<endl;;
      npnt++;
    }

    if(isnan(checknan)==false && Flip==true){
      aoutD<<npnt<<" "<<x1-xn<<" "<<zn<<endl;;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;      
    }  
  }

  ////Map out the 1st part of the reflected ray
  zn=z1;
  npnt=0;
  for(int i=0;i<dmax;i++){
    params6a = {A_ice, GetB(zn), -GetC(zn), lvalueR};
    params6b = {A_ice, GetB(z0), -GetC(z0), lvalueR};
    params6c = {A_ice, GetB(0.0000001), -GetC(0.0000001), lvalueR};
    xn=(fDnfR(-zn,&params6a)-fDnfR(-z0,&params6b)+2*fabs(fDnfR(0.0000001,&params6c)-fDnfR(-z0,&params6b)));
    checknan=fDnfR(-zn,&params6a);
    if(isnan(checknan)==false && zn<=0 && Flip==false){
      aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(isnan(checknan)==false && zn<=0 && Flip==true){
      aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
      
    zn=zn+h;
    if(zn>0){
      i=dmax+2;      
    }
  }

  ////Map out the 2nd part of the reflected ray
  zn=-0.0000001;
  for(int i=0;i<dmax;i++){  
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueR};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueR};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
    checknan=fDnfR(zn,&params6a);
    if(isnan(checknan)==false && Flip==false){
      aoutR<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }
      
    if(isnan(checknan)==false && Flip==true){
      aoutR<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }

    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;
    }
  }

  ////Map out the 1st part of the refracted ray
  zn=z1;
  npnt=0;
    
  for(int i=0;i<dmax;i++){
    params6a = {A_ice, GetB(zn), -GetC(zn), lvalueRa};
    params6b = {A_ice, GetB(z0), -GetC(z0), lvalueRa};
    params6c = {A_ice, GetB(zmax), -GetC(zmax), lvalueRa};
    xn=(fDnfR(-zn,&params6a)-fDnfR(-z0,&params6b)+2*fabs(fDnfR(zmax,&params6c)-fDnfR(-z0,&params6b)));
    checknan=fDnfR(-zn,&params6a);
    if(isnan(checknan)==false && zn<=0 && Flip==false){
      aoutRa<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(isnan(checknan)==false && zn<=0 && Flip==true){
      aoutRa<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
    
    zn=zn+h;
    if(zn>-zmax){
      i=dmax+2;      
    }
  }

  ////Map out the 2nd part of the refracted ray
  zn=-zmax;
  for(int i=0;i<dmax;i++){  
    params6a = {A_ice, GetB(zn), GetC(zn), lvalueRa};
    params6b = {A_ice, GetB(z0), GetC(z0), lvalueRa};
    xn=fDnfR(zn,&params6a)-fDnfR(z0,&params6b);
    checknan=fDnfR(zn,&params6a);
    if(isnan(checknan)==false && Flip==false){
      aoutRa<<npnt<<" "<<xn<<" "<<zn<<endl;
      npnt++;
    }

    if(isnan(checknan)==false && Flip==true){
      aoutRa<<npnt<<" "<<x1-xn<<" "<<zn<<endl;
      npnt++;
    }
    
    zn=zn-h;
    if(zn<z0){
      zn=z0;
      i=dmax+2;
    }
  }

}

////This function is used to measure the amount of time the code takes to run
typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp (){
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

double *IceRayTracing(double x0, double z0, double x1, double z1){

  ////define a pointer to give back the output of raytracing 
  double *output=new double[11];

  ////Store the ray paths in text files
  bool StoreRayPaths=false;
  ////calculate the attenuation (not included yet!)
  bool attcal=false;
  bool Flip=false;
  
  double Txcor[2]={x0,z0};//Tx positions
  double Rxcor[2]={x1,z1};//Rx Positions

  ////figure out what is the lowest depth from amongst the Rx and Tx depth. This is later used in deciding on how many steps do we need to map out the ray paths.
  double lowerz=0;
  double upz=0;
  if(Rxcor[1]<Txcor[1]){
    lowerz=Rxcor[1];
    upz=Txcor[1];
  }
  if(Rxcor[1]>Txcor[1]){
    lowerz=Txcor[1];
    upz=Rxcor[1];
  }
  if(Rxcor[1]==Txcor[1]){
    lowerz=Txcor[1];
    upz=Rxcor[1];
  }
  lowerz=lowerz;

  ////My raytracer can only work the Tx is below the Rx. If the Tx is higher than the Rx than we need to flip the depths to allow for raytracing and then we will flip them back later at the end
  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }

  ////*********This part of the code will try to get the Direct ray between Rx and Tx.***********
  ////First we setup the fDa function that will be minimised to get the launch angle (or the L parameter) for the direct ray.
  gsl_function F1;
  struct fDanfRa_params params1= {A_ice, z0, x1, z1};
  F1.function = &fDa;
  F1.params = &params1;

  ////In my raytracing solution given in the function fDnfR the launch angle (or the L parameter) has limit placed on it by this part in the solution sqrt( n(z)^2 - L^2) . This sqrt cannot be negative for both z0 and z1 and this sets the upper limit in our minimisation to get the launch angle (or the L parameter). Here I am basically setting the upper limit as GSL requires that my function is well behaved on the upper and lower bounds I give it for minimisation. 
  double UpperLimitL=Getnz(z0)*sin(90*(pi/180.0));
  if(pow(Getnz(z1),2)-pow(UpperLimitL,2)<0){
    UpperLimitL=Getnz(z1);
  }

  ////Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fDa function.
  double lvalueD=FindFunctionRoot(F1,0.0000001,UpperLimitL);
  double langD=asin(lvalueD/Getnz(z0))*(180.0/pi);
  double checkzeroD=fDa(lvalueD,&params1);

  ////Get the propagation time for the direct ray using the ftimeD function after we have gotten the value of the L parameter.
  struct ftimeD_params params2a = {A_ice, GetB(z0), -GetC(z0), spedc,lvalueD};
  struct ftimeD_params params2b = {A_ice, GetB(z1), -GetC(z1), spedc,lvalueD};

  ////we do the subtraction because we are measuring the time taken between the Tx and Rx positions
  double timeD=+ftimeD(-z0,&params2a) - ftimeD(-z1,&params2b);

  ////Setup the function that will be used to calculate the angle of reception for all the rays
  gsl_function F5;
  struct fDnfR_params params5a = {A_ice, GetB(z1), -GetC(z1), lvalueD};
  double result, abserr;
  F5.function = &fDnfR;

  ////Calculate the recieve angle for direc rays by calculating the derivative of the function at the Rx position
  F5.params = &params5a;
  gsl_deriv_central (&F5, -z1, 1e-8, &result, &abserr);
  double RangD=atan(result)*(180.0/pi);

  ////When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line.
  if(z1==z0 && isnan(RangD)==true){
    RangD=180-langD;
  }
  
  ////This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation.
  if(z1!=z0 && isnan(RangD)==true){
    RangD=90;
  }
  
  ////*********This part of the code will try to get the Reflected ray between Rx and Tx.***********
  ////First we setup the fRa function that will be minimised to get the launch angle (or the L parameter) for the reflected ray.
  gsl_function F3;
  struct fDanfRa_params params3= {A_ice, z0, x1, z1};
  F3.function = &fRa;
  F3.params = &params3;

  ////Set the upper limit for the minimisation to get the value of the launch angle (or the L parameter).  In the reflected case we set the upper limit at depth=0 m . I do not go exactly to 0 m depth since my solution blows up at the peak point of the ray. So just to be cautious I stay close to it but do not go exactly to that point.
  UpperLimitL=Getnz(0.0000001);

  ////Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fRa function.
  double lvalueR=FindFunctionRoot(F3,0.0000001,UpperLimitL);
  double langR=asin(lvalueR/Getnz(z0))*(180.0/pi);
  double checkzeroR=fRa(lvalueR,&params3); 

  ////Get the propagation time for the reflected ray using the ftimeD function after we have gotten the value of the L parameter.
  struct ftimeD_params params3a = {A_ice, GetB(z0), GetC(z0), spedc,lvalueR};
  struct ftimeD_params params3b = {A_ice, GetB(z1), GetC(z1), spedc,lvalueR};
  struct ftimeD_params params3c = {A_ice, GetB(0.0000001), GetC(0.0000001), spedc,lvalueR};

  ////we do the subtraction because we are measuring the time taken between the Tx and Rx positions. In the reflected case we basically have two direct rays 1) from Tx to surface 2) from surface to Rx.
  double timeR= 2*ftimeD(-0.0000001,&params3c) - ftimeD(z0,&params3a) - ftimeD(z1,&params3b);

  ////Also get the time for the two individual direct rays separately
  double timeR1=ftimeD(-0.0000001,&params3c) - ftimeD(z0,&params3a);
  double timeR2=ftimeD(-0.0000001,&params3c) - ftimeD(z1,&params3b);

  if(Flip==true){
    double dumR=timeR2;
    timeR2=timeR1;
    timeR1=dumR;
  }
  timeR1=timeR1;
  timeR2=timeR2;

  ////Calculate the recieve angle for reflected ray by calculating the derivative of the function at the Rx position
  struct fDnfR_params params5b = {A_ice, GetB(z1), GetC(z1), lvalueR};
  F5.params = &params5b;
  gsl_deriv_central (&F5, z1, 1e-8, &result, &abserr);
  double RangR=180-atan(result)*(180.0/pi);

  ////When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line.
  if(z1==z0 && isnan(RangR)==true){
    RangR=180-langR;
  }

  ////This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation.
  if(z1!=z0 && isnan(RangR)==true){
    RangR=90;
  }
  
  ////*********This part of the code will try to get the Refracted ray between Rx and Tx.***********
  ////Set up all the variables that will be used to get the parameters for the refracted ray
  double lvalueRa=0;
  double langRa=0;
  double checkzeroRa=1000;

  double timeRa=0;
  double timeRa1=0;
  double timeRa2=0;
  double raytime=0;
  double RangRa=0;
  
  double zmax=10;

  ////Set the upper limit for the minimisation to get the value of the launch angle (or the L parameter).  In the refracted case we set the upper limit at depth of z1 which is what we also do for the direct ray case.
  UpperLimitL=Getnz(z1);

  ////This if condition makes sure that we only try to find a refracted ray if we don't get two possible ray paths from the direct and reflected case. This saves us alot of time since we know that between each Tx and Rx position we only expect 2 rays.
  if(fabs(checkzeroR)>0.5 || fabs(checkzeroD)>0.5){
    
    ////First we setup the fRa function that will be minimised to get the launch angle (or the L parameter) for the refracted ray.
    gsl_function F4;
    struct fDanfRa_params params4= {A_ice, z0, x1, z1};
    F4.function = &fRaa;
    F4.params = &params4;

    ////Do the minimisation and get the value of the L parameter and the launch angle and then verify to see that the value of L that we got was actually a root of fRaa function. The thing to note here is the lower limit of the minimisation function is set to the L value corresponding to the reflected ray launch angle. Since we know the refracted ray always has bigger launch angle the reflected ray this reduces our range and makes the function more efficient at finding the refracted ray launch angle.
    lvalueRa=FindFunctionRoot(F4,Getnz(z0)*sin(((langR)*(pi/180.0))),UpperLimitL);
    langRa=asin(lvalueRa/Getnz(z0))*(180.0/pi);
    checkzeroRa=fRaa(lvalueRa,&params4);

    ////If the above strategy did not work then we start decreasing the reflected ray launch angle in steps of 5 degree and increase our range for minimisation to find the launch angle (or the L parameter). Sometimes the refracted and reflected rays come out to be the same in that case also I forced my solution to try harder by changing the minimisation range.
    double iangstep=5;
    while( (isnan(checkzeroRa)==true || fabs(checkzeroRa)>0.5 || fabs(lvalueRa-lvalueR)<pow(10,-9)) && langR>iangstep && iangstep<90){
      //cout<<"2nd try to get Refracted ray "<<isnan(checkzeroRa)<<" "<<fabs(checkzeroRa)<<endl;
      lvalueRa=FindFunctionRoot(F4,Getnz(z0)*sin(((langR-iangstep)*(pi/180.0))),UpperLimitL);
      langRa=asin(lvalueRa/Getnz(z0))*(180.0/pi);
      checkzeroRa=fRaa(lvalueRa,&params4);
      iangstep=iangstep+5;
    }///end the second attempt

    ////If we still did not find a refracted ray then set the check zero parameter to 1000 to make sure my code does not output this as a possible solution
    if(isnan(checkzeroRa)==true){
      checkzeroRa=1000;
    }

    ////If we did find a possible refracted ray then now we need to find the depth at which the ray turns back down without hitting the surface.
    zmax=GetZmax(A_ice,lvalueRa)+0.0000001;

    ////If the turning point depth also came out to be zero then now we are sure that there is no refracted ray
    if(zmax==0.0000001){
      checkzeroRa=1000;
    }

    ////Set parameters for ftimeD function to get the propagation time for the refracted ray
    struct ftimeD_params params4a = {A_ice, GetB(z0), GetC(z0), spedc,lvalueRa};
    struct ftimeD_params params4b = {A_ice, GetB(z1), GetC(z1), spedc,lvalueRa};
    struct ftimeD_params params4c = {A_ice, GetB(zmax), GetC(zmax), spedc,lvalueRa};

    ////This if condition checks if the function has not gone crazy and given us a turning point of the ray which is lower than both Tx and Rx and is shallower in depth than both
    if((z0<-zmax || zmax<-z1)){
      ////we do the subtraction because we are measuring the time taken between the Tx and Rx positions. In the refracted case we basically have two direct rays 1) from Tx to turning point 2) from turning point to Rx.
      raytime=2*ftimeD(zmax,&params4c) - ftimeD(z0,&params4a) - ftimeD(z1,&params4b);

      ////Also get the time for the two individual direct rays separately
      timeRa1=ftimeD(zmax,&params4c) - ftimeD(z0,&params4a);
      timeRa2=ftimeD(zmax,&params4c) - ftimeD(z1,&params4b);
      if(Flip==true){
	double dumRa=timeRa2;
	timeRa2=timeRa1;
	timeRa1=dumRa;
      }
      timeRa1=timeRa1;
      timeRa2=timeRa2;
    }
    timeRa=raytime;

    ////Calculate the recieve angle for refacted ray by calculating the derivative of the function at the Rx position
    struct fDnfR_params params5c = {A_ice, GetB(z1), GetC(z1), lvalueRa};
    F5.params = &params5c;
    gsl_deriv_central (&F5, z1, 1e-8, &result, &abserr);
    RangRa=180-atan(result)*(180.0/pi);

    ////When the Tx and Rx are at the same depth my function struggles to find a ray between them when they are very close to each other. In that case the ray is pretty much like a straight line.
    if(z1==z0 && isnan(RangRa)==true){
      RangRa=180-langRa;
    }

    ////This sometimes happens that when the Rx is very close to the peak point (or the turning point) of the ray then its hard to calculate the derivative around that area since the solution blows up around that area. therefore this is a good approximation.
    if(z1!=z0 && isnan(RangRa)==true){
      RangRa=90;
    }
    
  }///refracted if condition

  ////Fill in the output pointer after calculating all the results
  output[0]=langD;
  output[1]=langR;
  output[2]=langRa;
  output[3]=timeD;
  output[4]=timeR;
  output[5]=timeRa;
  output[6]=RangD;
  output[7]=RangR;
  output[8]=RangRa;

  ////If the flip case is true where we flipped Rx and Tx depths to trace rays then make sure everything is switched back before we give the output to the user.
  if(Flip==true){
    output[0]=180-RangD;
    output[1]=180-RangR;
    output[2]=180-RangRa;
    output[3]=timeD;
    output[4]=timeR;
    output[5]=timeRa;
    output[6]=180-langD;
    output[7]=180-langR;
    output[8]=180-langRa;
  }

  ////This part of the code can be used if the user wants to plot the individual ray paths. This part of the code prints out the individual ray paths in text files.
  if(StoreRayPaths==true){
    double lvalues[3];
    lvalues[0]=lvalueD;
    lvalues[1]=lvalueR;
    lvalues[2]=lvalueRa;

    double checkzeroes[3];
    checkzeroes[0]=checkzeroD;
    checkzeroes[1]=checkzeroR;
    checkzeroes[2]=checkzeroRa;
    
    PlotAndStoreRays(x0,z0,z1,x1,zmax,lvalues,checkzeroes,Flip);
  }

  dsw=0;
  ////If the Tx and Rx depth were switched then put them back to their original position
  if(Flip==true){
    dsw=z0;
    z0=z1;
    z1=dsw;
  }  
  
  ////print out all the output from the code
  //cout<<0<<" ,x0= "<<x0<<" ,z0= "<<z0<<" ,x1= "<<x1<<" ,z1= "<<z1<<" ,langRa= "<<output[2]<<" ,langR= "<<output[1]<<" ,langD= "<<output[0]<<" ,langD-langR= "<<output[0]-output[1]<<" ,langD-langRa= "<<output[0]-output[2]<<" ,RangRa= "<<output[8]<<" ,RangR= "<<output[7]<<" ,RangD= "<<output[6]<<" ,RangR-RangD= "<<output[7]-output[6]<<" ,RangRa-RangD= "<<output[8]-output[6]<<" ,timeRa= "<<output[5]<<" ,timeR= "<<output[4]<<" ,timeD= "<<output[3]<<" ,timeR-timeD= "<<output[4]-output[3]<<" ,timeRa-timeD= "<<output[5]-output[3]<<" ,lvalueRa "<<lvalueRa<<" ,lvalueR "<<lvalueR<<" "<<" ,lvalueD "<<lvalueD<<" ,checkzeroRa "<<checkzeroRa<<" ,checkzeroR "<<checkzeroR<<" ,checkzeroD "<<checkzeroD<<endl;

  ////fill in the output array part where you fill in the times for the two parts of the reflected or refracted rays
  if(fabs(checkzeroR)<0.5){
    output[9]=timeR1;
    output[10]=timeR2;
  }
  
  if(fabs(checkzeroRa)<0.5){
    output[9]=timeRa1;
    output[10]=timeRa2;
  }

  ////Set the recieve angle to be zero for a ray which did not give us a possible path between Tx and Rx. I use this as a flag to determine which two rays gave me possible ray paths.
  if(fabs(checkzeroD)>0.5){
    output[6]=0;
  }
  if(fabs(checkzeroR)>0.5){
    output[7]=0;
  }
  if(fabs(checkzeroRa)>0.5){
    output[8]=0;
  }
  
  return output;
}

int main(int argc, char ** argv){

  if(argc<7){
    cout<<"More Arguments needed!"<<endl;
    cout<<"Example run command: ./IceRayTracing 0 0 -100 0 100 -5"<<endl;
    cout<<"Here the first three arguments are Tx coordinates i.e. x_Tx=0 m, y_Tx=0 m, z_Tx=-100 m  and the second three arguments are Rx coordinates i.e. x_Rx=0 m, y_Rx=100 m, z_Rx=-5 m "<<endl;
    return 0;
  }
  if(argc==7){
    cout<<"Tx set at X="<<atof(argv[1])<<" m, Y="<<atof(argv[2])<<" m, Z="<<atof(argv[3])<<" m"<<endl;
    cout<<"Rx set at X="<<atof(argv[4])<<" m, Y="<<atof(argv[5])<<" m, Z="<<atof(argv[6])<<" m"<<endl;
  } 
  if(argc>7){
    cout<<"More Arguments than needed!"<<endl;
    cout<<"Example run command: ./IceRayTracing 0 0 -100 0 100 -5"<<endl;
    cout<<"Here the first three arguments are Tx coordinates i.e. x_Tx=0 m, y_Tx=0 m, z_Tx=-100 m  and the second three arguments are Rx coordinates i.e. x_Rx=0 m, y_Rx=100 m, z_Rx=-5 m "<<endl;
    return 0;
  }

  double TxCor[3]={atof(argv[1]),atof(argv[2]),atof(argv[3])};
  double RxCor[3]={atof(argv[4]),atof(argv[5]),atof(argv[6])};
  
  ////An example for Direct and Reflected Rays
  // double TxCor[3]={0,0,-100};
  // double RxCor[3]={0,100,-5};

  ////An example for Direct and Refracted Rays
  // double TxCor[3]={0,0,-67.5489};
  // double RxCor[3]={0,1338.3,-800};
 
  ////For recording how much time the process took
  timestamp_t t0 = get_timestamp();
  
  double x0=0;
  double z0=TxCor[2];
  double x1=sqrt(pow(TxCor[0]-RxCor[0],2)+pow(TxCor[1]-RxCor[1],2));
  double z1=RxCor[2];

  double * getresults=IceRayTracing(x0,z0,x1,z1);

  cout<<" "<<endl;
  cout<<"x0="<<x0<<" m , z0="<<z0<<" m ,x1="<<x1<<" m ,z1="<<z1<<" m "<<endl;
  if(getresults[6]!=0 && getresults[7]!=0 && getresults[8]!=0){
    cout<<"No Possible Ray Paths between Tx and Rx!!! :("<<endl;
  }
  
  if(getresults[6]!=0 && getresults[8]!=0){ 
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
      
  if(getresults[6]!=0 && getresults[7]!=0){
    cout<<"Direct and Reflected: dt(D,R)="<<(getresults[4]-getresults[3])*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Direct Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[6]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[3]*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Reflected Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[7]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;   
  }
      
  if(getresults[8]!=0 && getresults[7]!=0 && getresults[6]==0){
    cout<<"Refracted and Reflected: dt(D,R)="<<(getresults[4]-getresults[5])*pow(10,9)<<" ns "<<endl;
    cout<<"*******For the Reflected Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[7]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;
    cout<<"*******For the Refracted Ray********"<<endl;
    cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
    cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
    cout<<"Propogation Time: "<<getresults[5]*pow(10,9)<<" ns"<<endl;
  }
  
  delete []getresults;  

  timestamp_t t1 = get_timestamp();
  double secs = (t1 - t0) / 1000000.0L;
  cout<<"total time taken by the script: "<<secs<<" s"<<endl;

  return 0;
}
