#include <stdio.h>
#include <math.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_spline.h>
#include <sys/time.h>

using namespace std;

const double pi=4.0*atan(1.0); /**< Gives back value of Pi */
const double spedc=299792458.0; /**< Speed of Light in m/s */

const double A_ice=1.78;

double GetB(double z){
  double zabs=fabs(z);
  double B=0;

  //B=-0.375026*(1-(1.0/(1+exp(-(zabs-14.58)/8))))-0.60946*(1.0/(1+exp(-(zabs-14.58)/8)));
  //B=-0.36375*(1-(1.0/(1+exp(-(zabs-14.58)/8))))-0.60946*(1.0/(1+exp(-(zabs-14.58)/8)));
  B=-0.43;
  return B;
}

double GetC(double z){
  double zabs=fabs(z);
  double C=0;
  
  //C=+0.0196219*(1-(1.0/(1+exp(-(zabs-14.58)/8))))+0.017177*(1.0/(1+exp(-(zabs-14.58)/8)));
  //C=+0.0215883*(1-(1.0/(1+exp(-(zabs-14.58)/8))))+0.017177*(1.0/(1+exp(-(zabs-14.58)/8)));
  C=0.0132;
  return C;
}

double Getnz(double z){
  z=fabs(z);
  return A_ice+GetB(z)*exp(-GetC(z)*z);
}

////Use GSL minimiser which uses Brent's Method to find root for a given function
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

struct Minnz_params { double a,l; };
double GetMinnz(double x,void *params){
  struct Minnz_params *p= (struct Minnz_params *) params;
  double A = p->a;
  double L = p->l;
  return A+GetB(x)*exp(-GetC(x)*x)-L;
}

double GetZmax(double A, double L){
  gsl_function F1;
  struct Minnz_params params1= {A,L};
  F1.function = &GetMinnz;
  F1.params = &params1;
  double zmax=FindFunctionRoot(F1,0.0,5000);
  return zmax;
}


////Analytical solution describing the ray path in ice
struct fDnfR_params { double a, b, c, l; };
double fDnfR(double x,void *params){
  
  struct fDnfR_params *p= (struct fDnfR_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double L = p->l;
  
  return (L/C)*(1.0/sqrt(A*A-L*L))*(C*x-log(A*Getnz(x)-L*L+sqrt(A*A-L*L)*sqrt(pow(Getnz(x),2)-L*L)));;
}

////Analytical solution describing the ray path in ice
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

////Define the function that will be minimised to calculate the angle of reciept (from the vertical) on the antenna and the hit point of the ray on the ice surface given a ray incident angle
struct fdxdz_params { double a, b, c, lang, z0,z1; };
double fdxdz(double x,void *params){
  
  struct fdxdz_params *p= (struct fdxdz_params *) params;
  double A = p->a;
  double B = p->b;
  double C = p->c;
  double Lang = p->lang;
  double Z0 = p->z0;
  double Z1 = p->z1;

  return tan(asin( (Getnz(Z0)*sin(x))/Getnz(Z1) ) ) - tan(Lang);
}

////This function is used to measure the amount of time the code takes
typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp (){
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

double *IceRayTracing(double x0, double z0, double x1, double z1){
  
  double *output=new double[11];

  ////Store the ray paths in text files
  bool StoreRayPaths=false;
  ////calculate the attenuation (not included yet!)
  bool attcal=false;
  bool Flip=false;
  
  double Txcor[2]={x0,z0};//Tx positions
  double Rxcor[2]={x1,z1};//Rx Positions

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

  double dsw=z0;
  if(z0>z1){
    z0=z1;
    z1=dsw;
    Flip=true;
  }

  gsl_function F1;
  struct fDanfRa_params params1= {A_ice, z0, x1, z1};
  F1.function = &fDa;
  F1.params = &params1;

  double UpperLimitL=Getnz(z0)*sin(90*(pi/180.0));
  if(pow(Getnz(z1),2)-pow(UpperLimitL,2)<0){
    UpperLimitL=Getnz(z1);
  }

  double lvalueD=FindFunctionRoot(F1,0.0000001,UpperLimitL);
  double langD=asin(lvalueD/Getnz(z0))*(180.0/pi);
  double checkzeroD=fDa(lvalueD,&params1);
  
  struct ftimeD_params params2a = {A_ice, GetB(z0), -GetC(z0), spedc,lvalueD};
  struct ftimeD_params params2b = {A_ice, GetB(z1), -GetC(z1), spedc,lvalueD};
  double timeD=+ftimeD(-z0,&params2a) - ftimeD(-z1,&params2b);
  
  gsl_function F3;
  struct fDanfRa_params params3= {A_ice, z0, x1, z1};
  F3.function = &fRa;
  F3.params = &params3;

  UpperLimitL=Getnz(0.0000001);
  
  double lvalueR=FindFunctionRoot(F3,0.0000001,UpperLimitL);
  double langR=asin(lvalueR/Getnz(z0))*(180.0/pi);
  double checkzeroR=fRa(lvalueR,&params3); 
  
  struct ftimeD_params params3a = {A_ice, GetB(z0), GetC(z0), spedc,lvalueR};
  struct ftimeD_params params3b = {A_ice, GetB(z1), GetC(z1), spedc,lvalueR};
  struct ftimeD_params params3c = {A_ice, GetB(0.0000001), GetC(0.0000001), spedc,lvalueR};
  
  double timeR= 2*ftimeD(-0.0000001,&params3c) - ftimeD(z0,&params3a) - ftimeD(z1,&params3b);
  double timeR1=ftimeD(-0.0000001,&params3c) - ftimeD(z0,&params3a);
  double timeR2=ftimeD(-0.0000001,&params3c) - ftimeD(z1,&params3b);
  
  if(Flip==true){
    double dumR=timeR2;
    timeR2=timeR1;
    timeR1=dumR;
  }
  timeR1=timeR1*spedc;
  timeR2=timeR2*spedc;

  double lvalueRa=0;
  double langRa=0;
  double checkzeroRa=1000;

  double timeRa=0;
  double timeRa1=0;
  double timeRa2=0;
  double raytime=0;

  double zmax=10;

  UpperLimitL=Getnz(z1);
  
  if(fabs(checkzeroR)>0.5 || fabs(checkzeroD)>0.5){
    gsl_function F4;
    struct fDanfRa_params params4= {A_ice, z0, x1, z1};
    F4.function = &fRaa;
    F4.params = &params4;
    lvalueRa=FindFunctionRoot(F4,Getnz(z0)*sin(((langR)*(pi/180.0))),UpperLimitL);
    langRa=asin(lvalueRa/Getnz(z0))*(180.0/pi);
    checkzeroRa=fRaa(lvalueRa,&params4);
    
    double iangstep=5;
    while( (isnan(checkzeroRa)==true || fabs(checkzeroRa)>0.5 || fabs(lvalueRa-lvalueR)<pow(10,-9)) && langR>iangstep && iangstep<90){
      //cout<<"2nd try to get Refracted ray "<<isnan(checkzeroRa)<<" "<<fabs(checkzeroRa)<<endl;
      lvalueRa=FindFunctionRoot(F4,Getnz(z0)*sin(((langR-iangstep)*(pi/180.0))),UpperLimitL);
      langRa=asin(lvalueRa/Getnz(z0))*(180.0/pi);
      checkzeroRa=fRaa(lvalueRa,&params4);
      iangstep=iangstep+5;
    }///end the second attempt    
    
    if(isnan(checkzeroRa)==true){
      checkzeroRa=1000;
    }
    
    zmax=GetZmax(A_ice,lvalueRa)+0.0000001;
    //cout<<"zmax is "<<zmax<<endl;

    if(zmax==0.0000001){
      checkzeroRa=1000;
    }

    struct ftimeD_params params4a = {A_ice, GetB(z0), GetC(z0), spedc,lvalueRa};
    struct ftimeD_params params4b = {A_ice, GetB(z1), GetC(z1), spedc,lvalueRa};
    struct ftimeD_params params4c = {A_ice, GetB(zmax), GetC(zmax), spedc,lvalueRa};
    
    if((z0<-zmax || zmax<-z1)){
      raytime=2*ftimeD(zmax,&params4c) - ftimeD(z0,&params4a) - ftimeD(z1,&params4b);
      timeRa1=ftimeD(zmax,&params4c) - ftimeD(z0,&params4a);
      timeRa2=ftimeD(zmax,&params4c) - ftimeD(z1,&params4b);
      if(Flip==true){
	double dumRa=timeRa2;
	timeRa2=timeRa1;
	timeRa1=dumRa;
      }
      timeRa1=timeRa1*spedc;
      timeRa2=timeRa2*spedc;
    }
    timeRa=raytime;
  }///refracted if condition
  
  //calculate the recive angles
  gsl_function F5;
  struct fDnfR_params params5a = {A_ice, GetB(z1), -GetC(z1), lvalueD};
  struct fDnfR_params params5b = {A_ice, GetB(z1), GetC(z1), lvalueR};
  struct fDnfR_params params5c = {A_ice, GetB(z1), GetC(z1), lvalueRa};

  double result, abserr;
  F5.function = &fDnfR;
  F5.params = &params5a;
  gsl_deriv_central (&F5, -z1, 1e-8, &result, &abserr);
  double RangD=atan(result)*(180.0/pi);
  F5.params = &params5b;
  gsl_deriv_central (&F5, z1, 1e-8, &result, &abserr);
  double RangR=180-atan(result)*(180.0/pi);
  F5.params = &params5c;
  gsl_deriv_central (&F5, z1, 1e-8, &result, &abserr);
  double RangRa=180-atan(result)*(180.0/pi);

  if(z1==z0 && isnan(RangRa)==true){
    RangRa=180-langRa;
  }
  if(z1==z0 && isnan(RangR)==true){
    RangR=180-langR;
  }
  if(z1==z0 && isnan(RangD)==true){
    RangD=180-langD;
  }

  if(z1!=z0 && isnan(RangRa)==true){
    RangRa=90;
  }
  if(z1!=z0 && isnan(RangR)==true){
    RangR=90;
  }
  if(z1!=z0 && isnan(RangD)==true){
    RangD=90;
  }
  
  output[0]=langD;
  output[1]=langR;
  output[2]=langRa;
  output[3]=timeD;
  output[4]=timeR;
  output[5]=timeRa;
  output[6]=RangD;
  output[7]=RangR;
  output[8]=RangRa;
  
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

  if(StoreRayPaths==true){

    ofstream aoutD("DirectRay.txt");
    ofstream aoutR("ReflectedRay.txt");
    ofstream aoutRa("RefractedRay.txt");
    
    double h=0.1;
    int dmax=-round(lowerz/h);
    dmax=100000;
    double zn=z1;
    double xn=0;
    
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

    if(Flip==true){
      dsw=z0;
      z0=z1;
      z1=dsw;
    }  
  }
  
  //cout<<0<<" ,x0= "<<x0<<" ,z0= "<<z0<<" ,x1= "<<x1<<" ,z1= "<<z1<<" ,langRa= "<<output[2]<<" ,langR= "<<output[1]<<" ,langD= "<<output[0]<<" ,langD-langR= "<<output[0]-output[1]<<" ,langD-langRa= "<<output[0]-output[2]<<" ,RangRa= "<<output[8]<<" ,RangR= "<<output[7]<<" ,RangD= "<<output[6]<<" ,RangR-RangD= "<<output[7]-output[6]<<" ,RangRa-RangD= "<<output[8]-output[6]<<" ,timeRa= "<<output[5]<<" ,timeR= "<<output[4]<<" ,timeD= "<<output[3]<<" ,timeR-timeD= "<<output[4]-output[3]<<" ,timeRa-timeD= "<<output[5]-output[3]<<" ,lvalueRa "<<lvalueRa<<" ,lvalueR "<<lvalueR<<" "<<" ,lvalueD "<<lvalueD<<" ,checkzeroRa "<<checkzeroRa<<" ,checkzeroR "<<checkzeroR<<" ,checkzeroD "<<checkzeroD<<endl;

  if(fabs(checkzeroR)<0.5){
    output[9]=timeR1;
    output[10]=timeR2;
  }
  
  if(fabs(checkzeroRa)<0.5){
    output[9]=timeRa1;
    output[10]=timeRa2;
  }
  
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
  if(getresults[6]!=0 && getresults[8]!=0){ 
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
      
  if(getresults[6]!=0 && getresults[7]!=0){
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
  }
      
  if(getresults[8]!=0 && getresults[7]!=0 && getresults[6]==0){
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
  }
  
  delete []getresults;  

  timestamp_t t1 = get_timestamp();
  double secs = (t1 - t0) / 1000000.0L;
  cout<<"total time taken by the script: "<<secs<<" s"<<endl;

  return 0;
}
