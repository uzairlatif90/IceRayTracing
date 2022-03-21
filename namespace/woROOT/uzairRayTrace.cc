#include "IceRayTracing.hh"

/* This is an example script to call and play with Uzair's Ray Tracer */

int main(int argc, char**argv){
  
  if(argc<7){
    cout<<"More Arguments needed!"<<endl;
    cout<<"Example run command: ./uzairRayTrace 0 0 -100 0 100 -5"<<endl;
    cout<<"Here the first three arguments are Tx coordinates i.e. x_Tx=0 m, y_Tx=0 m, z_Tx=-100 m  and the second three arguments are Rx coordinates i.e. x_Rx=0 m, y_Rx=100 m, z_Rx=-5 m "<<endl;
    return 0;
  }
  if(argc==7){
    cout<<"Tx set at X="<<atof(argv[1])<<" m, Y="<<atof(argv[2])<<" m, Z="<<atof(argv[3])<<" m"<<endl;
    cout<<"Rx set at X="<<atof(argv[4])<<" m, Y="<<atof(argv[5])<<" m, Z="<<atof(argv[6])<<" m"<<endl;
  } 
  if(argc>7){
    cout<<"More Arguments than needed!"<<endl;
    cout<<"Example run command: ./uzairRayTrace 0 0 -100 0 100 -5"<<endl;
    cout<<"Here the first three arguments are Tx coordinates i.e. x_Tx=0 m, y_Tx=0 m, z_Tx=-100 m  and the second three arguments are Rx coordinates i.e. x_Rx=0 m, y_Rx=100 m, z_Rx=-5 m "<<endl;
    return 0;
  }

  double TxCor[3]={atof(argv[1]),atof(argv[2]),atof(argv[3])};
  double RxCor[3]={atof(argv[4]),atof(argv[5]),atof(argv[6])};

  
  ////For recording how much time the process took
  auto t1b = std::chrono::high_resolution_clock::now();  
  
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

  double * getresults=IceRayTracing::IceRayTracing(x0,z0,x1,z1);

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

  delete []getresults;
  
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

  return 0;
  
}

