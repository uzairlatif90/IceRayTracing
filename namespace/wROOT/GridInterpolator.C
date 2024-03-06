#include "IceRayTracing.cc"

void GridInterpolator(){
  
  double zT[1]={-100};  
  double ShowerHitDistance=50;

  int TotalAntennas=1;
  IceRayTracing::GridPositionXb.resize(TotalAntennas);
  IceRayTracing::GridPositionZb.resize(TotalAntennas);
  IceRayTracing::GridZValueb.resize(TotalAntennas);

  int AntNum=0;
      
  double xR=55;
  double zR=-10;
  int rtParameter=0; ///0 is for D ray optical time, 1 is for D ray geometric path length,  2 is for D launch angle, 3 is for D recieve angle, 4 is D for ray attenuation
                     ///5 is for R ray optical time, 6 is for R ray geometric path length,  7 is for R launch angle, 8 is for R recieve angle, 9 is R for ray attenuation
  double TimeRay[2]={0,0};
  double PathRay[2]={0,0};
  double LaunchAngle[2]={0,0};
  double RecieveAngle[2]={0,0};
  int IgnoreCh[2]={0,0};
  double IncidenceAngleInIce[2]={0,0};
  double A0=1;
  double frequency=0.1;
  double AttRay[2]={0,0};

  IceRayTracing::MakeTable(ShowerHitDistance,zR,zT[AntNum],AntNum);

  IceRayTracing::GetRayTracingSolutions(zR, xR, zT[AntNum], TimeRay, PathRay, LaunchAngle, RecieveAngle, IgnoreCh, IncidenceAngleInIce, A0, frequency, AttRay);

  vector <double> output;
  double rtresult=0,interresult=0;                                                           
  if(IgnoreCh[0]!=0){
    output.push_back(TimeRay[0]);
    output.push_back(PathRay[0]);
    output.push_back(LaunchAngle[0]);
    output.push_back(RecieveAngle[0]);
    output.push_back(AttRay[0]);
  }else{
    output.push_back(-1000);
    output.push_back(-1000);
    output.push_back(-1000);
    output.push_back(-1000);
    output.push_back(-1000);
  }
  
  if(IgnoreCh[1]!=0){
    output.push_back(TimeRay[1]);
    output.push_back(PathRay[1]);
    output.push_back(LaunchAngle[1]);
    output.push_back(RecieveAngle[1]);
    output.push_back(AttRay[1]);
  }else{
    output.push_back(-1000);
    output.push_back(-1000);
    output.push_back(-1000);
    output.push_back(-1000);
    output.push_back(-1000);
  }

  auto t1b = std::chrono::high_resolution_clock::now();                                      
  double NewZValue=IceRayTracing::GetInterpolatedValue(xR, zR, rtParameter,AntNum);
  auto t2b = std::chrono::high_resolution_clock::now();                                      
  double Duration = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b - t1b ).count();
  //cout<<"total time taken to interpolate: "<<Duration<<" nano s "<<endl;                     
  
  rtresult=output[rtParameter];
  interresult=NewZValue;
      
  interresult=NewZValue;
  cout<<"True Value ="<<rtresult<<" ,Interpolated Value ="<<interresult<<" ,True - Interpolated value "<<rtresult-interresult<<endl;
  cout<<"total time taken to interpolate: "<<Duration<<" nano s "<<endl;
  
}
