#include "AirToIceRayTracing.cc"

void TraceAirToIce(){
  
  ////All variables are in m here
  double AntennaDepth=-100;////Depth of antenna in the ice
  double IceLayerHeight=3000;////Height where the ice layer starts off
  double AntennaNumber=0;
  double AirTxHeight=3000.1;
  double HorizontalDistance=20;
  
  double opticalPathLengthInIce;
  double opticalPathLengthInAir;
  double geometricalPathLengthInIce;
  double geometricalPathLengthInAir;
  double launchAngle;
  double horidist2interpnt;
  double transmissionCoefficientS;
  double transmissionCoefficientP;
  double RecievedAngleInIce;
  
  bool CheckSol=false;////check if solution exists or not
  
  AirToIceRayTracing::MakeAtmosphere();

  AirToIceRayTracing::A_const=AirToIceRayTracing::Getnz_air(IceLayerHeight);
  AirToIceRayTracing::UseConstantRefractiveIndex=true;
  AirToIceRayTracing::A_air=AirToIceRayTracing::A_const;
  
  CheckSol=AirToIceRayTracing::GetRayTracingSolution(AirTxHeight*100, HorizontalDistance*100 ,AntennaDepth*100, IceLayerHeight*100,opticalPathLengthInIce, opticalPathLengthInAir, geometricalPathLengthInIce, geometricalPathLengthInAir, launchAngle, horidist2interpnt,transmissionCoefficientS,transmissionCoefficientP,RecievedAngleInIce);

  if(CheckSol==true){
    cout<<" We have a solution!!!"<<endl;
    cout<<"AntennaNumber: "<< AntennaNumber<<endl;
    cout<<"AirTxHeight: "<<AirTxHeight<<endl;
    cout<<"HorizontalDistance: "<<HorizontalDistance<<endl;
    cout<<"geometricalPathLengthInIce: "<<geometricalPathLengthInIce/100<<endl;
    cout<<"geometricalPathLengthInAir: "<<geometricalPathLengthInAir/100<<endl;
    cout<<"launchAngle: "<<launchAngle*(180/AirToIceRayTracing::pi)<<endl;
    cout<<"horidist2interpnt: "<<horidist2interpnt/100<<endl;
    cout<<"transmissionCoefficientS: "<<transmissionCoefficientS<<endl;
    cout<<"transmissionCoefficientP: "<<transmissionCoefficientP<<endl;
    cout<<"RecievedAngleInIce: "<<RecievedAngleInIce<<endl;
  }else{
    cout<<" We do NOT have a solution!!!"<<endl;
  }
  
}
