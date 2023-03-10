#include "../SimonCR/AirToIceRayTracing.cc"

void GetFocusingFactorAir(double zT, double xR, double zR, double &focusing, double focusparts[3]){

  double A0=1;
  double frequency=0.1;//Tx frequency in GHz  

  double distance[2]={0,0};
  double recAng[2]={0,0};
  double lauAng[2]={0,0};
  double recPos[2]={0,0};
  double nTx= AirToIceRayTracing::Getnz_ice(zT);  // emitter
  double nRx= AirToIceRayTracing::Getnz_air(zR); // receiver   
 
  double opticalPathLengthInIce;
  double opticalPathLengthInAir;
  double geometricalPathLengthInIce;
  double geometricalPathLengthInAir;
  double launchAngle;
  double horidist2interpnt;
  double transmissionCoefficientS;
  double transmissionCoefficientP;
  double RecievedAngleInIce;
  double AngleOfIncidenceOnIce;
  
  ////All variables are in m here
  double AntennaDepth=zT;////Depth of antenna in the ice
  double IceLayerHeight=3000;////Height where the ice layer starts off
  double AirTxHeight=3000+zR;
  double HorizontalDistance=xR; 
  bool CheckSol1=false;////check if solution exists or not
  
  //cout<<"check parameters "<<AntennaDepth<<" "<<IceLayerHeight<<" "<<AirTxHeight<<" "<<HorizontalDistance<<endl;
    
  CheckSol1=AirToIceRayTracing::GetRayTracingSolution(AirTxHeight*100, HorizontalDistance*100 ,AntennaDepth*100, IceLayerHeight*100,opticalPathLengthInIce, opticalPathLengthInAir, geometricalPathLengthInIce, geometricalPathLengthInAir, launchAngle, horidist2interpnt,AngleOfIncidenceOnIce,RecievedAngleInIce);
  
  distance[0]=(geometricalPathLengthInIce+ geometricalPathLengthInAir)/100;
  recAng[0]=RecievedAngleInIce*(AirToIceRayTracing::pi/180);
  lauAng[0]=launchAngle*(AirToIceRayTracing::pi/180);
  recPos[0]=zR;

  AirTxHeight=3000+zR+0.01;
  
  bool CheckSol2=false;////check if solution exists or not

  CheckSol2=AirToIceRayTracing::GetRayTracingSolution(AirTxHeight*100, HorizontalDistance*100 ,AntennaDepth*100, IceLayerHeight*100,opticalPathLengthInIce, opticalPathLengthInAir, geometricalPathLengthInIce, geometricalPathLengthInAir, launchAngle, horidist2interpnt,AngleOfIncidenceOnIce,RecievedAngleInIce);

  distance[1]=(geometricalPathLengthInIce+ geometricalPathLengthInAir)/100;
  recAng[1]=RecievedAngleInIce*(AirToIceRayTracing::pi/180);
  lauAng[1]=launchAngle*(AirToIceRayTracing::pi/180);
  recPos[1]=zR+0.01;  
  
  if(CheckSol1!=false && CheckSol2!=false){
    focusing= sqrt( ((distance[0] / (sin(recAng[0]) * fabs( (recPos[1] - recPos[0]) / (lauAng[1] - lauAng[0]) ) ) ) * (nTx / nRx) ));
    focusparts[0]=sqrt(distance[0] / sin(recAng[0]));
    focusparts[1]=sqrt(1.0/fabs( (recPos[1] - recPos[0]) / (lauAng[1] - lauAng[0]) ));
    focusparts[2]=sqrt((nTx / nRx));
    //cout<<"we are here "<<focusing<<" "<<focusparts[0]<<" "<<focusparts[1]<<" "<<focusparts[2]<<endl;
  }
  
}


void FocusingFactorAir(){

  AirToIceRayTracing::MakeAtmosphere("../SimonCR/Atmosphere.dat");

  // AirToIceRayTracing::A_const=AirToIceRayTracing::Getnz_air(3000);
  // AirToIceRayTracing::UseConstantRefractiveIndex=true;
  // AirToIceRayTracing::A_air=AirToIceRayTracing::A_const;
  
  double xR=0;
  double zR=0;
  double zT=-50;
  
  double GridStartX=0.5;
  double GridStopX=100.5;

  double GridStartZ=+50.5-(100/2);
  double GridStopZ=+50.5+(100/2);

  double GridStepSizeX_O=0.5;
  double GridStepSizeZ_O=0.5;
  
  // int TotalStepsX_O=(100/GridStepSizeX_O)+1;
  // int TotalStepsZ_O=(100/GridStepSizeZ_O)+1;

  int TotalStepsX_O=(100/GridStepSizeX_O)+1;
  int TotalStepsZ_O=(100/GridStepSizeZ_O)+1;

  TGraph2D *gr2A=new TGraph2D();
  TGraph2D *gr2B=new TGraph2D();
  TGraph2D *gr2C=new TGraph2D();
  TGraph2D *gr2D=new TGraph2D();

  TGraph2D *gr3A=new TGraph2D();
  TGraph2D *gr3B=new TGraph2D();
  TGraph2D *gr3C=new TGraph2D();
  TGraph2D *gr3D=new TGraph2D();


  gr2A->SetTitle("Focusing Factor;Distance (m); Depth (m)");
  gr2B->SetTitle("#sqrt{R / sin(#theta_{recieve})};Distance (m); Depth (m)");
  gr2C->SetTitle("#sqrt{|#Delta#theta_{launch}/#Deltaz |};Distance (m); Depth (m)");
  gr2D->SetTitle("#sqrt{nTx/nRx};Distance (m); Depth (m)");
  
  int count=0;
  int count2=0;
  
  for(int ix=0;ix<TotalStepsX_O;ix++){
    for(int iz=0;iz<TotalStepsZ_O;iz++){

      double xR=GridStartX+GridStepSizeX_O*ix;
      double zR=GridStartZ+GridStepSizeZ_O*iz;

      double focusing=1;
      double focusparts[3]={1,1,1};
      GetFocusingFactorAir(zT, xR, zR, focusing, focusparts);
      if(std::isnan(focusing)==false && focusing!=1){
	gr2A->SetPoint(count, xR,zR, focusing);
	gr2B->SetPoint(count, xR,zR, focusparts[0]);
	gr2C->SetPoint(count, xR,zR, focusparts[1]);
	gr2D->SetPoint(count, xR,zR, focusparts[2]);
	count++;
      }
      
    }
  }
 
  TCanvas *c1=new TCanvas("c1","c1");
  c1->Divide(2,2);
  c1->cd(1);
  gr2A->Draw("cont4z");
  c1->cd(2);
  gr2B->Draw("cont4z");
  c1->cd(3);
  gr2C->Draw("cont4z");
  c1->cd(4);
  gr2D->Draw("cont4z");

  
}
