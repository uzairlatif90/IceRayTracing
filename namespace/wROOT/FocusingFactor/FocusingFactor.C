#include "../SimonCR/IceRayTracing.cc"

void GetFocusingFactor(double zT, double xR, double zR, double focusing[2], double focusparts[6]){

  double A0=1;
  double frequency=0.1;//Tx frequency in GHz  

  double distance[4]={0,0,0,0};
  double recAng[4]={0,0,0,0};
  double lauAng[4]={0,0,0,0};
  double recPos[4]={0,0,0,0};
  double nTx= IceRayTracing::Getnz(zT);  // emitter
  double nRx= IceRayTracing::Getnz(zR); // receiver  
  
  double TimeRay[2]={0,0};
  double PathRay[2]={0,0};
  double LaunchAngle[2]={0,0};
  double RecieveAngle[2]={0,0};
  int IgnoreCh[2]={0,0};
  double IncidenceAngleInIce[2]={0,0};
  double AttRay[2]={0,0};
  IceRayTracing::GetRayTracingSolutions(zR, xR, zT, TimeRay, PathRay, LaunchAngle, RecieveAngle, IgnoreCh, IncidenceAngleInIce, A0, frequency, AttRay);
  
  distance[0]=PathRay[0];
  recAng[0]=RecieveAngle[0]*(IceRayTracing::pi/180);
  lauAng[0]=LaunchAngle[0]*(IceRayTracing::pi/180);
  recPos[0]=zR;

  distance[1]=PathRay[1];
  recAng[1]=RecieveAngle[1]*(IceRayTracing::pi/180);
  lauAng[1]=LaunchAngle[1]*(IceRayTracing::pi/180);
  recPos[1]=zR;
  
  double TimeRayB[2]={0,0};
  double PathRayB[2]={0,0};
  double LaunchAngleB[2]={0,0};
  double RecieveAngleB[2]={0,0};
  int IgnoreChB[2]={0,0};
  double IncidenceAngleInIceB[2]={0,0};
  double AttRayB[2]={0,0};
  
  IceRayTracing::GetRayTracingSolutions(zR-0.01, xR, zT, TimeRayB, PathRayB, LaunchAngleB, RecieveAngleB, IgnoreChB, IncidenceAngleInIceB, A0, frequency, AttRayB);
  
  distance[2]=PathRayB[0];
  recAng[2]=RecieveAngleB[0]*(IceRayTracing::pi/180);
  lauAng[2]=LaunchAngleB[0]*(IceRayTracing::pi/180);
  recPos[2]=zR-0.01;

  distance[3]=PathRayB[1];
  recAng[3]=RecieveAngleB[1]*(IceRayTracing::pi/180);
  lauAng[3]=LaunchAngleB[1]*(IceRayTracing::pi/180);
  recPos[3]=zR-0.01;

  // cout<<" 1 "<<distance[0]<<" "<<recAng[0]*(180./IceRayTracing::pi)<<" "<<lauAng[0]*(180./IceRayTracing::pi)<<endl;
  // cout<<" 2 "<<distance[1]<<" "<<recAng[1]*(180./IceRayTracing::pi)<<" "<<lauAng[1]*(180./IceRayTracing::pi)<<endl;
  // cout<<" 3 "<<distance[2]<<" "<<recAng[2]<<" "<<lauAng[2]<<endl;
  // cout<<" 4 "<<distance[3]<<" "<<recAng[3]<<" "<<lauAng[3]<<endl;
  
  if(RecieveAngle[0]!=-1000 && RecieveAngleB[0]!=-1000){
    focusing[0]= sqrt( ((distance[0] / (sin(recAng[0]) * fabs( (recPos[2] - recPos[0]) / (lauAng[2] - lauAng[0]) ) ) ) * (nTx / nRx) ));
    focusparts[0]=sqrt(distance[0] / sin(recAng[0]));
    focusparts[1]=sqrt(1.0/fabs( (recPos[2] - recPos[0]) / (lauAng[2] - lauAng[0]) ));
    focusparts[2]=sqrt((nTx / nRx));
  }
  
  if(RecieveAngle[1]!=-1000 && RecieveAngleB[1]!=-1000){
    focusing[1]= sqrt( ((distance[1] /( sin(recAng[1]) * fabs( (recPos[3] - recPos[1]) / (lauAng[3] - lauAng[1]) ) ) ) * (nTx / nRx)));
    focusparts[3]=sqrt(distance[1] / sin(recAng[1]));
    focusparts[4]=sqrt(1.0/fabs( (recPos[3] - recPos[1]) / (lauAng[3] - lauAng[1]) ));
    focusparts[5]=sqrt((nTx / nRx));
  }
  
  //cout<<"focusing factor is "<<focusing[0]<<" "<<focusing[1]<<endl;
}

void FocusingFactor(){

  
  double GridStartX=300.5;
  double GridStopX=500.5;

  double GridStartZ=-150.1;//-50.1-(100/2);
  double GridStopZ=-50.1;//-50.1+(100/2);

  double GridStepSizeX_O=0.5;
  double GridStepSizeZ_O=0.5; 
 
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

  gr3A->SetTitle("Focusing Factor;Distance (m); Depth (m)");
  gr3B->SetTitle("#sqrt{R / sin(#theta_{recieve})};Distance (m); Depth (m)");
  gr3C->SetTitle("#sqrt{|#Delta#theta_{launch}/#Deltaz |};Distance (m); Depth (m)");
  gr3D->SetTitle("#sqrt{nTx/nRx};Distance (m); Depth (m)");
  
  int count=0;
  int count2=0;

  double xT=0;
  double zT=0;
  double zR=-50; 
  
  for(int ix=0;ix<TotalStepsX_O;ix++){
    for(int iz=0;iz<TotalStepsZ_O;iz++){

      double xT=GridStartX+GridStepSizeX_O*ix;
      double zT=GridStartZ+GridStepSizeZ_O*iz;

      double focusing[2]={1,1};
      double focusparts[6]={1,1,1,1,1,1};
  
      GetFocusingFactor(zT, xT, zR, focusing, focusparts);
      if(std::isnan(focusing[0])==false && focusing[0]!=1){
	gr2A->SetPoint(count, xT,zT, focusing[0]);
	gr2B->SetPoint(count, xT,zT, focusparts[0]);
	gr2C->SetPoint(count, xT,zT, focusparts[1]);
	gr2D->SetPoint(count, xT,zT, focusparts[2]);
	count++;
      }
      if(std::isnan(focusing[1])==false && focusing[1]!=1){
	gr3A->SetPoint(count2, xT,zT, focusing[1]);
	gr3B->SetPoint(count2, xT,zT, focusparts[3]);
	gr3C->SetPoint(count2, xT,zT, focusparts[4]);
	gr3D->SetPoint(count2, xT,zT, focusparts[5]);
	count2++;
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

  TCanvas *c2=new TCanvas("c2","c2");
  c2->Divide(2,2);
  c2->cd(1);
  gr3A->Draw("cont4z");
  c2->cd(2);
  gr3B->Draw("cont4z");
  c2->cd(3);
  gr3C->Draw("cont4z");
  c2->cd(4);
  gr3D->Draw("cont4z");
  
}
