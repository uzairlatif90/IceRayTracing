#include "IceRayTracing.cc"

void DrawComparisonPlots(){
                                                                                             
  double zT[2]={-100,-200};                                                                            
  double ShowerHitDistance=50;

  IceRayTracing::GridPositionXb.resize(2);
  IceRayTracing::GridPositionZb.resize(2);
  IceRayTracing::GridZValueb.resize(2);

  Int_t AntNum=0;
  IceRayTracing::MakeTable(ShowerHitDistance,zT[AntNum],AntNum);
  AntNum=1;
  IceRayTracing::MakeTable(ShowerHitDistance,zT[AntNum],AntNum);

  AntNum=1;
  IceRayTracing::GridStartX=IceRayTracing::GridPositionXb[AntNum][0];
  IceRayTracing::GridStartZ=IceRayTracing::GridPositionZb[AntNum][0];
 
  IceRayTracing::GridStopX=IceRayTracing::GridPositionXb[AntNum][IceRayTracing::GridPositionXb[AntNum].size()-1];
  IceRayTracing::GridStopZ=IceRayTracing::GridPositionZb[AntNum][IceRayTracing::GridPositionZb[AntNum].size()-1];
 
  double GridStepSizeXb=0.13;
  double GridStepSizeZb=0.13;

  double TotalStepsXb=IceRayTracing::GridWidthX/GridStepSizeXb+1;
  double TotalStepsZb=IceRayTracing::GridWidthZ/GridStepSizeZb+1;  
  
  TH2D *h2b=new TH2D("","",TotalStepsXb+200,IceRayTracing::GridStartX-1,IceRayTracing::GridStopX+1,TotalStepsZb+200,IceRayTracing::GridStartZ-1,IceRayTracing::GridStopZ+1);
  TH2D *h2c=new TH2D("","",TotalStepsXb+200,IceRayTracing::GridStartX-1,IceRayTracing::GridStopX+1,TotalStepsZb+200,IceRayTracing::GridStartZ-1,IceRayTracing::GridStopZ+1);
  TH2D *h2corr=new TH2D("","",TotalStepsXb+200,IceRayTracing::GridStartX-1,IceRayTracing::GridStopX+1,TotalStepsZb+200,IceRayTracing::GridStartZ-1,IceRayTracing::GridStopZ+1);
  TH1D * h1=new TH1D("","",200,0,500);
  TH1D * h1error=new TH1D("","",100,-2,2);

  double minb=1e9;
  double minc=1e9;
  
  for(int ix=0;ix<TotalStepsXb;ix++){
    for(int iz=0;iz<TotalStepsZb;iz++){
      
      double xR=IceRayTracing::GridStartX+GridStepSizeXb*ix;
      double zR=IceRayTracing::GridStartZ+GridStepSizeZb*iz;

      double TimeRay[2]={0,0};
      double PathRay[2]={0,0};
      double LaunchAngle[2]={0,0};
      double RecieveAngle[2]={0,0};
      int IgnoreCh[2]={0,0};
      double IncidenceAngleInIce[2]={0,0};
      double A0=1;
      double frequency=0.1;
      double AttRay[2]={0,0};
  
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
      int rtParameter=0; ///0 is for launch angle, 1 is for recieve angle,  2 is for distance, 3 \
  // // 			 ////For recording how much time the process took                                           
      auto t1b = std::chrono::high_resolution_clock::now();                                      
      double NewZValue=IceRayTracing::GetInterpolatedValue(xR, zR, rtParameter,AntNum);
      auto t2b = std::chrono::high_resolution_clock::now();                                      
      double Duration = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b - t1b ).count();
      //cout<<"total time taken to interpolate: "<<Duration<<" nano s "<<endl;                     
      
      rtresult=output[rtParameter];
      interresult=NewZValue;

      output.clear();
      
      // if(RTresults[0]>1){
      // 	rtresult=RTresults[0];

      if(rtresult!=-1000){
	h2b->Fill(xR,zR,rtresult);
      }
      
      if(rtresult<minb && rtresult!=-1000){
        minb=rtresult;
      }
      if(interresult<minc && interresult!=-1000){
        minc=interresult;
      }

      if(interresult!=-1000){
	h2c->Fill(xR,zR,interresult);
      }
      if(rtresult!=-1000 && interresult!=-1000){
	h2corr->Fill(xR,zR,rtresult-interresult);
	h1error->Fill(rtresult-interresult);
      }
    }
  }
	//}
      
  // double xR=50;
  // double zR=-15.6;
  
  h2b->GetZaxis()->SetRangeUser(minb,h2b->GetMaximum());
  h2c->GetZaxis()->SetRangeUser(minc,h2c->GetMaximum());
  
  TCanvas *c1=new TCanvas("c1","c1");
  c1->Divide(2,2);
  c1->cd(1);
  h2b->Draw("colz");
  c1->cd(2);
  h2c->Draw("colz");
  c1->cd(3);
  h2corr->Draw("colz");
  c1->cd(4);
  h1->Draw();
  
  TCanvas *c2=new TCanvas("c2","c2");
  c2->cd();
  h1error->Draw();
  
  
  // // cout<<"True Value ="<<rtresult<<" ,Interpolated Value ="<<interresult<<" ,True - Interpolated value "<<rtresult-interresult<<endl;                                                     
  
                                                                                             
}
