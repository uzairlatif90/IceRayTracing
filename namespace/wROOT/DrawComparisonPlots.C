#include "IceRayTracing.cc"

void DrawComparisonPlots(){
                                                                                             
  double zT[2]={-100,-200};                                                                            
  double ShowerHitDistance=100;
  double ShowerDepth=0;

  IceRayTracing::GridPositionXb.resize(2);
  IceRayTracing::GridPositionZb.resize(2);
  IceRayTracing::GridZValueb.resize(2);

  Int_t AntNum=0;
  IceRayTracing::MakeTable(ShowerHitDistance,ShowerDepth,zT[AntNum],AntNum);
  // AntNum=1;
  // IceRayTracing::MakeTable(ShowerHitDistance,ShowerDepth,zT[AntNum],AntNum);

  //AntNum=1;
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
  TH1D * h1error_dRR=new TH1D("","",200,0,100);

  TGraph2D* gr2b=new TGraph2D();
  TGraph2D* gr2c=new TGraph2D();
  TGraph2D* gr2corr=new TGraph2D();

  int count1=0;
  int count2=0;
  int count3=0;
  int count4=0;
  
  double minb=1e9;
  double minc=1e9;

  cout<<"start points are "<<IceRayTracing::GridStartX<<" "<<IceRayTracing::GridStartZ<<endl;
  
  for(int ix=0;ix<TotalStepsXb;ix++){
    for(int iz=0;iz<TotalStepsZb;iz++){
      
      double xR=IceRayTracing::GridStartX+GridStepSizeXb*ix;
      double zR=IceRayTracing::GridStartZ+GridStepSizeZb*iz;

      //cout<<xR<<" "<<zR<<endl;
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
      int rtParameter=1;
      ////For recording how much time the process took                                           
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
	gr2b->SetPoint(count1,xR,zR,rtresult);
	count1++;
      }
      
      if(rtresult<minb && rtresult!=-1000){
        minb=rtresult;
      }
      if(interresult<minc && interresult!=-1000){
        minc=interresult;
      }

      if(interresult!=-1000 && std::isnan(interresult)==false ){
	h2c->Fill(xR,zR,interresult);
	gr2c->SetPoint(count2,xR,zR,interresult);
	count2++;
      }
      if(rtresult!=-1000 && interresult!=-1000 && std::isnan(interresult)==false){
	h2corr->Fill(xR,zR,rtresult-interresult);
	h1error->Fill(rtresult-interresult);
	h1error_dRR->Fill((fabs(rtresult-interresult)/rtresult) *100);
	gr2corr->SetPoint(count3,xR,zR,((rtresult-interresult)/rtresult)*100);
	count3++;
      }
    }
  }

  h2b->GetZaxis()->SetRangeUser(minb,h2b->GetMaximum());
  h2c->GetZaxis()->SetRangeUser(minc,h2c->GetMaximum());

  gr2b->GetZaxis()->SetRangeUser(minb,gr2b->GetMaximum());
  gr2c->GetZaxis()->SetRangeUser(minc,gr2c->GetMaximum());

  gr2b->SetTitle("RayTrace results for OP_Ice; Horizontal Distance (m); Depth (m); OP_Ice (m)");
  gr2c->SetTitle("Interpolated results for OP_Ice; Horizontal Distance (m); Depth (m); OP_Ice (m)");
  gr2corr->SetTitle("Percentage Error for OP_Ice; Horizontal Distance (m); Depth (m); |#DeltaOP_Ice|/OP_Ice x 100");
  h1error->SetTitle("Absolute Error for OP_Ice; |#DeltaOP_Ice| (m) ; Cumulative Fraction of Tx Positions;");
  h1error_dRR->SetTitle("Percentage Error for OP_Ice; |#DeltaOP_Ice|/OP_Ice x 100 ; Cumulative Fraction of Tx Positions;");
  h1->SetTitle("Time taken to do interpolation; Duration (ns) ; Interpolation calls;");
  
   TCanvas *c2a=new TCanvas("c2a","c2a");
  c2a->cd();
  //c2a->SetLogz();
  c2a->SetGridx();
  c2a->SetGridy();
  gr2b->Draw("cont4z");
  c2a->SaveAs("OP_Ice_RT.png");
  
  TCanvas *c2b=new TCanvas("c2b","c2b");
  c2b->cd();
  //c2b->SetLogz();
  c2b->SetGridx();
  c2b->SetGridy();
  gr2c->Draw("cont4z");
  c2b->SaveAs("OP_Ice_Inter.png");

  TCanvas *c2c=new TCanvas("c2c","c2c");
  c2c->cd();
  ////c2c->SetLogz();
  c2c->SetGridx();
  c2c->SetGridy();
  gr2corr->Draw("cont4z");
  c2c->SaveAs("OP_Ice_RT_IT.png");

  TH1* h1error_cum = h1error->GetCumulative();
  TH1* h1error_dRR_cum = h1error_dRR->GetCumulative();

  h1error_cum->Scale(1./count3);
  h1error_dRR_cum->Scale(1./count3);

  TCanvas *c2d=new TCanvas("c2d","c2d");
  c2d->Divide(2,1);
  c2d->cd(1);
  ////c2d->cd(1)->SetLogy();
  c2d->cd(1)->SetGridx();
  c2d->cd(1)->SetGridy();
  h1error_cum->SetStats(kFALSE);
  h1error_cum->Draw("");

  c2d->cd(2);
  ////c2d->cd(2)->SetLogy();
  c2d->cd(2)->SetGridx();
  c2d->cd(2)->SetGridy();
  h1error_dRR_cum->SetStats(kFALSE);
  h1error_dRR_cum->Draw("");
  c2d->SaveAs("OP_Ice_RT_IT_hist.png");

  TCanvas *c2e=new TCanvas("c2e","c2e");
  c2e->cd();
  c2e->SetLogy();
  c2e->SetGridx();
  c2e->SetGridy();
  h1->Draw("");
  c2e->SaveAs("OP_Ice_Duration_hist.png");

  // TCanvas *c2f=new TCanvas("c2f","c2f");
  // c2f->cd();
  // //c2e->SetLogy();
  // c2f->SetGridx();
  // c2f->SetGridy();
  // //c2e->SaveAs("OP_Ice_Duration_hist.png");
  
  TCanvas *c1=new TCanvas("c1","c1");
  c1->Divide(2,2);
  c1->cd(1);
  //c1->cd(1)->SetLogz();
  c1->cd(1)->SetGridx();
  c1->cd(1)->SetGridy();
  gr2b->Draw("cont4z");
  c1->cd(2);
  //c1->cd(2)->SetLogz();
  c1->cd(2)->SetGridx();
  c1->cd(2)->SetGridy();
  gr2c->Draw("cont4z");
  c1->cd(3);
  ////c1->cd(3)->SetLogz();
  c1->cd(3)->SetGridx();
  c1->cd(3)->SetGridy();
  gr2corr->Draw("cont4z");
  c1->cd(4);
  ////c1->cd(4)->SetLogy();
  c1->cd(4)->SetGridx();
  c1->cd(4)->SetGridy();
  h1error_dRR_cum->Draw();
  c1->SaveAs("OP_Ice_All.png");

  // TCanvas *c3=new TCanvas("c3","c3");
  // c3->cd();
  // c3->SetLogy();
  // c3->SetGridx();
  // c3->SetGridy();
  // h1RTduration->Draw("");
                                                                                             
}
