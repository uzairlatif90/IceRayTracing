#include "IceRayTracing.cc"

void MultRay(){
  
  /*These are all dummy variables you do not need to adjust them*/
  Double_t x0=0;//Tx x positions 
  Double_t x1=1000;//Rx dist position. This is just a dummy variable for now
  Double_t launchangle=0;//variable defined for the for loop
  
  Int_t totray=2000;//90*2*2;//total number of rays
  //Double_t launchinterval=0.05;//launch angle step size in degrees
  Double_t launchinterval=0.25;//launch angle step size in degrees
  
  ////store ray paths
  //TGraph *grR[totray];
  TMultiGraph *mg=new TMultiGraph();

  TGraph *gr1=new TGraph();
  TGraph *gr2=new TGraph();
  
  for(Int_t iTx=1;iTx<101;iTx++){

  //Plot the ray solutions
  Double_t z0=-iTx;//Tx z position. The rays will start at this depth
  //Double_t zlim=-(iTx+10);//How deep do you want the rays to go into the ice
  Double_t z1=z0;//Rx z position.  This is just a dummy variable for now
  
  Int_t iang=0;
  Bool_t RefrReached=false;
 
  while(RefrReached==false){
    launchangle=iang*launchinterval;
    
    //Double_t lvalueR=(IceRayTracing::A_ice+IceRayTracing::GetB(z0)*exp(IceRayTracing::GetC(z0)*z0))*sin(launchangle*(IceRayTracing::pi/180));
    Double_t lvalueRa=(IceRayTracing::A_ice+IceRayTracing::GetB(z0)*exp(IceRayTracing::GetC(z0)*z0))*sin(launchangle*(IceRayTracing::pi/180));
    Double_t zn=z1;

    /* This function returns the zmax for Reflected/Refracted ray path in a TGraph*/
    Double_t zmax=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa)+0.0000001;
    
    if(zmax>1e-5){
      struct IceRayTracing::fDnfR_params params6a;
      struct IceRayTracing::fDnfR_params params6b;
      struct IceRayTracing::fDnfR_params params6c;
      struct IceRayTracing::fDnfR_params params6d;
      
      params6a = {IceRayTracing::A_ice,IceRayTracing::GetB(z0), -IceRayTracing::GetC(z0), lvalueRa};
      params6b = {IceRayTracing::A_ice,IceRayTracing::GetB(zmax), -IceRayTracing::GetC(zmax), lvalueRa};
      Double_t xn1=+IceRayTracing::fDnfR(-z0,&params6a)-IceRayTracing::fDnfR(zmax,&params6b);
      
      params6c = {IceRayTracing::A_ice,IceRayTracing::GetB(z0), -IceRayTracing::GetC(z0), lvalueRa};
      params6d = {IceRayTracing::A_ice,IceRayTracing::GetB(5), -IceRayTracing::GetC(5), lvalueRa};
      Double_t xn2=-IceRayTracing::fDnfR(zmax,&params6c)+IceRayTracing::fDnfR(5,&params6d)+xn1;
      
      gr1->SetPoint(iTx-1,xn1,iTx);
      gr2->SetPoint(iTx-1,xn2,iTx);
      
      /* Thisfunction returns the x and z values for the full Refracted ray path in a TGraph */
      //grR[iang]=IceRayTracing::GetFullRefractedRayPath(z0,x1,z1,zmax,lvalueRa,zlim);
      RefrReached=true;
      //mg->Add(grR[iang]);
      cout<<"zmax is "<<zmax<<" "<<xn1<<" "<<xn2<<" "<<iang<<endl;
    }// else{
    //   zn=z1;
    //   /* This function returns the x and z values for the full Reflected ray path in a TGraph */
    //   grR[iang]=IceRayTracing::GetFullReflectedRayPath(z0,x1,z1,lvalueR,zlim);
    // }
    //grR[iang]->SetLineColor(kBlue);
    
    
    iang=iang+1;
  }//while loop
  }//iTx loop

  // TString title="Depth vs Distance, Tx Depth=";
  // title+=z0;
  // title+=" m; Distance (m);Depth (m)";
  // mg->SetTitle(title); 

  TString title="";
  title+="; Horizon distance (m); Tx Depth (m)";
  mg->SetTitle(title);
  mg->Add(gr1);
  mg->Add(gr2);
  
  TCanvas *c1=new TCanvas("c1","c1");
  c1->cd();
  gr1->SetMarkerStyle(20);
  gr2->SetMarkerStyle(20);

  gr1->SetTitle("Peak Point of Ray");
  gr2->SetTitle("5 m Below Peak Point");
  
  //gr->SetMarkerColor(kRed);
  gr2->SetMarkerColor(kRed);

  mg->Draw("ALP");
  mg->GetXaxis()->SetNdivisions(20);
  mg->GetYaxis()->SetNdivisions(20);
  c1->SetGridx();
  c1->SetGridy();
  c1->BuildLegend();
  
}
