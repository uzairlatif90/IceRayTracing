#include "IceRayTracing.cc"

void SelectPoints(){
  
  //Plot the ray solutions
  Double_t z0=-10;//Tx z position. The rays will start at this depth
  Double_t zlim=-20;//How deep do you want the rays to go into the ice
  
  /*These are all dummy variables you do not need to adjust them*/
  Double_t x0=0;//Tx x positions 
  Double_t x1=1000;//Rx dist position. This is just a dummy variable for now
  Double_t z1=z0;//Rx z position.  This is just a dummy variable for now
  Double_t launchangle=0;//variable defined for the for loop
  
  Int_t totray=90*2*2;//total number of rays
  Double_t launchinterval=0.25;//launch angle step size in degrees

  ////store ray paths
  TGraph *grR[totray];
  TMultiGraph *mg=new TMultiGraph();
  
  for(Int_t iang=297;iang<298;iang++){
  //for(Int_t iang=295;iang<305;iang++){
    launchangle=iang*launchinterval;
    
    Double_t lvalueR=(IceRayTracing::A_ice+IceRayTracing::GetB(z0)*exp(IceRayTracing::GetC(z0)*z0))*sin(launchangle*(IceRayTracing::pi/180));
    Double_t lvalueRa=(IceRayTracing::A_ice+IceRayTracing::GetB(z0)*exp(IceRayTracing::GetC(z0)*z0))*sin(launchangle*(IceRayTracing::pi/180));
    Double_t zn=z1;

    /* This function returns the zmax for Reflected/Refracted ray path in a TGraph*/
    Double_t zmax=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa)+0.0000001;
    
    if(zmax>1e-5){
      cout<<iang<<endl;
      zn=z1;
      /* This function returns the x and z values for the full Refracted ray path in a TGraph */
      grR[iang]=IceRayTracing::GetFullRefractedRayPath(z0,x1,z1,zmax,lvalueRa,0);//,zlim);
    }else{
      zn=z1;
      /* This function returns the x and z values for the full Reflected ray path in a TGraph */
      grR[iang]=IceRayTracing::GetFullReflectedRayPath(z0,x1,z1,lvalueR);//,zlim);
    }
    grR[iang]->SetLineColor(kBlue);
    
    mg->Add(grR[iang]);
  }//iang loop

  TGraph *grnew1=new TGraph();
  TGraph *grnew2=new TGraph();

  grnew1->SetPoint(0,69.2356,-0.08665);
  grnew2->SetPoint(0,117.79,-5);

  grnew1->SetMarkerStyle(20);
  grnew2->SetMarkerStyle(20);

  //grnew1->SetMarkerStyle(20);
  grnew2->SetMarkerColor(kRed);

  
  mg->Add(grnew1);
  mg->Add(grnew2);
  
  TString title="Depth vs Distance, Tx Depth=";
  title+=z0;
  title+=" m; Distance (m);Depth (m)";
  mg->SetTitle(title); 
  
  TCanvas *c1=new TCanvas("c1","c1");
  c1->cd();
  mg->Draw("ALP");
  mg->GetXaxis()->SetNdivisions(20);
  c1->SetGridx();
  
}
