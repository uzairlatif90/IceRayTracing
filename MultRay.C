#include "IceRayTracing_wROOTplot.C"

void MultRay(Double_t z0, Double_t LaunchInterval){

  ////Plot the ray solutions
  //Double_t z0=-180;//Tx z position
  Double_t LowerPlotLimit=z0;//Set the lower limit for the rays
  //Double_t LaunchInterval=0.25;//launch angle step size in degrees
  Int_t TotalRays=90.0/LaunchInterval;//total number of rays
  
  ////TGraphs to store ray paths
  TGraph *grR[TotalRays];
  TMultiGraph *mg=new TMultiGraph();

  Double_t DummyVariable=0;
  Double_t launchangle=0;//variable defined for the for-loop
  
  for(Int_t iang=0;iang<TotalRays;iang++){    
    launchangle=iang*LaunchInterval;
    Double_t B=IceRayTracing::GetB(z0);
    Double_t C=IceRayTracing::GetC(z0);
    
    Double_t lvalueR=(IceRayTracing::A_ice+B*exp(C*z0))*sin(launchangle*(IceRayTracing::pi/180));
    Double_t lvalueRa=(IceRayTracing::A_ice+B*exp(C*z0))*sin(launchangle*(IceRayTracing::pi/180));
    Double_t zn=LowerPlotLimit;
 
    zn=LowerPlotLimit;
    /* This function returns the x and z values for the full Reflected ray path in a TGraph and also prints out the ray path in a text file */
    grR[iang]=IceRayTracing::GetFullReflectedRayPath(z0,DummyVariable,LowerPlotLimit,lvalueR);

    /* Setup the function that will be used to calculate the angle of reception for all the rays */
    gsl_function F5;
    struct IceRayTracing::fDnfR_params paramsIAngB = {IceRayTracing::A_ice, IceRayTracing::GetB(0.0000001), IceRayTracing::GetC(0.0000001), lvalueR};
    double result, abserr;
    F5.function = &IceRayTracing::fDnfR;

    /* Calculate the angle of incidence of the reflected ray at the surface ice. This will be used to calculate the Fresnel Coefficients. The angle is calculated by calculating the derivative of the ray path fucnction at the surface*/
    F5.function = &IceRayTracing::fDnfR; 
    F5.params = &paramsIAngB;
    gsl_deriv_central (&F5, -0.0000001, 1e-8, &result, &abserr);
    double IncidenceAngleInIce=atan(result)*(180.0/IceRayTracing::pi);
    
    zn=LowerPlotLimit;
    /* This function returns the zmax for Reflected/Refracted ray path in a TGraph*/
    Double_t zmax=IceRayTracing::GetZmax(IceRayTracing::A_ice,lvalueRa)+0.0000001;
    if(zmax>1e-6){
       /* This function returns the x and z values for the full Refracted ray path in a TGraph and also prints out the ray path in a text file */
      grR[iang]=IceRayTracing::GetFullRefractedRayPath(z0,DummyVariable,LowerPlotLimit,zmax,lvalueRa,0);
    }
    mg->Add(grR[iang]);
  }//iang loop

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
