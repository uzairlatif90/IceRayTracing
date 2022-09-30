#include "IceRayTracing.cc"

void MakeRxTxPlots(){
  const char* filename;
  	
  Double_t dist=0;
  Double_t x0=0;
  Double_t z0=0;
  Double_t x1=0;
  Double_t z1=0;
  Double_t t0=0;
  Double_t t0error=0;

  Double_t ShowerLength=10;
  Int_t NTxTot=ShowerLength;
  Int_t NTxStrt=0;
  Int_t NTxStop=ShowerLength;

  // NTxStrt=9;
  // NTxStop=13;
  
  TGraph *gr[4][NTxTot];
  TGraph *grB[4][NTxTot];

  TGraph *grRx=new TGraph();
  grRx->SetMarkerStyle(21);
  grRx->SetMarkerColor(6);
  grRx->SetMarkerSize(2);

  TGraph *grTx=new TGraph();
  grTx->SetMarkerStyle(21);
  grTx->SetMarkerColor(7);
  grTx->SetMarkerSize(2);
  
  TGraph *grShwr=new TGraph();
  grShwr->SetMarkerStyle(20);
  grShwr->SetMarkerColor(kRed);
  grShwr->SetMarkerSize(2);
  grShwr->SetLineColor(12);
  grShwr->SetLineWidth(4);

  TGraph *grShwrTx=new TGraph();
  grShwrTx->SetMarkerStyle(20);
  grShwrTx->SetMarkerColor(kBlue);
  grShwrTx->SetLineColor(12);
  grShwrTx->SetLineWidth(4);
	
  Double_t TxX=-20;
  Double_t TxZ=-8;
  Double_t RxX=20;
  Double_t RxZ=-8;
  Double_t ShwrX=0;
  Double_t ShwrZ=-5;
  Double_t Tx2RxDist=fabs(TxX-RxX);
  Double_t ShwrDistRx=fabs(ShwrX-RxX);
  Double_t ShwrDistTx=fabs(ShwrX-TxX);

  cout<<Tx2RxDist<<" "<<TxX<<" "<<ShwrDistRx<<" "<<ShwrDistTx<<endl;
  
  grRx->SetPoint(0,RxX,RxZ);
  grTx->SetPoint(0,TxX,TxZ);
  
  Int_t nidpRx=0;
  Int_t nidpTx=0;
  Double_t RotationAngle=0;
	
  TMultiGraph *mg=new TMultiGraph();
  TCanvas *c1=new TCanvas();
  nidpRx=0;
  nidpTx=0;
  RotationAngle=90;
	
  cout<<"Shower at angle="<<RotationAngle<<endl;
  for(Int_t idp=NTxStrt;idp<NTxStop;idp++){
	  
    x0=RxX;
    z0=RxZ;
    if(ShwrX<RxX){
      x1=ShwrDistRx-(ShowerLength/2)+idp;
    }
    if(ShwrX>RxX){
      x1=ShwrDistRx+(ShowerLength/2)-idp;
    }
    z1=ShwrZ;
      
    TVector3 v(x1-ShwrDistRx,z1-ShwrZ,0);
    if(ShwrX<RxX){
      v.RotateZ(-RotationAngle*(TMath::Pi()/180.0));
    }
    if(ShwrX>RxX){
      v.RotateZ(RotationAngle*(TMath::Pi()/180.0));
    }
    x1=v.X()+ShwrDistRx;
    z1=v.Y()+ShwrZ;
	  
    Double_t CascadePropTime=((fabs(z1)/cos((90-RotationAngle)*(TMath::Pi()/180.0)))/IceRayTracing::c_light_ms)*pow(10,9);
	  
    if(z1==0){
      z1=-0.01;
    }
    if(z0==0){
      z0=-0.01;
    }
    cout<<"Rx to Shwr: "<<idp<<" "<<0<<" "<<z0<<" "<<x1<<" "<<z1<<endl;
    Double_t *getres1=IceRayTracing::IceRayTracing(0,z0,x1,z1);
    gr[0][idp]=IceRayTracing::GetFullDirectRayPath(z0,x1,z1,getres1[19]);
    gr[1][idp]=IceRayTracing::GetFullReflectedRayPath(z0,x1,z1,getres1[20]);
    gr[2][idp]=IceRayTracing::GetFullRefractedRayPath(z0,x1,z1,getres1[23],getres1[21],1);
    gr[3][idp]=IceRayTracing::GetFullRefractedRayPath(z0,x1,z1,getres1[24],getres1[22],2);

    gr[0][idp]->SetLineColor(46);
    gr[1][idp]->SetLineColor(46);
    gr[2][idp]->SetLineColor(46);
    gr[3][idp]->SetLineColor(46);
	  
    gr[0][idp]->SetLineWidth(2);
    gr[1][idp]->SetLineWidth(2);
    gr[2][idp]->SetLineWidth(2);
    gr[3][idp]->SetLineWidth(2);
      
    for(Int_t ic=0;ic<4;ic++){
      for(Int_t inum=0;inum<gr[ic][idp]->GetN();inum++){
	Double_t x,y;
	gr[ic][idp]->GetPoint(inum,x,y);
	if(ShwrX<RxX){
	  gr[ic][idp]->SetPoint(inum,-x+RxX,y);
	}
	if(ShwrX>RxX){
	  gr[ic][idp]->SetPoint(inum,x+RxX,y);
	}	  
      }
    }
      
    if(ShwrX<RxX){
      x1=-x1+RxX;
    }
    if(ShwrX>RxX){
      x1=x1+RxX;
    }
      
      
    if(getres1[8]!=-1000){
      //cout<<"Direct to Rx: "<<getres1[3]*pow(10,9)<<" "<<CascadePropTime<<endl;
      mg->Add(gr[0][idp]);;
      grShwr->SetPoint(nidpRx,x1,z1);
      nidpRx++;
    }
    if(getres1[9]!=-1000){
      //cout<<"Reflected to Rx: "<<getres1[4]*pow(10,9)<<" "<<CascadePropTime<<endl;
      mg->Add(gr[1][idp]);;
      grShwr->SetPoint(nidpRx,x1,z1);
      nidpRx++;
    }
    if(getres1[10]!=-1000){
      //cout<<"Refracted to Rx: "<<getres1[5]*pow(10,9)<<" "<<CascadePropTime<<endl;
      mg->Add(gr[2][idp]);;
      grShwr->SetPoint(nidpRx,x1,z1);
      nidpRx++;
    }
    if(getres1[11]!=-1000){
      //cout<<"Refracted to Rx: "<<getres1[5]*pow(10,9)<<" "<<CascadePropTime<<endl;
      mg->Add(gr[3][idp]);;
      grShwr->SetPoint(nidpRx,x1,z1);
      nidpRx++;
    }
      
    x0=TxX;
    z0=TxZ;
    if(ShwrX<TxX){
      x1=ShwrDistTx-(ShowerLength/2)+idp;
    }
    if(ShwrX>TxX){
      x1=ShwrDistTx+(ShowerLength/2)-idp;
    }
    z1=ShwrZ;
    //cout<<"Before Tx to Shwr: "<<idp<<" "<<0<<" "<<z0<<" "<<x1<<" "<<z1<<endl;
      
    TVector3 vB(x1-ShwrDistTx,z1-ShwrZ,0);
    if(ShwrX<TxX){
      vB.RotateZ(-RotationAngle*(TMath::Pi()/180.0));
    }
    if(ShwrX>TxX){
      vB.RotateZ(RotationAngle*(TMath::Pi()/180.0));
    }
      
    x1=vB.X()+ShwrDistTx;
    z1=vB.Y()+ShwrZ;
    if(z1==0){
      z1=-0.01;
    }
    if(z0==0){
      z0=-0.01;
    }
      
    cout<<"Tx to Shwr: "<<idp<<" "<<0<<" "<<z0<<" "<<x1<<" "<<z1<<endl;
    Double_t *getres2=IceRayTracing::IceRayTracing(0,z0,x1,z1);
    grB[0][idp]=IceRayTracing::GetFullDirectRayPath(z0,x1,z1,getres2[19]);
    grB[1][idp]=IceRayTracing::GetFullReflectedRayPath(z0,x1,z1,getres2[20]);
    grB[2][idp]=IceRayTracing::GetFullRefractedRayPath(z0,x1,z1,getres2[23],getres2[21],1);
    grB[3][idp]=IceRayTracing::GetFullRefractedRayPath(z0,x1,z1,getres2[24],getres2[22],2);
    
    grB[0][idp]->SetLineColor(9);
    grB[1][idp]->SetLineColor(9);
    grB[2][idp]->SetLineColor(9);
    grB[3][idp]->SetLineColor(9);
      
    grB[0][idp]->SetLineWidth(2);
    grB[1][idp]->SetLineWidth(2);
    grB[2][idp]->SetLineWidth(2);
    grB[3][idp]->SetLineWidth(2);
      
    for(Int_t ic=0;ic<4;ic++){
      for(Int_t inum=0;inum<grB[ic][idp]->GetN();inum++){
	Double_t x,y;
	grB[ic][idp]->GetPoint(inum,x,y);
	if(ShwrX<TxX){
	  grB[ic][idp]->SetPoint(inum,-x+TxX,y);
	}
	if(ShwrX>TxX){
	  grB[ic][idp]->SetPoint(inum,x+TxX,y);
	}	  
      }
    }
      
    if(ShwrX<TxX){
      x1=-x1+TxX;
    }
    if(ShwrX>TxX){
      x1=x1+TxX;
    }
     
    if(getres2[8]!=-1000){
      //cout<<"Direct to Tx: "<<getres2[3]*pow(10,9)+CascadePropTime<<endl;
      mg->Add(grB[0][idp]);;
      grShwrTx->SetPoint(nidpTx,x1,z1);
      nidpTx++;
    }
    if(getres2[9]!=-1000){
      //cout<<"Reflected to Tx: "<<getres2[4]*pow(10,9)+CascadePropTime<<endl;
      mg->Add(grB[1][idp]);;
      grShwrTx->SetPoint(nidpTx,x1,z1);
      nidpTx++;
    }
    if(getres2[10]!=-1000){
      //cout<<"Refracted to Tx: "<<getres2[5]*pow(10,9)+CascadePropTime<<endl;
      mg->Add(grB[2][idp]);;
      grShwrTx->SetPoint(nidpTx,x1,z1);
      nidpTx++;
    }
    if(getres2[11]!=-1000){
      //cout<<"Refracted to Tx: "<<getres2[5]*pow(10,9)+CascadePropTime<<endl;
      mg->Add(grB[3][idp]);;
      grShwrTx->SetPoint(nidpTx,x1,z1);
      nidpTx++;
    }
      
    delete []getres2;
    delete []getres1;
      
  }////idp loop
  mg->Add(grShwr);
  mg->Add(grShwrTx);
  mg->Add(grRx);
  mg->Add(grTx);
	
  c1->cd();
  c1->SetGridx();
  c1->SetGridy();
  TString mgtitle=";Distance (m); Depth (m)";
  mg->SetTitle(mgtitle);
  mg->Draw("APL");
  TString file="./plots/RxX=";
  file+=0;
  file+="_RxDepth=";
  file+=(int)fabs(RxZ);
  file+="TxX=";
  file+=(int)TxX;
  file+="_TxDepth=";
  file+=(int)fabs(TxZ);
  file+=".png";
  filename=file.Data();
  c1->SaveAs(filename);	

}
