#include "IceRayTracing_wROOTplot.C"

void MakeTxShwrRxPlot(){
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
	
  Double_t TxX=20;
  Double_t TxZ=-5;
  Double_t RxX=0;
  Double_t RxZ=-5;
  Double_t ShwrX=10;
  Double_t ShwrZ=-3;
  Double_t Tx2RxDist=fabs(TxX-RxX);
  Double_t ShwrDistRx=fabs(ShwrX-RxX);
  Double_t ShwrDistTx=fabs(ShwrX-TxX);
  
  grRx->SetPoint(0,RxX,RxZ);
  grTx->SetPoint(0,TxX,TxZ);
  
  Int_t nidpRx=0;
  Int_t nidpTx=0;
  Double_t RotationAngle=40;
	
  TMultiGraph *mg=new TMultiGraph();

  vector <double> ShwrTimeTx[4];
  vector <double> ShwrTimeRx[4];
  vector <double> ShwrTime[4];
  vector <double> ShwrPntIndex;
  
  cout<<"Shower at angle="<<RotationAngle<<endl;
  for(Int_t idp=NTxStrt;idp<NTxStop;idp++){

    Bool_t checkfirsthit=false;
    Double_t *getres1;
    Double_t *getres2;
    Double_t Dtime[2];
    Double_t Rtime[2];
    
    x0=RxX;
    z0=RxZ;   

    if(ShwrX<RxX){
      x1=ShwrDistRx-(ShowerLength/2)+idp;
    }
    if(ShwrX>RxX){
      x1=ShwrDistRx+(ShowerLength/2)-idp;
    }
    if(ShwrX==RxX){
      x1=(ShowerLength/2)-idp;
    }
    z1=ShwrZ;

    TVector3 v(x1-ShwrDistRx,z1-ShwrZ,0);
    if(ShwrX<RxX){
      v.RotateZ(-RotationAngle*(TMath::Pi()/180.0));
    }
    if(ShwrX>RxX){
      v.RotateZ(RotationAngle*(TMath::Pi()/180.0));
    }
    if(ShwrX==RxX){
      v.RotateZ(RotationAngle*(TMath::Pi()/180.0));
    }
    
    x1=v.X()+ShwrDistRx;
    z1=v.Y()+ShwrZ;

    if(z1<0){
    
      Double_t CascadePropTime=((fabs(z1)/cos((90-RotationAngle)*(TMath::Pi()/180.0)))/IceRayTracing::c_light_ms)*pow(10,9);
	  
      if(z1==0){
	z1=-0.01;
      }
      if(z0==0){
	z0=-0.01;
      }

      dist=fabs(x1);
      getres1=IceRayTracing::IceRayTracing(0,z0,dist,z1,false);
      
      gr[0][idp]=IceRayTracing::GetFullDirectRayPath(z0,dist,z1,getres1[19]);
      gr[1][idp]=IceRayTracing::GetFullReflectedRayPath(z0,dist,z1,getres1[20]);

      if((getres1[9]==-1000 || getres1[8]==-1000) && getres1[10]!=-1000){  
	gr[2][idp]=IceRayTracing::GetFullRefractedRayPath(z0,x1,z1,getres1[23],getres1[21],1);
	if(getres1[11]!=-1000){ 
	  gr[3][idp]=IceRayTracing::GetFullRefractedRayPath(z0,x1,z1,getres1[24],getres1[22],2);
	}
      }

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
	    if(x1>0){
	      gr[ic][idp]->SetPoint(inum,-x+RxX,y);
	    }else{
	      gr[ic][idp]->SetPoint(inum,x+RxX,y);
	    }
	  }
	  if(ShwrX>RxX){
	    if(x1>0){
	      gr[ic][idp]->SetPoint(inum,x+RxX,y);
	    }else{
	      gr[ic][idp]->SetPoint(inum,-x+RxX,y);
	    }
	  }	  
	  if(ShwrX==RxX){
	    if(x1>0){
	      gr[ic][idp]->SetPoint(inum,x+RxX,y);
	    }else{
	      gr[ic][idp]->SetPoint(inum,-x+RxX,y);
	    }	    
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
	Dtime[0]=getres1[4]*pow(10,9);
	mg->Add(gr[0][idp]);;
	grShwr->SetPoint(nidpRx,x1,z1);
	nidpRx++;
      }
      if(getres1[9]!=-1000){
	Rtime[0]=getres1[5]*pow(10,9);
	mg->Add(gr[1][idp]);;
	grShwr->SetPoint(nidpRx,x1,z1);
	nidpRx++;
      }
      if(getres1[10]!=-1000){
	Dtime[0]=getres1[6]*pow(10,9);
	mg->Add(gr[2][idp]);;
	grShwr->SetPoint(nidpRx,x1,z1);
	nidpRx++;
      }
      if(getres1[11]!=-1000){
	Rtime[0]=getres1[7]*pow(10,9);
	mg->Add(gr[3][idp]);;
	grShwr->SetPoint(nidpRx,x1,z1);
	nidpRx++;
      }
      
      if(getres1[8]!=-1000 || getres1[9]!=-1000 || getres1[10]!=-1000 || getres1[11]!=-1000){
	checkfirsthit=true;
      }
    
    }

    x0=TxX;
    z0=TxZ;
    if(ShwrX<TxX){
      x1=ShwrDistTx-(ShowerLength/2)+idp;
    }
    if(ShwrX>TxX){
      x1=ShwrDistTx+(ShowerLength/2)-idp;
    }
    if(ShwrX==TxX){
      x1=(ShowerLength/2)-idp;
    }
    z1=ShwrZ;
      
    TVector3 vB(x1-ShwrDistTx,z1-ShwrZ,0);
    if(ShwrX<TxX){
      vB.RotateZ(-RotationAngle*(TMath::Pi()/180.0));
    }
    if(ShwrX>TxX){
      vB.RotateZ(RotationAngle*(TMath::Pi()/180.0));
    }
    if(ShwrX==TxX){
      vB.RotateZ(RotationAngle*(TMath::Pi()/180.0));
    }
    x1=(vB.X()+ShwrDistTx);
    z1=vB.Y()+ShwrZ;
    
    if(z1<0){
      if(z1==0){
	z1=-0.01;
      }
      if(z0==0){
	z0=-0.01;
      }

      dist=fabs(x1);      
      getres2=IceRayTracing::IceRayTracing(0,z0,dist,z1,false);
      
      grB[0][idp]=IceRayTracing::GetFullDirectRayPath(z0,dist,z1,getres2[19]);
      grB[1][idp]=IceRayTracing::GetFullReflectedRayPath(z0,dist,z1,getres2[20]);

      if((getres2[9]==-1000 || getres2[8]==-1000) && getres2[10]!=-1000){  
	grB[2][idp]=IceRayTracing::GetFullRefractedRayPath(z0,x1,z1,getres2[23],getres2[21],1);
	if(getres2[11]!=-1000){ 
	  grB[3][idp]=IceRayTracing::GetFullRefractedRayPath(z0,x1,z1,getres2[24],getres2[22],2);
	}
      }
      
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
	    if(x1>0){
	      grB[ic][idp]->SetPoint(inum,-x+TxX,y);
	    }else{
	      grB[ic][idp]->SetPoint(inum,x+TxX,y);
	    }
	  }
	  if(ShwrX>TxX){
	    if(x1>0){
	      grB[ic][idp]->SetPoint(inum,x+TxX,y);
	    }else{
	      grB[ic][idp]->SetPoint(inum,-x+TxX,y);
	    }
	  }
	  if(ShwrX==TxX){
	    if(x1>0){
	      grB[ic][idp]->SetPoint(inum,x+TxX,y);
	    }else{
	      grB[ic][idp]->SetPoint(inum,-x+TxX,y);
	    }	    
	  }
	}
      }
      
      if(ShwrX<TxX){
	x1=-x1+TxX;
      }
      if(ShwrX>TxX){
	x1=x1+TxX;
      }
      if(ShwrX==TxX){
	x1=x1+TxX;
      }

      if(getres2[8]!=-1000){
	Dtime[1]=getres2[4]*pow(10,9);
	mg->Add(grB[0][idp]);;
	grShwrTx->SetPoint(nidpTx,x1,z1);
	nidpTx++;
      }
      if(getres2[9]!=-1000){
	Rtime[1]=getres2[5]*pow(10,9);
	mg->Add(grB[1][idp]);;
	grShwrTx->SetPoint(nidpTx,x1,z1);
	nidpTx++;
      }
      if(getres2[10]!=-1000){
	Dtime[1]=getres2[6]*pow(10,9);
	mg->Add(grB[2][idp]);;
	grShwrTx->SetPoint(nidpTx,x1,z1);
	nidpTx++;
      }
      if(getres2[11]!=-1000){
	Rtime[1]=getres2[7]*pow(10,9);
	mg->Add(grB[3][idp]);;
	grShwrTx->SetPoint(nidpTx,x1,z1);
	nidpTx++;
      }
      
      if(getres2[8]!=-1000 || getres2[9]!=-1000 || getres2[10]!=-1000 || getres2[11]!=-1000){
	if(checkfirsthit==true){
	  ShwrPntIndex.push_back(idp);
	  ShwrTime[0].push_back(Dtime[0]+Dtime[1]);
	  ShwrTime[1].push_back(Rtime[0]+Rtime[1]);
	  ShwrTime[2].push_back(Dtime[0]+Rtime[1]);
	  ShwrTime[3].push_back(Rtime[0]+Dtime[1]);

	  ShwrTimeRx[0].push_back(Dtime[0]);
	  ShwrTimeRx[1].push_back(Rtime[0]);
	  ShwrTimeRx[2].push_back(Dtime[0]);
	  ShwrTimeRx[3].push_back(Rtime[0]);

	  ShwrTimeTx[0].push_back(Dtime[1]);
	  ShwrTimeTx[1].push_back(Rtime[1]);
	  ShwrTimeTx[2].push_back(Rtime[1]);
	  ShwrTimeTx[3].push_back(Dtime[1]);
	}
      }
     
    }

    // delete []getres1;
    // delete []getres2;
    
  }////idp loop
  mg->Add(grShwr);
  mg->Add(grShwrTx);
  mg->Add(grRx);
  mg->Add(grTx);

  TCanvas *c1=new TCanvas();
  c1->cd();
  c1->SetGridx();
  c1->SetGridy();
  TString mgtitle=";Distance (m); Depth (m)";
  mg->SetTitle(mgtitle);
  mg->Draw("APL");
  TString file="RxX=";
  file+=0;
  file+="_RxDepth=";
  file+=(int)fabs(RxZ);
  file+="TxX=";
  file+=(int)TxX;
  file+="_TxDepth=";
  file+=(int)fabs(TxZ);
  file+=".png";
  filename=file.Data();
  //c1->SaveAs(filename);  

  // TGraph * grRxTxTimeReal[4];
  // TMultiGraph *mgtest[4];  

  // TGraph * grTxTimeReal[4];
  // TGraph * grRxTimeReal[4];

  // int totalpoints=ShwrPntIndex.size();
  
  // for (int irxtx = 0; irxtx < 4; irxtx++){
  //   grRxTxTimeReal[irxtx]=new TGraph();
  //   mgtest[irxtx]=new TMultiGraph();
  //   grTxTimeReal[irxtx]=new TGraph();
  //   grRxTimeReal[irxtx]=new TGraph();
   
  //   for (int ixy = 0; ixy < totalpoints; ixy++){
  //     grRxTxTimeReal[irxtx]->SetPoint(ixy,ShwrPntIndex[ixy],(ShwrTime[irxtx][ixy]));
  //     grRxTimeReal[irxtx]->SetPoint(ixy,ShwrPntIndex[ixy],(ShwrTimeRx[irxtx][ixy]));
  //     grTxTimeReal[irxtx]->SetPoint(ixy,ShwrPntIndex[ixy],(ShwrTimeTx[irxtx][ixy]));
  //   }
    	
  //   grRxTxTimeReal[irxtx]->SetMarkerStyle(20);
  //   grRxTxTimeReal[irxtx]->SetMarkerSize(4);

  //   grRxTimeReal[irxtx]->SetMarkerStyle(20);
  //   grRxTimeReal[irxtx]->SetMarkerSize(4);

  //   grTxTimeReal[irxtx]->SetMarkerStyle(20);
  //   grTxTimeReal[irxtx]->SetMarkerSize(4);

  //   mgtest[irxtx]->Add(grRxTxTimeReal[irxtx]); 
  // }

  // TCanvas *cRxInter=new TCanvas("cRxInter","cRxInter");
  // cRxInter->Divide(2,2);
  // cRxInter->cd(1);
  // grRxTimeReal[0]->SetTitle("D(Rx);Index of the point on the shower;D(Rx) Ray propagation Time (ns);");
  // grRxTimeReal[0]->Draw("AP");

  // cRxInter->cd(2);
  // grRxTimeReal[1]->SetTitle("R(Rx);Index of the point on the shower;R(Rx) Ray propagation Time (ns);");
  // grRxTimeReal[1]->Draw("AP");
  
  // cRxInter->cd(3);
  // grTxTimeReal[0]->SetTitle("D(Tx);Index of the point on the shower;D(Tx) Ray propagation Time (ns);");
  // grTxTimeReal[0]->Draw("AP");

  // cRxInter->cd(4);
  // grTxTimeReal[1]->SetTitle("R(Tx);Index of the point on the shower;R(Tx) Ray propagation Time (ns);");
  // grTxTimeReal[1]->Draw("AP");

  
  // TCanvas *cRxTxInter=new TCanvas("cRxTxInter","cRxTxInter");
  // cRxTxInter->Divide(2,2);
  // cRxTxInter->cd(1);
  // mgtest[0]->SetTitle("D(Rx)+D(Tx);Index of the point on the shower;D(Rx)+D(Tx) Ray propagation Time (ns);");
  // mgtest[0]->Draw("AP");

  // cRxTxInter->cd(2);
  // mgtest[1]->SetTitle("R(Rx)+R(Tx);Index of the point on the shower;R(Rx)+R(Tx) Ray propagation Time (ns);");
  // mgtest[1]->Draw("AP");
  
  // cRxTxInter->cd(3);
  // mgtest[2]->SetTitle("D(Rx)+R(Tx);Index of the point on the shower;D(Rx)+R(Tx) Ray propagation Time (ns);");
  // mgtest[2]->Draw("AP");

  // cRxTxInter->cd(4);
  // mgtest[3]->SetTitle("R(Rx)+D(Tx);Index of the point on the shower;R(Rx)+D(Tx) Ray propagation Time (ns);");
  // mgtest[3]->Draw("AP");
  
}
