// Goal: produce ratio plot of Data/MC for dsigma/dpt vs pt. Eta binning integrated from -2.1 < eta < 2.1. Works for pp or PbPb by flipping kSample switch
// KEEP THESE SCRIPTS SMALL AND SINGLE PURPOSED!

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
#include "TH1D.h"

#include "../commonUtility.h"
#include "../jzWeight.h"
#include "../dataWeight.h"
// this makes a class out of ntuple so all struct values can be accessed -- there were problems accessing the final memebers of the struct directly
#include "../ntuples/tashka.C"
#endif
#include <TMath.h>
#include <TLatex.h>
#include <TPaletteAxis.h>
#include "unfoldingUtil.h"
#include "systematicsTool.h"

// fraction of file you want to process ie) 0.5 would process half the file
double statFrac = 1.;
// set this to either kPP or kPbPb before running
int kSample = kPP;

int main(int argc, char *argv[]) {
  
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  gStyle->SetOptTitle(0);
  int nCentBins;
  if (kSample == kPP) nCentBins = 1;
  else nCentBins = 7; 
  int nEtaBins = 6;  

  TString jz2;
  TString jz3;
  TString jz4;
  if (kSample == kPP) {
    jz2 = jz2PPString;
    jz3 = jz3PPString;
    jz4 = jz4PPString;
  }
  else {
    jz2 = jz2PbPbString;
    jz3 = jz3PbPbString;
    jz4 = jz4PbPbString;
  }
 
  cout << " Setting Tashkas for mc ..." << endl;
  TFile* fjz2 = new TFile(Form("../ntuples/%s",jz2.Data()));
  TTree* tr2 = (TTree*)fjz2->Get("tr");
  tashka t_jz2;

  TFile* fjz3 = new TFile(Form("../ntuples/%s",jz3.Data()));
  TTree* tr3 = (TTree*)fjz3->Get("tr");
  tashka t_jz3;

  TFile* fjz4 = new TFile(Form("../ntuples/%s",jz4.Data()));
  TTree* tr4 = (TTree*)fjz4->Get("tr");
  tashka t_jz4;

  int nPtBins;
  double ptBin[30];
  double binWidth[30];
  int nMBins;
  double mBin[40];
  //  double binWidth[60];
  // located in ../unfoldingUtil.h
  // Xbin -> pt
  getXbin(nPtBins, ptBin, 77);
  cout << " nPtBins = " << nPtBins << endl;
  // Ybin -> m2/pt2
  getYbin(nMBins, mBin, 772);

  TH1D* ptBinTemp = new TH1D("ptBinTemp","", nPtBins, ptBin);
  TH1D* mBinTemp = new TH1D("mBinTemp","", nMBins, mBin);

  for (int i = 0; i < nPtBins; ++i) {
    cout << " ptBin " << i << " = " <<  ptBin[i] <<endl;
    binWidth[i] = ptBin[i+1]-ptBin[i];
  }
  for (int i = 0; i < nMBins; ++i) {
    cout << " mBin " << i << " = " <<  mBin[i] <<endl;
    //    binWidth[i] = mBin[i+1]-mBin[i];
  }

  TH2D* hJMS[6];
  TH2D* hJMR[6];
  TH2D* hJES[6];
  TH2D* hJER[6];
  TH1D* hDist[6][20][30]; // reco_m2/pt2 / gen_m2/pT2
  TH1D* hDistJes[6][20][30]; // reco_Pt / gen_pT
  cout << " initializing hdist with etaBins " << nEtaBins << " cent bins " << nCentBins << " pt bins " << nPtBins << " mbins " << nMBins << endl;
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    // note: two extra bins initialized for underflow and overflow
    for ( int ptBin = 0 ; ptBin <= nPtBins+1   ;  ++ptBin) {
      for ( int mBin = 0 ; mBin <= nMBins+1 ; ++mBin) {
	char name[256];
	sprintf(name, "hDistM2Pt2_eta%d_pt%d_m%d",etaBin,ptBin,mBin);
	hDist[etaBin][ptBin][mBin] = new TH1D(name,name, 100,0,2);
	sprintf(name, "hDistPt_eta%d_pt%d_m%d",etaBin,ptBin,mBin);
	hDistJes[etaBin][ptBin][mBin] = new TH1D(name,name, 100,0,2);
      }
    }
  }

  //  TH1D* hPtWidth;
  int cent = 0;
  float eta = 0.;
  long long entries;
  TString fname;
  if ( kSample == kPbPb ) {
    fname = pbpbDataString; 
  }
  else if ( kSample == kPP) {
    fname = ppDataString;
  }
  cout << " setting JMS plots "<< endl;  
  for (int i=0; i < nEtaBins; ++i) { 
    //    for (int j=0; j < nCentBins; ++j) {
      char histName[256];
      sprintf(histName,"hJMS_eta%d_sample%d",i,kSample);
      hJMS[i] = new TH2D(histName,histName, nPtBins, ptBin, nMBins, mBin);
      sprintf(histName,"hJMR_eta%d_sample%d",i,kSample);
      hJMR[i] = new TH2D(histName,histName, nPtBins, ptBin, nMBins, mBin);
      sprintf(histName,"hJER_eta%d_sample%d",i,kSample);
      hJER[i] = new TH2D(histName,histName, nPtBins, ptBin, nMBins, mBin);
      sprintf(histName,"hJES_eta%d_sample%d",i,kSample);
      hJES[i] = new TH2D(histName,histName, nPtBins, ptBin, nMBins, mBin);
      //    }
  }

  TH1D* hFcalReweight;
  if ( kSample == kPbPb ) {
    TFile* fcal = new TFile("reweightFactors/FCal_HP_v_MB_weights.root");
    hFcalReweight = (TH1D*)fcal->Get("weight_MBov_to_HP");
  }
  
  TFile* fData = new TFile(Form("../ntuples/%s",fname.Data()));
  TTree* tree = (TTree*)fData->Get("tr");
  cout << " Setting Tashka for data ..." << endl;
  tashka tr;
  tr.Init(tree);
  entries = tree->GetEntries();

  cout << " Fill MC " << endl;
  for ( int ijz =2 ; ijz<=4 ; ijz++) {
    double jzNorm=0;
    if ( ijz==2)  {
      tr = t_jz2;
      tr.Init(tr2);
      entries = tr2->GetEntries();
      jzNorm = hi9EvtWgtJZ2;
      // jz reweighting global values set in ../jzWeight.h
      // note hi9EvtWgtJZ2 = (JZ2_cs * JZ2_fltEff / JZ2_wNorm) / 5856854
    }

    else if ( ijz==3)  {
      tr = t_jz3;
      tr.Init(tr3);
      entries = tr3->GetEntries();
      jzNorm = hi9EvtWgtJZ3;
    }
    else if ( ijz==4)  {
      tr = t_jz4;
      tr.Init(tr4);
      entries = tr4->GetEntries();
      jzNorm = hi9EvtWgtJZ4;
    }
    cout << "Scanning JZ"<<ijz<<" file.  Total events = " << entries << endl;
    int jetCount = 0;
    for (Int_t i= 0; i<entries ; i++) {
      //  way to cut processing off early ie if statFrac == 0.5, process half the file
      if (ijz == 4 && jetCount == 6209 ) cout << " jetCount " << jetCount << " pt " << tr.jets_ptCalib << " mass2Calib " << tr.jets_mass2Calib << " eta " << tr.jets_yCalib <<  endl;
      if ( i > entries * statFrac ) break;      
      tr.GetEntry(i);
      // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
      //      if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_cent, tr.jets_etaCalib, tr.jets_phiCalib, true) ) // isMC = true
      // passJesEvent has lower pt cut than passEvent to allow truth matching
      //note: THIS IS A HACK: genMass2 should be changed back to genPt with new ttrees
      if ( ! passJesEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib) ) // isMC = true
	continue;
      if (ijz == 4 && jetCount == 6209 ) cout << " passed Jes event " << endl;     
      double fcalWeight;
      if ( kSample==kPbPb) {
	//	cout << " jet fcal " << tr.jets_chMassRcSubt << " rcalreweight " << hFcalReweight->GetXaxis()->FindBin(tr.jets_chMassRcSubt) << endl;
	//	fcalWeight = hFcalReweight->GetBinContent(hFcalReweight->GetXaxis()->FindBin(tr.jets_fcal));
	fcalWeight = hFcalReweight->GetBinContent(hFcalReweight->GetXaxis()->FindBin(tr.jets_fcal));
      }
      else fcalWeight = 1.0;
      cent = tr.jets_cent;
      eta  = TMath::Abs(tr.jets_yCalib);
      int etaBin = -1;
      if (eta < 0.3) etaBin = 0;
      else if (eta >= 0.3 && eta < 0.8) etaBin = 1;
      else if (eta >= 0.8 && eta < 1.2) etaBin = 2;
      else if (eta >= 1.2 && eta < 1.6) etaBin = 3;
      else if (eta >= 1.6 && eta < 2.1) etaBin = 4;   
      else continue;
      int ptBin = ptBinTemp->FindBin(tr.jets_ptCalib);
      int mBin = mBinTemp->FindBin(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib));
      //      cout << " ptBin " << ptBin << " pt " << tr.jets_ptCalib << endl;
      if (eta < 2.1 && kSample == kPP) {
	hDist[5][ptBin][mBin]->Fill (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib), tr.jets_weight * jzNorm * fcalWeight);
	hDistJes[5][ptBin][mBin]->Fill( tr.jets_ptCalib / tr.jets_genPt, tr.jets_weight * jzNorm * fcalWeight);
      }
      hDist[etaBin][ptBin][mBin]->Fill(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib), tr.jets_weight * jzNorm * fcalWeight);      
      hDistJes[etaBin][ptBin][mBin]->Fill( tr.jets_ptCalib / tr.jets_genPt, tr.jets_weight * jzNorm * fcalWeight );

      ++jetCount;
    }

  }
  /*
  cout << " setting ptWidth " << endl;
  hPtWidth = new TH1D("ptWidth","ptWidth",nPtBins+2,ptBin);  
  for (int j=0; j < hDistJes[0][0][0][0]->GetNbinsX(); ++j) {
    hPtWidth->SetBinContent(j+1,binWidth[j]);
    cout << "pt bin " << j << " width " << binWidth[j] << endl;
  }
  */
  
  TCanvas *c0;
  TCanvas *c0b;
  //TFile* foutJES = new TFile(Form("Plots/JES/jes_kSample%d.root",kSample),"recreate");
    //  for (int centBin = 0; centBin < nCentBins; ++centBin) { 

  for ( int etaBin = 0; etaBin < nEtaBins; ++etaBin ) {
    
    c0 = new TCanvas("c0","c0");
    c0->Divide(5,4);
    c0b = new TCanvas("c0b","c0b");
    c0b->Divide(5,4);    
    for ( int ptBin  = 0 ; ptBin < nPtBins ; ptBin++) {
      for ( int mBin = 0 ; mBin < nMBins ; mBin++)  {
	//	  hDistJes[etaBin][centBin][ptBin][mBin]->Divide(hPtWidth);
	double theMean = hDist[etaBin][ptBin][mBin]->GetBinCenter( hDist[etaBin][ptBin][mBin]->GetMaximumBin() );
	double theRms = hDist[etaBin][ptBin][mBin]->GetRMS();

	hDist[etaBin][ptBin][mBin]->Fit("gaus", "","",theMean - 1*theRms, theMean + .5*theRms);
	if ( hDist[etaBin][ptBin][mBin]->GetFunction("gaus") == NULL ) {
	  hJMS[etaBin]->SetBinContent(ptBin+1,mBin+1,-100);
	  hJMS[etaBin]->SetBinError(ptBin+1,mBin+1,-100);
	  hJMR[etaBin]->SetBinContent(ptBin+1,mBin+1,-100);
	  hJMR[etaBin]->SetBinError(ptBin+1,mBin+1,-100);
	}
	else {
	  float theJMS = hDist[etaBin][ptBin][mBin]->GetFunction("gaus")->GetParameter(1);
	  float theJMSerr = hDist[etaBin][ptBin][mBin]->GetFunction("gaus")->GetParError(1);
	  float theJMR = hDist[etaBin][ptBin][mBin]->GetFunction("gaus")->GetParameter(2);
	  float theJMRerr = hDist[etaBin][ptBin][mBin]->GetFunction("gaus")->GetParError(2);
	  hJMS[etaBin]->SetBinContent(ptBin+1,mBin+1,theJMS);
	  hJMS[etaBin]->SetBinError(ptBin+1,mBin+1,theJMSerr);
	  hJMR[etaBin]->SetBinContent(ptBin+1,mBin+1,theJMR);
	  hJMR[etaBin]->SetBinError(ptBin+1,mBin+1,theJMRerr);
	}
	if ( hDistJes[etaBin][ptBin][mBin]->GetFunction("gaus") == NULL ) {
	  hJES[etaBin]->SetBinContent(ptBin+1,mBin+1,-100);
	  hJES[etaBin]->SetBinError(ptBin+1,mBin+1,-100);
	  hJER[etaBin]->SetBinContent(ptBin+1,mBin+1,-100);
	  hJER[etaBin]->SetBinError(ptBin+1,mBin+1,-100);
	}
	else {
	  float theJES = hDistJes[etaBin][ptBin][mBin]->GetFunction("gaus")->GetParameter(1);
	  float theJESerr = hDistJes[etaBin][ptBin][mBin]->GetFunction("gaus")->GetParError(1);
	  float theJER = hDistJes[etaBin][ptBin][mBin]->GetFunction("gaus")->GetParameter(2);
	  float theJERerr = hDistJes[etaBin][ptBin][mBin]->GetFunction("gaus")->GetParError(2);
	  hJES[etaBin]->SetBinContent(ptBin+1,mBin+1,theJES);
	  hJES[etaBin]->SetBinError(ptBin+1,mBin+1,theJESerr);
	  hJER[etaBin]->SetBinContent(ptBin+1,mBin+1,theJER);
	  hJER[etaBin]->SetBinError(ptBin+1,mBin+1,theJERerr);
	}
	hDist[etaBin][ptBin][mBin]->Write();
	hDistJes[etaBin][ptBin][mBin]->Write();

	c0->cd(etaBin+1);
	handsomeTH1(hDist[etaBin][ptBin][mBin],2);
	fixedFontHist(hDist[etaBin][ptBin][mBin],2.5,2,20);
	hDist[etaBin][ptBin][mBin]->Draw();
	hDist[etaBin][ptBin][mBin]->GetXaxis()->SetTitle("(m^{2}/#it{p}_{T}^{2})^{reco} / (m^{2}/#it{p}_{T}^{2})^{truth}");
	drawCentrality(kSample,-1,0.70,0.86,1,24);
	ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
	c0b->cd(etaBin+1);
	handsomeTH1(hDistJes[etaBin][ptBin][mBin],2);
	fixedFontHist(hDistJes[etaBin][ptBin][mBin],2.5,2,20);
	hDistJes[etaBin][ptBin][mBin]->Draw();
	hDistJes[etaBin][ptBin][mBin]->GetXaxis()->SetTitle("(#it{p}_{T}){reco} / (#it{p}_{T})^{truth}");
	drawCentrality(kSample,-1,0.70,0.86,1,24);
	ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
	
      }
    }
    
      char fileName[256];
      sprintf(fileName,"Plots/JES/JMSFit_eta%d_sample%d.png",etaBin,kSample);
      c0->SaveAs(fileName);
      sprintf(fileName,"Plots/JES/JESFit_eta%d_sample_%d.png",etaBin,kSample);
      c0b->SaveAs(fileName);
      delete c0; delete c0b;
    
  }  
    //  }
    foutJES->Close();
    gStyle->SetPalette(kBird);
    TCanvas *c1; TCanvas *c2;
    for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
      //    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      c1 = new TCanvas("c1","",800,400);
      c1->Divide(2);
      c1->cd(1);
      hJMS[etaBin]->GetXaxis()->SetTitle("Truth #it{p}_{T}");
      hJMS[etaBin]->GetYaxis()->SetTitle("Truth m^{2}/#it{p}_{T}^{2}");
      hJMS[etaBin]->Draw("colz");
      hJMS[etaBin]->SetContour(99);
      c1->cd(2);
      hJMR[etaBin]->GetXaxis()->SetTitle("Truth #it{p}_{T}");
      hJMR[etaBin]->GetYaxis()->SetTitle("Truth m^{2}/#it{p}_{T}^{2}");
      hJMR[etaBin]->Draw("colz");
      hJMR[etaBin]->SetContour(99);
      drawCentrality(kSample,-1,0.70,0.86,1,24);
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
      char fileName[256];
      sprintf(fileName,"Plots/JES/ptVsMassJmsJmr_eta%d__sample_%d.png",etaBin,kSample);
      c1->SaveAs(fileName);
      delete c1;
      
      c2 = new TCanvas("c2","",800,400);
      c2->Divide(2);
      c2->cd(1);
      hJES[etaBin]->GetXaxis()->SetTitle("Truth #it{p}_{T}");
      hJES[etaBin]->GetYaxis()->SetTitle("Truth m^{2}/#it{p}_{T}^{2}");
      hJES[etaBin]->Draw("colz");
      hJES[etaBin]->SetContour(99);
      c2->cd(2);
      hJER[etaBin]->GetXaxis()->SetTitle("Truth #it{p}_{T}");
      hJER[etaBin]->GetYaxis()->SetTitle("Truth m^{2}/#it{p}_{T}^{2}");
      hJER[etaBin]->Draw("colz");
      hJER[etaBin]->SetContour(99);
      drawCentrality(kSample,-1,0.70,0.86,1,24);
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
      sprintf(fileName,"Plots/JES/ptVsMassJesJer_eta%d__sample_%d.png",etaBin,kSample);
      c2->SaveAs(fileName);
      delete c2;
      /*
      for ( int ix = lowPtBin+1 ; ix<=highPtBin ; ix++) { 
	TH1D* hs = (TH1D*)hJES->ProjectionY(Form("hs_ix%d",ix),ix,ix);
	handsomeTH1(hs,vColor[ix-lowPtBin]);
	hs->SetXTitle("(m/p_{T})^{Truth}");
	hs->SetYTitle("JES");
	hs->SetAxisRange(0.001,0.239,"X");
	hs->SetAxisRange(0.8,1.4,"Y");
	hs->SetNdivisions(505,"X");
	hs->SetNdivisions(505,"Y");
	if ( ix == lowPtBin+1) hs->Draw();
	else hs->Draw("same");
	leg1s->AddEntry(hs, textBin(xBin,ix,"GeV").Data(),"pl" );
      }
      jumSun(0,1,0.36,1);
      drawCentrality(kSample, icent, 0.2,0.85,1,25);
      drawText("<p_{T}^{Reco}/p_{T}^{Truth}>", 0.2,0.78,1,25);
      leg1s->Draw();
 
      c2s->SaveAs(Form("jms/jes_kSample%d_i-3.pdf",kSample,icent));


      TCanvas* c3a = new TCanvas("c3a", "",500,500);

      for ( int ix = lowPtBin+1 ; ix<=highPtBin ; ix++) { 
	TH1D* hr = (TH1D*)hJER->ProjectionY(Form("hr_ix%d",ix),ix,ix);
	handsomeTH1(hr,vColor[ix-lowPtBin]);
	hr->SetXTitle("(m/p_{T})^{Truth}");
	hr->SetYTitle("JER");
	hr->SetAxisRange(0.001,0.239,"X");
	hr->SetAxisRange(0,0.3,"Y");
	hr->SetNdivisions(505,"X");
	hr->SetNdivisions(505,"Y");
	if ( ix == lowPtBin+1) hr->Draw();
	else hr->Draw("same");
      }
      drawCentrality(kSample, icent, 0.2,0.85,1,25);
      */      
      //    }
  }

  return 0;
}
