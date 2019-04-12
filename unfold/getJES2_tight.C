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
  gStyle->SetPalette(kBird);
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
  double ptBins[30];
  double binWidth[30];
  int nMBins;
  double mBins[40];
  int nMBins_mpt;
  double mBins_mpt[40];
  //  double binWidth[60];
  // located in ../unfoldingUtil.h
  // Xbin -> pt
  getXbin(nPtBins, ptBins, 77);
  cout << " nPtBins = " << nPtBins << endl;
  // Ybin -> m2/pt2
  getYbin(nMBins, mBins, 772);


  TH1D* ptBinTemp = new TH1D("ptBinTemp","", nPtBins, ptBins);
  TH1D* mBinTemp = new TH1D("mBinTemp","", nMBins, mBins);

  for (int i = 0; i < nPtBins; ++i) {
    cout << " ptBin " << i << " = " <<  ptBins[i] <<endl;
    binWidth[i] = ptBins[i+1]-ptBins[i];
  }
  for (int i = 0; i < nMBins; ++i) {
    cout << " mBin " << i << " = " <<  mBins[i] <<endl;
    //    binWidth[i] = mBin[i+1]-mBin[i];
  }

  TH1D* hDist[6][22][32]; // reco_m2/pt2 / gen_m2/pT2
  TH1D* hDistJes[6][22][32]; // reco_Pt / gen_pT
  cout << " initializing hdist with etaBins " << nEtaBins << " cent bins " << nCentBins << " pt bins " << nPtBins << " mbins " << nMBins << endl;
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    // note: two extra bins initialized for underflow and overflow
    for ( int ptBin = 0 ; ptBin <= nPtBins+1   ;  ++ptBin) {
      for ( int mBin = 0 ; mBin <= nMBins+1 ; ++mBin) {
	char name[256];
	sprintf(name, "hDistM2Pt2_eta%d_pt%d_m%d",etaBin,ptBin,mBin);
	hDist[etaBin][ptBin][mBin] = new TH1D(name,name, 100,0,3);
	sprintf(name, "hDistPt_eta%d_pt%d_m%d",etaBin,ptBin,mBin);
	hDistJes[etaBin][ptBin][mBin] = new TH1D(name,name, 100,0,2);
      }
    }
  }
  TH1D* hJMS_eta_1d = new TH1D("JMS_eta","JMS_eta",6,0,6);
  TH1D* hJMS_pt_1d = new TH1D("JMS_pt","JMS_pt",nPtBins,ptBins);
  TH1D* hJMS_m_1d = new TH1D("JMS_m","JMS_m",nMBins,mBins);

  TH1D* hJMR_eta_1d = new TH1D("JMR_eta","JMR_eta",6,0,6);
  TH1D* hJMR_pt_1d = new TH1D("JMR_pt","JMR_pt",nPtBins,ptBins);
  TH1D* hJMR_m_1d = new TH1D("JMR_m","JMR_m",nMBins,mBins);
  TH1D* hJES_eta_1d = new TH1D("JES_eta","JES_eta",6,0,6);
  TH1D* hJES_pt_1d = new TH1D("JES_pt","JES_pt",nPtBins,ptBins);
  TH1D* hJES_m_1d = new TH1D("JES_m","JES_m",nMBins,mBins);
  TH1D* hJER_eta_1d = new TH1D("JER_eta","JER_eta",6,0,6);
  TH1D* hJER_pt_1d = new TH1D("JER_pt","JER_pt",nPtBins,ptBins);
  TH1D* hJER_m_1d = new TH1D("JER_m","JER_m",nMBins,mBins);
  
  TH2D* hJMS[6];
  TH2D* hJMR[6];
  TH2D* hJES[6];
  TH2D* hJER[6];
  TH1D* hDist_eta[6]; // reco_m2/pt2 / gen_m2/pT2
  TH1D* hDistJes_eta[6]; // reco_Pt / gen_pT
  TH1D* hDist_pt[22]; // reco_m2/pt2 / gen_m2/pT2
  TH1D* hDistJes_pt[22]; // reco_Pt / gen_pT
  TH1D* hDist_m[32]; // reco_m2/pt2 / gen_m2/pT2
  TH1D* hDistJes_m[32]; // reco_Pt / gen_pT

  cout << " initializing hdist with etaBins " << nEtaBins << " cent bins " << nCentBins << " pt bins " << nPtBins << " mbins " << nMBins << endl;
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {

    char name[256];
    sprintf(name, "hDistM2Pt2_eta%d",etaBin);
    hDist_eta[etaBin] = new TH1D(name,name, 100,0,3);
    sprintf(name, "hDistPt_eta%d",etaBin);
    hDistJes_eta[etaBin] = new TH1D(name,name, 100,0,2);

    char histName[256];
    sprintf(histName,"hJMS_eta%d_sample%d",etaBin,kSample);
    hJMS[etaBin] = new TH2D(histName,histName, nPtBins, ptBins, nMBins, mBins);
    sprintf(histName,"hJMR_eta%d_sample%d",etaBin,kSample);
    hJMR[etaBin] = new TH2D(histName,histName, nPtBins, ptBins, nMBins, mBins);
    sprintf(histName,"hJER_eta%d_sample%d",etaBin,kSample);
    hJER[etaBin] = new TH2D(histName,histName, nPtBins, ptBins, nMBins, mBins);
    sprintf(histName,"hJES_eta%d_sample%d",etaBin,kSample);
    hJES[etaBin] = new TH2D(histName,histName, nPtBins, ptBins, nMBins, mBins);
  }
  for (int mBin = 0; mBin <= nMBins+1; ++mBin) {
    // note: two extra bins initialized for underflow and overflow    
    char name[256];
    sprintf(name, "hDistM2Pt2_m%d",mBin);
    hDist_m[mBin] = new TH1D(name,name, 100,0,3);
    sprintf(name, "hDistm_m%d",mBin);
    hDistJes_m[mBin] = new TH1D(name,name, 100,0,2);
  }

  for (int ptBin = 0; ptBin <= nPtBins+1; ++ptBin) {
    // note: two extra bins initialized for underflow and overflow    
    char name[256];
    sprintf(name, "hDistPt2_pt%d",ptBin);
    hDist_pt[ptBin] = new TH1D(name,name, 100,0,3);
    sprintf(name, "hDistPt_pt%d",ptBin);
    hDistJes_pt[ptBin] = new TH1D(name,name, 100,0,2);

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
      int ptBin = ptBinTemp->FindBin(tr.jets_genPt);
      int mBin = mBinTemp->FindBin((tr.jets_genMass*tr.jets_genMass)/(tr.jets_genPt*tr.jets_genPt));

      //      cout << " ptBin " << ptBin << " pt " << tr.jets_ptCalib << endl;
      if (eta < 2.1 && kSample == kPP) {
	hDist_eta[5]->Fill ((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)) / ( tr.jets_genMass2 / (tr.jets_genPt * tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight );
	hDistJes_eta[5]->Fill( tr.jets_ptCalib / tr.jets_genPt , tr.jets_weight * jzNorm * fcalWeight );

	hDist[5][ptBin][mBin]->Fill ((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)) / (  tr.jets_genMass2 / (tr.jets_genPt * tr.jets_genPt) ) , tr.jets_weight * jzNorm * fcalWeight );
	hDistJes[5][ptBin][mBin]->Fill( tr.jets_ptCalib / tr.jets_genPt, tr.jets_weight * jzNorm * fcalWeight );
      }
      hDist_eta[etaBin]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)) / (  tr.jets_genMass2 / (tr.jets_genPt * tr.jets_genPt) ), tr.jets_weight * jzNorm * fcalWeight );      
      hDistJes_eta[etaBin]->Fill( tr.jets_ptCalib / tr.jets_genPt , tr.jets_weight * jzNorm * fcalWeight  );
      
      hDist_pt[ptBin]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)) / (  tr.jets_genMass2 / (tr.jets_genPt * tr.jets_genPt) ) , tr.jets_weight * jzNorm * fcalWeight );      
      hDistJes_pt[ptBin]->Fill( tr.jets_ptCalib / tr.jets_genPt  , tr.jets_weight * jzNorm * fcalWeight );
      
      hDist_m[mBin]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)) / (  tr.jets_genMass2 / (tr.jets_genPt * tr.jets_genPt) ) , tr.jets_weight * jzNorm * fcalWeight );
      hDistJes_m[mBin]->Fill( tr.jets_ptCalib / tr.jets_genPt, tr.jets_weight * jzNorm * fcalWeight  );      
      
      hDist[etaBin][ptBin][mBin]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)) / (  tr.jets_genMass2 / (tr.jets_genPt * tr.jets_genPt) ) , tr.jets_weight * jzNorm * fcalWeight );      
      hDistJes[etaBin][ptBin][mBin]->Fill( tr.jets_ptCalib / tr.jets_genPt , tr.jets_weight * jzNorm * fcalWeight );
      ++jetCount;
    }

  }
  
    //TFile* foutJES = new TFile(Form("Plots/JES/jes_kSample%d.root",kSample),"recreate");
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(111);
 
  TCanvas *c0 = new TCanvas("c0","c0");
  c0->Divide(6,1);
 
  TCanvas *c_b = new TCanvas("c_b","c_b");
  c_b->Divide(6,1);
  cout << " looping through eta " << endl;
  for ( int etaBin = 0; etaBin < nEtaBins; ++etaBin ) {    
    
    c0->cd(etaBin+1);
    
    double theMean = hDist_eta[etaBin]->GetBinCenter( hDist_eta[etaBin]->GetMaximumBin() );
    double theRms = hDist_eta[etaBin]->GetRMS();
    cout << "eta " << etaBin << " mean " << theMean << endl;
    hDist_eta[etaBin]->Fit("gaus", "","",theMean - 1*theRms, theMean + .5*theRms);
   
    if ( hDist_eta[etaBin]->GetFunction("gaus") == NULL ) {

      hJMS_eta_1d->SetBinContent(etaBin+1,-100);
      hJMS_eta_1d->SetBinError(etaBin+1,-100);
      hJMR_eta_1d->SetBinContent(etaBin+1,-100);
      hJMR_eta_1d->SetBinError(etaBin+1,-100);
    }
    else {
      float theJMS = hDist_eta[etaBin]->GetFunction("gaus")->GetParameter(1);
      float theJMSerr = hDist_eta[etaBin]->GetFunction("gaus")->GetParError(1);
      float theJMR = hDist_eta[etaBin]->GetFunction("gaus")->GetParameter(2);
      float theJMRerr = hDist_eta[etaBin]->GetFunction("gaus")->GetParError(2);

      hJMS_eta_1d->SetBinContent(etaBin+1,theJMS);
      hJMS_eta_1d->SetBinError(etaBin+1,theJMSerr);
      hJMR_eta_1d->SetBinContent(etaBin+1,theJMR);
      hJMR_eta_1d->SetBinError(etaBin+1,theJMRerr);
      cout << "eta " << etaBin << " JMS " << theJMS << " JMR " << theJMR << endl;
    }
    
    handsomeTH1(hDist_eta[etaBin],1);
    fixedFontHist(hDist_eta[etaBin],2.5,2,20);
    hDist_eta[etaBin]->Draw();

    hDist_eta[etaBin]->GetXaxis()->SetTitle("(m^{2}/#it{p}_{T}^{2})^{reco} / (m^{2}/#it{p}_{T}^{2})^{truth}");
    if (etaBin == 1)    drawCentrality(kSample,-1,0.65,0.7,1,24);
    ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    cout << " drawing histo eta bin " << etaBin << endl;
    
    c_b->cd(etaBin+1);
    double theMeanJes = hDistJes_eta[etaBin]->GetBinCenter( hDistJes_eta[etaBin]->GetMaximumBin() );
    double theRmsJes = hDistJes_eta[etaBin]->GetRMS();
    cout << "eta " << etaBin << " meanJes " << theMeanJes << endl;
    hDistJes_eta[etaBin]->Fit("gaus", "","",theMeanJes - 1*theRmsJes, theMeanJes + .5*theRmsJes);
 
    if ( hDistJes_eta[etaBin]->GetFunction("gaus") == NULL ) {

      hJES_eta_1d->SetBinContent(etaBin+1,-100);
      hJES_eta_1d->SetBinError(etaBin+1,-100);
      hJER_eta_1d->SetBinContent(etaBin+1,-100);
      hJER_eta_1d->SetBinError(etaBin+1,-100);
    }
    else {
      float theJES = hDistJes_eta[etaBin]->GetFunction("gaus")->GetParameter(1);
      float theJESerr = hDistJes_eta[etaBin]->GetFunction("gaus")->GetParError(1);
      float theJER = hDistJes_eta[etaBin]->GetFunction("gaus")->GetParameter(2);
      float theJERerr = hDistJes_eta[etaBin]->GetFunction("gaus")->GetParError(2);
      cout << "eta " << etaBin << " JES " << theJES << " JER " << theJER << endl;
      hJES_eta_1d->SetBinContent(etaBin+1,theJES);
      hJES_eta_1d->SetBinError(etaBin+1,theJESerr);
      hJER_eta_1d->SetBinContent(etaBin+1,theJER);
      hJER_eta_1d->SetBinError(etaBin+1,theJERerr);
    }
    
    
    handsomeTH1(hDistJes_eta[etaBin],1);
    //    fixedFontHist(hDistJes_eta[etaBin],2.5,2,20);
    if (etaBin == 0) {
      drawCentrality(kSample,-1,0.70,0.86,1,24);
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    } 
    hDistJes_eta[etaBin]->Draw();
    hDistJes_eta[etaBin]->GetXaxis()->SetTitle("(#it{p}_{T}){reco} / (#it{p}_{T})^{truth}");

      
  }
  char fileName[256];
  
  sprintf(fileName,"Plots/JES/JMSFit_eta_sample%d.png",kSample);
  c0->SaveAs(fileName);

  sprintf(fileName,"Plots/JES/JESFit_eta_sample_%d.png",kSample);
  c_b->SaveAs(fileName);
  delete c0;
  delete c_b;
  
  //   return 0;
  
  TCanvas *c1;
  TCanvas *c1b;
  c1 = new TCanvas("c1","c1");
  c1->Divide(5,4);
  c1b = new TCanvas("c1b","c1b");
  c1b->Divide(5,4);        
  for ( int ptBin = 0; ptBin < nPtBins; ++ptBin ) {    
    c1->cd(ptBin+1);
    double theMean = hDist_pt[ptBin]->GetBinCenter( hDist_pt[ptBin]->GetMaximumBin() );
    double theRms = hDist_pt[ptBin]->GetRMS();

    hDist_pt[ptBin]->Fit("gaus", "","",theMean - 1*theRms, theMean + .5*theRms);

    c1b->cd(ptBin+1);
    double theMeanJes = hDistJes_pt[ptBin]->GetBinCenter( hDistJes_pt[ptBin]->GetMaximumBin() );
    double theRmsJes = hDistJes_pt[ptBin]->GetRMS();
    cout << "pt " << ptBin << " meanJes " << theMeanJes << endl;
    hDistJes_pt[ptBin]->Fit("gaus", "","",theMeanJes - 1*theRmsJes, theMeanJes + .5*theRmsJes);
    c1->cd(ptBin+1);
    if ( hDist_pt[ptBin]->GetFunction("gaus") == NULL ) {
      hJMS_pt_1d->SetBinContent(ptBin+1,-100);
      hJMS_pt_1d->SetBinError(ptBin+1,-100);
      hJMR_pt_1d->SetBinContent(ptBin+1,-100);
      hJMR_pt_1d->SetBinError(ptBin+1,-100);
    }
    else {
      float theJMS = hDist_pt[ptBin]->GetFunction("gaus")->GetParameter(1);
      float theJMSerr = hDist_pt[ptBin]->GetFunction("gaus")->GetParError(1);
      float theJMR = hDist_pt[ptBin]->GetFunction("gaus")->GetParameter(2);
      float theJMRerr = hDist_pt[ptBin]->GetFunction("gaus")->GetParError(2);
      hJMS_pt_1d->SetBinContent(ptBin+1,theJMS);
      hJMS_pt_1d->SetBinError(ptBin+1,theJMSerr);
      hJMR_pt_1d->SetBinContent(ptBin+1,theJMR);
      hJMR_pt_1d->SetBinError(ptBin+1,theJMRerr);
      cout << " pt " << ptBin << " JMR " << theJMR << " JMS " << theJMS << endl;
    }
    c1b->cd(ptBin+1);
    if ( hDistJes_pt[ptBin]->GetFunction("gaus") == NULL ) {
      hJES_pt_1d->SetBinContent(ptBin+1,-100);
      hJES_pt_1d->SetBinError(ptBin+1,-100);
      hJER_pt_1d->SetBinContent(ptBin+1,-100);
      hJER_pt_1d->SetBinError(ptBin+1,-100);
    }
    else {
      float theJES = hDistJes_pt[ptBin]->GetFunction("gaus")->GetParameter(1);
      float theJESerr = hDistJes_pt[ptBin]->GetFunction("gaus")->GetParError(1);
      float theJER = hDistJes_pt[ptBin]->GetFunction("gaus")->GetParameter(2);
      float theJERerr = hDistJes_pt[ptBin]->GetFunction("gaus")->GetParError(2);
      hJES_pt_1d->SetBinContent(ptBin+1,theJES);
      hJES_pt_1d->SetBinError(ptBin+1,theJESerr);
      hJER_pt_1d->SetBinContent(ptBin+1,theJER);
      hJER_pt_1d->SetBinError(ptBin+1,theJERerr);
      cout << " pt " << ptBin << " JER " << theJER << " JES " << theJES << endl;
    }

    c1->cd(ptBin+1);
    handsomeTH1(hDist_pt[ptBin],1);
    //    fixedFontHist(hDist_pt[ptBin],2.5,2,20);
    hDist_pt[ptBin]->Draw();

    hDist_pt[ptBin]->GetXaxis()->SetTitle("(m^{2}/#it{p}_{T}^{2})^{reco} / (m^{2}/#it{p}_{T}^{2})^{truth}");
    if (ptBin == 0) {
      drawCentrality(kSample,-1,0.70,0.86,1,24);
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    }
    c1b->cd(ptBin+1);
    handsomeTH1(hDistJes_pt[ptBin],1);
    //    fixedFontHist(hDistJes_pt[ptBin],2.5,2,20);
    hDistJes_pt[ptBin]->Draw();
    hDistJes_pt[ptBin]->GetXaxis()->SetTitle("(#it{p}_{T}){reco} / (#it{p}_{T})^{truth}");
    if (ptBin == 0) {
      drawCentrality(kSample,-1,0.70,0.86,1,24);
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    }
  }

  sprintf(fileName,"Plots/JES/JMSFit_pt_sample%d.png",kSample);
  c1->SaveAs(fileName);
  sprintf(fileName,"Plots/JES/JESFit_pt_sample_%d.png",kSample);
  c1b->SaveAs(fileName);
  delete c1; delete c1b;

  TCanvas *c2;
  TCanvas *c2b;
  c2 = new TCanvas("c2","c2");
  c2->Divide(6,5);
  c2b = new TCanvas("c2b","c2b");
  c2b->Divide(6,5);        
  for ( int mBin = 0; mBin < nMBins; ++mBin ) {    
    c2->cd(mBin+1); 
    double theMean = hDist_m[mBin]->GetBinCenter( hDist_m[mBin]->GetMaximumBin() );
    double theRms = hDist_m[mBin]->GetRMS();
    hDist_m[mBin]->Fit("gaus", "","",theMean - 1*theRms, theMean + .5*theRms);
    if ( hDist_m[mBin]->GetFunction("gaus") == NULL ) {
      hJMS_m_1d->SetBinContent(mBin+1,-100);
      hJMS_m_1d->SetBinError(mBin+1,-100);
      hJMR_m_1d->SetBinContent(mBin+1,-100);
      hJMR_m_1d->SetBinError(mBin+1,-100);
    }
    else {
      float theJMS = hDist_m[mBin]->GetFunction("gaus")->GetParameter(1);
      float theJMSerr = hDist_m[mBin]->GetFunction("gaus")->GetParError(1);
      float theJMR = hDist_m[mBin]->GetFunction("gaus")->GetParameter(2);
      float theJMRerr = hDist_m[mBin]->GetFunction("gaus")->GetParError(2);
      hJMS_m_1d->SetBinContent(mBin+1,theJMS);
      hJMS_m_1d->SetBinError(mBin+1,theJMSerr);
      hJMR_m_1d->SetBinContent(mBin+1,theJMR);
      hJMR_m_1d->SetBinError(mBin+1,theJMRerr);
    }
    c2b->cd(mBin+1); 
    double theMeanJes = hDist_m[mBin]->GetBinCenter( hDist_m[mBin]->GetMaximumBin() );
    double theRmsJes = hDist_m[mBin]->GetRMS();
    hDistJes_m[mBin]->Fit("gaus", "","",theMeanJes - 1*theRmsJes, theMeanJes + .5*theRmsJes);
    if ( hDistJes_m[mBin]->GetFunction("gaus") == NULL ) {
      hJES_m_1d->SetBinContent(mBin+1,-100);
      hJES_m_1d->SetBinError(mBin+1,-100);
      hJER_m_1d->SetBinContent(mBin+1,-100);
      hJER_m_1d->SetBinError(mBin+1,-100);
    }
    else {
      float theJES = hDistJes_m[mBin]->GetFunction("gaus")->GetParameter(1);
      float theJESerr = hDistJes_m[mBin]->GetFunction("gaus")->GetParError(1);
      float theJER = hDistJes_m[mBin]->GetFunction("gaus")->GetParameter(2);
      float theJERerr = hDistJes_m[mBin]->GetFunction("gaus")->GetParError(2);
      hJES_m_1d->SetBinContent(mBin+1,theJES);
      hJES_m_1d->SetBinError(mBin+1,theJESerr);
      hJER_m_1d->SetBinContent(mBin+1,theJER);
      hJER_m_1d->SetBinError(mBin+1,theJERerr);
    }
    
    c2->cd(mBin+1);
    handsomeTH1(hDist_m[mBin],1);
    //    fixedFontHist(hDist_m[mBin],2.5,2,20);
    hDist_m[mBin]->Draw();
    hDist_m[mBin]->GetXaxis()->SetTitle("(m^{2}/#it{p}_{T}^{2})^{reco} / (m^{2}/#it{p}_{T}^{2})^{truth}");
    if (mBin == 0) {
      drawCentrality(kSample,-1,0.70,0.86,1,24);
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    }
    c2b->cd(mBin+1);
    handsomeTH1(hDistJes_m[mBin],1);
    //    fixedFontHist(hDistJes_m[mBin],2.5,2,20);
    hDistJes_m[mBin]->Draw();
    hDistJes_m[mBin]->GetXaxis()->SetTitle("(#it{p}_{T}){reco} / (#it{p}_{T})^{truth}");
    if (mBin == 0) {
      drawCentrality(kSample,-1,0.70,0.86,1,24);
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    }
  }

  sprintf(fileName,"Plots/JES/JMSFit_m_sample%d.png",kSample);
  c2->SaveAs(fileName);
  sprintf(fileName,"Plots/JES/JESFit_m_sample_%d.png",kSample);
  c2b->SaveAs(fileName);
  delete c2; delete c2b;
  

  for ( int etaBin = 0; etaBin < nEtaBins; ++etaBin ) {    
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

	double theMeanJes = hDistJes[etaBin][ptBin][mBin]->GetBinCenter( hDistJes[etaBin][ptBin][mBin]->GetMaximumBin() );
	double theRmsJes = hDistJes[etaBin][ptBin][mBin]->GetRMS();

	hDistJes[etaBin][ptBin][mBin]->Fit("gaus", "","",theMeanJes - 1*theRmsJes, theMeanJes + .5*theRmsJes);
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
      }
    }    
  }  



  TCanvas *c4 = new TCanvas("c4","");
  handsomeTH1(hJMS_eta_1d,1);
  handsomeTH1(hJMR_eta_1d,2);
  //  fixedFontHist(hJMS_eta_1d,2.5,2,20);
  //  fixedFontHist(hJMR_eta_1d,2.5,2,20);
  hJMS_eta_1d->GetXaxis()->SetTitle("#eta Bin");
  hJMS_eta_1d->Draw();
  hJMS_eta_1d->GetYaxis()->SetRangeUser(0,1.4);
  hJMR_eta_1d->Draw("same");
  hJMR_eta_1d->GetYaxis()->SetRangeUser(0,1.4);
  drawCentrality(kSample,-1,0.70,0.86,1,24);
  ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
  TLegend *leg1 = new TLegend(0.5271854,0.4394175,0.9049639,0.6533743,NULL,"brNDC");      
  easyLeg(leg1,"");
  leg1->AddEntry(hJMS_eta_1d, "JMS","plfe");
  leg1->AddEntry(hJMR_eta_1d, "JMR","plfe");
  leg1->Draw();
  sprintf(fileName,"Plots/JES/JMS_eta_sample%d.png",kSample);
  c4->SaveAs(fileName);
  delete c4;
      
  TCanvas *c5 = new TCanvas("c5","");
  handsomeTH1(hJES_eta_1d,1);
  handsomeTH1(hJER_eta_1d,2);
  //  fixedFontHist(hJES_eta_1d,2.5,2,20);
  //fixedFontHist(hJER_eta_1d,2.5,2,20);
  hJES_eta_1d->GetXaxis()->SetTitle("#eta Bin");
  hJES_eta_1d->Draw();
  hJES_eta_1d->GetYaxis()->SetRangeUser(0,1.2);
  hJER_eta_1d->Draw("same");
  drawCentrality(kSample,-1,0.70,0.86,1,24);
  ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
  TLegend *leg2 = new TLegend(0.5271854,0.4394175,0.9049639,0.6533743,NULL,"brNDC");      
  easyLeg(leg2,"");
  leg2->AddEntry(hJES_eta_1d, "JES","plfe");
  leg2->AddEntry(hJER_eta_1d, "JER","plfe");
  leg2->Draw();
  sprintf(fileName,"Plots/JES/JES_eta_sample%d.png",kSample);
  c5->SaveAs(fileName);
  delete c5;

  TCanvas *c6 = new TCanvas("c6","");
  handsomeTH1(hJMS_pt_1d,1);
  handsomeTH1(hJMR_pt_1d,2);
  //fixedFontHist(hJMS_pt_1d,2.5,2,20);
  //  fixedFontHist(hJMR_pt_1d,2.5,2,20);
  hJMS_pt_1d->Draw();
  hJMS_pt_1d->GetYaxis()->SetRangeUser(0,1.4);
  hJMS_pt_1d->GetXaxis()->SetTitle("Truth #it{p}_{T} (GeV)");
  hJMR_pt_1d->Draw("same");
  //  hJMS_pt_1d->GetYaxis()->SetRangeUser(0,1.2);

  drawCentrality(kSample,-1,0.70,0.86,1,24);
  ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
  TLegend *leg6 = new TLegend(0.5271854,0.4394175,0.9049639,0.6533743,NULL,"brNDC");      
  easyLeg(leg6,"");
  leg6->AddEntry(hJMS_pt_1d, "JMS","plfe");
  leg6->AddEntry(hJMR_pt_1d, "JMR","plfe");
  leg6->Draw();
  sprintf(fileName,"Plots/JES/JMS_pt_sample%d.png",kSample);
  c6->SaveAs(fileName);
  delete c6;
      
  TCanvas *c7 = new TCanvas("c7","");
  handsomeTH1(hJES_pt_1d,1);
  handsomeTH1(hJER_pt_1d,2);
  //  fixedFontHist(hJES_pt_1d,2.5,2,20);
  //  fixedFontHist(hJER_pt_1d,2.5,2,20);
  hJES_pt_1d->GetXaxis()->SetTitle("Truth #it{p}_{T} (GeV)");
  hJES_pt_1d->Draw();
  hJES_pt_1d->GetYaxis()->SetRangeUser(0,1.2);
  hJER_pt_1d->Draw("same");
  drawCentrality(kSample,-1,0.70,0.86,1,24);
  ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
  TLegend *leg3 = new TLegend(0.5271854,0.4394175,0.9049639,0.6533743,NULL,"brNDC");      
  easyLeg(leg3,"");
  leg3->AddEntry(hJES_pt_1d, "JES","plfe");
  leg3->AddEntry(hJER_pt_1d, "JER","plfe");
  leg3->Draw();
  sprintf(fileName,"Plots/JES/JES_pt_sample%d.png",kSample);
  c7->SaveAs(fileName);
  delete c7;

  TCanvas *c8 = new TCanvas("c8","");
  handsomeTH1(hJMS_m_1d,1);
  handsomeTH1(hJMR_m_1d,2);
  //  fixedFontHist(hJMS_m_1d,2.5,2,20);
  //fixedFontHist(hJMR_m_1d,2.5,2,20);
  hJMS_m_1d->GetXaxis()->SetTitle("Truth m^{2}/#it{p}_{T}^{2}");
  hJMS_m_1d->Draw();
  hJMS_m_1d->GetYaxis()->SetRangeUser(0,1.4);
  hJMR_m_1d->Draw("same");
  hJMR_m_1d->GetYaxis()->SetRangeUser(0,1.4);
  drawCentrality(kSample,-1,0.70,0.86,1,24);
  ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
  TLegend *leg4 = new TLegend(0.5271854,0.4394175,0.9049639,0.6533743,NULL,"brNDC");      
  easyLeg(leg4,"");
  leg4->AddEntry(hJMS_m_1d, "JMS","plfe");
  leg4->AddEntry(hJMR_m_1d, "JMR","plfe");
  leg4->Draw();
  sprintf(fileName,"Plots/JES/JMS_m_sample%d.png",kSample);
  c8->SaveAs(fileName);
  delete c8;
  
  TCanvas *c9 = new TCanvas("c9","");
  handsomeTH1(hJES_m_1d,1);
  handsomeTH1(hJER_m_1d,2);
  //  fixedFontHist(hJES_m_1d,2.5,2,20);
  // fixedFontHist(hJER_m_1d,2.5,2,20);
  hJES_m_1d->GetXaxis()->SetTitle("Truth m^{2}/#it{p}_{T}^{2}");
  hJES_m_1d->Draw();
  hJES_m_1d->GetYaxis()->SetRangeUser(0,1.2);
  hJER_m_1d->Draw("same");
  drawCentrality(kSample,-1,0.70,0.86,1,24);
  ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
  TLegend *leg5 = new TLegend(0.5271854,0.4394175,0.9049639,0.6533743,NULL,"brNDC");      
  easyLeg(leg5,"");
  leg5->AddEntry(hJES_eta_1d, "JES","plfe");
  leg5->AddEntry(hJER_eta_1d, "JER","plfe");
  leg5->Draw();
  sprintf(fileName,"Plots/JES/JES_m_sample%d.png",kSample);
  c9->SaveAs(fileName);
  delete c9;

  TCanvas *c10; TCanvas *c11; TCanvas *c12; TCanvas *c13;
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    c10 = new TCanvas("c10","");
    //    handsomeTH2(hJMS[etaBin],2);
    hJMS[etaBin]->GetXaxis()->SetTitle("Truth #it{p}_{T}");
    hJMS[etaBin]->GetYaxis()->SetTitle("Truth m^{2}/#it{p}^{2}_{T}");
    hJMS[etaBin]->Draw("colz");
    hJMS[etaBin]->SetContour(99);
    hJMS[etaBin]->SetMaximum(1.4);
    hJMS[etaBin]->SetMinimum(0.8);
    drawCentrality(kSample,-1,0.70,0.86,1,24);
    ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    sprintf(fileName,"Plots/JES/JMS_MvsPt_sample%d.png",kSample);
    c10->SaveAs(fileName);
    delete c10;

    c11 = new TCanvas("c11","");
    //    handsomeTH2(hJMR[etaBin],2);
    hJMR[etaBin]->GetXaxis()->SetTitle("Truth #it{p}_{T}");
    hJMR[etaBin]->GetYaxis()->SetTitle("Truth m^{2}/#it{p}^{2}_{T}");
    hJMR[etaBin]->Draw("colz");
    hJMR[etaBin]->SetContour(99);
    hJMR[etaBin]->SetMaximum(0.8);
    hJMR[etaBin]->SetMinimum(0.);
    drawCentrality(kSample,-1,0.70,0.86,1,24);
    ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    sprintf(fileName,"Plots/JES/JMR_MvsPt_sample%d.png",kSample);
    c11->SaveAs(fileName);
    delete c11;

    c12 = new TCanvas("c12","");
    //    handsomeTH2(hJES[etaBin],2);
    hJES[etaBin]->GetXaxis()->SetTitle("Truth #it{p}_{T}");
    hJES[etaBin]->GetYaxis()->SetTitle("Truth m^{2}/#it{p}^{2}_{T}");
    hJES[etaBin]->Draw("colz");
    hJES[etaBin]->SetContour(99);
    hJES[etaBin]->SetMaximum(1.1);
    hJES[etaBin]->SetMinimum(0.8);
    drawCentrality(kSample,-1,0.70,0.86,1,24);
    ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    sprintf(fileName,"Plots/JES/JES_MvsPt_sample%d.png",kSample);
    c12->SaveAs(fileName);
    delete c12;

    c13 = new TCanvas("c13","");
    //    handsomeTH2(hJER[etaBin],2);
    hJER[etaBin]->GetXaxis()->SetTitle("Truth #it{p}_{T}");
    hJER[etaBin]->GetYaxis()->SetTitle("Truth m^{2}/#it{p}^{2}_{T}");
    hJER[etaBin]->Draw("colz");
    hJER[etaBin]->SetMaximum(0.2);
    hJER[etaBin]->SetMinimum(0.0);
    hJER[etaBin]->SetContour(99);
    drawCentrality(kSample,-1,0.70,0.86,1,24);
    ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    sprintf(fileName,"Plots/JES/JER_MvsPt_sample%d.png",kSample);
    c13->SaveAs(fileName);
    delete c13;
  }
  return 0;
}
