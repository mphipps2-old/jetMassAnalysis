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
  char fileName[256];
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
  getYbin(nMBins_mpt, mBins_mpt, 78);

  TH1D* ptBinTemp = new TH1D("ptBinTemp","", nPtBins, ptBins);
  TH1D* mBinTemp = new TH1D("mBinTemp","", nMBins, mBins);
  TH1D* mBinTemp_mpt = new TH1D("mBinTemp_mpt","", nMBins_mpt, mBins_mpt);

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
	hDist[etaBin][ptBin][mBin] = new TH1D(name,name, 100,-2,5);
	sprintf(name, "hDistPt_eta%d_pt%d_m%d",etaBin,ptBin,mBin);
	hDistJes[etaBin][ptBin][mBin] = new TH1D(name,name, 100,0,2);
      }
    }
  }

  TH1D* hJMS_pt1_m_1d = new TH1D("JMS_pt1_mpt","JMS_pt1_mpt",nMBins_mpt,mBins_mpt);
  TH1D* hJMS_pt2_m_1d = new TH1D("JMS_pt2_mpt","JMS_pt2_mpt",nMBins_mpt,mBins_mpt);
  TH1D* hJMS_pt3_m_1d = new TH1D("JMS_pt3_mpt","JMS_pt3_mpt",nMBins_mpt,mBins_mpt);

  TH1D* hJMR_pt1_m_1d = new TH1D("JMR_pt1_mpt","JMR_pt1_mpt",nMBins_mpt,mBins_mpt);
  TH1D* hJMR_pt2_m_1d = new TH1D("JMR_pt2_mpt","JMR_pt2_mpt",nMBins_mpt,mBins_mpt);
  TH1D* hJMR_pt3_m_1d = new TH1D("JMR_pt3_mpt","JMR_pt3_mpt",nMBins_mpt,mBins_mpt);
 
  TH1D* hDistRecoPt1_m[32]; // reco_m2/pt2, pt range 1
  TH1D* hDistRecoPt2_m[32]; // reco_m2/pt2, pt range 2
  TH1D* hDistRecoPt3_m[32]; // reco_m2/pt2, pt range s3

  // to match yongsuns binning
  for (int mBin_mpt = 0; mBin_mpt <= nMBins_mpt+1; ++mBin_mpt) {
    // note: two extra bins initialized for underflow and overflow
    if (mBin_mpt > 1) {
      char name[256];
      sprintf(name, "hDistRecoPt1_m%d",mBin_mpt);
      hDistRecoPt1_m[mBin_mpt] = new TH1D(name,name, 100,0,0.05);
      sprintf(name, "hDistRecoPt2_m%d",mBin_mpt);
      hDistRecoPt2_m[mBin_mpt] = new TH1D(name,name, 100,0.,0.05);
      sprintf(name, "hDistRecoP3_m%d",mBin_mpt);
      hDistRecoPt3_m[mBin_mpt] = new TH1D(name,name, 100,0.,0.05);
    }
    else if (mBin_mpt == 1) {
      char name[256];
      sprintf(name, "hDistRecoPt1_m%d",mBin_mpt);
      hDistRecoPt1_m[mBin_mpt] = new TH1D(name,name, 100,0,0.02);
      sprintf(name, "hDistRecoPt2_m%d",mBin_mpt);
      hDistRecoPt2_m[mBin_mpt] = new TH1D(name,name, 100,0.,0.02);
      sprintf(name, "hDistRecoP3_m%d",mBin_mpt);
      hDistRecoPt3_m[mBin_mpt] = new TH1D(name,name, 100,0.,0.02);
    }
      else  {
      char name[256];
      sprintf(name, "hDistRecoPt1_m%d",mBin_mpt);
      hDistRecoPt1_m[mBin_mpt] = new TH1D(name,name, 100,0,0.005);
      sprintf(name, "hDistRecoPt2_m%d",mBin_mpt);
      hDistRecoPt2_m[mBin_mpt] = new TH1D(name,name, 100,0.,0.005);
      sprintf(name, "hDistRecoP3_m%d",mBin_mpt);
      hDistRecoPt3_m[mBin_mpt] = new TH1D(name,name, 100,0.,0.005);
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
      int mBin_mpt = mBinTemp_mpt->FindBin(tr.jets_genMass/(tr.jets_genPt));
      //      cout << " ptBin " << ptBin << " pt " << tr.jets_ptCalib << endl;

      if (ptBin == 1 || ptBin == 2) {
	//	if (mBin_mpt < 3) cout << " pt " << tr.jets_ptCalib << " m2/pt2 " << tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib) << " m/pt bin " << mBin_mpt << endl;
	hDistRecoPt1_m[mBin_mpt]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)) );
	//		hDistRecoPt1_m[mBin_mpt]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight );
      }
      else if (ptBin == 5 || ptBin == 6) {
	hDistRecoPt2_m[mBin_mpt]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)) );
	//	hDistRecoPt2_m[mBin_mpt]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight );
      }
      else if (ptBin == 13 || ptBin == 14) {
	hDistRecoPt3_m[mBin_mpt]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)) );
	//	hDistRecoPt3_m[mBin_mpt]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight );
      }
      
      ++jetCount;
    }

  }
 
  TCanvas *c1;
  c1 = new TCanvas("c1","c1");
  c1->Divide(4,2);
   
  TCanvas *c2;
  c2 = new TCanvas("c2","c2");
  c2->Divide(4,2);
   
  TCanvas *c3;
  c3 = new TCanvas("c3","c3");
  c3->Divide(4,2);

  for ( int mBin_mpt = 0; mBin_mpt < nMBins_mpt; ++mBin_mpt ) {
    c1->cd(mBin_mpt+1);

   
    double theMeanJmsPt1 = hDistRecoPt1_m[mBin_mpt]->GetBinCenter( hDistRecoPt1_m[mBin_mpt]->GetMaximumBin() );
    double theRmsJmsPt1 = hDistRecoPt1_m[mBin_mpt]->GetRMS();
    hDistRecoPt1_m[mBin_mpt]->Fit("gaus", "","",theMeanJmsPt1 - 1*theRmsJmsPt1, theMeanJmsPt1 + .5*theRmsJmsPt1);
    if ( hDistRecoPt1_m[mBin_mpt]->GetFunction("gaus") == NULL ) {
      hJMS_pt1_m_1d->SetBinContent(mBin_mpt+1,-100);
      hJMS_pt1_m_1d->SetBinError(mBin_mpt+1,-100);
      hJMR_pt1_m_1d->SetBinContent(mBin_mpt+1,-100);
      hJMR_pt1_m_1d->SetBinError(mBin_mpt+1,-100);
    }
    else {
      float theJMS2 = hDistRecoPt1_m[mBin_mpt]->GetFunction("gaus")->GetParameter(1);
      float theJMSerr2 = hDistRecoPt1_m[mBin_mpt]->GetFunction("gaus")->GetParError(1);
      float theJMR2 = hDistRecoPt1_m[mBin_mpt]->GetFunction("gaus")->GetParameter(2);
      float theJMRerr2 = hDistRecoPt1_m[mBin_mpt]->GetFunction("gaus")->GetParError(2);
      hJMS_pt1_m_1d->SetBinContent(mBin_mpt+1,theJMS2);
      hJMS_pt1_m_1d->SetBinError(mBin_mpt+1,theJMSerr2);
      hJMR_pt1_m_1d->SetBinContent(mBin_mpt+1,theJMR2);
      hJMR_pt1_m_1d->SetBinError(mBin_mpt+1,theJMRerr2);
    }
    hDistRecoPt1_m[mBin_mpt]->Draw();
    c2->cd(mBin_mpt+1);
    double theMeanJmsPt2 = hDistRecoPt2_m[mBin_mpt]->GetBinCenter(hDistRecoPt2_m[mBin_mpt]->GetMaximumBin() );
    double theRmsJmsPt2 = hDistRecoPt2_m[mBin_mpt]->GetRMS();
    hDistRecoPt2_m[mBin_mpt]->Fit("gaus", "","",theMeanJmsPt2 - 1*theRmsJmsPt2, theMeanJmsPt2 + .5*theRmsJmsPt2);
    if ( hDistRecoPt2_m[mBin_mpt]->GetFunction("gaus") == NULL ) {
      hJMS_pt2_m_1d->SetBinContent(mBin_mpt+1,-100);
      hJMS_pt2_m_1d->SetBinError(mBin_mpt+1,-100);
      hJMR_pt2_m_1d->SetBinContent(mBin_mpt+1,-100);
      hJMR_pt2_m_1d->SetBinError(mBin_mpt+1,-100);
    }
    else {
      float theJMS2 = hDistRecoPt2_m[mBin_mpt]->GetFunction("gaus")->GetParameter(1);
      float theJMSerr2 = hDistRecoPt2_m[mBin_mpt]->GetFunction("gaus")->GetParError(1);
      float theJMR2 = hDistRecoPt2_m[mBin_mpt]->GetFunction("gaus")->GetParameter(2);
      float theJMRerr2 = hDistRecoPt2_m[mBin_mpt]->GetFunction("gaus")->GetParError(2);
      hJMS_pt2_m_1d->SetBinContent(mBin_mpt+1,theJMS2);
      hJMS_pt2_m_1d->SetBinError(mBin_mpt+1,theJMSerr2);
      hJMR_pt2_m_1d->SetBinContent(mBin_mpt+1,theJMR2);
      hJMR_pt2_m_1d->SetBinError(mBin_mpt+1,theJMRerr2);
    }
    hDistRecoPt2_m[mBin_mpt]->Draw();
    c3->cd(mBin_mpt+1);
    double theMeanJmsPt3 = hDistRecoPt3_m[mBin_mpt]->GetBinCenter( hDistRecoPt3_m[mBin_mpt]->GetMaximumBin() );
    double theRmsJmsPt3 = hDistRecoPt3_m[mBin_mpt]->GetRMS();
    hDistRecoPt3_m[mBin_mpt]->Fit("gaus", "","",theMeanJmsPt3 - 1*theRmsJmsPt3, theMeanJmsPt3 + .5*theRmsJmsPt3);
    if ( hDistRecoPt3_m[mBin_mpt]->GetFunction("gaus") == NULL ) {
      hJMS_pt3_m_1d->SetBinContent(mBin_mpt+1,-100);
      hJMS_pt3_m_1d->SetBinError(mBin_mpt+1,-100);
      hJMR_pt3_m_1d->SetBinContent(mBin_mpt+1,-100);
      hJMR_pt3_m_1d->SetBinError(mBin_mpt+1,-100);
    }
    else {
      float theJMS2 = hDistRecoPt3_m[mBin_mpt]->GetFunction("gaus")->GetParameter(1);
      float theJMSerr2 = hDistRecoPt3_m[mBin_mpt]->GetFunction("gaus")->GetParError(1);
      float theJMR2 = hDistRecoPt3_m[mBin_mpt]->GetFunction("gaus")->GetParameter(2);
      float theJMRerr2 = hDistRecoPt3_m[mBin_mpt]->GetFunction("gaus")->GetParError(2);
      hJMS_pt3_m_1d->SetBinContent(mBin_mpt+1,theJMS2);
      hJMS_pt3_m_1d->SetBinError(mBin_mpt+1,theJMSerr2);
      hJMR_pt3_m_1d->SetBinContent(mBin_mpt+1,theJMR2);
      hJMR_pt3_m_1d->SetBinError(mBin_mpt+1,theJMRerr2);
    }
    hDistRecoPt3_m[mBin_mpt]->Draw();
  }

  sprintf(fileName,"Plots/JES/JMSFitYongsun_pt1.png");
  c1->SaveAs(fileName);
  delete c1;
  sprintf(fileName,"Plots/JES/JMSFitYongsun_pt2.png");
  c2->SaveAs(fileName);
  delete c2;
  sprintf(fileName,"Plots/JES/JMSFitYongsun_pt3.png");
  c3->SaveAs(fileName);
  delete c3;
  
  TCanvas *c8b = new TCanvas("c8b","");
  handsomeTH1(hJMS_pt1_m_1d,1);
  handsomeTH1(hJMS_pt2_m_1d,2);
  handsomeTH1(hJMS_pt3_m_1d,4);
  hJMS_pt1_m_1d->GetXaxis()->SetTitle("Truth m/#it{p}_{T}");
  hJMS_pt1_m_1d->GetYaxis()->SetTitle("<(m/#it{p}_{T})^{2}_{Reco}>");
  hJMS_pt1_m_1d->Draw();
  hJMS_pt1_m_1d->GetXaxis()->SetRangeUser(0,0.24);
  hJMS_pt1_m_1d->GetYaxis()->SetRangeUser(0,0.05);
  hJMS_pt2_m_1d->Draw("same");
  hJMS_pt2_m_1d->SetMarkerStyle(25);
  hJMS_pt3_m_1d->Draw("same");
  hJMS_pt3_m_1d->SetMarkerStyle(28);
  drawCentrality(kSample,-1,0.70,0.86,1,24);
  ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
  TLegend *leg4b = new TLegend(0.1571854,0.6394175,0.6049639,0.8533743,NULL,"brNDC");      
  easyLeg(leg4b,"");
  leg4b->AddEntry(hJMS_pt1_m_1d, "100 < Truth Jet #it{p}_{T} < 125 GeV","plfe");
  leg4b->AddEntry(hJMS_pt2_m_1d, "158 < Truth Jet #it{p}_{T} < 199 GeV","plfe");
  leg4b->AddEntry(hJMS_pt3_m_1d, "398 < Truth Jet #it{p}_{T} < 500 GeV","plfe");
  leg4b->Draw();
  sprintf(fileName,"Plots/JES/JMS_ptn_m_sample%d.png",kSample);
  c8b->SaveAs(fileName);
  sprintf(fileName,"Plots/JES/JMS_ptn_m_sample%d.root",kSample);
  c8b->SaveAs(fileName);
  delete c8b;
  
  TCanvas *c8c = new TCanvas("c8c","");
  handsomeTH1(hJMR_pt1_m_1d,1);
  handsomeTH1(hJMR_pt2_m_1d,2);
  handsomeTH1(hJMR_pt3_m_1d,4);
  hJMR_pt1_m_1d->GetXaxis()->SetTitle("Truth m/#it{p}_{T}");
  hJMR_pt1_m_1d->GetYaxis()->SetTitle("#sigma[(m/#it{p}_{T})^{2}_{Reco}]");
  hJMR_pt1_m_1d->Draw();
  hJMR_pt1_m_1d->GetXaxis()->SetRangeUser(0,0.24);
  hJMR_pt1_m_1d->GetYaxis()->SetRangeUser(0,0.027);
  hJMR_pt2_m_1d->Draw("same");
  hJMR_pt2_m_1d->SetMarkerStyle(25);
  hJMR_pt3_m_1d->Draw("same");
  hJMR_pt3_m_1d->SetMarkerStyle(28);
  drawCentrality(kSample,-1,0.70,0.86,1,24);
  ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
  TLegend *leg4c = new TLegend(0.1571854,0.6394175,0.6049639,0.8533743,NULL,"brNDC");      
  easyLeg(leg4c,"");
  leg4c->AddEntry(hJMS_pt1_m_1d, "100 < Truth Jet #it{p}_{T} < 125 GeV","plfe");
  leg4c->AddEntry(hJMS_pt2_m_1d, "158 < Truth Jet #it{p}_{T} < 199 GeV","plfe");
  leg4c->AddEntry(hJMS_pt3_m_1d, "398 < Truth Jet #it{p}_{T} < 500 GeV","plfe");
  leg4c->Draw();
  sprintf(fileName,"Plots/JES/JMR_ptn_m_sample%d.png",kSample);
  c8c->SaveAs(fileName);
  sprintf(fileName,"Plots/JES/JMR_ptn_m_sample%d.root",kSample);
  c8c->SaveAs(fileName);
  delete c8c;
  
  for (int i = 0; i < 7; ++i ) {
    cout << " bin " << i << " JMS_pt1_m " << hJMS_pt1_m_1d->GetBinContent(i+1) << endl;
    cout << " bin " << i << " JMR_pt1_m " << hJMR_pt1_m_1d->GetBinContent(i+1) << endl;
  }      

  return 0;
}
