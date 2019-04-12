 // Goal: take MC from TTree to Spectra applying MC reweighting (pt, mass and pt vs mass)
// KEEP THESE SCRIPTS SMALL AND SINGLE PURPOSED!

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include </home/mike/RootUtils/AtlasStyle.C>
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
  SetAtlasStyle();
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  gStyle->SetOptTitle(0);
  int nCentIntervals;
  if (kSample == kPP) nCentIntervals = 1;
  else nCentIntervals = 7; 
  int nEtaIntervals = 6;  

  double eventWeight;
  int nPtBins;
  double ptBins[30];
  int nMBins;
  double mBins[40];
  getXbin(nPtBins, ptBins, 77);
  cout << " nXbins = " << nPtBins << endl;
  getYbin(nMBins, mBins, 772);
  cout << " nMBins = " << nMBins << endl;
  TH1D* hCrossSec_pt[6][7];
  TH1D* hCrossSec_m[6][7];
  TH2D* hCrossSec_m_pt[6][7];
  TH2D* hCrossSec_m_pt_lowPt[6][7];
  TH2D* hCrossSec_m_pt_lowPt_stats[6][7];
  int cent = 0;  
  float eta = 0.;
  long long entries;
  TString fname;

 
  for (int i=0; i < nEtaIntervals; ++i) { 
    for (int j=0; j < nCentIntervals; ++j) {
      char histName[256];
      sprintf(histName,"hCrossSecData_pt_eta%dcent%d",i,j);
      hCrossSec_pt[i][j]  = new TH1D(histName,histName,nPtBins,ptBins);
      sprintf(histName,"hCrossSecData_m_eta%dcent%d",i,j);
      hCrossSec_m[i][j]  = new TH1D(histName,histName,nMBins,mBins);
      sprintf(histName,"hCrossSecData_m_pt_eta%dcent%d",i,j);
      hCrossSec_m_pt[i][j]  = new TH2D(histName,histName,nPtBins,ptBins,nMBins, mBins);
    }
  }
  getXbin(nPtBins, ptBins, 80);
  cout << " nXbins in 2d histo = " << nPtBins << endl;
  for (int i=0; i < nEtaIntervals; ++i) { 
    for (int j=0; j < nCentIntervals; ++j) {
      char histName[256];
      sprintf(histName,"hCrossSecData_m_pt_eta%dcent%d_lowPt",i,j);
      hCrossSec_m_pt_lowPt[i][j]  = new TH2D(histName,histName,nPtBins,ptBins,nMBins, mBins);
      sprintf(histName,"hCrossSecData_m_pt_lowPt_stats_eta%dcent%d",i,j);
      hCrossSec_m_pt_lowPt_stats[i][j]  = new TH2D(histName,histName,nPtBins,ptBins,nMBins, mBins);
    }
  }
  
  cout << " Fill Data " << endl;
  if ( kSample == kPbPb ) {
    fname = pbpbDataString; 
  }
  else if ( kSample == kPP) {
    fname = ppDataString;
  }
  TFile* fData = new TFile(Form("../ntuples/%s",fname.Data()));
  TTree* tree = (TTree*)fData->Get("tr");
  cout << " Setting Tashka for data ..." << endl;
  tashka tr;
  tr.Init(tree);
  entries = tree->GetEntries();

  for (Int_t i= 0; i<entries ; i++) {
    tr.GetEntry(i);
    //  way to cut processing off early ie if statFrac == 0.5, process half the file
    if ( i > entries * statFrac) break;
    // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
    cent = tr.jets_cent;
    eta  = TMath::Abs(tr.jets_yCalib);

    int etaBin = -1;
    if (eta < 0.3) etaBin = 0;
    else if (eta >= 0.3 && eta < 0.8) etaBin = 1;
    else if (eta >= 0.8 && eta < 1.2) etaBin = 2;
    else if (eta >= 1.2 && eta < 1.6) etaBin = 3;
    else if (eta >= 1.6 && eta < 2.1) etaBin = 4;
    else continue;
    //      cout << " passed " << endl;
    if (kSample == kPP) {
      // 1 / pp lumi (nb)
      eventWeight = ppEvtWgt;
    }
    else if (kSample == kPbPb && cent == 0) {
      // 1/ ( TAA_cent * nEvents_cent) (nb)
      eventWeight = hi9EvtWgt_0;
    }
    else if (kSample == kPbPb && cent == 1) {
      eventWeight = hi9EvtWgt_1;
    }
    else if (kSample == kPbPb && cent == 2) {
      eventWeight = hi9EvtWgt_2;
    }
    else if (kSample == kPbPb && cent == 3) {
      eventWeight = hi9EvtWgt_3;
    }
    else if (kSample == kPbPb && cent == 4) {
      eventWeight = hi9EvtWgt_4;
    }
    else if (kSample == kPbPb && cent == 5) {
      eventWeight = hi9EvtWgt_5;
    }
    else if (kSample == kPbPb && cent == 6) {
      eventWeight = hi9EvtWgt_6;
    }

    if ( ! pass2dEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, true) ) // isMC = true    
	continue;
    hCrossSec_m_pt_lowPt[etaBin][cent]->Fill(tr.jets_ptCalib,(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), eventWeight);
    hCrossSec_m_pt_lowPt_stats[etaBin][cent]->Fill(tr.jets_ptCalib,(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)));
    if (kSample == kPP && eta < 2.1) {
      hCrossSec_m_pt_lowPt[5][0]->Fill(tr.jets_ptCalib, (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), eventWeight);
      hCrossSec_m_pt_lowPt_stats[5][0]->Fill(tr.jets_ptCalib, (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)));
    }

    if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, false) ) // isMC = false
      continue;
    hCrossSec_m_pt[etaBin][cent]->Fill(tr.jets_ptCalib,(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), eventWeight);
    hCrossSec_pt[etaBin][cent]->Fill( tr.jets_ptCalib, eventWeight );
    hCrossSec_m[etaBin][cent]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), eventWeight);
    if (kSample == kPP && eta < 2.1) {
      hCrossSec_pt[5][0]->Fill( tr.jets_ptCalib, eventWeight);
      hCrossSec_m[5][0]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), eventWeight);
      hCrossSec_m_pt[5][0]->Fill(tr.jets_ptCalib, (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), eventWeight);
    }
  }

  cout << " writing to file " << endl; 
  // write results to file
  TFile *fout = new TFile(Form("rawSpectra/RawSpectraData_sample%d.root",kSample),"RECREATE"); 
  for (int etaBin = 0; etaBin < nEtaIntervals; ++etaBin) {
    for (int centBin = 0; centBin < nCentIntervals; ++centBin) {  
	hCrossSec_m[etaBin][centBin]->Write();
	hCrossSec_pt[etaBin][centBin]->Write();
	hCrossSec_m_pt[etaBin][centBin]->Write();
	hCrossSec_m_pt_lowPt[etaBin][centBin]->Write();
	hCrossSec_m_pt_lowPt_stats[etaBin][centBin]->Write();
    }
  }

  fout->Close();  
  return 0;
}
