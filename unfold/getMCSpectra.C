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

  TH2D* hReweightingFactors; 
  TFile *reweightFile = new TFile(Form("ReweightingFactors/ReweightingFactors_kSample%d_DataOverMc.root",kSample));
  hReweightingFactors = (TH2D*) reweightFile->Get("ReweightingHisto_etaBin5_centBin0");
  
  int nPtBins;
  double ptBins[30];
  int nMBins;
  double mBins[40];
  getXbin(nPtBins, ptBins, 77);
  cout << " nXbins = " << nPtBins << endl;
  getYbin(nMBins, mBins, 772);
  cout << " nMBins = " << nMBins << endl;
  TH1D* hCrossSec_pt[6][7];
  TH1D* hCrossSec_pt_jz2[6][7];
  TH1D* hCrossSec_pt_jz2_weights[6][7];
  TH1D* hCrossSec_pt_jz2_stats[6][7];
  TH1D* hCrossSec_pt_jz3[6][7];
  TH1D* hCrossSec_pt_jz3_weights[6][7];
  TH1D* hCrossSec_pt_jz3_stats[6][7];
  TH1D* hCrossSec_m[6][7];
  TH2D* hCrossSec_m_pt[6][7];
  TH2D* hCrossSecTruth_m_pt[6][7];
  TH2D* hCrossSecTruth_secondHalfReweight_m_pt[6][7];
  TH2D* hCrossSec_m_pt_lowPt[6][7];
  TH2D* hCrossSec_m_pt_firstHalf[6][7];
  TH2D* hCrossSec_m_pt_secondHalfReweight[6][7];
  TH2D* hCrossSec_m_pt_lowPt_stats[6][7];
  int cent = 0;  
  float eta = 0.;
  long long entries;
  TString fname;
 
  for (int i=0; i < nEtaIntervals; ++i) { 
    for (int j=0; j < nCentIntervals; ++j) {
      char histName[256];
      sprintf(histName,"hCrossSecMc_pt_eta%dcent%d",i,j);
      hCrossSec_pt[i][j]  = new TH1D(histName,histName,nPtBins,ptBins);
      sprintf(histName,"hCrossSecMc_pt_jz2_eta%dcent%d",i,j);
      hCrossSec_pt_jz2[i][j]  = new TH1D(histName,histName,nPtBins,ptBins);
      sprintf(histName,"hCrossSecMc_pt_jz2Stats_eta%dcent%d",i,j);
      hCrossSec_pt_jz2_stats[i][j]  = new TH1D(histName,histName,nPtBins,ptBins);
      sprintf(histName,"hCrossSecMc_pt_jz3_eta%dcent%d",i,j);
      hCrossSec_pt_jz3[i][j]  = new TH1D(histName,histName,nPtBins,ptBins);
      sprintf(histName,"hCrossSecMc_pt_jz3Stats_eta%dcent%d",i,j);
      hCrossSec_pt_jz3_stats[i][j]  = new TH1D(histName,histName,nPtBins,ptBins);
      sprintf(histName,"hCrossSecMc_m_eta%dcent%d",i,j);
      hCrossSec_m[i][j]  = new TH1D(histName,histName,nMBins,mBins);
      sprintf(histName,"hCrossSecMc_m_pt_eta%dcent%d",i,j);
      hCrossSec_m_pt[i][j]  = new TH2D(histName,histName,nPtBins,ptBins,nMBins, mBins);
      sprintf(histName,"hCrossSecMcTruth_m_pt_eta%dcent%d",i,j);
      hCrossSecTruth_m_pt[i][j]  = new TH2D(histName,histName,nPtBins,ptBins,nMBins, mBins);
      sprintf(histName,"hCrossSecMcTruth_secondHalfReweight_m_pt_eta%dcent%d",i,j);
      hCrossSecTruth_secondHalfReweight_m_pt[i][j]  = new TH2D(histName,histName,nPtBins,ptBins,nMBins, mBins);
      sprintf(histName,"hCrossSecMc_m_pt_firstHalf_eta%dcent%d",i,j);
      hCrossSec_m_pt_firstHalf[i][j]  = new TH2D(histName,histName,nPtBins,ptBins,nMBins, mBins);
      sprintf(histName,"hCrossSecMc_m_pt_secondHalfReweight_eta%dcent%d",i,j);
      hCrossSec_m_pt_secondHalfReweight[i][j]  = new TH2D(histName,histName,nPtBins,ptBins,nMBins, mBins);
    }
  }
  
  getXbin(nPtBins, ptBins, 80);
  cout << " nXbins for 2d histos = " << nPtBins << endl;  
  for (int i=0; i < nEtaIntervals; ++i) { 
    for (int j=0; j < nCentIntervals; ++j) {
      char histName[256];

      sprintf(histName,"hCrossSecMc_m_pt_lowPt_eta%dcent%d",i,j);
      hCrossSec_m_pt_lowPt[i][j]  = new TH2D(histName,histName,nPtBins,ptBins,nMBins, mBins);
      sprintf(histName,"hCrossSecMc_m_pt_lowPt_stats_eta%dcent%d",i,j);
      hCrossSec_m_pt_lowPt_stats[i][j]  = new TH2D(histName,histName,nPtBins,ptBins,nMBins, mBins);
    }
  }

  TH1D* hFcalReweight;
  if ( kSample == kPbPb ) {
    TFile* fcal = new TFile("reweightFactors/FCal_HP_v_MB_weights.root");
    hFcalReweight = (TH1D*)fcal->Get("weight_MBov_to_HP");
  }
  
  cout << " Fill MC " << endl;

  tashka tr;

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
    for (Int_t i= 0; i<entries ; i++) {      
      //  way to cut processing off early ie if statFrac == 0.5, process half the file
      if ( i > entries * statFrac ) break;      
      tr.GetEntry(i);      
      
      double fcalWeight;
      if ( kSample==kPbPb) {
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

      double reweightFactor = 1.;
     
      //      cout << " genPt " << tr.jets_genPt << endl;
      if (tr.jets_genPt >= 100) {
	reweightFactor = hReweightingFactors->GetBinContent(hReweightingFactors->GetXaxis()->FindBin(tr.jets_genPt),hReweightingFactors->GetYaxis()->FindBin(tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)));
      }
      else {
	reweightFactor = hReweightingFactors->GetBinContent(hReweightingFactors->GetXaxis()->FindBin(100.),hReweightingFactors->GetYaxis()->FindBin(tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)));
	}
        
     
      //      cout << " reweightFactor " << reweightFactor << endl;
      //      if (reweightFactor < 0.2) reweightFactor = 0.2;
      //  if (reweightFactor > 2.0) reweightFactor = 2.0;

      //  if (ijz == 2 && tr.jets_ptCalib > 205.) continue; // 100 events
      //           if (ijz == 2 && tr.jets_ptCalib > 220.) continue; // 10 events
	     // if (ijz == 3 && tr.jets_ptCalib > 485.) continue; // 100 events
	//    if (ijz == 3 && tr.jets_ptCalib > 505.) continue; // 10 events
      // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32) -- looser pt cut 
      if ( ! pass2dEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, true) ) // isMC = true    
	continue;
      
      hCrossSec_m_pt_lowPt[etaBin][cent]->Fill(tr.jets_ptCalib,(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight);
      hCrossSec_m_pt_lowPt_stats[etaBin][cent]->Fill(tr.jets_ptCalib,(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)));
      if (kSample == kPP && eta < 2.1) {
	hCrossSec_m_pt_lowPt[5][0]->Fill(tr.jets_ptCalib, (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight);
	hCrossSec_m_pt_lowPt_stats[5][0]->Fill(tr.jets_ptCalib, (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)));
      }
      
      // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
      if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, true) ) // isMC = true    
	continue;

      hCrossSec_m_pt[etaBin][cent]->Fill(tr.jets_ptCalib,(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight);
      hCrossSecTruth_m_pt[etaBin][cent]->Fill(tr.jets_genPt,(tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight);
      if (i < entries/2) hCrossSec_m_pt_firstHalf[etaBin][cent]->Fill(tr.jets_ptCalib,(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight);
      else {
	hCrossSec_m_pt_secondHalfReweight[etaBin][cent]->Fill(tr.jets_ptCalib,(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight * reweightFactor);
	hCrossSecTruth_secondHalfReweight_m_pt[etaBin][cent]->Fill(tr.jets_genPt,(tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight * reweightFactor);
      }
      hCrossSec_pt[etaBin][cent]->Fill( tr.jets_ptCalib, tr.jets_weight * jzNorm * fcalWeight );
      if (ijz == 2) {
	hCrossSec_pt_jz2[etaBin][cent]->Fill( tr.jets_ptCalib, tr.jets_weight * jzNorm * fcalWeight );
	hCrossSec_pt_jz2_stats[etaBin][cent]->Fill( tr.jets_ptCalib);
      }
      if (ijz == 3) {
	hCrossSec_pt_jz3[etaBin][cent]->Fill( tr.jets_ptCalib, tr.jets_weight * jzNorm * fcalWeight );
	hCrossSec_pt_jz3_stats[etaBin][cent]->Fill( tr.jets_ptCalib);
      }
      hCrossSec_m[etaBin][cent]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight);
      if (kSample == kPP && eta < 2.1) {
	if (i < entries/2) hCrossSec_m_pt_firstHalf[5][0]->Fill(tr.jets_ptCalib,(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight);
	else {
	  hCrossSec_m_pt_secondHalfReweight[5][0]->Fill(tr.jets_ptCalib,(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight * reweightFactor);
	  hCrossSecTruth_secondHalfReweight_m_pt[5][0]->Fill(tr.jets_genPt,(tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight * reweightFactor);
	}
	hCrossSec_m_pt[5][0]->Fill(tr.jets_ptCalib,(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight);
	hCrossSecTruth_m_pt[5][0]->Fill(tr.jets_genPt,(tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight);
	hCrossSec_pt[5][0]->Fill( tr.jets_ptCalib, tr.jets_weight * jzNorm * fcalWeight);
	if (ijz == 2) {
	  hCrossSec_pt_jz2[5][0]->Fill( tr.jets_ptCalib, tr.jets_weight * jzNorm * fcalWeight);
	  hCrossSec_pt_jz2_stats[5][0]->Fill( tr.jets_ptCalib);
	}
	if (ijz == 3) {
	  hCrossSec_pt_jz3[5][0]->Fill( tr.jets_ptCalib, tr.jets_weight * jzNorm * fcalWeight);
	  hCrossSec_pt_jz3_stats[5][0]->Fill( tr.jets_ptCalib);
	}
	hCrossSec_m[5][0]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight);
      }
      
    }
  }
  cout << " writing to file " << endl; 
  // write results to file
  //    TFile *fout = new TFile(Form("rawSpectra/RawSpectraMC_sample%d_ptCut10Events.root",kSample),"RECREATE");
  // TFile *fout = new TFile(Form("rawSpectra/RawSpectraMC_sample%d_ptCut100Events.root",kSample),"RECREATE"); 
  TFile *fout = new TFile(Form("rawSpectra/RawSpectraMC_sample%d.root",kSample),"RECREATE"); 
  for (int etaBin = 0; etaBin < nEtaIntervals; ++etaBin) {
    for (int centBin = 0; centBin < nCentIntervals; ++centBin) {
      hCrossSec_pt_jz2_weights[etaBin][centBin] = (TH1D*) hCrossSec_pt_jz2[etaBin][centBin]->Clone(Form("jz2Weights_centBin%d_etaBin%d",centBin,etaBin));
	hCrossSec_pt_jz2_weights[etaBin][centBin]->Divide(hCrossSec_pt_jz2_stats[etaBin][centBin]);
	hCrossSec_pt_jz3_weights[etaBin][centBin] = (TH1D*) hCrossSec_pt_jz3[etaBin][centBin]->Clone(Form("jz3Weights_centBin%d_etaBin%d",centBin,etaBin));
	hCrossSec_pt_jz3_weights[etaBin][centBin]->Divide(hCrossSec_pt_jz3_stats[etaBin][centBin]);
	hCrossSec_pt_jz2_weights[etaBin][centBin]->Write();
	hCrossSec_pt_jz3_weights[etaBin][centBin]->Write();
	hCrossSec_m[etaBin][centBin]->Write();
	hCrossSec_pt[etaBin][centBin]->Write();
	hCrossSec_pt_jz2[etaBin][centBin]->Write();
	hCrossSec_pt_jz2_stats[etaBin][centBin]->Write();
	hCrossSec_pt_jz3[etaBin][centBin]->Write();
	hCrossSec_pt_jz3_stats[etaBin][centBin]->Write();
	hCrossSec_m_pt[etaBin][centBin]->Write();
	hCrossSecTruth_m_pt[etaBin][centBin]->Write();
	hCrossSecTruth_secondHalfReweight_m_pt[etaBin][centBin]->Write();
	hCrossSec_m_pt_lowPt[etaBin][centBin]->Write();
	hCrossSec_m_pt_firstHalf[etaBin][centBin]->Write();
	hCrossSec_m_pt_secondHalfReweight[etaBin][centBin]->Write();
	hCrossSec_m_pt_lowPt_stats[etaBin][centBin]->Write();
    }
  }

  fout->Close();  
  return 0;
}
