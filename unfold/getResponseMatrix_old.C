// Goal: produce ratio plot of Data/MC for dsigma/dpt vs pt. Eta binning integrated from -2.1 < eta < 2.1. Works for pp or PbPb by flipping kSample switch
// KEEP THESE SCRIPTS SMALL AND SINGLE PURPOSED!

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include </home/mike/RootUtils/AtlasStyle.C>
//#include <AtlasStyle.h>
#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
#include <vector>
#include <string>
#include "TH1D.h"

#include "/home/mike/Desktop/JetAnalysis/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "/home/mike/Desktop/JetAnalysis/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
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
  gStyle->SetPalette(kBird);
  int nCentBins;
  if (kSample == kPP) nCentBins = 1;
  else nCentBins = 7; 
  int nEtaBins = 6;
  // Number of iterations -- note this number already optimized in previous script
  
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
  int nMBins;
  double mBins[40];

  //  double binWidth[60];
  // located in ../unfoldingUtil.h
  // Xbin -> pt
  getXbin(nPtBins, ptBins, 80);
  cout << " nPtBins = " << nPtBins << endl;
  // Ybin -> m2/pt2
  getYbin(nMBins, mBins, 772);

  long long entries;
  TString fname;

  TH1D* hFcalReweight;
  if ( kSample == kPbPb ) {
    TFile* fcal = new TFile("reweightFactors/FCal_HP_v_MB_weights.root");
    hFcalReweight = (TH1D*)fcal->Get("weight_MBov_to_HP");
  }
  TH2D* hReweightingFactors; 
  TFile *reweightFile = new TFile(Form("ReweightingFactors/ReweightingFactors_kSample%d_DataOverMc.root",kSample));
  hReweightingFactors = (TH2D*) reweightFile->Get("ReweightingHisto_etaBin5_centBin0");
  //  hReweightingFactors->Smooth();
  

  RooUnfoldResponse* rooUnfoldRes[6][7]; // eta, cent
  TH2D *hMcReco[6][7];
  TH2D *hMcTruth[6][7];
  TH2D *hResponse_pt[6][7];
  TH2D *hResponse_m[6][7];
  TH2D *hResMatrix[6][7];

  std::vector<std::vector<std::vector<std::vector<double> > > > response4d;
  
  for (int eta = 0; eta < nEtaBins; ++eta) {
    for (int cent = 0; cent < nCentBins; ++cent) {
      hMcReco[eta][cent] = new TH2D(Form("hMcReco_pt_eta%d_cent%d",eta,cent), ";Reco #it{p}_{T} (GeV);Reco m^{2}/#it{p}_{T}^{2}",nPtBins,ptBins,nMBins,mBins);
      hMcTruth[eta][cent] = new TH2D(Form("hMcTruth_m_eta%d_cent%d",eta,cent), ";Reco #it{p}_{T};Reco m^{2}/#it{p}_{T}^{2}",nPtBins,ptBins,nMBins,mBins);	
      hResponse_pt[eta][cent] = new TH2D(Form("hResponse_pt_eta%d_cent%d",eta,cent), ";Reco #it{p}_{T} (GeV);Truth #it{p}_{T} (GeV)",nPtBins,ptBins,nPtBins,ptBins);
      hResponse_m[eta][cent] = new TH2D(Form("hResponse_m_eta%d_cent%d",eta,cent), ";Reco m^{2}/#it{p}_{T}^{2};Truth m^{2}/#it{p}_{T}^{2}",nMBins,mBins,nMBins,mBins);
      rooUnfoldRes[eta][cent] = new RooUnfoldResponse( hMcReco[eta][cent], hMcTruth[eta][cent] );
      rooUnfoldRes[eta][cent]->SetName(Form("responseMatrix_eta%d_cent%d",eta,cent));
    }
  }
  
  tashka tr;
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
    for (Int_t i= 0; i<entries ; i++) {
      //  way to cut processing off early ie if statFrac == 0.5, process half the file
      if ( i > entries * statFrac ) break;      
      tr.GetEntry(i);
      // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
      if ( ! passJesEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib) ) // isMC = true
	continue;

      double fcalWeight;
      if ( kSample==kPbPb) {
	fcalWeight = hFcalReweight->GetBinContent(hFcalReweight->GetXaxis()->FindBin(tr.jets_fcal));
      }

      else fcalWeight = 1.0;
      // double reweightFactor = hReweightingFactors->GetBinContent(hReweightingFactors->GetXaxis()->FindBin(tr.jets_ptCalib),hReweightingFactors->GetYaxis()->FindBin(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)));
      double reweightFactor = 1.;
      //      cout << " genPt " << tr.jets_genPt << endl;
      if (tr.jets_genPt >= 100) {
	reweightFactor = hReweightingFactors->GetBinContent(hReweightingFactors->GetXaxis()->FindBin(tr.jets_genPt),hReweightingFactors->GetYaxis()->FindBin(tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)));
      }
      else {
	reweightFactor = hReweightingFactors->GetBinContent(hReweightingFactors->GetXaxis()->FindBin(100.),hReweightingFactors->GetYaxis()->FindBin(tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)));
	}
        
      //      if (reweightFactor == 0 || reweightFactor != reweightFactor || reweightFactor > 1000) continue;
      //      cout << " reweightFactor " << reweightFactor << endl;
      if (reweightFactor < 0.4) reweightFactor = 0.4;
      if (reweightFactor > 1.5) reweightFactor = 1.5;
      // data/MC
      //      reweightFactor = 1. / reweightFactor;
      
      int cent = tr.jets_cent;
      float eta  = TMath::Abs(tr.jets_yCalib);
      int etaBin = -1;
      if (eta < 0.3) etaBin = 0;
      else if (eta >= 0.3 && eta < 0.8) etaBin = 1;
      else if (eta >= 0.8 && eta < 1.2) etaBin = 2;
      else if (eta >= 1.2 && eta < 1.6) etaBin = 3;
      else if (eta >= 1.6 && eta < 2.1) etaBin = 4;   
      else continue;
      // note: now we're using the full dataset for the unfolding
      if (eta < 2.1 && kSample == kPP) {
	hResponse_pt[5][cent]->Fill(tr.jets_ptCalib, tr.jets_genPt, tr.jets_weight * jzNorm * fcalWeight *reweightFactor);	
	hResponse_m[5][cent]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), (tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight * reweightFactor);
	rooUnfoldRes[5][cent]->Fill(tr.jets_ptCalib, (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_genPt, (tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight * reweightFactor);	  
      }
      hResponse_pt[etaBin][cent]->Fill(tr.jets_ptCalib, tr.jets_genPt, tr.jets_weight * jzNorm * fcalWeight *reweightFactor);	
      hResponse_m[etaBin][cent]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), (tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight * reweightFactor);
      rooUnfoldRes[etaBin][cent]->Fill(tr.jets_ptCalib, (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_genPt, (tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight * reweightFactor);      
    }
  }
  
  // write results to file
  TFile *fout = new TFile(Form("unfSpectra/rooUnfoldResponse_ReweightingMin0P4Max1P5_sample%d.root",kSample),"RECREATE");
 
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      hResMatrix[etaBin][centBin] = (TH2D*) rooUnfoldRes[etaBin][centBin]->Hresponse();
      hResponse_pt[etaBin][centBin]->Write();
      hResponse_m[etaBin][centBin]->Write();
      hResMatrix[etaBin][centBin]->Write();
      rooUnfoldRes[etaBin][centBin]->Write(); 
    }
  }
  
  fout->Close();
  
  return 0;
}
