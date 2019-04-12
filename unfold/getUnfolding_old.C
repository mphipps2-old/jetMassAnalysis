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
  int nUnfIterations = 2;
  
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

  //  double binWidth[60];
  // located in ../unfoldingUtil.h
  // Xbin -> pt
  getXbin(nPtBins, ptBins, 77);
  cout << " nPtBins = " << nPtBins << endl;
  // Ybin -> m2/pt2
  getYbin(nMBins, mBins, 772);

  for (int i = 0; i < nPtBins; ++i) {
    cout << " ptBin " << i << " = " <<  ptBins[i] <<endl;
    binWidth[i] = ptBins[i+1]-ptBins[i];
  }
  for (int i = 0; i < nMBins; ++i) {
    cout << " mBin " << i << " = " <<  mBins[i] <<endl;
    //    binWidth[i] = mBin[i+1]-mBin[i];
  }


  double eventWeight = 1.;
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
  TH2D* hReweightingFactors; 
  TFile *reweightFile = new TFile("ReweightingFactors/ReweightingFactors_kSample0_etaBin5.root");
  hReweightingFactors = (TH2D*) reweightFile->Get("ReweightingHisto");
  hReweightingFactors->Smooth();
  
  TH2D *hDataReco[6][7];  
  TFile* fData = new TFile(Form("../ntuples/%s",fname.Data()));
  TTree* tree = (TTree*)fData->Get("tr");
  cout << " Setting Tashka for data ..." << endl;
  tashka tr;
  tr.Init(tree);
  entries = tree->GetEntries();

  RooUnfoldResponse* rooUnfoldRes[6][7]; // eta, cent
  TH2D *hMcReco[6][7];
  TH2D *hMcTruth[6][7];
  TH2D *hResMatrix[6][7];
  TH2D* hDataUnf[6][7]; // unfolding iter
  for (int eta = 0; eta < nEtaBins; ++eta) {
    for (int cent = 0; cent < nCentBins; ++cent) {
	
      hMcReco[eta][cent] = new TH2D(Form("hMcRecoPtVsM2Pt2_eta%d_cent%d_firstHalf",eta,cent), ";Reco #it{p}_{T} (GeV);Reco (m/p_{T})^{2}",nPtBins,ptBins,nMBins,mBins);
      hMcTruth[eta][cent] = new TH2D(Form("hTruthPtVsM2Pt2_eta%d_cent%d_firstHalf",eta,cent), ";#it{p}_{T} (GeV);m^{2}/#it{p}_{T}^{2}",nPtBins,ptBins,nMBins,mBins);
      hDataReco[eta][cent] = new TH2D(Form("hDataRecoPtVsM2Pt2_eta%d_cent%d",eta,cent), ";Reco #it{p}_{T} (GeV);Reco (m/p_{T})^{2}",nPtBins,ptBins,nMBins,mBins);
      rooUnfoldRes[eta][cent] = new RooUnfoldResponse( hMcReco[eta][cent], hMcTruth[eta][cent] );
      rooUnfoldRes[eta][cent]->SetName(Form("responseMatrix_eta%d_cent%d",eta,cent));
      hResMatrix[eta][cent] = new TH2D(Form("h2dResponseMatrix_eta%d_cent%d",eta,cent), ";Bin # of reco ({p}_{T}, (m/#it{p}_{T})^{2});Bin # of truth ({p}_{T}, (m/#it{p}_{T})^{2})",nPtBins,ptBins,nMBins,mBins);
    }
  }
  
  for (Int_t i= 0; i<entries ; i++) {
    tr.GetEntry(i);
    //  way to cut processing off early ie if statFrac == 0.5, process half the file
    if ( i > entries * statFrac) break;
    // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
    if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, false) ) // isMC = false
      continue;
    if (kSample == kPP) {
      // 1 / pp lumi (nb)
      eventWeight = ppEvtWgt;
    }
    else if (kSample == kPbPb && tr.jets_cent == 0) {
      eventWeight = hi9EvtWgt_0;
    }
    else if (kSample == kPbPb && tr.jets_cent == 1) {
      eventWeight = hi9EvtWgt_1;
    }
    else if (kSample == kPbPb && tr.jets_cent == 2) {
      eventWeight = hi9EvtWgt_2;
    }
    else if (kSample == kPbPb && tr.jets_cent == 3) {
      eventWeight = hi9EvtWgt_3;
    }
    else if (kSample == kPbPb && tr.jets_cent == 4) {
      eventWeight = hi9EvtWgt_4;
    }
    else if (kSample == kPbPb && tr.jets_cent == 5) {
      eventWeight = hi9EvtWgt_5;
    }
    else if (kSample == kPbPb && tr.jets_cent == 6) {
      // 1/ ( TAA_cent * nEvents_cent) (nb)
      eventWeight = hi9EvtWgt_6;
    }
    int cent = tr.jets_cent;
    float eta  = TMath::Abs(tr.jets_yCalib);
    int etaBin = -1;
    if (eta < 0.3) etaBin = 0;
    else if (eta >= 0.3 && eta < 0.8) etaBin = 1;
    else if (eta >= 0.8 && eta < 1.2) etaBin = 2;
    else if (eta >= 1.2 && eta < 1.6) etaBin = 3;
    else if (eta >= 1.6 && eta < 2.1) etaBin = 4;   
    else continue;
    //    cout << " scalingn by eventweight: " << eventWeight << endl;
    if (eta < 2.1 && kSample == kPP) {
      hDataReco[5][cent]->Fill( tr.jets_ptCalib, tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib), eventWeight); 
    }
    hDataReco[etaBin][cent]->Fill( tr.jets_ptCalib, tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib), eventWeight);
  }  

  
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
      double reweightFactor = hReweightingFactors->GetBinContent(hReweightingFactors->GetXaxis()->FindBin(tr.jets_ptCalib),hReweightingFactors->GetYaxis()->FindBin(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)));
      //    double reweightFactor = 1.;
      
           if (tr.jets_genPt >= 90) reweightFactor = hReweightingFactors->GetBinContent(hReweightingFactors->GetXaxis()->FindBin(tr.jets_genPt),hReweightingFactors->GetYaxis()->FindBin(tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)));
      else {
	reweightFactor = hReweightingFactors->GetBinContent(hReweightingFactors->GetXaxis()->FindBin(90.),hReweightingFactors->GetYaxis()->FindBin(tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)));
	//	cout << " truthPt " << tr.jets_genPt << " recoPt " << tr.jets_ptCalib << " truthM2 " << tr.jets_genMass2 << 
	}
            
	   //      if (reweightFactor == 0 || reweightFactor != reweightFactor || reweightFactor > 1000) continue;
      if (reweightFactor < 0.3) reweightFactor = 0.3;
      if (reweightFactor > 1.5) reweightFactor = 1.5;
      // data/MC
      reweightFactor = 1. / reweightFactor;

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
	rooUnfoldRes[5][cent]->Fill(tr.jets_ptCalib, (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_genPt, (tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight * reweightFactor);	  
      }
      rooUnfoldRes[etaBin][cent]->Fill(tr.jets_ptCalib, (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_genPt, (tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight * reweightFactor);      
    }
  }
  /*
  cout << " get response matrices " << endl;  
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      hResMatrix[etaBin][centBin] = (TH2D*) rooUnfoldRes[etaBin][centBin]->Hresponse();
    }
  }
  */
  cout << " unfold " << endl;
  
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for ( int centBin=0 ; centBin < nCentBins; ++centBin) {
      if ( (kSample == kPP) && ( centBin != 0 ) )      continue;
      // Data
      RooUnfoldBayes unfoldData (rooUnfoldRes[etaBin][centBin], hDataReco[etaBin][centBin], nUnfIterations);    
      hDataUnf[etaBin][centBin]= (TH2D*)unfoldData.Hreco();
      hDataUnf[etaBin][centBin]->SetName( Form("hDataUnf_eta%d_cent%d",etaBin,centBin));
    }
  }

   
  cout << " done unfolding " << endl;
  TH1D *hDataUnf_Pt[6][7][30];
  TH1D *hDataUnf_Pt_Integrated[6][7];
  TH1D *hDataUnf_M[6][7][20];
  TH1D *hDataUnf_M_Integrated[6][7];
 
  // write results to file
  TFile *fout = new TFile(Form("unfSpectra/UnfoldingSpectra_sample%d_withWeights.root",kSample),"RECREATE");
 
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      for (int mBin = 0; mBin < nMBins; ++mBin) {
	hDataUnf_Pt[etaBin][centBin][mBin] = (TH1D*)hDataUnf[etaBin][centBin]->ProjectionX(Form("hdataUnf_Pt_eta%d_icent%d_m%d",etaBin,centBin,mBin),mBin+2,mBin+2);	 
	TH1ScaleByWidth(hDataUnf_Pt[etaBin][centBin][mBin]);
	TH1ScaleByRapidity(hDataUnf_Pt[etaBin][centBin][mBin], etaBin);
	hDataUnf_Pt[etaBin][centBin][mBin]->Write();
      }
    }
  }


  
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) { 
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      for (int ptBin = 0; ptBin < nPtBins; ++ptBin) {
	hDataUnf_M[etaBin][centBin][ptBin] = (TH1D*)hDataUnf[etaBin][centBin]->ProjectionY(Form("hdataUnf_M_eta%d_icent%d_pt%d",etaBin,centBin,ptBin),ptBin+2,ptBin+2);
	TH1ScaleByWidth(hDataUnf_M[etaBin][centBin][ptBin]);
	TH1ScaleByRapidity(hDataUnf_M[etaBin][centBin][ptBin],etaBin);
    	hDataUnf_M[etaBin][centBin][ptBin]->Write();
      }
    }
  }

  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {     
      hDataUnf_Pt_Integrated[etaBin][centBin] = (TH1D*)hDataUnf[etaBin][centBin]->ProjectionX(Form("hDataUnfPt_eta%d_cent%d",etaBin,centBin),2,nMBins+1);
      hDataUnf_M_Integrated[etaBin][centBin] = (TH1D*)hDataUnf[etaBin][centBin]->ProjectionY(Form("hDataUnfM_eta%d_cent%d",etaBin,centBin),2,nPtBins+1);
      TH1ScaleByWidth(hDataUnf_M_Integrated[etaBin][centBin]);
      TH1ScaleByRapidity(hDataUnf_M_Integrated[etaBin][centBin],etaBin);
      TH1ScaleByWidth(hDataUnf_Pt_Integrated[etaBin][centBin]);
      TH1ScaleByRapidity(hDataUnf_Pt_Integrated[etaBin][centBin],etaBin);
      hDataUnf_Pt_Integrated[etaBin][centBin]->Write();
      hDataUnf_M_Integrated[etaBin][centBin]->Write();
    }
  }
  
  fout->Close();
  
  return 0;
}
