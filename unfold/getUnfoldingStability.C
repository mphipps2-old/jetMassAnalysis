// Goal: produce ratio plot of Data/MC for dsigma/dpt vs pt. Eta binning integrated from -2.1 < eta < 2.1. Works for pp or PbPb by flipping kSample switch
// KEEP THESE SCRIPTS SMALL AND SINGLE PURPOSED!

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
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
  
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kRainBow);
  int nCentBins;
  if (kSample == kPP) nCentBins = 1;
  else nCentBins = 7; 
  int nEtaBins = 6;  
  int nUnfIterations = 3;
  
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
    if (mBin < 6)    hDist_m[mBin] = new TH1D(name,name, 100,0,9);
    else hDist_m[mBin] = new TH1D(name,name, 100,0,2);
    sprintf(name, "hDistm_m%d",mBin);
    hDistJes_m[mBin] = new TH1D(name,name, 50,0.7,1.3);

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
  
  RooUnfoldResponse* hRooUnfold[6][7]; // eta, cent
  TH2D *hMcReco[6][7];
  TH2D *hTruth[6][7];
  TH2D *hDataReco[6][7];
  TH2D *hResMatrix[6][7];
  TH2D* hMcUnf[7][20]; // unfolding iter
  TH2D* hDataUnf[7][20]; // unfolding iter
  for (int eta = 0; eta < nEtaBins; ++eta) {
    for (int cent = 0; cent < nCentBins; ++cent) {
      hMcReco[eta][cent] = new TH2D(Form("hMcRecoPtVsM2Pt2_eta%d_cent%d",eta,cent), ";Reco #it{p}_{T} (GeV);Reco (m/p_{T})^{2}",nPtBins,ptBin,nMBins,mBin;)
      hTruth[eta][cent] = new TH2D(Form("hTruthPtVsM2Pt2_eta%d_cent%d",eta,cent), ";Truth #it{p}_{T} (GeV);Truth (m/#it{p}_{T})^{2}",nPtBins,ptBin,nMBins,mBin);
      hDataReco[eta][cent] = new TH2D(Form("hDataRecoPtVsM2Pt2_eta%d_cent%d",eta,cent), ";Reco #it{p}_{T} (GeV);Reco (m/p_{T})^{2}",nPtBins,ptBin,nMBins,mBin;)
      rooUnfoldRes[eta][cent] = new RooUnfoldResponse( hReco[eta][cent], hTruth[eta][cent] );
      rooUnfoldRes[eta][cent]->SetName(Form("responseMatrix_eta%d_cent%d",eta,cent));
      hResMatrix[eta][cent] = new TH2D(Form("h2dResponseMatrix_eta%d_cent%d",eta,cent), ";Bin # of reco ({p}_{T}, (m/#it{p}_{T})^{2});Bin # of truth ({p}_{T}, (m/#it{p}_{T})^{2})",nPtBins,ptBin,nMBins,mBin);
    }
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
      if (ijz == 4 && jetCount == 6209 ) cout << " jetCount " << jetCount << " pt " << tr.jets_ptCalib << " mass2Calib " << tr.jets_mass2Calib << " eta " << tr.jets_yCalib <<  endl;
      if ( i > entries * statFrac ) break;      
      tr.GetEntry(i);
      // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
      if ( ! passJesEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib) ) // isMC = true
	continue;
      if (ijz == 4 && jetCount == 6209 ) cout << " passed Jes event " << endl;     
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

      if (eta < 2.1 && kSample == kPP) {
xo	hMcReco[5][cent]->Fill(tr.jets_ptCalib, tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib), tr.jets_weight * jzNorm * fcalWeight);
	hMcTruth[5][cent]->Fill(tr.jets_genCalib, (tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt)), tr.jets_weight * jzNorm * fcalWeight);
      }
      hMcReco[eta][cent]->Fill(tr.jets_ptCalib, (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight);
      hMcTruth[eta][cent]->Fill(tr.jets_genCalib, tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt), tr.jets_weight * jzNorm * fcalWeight);      
    }
  }


  for (Int_t i= 0; i<entries ; i++) {
    tr.GetEntry(i);
    //  way to cut processing off early ie if statFrac == 0.5, process half the file
    if ( i > entries * statFrac) break;
    // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
    //    if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_cent, tr.jets_etaCalib, tr.jets_phiCalib, false) ) // isMC = false
    // NOTE: CHANGE BACK WITH NEW TTREE!!!
    if ( ! passEvent(tr.jets_cent, tr.jets_genMass2, tr.jets_ptCalib, tr.jets_cent, tr.jets_etaCalib, tr.jets_phiCalib, false) ) // isMC = false
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
    cent = tr.jets_cent;
    eta  = TMath::Abs(tr.jets_yCalib);
    int etaBin = -1;
    if (eta < 0.3) etaBin = 0;
    else if (eta >= 0.3 && eta < 0.8) etaBin = 1;
    else if (eta >= 0.8 && eta < 1.2) etaBin = 2;
    else if (eta >= 1.2 && eta < 1.6) etaBin = 3;
    else if (eta >= 1.6 && eta < 2.1) etaBin = 4;   
    else continue;
    if (eta < 2.1 && kSample == kPP) {
      hDataReco[5][cent]->Fill( tr.jets_ptCalib, tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib), eventWeight);
    }
    hDataReco[eta][cent]->Fill( tr.jets_ptCalib, tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib), eventWeight);
  }


  
  
  for (int etaBin = 0; etaBin < nEtaBins; +etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      hResMatrix[etaBin][centBin] = (TH2D*) rooUnfoldRes[etaBin][centBin]->Hresponse();
    }
  }

  for (int etaBin = 0; etaBin < nEtaBins; +etaBin) {
    for ( int centBin=0 ; centBin < nCentBins; ++centBin) {
    if ( (kSample == kPP) && ( centBin != 0 ) )      continue;
    // loop through unfolding iterations
    for ( int it = 0 ; it < nUnfIterations; it++) {  
      RooUnfoldBayes unfoldMc (rooUnfoldRes[etaBin][centBin], hMcReco[etaBin][centBin], it);    
      hMcUnf[etaBin][centBin][it] = (TH2D*)unfoldMc.Hreco();
      hMcUnf[etaBin][centBin][it]->SetName( Form("hMcUnf_eta%d_cent%d_iter%d",etaBin,centBin,it);
      // Data
      RooUnfoldBayes unfoldData (rooUnfoldRes[etaBin][centBin], hDataReco[etaBin][centBin], it);    
      hDataUnf[etaBin][centBin][it]= (TH2D*)unfoldData.Hreco();
      hDataUnf[etaBin][centBin][it]->SetName( Form("hDataUnf_eta%d_cent%d_iter%d",etaBin,centBin,it));
    }
  }

  for (int etaBin = 0; etaBin < nEtaBins; +etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      for (int ptBin = 0; ptBin < nPtBins; ++ptBin) {
	hMcTruth_Pt[etaBin][centBin][ptBin] = (TH1D*)hMcTruth[etaBin][centBin]->ProjectionY(Form("hMcTruthPt_eta%d_cent%d_pt%d",etaBin,centBin,ptBin),ptBin,ptBin);
	hMcReco_Pt[etaBin][centBin][ptBin] = (TH1D*)hMcRaw[etaBin][centBin]->ProjectionY(Form("hMcRecoPt_eta%d_cent%d_pt%d",etaBin,centBin,ptBin),ptBin,ptBin);	
	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  hMcUnf_Pt[etaBin][centBin][ptBin][it] = (TH1D*)hMcUnf[etaBin][centBin][it]->ProjectionY(Form("hMcUnfPt_eta%d_cent%d_pt%d_iter%d",etaBin,centBin,ptBin,it],ptBin,ptBin);
	}

	TH1ScaleByWidth(hMcTruth_Pt[etaBin][centBin][ptBin]);
	TH1ScaleByWidth(hMcRaw_Pt[etaBin][centBin][ptBin]);
	TH1ScaleByRapidity(hMcTruth_Pt[etaBin][centBin][ptBin],etaBin);
	TH1ScaleByRapidity(hMcRaw_Pt[etaBin][centBin][ptBin],etaBin);
	
	for ( int it = 0 ; it< nUnfIterations ; it++) {      
	  TH1ScaleByWidth(hMcUnf_Pt[etaBin][centBin][ptBin][it]);
	  TH1ScaleByRapidity(hMcUnf_Pt[etaBin][centBin][ptBin][it], etaBin);
	}

	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  hDataUnf_Pt[etaBin][centBin][ptBin][it] = (TH1D*)hDataUnf[etaBin][centBin][it]->ProjectionY(Form("hdataUnf_Pt_icent%d_pt%d_iter%d",etaBin,centBin,ptBin,it,ptBin,ptBin) );
	}
	hDataReco_Pt[etaBin][centBin][ptBin] = (TH1D*)hDataReco[etaBin][centBin]->ProjectionY(Form("hdataRaw_Pt_eta%d_cent%d_pt%d",etaBin,centBin,ptBin),ptBin,ptBin);
      
	TH1ScaleByWidth(hDataRaw_Pt[etaBin][centBin][ptBin]);
	TH1ScaleByRapidity(hDataRaw_Pt[etaBin][centBin][ptBin], etaBin);
	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  TH1ScaleByWidth(hDataUnf_Pt[etaBin][centBin][ptBin][it]);
	  TH1ScaleByRapidity(hDataUnf_Pt[etaBin][centBin][ptBin][it], etaBin);
	}	       	
      }
    }
  }
  TCanvas *c1; 
  for (int etaBin = 0; etaBin < nEtaBins; +etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      for (int mBin = 0; mBin < nMBins; ++mBin) {
	hMcTruth_M[etaBin][centBin][mBin] = (TH1D*)hMcTruth[etaBin][centBin]->ProjectionY(Form("hMcTruthM_eta%d_cent%d_m%d",etaBin,centBin,mBin),mBin,mBin);
	hMcReco_M[etaBin][centBin][mBin] = (TH1D*)hMcRaw[etaBin][centBin]->ProjectionY(Form("hMcRecoM_eta%d_cent%d_m%d",etaBin,centBin,mBin),mBin,mBin);	
	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  hMcUnf_M[etaBin][centBin][mBin][it] = (TH1D*)hMcUnf[etaBin][centBin]->ProjectionY(Form("hMcUnfM_eta%d_cent%d_m%d_iter%d",etaBin,centBin,mBin,it],mBin,mBin);
	}

	TH1ScaleByWidth(hMcTruth_M[etaBin][centBin][mBin]);
	TH1ScaleByWidth(hMcRaw_M[etaBin][centBin][mBin]);
	TH1ScaleByRapidity(hMcTruth_M[etaBin][centBin][mBin],etaBin);
	TH1ScaleByRapidity(hMcRaw_M[etaBin][centBin][mBin],etaBin);
	for ( int it = 0 ; it< nUnfIterations ; it++) {      
	  TH1ScaleByWidth(hMcUnf_M[etaBin][centBin][mBin][it]);
	  TH1ScaleByRapidity(hMcUnf_M[etaBin][centBin][mBin][it], etaBin);
	  char name[256];
	  sprintf(name,"c1_iter%d",it);
	  c1 = new TCanvas(name,name,500,500);
	  hMcUnf_M[etaBin][centBin][mBin][it]->Draw();
	  c1->SaveAs("");
	}

	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  hDataUnf_M[etaBin][centBin][mBin][it] = (TH1D*)hDataUnf[etaBin][centBin][it]->ProjectionY(Form("hdataUnf_M_icent%d_m%d_iter%d",etaBin,centBin,mBin,it,mBin,mBin) );
	}
	hDataReco_M[etaBin][centBin][mBin] = (TH1D*)hDataReco[etaBin][centBin]->ProjectionY(Form("hdataRaw_M_eta%d_cent%d_m%d",etaBin,centBin,mBin),mBin,mBin);
      
	TH1ScaleByWidth(hDataRaw_M[etaBin][centBin][mBin]);
	TH1ScaleByRapidity(hDataRaw_M[etaBin][centBin][mBin], etaBin);
	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  TH1ScaleByWidth(hDataUnf_M[etaBin][centBin][mBin][it]);
	  TH1ScaleByRapidity(hDataUnf_M[etaBin][centBin][mBin][it], etaBin);	
	}	       	
      }
    }
  }




  TH1D* hStat[30];
  TH1D* hDevi[30];
  for ( int ix = lowPtBin ; ix<= highPtBin ; ix++)  { 
    hStat[ix] = new TH1D(Form("hstat_ix%d",ix),";Number of iterations;Relative uncertainty",29,0.5,29.5);
    hDevi[ix] = (TH1D*)hStat[ix]->Clone(Form("hdevi_ix%d",ix));
  }

  
  TCanvas* c1=  new TCanvas("c1","",1400,400);
  c1->Divide(nPtPannels,1,0,0);
  
  for ( int ipt = lowPtBin ; ipt<= highPtBin ; ipt++)  {

    c1->cd(ipt - lowPtBin + 1);
    
    vector<valErr> vPair;
    valErr nullVal;   nullVal.val = 0 ;  nullVal.err = 0 ;
    vPair.push_back(nullVal);
    
    for (int in = 1; in <= nIter+1 ; in++) {
      vPair.push_back( getDATApoint(kSample, icent, etaBin, ipt, in, matRwt, specRwt, massBin) ) ; 
    }
    for (int in = 2; in <= nIter ; in++) {
      hStat[ipt]->SetBinContent( in, vPair.at(in).err / vPair.at(in).val) ;
      hDevi[ipt]->SetBinContent( in, fabs(vPair.at(in-1).val - vPair.at(in).val)/vPair.at(in).val ) ;
      hStat[ipt]->SetBinError( in, 0.0000001);
      hDevi[ipt]->SetBinError( in, 0.0000001);
    }    
    handsomeTH1(hStat[ipt],4,1,26);
    handsomeTH1(hDevi[ipt],2,1,32);

    //    int fScale = 5;
    //    if ( ipt > lowPtBin + 2)  hStat[ipt]->Scale(1./fScale);
    if ( icent <=3)  hStat[ipt]->SetAxisRange(-0.01,0.03,"Y");
    else if ( icent ==4)  hStat[ipt]->SetAxisRange(-0.50,1.21,"Y");
    else if ( icent ==5)  hStat[ipt]->SetAxisRange(-0.50,1.21,"Y");
    else if ( icent ==6)  hStat[ipt]->SetAxisRange(-1.0,5.1,"Y");
    hStat[ipt]->Draw();
    hDevi[ipt]->Draw("same");


    if ( ipt == lowPtBin ) {
      drawCentrality(kSample, icent, 0.38,0.86,1,24);
      TLegend *leg1 = new TLegend(0.3830109,0.5370667,0.9199772,0.74624,NULL,"brNDC");
      easyLeg(leg1,Form("%.2f < m/p_{T} < %.2f",(float)yBin[massBin-1], (float)yBin[massBin]),0.07);
      leg1->AddEntry(hStat[ipt], "Stat. uncertainty","p");
      leg1->AddEntry(hDevi[ipt], "|y_{N} - y_{N-1}| / y_{N}","p");
      leg1->Draw();
      ATLASLabel(0.2,0.93,"Internal",0.08,0.28);
    }
    drawBin(xBin,ipt,"GeV",0.33 + (0.05* (ipt==lowPtBin)), 0.78,1,20);
    //    if ( ipt > lowPtBin + 2) 
      //      drawText(Form("Stat. Unc scaled by 1/%d",fScale), 0.16 + (0.05* (ipt==lowPtBin)), 0.3,1,16);
    jumSun(0,0,49,0);
  }
  c1->SaveAs(Form("choiceIter/data_coll%d_icent%d_matrixRwt%d_spectraRwt%d_massBin%d.pdf",kSample,icent,(int)matRwt, (int)specRwt,massBin));

  
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(111);
 

  return 0;
}
