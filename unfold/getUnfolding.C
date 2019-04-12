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
  TFile *f1 = TFile::Open(Form("rawSpectra/RawSpectraData_sample%d.root",kSample)); 
  TFile *f2 = TFile::Open(Form("unfSpectra/rooUnfoldResponse_ReweightingMin0P4Max1P5_sample%d.root",kSample)); 
  int nPtBins;
  double ptBins[30];
  int nMBins;
  double mBins[40];

  //  double binWidth[60];
  // located in ../unfoldingUtil.h
  // Xbin -> pt
  getXbin(nPtBins, ptBins, 80);  
  // Ybin -> m2/pt2
  getYbin(nMBins, mBins, 772);

  long long entries;
  TString fname;

  RooUnfoldResponse* rooUnfoldRes[6][7]; // eta, cent
  TH2D* hDataReco[6][7]; 
  TH2D* hDataUnf[6][7]; // unfolding iter
  //  rooUnfoldRes[0][centBin] = (RooUnfoldResponse*) f2->Get("responseMatrix_eta0cent0");
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {	
      hDataReco[etaBin][centBin] = (TH2D*) f1->Get(Form("hCrossSecData_m_pt_eta%dcent%d",etaBin,centBin));
      rooUnfoldRes[etaBin][centBin] = (RooUnfoldResponse*) f2->Get(Form("responseMatrix_eta%d_cent%d",etaBin,centBin));
      cout << " eta " << etaBin << " cent " << centBin << endl;
    }
  }
   
  cout << " unfold " << endl;
  
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for ( int centBin=0 ; centBin < nCentBins; ++centBin) {
      if ( (kSample == kPP) && ( centBin != 0 ) )      continue;
      // Data
      cout << " etaBin " << etaBin << " centBin " << centBin << endl;
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
  TFile *fout = new TFile(Form("unfSpectra/UnfoldingSpectra_sample%d_WeightsMin0P4Max1P5.root",kSample),"RECREATE");
 
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      for (int mBin = 0; mBin < nMBins; ++mBin) {
	hDataUnf_Pt[etaBin][centBin][mBin] = (TH1D*)hDataUnf[etaBin][centBin]->ProjectionX(Form("hdataUnf_Pt_eta%d_icent%d_m%d",etaBin,centBin,mBin),mBin+2,mBin+2);	 
	hDataUnf_Pt[etaBin][centBin][mBin]->Write();
      }
    }
  }


  
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) { 
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      for (int ptBin = 0; ptBin < nPtBins; ++ptBin) {
	hDataUnf_M[etaBin][centBin][ptBin] = (TH1D*)hDataUnf[etaBin][centBin]->ProjectionY(Form("hdataUnf_M_eta%d_icent%d_pt%d",etaBin,centBin,ptBin),ptBin+2,ptBin+2);
    	hDataUnf_M[etaBin][centBin][ptBin]->Write();
      }
    }
  }

  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {     
      hDataUnf_Pt_Integrated[etaBin][centBin] = (TH1D*)hDataUnf[etaBin][centBin]->ProjectionX(Form("hDataUnfPt_eta%d_cent%d",etaBin,centBin),2,nMBins+1);
      hDataUnf_M_Integrated[etaBin][centBin] = (TH1D*)hDataUnf[etaBin][centBin]->ProjectionY(Form("hDataUnfM_eta%d_cent%d",etaBin,centBin),2,nPtBins+1);
      hDataUnf_Pt_Integrated[etaBin][centBin]->Write();
      hDataUnf_M_Integrated[etaBin][centBin]->Write();
    }
  }
  
  fout->Close();
  
  return 0;
}
