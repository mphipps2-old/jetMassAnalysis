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
#include "TFile.h"
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
  gStyle->SetPalette(kRainBow);
  int nCentBins;
  if (kSample == kPP) nCentBins = 1;
  else nCentBins = 7; 
  int nEtaBins = 6;  
  int nUnfIterations = 6;

  int nPtBins;
  double ptBins[30];
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

  TH1D* hCrossSecData_pt[6][7];
  TH1D* hCrossSecMC_pt[6][7];
  TH1D* hMcOverData_pt[6][7];
  TH1D* hCrossSecData_m[6][7];
  TH1D* hCrossSecMC_m[6][7];
  TH1D* hMcOverData_m[6][7];
  TH2D* hCrossSecData[6][7];
  TH2D* hCrossSecData_stats[6][7];
  TH2D* hCrossSecMC[6][7];
  TH2D* hCrossSecMC_stats[6][7];
  TH2D* hMcOverData[6][7];
  TH2D* hMcOverData_stats[6][7];
  TFile *f1 = TFile::Open(Form("rawSpectra/RawSpectraMC_sample%d.root",kSample));
  // TFile *f1 = TFile::Open(Form("rawSpectra/RawSpectraMC_sample%d_ptCut10Events.root",kSample));
  //   TFile *f1 = TFile::Open(Form("rawSpectra/RawSpectraMC_sample%d_ptCut100Events.root",kSample));
  TFile *f2 = TFile::Open(Form("rawSpectra/RawSpectraData_sample%d.root",kSample));
  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      hCrossSecMC_pt[etaBin][centBin] = (TH1D*) f1->Get(Form("hCrossSecMc_pt_eta%dcent%d",etaBin,centBin));
      hCrossSecMC_m[etaBin][centBin] = (TH1D*) f1->Get(Form("hCrossSecMc_m_eta%dcent%d",etaBin,centBin));
      hCrossSecMC[etaBin][centBin] = (TH2D*) f1->Get(Form("hCrossSecMc_m_pt_eta%dcent%d",etaBin,centBin));
      hCrossSecMC_stats[etaBin][centBin] = (TH2D*) f1->Get(Form("hCrossSecMc_m_pt_lowPt_stats_eta%dcent%d",etaBin,centBin));
      hCrossSecData_pt[etaBin][centBin] = (TH1D*) f2->Get(Form("hCrossSecData_pt_eta%dcent%d",etaBin,centBin));
      hCrossSecData_m[etaBin][centBin] = (TH1D*) f2->Get(Form("hCrossSecData_m_eta%dcent%d",etaBin,centBin));
      hCrossSecData[etaBin][centBin] = (TH2D*) f2->Get(Form("hCrossSecData_m_pt_eta%dcent%d",etaBin,centBin));
      hCrossSecData_stats[etaBin][centBin] = (TH2D*) f2->Get(Form("hCrossSecData_m_pt_lowPt_stats_eta%dcent%d",etaBin,centBin));
    }   
  }
  
  
  gStyle->SetPalette(kBird);
  char outName[256];
  sprintf(outName,"ReweightingFactors/ReweightingFactors_kSample%d_DataOverMc.root",kSample);
  
  TFile *fout = new TFile(outName,"recreate");

  for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      
      hMcOverData[etaBin][centBin] = (TH2D*) hCrossSecData[etaBin][centBin]->Clone(Form("ReweightingHisto_etaBin%d_centBin%d",etaBin,centBin));
      hMcOverData[etaBin][centBin]->Divide(hCrossSecMC[etaBin][centBin]);
      hMcOverData[etaBin][centBin]->GetXaxis()->SetTitle("#it{p}_{T}");      
      hMcOverData[etaBin][centBin]->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");
      hMcOverData[etaBin][centBin]->GetYaxis()->SetTitleOffset(1.55);
      hMcOverData[etaBin][centBin]->SetContour(99);
      drawCentrality(kSample,0,0.70,0.86,1,24);
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
      hMcOverData[etaBin][centBin]->Write("",TObject::kOverwrite);

      hMcOverData_stats[etaBin][centBin] = (TH2D*) hCrossSecData_stats[etaBin][centBin]->Clone(Form("ReweightingHisto_stats_etaBin%d_centBin%d",etaBin,centBin));
      hMcOverData_stats[etaBin][centBin]->Divide(hCrossSecMC_stats[etaBin][centBin]);
      hMcOverData_stats[etaBin][centBin]->GetXaxis()->SetTitle("#it{p}_{T}");      
      hMcOverData_stats[etaBin][centBin]->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");
      hMcOverData_stats[etaBin][centBin]->GetYaxis()->SetTitleOffset(1.55);
      hMcOverData_stats[etaBin][centBin]->SetContour(99);
      drawCentrality(kSample,0,0.70,0.86,1,24);
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
      hMcOverData_stats[etaBin][centBin]->Write("",TObject::kOverwrite);
   
      hMcOverData_m[etaBin][centBin] = (TH1D*) hCrossSecData_m[etaBin][centBin]->Clone(Form("ReweightingHisto_M_etaBin%d_centBin%d",etaBin,centBin));
      hMcOverData_m[etaBin][centBin]->Divide(hCrossSecMC_m[etaBin][centBin]);
      hMcOverData_m[etaBin][centBin]->GetXaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");
      hMcOverData_m[etaBin][centBin]->GetYaxis()->SetTitle("Ratio (Data/MC)");
      hMcOverData_m[etaBin][centBin]->GetYaxis()->SetTitleOffset(1.55);
      drawCentrality(kSample,0,0.70,0.86,1,24);
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
      hMcOverData_m[etaBin][centBin]->Write("",TObject::kOverwrite);
      /*
	std::string label = "()";
	char name2[256];
	sprintf(name2,"Title:Reweighting");
	c6->Print(Form("Unfolding/Closure/Reweighting_sample%d_icent%d_etaBin%d.pdf",kSample,0,5), name2); 
      */


      hMcOverData_pt[etaBin][centBin] = (TH1D*) hCrossSecData_pt[etaBin][centBin]->Clone(Form("ReweightingHisto_Pt_etaBin%d_centBin%d",etaBin,centBin));
      hMcOverData_pt[etaBin][centBin]->Divide(hCrossSecMC_pt[etaBin][centBin]);
      hMcOverData_pt[etaBin][centBin]->GetXaxis()->SetTitle("#it{p}_{T}");
      hMcOverData_pt[etaBin][centBin]->GetYaxis()->SetTitle("Ratio (Data/MC)");
      drawCentrality(kSample,0,0.70,0.86,1,24);
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
      hMcOverData_pt[etaBin][centBin]->Write("",TObject::kOverwrite);
    }
  }
  fout->Close();

  return 0;
}
