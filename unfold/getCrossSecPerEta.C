// Goal: produce ratio plot of Data/MC for dsigma/dpt vs pt. Eta binning integrated from -2.1 < eta < 2.1. Works for pp or PbPb by flipping kSample switch
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
  int nPtBins;
  double ptBins[30];
  int nMBins;
  double mBins[40];
  getXbin(nPtBins, ptBins, 77);
  cout << " nXbins = " << nPtBins << endl;
  getYbin(nMBins, mBins, 772);
  cout << " nMBins = " << nMBins << endl;
  TH1D* hCrossSecData_pt[6][7];
  TH1D* hCrossSecMC_pt[6][7];
  TH1D* hCrossSecRatio_pt[6][7];
  TH1D* hCrossSecData_m[6][7];
  TH1D* hCrossSecMC_m[6][7];
  TH1D* hCrossSecRatio_m[6][7];
  TString fname;
  double histoOffset;
  
  TFile *f1 = TFile::Open(Form("rawSpectra/RawSpectraMC_sample%d.root",kSample));
  TFile *f2 = TFile::Open(Form("rawSpectra/RawSpectraData_sample%d.root",kSample));
   for (int etaBin = 0; etaBin < nEtaIntervals; ++etaBin) {
     for (int centBin = 0; centBin < nCentIntervals; ++centBin) {
       hCrossSecMC_pt[etaBin][centBin] = (TH1D*) f1->Get(Form("hCrossSecMc_pt_eta%dcent%d",etaBin,centBin));
       hCrossSecMC_m[etaBin][centBin] = (TH1D*) f1->Get(Form("hCrossSecMc_m_eta%dcent%d",etaBin,centBin));
       hCrossSecData_pt[etaBin][centBin] = (TH1D*) f2->Get(Form("hCrossSecData_pt_eta%dcent%d",etaBin,centBin));
       hCrossSecData_m[etaBin][centBin] = (TH1D*) f2->Get(Form("hCrossSecData_m_eta%dcent%d",etaBin,centBin));
    }   
  }
  if (kSample == kPbPb) {
    TCanvas *c1;
    TCanvas *c2;
    for (int i=0; i < nEtaIntervals; ++i) {
      c1 = new TCanvas("c1","",500,500);
      makeEfficiencyCanvas(c1,1, 0.0, 0.01, 0.2, 0.25, 0.01);
      c2 = new TCanvas("c2","",500,500);
      makeEfficiencyCanvas(c2,1, 0.0, 0.01, 0.2, 0.25, 0.01);
      for (int j=0; j < nCentIntervals; ++j) {
	c1->cd(1);
	TH1ScaleByWidth(hCrossSecData_pt[i][j]);
	TH1ScaleByWidth(hCrossSecMC_pt[i][j]);
	TH1ScaleByRapidity(hCrossSecData_pt[i][j],i);
	TH1ScaleByRapidity(hCrossSecMC_pt[i][j],i);
	TH1ScaleByWidth(hCrossSecData_m[i][j]);
	TH1ScaleByWidth(hCrossSecMC_m[i][j]);
	TH1ScaleByRapidity(hCrossSecData_m[i][j],i);
	TH1ScaleByRapidity(hCrossSecMC_m[i][j],i);

	// defined in commonUtility: 1 is black; 2 is red    
	if (j == 0) {
	  handsomeTH1(hCrossSecData_pt[i][j],2);
	  handsomeTH1(hCrossSecMC_pt[i][j],2);
	  handsomeTH1(hCrossSecData_m[i][j],2);
	  handsomeTH1(hCrossSecMC_m[i][j],2);	  
	}
	else if (j == 1) {
	  handsomeTH1(hCrossSecData_pt[i][j],3);
	  handsomeTH1(hCrossSecMC_pt[i][j],3);
	  handsomeTH1(hCrossSecData_m[i][j],3);
	  handsomeTH1(hCrossSecMC_m[i][j],3);
	  hCrossSecData_pt[i][j]->Scale(1e+2);
	  hCrossSecMC_pt[i][j]->Scale(1e+2);
	  hCrossSecData_m[i][j]->Scale(1e+2);
	  hCrossSecMC_m[i][j]->Scale(1e+2);
	}
	else if (j == 2) {
	  handsomeTH1(hCrossSecData_pt[i][j],4);
	  handsomeTH1(hCrossSecMC_pt[i][j],4);
	  handsomeTH1(hCrossSecData_m[i][j],4);
	  handsomeTH1(hCrossSecMC_m[i][j],4);
	  hCrossSecData_pt[i][j]->Scale(1e+4);
	  hCrossSecMC_pt[i][j]->Scale(1e+4);
	  hCrossSecData_m[i][j]->Scale(1e+4);
	  hCrossSecMC_m[i][j]->Scale(1e+4);
	}
	else if (j == 3) {
	  handsomeTH1(hCrossSecData_pt[i][j],6);
	  handsomeTH1(hCrossSecMC_pt[i][j],6);
	  handsomeTH1(hCrossSecData_m[i][j],6);
	  handsomeTH1(hCrossSecMC_m[i][j],6);
	  hCrossSecData_pt[i][j]->Scale(1e+6);
	  hCrossSecMC_pt[i][j]->Scale(1e+6);
	  hCrossSecData_m[i][j]->Scale(1e+6);
	  hCrossSecMC_m[i][j]->Scale(1e+6);
	}
	else if (j == 4) {
	  handsomeTH1(hCrossSecData_pt[i][j],7);
	  handsomeTH1(hCrossSecMC_pt[i][j],7);
	  handsomeTH1(hCrossSecData_m[i][j],7);
	  handsomeTH1(hCrossSecMC_m[i][j],7);
	  hCrossSecData_pt[i][j]->Scale(1e+8);
	  hCrossSecMC_pt[i][j]->Scale(1e+8);
	  hCrossSecData_m[i][j]->Scale(1e+8);
	  hCrossSecMC_m[i][j]->Scale(1e+8);
	}
	else if (j == 5) {
	  handsomeTH1(hCrossSecData_pt[i][j],9);
	  handsomeTH1(hCrossSecMC_pt[i][j],9);
	  handsomeTH1(hCrossSecData_m[i][j],9);
	  handsomeTH1(hCrossSecMC_m[i][j],9);
	  hCrossSecData_pt[i][j]->Scale(1e+10);
	  hCrossSecMC_pt[i][j]->Scale(1e+10);
	  hCrossSecData_m[i][j]->Scale(1e+10);
	  hCrossSecMC_m[i][j]->Scale(1e+10);
	}
	else if (j == 6) {
	  handsomeTH1(hCrossSecData_pt[i][j],30);
	  handsomeTH1(hCrossSecMC_pt[i][j],30);
	  handsomeTH1(hCrossSecData_m[i][j],30);
	  handsomeTH1(hCrossSecMC_m[i][j],30);
	  hCrossSecData_pt[i][j]->Scale(1e+12);
	  hCrossSecMC_pt[i][j]->Scale(1e+12);
	  hCrossSecData_m[i][j]->Scale(1e+12);
	  hCrossSecMC_m[i][j]->Scale(1e+12);
	}
	c1->cd(1);
       	hCrossSecMC_pt[i][j]->SetMarkerStyle(kOpenCircle);
	fixedFontHist(hCrossSecData_pt[i][j],2.5,2,20);
	fixedFontHist(hCrossSecMC_pt[i][j],2.5,2,20);
	hCrossSecData_pt[i][j]->SetAxisRange(5e-9,5e17,"Y");
	hCrossSecData_pt[i][j]->SetYTitle("#frac{1}{#LT T_{AA} #GT} #frac{1}{N_{evt}} #frac{d#sigma}{dp_{T}} [nb/GeV]");
	if (j == 0) hCrossSecData_pt[i][j]->Draw();
	else hCrossSecData_pt[i][j]->Draw("same");
	hCrossSecMC_pt[i][j]->Draw("same");
	//	cout << "data eta " << i << hCrossSecData_pt[i][j]->GetEntries() << endl;
	c1->cd(2);
	char ratioName[256];
	sprintf(ratioName,"hCrossSecRatio_pt_eta%d_cent%d",i,j);
	hCrossSecRatio_pt[i][j] = (TH1D*)hCrossSecData_pt[i][j]->Clone(ratioName);
	hCrossSecRatio_pt[i][j]->Divide(hCrossSecMC_pt[i][j]);
	hCrossSecRatio_pt[i][j]->SetXTitle("#it{p}_{T}^{Reco}");
	hCrossSecRatio_pt[i][j]->SetYTitle("Data/MC");
	hCrossSecRatio_pt[i][j]->SetNdivisions(505,"X");
	fixedFontHist(hCrossSecRatio_pt[i][j],3,2,20);
	if (j == 0) hCrossSecRatio_pt[i][j]->Draw();
	else hCrossSecRatio_pt[i][j]->Draw("same");
	
	c2->cd(1);
	hCrossSecMC_m[i][j]->SetMarkerStyle(kOpenCircle);
	fixedFontHist(hCrossSecData_m[i][j],2.5,2,20);
	fixedFontHist(hCrossSecMC_m[i][j],2.5,2,20);
	hCrossSecData_m[i][j]->SetAxisRange(5e-9,5e17,"Y");
	hCrossSecData_m[i][j]->SetYTitle("#frac{1}{#LT T_{AA} #GT} #frac{1}{N_{evt}} #frac{d#sigma}{dp_{T}} [nb/GeV]");
	if (j == 0) hCrossSecData_m[i][j]->Draw();
	else hCrossSecData_m[i][j]->Draw("same");
	hCrossSecMC_m[i][j]->Draw("same");
	//	cout << "data eta " << i << hCrossSecData_m[i][j]->GetEntries() << endl;
	c2->cd(2);

	sprintf(ratioName,"hCrossSecRatio_m_eta%d_cent%d",i,j);
	hCrossSecRatio_m[i][j] = (TH1D*)hCrossSecData_m[i][j]->Clone(ratioName);
	hCrossSecRatio_m[i][j]->Divide(hCrossSecMC_m[i][j]);
	hCrossSecRatio_m[i][j]->SetXTitle("m^{2} / #it{p}_{T}^{2}");
	hCrossSecRatio_m[i][j]->SetYTitle("Data/MC");
	hCrossSecRatio_m[i][j]->SetNdivisions(505,"X");
	fixedFontHist(hCrossSecRatio_m[i][j],3,2,20);
	if (j == 0) hCrossSecRatio_m[i][j]->Draw();
	else hCrossSecRatio_m[i][j]->Draw("same");
      }

      //      gPad->SetLogy();
      c1->cd(1);
      gPad->SetLogy();
      drawCentrality(kSample,-1,0.70,0.86,1,24);
      TLegend *leg1 = new TLegend(0.3271854,0.6394175,0.9049639,0.8533743,NULL,"brNDC");      
      easyLeg(leg1,"");
      leg1->SetNColumns(2);
      leg1->AddEntry(hCrossSecData_pt[i][6], "Data: 60-80% (x10^{12})","plfe");
      leg1->AddEntry(hCrossSecMC_pt[i][6], "MC: 60-80% (x10^{12})","plfe");
      leg1->AddEntry(hCrossSecData_pt[i][5], "Data: 50-60% (x10^{10})","plfe");
      leg1->AddEntry(hCrossSecMC_pt[i][5], "MC: 50-60% (x10^{10})","plfe");
      leg1->AddEntry(hCrossSecData_pt[i][4], "Data: 40-50% (x10^{8})","plfe");
      leg1->AddEntry(hCrossSecMC_pt[i][4], "MC: 40-50% (x10^{8})","plfe");
      leg1->AddEntry(hCrossSecData_pt[i][3], "Data: 30-40% (x10^{6})","plfe");
      leg1->AddEntry(hCrossSecMC_pt[i][3], "MC: 30-40% (x10^{6})","plfe");
      leg1->AddEntry(hCrossSecData_pt[i][2], "Data: 20-30% (x10^{4})","plfe");
      leg1->AddEntry(hCrossSecMC_pt[i][2], "MC: 20-30% (x10^{4})","plfe");
      leg1->AddEntry(hCrossSecData_pt[i][1], "Data: 10-20% (x10^{2})","plfe");
      leg1->AddEntry(hCrossSecMC_pt[i][1], "MC: 10-20% (x10^{2})","plfe");      
      leg1->AddEntry(hCrossSecData_pt[i][0], "Data: 0-10%","plfe");
      leg1->AddEntry(hCrossSecMC_pt[i][0], "MC: 0-10%","plfe");
      leg1->Draw();
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    
      c1->SaveAs(Form("Plots/CrossSection/CrossSection_kSample%d_eta%d_NoReweighting.png",kSample,i));
      c1->SaveAs(Form("Plots/CrossSection/CrossSection_kSample%d_eta%d_NoReweighting.root",kSample,i));
      delete c1;

      c2->cd(1);
      gPad->SetLogy();
      drawCentrality(kSample,-1,0.70,0.86,1,24);
      TLegend *leg2 = new TLegend(0.3271854,0.6394175,0.9049639,0.8533743,NULL,"brNDC");      
      easyLeg(leg2,"");
      leg2->SetNColumns(2);
      leg2->AddEntry(hCrossSecData_m[i][6], "Data: 60-80% (x10^{12})","plfe");
      leg2->AddEntry(hCrossSecMC_m[i][6], "MC: 60-80% (x10^{12})","plfe");
      leg2->AddEntry(hCrossSecData_m[i][5], "Data: 50-60% (x10^{10})","plfe");
      leg2->AddEntry(hCrossSecMC_m[i][5], "MC: 50-60% (x10^{10})","plfe");
      leg2->AddEntry(hCrossSecData_m[i][4], "Data: 40-50% (x10^{8})","plfe");
      leg2->AddEntry(hCrossSecMC_m[i][4], "MC: 40-50% (x10^{8})","plfe");
      leg2->AddEntry(hCrossSecData_m[i][3], "Data: 30-40% (x10^{6})","plfe");
      leg2->AddEntry(hCrossSecMC_m[i][3], "MC: 30-40% (x10^{6})","plfe");
      leg2->AddEntry(hCrossSecData_m[i][2], "Data: 20-30% (x10^{4})","plfe");
      leg2->AddEntry(hCrossSecMC_m[i][2], "MC: 20-30% (x10^{4})","plfe");
      leg2->AddEntry(hCrossSecData_m[i][1], "Data: 10-20% (x10^{2})","plfe");
      leg2->AddEntry(hCrossSecMC_m[i][1], "MC: 10-20% (x10^{2})","plfe");      
      leg2->AddEntry(hCrossSecData_m[i][0], "Data: 0-10%","plfe");
      leg2->AddEntry(hCrossSecMC_m[i][0], "MC: 0-10%","plfe");


      leg2->Draw();
      ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    
      c2->SaveAs(Form("Plots/CrossSection/CrossSection_m2pt2_kSample%d_eta%d_NoReweighting.png",kSample,i));
      c2->SaveAs(Form("Plots/CrossSection/CrossSection_m2pt2_kSample%d_eta%d_NoReweighting.root",kSample,i));
      delete c2;
    } 
  }
  else if (kSample == kPP) {    
    TCanvas *c1 = new TCanvas("c1","",500,500);
    makeEfficiencyCanvas(c1,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    TCanvas *c2 = new TCanvas("c2","",500,500);
    makeEfficiencyCanvas(c2,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c2->cd(1);
    gPad->SetLogy();
    c1->cd(1);
    //    gPad->SetLogx();
    gPad->SetLogy();
    c1->cd(2);

    //    gPad->SetLogx();
    for (int i=0; i < nEtaIntervals; ++i) {
      c1->cd(1);

      TH1ScaleByWidth(hCrossSecData_pt[i][0]);
      TH1ScaleByWidth(hCrossSecMC_pt[i][0]);
      TH1ScaleByWidth(hCrossSecData_m[i][0]);
      TH1ScaleByWidth(hCrossSecMC_m[i][0]);
      TH1ScaleByRapidity(hCrossSecData_pt[i][0], i);
      TH1ScaleByRapidity(hCrossSecMC_pt[i][0], i);
      TH1ScaleByRapidity(hCrossSecData_m[i][0], i);
      TH1ScaleByRapidity(hCrossSecMC_m[i][0], i);
 
      // defined in commonUtility: 1 is black; 2 is red
      if (i == 5) {
	handsomeTH1(hCrossSecData_pt[i][0],kMagenta+1);
	handsomeTH1(hCrossSecMC_pt[i][0],kMagenta+1);
	handsomeTH1(hCrossSecData_m[i][0],kMagenta+1);
	handsomeTH1(hCrossSecMC_m[i][0],kMagenta+1);
	hCrossSecMC_m[i][0]->SetMarkerStyle(kOpenCircle);
	hCrossSecMC_pt[i][0]->SetMarkerStyle(kOpenCircle);
	hCrossSecData_m[i][0]->SetMarkerStyle(kFullCircle);
	hCrossSecData_pt[i][0]->SetMarkerStyle(kFullCircle);
	hCrossSecData_pt[i][0]->Scale(1e+10);
	hCrossSecMC_pt[i][0]->Scale(1e+10);
	hCrossSecData_m[i][0]->Scale(1e+10);
	hCrossSecMC_m[i][0]->Scale(1e+10);
      }
      
      else if (i == 0) {
	handsomeTH1(hCrossSecData_pt[i][0],kGreen+1);
	handsomeTH1(hCrossSecMC_pt[i][0],kGreen+1);
	handsomeTH1(hCrossSecData_m[i][0],kGreen+1);
	handsomeTH1(hCrossSecMC_m[i][0],kGreen+1);
	hCrossSecMC_m[i][0]->SetMarkerStyle(kOpenCross);
	hCrossSecMC_pt[i][0]->SetMarkerStyle(kOpenCross);
	hCrossSecData_m[i][0]->SetMarkerStyle(kFullCross);
	hCrossSecData_pt[i][0]->SetMarkerStyle(kFullCross);
	hCrossSecData_pt[i][0]->Scale(1e+8);
	hCrossSecMC_pt[i][0]->Scale(1e+8);
	hCrossSecData_m[i][0]->Scale(1e+8);
	hCrossSecMC_m[i][0]->Scale(1e+8);
      }
      else if (i == 1) {
	handsomeTH1(hCrossSecData_pt[i][0],kBlue-1);
	handsomeTH1(hCrossSecMC_pt[i][0],kBlue-1);
	handsomeTH1(hCrossSecData_m[i][0],kBlue-1);
	handsomeTH1(hCrossSecMC_m[i][0],kBlue-1);
	hCrossSecMC_m[i][0]->SetMarkerStyle(kOpenDiamond);
	hCrossSecMC_pt[i][0]->SetMarkerStyle(kOpenDiamond);
	hCrossSecData_m[i][0]->SetMarkerStyle(kFullDiamond);
	hCrossSecData_pt[i][0]->SetMarkerStyle(kFullDiamond);
	hCrossSecData_pt[i][0]->Scale(1e+6);
	hCrossSecMC_pt[i][0]->Scale(1e+6);
	hCrossSecData_m[i][0]->Scale(1e+6);
	hCrossSecMC_m[i][0]->Scale(1e+6);
      }
      else if (i == 2) {
	handsomeTH1(hCrossSecData_pt[i][0],kRed+1);
	handsomeTH1(hCrossSecMC_pt[i][0],kRed+1);
	handsomeTH1(hCrossSecData_m[i][0],kRed+1);
	handsomeTH1(hCrossSecMC_m[i][0],kRed+1);
	hCrossSecMC_m[i][0]->SetMarkerStyle(kOpenSquare);
	hCrossSecMC_pt[i][0]->SetMarkerStyle(kOpenSquare);
	hCrossSecData_m[i][0]->SetMarkerStyle(kFullSquare);
	hCrossSecData_pt[i][0]->SetMarkerStyle(kFullSquare);
	hCrossSecData_pt[i][0]->Scale(1e+4);
	hCrossSecMC_pt[i][0]->Scale(1e+4);
	hCrossSecData_m[i][0]->Scale(1e+4);
	hCrossSecMC_m[i][0]->Scale(1e+4);
      }
      else if (i == 3) {
	handsomeTH1(hCrossSecData_pt[i][0],kOrange+5);
	handsomeTH1(hCrossSecMC_pt[i][0],kOrange+5);
	handsomeTH1(hCrossSecData_m[i][0],kOrange+5);
	handsomeTH1(hCrossSecMC_m[i][0],kOrange+5);
	hCrossSecMC_m[i][0]->SetMarkerStyle(kOpenStar);
	hCrossSecMC_pt[i][0]->SetMarkerStyle(kOpenStar);
	hCrossSecData_m[i][0]->SetMarkerStyle(kFullStar);
	hCrossSecData_pt[i][0]->SetMarkerStyle(kFullStar);
	hCrossSecData_pt[i][0]->Scale(1e+2);
	hCrossSecMC_pt[i][0]->Scale(1e+2);
	hCrossSecData_m[i][0]->Scale(1e+2);
	hCrossSecMC_m[i][0]->Scale(1e+2);
      }
      else if (i == 4) {
	handsomeTH1(hCrossSecData_pt[i][0],kPink+1);
	handsomeTH1(hCrossSecMC_pt[i][0],kPink+1);
	handsomeTH1(hCrossSecData_m[i][0],kPink+1);
	handsomeTH1(hCrossSecMC_m[i][0],kPink+1);
	hCrossSecMC_m[i][0]->SetMarkerStyle(kOpenCrossX);
	hCrossSecMC_pt[i][0]->SetMarkerStyle(kOpenCrossX);
	hCrossSecData_m[i][0]->SetMarkerStyle(kFullCrossX);	
	hCrossSecData_pt[i][0]->SetMarkerStyle(kFullCrossX);
	hCrossSecData_pt[i][0]->Scale(1);
	hCrossSecMC_pt[i][0]->Scale(1);
	hCrossSecData_m[i][0]->Scale(1);
	hCrossSecMC_m[i][0]->Scale(1);
      }

      c1->cd(1);
      fixedFontHist(hCrossSecData_pt[i][0],2.5,2,20);
      fixedFontHist(hCrossSecMC_pt[i][0],2.5,2,20);

      hCrossSecData_pt[i][0]->SetAxisRange(9e-9,5e20,"Y");
      hCrossSecData_pt[i][0]->SetAxisRange(125,1000,"X");
      hCrossSecData_pt[i][0]->SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T}dy} [nb/GeV]");
      //      hCrossSecData_pt[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
      hCrossSecMC_pt[i][0]->SetAxisRange(9e-9,5e20,"Y");
      hCrossSecMC_pt[i][0]->SetAxisRange(125,1000,"X");
      hCrossSecMC_pt[i][0]->SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T}dy} [nb/GeV]");
      //      hCrossSecMC_pt[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
      if (i == 0) hCrossSecData_pt[i][0]->Draw();
      else hCrossSecData_pt[i][0]->Draw("same");
      hCrossSecMC_pt[i][0]->Draw("same");
      c1->cd(2);
      char ratioName[256];
      sprintf(ratioName,"hCrossSecRatio_pt_eta%d_cent%d",i,0);
      hCrossSecRatio_pt[i][0] = (TH1D*)hCrossSecData_pt[i][0]->Clone(ratioName);
      hCrossSecRatio_pt[i][0]->Divide(hCrossSecMC_pt[i][0]);
      hCrossSecRatio_pt[i][0]->SetXTitle("#it{p}_{T}^{Reco}");
      hCrossSecRatio_pt[i][0]->SetYTitle("Data/MC");
      //      hCrossSecRatio_pt[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
      hCrossSecRatio_pt[i][0]->SetNdivisions(505,"X");
      hCrossSecRatio_pt[i][0]->SetAxisRange(0.38,0.82,"Y");
      fixedFontHist(hCrossSecRatio_pt[i][0],3,2,20);
      if (i == 0) hCrossSecRatio_pt[i][0]->Draw();
      else hCrossSecRatio_pt[i][0]->Draw("same");

      
      c2->cd(1);
      hCrossSecMC_m[i][0]->SetMarkerStyle(kOpenCircle);
      fixedFontHist(hCrossSecData_m[i][0],2.5,2,20);
      fixedFontHist(hCrossSecMC_m[i][0],2.5,2,20);
      //      hCrossSecData_m[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
       hCrossSecData_m[i][0]->SetAxisRange(9e-9,5e26,"Y");
      hCrossSecData_m[i][0]->SetYTitle("#frac{d^{2}#sigma}{d(m/#it{p}_{T})^{2}dy} [nb]");
      //      hCrossSecMC_m[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
      hCrossSecMC_m[i][0]->SetAxisRange(9e-9,5e26,"Y");
      hCrossSecMC_m[i][0]->SetYTitle("#frac{d^{2}#sigma}{d(m/#it{p}_{T})^{2}}_{T}dy} [nb]"); 
      if (i == 0) hCrossSecData_m[i][0]->Draw();
      else hCrossSecData_m[i][0]->Draw("same");
      hCrossSecMC_m[i][0]->Draw("same");
      c2->cd(2);
      sprintf(ratioName,"hCrossSecRatio_m_eta%d_cent%d",i,0);
      hCrossSecRatio_m[i][0] = (TH1D*)hCrossSecData_m[i][0]->Clone(ratioName);
      hCrossSecRatio_m[i][0]->Divide(hCrossSecMC_m[i][0]);
      hCrossSecRatio_m[i][0]->SetXTitle("m^{2} / #it{p}_{T}^{2}");
      hCrossSecRatio_m[i][0]->SetYTitle("Data/MC");
      hCrossSecRatio_m[i][0]->SetNdivisions(505,"X");
      //      hCrossSecRatio_m[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
      hCrossSecRatio_m[i][0]->SetAxisRange(0,5,"Y");
      fixedFontHist(hCrossSecRatio_m[i][0],3,2,20);
      if (i == 0) hCrossSecRatio_m[i][0]->Draw();
      else hCrossSecRatio_m[i][0]->Draw("same");

    }
    //    c1->cd(2);
    //    gPad->SetLogy();
    c1->cd(1);

    drawCentrality(kSample,0,0.77,0.94,1,24); 
    TLegend *leg1 = new TLegend(0.3571854,0.6194175,0.9849639,0.9333743,NULL,"brNDC");
    easyLeg(leg1,"");
    leg1->SetNColumns(2);
    leg1->AddEntry(hCrossSecData_pt[5][0], "Data: |y| < 2.1 (x10^{10})","plfe");
    leg1->AddEntry(hCrossSecMC_pt[5][0], "MC: |y| < 2.1 (x10^{10})","plfe");
    leg1->AddEntry(hCrossSecData_pt[0][0], "Data: |y| < 0.3 (x10^{8})","plfe");
    leg1->AddEntry(hCrossSecMC_pt[0][0], "MC: |y| < 0.3 (x10^{8})","plfe");    
    leg1->AddEntry(hCrossSecData_pt[1][0], "Data: 0.3 #leq |y| < 0.8 (x10^{6})","plfe");
    leg1->AddEntry(hCrossSecMC_pt[1][0], "MC: 0.3 #leq |y| < 0.8 (x10^{6})","plfe");    
    leg1->AddEntry(hCrossSecData_pt[2][0], "Data: 0.8 #leq |y| < 1.2 (x10^{4})","plfe");
    leg1->AddEntry(hCrossSecMC_pt[2][0], "MC: 0.8 #leq |y| < 1.2 (x10^{4})","plfe");    
    leg1->AddEntry(hCrossSecData_pt[3][0], "Data: 1.2 #leq |y| < 1.6 (x10^{2})","plfe");
    leg1->AddEntry(hCrossSecMC_pt[3][0], "MC: 1.2 #leq |y| < 1.6 (x10^{2})","plfe");    
    leg1->AddEntry(hCrossSecData_pt[4][0], "Data: 1.6 #leq |y| < 2.1","plfe");
    leg1->AddEntry(hCrossSecMC_pt[4][0], "MC: 1.6 #leq |y| < 2.1","plfe");
     
    leg1->Draw();
    ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    c1->cd(2);


    c1->SaveAs(Form("Plots/CrossSection/CrossSection_kSample%d_eta_NoReweighting_ptCut10Events.png",kSample));
    c1->SaveAs(Form("Plots/CrossSection/CrossSection_kSample%d_eta_NoReweighting_ptCut10Events.root",kSample));
    delete c1;


    c2->cd(1);
    //    c2->SetLogx();
    //    c2->SetLogy();
    drawCentrality(kSample,0,0.77,0.94,1,24); 
    TLegend *leg2 = new TLegend(0.3571854,0.6194175,0.9849639,0.9333743,NULL,"brNDC");
    easyLeg(leg2,"");
    leg2->SetNColumns(2);
    leg2->AddEntry(hCrossSecData_m[5][0], "Data: |y| < 2.1 (x10^{10})","plfe");
    leg2->AddEntry(hCrossSecMC_m[5][0], "MC: |y| < 2.1 (x10^{10})","plfe");     
    leg2->AddEntry(hCrossSecData_m[0][0], "Data: |y| < 0.3 (x10^{8})","plfe");
    leg2->AddEntry(hCrossSecMC_m[0][0], "MC: |y| < 0.3 (x10^{8})","plfe");    
    leg2->AddEntry(hCrossSecData_m[1][0], "Data: 0.3 #leq |y| < 0.8 (x10^{6})","plfe");
    leg2->AddEntry(hCrossSecMC_m[1][0], "MC: 0.3 #leq |y| < 0.8 (x10^{6})","plfe");    
    leg2->AddEntry(hCrossSecData_m[2][0], "Data: 0.8 #leq |y| < 1.2 (x10^{4})","plfe");
    leg2->AddEntry(hCrossSecMC_m[2][0], "MC: 0.8 #leq |y| < 1.2 (x10^{4})","plfe");    
    leg2->AddEntry(hCrossSecData_m[3][0], "Data: 1.2 #leq |y| < 1.6 (x10^{2})","plfe");
    leg2->AddEntry(hCrossSecMC_m[3][0], "MC: 1.2 #leq |y| < 1.6 (x10^{2})","plfe");    
    leg2->AddEntry(hCrossSecData_m[4][0], "Data: 1.6 #leq |y| < 2.1","plfe");
    leg2->AddEntry(hCrossSecMC_m[4][0], "MC: 1.6 #leq |y| < 2.1","plfe");
    leg2->Draw();
    ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    c2->SaveAs(Form("Plots/CrossSection/CrossSection_m2pt2_kSample%d_eta_NoReweighting_ptCut10Events.png",kSample));
    c2->SaveAs(Form("Plots/CrossSection/CrossSection_m2pt2_kSample%d_eta_NoReweighting_ptCut10Events.root",kSample));
    delete c2;
  }

  return 0;
}
