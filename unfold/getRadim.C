// Goal: compare my unfolded result vs radims. Works for pp or PbPb by flipping kSample switch
// KEEP THESE SCRIPTS SMALL AND SINGLE PURPOSED!

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include </home/mike/RootUtils/AtlasStyle.C>
#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"

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
  gStyle->SetErrorX(0.5);
  gStyle->SetOptTitle(0);

  TFile *f1 = TFile::Open("unfSpectra/HEPData-ins1673184-v1-root.root");
  TH1D *hRadim[5];
  //  TH1D *hRadim_e1plus[5]; 
  //  TH1D *hRadim_e1minus[5];
  TH1D *hRadim_e2[5];
  //  TFile *f2 = TFile::Open("unfSpectra/UnfoldingSpectra_sample0.root");
  //  TFile *f2 = TFile::Open("unfSpectra/UnfoldingSpectra_sample0_withWeights.root");
  //TFile *f2 = TFile::Open("unfSpectra/UnfoldingSpectra_sample0_withWeights.root");
  //  TFile *f2 = TFile::Open("unfSpectra/UnfoldingSpectra_sample0_WeightsNoCuts.root");
   TFile *f2 = TFile::Open("unfSpectra/UnfoldingSpectra_sample0_NoWeights.root");
  TFile *f3 = TFile::Open(Form("rawSpectra/RawSpectraData_sample%d.root",kSample));

  TH1D *hMe[5];
  TH1D *hMeNew[5]; // rebin to match Radim's
  TH1D *hRadimNew[5]; // rebin to match what I'm capable of: 12 bins. [63,79],[79,100],[100,126],[126,158],[158,200],[200,251],[251,316],[316,398],[398,501],[501,631],[631,800],[800,1000]
  TH1D *hMyRaw[5];
  
  int nEtaIntervals = 5;
  
  for (int etaBin = 0; etaBin < 5; ++etaBin) {
    hMyRaw[etaBin] = (TH1D*) f3->Get(Form("hCrossSecData_pt_eta%dcent0",etaBin));
    if (etaBin == 0) {
      hRadim[etaBin] = (TH1D*) f1->Get("Table 5/Hist1D_y1");
      //      hRadim_e1plus[etaBin] = (TH1D*) f1->Get("Table 5/Hist1D_y1_e1minus");
      //      hRadim_e1minus[etaBin] = (TH1D*) f1->Get("Table 5/Hist1D_y1_e1minus");
      hRadim_e2[etaBin] = (TH1D*) f1->Get("Table 5/Hist1D_y1_e2");
      hMe[etaBin] = (TH1D*) f2->Get("hDataUnfPt_eta0_cent0");
    }
    else if (etaBin == 1) {
      hRadim[etaBin] = (TH1D*) f1->Get("Table 6/Hist1D_y1");
      hRadim_e2[etaBin] = (TH1D*) f1->Get("Table 6/Hist1D_y1_e2");
      hMe[etaBin] = (TH1D*) f2->Get("hDataUnfPt_eta1_cent0");
    }
    else if (etaBin == 2) {
      hRadim[etaBin] = (TH1D*) f1->Get("Table 7/Hist1D_y1");
      hRadim_e2[etaBin] = (TH1D*) f1->Get("Table 7/Hist1D_y1_e2");
      hMe[etaBin] = (TH1D*) f2->Get("hDataUnfPt_eta2_cent0");
    }
    else if (etaBin == 3) {
      hRadim[etaBin] = (TH1D*) f1->Get("Table 8/Hist1D_y1");
      hRadim_e2[etaBin] = (TH1D*) f1->Get("Table 8/Hist1D_y1_e2");
      hMe[etaBin] = (TH1D*) f2->Get("hDataUnfPt_eta3_cent0");
    }
    else if (etaBin == 4) {
      hRadim[etaBin] = (TH1D*) f1->Get("Table 9/Hist1D_y1");
      hRadim_e2[etaBin] = (TH1D*) f1->Get("Table 9/Hist1D_y1_e2");
      hMe[etaBin] = (TH1D*) f2->Get("hDataUnfPt_eta4_cent0");
    }

    
  }

  TH1D *hRatio_radim[5];
  TH1D *hRatio_me[5];
  int nPtBins;
  double ptBins[30];
  getXbin(nPtBins, ptBins, 81);
  for (int i = 0; i < 5; ++i ) {
    //    int binNum = 1;
    char name[256];
    sprintf(name,"Radim%d",i);
    hRadimNew[i] = new TH1D(name,name,nPtBins,ptBins);
    sprintf(name,"Radim%d_e",i);
    //    hRadimNew_e1[i] = new TH1D(name,name,nPtBins,ptBins); 
    for (int j = 2; j < hRadim[0]->GetNbinsX(); ++j ) {
      //      if (i == 0) cout << "radim bin " << j << " lower edge " << hRadim[0]->GetBinLowEdge(1+j) << " bin width  " << hRadim[0]->GetBinWidth(j+1) << " bin content " << hRadim[i]->GetBinContent(j+1) << endl;
      //   if (j%2 == 1) ++binNum; 
      hRadimNew[i]->SetBinContent(j-1,hRadim[i]->GetBinContent(j+1));
      hRadimNew[i]->SetBinError(j-1,hRadim_e2[i]->GetBinContent(j+1));
    }
  }
  for (int i = 0; i < 5; ++i ) {
    int binNum = 0;
    char name[256];
    sprintf(name,"Me%d",i);
    hMeNew[i] = new TH1D(name,name,nPtBins,ptBins);
    double binContent = 0.;
    for (int j = 0; j < hMe[0]->GetNbinsX(); ++j ) {

      binContent += hMe[i]->GetBinContent(j+1);
      if (j%2 == 1) {
	//	if (i == 0) cout << "my bin " << j <<  " new bin " << binNum << " binContent " << binContent << endl;
	hMeNew[i]->SetBinContent(binNum+1,binContent);
	hMeNew[i]->SetBinError(binNum+1,0.000001);
	++binNum;
	binContent = 0;
      }
    }
  }
  for (int j = 0; j < hRadim[0]->GetNbinsX(); ++j ) {
    cout << "radim bin " << j << " lower edge " << hRadim[0]->GetBinLowEdge(1+j) << " bin width  " << hRadim[0]->GetBinWidth(j+1) << " bin content " << hRadim[0]->GetBinContent(1+j) << endl;
  }
  for (int j = 0; j < hRadimNew[0]->GetNbinsX(); ++j ) {
    cout << "radimNew bin " << j << " lower edge " << hRadimNew[0]->GetBinLowEdge(1+j) << " bin width  " << hRadimNew[0]->GetBinWidth(j+1) << " bin content " << hRadimNew[0]->GetBinContent(1+j) << endl;
  }
  for (int j = 0; j < hMe[0]->GetNbinsX(); ++j ) {
    cout << "me bin " << j << " lower edge " << hMe[0]->GetBinLowEdge(1+j) << " bin width  " << hMe[0]->GetBinWidth(j+1) << " bin content " << hMe[0]->GetBinContent(1+j) << endl;
  }  
  for (int j = 0; j < hMeNew[0]->GetNbinsX(); ++j ) {
    cout << "meNew bin " << j << " lower edge " << hMeNew[0]->GetBinLowEdge(1+j) << " bin width  " << hMeNew[0]->GetBinWidth(j+1) << " bin content " << hMeNew[0]->GetBinContent(1+j) << endl;    
  }

  for (int i=0; i < nEtaIntervals; ++i) {
    TH1ScaleByWidth(hMeNew[i]);
    TH1ScaleByRapidity(hMeNew[i],i);
    TH1ScaleByAcceptance(hMeNew[i],i);
    TH1ScaleByWidth(hMyRaw[i]);
    TH1ScaleByRapidity(hMyRaw[i],i);
    TH1ScaleByAcceptance(hMyRaw[i],i);

    if (i == 0) {
      handsomeTH1(hMeNew[i],kMagenta+1);
      handsomeTH1(hRadimNew[i],kMagenta+1);
      hMeNew[i]->SetMarkerStyle(kFullCross);
      hRadimNew[i]->SetMarkerStyle(kFullCross);
      handsomeTH1(hMyRaw[i],kMagenta+1);
      hMyRaw[i]->SetMarkerStyle(kFullCross);
    }
    else if (i == 1) {
      handsomeTH1(hMeNew[i],kGreen+1);
      handsomeTH1(hRadimNew[i],kGreen+1);
      hRadimNew[i]->SetMarkerStyle(kFullDiamond);
      hMeNew[i]->SetMarkerStyle(kFullDiamond);
      handsomeTH1(hMyRaw[i],kGreen+1);
      hMyRaw[i]->SetMarkerStyle(kFullDiamond);
    }
    else if (i == 2) {
      handsomeTH1(hMeNew[i],kBlue-1);
      handsomeTH1(hRadimNew[i],kBlue-1);
      hRadimNew[i]->SetMarkerStyle(kFullSquare);
      hMeNew[i]->SetMarkerStyle(kFullSquare);
      handsomeTH1(hMyRaw[i],kBlue-1);
      hMyRaw[i]->SetMarkerStyle(kFullSquare);
    }
    else if (i == 3) {
      handsomeTH1(hMeNew[i],kRed+1);
      handsomeTH1(hRadimNew[i],kRed+1);
      hRadimNew[i]->SetMarkerStyle(kFullStar);
      hMeNew[i]->SetMarkerStyle(kFullStar);
      handsomeTH1(hMyRaw[i],kRed+1);
      hMyRaw[i]->SetMarkerStyle(kFullStar);

    }
    else if (i == 4) {
      handsomeTH1(hMeNew[i],kOrange+5);
      handsomeTH1(hRadimNew[i],kOrange+5);
      hRadimNew[i]->SetMarkerStyle(kFullCrossX);
      hMeNew[i]->SetMarkerStyle(kFullCrossX);
      handsomeTH1(hMyRaw[i],kOrange+5);
      hMyRaw[i]->SetMarkerStyle(kFullCrossX);
    }
    fixedFontHist(hMeNew[i],1.1,2.3,20);
    fixedFontHist(hRadimNew[i],1.1,2.3,20);
    fixedFontHist(hMyRaw[i],1.1,2.3,20);
  }
  
  TCanvas *c0 = new TCanvas("c0","",1100,500);
  c0->Divide(3,1);
  for (int i=0; i < nEtaIntervals; ++i) {
    c0->cd(1);
    gPad->SetLogy();
    gStyle->SetErrorX(0.5);
    hMeNew[i]->GetXaxis()->SetRangeUser(126,1000);
    hMeNew[i]->SetXTitle("#it{p}_{T}");
    hMeNew[i]->SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T}dy} [nb/GeV]");
    if (i == 0) {
      hMeNew[i]->Draw("P");
    }
    else hMeNew[i]->Draw("sameP");

    c0->cd(2);
    gPad->SetLogy();
    hRadimNew[i]->GetXaxis()->SetRangeUser(126,1000);
    hRadimNew[i]->SetXTitle("#it{p}_{T}");
    hRadimNew[i]->SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T}dy} [nb/GeV]");
    if (i == 0) {
      hRadimNew[i]->Draw("P");
    }
    else hRadimNew[i]->Draw("sameP");
  
    c0->cd(3);
    gPad->SetLogy();
    hMyRaw[i]->GetXaxis()->SetRangeUser(126,1000);
    hMyRaw[i]->SetXTitle("#it{p}_{T}");
    hMyRaw[i]->SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T}dy} [nb/GeV]");
    if (i == 0) {
      hMyRaw[i]->Draw("P");
    }
    else hMyRaw[i]->Draw("sameP");
  }
  
  TLegend *leg1 = new TLegend(0.3571854,0.6194175,0.8649639,0.8733743,NULL,"brNDC");
  easyLeg(leg1,"My Unfolded");
  c0->cd(1); 
  leg1->AddEntry(hMeNew[0], "|y| < 0.3","plfe");    
  leg1->AddEntry(hMeNew[1], "0.3 #leq |y| < 0.8","plfe");    
  leg1->AddEntry(hMeNew[2], "0.8 #leq |y| < 1.2","plfe");    
  leg1->AddEntry(hMeNew[3], "1.2 #leq |y| < 1.6","plfe");    
  leg1->AddEntry(hMeNew[4], "1.6 #leq |y| < 2.1","plfe");  
  leg1->Draw();
  ATLASLabel(0.25, 0.88, "Internal",0.07,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
  drawCentrality(kSample,0,0.77,0.88,1,24);
  
  c0->cd(2);
  TLegend *leg2 = new TLegend(0.3571854,0.6194175,0.8649639,0.8733743,NULL,"brNDC");
  easyLeg(leg2,"Radim Unfolded");
  leg2->AddEntry(hRadimNew[0], "|y| < 0.3","plfe");
   leg2->AddEntry(hRadimNew[1], "0.3 #leq |y| < 0.8","plfe");
  leg2->AddEntry(hRadimNew[2], "0.8 #leq |y| < 1.2","plfe");
  leg2->AddEntry(hRadimNew[3], "1.2 #leq |y| < 1.6","plfe");
  leg2->AddEntry(hRadimNew[4], "1.6 #leq |y| < 2.1","plfe");
  leg2->Draw();

  c0->cd(3);
  TLegend *leg3 = new TLegend(0.3571854,0.6194175,0.8649639,0.8733743,NULL,"brNDC");
  easyLeg(leg3,"My Raw");
  leg3->AddEntry(hRadimNew[0], "|y| < 0.3","plfe");
  leg3->AddEntry(hRadimNew[1], "0.3 #leq |y| < 0.8","plfe");
  leg3->AddEntry(hRadimNew[2], "0.8 #leq |y| < 1.2","plfe");
  leg3->AddEntry(hRadimNew[3], "1.2 #leq |y| < 1.6","plfe");
  leg3->AddEntry(hRadimNew[4], "1.6 #leq |y| < 2.1","plfe");
  leg3->Draw();



  //    c0->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparisonNoYScale_kSample%d_eta_ReweightingNoCuts.png",kSample));
  //  c0->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparisonNoYScale_kSample%d_eta_ReweightingNoCuts.root",kSample));
  c0->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparisonNoYScale_kSample%d_eta_NoReweighting.png",kSample));
  c0->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparisonNoYScale_kSample%d_eta_NoReweighting.root",kSample));
  delete c0;

  
  TCanvas *c1 = new TCanvas("c1","",500,500);
  makeEfficiencyCanvas(c1,1, 0.0, 0.01, 0.2, 0.25, 0.01);
  //  gStyle->SetErrorX(0);  
  for (int i=0; i < nEtaIntervals; ++i) {
    // defined in commonUtility: 1 is black; 2 is red
    c1->cd(1);
    gPad->SetLogy();
    gStyle->SetErrorX(0.5);
    if (i == 0) {
      hMeNew[i]->Scale(1e+10);
      hRadimNew[i]->Scale(1e+10);

      handsomeTH1(hMeNew[i],kMagenta+1);
      handsomeTH1(hRadimNew[i],kMagenta+1);
      hMeNew[i]->SetMarkerStyle(kFullCross);
      hRadimNew[i]->SetMarkerStyle(kOpenCross);
    }
    else if (i == 1) {
      hMeNew[i]->Scale(1e+8);
      hRadimNew[i]->Scale(1e+8);

      handsomeTH1(hMeNew[i],kGreen+1);
      handsomeTH1(hRadimNew[i],kGreen+1);
      hRadimNew[i]->SetMarkerStyle(kFullDiamond);
      hMeNew[i]->SetMarkerStyle(kOpenDiamond);
    }
    else if (i == 2) {
      hMeNew[i]->Scale(1e+6);
      hRadimNew[i]->Scale(1e+6);

      handsomeTH1(hMeNew[i],kBlue-1);
      handsomeTH1(hRadimNew[i],kBlue-1);
      hRadimNew[i]->SetMarkerStyle(kFullSquare);
      hMeNew[i]->SetMarkerStyle(kOpenSquare);
    }
    else if (i == 3) {
      hMeNew[i]->Scale(1e+4);
      hRadimNew[i]->Scale(1e+4);

      handsomeTH1(hMeNew[i],kRed+1);
      handsomeTH1(hRadimNew[i],kRed+1);
      hRadimNew[i]->SetMarkerStyle(kFullStar);
      hMeNew[i]->SetMarkerStyle(kOpenStar);
    }
    else if (i == 4) {
      hMeNew[i]->Scale(1e+2);
      hRadimNew[i]->Scale(1e+2);
      
      handsomeTH1(hMeNew[i],kOrange+5);
      handsomeTH1(hRadimNew[i],kOrange+5);
      hRadimNew[i]->SetMarkerStyle(kFullCrossX);
      hMeNew[i]->SetMarkerStyle(kOpenCrossX);
    }

    /*
    for (int bin = 1; bin < hRadimNew[i]->GetXaxis()->GetNbins(); ++bin) {
      //      double binError = hRadimNew_e2[i]->GetBinContent(bin);
      cout << " bin " << bin << " i " << i << " binError " << binError << endl;
       
      //      hRadimNew[i]->SetBinError(bin, binError);
    }
    */
 
    hMeNew[i]->SetAxisRange(9e-6,5e20,"Y");
    hMeNew[i]->GetXaxis()->SetRangeUser(126,1000);
    hMeNew[i]->SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T}dy} [nb/GeV]");
    //      hMeNew[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
    //    hRadimNew[i]->SetAxisRange(9e-6,5e20,"Y");
    hRadimNew[i]->GetXaxis()->SetRangeUser(126,1000);
    //    hRadimNew[i]->GetYaxis()->SetRange(9e-6,5e20);
    //    fit[i]->SetRange(100,1000);
    hRadimNew[i]->SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T}dy} [nb/GeV]");
    //      hCrossSec_mc[i]->GetXaxis()->SetMoreLogLabels(kTRUE);

    if (i == 0) {
      hMeNew[i]->Draw("P");
      //      fit[i]->Draw();
    }
    else hMeNew[i]->Draw("sameP");
    hRadimNew[i]->Draw("sameP");
    //    fit[i]->Draw("same");    
    c1->cd(2);
    gStyle->SetErrorX(0.5);
       gPad->SetLogy(0);
    char ratioName[256];
    sprintf(ratioName,"hRatio_radim_eta%d",i);
    //    hRatio_radim[i] = (TH1D*)hRadimNew[i]->Clone(ratioName);
    hRatio_radim[i] = (TH1D*)hMeNew[i]->Clone(ratioName);

    hRatio_radim[i]->Divide(hRadimNew[i]);

    fixedFontHist(hRatio_radim[i],2.3,3.3,13);
    if (i == 0) {

      handsomeTH1(hRatio_radim[i],kMagenta+1,1,kFullCross);
    }
    else if (i == 1) {
      handsomeTH1(hRatio_radim[i],kGreen+1,1,kFullDiamond);
    }
    else if (i == 2) {
      handsomeTH1(hRatio_radim[i],kBlue-1,1, kFullSquare);
    }
    else if (i == 3) {
      handsomeTH1(hRatio_radim[i],kRed+1,1,kFullStar);
    }
    else if (i == 4) {      
      handsomeTH1(hRatio_radim[i],kOrange+5,1,kFullCrossX);
    }
    
    hRatio_radim[i]->SetXTitle("#it{p}_{T} (GeV)");
    hRatio_radim[i]->SetYTitle("Ratio (Me/Radim)");
    hRatio_radim[i]->GetXaxis()->SetRangeUser(126,1000);
    hRatio_radim[i]->SetAxisRange(0.75,1.25,"Y");
    hRatio_radim[i]->SetNdivisions(505,"X");
    if (i == 0) hRatio_radim[i]->Draw("P");
    else hRatio_radim[i]->Draw("sameP");
   
  }
  

  //    c1->cd(2);
  //    gPad->SetLogy();
  c1->cd(1);

  drawCentrality(kSample,0,0.77,0.88,1,24); 
  TLegend *leg0 = new TLegend(0.2371854,0.6194175,0.8049639,0.9233743,NULL,"brNDC");
  easyLeg(leg0,"");
  leg0->SetNColumns(2);
  
  leg0->AddEntry(hRadimNew[0], "Radim: |y| < 0.3 (x10^{10})","plfe");
  leg0->AddEntry(hMeNew[0], "Me: |y| < 0.3 (x10^{10})","plfe");    
  leg0->AddEntry(hRadimNew[1], "Radim: 0.3 #leq |y| < 0.8 (x10^{8})","plfe");
  leg0->AddEntry(hMeNew[1], "Me: 0.3 #leq |y| < 0.8 (x10^{8})","plfe");    
  leg0->AddEntry(hRadimNew[2], "Radim: 0.8 #leq |y| < 1.2 (x10^{6})","plfe");
  leg0->AddEntry(hMeNew[2], "Me: 0.8 #leq |y| < 1.2 (x10^{6})","plfe");    
  leg0->AddEntry(hRadimNew[3], "Radim: 1.2 #leq |y| < 1.6 (x10^{4})","plfe");
  leg0->AddEntry(hMeNew[3], "Me: 1.2 #leq |y| < 1.6 (x10^{4})","plfe");    
  leg0->AddEntry(hRadimNew[4], "Radim: 1.6 #leq |y| < 2.1 (x10^{2})","plfe");
  leg0->AddEntry(hMeNew[4], "Me: 1.6 #leq |y| < 2.1 (x10^{2})","plfe");
  
  leg0->Draw();
  ATLASLabel(0.25, 0.88, "Internal",0.07,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);


  //  c1->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparison_kSample%d_eta_ReweightingNoCuts.png",kSample));
  //  c1->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparison_kSample%d_eta_ReweightingNoCuts.root",kSample));
  c1->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparison_kSample%d_eta_NoReweighting.png",kSample));
  c1->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparison_kSample%d_eta_NoReweighting.root",kSample));
  delete c1;
  
  return 0;

}
