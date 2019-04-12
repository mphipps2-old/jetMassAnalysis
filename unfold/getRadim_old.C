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
  gStyle->SetErrorX(0);
  gStyle->SetOptTitle(0);

  TFile *f1 = TFile::Open("unfSpectra/HEPData-ins1673184-v1-root.root");
  TH1D *hRadim[5];
  //  TH1D *hRadim_e1plus[5]; 
  //  TH1D *hRadim_e1minus[5];
  TH1D *hRadim_e2[5];
  //  TFile *f2 = TFile::Open("unfSpectra/UnfoldingSpectra_sample0.root");
  //  TFile *f2 = TFile::Open("unfSpectra/UnfoldingSpectra_sample0_withWeights.root");
  //TFile *f2 = TFile::Open("unfSpectra/UnfoldingSpectra_sample0_withWeights.root");
  TFile *f2 = TFile::Open("unfSpectra/UnfoldingSpectra_sample0_WeightsNoCuts.root");
  TH1D *hMe[5];
  TH1D *hMeNew[5]; // rebin to match Radim's
  TH1D *hRadimNew[5]; // rebin to match what I'm capable of: 12 bins. [63,79],[79,100],[100,126],[126,158],[158,200],[200,251],[251,316],[316,398],[398,501],[501,631],[631,800],[800,1000]
  TH1D *hRadimNew_e1[5]; 
  int nEtaIntervals = 5;
  
  for (int etaBin = 0; etaBin < 5; ++ etaBin) {
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
    cout << " first histo got " << endl;
    cout << " hRadim " << hRadim[etaBin] << endl;
    
  }
  cout << " making canvas " << endl;
  TF1 *fit[5];
  TCanvas *c1 = new TCanvas("c1","",800,500);
  //  makeEfficiencyCanvas(c1,1, 0.0, 0.01, 0.2, 0.25, 0.01);
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
    hRadimNew_e1[i] = new TH1D(name,name,nPtBins,ptBins); 
    for (int j = 2; j < hRadim[0]->GetNbinsX(); ++j ) {
      cout << "radim bin " << j << " lower edge " << hRadim[0]->GetBinLowEdge(1+j) << " bin width  " << hRadim[0]->GetBinWidth(j+1) << endl;
      //   if (j%2 == 1) ++binNum; 
      hRadimNew[i]->SetBinContent(j+1,hRadim[i]->GetBinContent(j+1));
    }
  }
  for (int i = 0; i < 5; ++i ) {
    int binNum = 0;
    char name[256];
    sprintf(name,"Me%d",i);
    hMeNew[i] = new TH1D(name,name,nPtBins,ptBins); 
    for (int j = 0; j < hMe[0]->GetNbinsX(); ++j ) {
      cout << "my bin " << j << " new bin " << binNum << endl;
      hMeNew[i]->SetBinContent(binNum+1,hMe[i]->GetBinContent(j+1));
      if (j%2 == 0) ++binNum; 
    }
  }
  
  
  //    gPad->SetLogx();
  for (int i=0; i < nEtaIntervals; ++i) {
    cout << " eta " << i << endl;
    // defined in commonUtility: 1 is black; 2 is red
    c1->cd(1);
    gPad->SetLogy();
    TH1ScaleByWidth(hMeNew[i]);
    TH1ScaleByRapidity(hMeNew[i],i);
    if (i == 0) {
      cout << " setting marker style " << endl;
      handsomeTH1(hMeNew[i],kMagenta+1);
      handsomeTH1(hRadimNew[i],kMagenta+1);
      hMeNew[i]->SetMarkerStyle(kFullCross);
      cout << " done setting my marker style " << endl;
      hRadimNew[i]->SetMarkerStyle(kOpenCross);
      fit[0] = new TF1("fit0","1.15043e+22/(TMath::Power(x,5.73810))");
      fit[0]->SetLineColor(kMagenta+1);
      cout << " done setting marker style " << endl;
    }
    else if (i == 1) {
      handsomeTH1(hMeNew[i],kGreen+1);
      handsomeTH1(hRadimNew[i],kGreen+1);
      hRadimNew[i]->SetMarkerStyle(kFullDiamond);
      hMeNew[i]->SetMarkerStyle(kOpenDiamond);
      fit[1] = new TF1("fit1","1.18562e+20/(TMath::Power(x,5.75677))");
      fit[1]->SetLineColor(kGreen+1);
    }
    else if (i == 2) {
      handsomeTH1(hMeNew[i],kBlue-1);
      handsomeTH1(hRadimNew[i],kBlue-1);
      hRadimNew[i]->SetMarkerStyle(kFullSquare);
      hMeNew[i]->SetMarkerStyle(kOpenSquare);
      fit[2] = new TF1("fit2","1.61646e+18/(TMath::Power(x,5.85756))");
      fit[2]->SetLineColor(kBlue-1);
    }
    else if (i == 3) {
      handsomeTH1(hMeNew[i],kRed+1);
      handsomeTH1(hRadimNew[i],kRed+1);
      hRadimNew[i]->SetMarkerStyle(kFullStar);
      hMeNew[i]->SetMarkerStyle(kOpenStar);
      fit[3]= new TF1("fit3","2.95101e+16/(TMath::Power(x,6.03325))");
      fit[3]->SetLineColor(kRed+1);
    }
    else if (i == 4) {
      handsomeTH1(hMeNew[i],kOrange+5);
      handsomeTH1(hRadimNew[i],kOrange+5);
      hRadimNew[i]->SetMarkerStyle(kFullCrossX);
      hMeNew[i]->SetMarkerStyle(kOpenCrossX);
      fit[4] = new TF1("fit4","1.76916e+15/(TMath::Power(x,6.46588))");
      fit[4]->SetLineColor(kOrange+5);
    }
    /*
    hMeNew[i]->GetYaxis()->SetTitleOffset(1.1);
    hMeNew[i]->GetXaxis()->SetTitleOffset(1.1);
    hRadimNew[i]->GetYaxis()->SetTitleOffset(1.1);
    hRadimNew[i]->GetXaxis()->SetTitleOffset(1.1);
    */
    fixedFontHist(hMeNew[i],1.1,1.2,20);
    fixedFontHist(hRadimNew[i],1.1,1.2,20);
    if (i == 0) {
      hMeNew[i]->Scale(1e+10);
      hRadimNew[i]->Scale(1e+10);
    }
    else if (i == 1) {
      hMeNew[i]->Scale(1e+8);
      hRadimNew[i]->Scale(1e+8);
    }
    else if (i == 2) {
      hMeNew[i]->Scale(1e+6);
      hRadimNew[i]->Scale(1e+6);
    }
    else if (i == 3) {
      hMeNew[i]->Scale(1e+4);
      hRadimNew[i]->Scale(1e+4);
    }
    else if (i == 4) {
      hMeNew[i]->Scale(1e+2);
      hRadimNew[i]->Scale(1e+2);
    }
    /*
    for (int bin = 1; bin < hRadimNew[i]->GetXaxis()->GetNbins(); ++bin) {
      //      double binError = hRadimNew_e2[i]->GetBinContent(bin);
      cout << " bin " << bin << " i " << i << " binError " << binError << endl;
      
      //      hRadimNew[i]->SetBinError(bin, binError);
    }
    */
    cout << " setting axis range " << endl;
    hMeNew[i]->SetAxisRange(9e-6,5e20,"Y");
    hMeNew[i]->SetAxisRange(125,1000,"X");
    hMeNew[i]->SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T}dy} [nb/GeV]");
    //      hMeNew[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
    hRadimNew[i]->SetAxisRange(9e-6,5e20,"Y");
    hRadimNew[i]->SetAxisRange(125,1000,"X");
    //    hRadimNew[i]->GetYaxis()->SetRange(9e-6,5e20);
    fit[i]->SetRange(100,1000);
    hRadimNew[i]->SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T}dy} [nb/GeV]");
    //      hCrossSec_mc[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
    cout << " drawing " << endl;
    if (i == 0) {
      hMeNew[i]->Draw();
      //      fit[i]->Draw();
    }
    else hMeNew[i]->Draw("same");
    hRadimNew[i]->Draw("same");
    //    c1->cd(2);
        fit[i]->Draw("same");
  }
  

  //    c1->cd(2);
  //    gPad->SetLogy();
  c1->cd(1);

  drawCentrality(kSample,0,0.77,0.88,1,24); 
  TLegend *leg1 = new TLegend(0.2371854,0.6194175,0.8049639,0.9233743,NULL,"brNDC");
  easyLeg(leg1,"");
  leg1->SetNColumns(2);
  
  leg1->AddEntry(hRadimNew[0], "Radim: |y| < 0.3 (x10^{10})","plfe");
  leg1->AddEntry(hMeNew[0], "Me: |y| < 0.3 (x10^{10})","plfe");    
  leg1->AddEntry(hRadimNew[1], "Radim: 0.3 #leq |y| < 0.8 (x10^{8})","plfe");
  leg1->AddEntry(hMeNew[1], "Me: 0.3 #leq |y| < 0.8 (x10^{8})","plfe");    
  leg1->AddEntry(hRadimNew[2], "Radim: 0.8 #leq |y| < 1.2 (x10^{6})","plfe");
  leg1->AddEntry(hMeNew[2], "Me: 0.8 #leq |y| < 1.2 (x10^{6})","plfe");    
  leg1->AddEntry(hRadimNew[3], "Radim: 1.2 #leq |y| < 1.6 (x10^{4})","plfe");
  leg1->AddEntry(hMeNew[3], "Me: 1.2 #leq |y| < 1.6 (x10^{4})","plfe");    
  leg1->AddEntry(hRadimNew[4], "Radim: 1.6 #leq |y| < 2.1 (x10^{2})","plfe");
  leg1->AddEntry(hMeNew[4], "Me: 1.6 #leq |y| < 2.1 (x10^{2})","plfe");
  
  leg1->Draw();
  ATLASLabel(0.25, 0.88, "Internal",0.07,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);


  c1->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparison_kSample%d_eta_ReweightingNoCuts.png",kSample));
  c1->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparison_kSample%d_eta_ReweightingNoCuts.root",kSample));
  delete c1;


  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->Divide(3,2);
  TLegend *leg2[5];
  
  for (int i=0; i < nEtaIntervals; ++i) {
    c2->cd(i+1);
    gPad->SetLogy(0);
    char ratioName[256];
    sprintf(ratioName,"hRatio_radim_eta%d",i);
    hRatio_radim[i] = (TH1D*)hRadimNew[i]->Clone(ratioName);
    hRatio_radim[i]->Divide(fit[i]);
    sprintf(ratioName,"hRatio_me_eta%d",i);
    hRatio_me[i] = (TH1D*)hMeNew[i]->Clone(ratioName);
    hRatio_me[i]->Divide(fit[i]);
    fixedFontHist(hRatio_radim[i],2.3,3.3,13);
    fixedFontHist(hRatio_me[i],2.3,3.3,13);
    if (i == 0) {
      cout << " setting marker style " << endl;
      handsomeTH1(hRatio_me[i],kMagenta+1,1,kFullCross);
      handsomeTH1(hRatio_radim[i],kMagenta+1,1,kOpenCross);
    }
    else if (i == 1) {
      handsomeTH1(hRatio_me[i],kGreen+1,1,kFullDiamond);
      handsomeTH1(hRatio_radim[i],kGreen+1,1,kOpenDiamond);
    }
    else if (i == 2) {
      handsomeTH1(hRatio_me[i],kBlue-1,1,kFullSquare);
      handsomeTH1(hRatio_radim[i],kBlue-1,1, kOpenSquare);
    }
    else if (i == 3) {
      handsomeTH1(hRatio_me[i],kRed+1,1,kFullStar);
      handsomeTH1(hRatio_radim[i],kRed+1,1,kOpenStar);
    }
    else if (i == 4) {
      handsomeTH1(hRatio_me[i],kOrange+5,1,kFullCrossX);
      handsomeTH1(hRatio_radim[i],kOrange+5,1,kOpenCrossX);
    }
    
    hRatio_radim[i]->SetXTitle("#it{p}_{T}^{Reco}");
    hRatio_radim[i]->SetYTitle("Ratio");
    //    hRatio_me[i]->GetYaxis()->SetTitleOffset(1.4);
    //    hRatio_radim[i]->GetYaxis()->SetTitleOffset(1.4);
    hRatio_radim[i]->SetAxisRange(125,1000,"X");
    hRatio_radim[i]->SetNdivisions(505,"X");
    if (i == 0) hRatio_radim[i]->Draw();
    else hRatio_radim[i]->Draw("same");
    hRatio_me[i]->Draw("same");
    
    char legTitle[256];
    if (i == 0)  sprintf(legTitle,"|y| < 0.3");
    else if (i == 1)  sprintf(legTitle,"0.3 #leq |y| < 0.8");
    else if (i == 2)  sprintf(legTitle,"0.8 #leq |y| < 1.2");
    else if (i == 3)  sprintf(legTitle,"1.2 #leq |y| < 1.6");
    else if (i == 4)  sprintf(legTitle,"1.6 #leq |y| < 2.1");
    leg2[i] = new TLegend(0.4471854,0.6194175,0.8649639,0.9233743,NULL,"brNDC");
    easyLeg(leg2[i],legTitle);
    leg2[i]->AddEntry(hRatio_radim[i], "Radim/Fit","plfe");
    leg2[i]->AddEntry(hRatio_me[i], "Me/Fit","plfe");
    leg2[i]->Draw();
  }

  c2->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparisonRatio_kSample%d_eta_ReweightingNoCuts.png",kSample));
  c2->SaveAs(Form("Plots/CrossSection/CrossSectionRadimComparisonRatio_kSample%d_eta_ReweightingNoCuts.root",kSample));
  delete c2;
  
  return 0;

}
