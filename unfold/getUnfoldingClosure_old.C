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
  int nUnfIterations = 4;


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
  
  //  RooUnfoldResponse* hRooUnfold[6][7]; // eta, cent
  RooUnfoldResponse* rooUnfoldRes[6][7]; // eta, cent
  TH2D *hResponseMatrix_pt[6][7];
  TH2D *hResponseMatrix_m[6][7];
  TH2D *hMcReco1[6][7];
  TH2D *hMcReco2[6][7];
  TH2D *hMcTruth[6][7];
  TH2D *hResMatrix[6][7];
  TH2D* hMcUnf[6][7][20]; // unfolding iter
  TH2D* hDataUnf[6][7][20]; // unfolding iter
 
  
  for (int etaBin = 5; etaBin < nEtaBins; ++etaBin) {
    for ( int centBin=0 ; centBin < nCentBins; ++centBin) {
      if ( (kSample == kPP) && ( centBin != 0 ) )      continue;
      // loop through unfolding iterations
      for ( int it = 0 ; it < nUnfIterations; it++) {  
	RooUnfoldBayes unfoldMc (rooUnfoldRes[etaBin][centBin], hMcReco2[etaBin][centBin], it);    
	hMcUnf[etaBin][centBin][it] = (TH2D*)unfoldMc.Hreco();
	hMcUnf[etaBin][centBin][it]->SetName( Form("hMcUnf_eta%d_cent%d_iter%d",etaBin,centBin,it));
	// Data
	RooUnfoldBayes unfoldData (rooUnfoldRes[etaBin][centBin], hDataReco[etaBin][centBin], it);    
	hDataUnf[etaBin][centBin][it]= (TH2D*)unfoldData.Hreco();
	hDataUnf[etaBin][centBin][it]->SetName( Form("hDataUnf_eta%d_cent%d_iter%d",etaBin,centBin,it));
      }
    }
  }


  
  TCanvas *c4 = new TCanvas("c4","c4");
   
   hDataUnf[5][0][1]->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");      
   hDataUnf[5][0][1]->GetXaxis()->SetTitle("#it{p}_{T}");
   hDataUnf[5][0][1]->GetYaxis()->SetTitleOffset(1.55);
   hDataUnf[5][0][1]->SetContour(99);
   hDataUnf[5][0][1]->Draw("colz");
   c4->SetLogz();
   //   hDataUnf[5][0][1]->SetMinimum(1e-5);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   label = "";
   sprintf(name2,"Title:dataUnf");
   c4->Print(Form("Unfolding/Closure/Closure_diagnostics_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, label.c_str()), name2);

   TCanvas *c5 = new TCanvas("c5","c5");
   hMcUnf[5][0][1]->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");      
   hMcUnf[5][0][1]->GetXaxis()->SetTitle("#it{p}_{T}");
   hMcUnf[5][0][1]->GetYaxis()->SetTitleOffset(1.55);
   hMcUnf[5][0][1]->SetContour(99);
   hMcUnf[5][0][1]->Draw("colz");
   c5->SetLogz();
   //   hMcUnf[5][0][1]->SetMinimum(1e-5);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   label = ")";
   sprintf(name2,"Title:mcUnf");
   c5->Print(Form("Unfolding/Closure/Closure_diagnostics_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, label.c_str()), name2);
   
  TH1D *hMcTruth_Pt[6][7][30];
  TH1D *hMcReco_Pt[6][7][30];
  TH1D *hMcUnf_Pt[6][7][30][20]; // nUnfIterations
  TH1D *hMcRatio_Pt[6][7][30][20];
  TH1D *hMcTruth_M[6][7][20];
  TH1D *hMcReco_M[6][7][20];
  TH1D *hMcUnf_M[6][7][20][20];
  TH1D *hMcRatio_M[6][7][20][20];
  TH1D *hDataReco_Pt[6][7][30];
  TH1D *hDataUnf_Pt[6][7][30][20];
  TH1D *hDataRatio_Pt[6][7][30][20];
  TH1D *hDataReco_M[6][7][20];
  TH1D *hDataUnf_M[6][7][20][20];
  TH1D *hDataRatio_M[6][7][20][20];
    
  for (int etaBin = 5; etaBin < nEtaBins; ++etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      for (int mBin = 0; mBin < nMBins; ++mBin) {
	hMcTruth_Pt[etaBin][centBin][mBin] = (TH1D*)hMcTruth[etaBin][centBin]->ProjectionX(Form("hMcTruthPt_eta%d_cent%d_m%d",etaBin,centBin,mBin),mBin+2,mBin+2);
	hMcReco_Pt[etaBin][centBin][mBin] = (TH1D*)hMcReco2[etaBin][centBin]->ProjectionX(Form("hMcRecoPt_eta%d_cent%d_m%d",etaBin,centBin,mBin),mBin+2,mBin+2);	
	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  hMcUnf_Pt[etaBin][centBin][mBin][it] = (TH1D*)hMcUnf[etaBin][centBin][it]->ProjectionX(Form("hMcUnfPt_eta%d_cent%d_m%d_iter%d",etaBin,centBin,mBin,it),mBin+2,mBin+2);
	}

	TH1ScaleByWidth(hMcTruth_Pt[etaBin][centBin][mBin]);
	TH1ScaleByWidth(hMcReco_Pt[etaBin][centBin][mBin]);
	TH1ScaleByRapidity(hMcTruth_Pt[etaBin][centBin][mBin],etaBin);
	TH1ScaleByRapidity(hMcReco_Pt[etaBin][centBin][mBin],etaBin);
	
	for ( int it = 0 ; it< nUnfIterations ; it++) {      
	   TH1ScaleByWidth(hMcUnf_Pt[etaBin][centBin][mBin][it]);
	   TH1ScaleByRapidity(hMcUnf_Pt[etaBin][centBin][mBin][it], etaBin);
	}

	hDataReco_Pt[etaBin][centBin][mBin] = (TH1D*)hDataReco[etaBin][centBin]->ProjectionX(Form("hdataRaw_Pt_eta%d_cent%d_m%d",etaBin,centBin,mBin),mBin+2,mBin+2);		
	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  hDataUnf_Pt[etaBin][centBin][mBin][it] = (TH1D*)hDataUnf[etaBin][centBin][it]->ProjectionX(Form("hdataUnf_Pt_eta%d_icent%d_m%d_iter%d",etaBin,centBin,mBin,it),mBin+2,mBin+2);
	}      
	TH1ScaleByWidth(hDataReco_Pt[etaBin][centBin][mBin]);
	TH1ScaleByRapidity(hDataReco_Pt[etaBin][centBin][mBin], etaBin);
	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  TH1ScaleByWidth(hDataUnf_Pt[etaBin][centBin][mBin][it]);
	  TH1ScaleByRapidity(hDataUnf_Pt[etaBin][centBin][mBin][it], etaBin);
	}	       	
      }
    }
  }


  
  for (int etaBin = 5; etaBin < nEtaBins; ++etaBin) { 
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      for (int ptBin = 0; ptBin < nPtBins; ++ptBin) {
	hMcTruth_M[etaBin][centBin][ptBin] = (TH1D*)hMcTruth[etaBin][centBin]->ProjectionY(Form("hMcTruthM_eta%d_cent%d_pt%d",etaBin,centBin,ptBin),ptBin+2,ptBin+2);
	hMcReco_M[etaBin][centBin][ptBin] = (TH1D*)hMcReco2[etaBin][centBin]->ProjectionY(Form("hMcRecoM_eta%d_cent%d_pt%d",etaBin,centBin,ptBin),ptBin+2,ptBin+2);	
	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  hMcUnf_M[etaBin][centBin][ptBin][it] = (TH1D*)hMcUnf[etaBin][centBin][it]->ProjectionY(Form("hMcUnfM_eta%d_cent%d_pt%d_iter%d",etaBin,centBin,ptBin,it),ptBin+2,ptBin+2);
	  //	  TH1ScaleByWidth(hMcUnf_M[etaBin][centBin][ptBin][it]);
	  //	  TH1ScaleByRapidity(hMcUnf_M[etaBin][centBin][ptBin][it], etaBin);
	}
	//	TH1ScaleByWidth(hMcTruth_M[etaBin][centBin][ptBin]);
	//	TH1ScaleByWidth(hMcReco_M[etaBin][centBin][ptBin]);
	//	TH1ScaleByRapidity(hMcTruth_M[etaBin][centBin][ptBin],etaBin);
	//	TH1ScaleByRapidity(hMcReco_M[etaBin][centBin][ptBin],etaBin);
	cout << " eta " << etaBin << " cent " << centBin << " ptBin " << ptBin << endl;
	hDataReco_M[etaBin][centBin][ptBin] = (TH1D*)hDataReco[etaBin][centBin]->ProjectionY(Form("hdataRaw_M_eta%d_cent%d_pt%d",etaBin,centBin,ptBin),ptBin+2,ptBin+2);      
	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  hDataUnf_M[etaBin][centBin][ptBin][it] = (TH1D*)hDataUnf[etaBin][centBin][it]->ProjectionY(Form("hdataUnf_M_eta%d_icent%d_pt%d_iter%d",etaBin,centBin,ptBin,it),ptBin+2,ptBin+2);
	  //	  TH1ScaleByWidth(hDataUnf_M[etaBin][centBin][ptBin][it]);
	  //	  TH1ScaleByRapidity(hDataUnf_M[etaBin][centBin][ptBin][it], etaBin);	
	}
	//	TH1ScaleByWidth(hDataReco_M[etaBin][centBin][ptBin]);
	//	TH1ScaleByRapidity(hDataReco_M[etaBin][centBin][ptBin], etaBin);
      }
    }
  }


  TCanvas* c1a=  new TCanvas("c1a","",1200,550);
  TCanvas* c1b=  new TCanvas("c1b","",1200,550);
  TCanvas* c1c=  new TCanvas("c1c","",1200,550);
  TCanvas* c1d=  new TCanvas("c1d","",1200,550);
  TCanvas* c1e=  new TCanvas("c1e","",1200,550);

  makeEfficiencyCanvas(c1a,5+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c1b,5+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c1c,5+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c1d,5+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  //  makeEfficiencyCanvas(c1e,nPtBins+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  // place legend in final pannel
  c1a->cd(5+1);
  TLegend * leg1 = new TLegend(0.0,0.00,1.0,1.0,NULL,"brNDC");
  //      easyLeg(leg1,"Unfolded",0.06);
  easyLeg(leg1,"MC closure test",0.12);
  leg1->AddEntry(hMcTruth_M[5][0][0], "Truth","l");
  leg1->AddEntry(hMcReco_M[5][0][0], "Raw (Reco)","l");
  for (int in = 0; in < nUnfIterations ; in++)  {
    if ( in == 0 ) leg1->AddEntry(hMcUnf_M[5][0][0][in], Form("%dst iter.",in+1));
    else if ( in == 1 ) leg1->AddEntry(hMcUnf_M[5][0][0][in], Form("%dnd iter.",in+1));
    else if ( in == 2 ) leg1->AddEntry(hMcUnf_M[5][0][0][in], Form("%drd iter.",in+1));
    else  leg1->AddEntry(hMcUnf_M[5][0][0][in], Form("%dth iter.",in+1));
  }
  leg1->Draw();
  c1b->cd(5+1);
  leg1->Draw();
  c1c->cd(5+1);
  leg1->Draw();
  c1d->cd(5+1);
  leg1->Draw();
  //  int underflow = 2
  
  for ( int ipt = 0 ; ipt < nPtBins-1; ipt++)  {

    if (ipt < 5) c1a->cd(ipt%5 + 1);
    else if (ipt >= 5 && ipt < 10) c1b->cd(ipt%5 + 1);
    else if (ipt >= 10 && ipt < 15) c1c->cd(ipt%5 + 1);
    else if (ipt >= 15 && ipt < 20) c1d->cd(ipt%5 + 1);

    for (int iter = 0; iter < nUnfIterations ; iter++)  { 
      hMcTruth_M[5][0][ipt]->SetXTitle("m^{2}/#it{p}_{T}^{2}");
      hMcTruth_M[5][0][ipt]->SetXTitle("m/p_{T}");
      //      if ( hMcTruth_M[5][0][ipt]->Integral()>0) cleverRangeLog(hMcTruth_M[5][0][ipt],100,0.000001);
      hMcTruth_M[5][0][ipt]->SetYTitle("Entries");
      handsomeTH1(hMcTruth_M[5][0][ipt],1);
      handsomeTH1(hMcReco_M[5][0][ipt],1);
      handsomeTH1(hMcUnf_M[5][0][ipt][iter],color[iter]);      
      hMcTruth_M[5][0][ipt]->SetLineStyle(6);
      hMcUnf_M[5][0][ipt][iter]->SetMarkerStyle(mstyle[iter]);

      if ( iter == 0 )  {
	hMcTruth_M[5][0][ipt]->Draw("hist");
	hMcReco_M[5][0][ipt]->Draw("hist same");
      }
      hMcUnf_M[5][0][ipt][iter]->Draw("same e");
      gPad->SetLogy();
      //      gPad->SetLogx();
      
    }
    if ( ipt%5 == 0 )    {
      drawCentrality(kSample, 0, 0.25,0.8,1,24);
      ATLASLabel(0.25,0.92,"Internal",0.085,0.28);
    }
    // fix drawBins+
    //    drawBin(ptBins,ipt+1,"GeV",0.1 + (3* (ipt==0)), 0.3,1,18);
    drawBin(ptBins,ipt+2,"GeV",0.2, 0.3,1,18);
    //    drawBin(ptBins,ipt+1,"GeV",0.1, 0.3,1,18);

    if (ipt < 5) c1a->cd(ipt%5 + 2 + 5);
    else if (ipt >= 5 && ipt < 10) c1b->cd(ipt%5 + 2 + 5);
    else if (ipt >= 10 && ipt < 15) c1c->cd(ipt%5 + 2 + 5);
    else if (ipt >= 15 && ipt < 20) c1d->cd(ipt%5 + 2 + 5);
    for (int iter = 0; iter < nUnfIterations ; iter++)  {
      hMcRatio_M[5][0][ipt][iter] = (TH1D*)hMcUnf_M[5][0][ipt][iter]->Clone(Form("mcRatioSq_ix%d_iter%d",ipt,iter));
      hMcRatio_M[5][0][ipt][iter]->Divide(hMcTruth_M[5][0][ipt]);
      hMcRatio_M[5][0][ipt][iter]->SetAxisRange(.45,1.55,"Y");
      //      hMcRatio_M[5][0][ipt][iter]->SetAxisRange(0.07,0.3,"X");
      hMcRatio_M[5][0][ipt][iter]->SetYTitle("Ratio");
      hMcRatio_M[5][0][ipt][iter]->SetNdivisions(505,"X");
      hMcRatio_M[5][0][ipt][iter]->SetNdivisions(505,"Y");
      hMcRatio_M[5][0][ipt][iter]->SetTitleSize(.15,"X");
      hMcRatio_M[5][0][ipt][iter]->SetTitleSize(0.1,"Y");
      hMcRatio_M[5][0][ipt][iter]->SetTitleOffset(.5,"X");
      hMcRatio_M[5][0][ipt][iter]->SetTitleOffset(.7,"Y");
      fixedFontHist(hMcRatio_M[5][0][ipt][iter],2.5,2.5,20);
      //      hMcRatio_M[5][0][ipt][iter]->GetXaxis()->SetLabelOffset(0.02);
      if ( iter==0)  hMcRatio_M[5][0][ipt][iter]->Draw();
      else  hMcRatio_M[5][0][ipt][iter]->Draw("same");
    }
    if ( ipt == 0) 
      drawText("Unfolded/Truth", 0.22, 0.9,  1);

    //    gPad->SetLogx();

  }


  c1a->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:mcM_pt0-4"));
  //  c1a->SaveAs(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d_mcM_pt0-4_Reweighting.root",kSample,0,5));
  //    c1a->SaveAs(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d_mcM_pt0-4_NoReweighting.root",kSample,0,5));
  pdf_label = "";
  c1b->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:mcM_pt5-9"));
  //  c1b->SaveAs(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d_mcM_pt5-9_Reweighting.root",kSample,0,5));
  //  c1b->SaveAs(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d_mcM_pt5-9_NoReweighting.root",kSample,0,5));
  c1c->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:mcM_pt10-14"));
  //  c1c->SaveAs(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d_mcM_pt10-14_Reweighting.root",kSample,0,5));
  //c1c->SaveAs(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d_mcM_pt10-14_NoReweighting.root",kSample,0,5));
  c1d->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:mcM_pt15-19"));
  //  c1d->SaveAs(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d_mcM_pt14-19_Reweighting.root",kSample,0,5));
  //  c1d->SaveAs(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d_mcM_pt14-19_NoReweighting.root",kSample,0,5));

  TCanvas* c2 =  new TCanvas("c2","",1200,550);
  makeEfficiencyCanvas(c2,nPtBins-1+1, 0.02, 0.01, 0.2, 0.3, 0.01);

  TCanvas* c2a=  new TCanvas("c2a","",1200,550);
  TCanvas* c2b=  new TCanvas("c2b","",1200,550);
  TCanvas* c2c=  new TCanvas("c2c","",1200,550);
  TCanvas* c2d=  new TCanvas("c2d","",1200,550);
  TCanvas* c2e=  new TCanvas("c2e","",1200,550);

  makeEfficiencyCanvas(c2a,5+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c2b,5+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c2c,5+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c2d,5+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  cout << " making leg2 " << endl;  
  c2a->cd(5+1);
  TLegend *leg2 = new TLegend(0.002212389,0.001490066,0.9977876,1,NULL,"brNDC");
  easyLeg(leg2,"Data",0.12);
  leg2->AddEntry( hDataReco_M[5][0][0], "Raw","l");
    
  for (int iter = 0; iter < nUnfIterations ; iter++)  {
    if ( iter == 0 ) leg2->AddEntry(hDataUnf_M[5][0][0][iter], Form("%dst iter.",iter+1));
    else if ( iter == 1 ) leg2->AddEntry(hDataUnf_M[5][0][0][iter], Form("%dnd iter.",iter+1));
    else if ( iter == 2 ) leg2->AddEntry(hDataUnf_M[5][0][0][iter], Form("%drd iter.",iter+1));
    else  leg2->AddEntry(hDataUnf_M[5][0][0][iter], Form("%dth iter.",iter+1));
  }
  //    leg2->AddEntry( hMcTruth_M[5][0][0][0], "PYTHIA","l");
  leg2->Draw();

  c2b->cd(5+1);
  leg2->Draw();
  c2c->cd(5+1);
  leg2->Draw();
  c2d->cd(5+1);
  leg2->Draw();

  for ( int ipt = 0 ; ipt< nPtBins-1 ; ipt++)  {
    //    cout << " ipt " << ipt << endl;
    if (ipt < 5) c2a->cd(ipt%5 + 1);
    else if (ipt >= 5 && ipt < 10) c2b->cd(ipt%5 + 1);
    else if (ipt >= 10 && ipt < 15) c2c->cd(ipt%5 + 1);
    else if (ipt >= 15 && ipt < 20) c2d->cd(ipt%5 + 1);
    
    for (int iter = 0; iter < nUnfIterations ; iter++)  {
      hDataReco_M[5][0][ipt]->SetXTitle("m^{2}/#it{p}_{T}^{2}");
      //if ( hDataReco_M[5][0][ipt]->Integral()>0) cleverRangeLog(hDataReco_M[5][0][ipt],100,0.000001);
      hDataReco_M[5][0][ipt]->SetYTitle("Entries");
      handsomeTH1(hDataReco_M[5][0][ipt],1);
      handsomeTH1(hDataUnf_M[5][0][ipt][iter],color[iter]);
      hDataUnf_M[5][0][ipt][iter]->SetMarkerStyle(mstyle[iter]);
      if ( iter == 0 )  {
	hDataReco_M[5][0][ipt]->Draw("hist");
      }
      hDataUnf_M[5][0][ipt][iter]->Draw("same e");
	
      gPad->SetLogy();
      //      gPad->SetLogx();
	
    }
    //      hMcTruth_M[ipt][0]->Draw("same hist");
    //    drawBin(ptBins,ipt+1,"GeV",0.2 + (3* (ipt==0)), 0.7,1,18);
    drawBin(ptBins,ipt+2,"GeV",0.2, 0.3,1,18);
    // drawBin(ptBins,ipt+1,"GeV",0.1 , 0.3,1,18);
    if (ipt < 5) c2a->cd(ipt%5 + 2 + 5);
    else if (ipt >= 5 && ipt < 10) c2b->cd(ipt%5 + 2 + 5);
    else if (ipt >= 10 && ipt < 15) c2c->cd(ipt%5 + 2 + 5);
    else if (ipt >= 15 && ipt < 20) c2d->cd(ipt%5 + 2 + 5);
    bool drawFirst=true; 
    for (int iter = 0; iter < nUnfIterations ; iter++)  {
      //	if ( iter ==refId  ) continue;
      hDataRatio_M[5][0][ipt][iter] = (TH1D*)hDataUnf_M[5][0][ipt][iter]->Clone(Form("dataRatioSq_ix%d_in%d",ipt,iter));
      hDataRatio_M[5][0][ipt][iter]->Divide(hDataUnf_M[5][0][ipt][2]);

      hDataRatio_M[5][0][ipt][iter]->SetAxisRange(0.45,1.55,"Y");
      //      hDataRatio_M[5][0][ipt][iter]->SetAxisRange(0.07,0.3,"X");
      //      if ( optY==1)  hDataRatio_M[5][0][ipt][iter]->SetAxisRange(0.00,100,"X");
      hDataRatio_M[5][0][ipt][iter]->SetYTitle("Ratio");
      hDataRatio_M[5][0][ipt][iter]->SetNdivisions(505,"X");
      hDataRatio_M[5][0][ipt][iter]->SetNdivisions(505,"Y");
      fixedFontHist(hDataRatio_M[5][0][ipt][iter],2.5,2.5,20);
      if ( drawFirst)  { 	
	hDataRatio_M[5][0][ipt][iter]->Draw();

	drawFirst=false;  
      }
      else   hDataRatio_M[5][0][ipt][iter]->Draw("same");
      //	if ( optY == 2)  jumSun(0,1,0.3,1);
	
    }
      
      
    //    gPad->SetLogx();
      
  }
    

  c2a->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:dataM_pt0-4"));
  pdf_label = "";
  c2b->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:dataM_pt5-9"));
  c2c->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:dataM_pt10-14"));
  c2d->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:dataM_pt15-19"));  


  TCanvas* c3a=  new TCanvas("c3a","",1200,550);
  TCanvas* c3b=  new TCanvas("c3b","",1200,550);
  TCanvas* c3c=  new TCanvas("c3c","",1200,550);
  TCanvas* c3d=  new TCanvas("c3d","",1200,550);
  TCanvas* c3e=  new TCanvas("c3e","",1200,550);

  makeEfficiencyCanvas(c3a,6+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c3b,6+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c3c,6+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c3d,6+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c3e,6+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  
  //  makeEfficiencyCanvas(c3e,nPtBins-1+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  // place legend in final pannel
  c3a->cd(6+1);
  TLegend * leg3 = new TLegend(0.0,0.00,1.0,1.0,NULL,"brNDC");
  //      easyLeg(leg3,"Unfolded",0.06);
  easyLeg(leg3,"MC closure test",0.12);
  leg3->AddEntry(hMcTruth_M[5][0][0], "Truth","l");
  leg3->AddEntry(hMcReco_M[5][0][0], "Raw (Reco)","l");
  for (int in = 0; in < nUnfIterations ; in++)  {
    if ( in == 0 ) leg3->AddEntry(hMcUnf_M[5][0][0][in], Form("%dst iter.",in+1));
    else if ( in == 1 ) leg3->AddEntry(hMcUnf_M[5][0][0][in], Form("%dnd iter.",in+1));
    else if ( in == 2 ) leg3->AddEntry(hMcUnf_M[5][0][0][in], Form("%drd iter.",in+1));
    else  leg3->AddEntry(hMcUnf_M[5][0][0][in], Form("%dth iter.",in+1));
  }
  leg3->Draw();
  c3b->cd(6+1);
  leg3->Draw();
  c3c->cd(6+1);
  leg3->Draw();
  c3d->cd(6+1);
  leg3->Draw();
  c3e->cd(6+1);
  leg3->Draw();
  
  for ( int im = 0 ; im < nMBins-1 ; im++)  {

    if (im < 6) c3a->cd(im%6 + 1);
    else if (im >= 6 && im < 12) c3b->cd(im%6 + 1);
    else if (im >= 12 && im < 18) c3c->cd(im%6 + 1);
    else if (im >= 18 && im < 24) c3d->cd(im%6 + 1);
    else if (im >= 24 && im < 30) c3e->cd(im%6 + 1);

    for (int iter = 0; iter < nUnfIterations ; iter++)  { 
      hMcTruth_Pt[5][0][im]->SetXTitle("#it{p}_{T}");
      //      if ( hMcTruth_Pt[5][0][im]->Integral()>0) cleverRangeLog(hMcTruth_Pt[5][0][im],100,0.000001);
      hMcTruth_Pt[5][0][im]->SetYTitle("Entries");
      handsomeTH1(hMcTruth_Pt[5][0][im],1);
      handsomeTH1(hMcReco_Pt[5][0][im],1);
      handsomeTH1(hMcUnf_Pt[5][0][im][iter],color[iter]);      
      hMcTruth_Pt[5][0][im]->SetLineStyle(6);
      hMcUnf_Pt[5][0][im][iter]->SetMarkerStyle(mstyle[iter]);

      if ( iter == 0 )  {
	hMcTruth_Pt[5][0][im]->Draw("hist");
	hMcReco_Pt[5][0][im]->Draw("hist same");
      }
      hMcUnf_Pt[5][0][im][iter]->Draw("same e");
      gPad->SetLogy();
      //      gPad->SetLogx();
      
    }
    if ( im == 0 )    {
      drawCentrality(kSample, 0, 0.25,0.8,1,24);
      ATLASLabel(0.25,0.92,"Internal",0.085,0.28);
    }
    //    drawBin(mBins,im+1,"GeV",0.2 + (3* (im==0)), 0.7,1,18);
    drawBin(mBins,im+2,"GeV",0.2, 0.7,1,18);
    // drawBin(mBins,im+1,"GeV",0.2 + (0.05* (im==0)), 0.7,1,18);
    if (im < 6) c3a->cd(im%6 + 2 + 6);
    else if (im >= 6 && im < 12) c3b->cd(im%6 + 2 + 6);
    else if (im >= 12 && im < 18) c3c->cd(im%6 + 2 + 6);
    else if (im >= 18 && im < 24) c3d->cd(im%6 + 2 + 6);
    else if (im >= 24 && im < 30) c3e->cd(im%6 + 2 + 6);
    for (int iter = 0; iter < nUnfIterations ; iter++)  {
      hMcRatio_Pt[5][0][im][iter] = (TH1D*)hMcUnf_Pt[5][0][im][iter]->Clone(Form("mcRatioSq_ix%d_iter%d",im,iter));
      hMcRatio_Pt[5][0][im][iter]->Divide(hMcTruth_Pt[5][0][im]);
      hMcRatio_Pt[5][0][im][iter]->SetAxisRange(.5,1.5,"Y");
      //      hMcRatio_Pt[5][0][im][iter]->SetAxisRange(90,1000,"X");
      hMcRatio_Pt[5][0][im][iter]->SetYTitle("Ratio");
      hMcRatio_Pt[5][0][im][iter]->SetNdivisions(505,"X");
      hMcRatio_Pt[5][0][im][iter]->SetNdivisions(505,"Y");
      hMcRatio_Pt[5][0][im][iter]->SetTitleSize(.15,"X");
      hMcRatio_Pt[5][0][im][iter]->SetTitleSize(0.1,"Y");
      hMcRatio_Pt[5][0][im][iter]->SetTitleOffset(.5,"X");
      hMcRatio_Pt[5][0][im][iter]->SetTitleOffset(.7,"Y");
      fixedFontHist(hMcRatio_Pt[5][0][im][iter],2.5,2.5,20);
  
      //      hMcRatio_Pt[5][0][im][iter]->GetXaxis()->SetLabelOffset(0.02);
      if ( iter==0)  hMcRatio_Pt[5][0][im][iter]->Draw();
      else  hMcRatio_Pt[5][0][im][iter]->Draw("same");
    }
    if ( im == 0) 
      drawText("Unfolded/Truth", 0.22, 0.9,  1);

    //    gPad->SetLogx();

  }

  c3a->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:mcPt_m0-5"));
  pdf_label = "";
  c3b->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:mcPt_m6-11"));
  c3c->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:mcPt_m12-17"));
  c3d->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:mcPt_m18-23"));
  c3e->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:mcPt_m24-29"));  
  

  TCanvas* c4a=  new TCanvas("c4a","",1200,550);
  TCanvas* c4b=  new TCanvas("c4b","",1200,550);
  TCanvas* c4c=  new TCanvas("c4c","",1200,550);
  TCanvas* c4d=  new TCanvas("c4d","",1200,550);
  TCanvas* c4e=  new TCanvas("c4e","",1200,550);

  makeEfficiencyCanvas(c4a,6+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c4b,6+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c4c,6+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c4d,6+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  makeEfficiencyCanvas(c4e,6+1, 0.02, 0.01, 0.2, 0.3, 0.01);
  cout << " making leg2 " << endl;  
  c4a->cd(6+1);
  TLegend *leg4 = new TLegend(0.002212389,0.001490066,0.9977876,1,NULL,"brNDC");
  easyLeg(leg4,"Data",0.12);
  leg4->AddEntry( hDataReco_Pt[5][0][0], "Raw","l");
    
  for (int iter = 0; iter < nUnfIterations ; iter++)  {
    if ( iter == 0 ) leg4->AddEntry(hDataUnf_Pt[5][0][0][iter], Form("%dst iter.",iter+1));
    else if ( iter == 1 ) leg4->AddEntry(hDataUnf_Pt[5][0][0][iter], Form("%dnd iter.",iter+1));
    else if ( iter == 2 ) leg4->AddEntry(hDataUnf_Pt[5][0][0][iter], Form("%drd iter.",iter+1));
    else  leg4->AddEntry(hDataUnf_Pt[5][0][0][iter], Form("%dth iter.",iter+1));
  }
  //    leg4->AddEntry( hMcTruth_Pt[5][0][0][0], "PYTHIA","l");
  leg4->Draw();

  c4b->cd(6+1);
  leg4->Draw();
  c4c->cd(6+1);
  leg4->Draw();
  c4d->cd(6+1);
  leg4->Draw();
  c4e->cd(6+1);
  leg4->Draw();

  for ( int im = 0 ; im< nMBins-1 ; im++)  {
    //    cout << " im " << im << endl;
    if (im < 6) c4a->cd(im%6 + 1);
    else if (im >= 6 && im < 12) c4b->cd(im%6 + 1);
    else if (im >= 12 && im < 18) c4c->cd(im%6 + 1);
    else if (im >= 18 && im < 24) c4d->cd(im%6 + 1);
    else if (im >= 24 && im < 30) c4e->cd(im%6 + 1);
    
    for (int iter = 0; iter < nUnfIterations ; iter++)  {
      hDataReco_Pt[5][0][im]->SetXTitle("#it{p}_{T}");
      //if ( hDataReco_Pt[5][0][im]->Integral()>0) cleverRangeLog(hDataReco_Pt[5][0][im],100,0.000001);
      hDataReco_Pt[5][0][im]->SetYTitle("Entries");
      handsomeTH1(hDataReco_Pt[5][0][im],1);
      handsomeTH1(hDataUnf_Pt[5][0][im][iter],color[iter]);
      hDataUnf_Pt[5][0][im][iter]->SetMarkerStyle(mstyle[iter]);
      if ( iter == 0 )  {
	hDataReco_Pt[5][0][im]->Draw("hist");
      }
      hDataUnf_Pt[5][0][im][iter]->Draw("same e");
	
      gPad->SetLogy();
      //      gPad->SetLogx();
	
    }
    //      hMcTruth_Pt[im][0]->Draw("same hist");
    drawBin(mBins,im+2,"GeV",0.2, 0.7,1,18);
    //    drawBin(mBins,im+1,"GeV",0.2 + (3* (im==0)), 0.7,1,18);
    //drawBin(mBins,im+1,"GeV",0.2 + (0.05* (im==0)), 0.7,1,18);
    if (im < 6) c4a->cd(im%6 + 2 + 6);
    else if (im >= 6 && im < 12) c4b->cd(im%6 + 2 + 6);
    else if (im >= 12 && im < 18) c4c->cd(im%6 + 2 + 6);
    else if (im >= 18 && im < 24) c4d->cd(im%6 + 2 + 6);
    else if (im >= 24 && im < 30) c4e->cd(im%6 + 2 + 6);
    bool drawFirst=true; 
    for (int iter = 0; iter < nUnfIterations ; iter++)  {
      //	if ( iter ==refId  ) continue;
      hDataRatio_Pt[5][0][im][iter] = (TH1D*)hDataUnf_Pt[5][0][im][iter]->Clone(Form("dataRatioSq_ix%d_in%d",im,iter));
      hDataRatio_Pt[5][0][im][iter]->Divide(hDataUnf_Pt[5][0][im][2]);

      hDataRatio_Pt[5][0][im][iter]->SetAxisRange(0.5,1.5,"Y");
      //      hDataRatio_Pt[5][0][im][iter]->SetAxisRange(90,1000,"X");
      //      if ( optY==1)  hDataRatio_Pt[5][0][im][iter]->SetAxisRange(0.00,100,"X");
      hDataRatio_Pt[5][0][im][iter]->SetYTitle("Ratio");
      hDataRatio_Pt[5][0][im][iter]->SetNdivisions(505,"X");
      hDataRatio_Pt[5][0][im][iter]->SetNdivisions(505,"Y");
      fixedFontHist(hDataRatio_Pt[5][0][im][iter],2.5,2.5,20);
      if ( drawFirst)  { 	
	hDataRatio_Pt[5][0][im][iter]->Draw();

	drawFirst=false;  
      }
      else   hDataRatio_Pt[5][0][im][iter]->Draw("same");
      //	if ( optY == 2)  jumSun(0,1,0.3,1);
	
    }
      
      
    //    gPad->SetLogx();
      
  }
    

  c4a->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:dataPt_m0-5"));
  pdf_label = "";
  c4b->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:dataPt_m6-11"));
  c4c->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:dataPt_m12-17"));
  c4d->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:dataPt_m18-23"));
  pdf_label = ")";
  c4e->Print(Form("Unfolding/Closure/Closure_sample%d_icent%d_etaBin%d.pdf%s",kSample,0,5, pdf_label.c_str()), Form("Title:dataPt_m24-29"));  
  


return 0;
}
