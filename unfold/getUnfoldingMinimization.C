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
  gStyle->SetPalette(kRainBow);
  int nCentBins;
  if (kSample == kPP) nCentBins = 1;
  else nCentBins = 7; 
  int nEtaBins = 6;  
  int nUnfIterations = 9;
  


  TFile *f1 = TFile::Open(Form("rawSpectra/RawSpectraData_sample%d.root",kSample)); 
  //  TFile *f2 = TFile::Open(Form("unfSpectra/rooUnfoldResponse_ReweightingMin0P4Max1P5_sample%d.root",kSample));

  TFile *f2 = TFile::Open(Form("unfSpectra/rooUnfoldResponse_Reweighting_sample%d_ptCut100Events.root",kSample));
  //  TFile *f2 = TFile::Open(Form("unfSpectra/rooUnfoldResponse_NoReweighting_sample%d.root",kSample));
  //  TFile *f2 = TFile::Open(Form("unfSpectra/rooUnfoldResponse_ReweightMin0P2Max2P0_sample%d.root",kSample));
  //TFile *f2 = TFile::Open(Form("unfSpectra/rooUnfoldResponse_ReweightMin0P4Max1P5_sample%d.root",kSample));
 //TFile *f2 = TFile::Open(Form("unfSpectra/rooUnfoldResponse_ReweightNoLimit_sample%d.root",kSample));
    //    TFile *f2 = TFile::Open(Form("unfSpectra/rooUnfoldResponse_ReweightingMin0P4Max1P5_sample%d_ptCut10Events.root",kSample)); 
  
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

  long long entries;
  TString fname;
  
  RooUnfoldResponse* rooUnfoldRes[6][7]; // eta, cent
  TH2D* hDataReco[6][7]; 
  TH2D* hDataUnf[6][7][20]; // unfolding iter
  std::string pdf_label;
  
   for (int etaBin = 0; etaBin < nEtaBins; ++etaBin) {
     for (int centBin = 0; centBin < nCentBins; ++centBin) {
       hDataReco[etaBin][centBin] = (TH2D*) f1->Get(Form("hCrossSecData_m_pt_eta%dcent%d",etaBin,centBin));
       rooUnfoldRes[etaBin][centBin] = (RooUnfoldResponse*) f2->Get(Form("responseMatrix_eta%d_cent%d",etaBin,centBin));

    }   
  }
  for (int etaBin = 5; etaBin < nEtaBins; ++etaBin) {
    for ( int centBin=0 ; centBin < nCentBins; ++centBin) {
      if ( (kSample == kPP) && ( centBin != 0 ) )      continue;
      // loop through unfolding iterations
      for ( int it = 0 ; it < nUnfIterations; it++) {  

	// Data
	RooUnfoldBayes unfoldData (rooUnfoldRes[etaBin][centBin], hDataReco[etaBin][centBin], it);    
	hDataUnf[etaBin][centBin][it]= (TH2D*)unfoldData.Hreco();
	hDataUnf[etaBin][centBin][it]->SetName( Form("hDataUnf_eta%d_cent%d_iter%d",etaBin,centBin,it));
      }
    }
  }
  cout << " done unfolding " << endl;

  TH1D *hDataReco_Pt[6][7][30];
  TH1D *hDataUnf_Pt[6][7][30][20];
  TH1D *hDataRatio_Pt[6][7][30][20];
  TH1D *hDataReco_M[6][7][20];
  TH1D *hDataUnf_M[6][7][20][20];
  TH1D *hDataRatio_M[6][7][20][20];
    
  for (int etaBin = 5; etaBin < nEtaBins; ++etaBin) {
    for (int centBin = 0; centBin < nCentBins; ++centBin) {
      for (int mBin = 0; mBin < nMBins; ++mBin) {

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

	hDataReco_M[etaBin][centBin][ptBin] = (TH1D*)hDataReco[etaBin][centBin]->ProjectionY(Form("hdataRaw_M_eta%d_cent%d_pt%d",etaBin,centBin,ptBin),ptBin+2,ptBin+2);      
	for ( int it = 0 ; it< nUnfIterations ; it++) {
	  hDataUnf_M[etaBin][centBin][ptBin][it] = (TH1D*)hDataUnf[etaBin][centBin][it]->ProjectionY(Form("hdataUnf_M_eta%d_icent%d_pt%d_iter%d",etaBin,centBin,ptBin,it),ptBin+2,ptBin+2);
	  TH1ScaleByWidth(hDataUnf_M[etaBin][centBin][ptBin][it]);
	  TH1ScaleByRapidity(hDataUnf_M[etaBin][centBin][ptBin][it], etaBin);	
	}
	TH1ScaleByWidth(hDataReco_M[etaBin][centBin][ptBin]);
	TH1ScaleByRapidity(hDataReco_M[etaBin][centBin][ptBin], etaBin);
      }
    }
  }



  pdf_label = "(";
  
  TH1D* hStat[30];
  TH1D* hDevi[30];
  for ( int ptBin = 0 ; ptBin < nPtBins ; ++ptBin)  {
    hDevi[ptBin] = new TH1D(Form("hdevi_pt%d",ptBin),";Number of iterations;Relative uncertainty",10,0.5,10.5);
    hStat[ptBin] = (TH1D*)hDevi[ptBin]->Clone(Form("hstat_ptBin%d",ptBin));

  }
  
  TCanvas* c1[nMBins-1];
  Double_t sumSigma[40];
  Double_t sumDev[40];
  Double_t sumStat[40];
  int sigmaCount[40];
  for (int i = 0; i < 40; ++i) {
    sumSigma[i] = 0;
    sumStat[i] = 0;
    sumDev[i] = 0;
    sigmaCount[i] = 0;
  }
  for ( int im = 0 ; im < nMBins-1 ; im++)  {
    char cName[256];
    sprintf(cName,"c_%d",im);
    c1[im] =  new TCanvas(cName,cName,1200,600);
    c1[im]->Divide(6,3);
    for ( int ipt = 0; ipt< nPtBins-2 ; ipt++)  {
      std::vector<Double_t> val;
      std::vector<Double_t> err;

      c1[im]->cd(ipt + 1);
      for (int in = 0; in < nUnfIterations ; in++) {
	// need at least one underflow bin in mass and pt
	//	val.push_back(hDataUnf[5][0][in]->GetBinContent( ipt+2,im+2));
	val.push_back(hDataUnf[5][0][in]->GetBinContent( ipt+3,im+2));
	//	err.push_back(hDataUnf[5][0][in]->GetBinError( ipt+2,im+2));
	err.push_back(hDataUnf[5][0][in]->GetBinError( ipt+3,im+2));
      }
      for (int in = 2; in < nUnfIterations ; in++) {
	Double_t relStatSigma = err.at(in) / val.at(in);
	Double_t relMigration = fabs(val.at(in) - val.at(in-1))/val.at(in);
	//	cout << " im " << im << " ipt " << ipt << " iter " << in << " val " << val[in] << " err " << err[in] << " relStatSigma " << relStatSigma << " relMigration " << relMigration << endl; 
	hStat[ipt]->SetBinContent( in-1, relStatSigma) ;
	hDevi[ipt]->SetBinContent( in-1, relMigration);
	hStat[ipt]->SetBinError( in-1, 0.0000001);
	hDevi[ipt]->SetBinError( in-1, 0.0000001);
	Double_t sigma = TMath::Sqrt((relStatSigma*relStatSigma) + (relMigration*relMigration)); 
	if (sigma == sigma) {
	  sumSigma[in-1] += sigma;
	  sumDev[in-1] += relMigration;
	  sumStat[in-1] += relStatSigma;
	  ++sigmaCount[in-1];
	}
	//	cout << " total sigma " << sumSigma[in-1] << " sigma " << TMath::Sqrt((relStatSigma*relStatSigma) + (relMigration*relMigration)) << endl;
      }

      handsomeTH1(hStat[ipt],4,1,26);
      handsomeTH1(hDevi[ipt],2,1,32);
      fixedFontHist(hStat[ipt],2.2,3.6,9);
      fixedFontHist(hDevi[ipt],2.2,3.6,9);
      hStat[ipt]->SetAxisRange(0.,0.03,"Y");
      //    int fScale = 5;
      //    if ( ipt > lowPtBin + 2)  hStat[ipt]->Scale(1./fScale);
      /*
      if ( icent <=3)  hStat[ipt]->SetAxisRange(-0.01,0.03,"Y");
      else if ( icent ==4)  hStat[ipt]->SetAxisRange(-0.50,1.21,"Y");
      else if ( icent ==5)  hStat[ipt]->SetAxisRange(-0.50,1.21,"Y");
      else if ( icent ==6)  hStat[ipt]->SetAxisRange(-1.0,5.1,"Y");
      */
      hDevi[ipt]->Draw();
      hStat[ipt]->Draw("same");



      if ( ipt == 1 ) {
	drawCentrality(kSample, 0, 0.33,0.60,1,18);
	ATLASLabel(0.33,0.68,"Internal",0.08,0.26);
      }
      else if ( ipt == 0 ) {
	TLegend *leg1 = new TLegend(0.3830109,0.5570667,0.9199772,0.74624,NULL,"brNDC");
	//	easyLeg(leg1,Form("%.2f < m^{2}/#it{p}_{T}^{2} < %.2f",(float)mBins[im], (float)mBins[im+1]),0.07);
	easyLeg(leg1);
	leg1->AddEntry(hDevi[ipt], "|y_{N} - y_{N-1}| / y_{N}","p");
	leg1->AddEntry(hStat[ipt], "Stat. uncertainty","p");
	leg1->Draw();
	drawBinMpt(mBins,im+2,"",0.33 + (0.02* (ipt==0)), 0.78,1,11);
      }
      drawBinPt(ptBins,ipt+3,"GeV",0.33 + (0.02* (ipt==0)), 0.86,1,11);
      //    if ( ipt > lowPtBin + 2) 
      //      drawText(Form("Stat. Unc scaled by 1/%d",fScale), 0.16 + (0.05* (ipt==lowPtBin)), 0.3,1,16);
      //      jumSun(0,0,49,0);
    }
    char name[256];
    sprintf(name,"Title:Iteration_optimization_m%d",im);
    if (im == 0) pdf_label = "(";
    else pdf_label = "";
    c1[im]->Print(Form("Unfolding/Closure/IterationMin_sample%d_icent%d_etaBin%dReweight_ptCut100Events.pdf%s",kSample,0,5, pdf_label.c_str()), name);
  }

  TCanvas *c3 = new TCanvas("c3","c3",500,500);
  TH1D *hSigmaSum = new TH1D("sigmaSum", ";Number of Iterations;#sigma",30,0,29);
  TH1D *hDevSum = new TH1D("devSum", ";Number of Iterations;#sigma",30,0,29);
  TH1D *hStatSum = new TH1D("statSum", ";Number of Iterations;#sigma",30,0,29);
  for (int i = 0; i < nUnfIterations; ++i) {
    if (sigmaCount[i] == 0)     hSigmaSum->SetBinContent(i+1, 0);
    else    hSigmaSum->SetBinContent(i+1, sumSigma[i]/sigmaCount[i]);
    if (sigmaCount[i] == 0)     hDevSum->SetBinContent(i+1, 0);
    else    hDevSum->SetBinContent(i+1, sumDev[i]/sigmaCount[i]);
    if (sigmaCount[i] == 0)     hStatSum->SetBinContent(i+1, 0);
    else    hStatSum->SetBinContent(i+1, sumStat[i]/sigmaCount[i]);
    cout << " iteration " << i << " total sigma " << sumSigma[i]/sigmaCount[i] << endl;  
  }
  handsomeTH1(hSigmaSum,kBlack,1,20);
  handsomeTH1(hDevSum,kRed+1,1,kFullStar);
  handsomeTH1(hStatSum,kGreen+1,1,kFullDiamond);
  hSigmaSum->Draw("P");
  hDevSum->Draw("SameP");
  hStatSum->Draw("SameP");
  drawCentrality(kSample,-1,0.65,0.75,1,24);
  ATLASLabel(0.25, 0.80, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
  TLegend *leg10 = new TLegend(0.3830109,0.5570667,0.9199772,0.74624,NULL,"brNDC");
  //	easyLeg(leg1,Form("%.2f < m^{2}/#it{p}_{T}^{2} < %.2f",(float)mBins[im], (float)mBins[im+1]),0.07);
  easyLeg(leg10);
  leg10->AddEntry(hDevSum, "|y_{N} - y_{N-1}| / y_{N}","p");
  leg10->AddEntry(hStatSum, "Stat. uncertainty","p");
  leg10->AddEntry(hSigmaSum, "#sigma_{Total}","p");
  leg10->Draw();
  pdf_label = "";
  char name[256];
  sprintf(name,"Title:SummedSigma");
  c3->Print(Form("Unfolding/Closure/IterationMin_sample%d_icent%d_etaBin%dReweight_ptCut100Events.pdf%s",kSample,0,5, pdf_label.c_str()), name);
  c3->SaveAs(Form("Unfolding/Closure/IterationMin_sample%d_icent%d_etaBin%dReweight_ptCut100Events.root",kSample,0,5));

  TH1D* hStatM[30];
  TH1D* hDeviM[30];

  for ( int mBin = 0 ; mBin < nMBins ; ++mBin)  {
    hDeviM[mBin] = new TH1D(Form("hdevi_m%d",mBin),";Number of iterations;Relative uncertainty",10,0.5,10.5);
    hStatM[mBin] = (TH1D*)hDeviM[mBin]->Clone(Form("hstat_mBin%d",mBin));

  }
  
  Double_t sumSigmaM[40];
  Double_t sumDevM[40];
  Double_t sumStatM[40];
  int sigmaCountM[40];
  TCanvas* c4[nPtBins-2];  
  for (int i = 0; i < 40; ++i) {
    sumSigmaM[i] = 0;
    sigmaCountM[i] = 0;
  }
  for ( int ipt = 0 ; ipt < nPtBins-2; ipt++)  {
    char cName[256];
    sprintf(cName,"c4_%d",ipt);
    c4[ipt] =  new TCanvas(cName,cName,1200,600);
    c4[ipt]->Divide(6,5);
    for ( int im = 0; im < nMBins-1 ; im++)  {

      std::vector<Double_t> val;
      std::vector<Double_t> err;

      c4[ipt]->cd(ipt + 1);
      for (int in = 0; in < nUnfIterations ; in++) {
	// need at least one underflow bin in mass and pt
	//	val.push_back(hDataUnf[5][0][in]->GetBinContent( ipt+2,im+2));
	val.push_back(hDataUnf[5][0][in]->GetBinContent( ipt+3,im+2));
	//	err.push_back(hDataUnf[5][0][in]->GetBinError( ipt+2,im+2));
	err.push_back(hDataUnf[5][0][in]->GetBinError( ipt+3,im+2));
      }
      for (int in = 2; in < nUnfIterations ; in++) {
	Double_t relStatSigma = err.at(in) / val.at(in);
	Double_t relMigration = fabs(val.at(in) - val.at(in-1))/val.at(in);
	//	cout << " im " << im << " ipt " << ipt << " iter " << in << " val " << val[in] << " err " << err[in] << " relStatSigma " << relStatSigma << " relMigration " << relMigration << endl; 
	hStatM[im]->SetBinContent( in-1, relStatSigma) ;
	hDeviM[im]->SetBinContent( in-1, relMigration);
	hStatM[im]->SetBinError( in-1, 0.0000001);
	hDeviM[im]->SetBinError( in-1, 0.0000001);
	Double_t sigma = TMath::Sqrt((relStatSigma*relStatSigma) + (relMigration*relMigration)); 
	if (sigma == sigma) {
	  sumSigmaM[in-1] += sigma;
	  sumDevM[in-1] += relMigration;
	  sumStatM[in-1] += relStatSigma;
	  ++sigmaCountM[in-1];
	}
	//	cout << " total sigma " << sumSigmaM[in-1] << " sigma " << TMath::Sqrt((relStatSigma*relStatSigma) + (relMigration*relMigration)) << endl;
      }

      handsomeTH1(hStatM[im],4,1,26);
      handsomeTH1(hDeviM[im],2,1,32);
      fixedFontHist(hStatM[im],2.2,3.6,9);
      fixedFontHist(hDeviM[im],2.2,3.6,9);
      //      hStatM[im]->SetAxisRange(0.,0.03,"Y");
      //    int fScale = 5;
      //    if ( ipt > lowPtBin + 2)  hStatM[im]->Scale(1./fScale);
      /*
      if ( icent <=3)  hStatM[im]->SetAxisRange(-0.01,0.03,"Y");
      else if ( icent ==4)  hStatM[im]->SetAxisRange(-0.50,1.21,"Y");
      else if ( icent ==5)  hStatM[im]->SetAxisRange(-0.50,1.21,"Y");
      else if ( icent ==6)  hStatM[im]->SetAxisRange(-1.0,5.1,"Y");
      */
      hDeviM[im]->Draw();
      hStatM[im]->Draw("same");



      if ( ipt == 1 ) {
	drawCentrality(kSample, 0, 0.33,0.60,1,18);
	ATLASLabel(0.33,0.68,"Internal",0.08,0.26);
      }
      else if ( ipt == 0 ) {
	TLegend *leg1 = new TLegend(0.3830109,0.5570667,0.9199772,0.74624,NULL,"brNDC");
	//	easyLeg(leg1,Form("%.2f < m^{2}/#it{p}_{T}^{2} < %.2f",(float)mBins[im], (float)mBins[im+1]),0.07);
	easyLeg(leg1);
	leg1->AddEntry(hDeviM[im], "|y_{N} - y_{N-1}| / y_{N}","p");
	leg1->AddEntry(hStatM[im], "Stat. uncertainty","p");
	leg1->Draw();
	drawBinPt(ptBins,ipt+3,"GeV",0.33 + (0.02* (ipt==0)), 0.78,1,11);
      }
	drawBinMpt(mBins,im+2,"",0.33 + (0.02* (ipt==0)), 0.86,1,11);
      //    if ( ipt > lowPtBin + 2) 
      //      drawText(Form("Stat. Unc scaled by 1/%d",fScale), 0.16 + (0.05* (ipt==lowPtBin)), 0.3,1,16);
      //      jumSun(0,0,49,0);
    }
    char name[256];
    sprintf(name,"Title:Iteration_optimization_pt%d",ipt);
    pdf_label = "";
    c4[ipt]->Print(Form("Unfolding/Closure/IterationMin_sample%d_icent%d_etaBin%dReweight_ptCut100Events.pdf%s",kSample,0,5, pdf_label.c_str()), name);
  }

  TCanvas *c5 = new TCanvas("c5","c5",500,500);
  TH1D *hSigmaSumM = new TH1D("sigmaSumM", ";Number of Iterations;#sigma_{Total}",30,0,29);
  TH1D *hDevSumM = new TH1D("devSumM", ";Number of Iterations;#sigma",30,0,29);
  TH1D *hStatSumM = new TH1D("statSumM", ";Number of Iterations;#sigma",30,0,29);
  for (int i = 0; i < nUnfIterations; ++i) {
    if (sigmaCountM[i] == 0)     hSigmaSumM->SetBinContent(i+1, 0);
    else    hSigmaSumM->SetBinContent(i+1, sumSigmaM[i]/sigmaCountM[i]);
    if (sigmaCountM[i] == 0)     hDevSumM->SetBinContent(i+1, 0);
    else    hDevSumM->SetBinContent(i+1, sumDevM[i]/sigmaCountM[i]);
    if (sigmaCountM[i] == 0)     hStatSumM->SetBinContent(i+1, 0);
    else    hStatSumM->SetBinContent(i+1, sumStatM[i]/sigmaCountM[i]);
    //    cout << " iteration " << i << " total sigma " << sumSigmaM[i]/sigmaCountM[i] << endl;  
  }
  handsomeTH1(hSigmaSumM,1,1,20);
  handsomeTH1(hDevSumM,kRed+1,1,kFullStar);
  handsomeTH1(hStatSumM,kGreen+1,1,kFullDiamond);
  hSigmaSumM->Draw("P");
  hDevSumM->Draw("SameP");
  hStatSumM->Draw("SameP");
  drawCentrality(kSample,-1,0.65,0.75,1,24);
  ATLASLabel(0.25, 0.80, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
  TLegend *leg11 = new TLegend(0.3830109,0.5570667,0.9199772,0.74624,NULL,"brNDC");
  //	easyLeg(leg1,Form("%.2f < m^{2}/#it{p}_{T}^{2} < %.2f",(float)mBins[im], (float)mBins[im+1]),0.07);
  easyLeg(leg11);
  leg11->AddEntry(hDevSumM, "|y_{N} - y_{N-1}| / y_{N}","p");
  leg11->AddEntry(hStatSumM, "Stat. uncertainty","p");
  leg11->AddEntry(hSigmaSumM, "#sigma_{Total}","p");
  leg11->Draw();
  pdf_label = ")";
  sprintf(name,"Title:SummedSigma");
  c5->Print(Form("Unfolding/Closure/IterationMin_sample%d_icent%d_etaBin%dReweight_ptCut100Events.pdf%s",kSample,0,5, pdf_label.c_str()), name);
  c5->SaveAs(Form("Unfolding/Closure/IterationMin_sample%d_icent%d_etaBin%dReweight_ptCut100Events.root",kSample,0,5));
  
  return 0;
}
 
