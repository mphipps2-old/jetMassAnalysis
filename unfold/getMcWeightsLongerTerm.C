// This file derives the pt and mass spectra reweighting factors. This is needed b/c the MC spectra are significantly different from data. These are derived by dividing the pt and m/pt distributions for each pp and PbPb centrailty bin

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

#include <TRandom3.h>
#include "TH1D.h"
#include "TList.h"
#include "../getSdHists.C"
//#include "../ntupleDefinition.h"
//#include "../ntupleDefinition.h"
#include "../commonUtility.h"
#include "../jzWeight.h"
// this makes a class out of ntuple so all struct values can be accessed -- there were problems accessing the final memebers of the struct directly
//#include "../ntuples/tashka.h"
#include "../ntuples/tashka.C"
#endif

#include "unfoldingUtil.h"

#include "../JssUtils.h"
#include <TPaletteAxis.h>
#include "systematicsTool.h"

double statFrac = 1.0;
bool doUnfData = true ;

int lowPtBin = 1;  int highPtBin = 13;
//int nPtPannels = highPtBin-lowPtBin+1;

bool isTooSmall(TH2D* hEntries=0, int recoVarX=0, int recoVarY=0, int minEntries=10);
void getMCspectra(int kSample=kPP, int icent=0, int etaBin=0, TH2D* hmcRaw=0, TH2D* hmcRaw_noWeight=0, TH2D* hmcRaw_jetWeight=0, TH2D* hmcRaw_jzNorm=0, TH2D* hmcRaw_fcalWeight=0, TH2D* hmcTruth=0,/* TH2D* hmcRawGlobal=0, TH2D* hmcRawGlobalTruth=0,*/ TF1* ptScale=0, int nSys=-1);
void getDATAspectra(int kSample=kPP, int icent=0, int etaBin=0, TH2D* hdataRaw=0,/*TH2D* hdataRawGlobal=0,*/ int nSys=0);

TH1D* getVariedHist(TH1D* hin=0, double variation=0);
//bool isTooSmall(TH2D* hEntries=0, int recoVarX=0, int recoVarY=0, int minEntries=10);

float flucCut = 1.;
void removeFluc2(TH2* h=0);
// unfoldingUtil.h defines this binning
// kSample: kPP or kPbPb
// icent: if kSample == kPP, icent = 0; else icent >=0 and <= 6
// etaBin: 0 -> |y| < 0.3; 1 -> |y| >= 0.3 && |y| < 1.2; 2 -> |y| >= 1.2 |y| < 2.1
// default optX: 77 -- sets ptBinning 
// default optY: 772 -- sets mass binning and retrieves m2/pt2

// declare eta wide global histos here and pt and m2/pt2 containers here
TH2D* hmcRawGlobal;
TH2D* hmcRawGlobalTruth;
TH2D* hmcRawGlobalCorr;
TH2D* hmcRawGlobalTruthCorr;
TH2D* hdataRawGlobal;
TList *listmcGlobal;
TList *listmcGlobalTruth;
TList *listdataGlobal;

void getMcWeights(int kSample = kPP, int icent=0, int etaBin=3, int nSys=-1) { 
  cout << "kSample " << kSample << " icent " << icent << " etaBin " << etaBin << " nSys " << nSys << endl;

  TH1::SetDefaultSumw2();
  
  int nXbins;
  double xBin[30];
  getXbin(nXbins, xBin, 77);
  cout << " nXbins = " << nXbins << endl;
  cout << " xBin = " << xBin[0] << ",   " << xBin[1] << ",   " <<xBin[2] << ", ..." <<endl;


  int nYbins ;
  double yBin[200] ;

  getYbin(nYbins, yBin, 772);
  cout << " nYbins = " << nYbins << endl;
  cout << " yBin = " << yBin[0] << ",   " << yBin[1] << ",   " <<yBin[2] << ", ..." <<endl;
  
  TH2D* hTemp = new TH2D("hptTemp","", nXbins, xBin, nYbins, yBin);
  int i = icent;
  
  if (etaBin == 0) {
    hmcRawGlobal = (TH2D*)hTemp->Clone(Form("hmcRawGlobal_kSample%d_icent%d",kSample,i));
    hmcRawGlobalTruth = (TH2D*)hTemp->Clone(Form("hmcRawGlobalTruth_kSample%d_icent%d",kSample,i));
    hdataRawGlobal = (TH2D*)hTemp->Clone(Form("hdataRawGlobal_kSample%d_icent%d",kSample,i));
    hmcRawGlobalTruthCorr = (TH2D*)hTemp->Clone(Form("hmcRawGlobalTruthCorr_kSample%d_icent%d",kSample,i));
    hmcRawGlobalCorr = (TH2D*)hTemp->Clone(Form("hdataRawGlobalCorr_kSample%d_icent%d",kSample,i));
    listmcGlobal = new TList;
    listmcGlobalTruth = new TList;
    listdataGlobal = new TList;
  }
  // MC 
  TH2D* hmcRaw   = (TH2D*)hTemp->Clone(Form("hmcRaw_kSample%d_icent%d",kSample,i));
  TH2D* hmcRaw_NoWeight   = (TH2D*)hTemp->Clone(Form("hmcRaw_kSample%d_icent%d_NoWeight",kSample,i));
  TH2D* hmcRaw_jzNorm   = (TH2D*)hTemp->Clone(Form("hmcRaw_kSample%d_icent%d_jzNorm",kSample,i));
  TH2D* hmcRaw_jetWeight   = (TH2D*)hTemp->Clone(Form("hmcRaw_kSample%d_icent%d_jetWeight",kSample,i));
  TH2D* hmcRaw_fcalWeight   = (TH2D*)hTemp->Clone(Form("hmcRaw_kSample%d_icent%d_fcalWeight",kSample,i));
  TH2D* hmcTruth = (TH2D*)hTemp->Clone(Form("hmcTruth_kSample%d_icent%d",kSample,i));
  TH2D* hdataRaw = (TH2D*)hTemp->Clone(Form("hdataRaw_kSample%d_icent%d",kSample,i));
  
  cout << " etaBin " << etaBin << " cent " << icent << endl;

  getMCspectra   ( kSample, icent, etaBin, hmcRaw, hmcRaw_NoWeight, hmcRaw_jetWeight, hmcRaw_jzNorm, hmcRaw_fcalWeight, hmcTruth, /*hmcRawGlobal, hmcRawGlobalTruth,*/ 0, nSys);
  getDATAspectra ( kSample, icent, etaBin, hdataRaw,/* hdataRawGlobal,*/ nSys);

  
  
  listmcGlobal->Add(hmcRaw);
  listmcGlobalTruth->Add(hmcTruth);
  listdataGlobal->Add(hdataRaw);
  
  /*
  hmcRawGlobal->Add(hmcRaw);
  hmcRawGlobalTruth->Add(hmcTruth);
  hdataRawGlobal->Add(hdataRaw);
  
  std::cout << " hmcRawGlobal entries " << hmcRawGlobal->GetEntries() << std::endl;
  std::cout << " hmcRawGlobalTruth entries " << hmcRawGlobalTruth->GetEntries() << std::endl;
  std::cout << " hmcdataGlobal entries " << hdataRawGlobal->GetEntries() << std::endl;
  */
  
  TH2D* hRatioRaw = (TH2D*)hdataRaw->Clone(Form("hRatioRaw_kSample%d_icent%d",kSample,i));
  hRatioRaw->Divide(hmcRaw);
  removeFluc2(hRatioRaw);
  /*
  TH2D* hRatioSmooth2 = (TH2D*)hRatioRaw->Clone(Form("hRatioSmooth2_kSample%d_icent%d",kSample,i));
  hRatioSmooth2->Smooth();

  TH2D* hmcRawSmooth = (TH2D*)hmcRaw->Clone(Form("%s_smooth",hmcRaw->GetName()));
  TH2D* hdataRawSmooth = (TH2D*)hdataRaw->Clone(Form("%s_smooth",hdataRaw->GetName()));
  
  hmcRawSmooth->Smooth();
  hdataRawSmooth->Smooth();
  TH2D* hRatioSmoothRaw = (TH2D*)hdataRawSmooth->Clone(Form("hRatioSmoothRaw_kSample%d_icent%d",kSample,i));
  hRatioSmoothRaw->Divide(hmcRawSmooth);
  
  TH2D* hRatioSmooth = (TH2D*)hRatioSmoothRaw->Clone(Form("hRatioSmooth_kSample%d_icent%d",kSample,i));
  */
    // all reweighing applied
  
  TCanvas* c1=  new TCanvas("c1","",500,500);
  makeEfficiencyCanvas(c1,1, 0.0, 0.01, 0.2, 0.25, 0.01);
  c1->cd(1);
  TH1D* ptmc = (TH1D*)hmcRaw->ProjectionX("ptmc");
  TH1D* ptdata = (TH1D*)hdataRaw->ProjectionX("ptdata");
  handsomeTH1(ptmc,1);
  handsomeTH1(ptdata,2);
  ptmc->SetAxisRange(1e-6,1e8,"Y");
  ptmc->SetYTitle("dN/dp_{T}^{Reco}");
  ptmc->Draw(); 
  ptdata->Draw("same");
  gPad->SetLogy();
  drawCentrality(kSample, icent, 0.6,0.86,1,24);

  TLegend *leg1 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

  easyLeg(leg1,"Mass integrated");
  leg1->AddEntry(ptdata, "Data","pl");
  leg1->AddEntry(ptmc, "MC","pl");
  leg1->Draw();
  ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

  c1->cd(2);
  TH1D* hptRatio = (TH1D*)ptdata->Clone("hptRatio");
  hptRatio->Divide(ptmc);
  if ( kSample == 0) hptRatio->SetAxisRange(0,1e5,"Y");

  hptRatio->SetXTitle("p_{T} (GeV)");
  hptRatio->SetYTitle("Data/MC");
  fixedFontHist(hptRatio,3,2,20);
  hptRatio->SetNdivisions(505,"X");
  hptRatio->Draw();
  
  c1->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d.png",kSample,icent,etaBin));
  c1->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d.root",kSample,icent,etaBin));
    
  TF1* fit;
  
  if ( kSample == kPP ) {
    hptRatio->Fit("pol2","","",100,900);
    fit = new TF1("myFit","[0]+ [1]*x + [2]*x*x",100,900);
    fit->SetParameter(0, hptRatio->GetFunction("pol2")->GetParameter(0));
    fit->SetParameter(1, hptRatio->GetFunction("pol2")->GetParameter(1));
    fit->SetParameter(2, hptRatio->GetFunction("pol2")->GetParameter(2));
  }
  else {
    hptRatio->Fit("pol1","","",100,700);
    fit = new TF1("myFit","[0]+ [1]*x",100,900);
    fit->SetParameter(0, hptRatio->GetFunction("pol1")->GetParameter(0));
    fit->SetParameter(1, hptRatio->GetFunction("pol1")->GetParameter(1));
  }
  fit->SetName("fitFunction");

  c1->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_fit.png",kSample,icent,etaBin));
  delete c1;
  
  // no reweighing
  TCanvas* c110=  new TCanvas("c110","",500,500);
  makeEfficiencyCanvas(c110,1, 0.0, 0.01, 0.2, 0.25, 0.01);
  c110->cd(1);
  TH1D* ptmc_NoWeight = (TH1D*)hmcRaw_NoWeight->ProjectionX("ptmc_NoWeight");
  handsomeTH1(ptmc_NoWeight,1);
  ptmc_NoWeight->SetAxisRange(1e-6,1e8,"Y");
  ptmc_NoWeight->SetYTitle("dN/dp_{T}^{Reco}");
  ptmc_NoWeight->Draw(); 
  ptdata->Draw("same");
  gPad->SetLogy();
  drawCentrality(kSample, icent, 0.6,0.86,1,24);

  TLegend *leg20 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");
  easyLeg(leg20,"Mass integrated");
  leg20->AddEntry(ptdata, "Data","pl");
  leg20->AddEntry(ptmc_NoWeight, "MC","pl");
  leg20->Draw();
  ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

  c110->cd(2);
  TH1D* hptRatio_2 = (TH1D*)ptdata->Clone("hptRatio2");
  hptRatio_2->Divide(ptmc_NoWeight);
  if ( kSample == 0) hptRatio_2->SetAxisRange(0,1e5,"Y");

  hptRatio_2->SetXTitle("p_{T} (GeV)");
  hptRatio_2->SetYTitle("Data/MC");
  fixedFontHist(hptRatio_2,3,2,20);
  hptRatio_2->SetNdivisions(505,"X");
  hptRatio_2->Draw();
  
  //  c1->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_yongsuns.png",kSample,icent,etaBin));
  c110->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_NoWeight.png",kSample,icent,etaBin));
  c110->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_NoWeiht.root",kSample,icent,etaBin));
  delete c110;
  
  // only jet weights
  TCanvas* c111=  new TCanvas("c111","",500,500);
  makeEfficiencyCanvas(c111,1, 0.0, 0.01, 0.2, 0.25, 0.01);
  c111->cd(1);
  TH1D* ptmc_jetWeight = (TH1D*)hmcRaw_jetWeight->ProjectionX("ptmc_jetWeight");
  handsomeTH1(ptmc_jetWeight,1);
  ptmc_jetWeight->SetAxisRange(1e-2,1e13,"Y");
  ptmc_jetWeight->SetYTitle("dN/dp_{T}^{Reco}");
  ptmc_jetWeight->Draw(); 
  ptdata->Draw("same");
  gPad->SetLogy();
  drawCentrality(kSample, icent, 0.6,0.86,1,24);

  TLegend *leg30 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

  easyLeg(leg30,"Mass integrated");
  leg30->AddEntry(ptdata, "Data","pl");
  leg30->AddEntry(ptmc_jetWeight, "MC","pl");
  leg30->Draw();
  ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

  c111->cd(2);
  TH1D* hptRatio_3 = (TH1D*)ptdata->Clone("hptRatio_3");
  hptRatio_3->Divide(ptmc_jetWeight);
  if ( kSample == 0) hptRatio_3->SetAxisRange(0,1e5,"Y");

  hptRatio_3->SetXTitle("p_{T} (GeV)");
  hptRatio_3->SetYTitle("Data/MC");
  fixedFontHist(hptRatio_3,3,2,20);
  hptRatio_3->SetNdivisions(505,"X");
  hptRatio_3->Draw();
  
  //  c1->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_yongsuns.png",kSample,icent,etaBin));
  c111->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_jetWeight.png",kSample,icent,etaBin));
  c111->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_jetWeight.root",kSample,icent,etaBin));
  delete c111;
  // only jzNorm
  TCanvas* c112=  new TCanvas("c112","",500,500);
  makeEfficiencyCanvas(c112,1, 0.0, 0.01, 0.2, 0.25, 0.01);
  c112->cd(1);
  TH1D* ptmc_jzNorm = (TH1D*)hmcRaw_jzNorm->ProjectionX("ptmc_jzNorm");
  handsomeTH1(ptmc_jzNorm,1);
  ptmc_jzNorm->SetAxisRange(1e-6,1e8,"Y");
  ptmc_jzNorm->SetYTitle("dN/dp_{T}^{Reco}");
  ptmc_jzNorm->Draw(); 
  ptdata->Draw("same");
  gPad->SetLogy();
  drawCentrality(kSample, icent, 0.6,0.86,1,24);

  TLegend *leg4 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

  easyLeg(leg4,"Mass integrated");
  leg4->AddEntry(ptdata, "Data","pl");
  leg4->AddEntry(ptmc_jzNorm, "MC","pl");
  leg4->Draw();
  ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

  c112->cd(2);
  TH1D* hptRatio_4 = (TH1D*)ptdata->Clone("hptRatio_4");
  hptRatio_4->Divide(ptmc_jzNorm);
  if ( kSample == 0) hptRatio_4->SetAxisRange(0,1e5,"Y");

  hptRatio_4->SetXTitle("p_{T} (GeV)");
  hptRatio_4->SetYTitle("Data/MC");
  fixedFontHist(hptRatio_4,3,2,20);
  hptRatio_4->SetNdivisions(505,"X");
  hptRatio_4->Draw();
  
  //  c112->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_yongsuns.png",kSample,icent,etaBin));
    c112->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_jzNorm.png",kSample,icent,etaBin));
    c112->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_jzNorm.root",kSample,icent,etaBin));
    delete c112;

  // no fcal weight
  TCanvas* c113=  new TCanvas("c113","",500,500);
  makeEfficiencyCanvas(c113,1, 0.0, 0.01, 0.2, 0.25, 0.01);
  c113->cd(1);
  TH1D* ptmc_fcal = (TH1D*)hmcRaw_fcalWeight->ProjectionX("ptmc_fcal");
  handsomeTH1(ptmc_fcal,1);
  ptmc_fcal->SetAxisRange(1e-6,1e8,"Y");
  ptmc_fcal->SetYTitle("dN/dp_{T}^{Reco}");
  ptmc_fcal->Draw(); 
  ptdata->Draw("same");
  gPad->SetLogy();
  drawCentrality(kSample, icent, 0.6,0.86,1,24);

  TLegend *leg5 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

  easyLeg(leg5,"Mass integrated");
  leg5->AddEntry(ptdata, "Data","pl");
  leg5->AddEntry(ptmc_fcal, "MC","pl");
  leg5->Draw();
  ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

  c113->cd(2); 
  TH1D* hptRatio_5 = (TH1D*)ptdata->Clone("hptRatio_5");
  hptRatio_5->Divide(ptmc_fcal);
  if ( kSample == 0) hptRatio_5->SetAxisRange(0,1e5,"Y");

  hptRatio_5->SetXTitle("p_{T} (GeV)");
  hptRatio_5->SetYTitle("Data/MC");
  fixedFontHist(hptRatio_5,3,2,20);
  hptRatio_5->SetNdivisions(505,"X");
  hptRatio_5->Draw();
  
    c113->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_NofcalWeight.png",kSample,icent,etaBin));
    c113->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_etaBin%d_NofcalWeight.root",kSample,icent,etaBin));        

    delete c113;
  if (etaBin == 2) { 
    // all reweighing applied: all eta bins combined: reco  
    TCanvas* c35=  new TCanvas("c35","",500,500);
    makeEfficiencyCanvas(c1,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c35->cd(1); 
    hmcRawGlobal->Reset();
    hmcRawGlobal->Merge(listmcGlobal);
    hdataRawGlobal->Reset();
    hdataRawGlobal->Merge(listdataGlobal);
    
    TH1D* ptmcGlobal = (TH1D*)hmcRawGlobal->ProjectionX("ptmcGlobal");
    TH1D* ptdataGlobal = (TH1D*)hdataRawGlobal->ProjectionX("ptdataGlobal");
    handsomeTH1(ptmcGlobal,1);
    handsomeTH1(ptdataGlobal,2);
    ptmcGlobal->SetAxisRange(1e-6,1e8,"Y");
    ptmcGlobal->SetYTitle("dN/dp_{T}^{Reco}");
    ptmcGlobal->Draw();
    ptdataGlobal->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.6,0.86,1,24);

    //    std::cout << "ptmcGlobal entries " << ptmcGlobal->GetEntries() << " ptdataGlobal entries " << ptdataGlobal->GetEntries() << std::endl;
    TLegend *leg35 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

    easyLeg(leg35,"Mass integrated");
    leg35->AddEntry(ptdata, "Data","pl");
    leg35->AddEntry(ptmc, "MC","pl");
    leg35->Draw();
    ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

    c35->cd(2);
    TH1D* hptRatioGlobal = (TH1D*)ptdataGlobal->Clone("hptRatioGlobal");
    hptRatioGlobal->Divide(ptmcGlobal);
    if ( kSample == 0) hptRatioGlobal->SetAxisRange(0,1e5,"Y");

    hptRatioGlobal->SetXTitle("p_{T} (GeV)");
    hptRatioGlobal->SetYTitle("Data/MC");
    fixedFontHist(hptRatioGlobal,3,2,20);
    hptRatioGlobal->SetNdivisions(505,"X");
    hptRatioGlobal->Draw();
  
    c35->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d.png",kSample,icent));
    c35->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d.root",kSample,icent));
    delete c35;
        // all reweighing applied: all eta bins combined: truth
    TCanvas* c36=  new TCanvas("c36","",500,500);
    makeEfficiencyCanvas(c36,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c36->cd(1);
    hmcRawGlobalTruth->Reset();
    hmcRawGlobalTruth->Merge(listmcGlobalTruth);
    
    TH1D* ptmcGlobalTruth = (TH1D*)hmcRawGlobalTruth->ProjectionX("ptmcGlobalTruth");
    handsomeTH1(ptmcGlobalTruth,1);
    ptmcGlobalTruth->SetAxisRange(1e-6,1e8,"Y");
    ptmcGlobalTruth->SetYTitle("dN/dp_{T}^{Reco}");
    ptmcGlobalTruth->Draw(); 
    ptdataGlobal->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.6,0.86,1,24);

    TLegend *leg36 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

    easyLeg(leg36,"Mass integrated");
    leg36->AddEntry(ptdataGlobal, "Data","pl");
    leg36->AddEntry(ptmcGlobalTruth, "MC","pl");
    leg36->Draw();
    ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

    c36->cd(2);
    TH1D* hptRatioGlobalTruth = (TH1D*)ptdataGlobal->Clone("hptRatioGlobalTruth");
    hptRatioGlobalTruth->Divide(ptmcGlobalTruth);
    if ( kSample == 0) hptRatioGlobalTruth->SetAxisRange(0,1e5,"Y");

    hptRatioGlobalTruth->SetXTitle("p_{T} (GeV)");
    hptRatioGlobalTruth->SetYTitle("Data/MC");
    fixedFontHist(hptRatioGlobalTruth,3,2,20);
    hptRatioGlobalTruth->SetNdivisions(505,"X");
    hptRatioGlobalTruth->Draw();
  
    c36->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_Truth.png",kSample,icent));
    c36->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_Truth.root",kSample,icent));
    delete c36;

    TCanvas* c37=  new TCanvas("c37","",500,500);
    makeEfficiencyCanvas(c37,1, 0.05, 0.01, 0.1, 0.3, 0.01);
    c37->cd(1);
    TH1D* hmcRawM2Pt2Global = (TH1D*)hmcRawGlobal->ProjectionY("hmcRawM2Pt2Global");
    TH1D* hdataRawM2Pt2Global = (TH1D*)hdataRawGlobal->ProjectionY("hdataRawM2Pt2Global");

    handsomeTH1(hmcRawM2Pt2Global,1);
    handsomeTH1(hdataRawM2Pt2Global,2);
    fixedFontHist(hmcRawM2Pt2Global,2.5,2,20);
    hmcRawM2Pt2Global->Draw();
    hdataRawM2Pt2Global->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.70,0.86,1,24);
    TLegend *leg37 = new TLegend(0.1871854,0.6894175,0.7049639,0.9033743,NULL,"brNDC");
    easyLeg(leg37,"After p_{T} weighting");
    leg37->AddEntry(ptdataGlobal, "Data","pl");
    leg37->AddEntry(ptmcGlobalTruth, "MC","pl");
    leg37->Draw();
    c37->cd(2);
    TH1D* hmptRatioPtGlobal = (TH1D*)hdataRawM2Pt2Global->Clone("hmptRatioPtGlobal");
    hmptRatioPtGlobal->Divide(hmcRawM2Pt2Global);
    TH1D* hmptRatioPtsmoothGlobal = (TH1D*)hmptRatioPtGlobal->Clone("hmptRatioPtsmoothGlobal");

    hmptRatioPtsmoothGlobal->Smooth();
    hmptRatioPtGlobal->SetAxisRange(0,5,"Y");
    hmptRatioPtGlobal->SetXTitle("(m/p_{T}^{Reco})^{2}");
    hmptRatioPtGlobal->SetYTitle("Data/MC");
    hmptRatioPtGlobal->SetNdivisions(505,"X");
    fixedFontHist(hmptRatioPtGlobal,3,2,20);
    hmptRatioPtGlobal->Draw();
    hmptRatioPtsmoothGlobal->Draw("same hist");
    jumSun(-.2,1,0.35,1);
  

    c37->SaveAs(Form("reweightFactors/MassIntreweighting_kSample%d_icent%d.png",kSample,icent));
    c37->SaveAs(Form("reweightFactors/MassIntreweighting_kSample%d_icent%d.root",kSample,icent));
    delete c37;
    
    TCanvas* c38=  new TCanvas("c38","",500,500);
    makeEfficiencyCanvas(c38,1, 0.05, 0.01, 0.1, 0.3, 0.01);
    c38->cd(1);
    TH1D* hmcRawM2Pt2GlobalTruth = (TH1D*)hmcRawGlobalTruth->ProjectionY("hmcRawM2Pt2GlobalTruth");
    //    TH1D* hdataRawM2Pt2Global = (TH1D*)hdataRawGlobal->ProjectionY("hdataRawM2Pt2Global");

    handsomeTH1(hmcRawM2Pt2GlobalTruth,1);
    //    handsomeTH1(hdataRawM2Pt2Global,2);
    fixedFontHist(hmcRawM2Pt2GlobalTruth,2.5,2,20);
    hmcRawM2Pt2GlobalTruth->Draw();
    hdataRawM2Pt2Global->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.70,0.86,1,24);
    TLegend *leg38 = new TLegend(0.1871854,0.6894175,0.7049639,0.9033843,NULL,"brNDC");
    easyLeg(leg38,"After p_{T} weighting");
    leg38->AddEntry(ptdataGlobal, "Data","pl");
    leg38->AddEntry(ptmcGlobalTruth, "MC","pl");
    leg38->Draw();
    c38->cd(2);
    TH1D* hmptRatioPtGlobalTruth = (TH1D*)hdataRawM2Pt2Global->Clone("hmptRatioPtGlobalTruth");
    hmptRatioPtGlobalTruth->Divide(hmcRawM2Pt2GlobalTruth);
    TH1D* hmptRatioPtsmoothGlobalTruth = (TH1D*)hmptRatioPtGlobalTruth->Clone("hmptRatioPtsmoothGlobalTruth");

    hmptRatioPtsmoothGlobalTruth->Smooth();
    hmptRatioPtGlobalTruth->SetAxisRange(0,5,"Y");
    hmptRatioPtGlobalTruth->SetXTitle("(m/p_{T}^{Reco})^{2}");
    hmptRatioPtGlobalTruth->SetYTitle("Data/MC");
    hmptRatioPtGlobalTruth->SetNdivisions(505,"X");
    fixedFontHist(hmptRatioPtGlobalTruth,3,2,20);
    hmptRatioPtGlobalTruth->Draw();
    hmptRatioPtsmoothGlobalTruth->Draw("same hist");
    jumSun(-.2,1,0.35,1);
  

    c38->SaveAs(Form("reweightFactors/MassIntreweighting_Truth_kSample%d_icent%d.png",kSample,icent));
    c38->SaveAs(Form("reweightFactors/MassIntreweighting_Truth_kSample%d_icent%d.root",kSample,icent));
    delete c38;
    
  }
    

  
  /*  
  TH2D* hmcPtCorr   = (TH2D*)hTemp->Clone(Form("hmcRawPtCorr_kSample%d_icent%d",kSample,i));
  TH2D* hmcPtCorr_NoWeight   = (TH2D*)hTemp->Clone(Form("hmcRawPtCorr_kSample%d_icent%d_NoWeight",kSample,i));
  TH2D* hmcPtCorr_jetWeight   = (TH2D*)hTemp->Clone(Form("hmcRawPtCorr_kSample%d_icent%d_jetWeight",kSample,i));
  TH2D* hmcPtCorr_jzNorm   = (TH2D*)hTemp->Clone(Form("hmcRawPtCorr_kSample%d_icent%d_jzNorm",kSample,i));
    TH2D* hmcPtCorr_fcalWeight  = (TH2D*)hTemp->Clone(Form("hmcRawPtCorr_kSample%d_icent%d_fcalWeight",kSample,i));
  TH2D* hmcTruthPtCorr = (TH2D*)hTemp->Clone(Form("hmcTruthPtCorr_kSample%d_icent%d",kSample,i));


  getMCspectra   ( kSample, icent, etaBin,  hmcPtCorr, hmcPtCorr_NoWeight, hmcPtCorr_jetWeight, hmcPtCorr_jzNorm, hmcPtCorr_fcalWeight, hmcTruthPtCorr, hmcRawGlobalCorr, hmcRawGlobalTruthCorr, fit, nSys);
  
  TCanvas* c1ptCorr=  new TCanvas("c1ptCorr","",500,500);
  makeEfficiencyCanvas(c1ptCorr,1, 0.05, 0.01, 0.1, 0.3, 0.01);
  c1ptCorr->cd(1);
  TH1D* ptmcCorr = (TH1D*)hmcPtCorr->ProjectionX("ptmcCorr");
  handsomeTH1(ptmcCorr,1);
  cleverRangeLog(ptmcCorr,10,1e-7);
  ptmcCorr->Draw();
  ptdata->Draw("same");
  gPad->SetLogy();
  c1ptCorr->cd(2);
  TH1D* hptRatioPtCorr = (TH1D*)ptdata->Clone("hptRatioPtcorr");
  hptRatioPtCorr->Divide(ptmcCorr);
  hptRatioPtCorr->Draw();  
  


  Canvas* c15=  new TCanvas("c15","",500,500);
  makeEfficiencyCanvas(c15,1, 0.05, 0.01, 0.1, 0.3, 0.01);
  c15->cd(1);
  TH1D* h1mcAllPt = (TH1D*)hmcPtCorr->ProjectionY("hmcPtCorrAllPt");
  TH1D* h1dataAllPt = (TH1D*)hdataRaw->ProjectionY("hdataRawAllPt");
  float tempMax = cleverRange(h1mcAllPt,10,1);
  //  float tempMax = cleverRange(h1mcAllPt,100,1e-8);
  //  h1mcAllPt->SetAxisRange(1,tempMax*10,"Y");
  handsomeTH1(h1mcAllPt,1);
  handsomeTH1(h1dataAllPt,2);
  fixedFontHist(h1mcAllPt,2.5,2,20);
  h1mcAllPt->Draw();
  h1dataAllPt->Draw("same");
  gPad->SetLogy();
  drawCentrality(kSample, icent, 0.70,0.86,1,24);
  TLegend *leg2 = new TLegend(0.1871854,0.6894175,0.7049639,0.9033743,NULL,"brNDC");
  easyLeg(leg2,"After p_{T} weighting");
  leg2->AddEntry(ptdata, "Data","pl");
  leg2->AddEntry(ptmc, "MC","pl");
  leg2->Draw();
  c15->cd(2);
  TH1D* hmptRatioPtAll = (TH1D*)h1dataAllPt->Clone("hmptRatioPtAll");
  hmptRatioPtAll->Divide(h1mcAllPt);
  TH1D* hmptRatioPtAllsmooth = (TH1D*)hmptRatioPtAll->Clone("hmptRatioPtAllsmooth");

  hmptRatioPtAllsmooth->Smooth();
  hmptRatioPtAll->SetAxisRange(0,5,"Y");
  hmptRatioPtAll->SetXTitle("(m/p_{T}^{Reco})^{2}");
  hmptRatioPtAll->SetYTitle("Data/MC");
  hmptRatioPtAll->SetNdivisions(505,"X");
  fixedFontHist(hmptRatioPtAll,3,2,20);
  hmptRatioPtAll->Draw();
  hmptRatioPtAllsmooth->Draw("same hist");
  jumSun(-.2,1,0.35,1);
  

  c15->SaveAs(Form("reweightFactors/MassIntreweighting_etaBin%d_kSample%d_icent%d.png",etaBin,kSample,icent));
  c15->SaveAs(Form("reweightFactors/MassIntreweighting_etaBin%d_kSample%d_icent%d.root",etaBin,kSample,icent));

   
  TCanvas* c16 = new TCanvas("c16","",500,500);

  TH1D* hsmooVarP = getVariedHist(hmptRatioPtAllsmooth,0.5);
  TH1D* hsmooVarM = getVariedHist(hmptRatioPtAllsmooth,-0.5);

  handsomeTH1( hsmooVarP,4);
  handsomeTH1( hsmooVarM,kGreen+4);
  hsmooVarP->SetLineWidth(2);
  hsmooVarM->SetLineWidth(2);
  hsmooVarP->SetLineStyle(9);
  hsmooVarM->SetLineStyle(6);

  cleverRange(hmptRatioPtAllsmooth,3);
  hmptRatioPtAllsmooth->SetAxisRange(-0.07,0.3,"X");
  if ( kSample == 1) hmptRatioPtAllsmooth->SetAxisRange(-0.18,0.3,"X");
  jumSun(-.2,1,0.35,1);
  hmptRatioPtAllsmooth->SetXTitle("(m/p_{T}^{Reco})^{2}");
  hmptRatioPtAllsmooth->SetYTitle("Data/MC");
  hmptRatioPtAllsmooth->Draw("hist");
  jumSun(-0.18,1,0.3,1);
  hmptRatioPtAll->Draw("same");
  hsmooVarP->Draw("same hist");
  hsmooVarM->Draw("same hist");
  drawCentrality(kSample, icent, 0.50,0.86,1,24);
  
  TLegend *leg3 = new TLegend(0.523618,0.6186667,0.959799,0.832,NULL,"brNDC");
  easyLeg(leg3,"Reweight factor");
  leg3->AddEntry(hsmooVarP, "Varied by +50%","l");
  leg3->AddEntry(hmptRatioPtAll, "Nominal","l");
  leg3->AddEntry(hsmooVarM, "Varied by -50%","l");
  leg3->Draw();

  ATLASLabel(0.18, 0.88, "Internal",0.05);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

  c16->SaveAs(Form("reweightFactors/MassIntreweighting_etaBin%d_kSample%d_icent%d_ratioOnly.png",etaBin,kSample,icent));
 
  //  c16->SaveAs(Form("reweightFactors/MassIntreweighting_etaBin%d_kSample%d_icent%d_ratioOnly_yongsuns.png",etaBin,kSample,icent));  
  
 

  TH1D* hRatio[20];
  TH1D* hRatioSmoPtBin[20];
  TCanvas* c2=  new TCanvas("c2","",1200,600);
  c2->Divide(5,2);
  for ( int ix = 1 ; ix<= nXbins ; ix++) {
    c2->cd(ix);
    TH1D* h1mc = (TH1D*)hmcPtCorr->ProjectionY(Form("hmcPtCorr_%d",ix),ix,ix);
    TH1D* h1data = (TH1D*)hdataRaw->ProjectionY(Form("hdataRaw_%d",ix),ix,ix);
    //    scaleInt(h1data);
    //    scaleInt(h1mc);
    handsomeTH1(h1mc,1);
    handsomeTH1(h1data,2);
    cleverRangeLog(h1mc,10,1e-7);
    h1mc->Draw();
    h1data->Draw("same");
    gPad->SetLogy();

    if ( ix==1)    drawCentrality(kSample, icent, 0.60,0.86,1,24);
    drawBin(xBin,ix,"GeV",0.4,0.78,49,16);
    //    drawBin(xBin,ix,"GeV",0.16 + (0.05* (ix==lowPtBin)), 0.78,1,16);
    
    hRatio[ix] = (TH1D*)h1data->Clone(Form("hRatio_ix%d",ix));
    hRatio[ix]->Divide(h1mc);
    hRatioSmoPtBin[ix] = (TH1D*)hRatio[ix]->Clone(Form("hRatioSmoPtBin_ix%d",ix));
    hRatioSmoPtBin[ix]->Smooth();
  }
  
  TCanvas* c3=  new TCanvas("c3","",1200,600);
  c3->Divide(5,2);
  for ( int ix = 1 ; ix<= nXbins ; ix++) {
    c3->cd(ix);
    hRatio[ix]->SetYTitle("Data/MC");
    //    cleverRangeLog(hRatio[ix],1000,1e-9);
    hRatio[ix]->SetAxisRange(0,5,"Y");
    hRatio[ix]->Draw();
    hRatioSmoPtBin[ix]->Draw("same hist");
    if ( ix==1)    drawCentrality(kSample, icent, 0.60,0.86,1,24);
    drawBin(xBin,ix,"GeV",0.4,0.78,49,16);
    jumSun(-0.5,1,0.5,1);
    //    drawBin(xBin,ix,"GeV",0.16 + (0.05* (ix==lowPtBin)), 0.78,1,16);
    //    hRatio[ix]->Fit("gaus");
    //    gPad->SetLogy();
  }


  TH2D* hFacRatio = new TH2D(Form("factorizedRatio_kSample%d_icent%d",kSample,i),"", nXbins, xBin, nYbins, yBin);
  TH2D* hFacRatioVarP = new TH2D(Form("factorizedRatioVarP_kSample%d_icent%d",kSample,i),"", nXbins, xBin, nYbins, yBin);
  TH2D* hFacRatioVarM = new TH2D(Form("factorizedRatioVarM_kSample%d_icent%d",kSample,i),"", nXbins, xBin, nYbins, yBin);
  TH2D* hFacRatio2 = new TH2D(Form("factorizedRatio2_kSample%d_icent%d",kSample,i),"", nXbins, xBin, nYbins, yBin);
  TH1D* hTemp1x  = (TH1D*)hFacRatio->ProjectionX("htemp1x");
  TH1D* hTemp1y  = (TH1D*)hFacRatio->ProjectionY("htemp1y");
  for ( int ii = 1 ; ii<=hFacRatio->GetNbinsX() ; ii++) {
    for ( int jj = 1 ; jj<=hFacRatio->GetNbinsY() ; jj++) {
      double recoPt = hTemp1x->GetBinCenter(ii);
      double ptWeight = fit->Eval(recoPt);

      double recoM  = hTemp1y->GetBinCenter(jj);

      int theYBin = hmptRatioPtAll->FindBin( recoM);
      double mWeight = hmptRatioPtAll->GetBinContent(theYBin);
      hFacRatio->SetBinContent(ii,jj, ptWeight * mWeight);

      double mWeightVarP = hsmooVarP->GetBinContent(theYBin);
      hFacRatioVarP->SetBinContent(ii,jj, ptWeight * mWeightVarP);

      double mWeightVarM = hsmooVarM->GetBinContent(theYBin);
      hFacRatioVarM->SetBinContent(ii,jj, ptWeight * mWeightVarM);

      int theYBin2 = hRatioSmoPtBin[ii]->FindBin(recoM);
      double mWeight2 = hRatioSmoPtBin[ii]->GetBinContent(theYBin2);
      hFacRatio2->SetBinContent(ii,jj, ptWeight * mWeight2);
    }
  }

  
  TFile * fout;
  if ( nSys < 0 ) 
    //    fout = new TFile(Form("reweightFactors/reweightingFactor_etaBin%d_factorized_v60_yongsuns.root`",etaBin, (float)flucCut),"update");
    fout = new TFile(Form("reweightFactors/reweightingFactor_kSample%d_etaBin%d_factorized_v60.root",kSample,etaBin, (float)flucCut),"update");
  else
    fout = new TFile(Form("reweightFactors/reweightingFactor_kSample%d_etaBin%d_factorized_v60_nSys%d.root",kSample,etaBin, (float)flucCut,nSys),"update");

  hmcPtCorr->Write("",TObject::kOverwrite);
  hmcTruth->Write("",TObject::kOverwrite);
  hdataRaw->Write("",TObject::kOverwrite);
  hmcRawSmooth->Write("",TObject::kOverwrite);
  hdataRawSmooth->Write("",TObject::kOverwrite);
  hdataRawSmooth->Write("",TObject::kOverwrite);
  hRatioSmoothRaw->Write("",TObject::kOverwrite);
  hRatioSmooth->Write("",TObject::kOverwrite);
  hRatioRaw->Write("",TObject::kOverwrite);
  hRatioSmooth2->Write("",TObject::kOverwrite);
  hFacRatio->Write("",TObject::kOverwrite);
  hFacRatioVarP->Write("",TObject::kOverwrite);
  hFacRatioVarM->Write("",TObject::kOverwrite);
  hFacRatio2->Write("",TObject::kOverwrite);

  fout->Close();
  */
  delete hTemp;
}


void getMCspectra(int kSample, int icent, int etaBin,  TH2D* hmcRaw, TH2D* hmcRaw_noWeight, TH2D* hmcRaw_jetWeight, TH2D* hmcRaw_jzNorm, TH2D* hmcRaw_fcalWeight, TH2D* hmcTruth, /*TH2D* hmcRawGlobal, TH2D* hmcRawGlobalTruth,*/ TF1* ptScale, int nSys) {
  
  TRandom3 genRandom;
  genRandom.SetSeed(200);

  TH1::SetDefaultSumw2();
  /*
  if (etaBin == 0) {
    hmcRawGlobal->Reset();
    hmcRawGlobalTruth->Reset();
  }
  */
  hmcRaw->Reset();
  hmcRaw_noWeight->Reset();
  hmcRaw_jetWeight->Reset();
  hmcRaw_jzNorm->Reset();
  hmcRaw_fcalWeight->Reset();
  hmcTruth->Reset();

  TF1* fjmscal[30];
  if( nSys == 300) {   
    TFile* fin = new TFile(Form("fJMScalibration_kSample%d_icent%d_mPt2.root",kSample,icent));
    for ( int ix = lowPtBin ; ix<=highPtBin ; ix++) {
      fjmscal[ix] = (TF1*)fin->Get(Form("f1_kSample%d_icent%d_ix%d",kSample,icent,ix));
    }
  }

  int nXbinsCal;
  double xBinCal[30];
  getXbin(nXbinsCal, xBinCal, 1);
  TH1D* xBinTemp = new TH1D("xBinTemp","", nXbinsCal, xBinCal);


  TString jz2;
  TString jz3;
  TString jz4;
  if ( kSample == kPbPb ) {
    if ( (nSys >= 0)&&(nSys<200)) {
      jz2 = jz2PbPbStringSys;
      jz3 = jz3PbPbStringSys;
      jz4 = jz4PbPbStringSys;
    }
    else {
      jz2 = jz2PbPbString;
      jz3 = jz3PbPbString;
      jz4 = jz4PbPbString;
    }
  }
  else if ( kSample == kPP ) {
    if ( (nSys >= 0)&&(nSys<200)) {
      jz2 = jz2PPStringSys;
      jz3 = jz3PPStringSys;
      jz4 = jz4PPStringSys;
    }
    else {
      jz2 = jz2PPString;
      jz3 = jz3PPString;
      jz4 = jz4PPString;
    }
  }
  
  RtrkProvider rtrkProv;
  rtrkProv.Setup(kSample, icent);



  TH1D* hFcalReweight;
  if ( kSample == kPbPb ) {
    TFile* fcal = new TFile("reweightFactors/FCal_HP_v_MB_weights.root");
    //TFile* fcal = new TFile("reweightFactors/fcalWeight_PbPb_v1.root");
    //    hFcalReweight = (TH1D*)fcal->Get("FCal_HP_v_MBOV_weights");
    hFcalReweight = (TH1D*)fcal->Get("weight_MBov_to_HP");
  }
  
  //  jetSubStr  myJetMc;
  //    jetSubStr  *myJetMc = malloc(sizeof (jetSubStr))ubStr();
  //  TBranch  *b_myJetSubMc;


  float ptSys;
  TBranch *b_ptSys ;
  TString jetSysName = getPtSysName(nSys);
  cout << " jetSysName = " << jetSysName << endl;

  cout << " Setting tree branch address..." << endl;
  TFile* fjz2 = new TFile(Form("../ntuples/%s",jz2.Data()));
  TTree* tr2 = (TTree*)fjz2->Get("tr");
  tashka t_jz2;
  //  tr2->SetBranchAddress("jets", &myJetMc.cent, &b_myJetSubMc);
  if ( (nSys>=0) && (nSys<200) )
    tr2->SetBranchAddress(jetSysName.Data(), &ptSys, &b_ptSys);

  TFile* fjz3 = new TFile(Form("../ntuples/%s",jz3.Data()));
  TTree* tr3 = (TTree*)fjz3->Get("tr");
  tashka t_jz3;
  //  tr3->SetBranchAddress("jets", &(myJetMc.cent), &b_myJetSubMc);
  if ( (nSys>=0) && (nSys<200) )
    tr3->SetBranchAddress(jetSysName.Data(), &ptSys, &b_ptSys);

  TFile* fjz4 = new TFile(Form("../ntuples/%s",jz4.Data()));
  TTree* tr4 = (TTree*)fjz4->Get("tr");
  tashka t_jz4;

  TH1F *ptSpectra = new TH1F("ptSpectra", "ptSpectra",1000,0 ,800 );
  TH1F *etaDist = new TH1F("etaDist","etaDist",1000,-2.5,2.5);
  TH2F *etaPhiMap = new TH2F("etaPhiMap","etaPhiMap",1000,-2.5,2.5,1000,0,TMath::Pi());
  TH2F *h2_ptWeight = new TH2F("ptWeight","ptWeight",1000,0,800,1000,0,0.5);
  
  //  tr4->SetBranchAddress("jets", &(myJetMc.cent), &b_myJetSubMc);
  if ( (nSys>=0) && (nSys<200) )
    tr4->SetBranchAddress(jetSysName.Data(), &ptSys, &b_ptSys);

  long long entries;
  for ( int ijz =2 ; ijz<=4 ; ijz++) {
    tashka tr;
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
    std::cout << "jz " << ijz << " norm " << jzNorm << " nEvents_jz4 " << nEvents_jz4 << std::endl;
    cout << "Scanning JZ"<<ijz<<" file.  Total events = " << entries << endl;
    int jetCount = 0;
    for (Int_t i= 0; i<entries ; i++) {
      if ( i > entries * statFrac ) break;
      if (i%2==0)  continue;

      tr.GetEntry(i);

      if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_rawPt, icent, true) ) // isMC = true
        continue;
      if ( !passEtaCut(tr.jets_eta, etaBin) )
        continue;
      
      double genX;
      double genY;
      double recoX  = -10000; 
      double recoY  = -10000; 
      //      cout << " fcalet " << tr.jets_fcal << " maxTrkPt " << tr.jets_maxTrkPt << " drTrkJetBkg " << tr.jets_drTrkJetBkg << " recoMass " << tr.jets_massRaw2 << " recoChPt " << tr.jets_chPtCSubt << " matchDr " << tr.jets_dr << " genNch " << tr.jets_genNch << " genChSdPt " << tr.jets_genChSdPt << " genChSdMass " << tr.jets_genChSdMass <<  " recoChPtRcSubt " << tr.jets_chPtRcSubt << " recoChMassRcSubt " << tr.jets_chMassRcSubt << endl;
      
      genX = tr.jets_genPt;
      genY = tr.jets_genMass2 /  (tr.jets_genPt * tr.jets_genPt);
      //      genY = tr.jets_genMass /  tr.jets_genPt;      
      recoX = tr.jets_rawPt;
      //      recoY = tr.jets_massRaw2 / tr.jets_rawPt;
      recoY = tr.jets_massRaw2 / (tr.jets_rawPt * tr.jets_rawPt);

      if ( genY < 0 ) genY = 0.000001;
      
      double fcalWeight = 1.0;
      if ( kSample==kPbPb) {
	fcalWeight = hFcalReweight->GetBinContent(hFcalReweight->GetXaxis()->FindBin(tr.jets_fcal));
	//	cout <<" fcal, weight = "<<tr.jets_fcal<<", "<<fcalWeight << " myjetMc.recoPt " << tr.jets_rawPt<<endl;
      }
      
      double ptWeight = 1;
      if ( ptScale != 0) {
	ptWeight = ptScale->Eval(recoX);
      }
      
      if ( (nSys>=0) && (nSys<200) ) {
	double extraPtScale = ptSys / tr.jets_rawPt ;
	recoX = recoX * extraPtScale ; //pt
	tr.jets_rawPt = ptSys;  // New pT!!!
	tr.jets_massRaw2 = tr.jets_massRaw2 * extraPtScale ; // new mass so that m/pT is invariant.
      }
      
      if (nSys==200) { // JMR
	// smear by 20% the recoY
	double theCenter = genY * getJMSscale( kSample, icent, tr.jets_rawPt);
	double recoDev = recoY - theCenter; 
	double theResol = getJMRsigma( kSample, icent, tr.jets_rawPt);
        double jmrUnc = getJmrUnc( kSample, icent, tr.jets_rawPt);
        double theVariation = sqrt ( (1 +jmrUnc)*(1+jmrUnc) - 1 );
	recoY = theCenter + recoDev * genRandom.Gaus(1, theVariation * theResol);  
	//	recoY = theCenter + recoDev * genRandom.Gaus(1, 0.66 * theResol);  
      }
      if (nSys==201) { // JMR HI
	// smear by 20% the recoY
	double theCenter = genY * getJMSscale( kSample, icent, tr.jets_rawPt);
	double recoDev = recoY - theCenter; 
	double theResol = getJMRsigma( kSample, icent, tr.jets_rawPt);
        double jmrUnc = getJmrUncHI( kSample, icent, tr.jets_rawPt);
        double theVariation = sqrt ( (1 +jmrUnc)*(1+jmrUnc) - 1 );
	recoY = theCenter + recoDev * genRandom.Gaus(1, theVariation * theResol);  
	//	recoY = theCenter + recoDev * genRandom.Gaus(1, 0.66 * theResol);  
      }
      /* commented out temporarily -- if using these nsys, fix the myJetMc section)
      if (nSys==210) { // JMS
	double theRtrk = rtrkProv.getR( myJetMc);
	recoY = recoY * theRtrk;
      }
      if (nSys==211) { // JMS
	recoY = recoY * 1.008;
      }
      if (nSys==213) { // JMS by Herwig
	double theRtrk = rtrkProv.getRPyHer(myJetMc);
	recoY = recoY / theRtrk;
      }
      if (nSys==217) { // JMS by Herwig
	recoY = recoY * 1.05;
      }
      if (nSys==300) { // JMS calibration
	int xx = xBinTemp->FindBin( tr.jets_rawPt);
	if ( xx > 11 )  xx = 11;
	if ( xx < 5 )  xx = 5;
	float mptVal = recoY;
	if (mptVal <0 )  mptVal = 0;

	double theFac = 1;
	if ( recoY < 0 )
	  theFac = 1;
	else if ( recoY > 0.25 )
	  theFac = fjmscal[xx]->Eval(0.25);
	else
	  theFac = fjmscal[xx]->Eval(recoY);

	recoY = recoY / theFac;
      }
      */
      //      std::cout << " recoX " << recoX << " recoY " << recoY << " jets_weight " << tr.jets_weight << " jzNorm " <<  jzNorm << "  fcalWeight " << fcalWeight << " ptWeight " << ptWeight << std::endl;
      //recoX: jet pt
      //recoY: m^2/pt^2
      //      hmcRawGlobal->Fill( recoX, recoY, tr.jets_weight * jzNorm * fcalWeight * ptWeight);
      //      hmcRawGlobalTruth->Fill( genX, genY, tr.jets_weight * jzNorm * fcalWeight * ptWeight);
      hmcRaw->Fill( recoX, recoY, tr.jets_weight * jzNorm * fcalWeight * ptWeight);
      //      if (fcalWeight != 1 || ptWeight != 1)  std::cout << "FCAL OR PTWEIGHT != 1!!!!!!   fcal " << fcalWeight << " ptWeight " << ptWeight << std::endl; 
      hmcRaw_noWeight->Fill( recoX, recoY);
      hmcRaw_jetWeight->Fill( recoX, recoY, tr.jets_weight);
      hmcRaw_fcalWeight->Fill( recoX, recoY, tr.jets_weight * jzNorm * ptWeight);
      hmcRaw_jzNorm->Fill( recoX, recoY, jzNorm);	    
      hmcTruth->Fill( genX, genY, tr.jets_weight * jzNorm * fcalWeight * ptWeight);

      etaDist->Fill(tr.jets_y, tr.jets_weight * jzNorm * fcalWeight * ptWeight);
      etaPhiMap->Fill(tr.jets_y, tr.jets_phi, tr.jets_weight * jzNorm * fcalWeight * ptWeight);
      //      cout << " recox " << recoX << " recoY " << recoY << " jet weight " << tr.jets_weight << " jzNorm " << jzNorm << " fcalWeight " << fcalWeight << " ptWeight " << ptWeight << endl;
      //      cout << " total weight " << tr.jets_weight * jzNorm * fcalWeight * ptWeight << endl;
      ++jetCount;
    }
    cout << "jz " << ijz << " jetCount " << jetCount << endl;
    //    ptSpectra->Print("ALL");
    if (ijz == 4) {
      //      cout << " c2" << endl;
      TCanvas *c2 = new TCanvas();
      etaDist->Draw();
      c2->SaveAs(Form("reweightFactors/etaDist_kSample%d.root",kSample));
      delete c2;
      TCanvas *c3 = new TCanvas();
      etaPhiMap->Draw("colz");
      c3->SaveAs(Form("reweightFactors/etaPhiMap_kSample%d.root",kSample));
      delete c3;

    }
  }
  //  std::cout << " hmcRawGlobal entriesss " << hmcRawGlobal->GetEntries() << std::endl;
 
  delete xBinTemp;  
}

void getDATAspectra(int kSample, int icent, int etaBin, TH2D* hdataRaw, /*TH2D* hdataRawGlobal,*/ int nSys) { 

  TH1::SetDefaultSumw2();
  hdataRaw->Reset();
  /*
  if (etaBin == 0) {
    hdataRawGlobal->Reset();
  }
  */
  TString fname;
  if ( kSample == kPbPb ) {
    fname = pbpbDataString; 
  }
  else if ( kSample == kPP) {
    fname = ppDataString;
  }

  TF1* fjmscal[30];
  if( nSys == 300) {
    TFile* fin = new TFile(Form("fJMScalibration_kSample%d_icent%d_num.root",kSample,icent));
    for ( int ix = lowPtBin ; ix<=highPtBin ; ix++) {
      fjmscal[ix] = (TF1*)fin->Get(Form("f1_kSample%d_icent%d_ix%d",kSample,icent,ix));
    }
  }
  int nXbinsCal;
  double xBinCal[30];
  getXbin(nXbinsCal, xBinCal, 1);
  TH1D* xBinTemp = new TH1D("xBinTemp","", nXbinsCal, xBinCal);
  
  
  TFile* fData = new TFile(Form("../ntuples/%s",fname.Data()));
  TTree* tree = (TTree*)fData->Get("tr");
  tashka tr;
  tr.Init(tree);
  //  jetSubStr myJet;
  //  TBranch       *b_myJet;
  //  tr->SetBranchAddress("jets", &(myJet.cent), &b_myJet);
  if ( kSample == kPP )  cout << " pp " ;
  else if ( kSample == kPbPb) cout << " PbPb " ;

  long long entries = tree->GetEntries();
  //  cout << "data entries = " << entries << endl;
  for (Int_t i= 0; i<entries ; i++) {
    tr.GetEntry(i);
    //    cout << " i " << i << " entries " << entries << " statFrac " << statFrac << endl;
    if ( i > entries * statFrac) break;
    //    cout << " cent " << tr.jets_cent << " pt  "<< tr.jets_genPt << endl;
    if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_rawPt, icent, false) ) // isMC = false
      continue;
    if ( !passEtaCut(tr.jets_eta, etaBin) )
      continue;
    
    double recoX  = -10000;
    double recoY  = -10000;
    recoX = tr.jets_rawPt;
    //    recoY = tr.jets_massRaw2 / tr.jets_rawPt;
    recoY = tr.jets_massRaw2 / (tr.jets_rawPt * tr.jets_rawPt);
    
    if (nSys==300) { // JMS calibration
      int xx = xBinTemp->FindBin( tr.jets_rawPt);
      if ( xx > 11 )  xx = 11;
      if ( xx < 5 )  xx = 5;
      float mptVal = recoY;
      if (mptVal <0 )  mptVal = 0;
      double theFac = 1;
      if ( recoY < 0 )
	theFac = 1; 
      else if ( recoY > 0.25 )
	theFac = fjmscal[xx]->Eval(0.25);
      else
	theFac = fjmscal[xx]->Eval(recoY);
      recoY = recoY / theFac;
    }

    //      cout << " data recoX " << recoX << " recoY " << recoY << endl;
    hdataRaw->Fill ( recoX, recoY);
    //    hdataRawGlobal->Fill ( recoX, recoY);
  }
 delete xBinTemp;
  //  std::cout << " hdataRawGlobal entriesss " << hdataRawGlobal->GetEntries() << std::endl;
}


bool isTooSmall(TH2D* hEntries, int recoVarX, int recoVarY, int minEntries) {
  int theBin = hEntries->FindBin(recoVarX, recoVarY);
  if (  hEntries->GetBinContent(theBin) < minEntries )
    return true;
  
  return false;

}


TH1D* getVariedHist(TH1D* hin, double variation)  {
  
  TH1D* hout = (TH1D*)hin->Clone(Form("%s_var%.2f",hin->GetName(),variation));
  hout->Reset();
  for ( int i = 1 ; i<= hin->GetNbinsX() ; i++) {
    double yin = hin->GetBinContent(i);
    double newY = (yin - 1 ) * (1+variation) + 1; 
    hout->SetBinContent(i, newY);
  }
  
  return hout;
  
  
}

void removeFluc2(TH2* h) {
  for ( int i =1 ;  i<=h->GetNbinsX() ; i++) {
    for ( int j =1 ;  j<=h->GetNbinsY() ; j++) {
      double val  = h->GetBinContent(i,j);
      double error  = h->GetBinError(i,j);
      if ( error > val * flucCut )   {
	h->SetBinContent(i,j,1);
	h->SetBinError(i,j, error/val);
      }
    }
  }
}

#ifndef __CINT__
int main (int argc, char *argv[]) {
  for (int i = 0; i < 4; ++i ) {
    cout << " argv " << i << " " << argv[i] << endl;
  }
  hmcRawGlobal->Reset();
  hmcRawGlobalTruth->Reset();
  hmcGlobalCorr->Reset();
  hmcGlobalTruthCorr->Reset();
  hdataRawGlobal->Reset();
  listmcGlobal->Reset();
  listmcGlobalTruth->Reset();
  listmcGlobalCorr->Reset();
  listmcGlobalTruthCorr->Reset();
  listdataGlobal->Reset();
  
  for (int etaBin = 0; etaBin<3; etaBin++){
    int kSample = atoi(argv[1]);
    int icent = atoi(argv[2]);
    //        int etaBin = atoi(argv[3]);
    int nSys = atoi(argv[3]);
    //	    int nSys = atoi(argv[4]);
    getMcWeights(kSample,icent,etaBin,nSys);
    //  }
  return 0;
}  // Main program when run stand-alone
#endif
