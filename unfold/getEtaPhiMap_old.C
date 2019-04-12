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
#include "../dataWeight.h"
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
void getMCspectra(int kSample=kPP, int icent=0, int etaBin=0, TH2D* hmcRaw=0, TH2D* hmcTruth=0, TH2D* hm2pt=0, TH2D* hm2ptTruth=0, TH2D* hmcRawGlobal_jzNorm=0, TH2D* hmcRawGlobal_fcalWeight=0, TH2D* hmcRawGlobal_jetWeight=0, TH2D* hmcRawGlobal_noWeight=0, TF1* ptScale=0, int nSys=-1);
void getDATAspectra(int kSample=kPP, int icent=0, int etaBin=0, TH2D* hdataRaw=0,int nSys=0);

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
TH2D* hmcRawGlobal_jzNorm;
TH2D* hmcRawGlobal_fcalWeight;
TH2D* hmcRawGlobal_jetWeight;
TH2D* hmcRawGlobal_noWeight;
TH2D* hmcRawGlobalTruth;
TH2D* hmcGlobalCorr;
TH2D* hmcGlobalTruthCorr;
TH2D* hdataRawGlobal;
TH2D* hm2ptRawGlobal;
TH2D* hm2ptRawGlobalTruth;


void getEtaPhiMap() { 

  gStyle->SetOptStat(0);
  //Variable binning
  TH2D* etaPhi = new TH2D("emp","", nXbins, xBin, nYbins, yBin);
  //Fixed binning
  //  TH2D* hTemp = new TH2D("hptTemp","", nXbins,100,1000, nYbins, yBin);
  std::cout << " created variable binning hist " << std::endl;
  int i = icent;
    
  if (etaBin == 0) {
    hmcRawGlobal = (TH2D*)hTemp->Clone(Form("hmcRawGlobal_kSample%d_icent%d",kSample,i));
    hmcRawGlobal_jzNorm = (TH2D*)hTemp->Clone(Form("hmcRawGlobal_jzNorm_kSample%d_icent%d",kSample,i));
    hmcRawGlobal_fcalWeight = (TH2D*)hTemp->Clone(Form("hmcRawGlobal_fcalWeight_kSample%d_icent%d",kSample,i));
    hmcRawGlobal_jetWeight = (TH2D*)hTemp->Clone(Form("hmcRawGlobal_jetWeight_kSample%d_icent%d",kSample,i));
    hmcRawGlobal_noWeight = (TH2D*)hTemp->Clone(Form("hmcRawGlobal_noWeight_kSample%d_icent%d",kSample,i));
    hmcRawGlobalTruth = (TH2D*)hTemp->Clone(Form("hmcRawGlobalTruth_kSample%d_icent%d",kSample,i));
    hmcGlobalCorr = (TH2D*)hTemp->Clone(Form("hmcGlobalCorr_kSample%d_icent%d",kSample,i));
    hmcGlobalTruthCorr = (TH2D*)hTemp->Clone(Form("hmcGlobalTruthCorr_kSample%d_icent%d",kSample,i));
    hdataRawGlobal = (TH2D*)hTemp->Clone(Form("hdataRawGlobal_kSample%d_icent%d",kSample,i));
    // pp: 0-400, 0-300
    hm2ptRawGlobal = new TH2D("hm2ptRawGlobal", "hm2pt",500,-300,2000,300,0,800);
    //PbPb: 0-600, 0-300
    hm2ptRawGlobalTruth = new TH2D("hm2ptRawGlobalTruth", "hm2ptTruth",500,-300,2000,300,0,800);

    
  }
  // MC 
  TH2D* hmcRaw   = (TH2D*)hTemp->Clone(Form("hmcRaw_kSample%d_icent%d",kSample,i));
  TH2D* hmcRaw_jzNorm   = (TH2D*)hTemp->Clone(Form("hmcRaw_jzNorm_kSample%d_icent%d",kSample,i));
  TH2D* hmcRaw_fcalWeight   = (TH2D*)hTemp->Clone(Form("hmcRaw_fcalWeight_kSample%d_icent%d",kSample,i));
  TH2D* hmcRaw_jetWeight   = (TH2D*)hTemp->Clone(Form("hmcRaw_jetWeight_kSample%d_icent%d",kSample,i));
  TH2D* hmcRaw_noWeight   = (TH2D*)hTemp->Clone(Form("hmcRaw_noWeight_kSample%d_icent%d",kSample,i));
  TH2D* hmcTruth = (TH2D*)hTemp->Clone(Form("hmcTruth_kSample%d_icent%d",kSample,i));
  TH2D* hdataRaw = (TH2D*)hTemp->Clone(Form("hdataRaw_kSample%d_icent%d",kSample,i));
  TH2D* hm2pt   = (TH2D*)hm2ptRawGlobal->Clone(Form("hm2pt_kSample%d_icent%d",kSample,i));
  TH2D* hm2ptTruth = (TH2D*)hm2ptRawGlobal->Clone(Form("hm2ptTruth_kSample%d_icent%d",kSample,i));
  
  cout << " etaBin " << etaBin << " cent " << icent << endl;

  getMCspectra   ( kSample, icent, etaBin, hmcRaw, hmcTruth, hm2pt, hm2ptTruth, hmcRaw_jzNorm, hmcRaw_fcalWeight, hmcRaw_jetWeight, hmcRaw_noWeight, 0, nSys);
  getDATAspectra ( kSample, icent, etaBin, hdataRaw,nSys);

  std::cout << " hmcRaw entries " << hmcRaw->GetEntries() << std::endl;
  std::cout << " hmcTruth entries " << hmcTruth->GetEntries() << std::endl;
  std::cout << " hdata entries " << hdataRaw->GetEntries() << std::endl;
  std::cout << " hm2pt entries " << hm2pt->GetEntries() << std::endl;
  std::cout << " hm2pt truth entries " << hm2ptTruth->GetEntries() << std::endl;
  listmcGlobal->Add(hmcRaw);
  listmcGlobal_jzNorm->Add(hmcRaw_jzNorm);
  listmcGlobal_fcalWeight->Add(hmcRaw_fcalWeight);
  listmcGlobal_jetWeight->Add(hmcRaw_jetWeight);
  listmcGlobal_noWeight->Add(hmcRaw_noWeight);
  listmcGlobalTruth->Add(hmcTruth);
  listdataGlobal->Add(hdataRaw);
  listm2ptGlobal->Add(hm2pt);
  listm2ptGlobalTruth->Add(hm2ptTruth);
      // all reweighing applied
  
  TCanvas* c1=  new TCanvas("c1","",500,500);
  makeEfficiencyCanvas(c1,1, 0.0, 0.01, 0.2, 0.25, 0.01);
  c1->cd(1);
  TH1D* ptmc = (TH1D*)hmcRaw->ProjectionX("ptmc");
  TH1D* ptdata = (TH1D*)hdataRaw->ProjectionX("ptdata");
  for (int i = 0; i < nXbins; ++i) {
    std::cout << " bin " << i << " lower bin edge " << ptmc->GetBinLowEdge(i+1) << " upper bin edge " << ptmc->GetBinLowEdge(i+1) + ptmc->GetBinWidth(i+1) << std::endl;
  }
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
  if ( kSample == 0) hptRatio->SetAxisRange(0,5e4,"Y");

  hptRatio->SetXTitle("p_{T} (GeV)");
  hptRatio->SetYTitle("Data/MC");
  fixedFontHist(hptRatio,3,2,20);
  hptRatio->SetNdivisions(505,"X");
  hptRatio->Draw();
  
  c1->SaveAs(Form("reweightFactors/Plots/pTreweighting_kSample%d_icent%d_etaBin%d.png",kSample,icent,etaBin));
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

  c1->SaveAs(Form("reweightFactors/Plots/pTreweighting_kSample%d_icent%d_etaBin%d_fit.png",kSample,icent,etaBin));


      // all reweighing applied: Truth
  
  TCanvas* c111=  new TCanvas("c111","",500,500);
  makeEfficiencyCanvas(c111,1, 0.0, 0.01, 0.2, 0.25, 0.01);
  c111->cd(1);
  TH1D* ptmcTruth = (TH1D*)hmcTruth->ProjectionX("ptmcTruth");
  handsomeTH1(ptmcTruth,1);
  ptmcTruth->SetAxisRange(1e-6,1e8,"Y");
  ptmcTruth->SetYTitle("dN/dp_{T}^{Reco}");
  ptmcTruth->Draw(); 
  ptdata->Draw("same");
  gPad->SetLogy();
  drawCentrality(kSample, icent, 0.6,0.86,1,24);

  TLegend *leg111 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

  easyLeg(leg111,"Mass integrated");
  leg111->AddEntry(ptdata, "Data","pl");
  leg111->AddEntry(ptmc, "MC","pl");
  leg111->Draw();
  ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

  c111->cd(2);
  TH1D* hptTruthRatio = (TH1D*)ptdata->Clone("hptTruthRatio");
  hptTruthRatio->Divide(ptmcTruth);
  if ( kSample == 0) hptRatio->SetAxisRange(0,1e4,"Y");

  hptTruthRatio->SetXTitle("p_{T} (GeV)");
  hptTruthRatio->SetYTitle("Data/MC");
  fixedFontHist(hptTruthRatio,3,2,20);
  hptTruthRatio->SetNdivisions(505,"X");
  hptTruthRatio->Draw();
  
  c111->SaveAs(Form("reweightFactors/Plots/pTreweighting_Truth_kSample%d_icent%d_etaBin%d.png",kSample,icent,etaBin));
  c111->SaveAs(Form("reweightFactors/pTreweighting_Truth_kSample%d_icent%d_etaBin%d.root",kSample,icent,etaBin));

  


   

  TH2D* hmcPtCorr   = (TH2D*)hTemp->Clone(Form("hmcRawPtCorr_kSample%d_icent%d",kSample,i));
  TH2D* hmcPtCorr_jzNorm   = (TH2D*)hTemp->Clone(Form("hmcRawPtCorr_jzNorm_kSample%d_icent%d",kSample,i));
  TH2D* hmcPtCorr_fcalWeight   = (TH2D*)hTemp->Clone(Form("hmcRawPtCorr_fcalWeight_kSample%d_icent%d",kSample,i));
  TH2D* hmcPtCorr_jetWeight   = (TH2D*)hTemp->Clone(Form("hmcRawPtCorr_jetWeight_kSample%d_icent%d",kSample,i));
  TH2D* hmcPtCorr_noWeight   = (TH2D*)hTemp->Clone(Form("hmcRawPtCorr_noWeight_kSample%d_icent%d",kSample,i));
  TH2D* hmcTruthPtCorr = (TH2D*)hTemp->Clone(Form("hmcTruthPtCorr_kSample%d_icent%d",kSample,i));
  TH2D* hm2ptCorr   = (TH2D*)hTemp->Clone(Form("hm2ptCorr_kSample%d_icent%d",kSample,i));
  TH2D* hm2ptCorrTruth = (TH2D*)hTemp->Clone(Form("hm2ptCorrTruth_kSample%d_icent%d",kSample,i));


  getMCspectra   ( kSample, icent, etaBin,  hmcPtCorr, hmcTruthPtCorr, hm2ptCorr, hm2ptCorrTruth,hmcPtCorr_jzNorm, hmcPtCorr_fcalWeight, hmcPtCorr_jetWeight, hmcPtCorr_noWeight,fit, nSys);
  listmcGlobalCorr->Add(hmcPtCorr);
  listmcGlobalTruthCorr->Add(hmcTruthPtCorr);

  
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
  c1ptCorr->SaveAs(Form("reweightFactors/Plots/pTreweightingCorr_kSample%d_icent%d.png",kSample,icent));


  TCanvas* c15=  new TCanvas("c15","",500,500);
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
  TLegend *leg2 = new TLegend(0.3071854,0.6894175,0.80049639,0.9033743,NULL,"brNDC");
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
  

  c15->SaveAs(Form("reweightFactors/Plots/MassIntreweighting_etaBin%d_kSample%d_icent%d.png",etaBin,kSample,icent));
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

  ATLASLabel(0.18, 0.88, "Internal",0.05);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);v

  c16->SaveAs(Form("reweightFactors/Plots/MassIntreweighting_etaBin%d_kSample%d_icent%d_ratioOnly.png",etaBin,kSample,icent));
 
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

  if (etaBin == 2) { 
    // all reweighing applied: all eta bins combined: reco

        hmcRawGlobal->Reset();
    hmcRawGlobal->Merge(listmcGlobal);
    hmcRawGlobal_jzNorm->Reset();
    hmcRawGlobal_jzNorm->Merge(listmcGlobal_jzNorm);
    hmcRawGlobal_jetWeight->Reset();
    hmcRawGlobal_jetWeight->Merge(listmcGlobal_jetWeight);
    hmcRawGlobal_fcalWeight->Reset();
    hmcRawGlobal_fcalWeight->Merge(listmcGlobal_fcalWeight);
    hmcRawGlobal_noWeight->Reset();
    hmcRawGlobal_noWeight->Merge(listmcGlobal_noWeight);
    hmcRawGlobalTruth->Reset();
    hmcRawGlobalTruth->Merge(listmcGlobalTruth);
    hdataRawGlobal->Reset();
    hdataRawGlobal->Merge(listdataGlobal);
    hmcGlobalCorr->Reset();
    hmcGlobalCorr->Merge(listmcGlobalCorr);
    hmcGlobalTruthCorr->Reset();
    hmcGlobalTruthCorr->Merge(listmcGlobalTruthCorr);

    
    TH1D* ptmcGlobal = (TH1D*)hmcRawGlobal->ProjectionX("ptmcGlobal");
    TH1D* ptmcGlobal_fcalWeight = (TH1D*)hmcRawGlobal_fcalWeight->ProjectionX("ptmcGlobal_fcalWeight");
    TH1D* ptmcGlobal_jzNorm = (TH1D*)hmcRawGlobal_jzNorm->ProjectionX("ptmcGlobal_jzNorm");
    TH1D* ptmcGlobal_jetWeight = (TH1D*)hmcRawGlobal_jetWeight->ProjectionX("ptmcGlobal_jetWeight");
    TH1D* ptmcGlobal_noWeight = (TH1D*)hmcRawGlobal_noWeight->ProjectionX("ptmcGlobal_noWeight");
    TH1D* ptdataGlobal = (TH1D*)hdataRawGlobal->ProjectionX("ptdataGlobal");
    TH1D* ptmcGlobalTruth = (TH1D*)hmcRawGlobalTruth->ProjectionX("ptmcGlobalTruth");
    TH1D* hmcRawM2Pt2GlobalTruth = (TH1D*)hmcRawGlobalTruth->ProjectionY("hmcRawM2Pt2GlobalTruth");
    TH1D* hmcRawM2Pt2Global = (TH1D*)hmcRawGlobal->ProjectionY("hmcRawM2Pt2Global");
    TH1D* hdataRawM2Pt2Global = (TH1D*)hdataRawGlobal->ProjectionY("hdataRawM2Pt2Global");
    TH1D* hmcM2Pt2GlobalTruthCorr = (TH1D*)hmcGlobalTruthCorr->ProjectionY("hmcM2Pt2GlobalTruthCorr");
    TH1D* hmcM2Pt2GlobalCorr = (TH1D*)hmcGlobalCorr->ProjectionY("hmcM2Pt2GlobalCorr");
    TH1D* ptmcGlobalCorr = (TH1D*)hmcGlobalCorr->ProjectionX("ptmcGlobalCorr");

    
    std::cout << "entries in mc mpt " << hmcRawM2Pt2Global->GetEntries() << std::endl;
    std::cout << "entries in ptmcGlobal no weight " << ptmcGlobal_noWeight->GetEntries() << std::endl;
    //rerun tomorrow and test from here
    std::cout << "entries in ptmcGlobal jet weight " << ptmcGlobal_jetWeight->GetEntries() << std::endl;
    std::cout << "entries in data mpt " << hdataRawM2Pt2Global->GetEntries() << std::endl;
    std::cout << "entries in mcTruth mpt " << hmcRawM2Pt2GlobalTruth->GetEntries() << std::endl;
    std::cout << "entries in mc pt " << ptmcGlobal->GetEntries() << std::endl;
    std::cout << "entries in data pt " << ptdataGlobal->GetEntries() << std::endl;
    std::cout << "entries in mcTruth pt " << ptmcGlobalTruth->GetEntries() << std::endl;

    TCanvas* c35=  new TCanvas("c35","",500,500);
    makeEfficiencyCanvas(c35,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c35->cd(1); 
    //    ptmcGlobal->Print("ALL");
    //    hmcRawM2Pt2Global->Print("ALL");
    
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
    leg35->AddEntry(ptdataGlobal, "Data","pl");
    leg35->AddEntry(ptmcGlobal, "MC Reco","pl");
    leg35->Draw();
    ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

    c35->cd(2);
    TH1D* hptRatioGlobal = (TH1D*)ptdataGlobal->Clone("hptRatioGlobal");
    hptRatioGlobal->Divide(ptmcGlobal);
    if ( kSample == 0) hptRatioGlobal->SetAxisRange(0,4e4,"Y");

    hptRatioGlobal->SetXTitle("p_{T} (GeV)");
    hptRatioGlobal->SetYTitle("Data/MC");
    fixedFontHist(hptRatioGlobal,3,2,20);
    hptRatioGlobal->SetNdivisions(505,"X");
    hptRatioGlobal->Draw();
  
    c35->SaveAs(Form("reweightFactors/Plots/pTreweighting_kSample%d_icent%d.png",kSample,icent));
    c35->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d.root",kSample,icent));

        // all reweighing applied: all eta bins combined: truth
    TCanvas* c36=  new TCanvas("c36","",500,500);
    makeEfficiencyCanvas(c36,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c36->cd(1);
    
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
    leg36->AddEntry(ptmcGlobalTruth, "MC Truth","pl");
    leg36->Draw();
    ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

    c36->cd(2);
    TH1D* hptRatioGlobalTruth = (TH1D*)ptdataGlobal->Clone("hptRatioGlobalTruth");
    hptRatioGlobalTruth->Divide(ptmcGlobalTruth);
    if ( kSample == 0) hptRatioGlobalTruth->SetAxisRange(0,1e5,"Y");

    hptRatioGlobalTruth->SetXTitle("p_{T}^{Reco} (GeV)");
    hptRatioGlobalTruth->SetYTitle("Data/MC");
    fixedFontHist(hptRatioGlobalTruth,3,2,20);
    hptRatioGlobalTruth->SetNdivisions(505,"X");
    hptRatioGlobalTruth->Draw();
  
    c36->SaveAs(Form("reweightFactors/Plots/pTreweighting_kSample%d_icent%d_Truth.png",kSample,icent));
    c36->SaveAs(Form("reweightFactors/pTreweighting_kSample%d_icent%d_Truth.root",kSample,icent));

    TCanvas* c37=  new TCanvas("c37","",500,500);
    makeEfficiencyCanvas(c37,1, 0.05, 0.01, 0.1, 0.3, 0.01);
    c37->cd(1);
    
    handsomeTH1(hmcRawM2Pt2Global,1);
    handsomeTH1(hdataRawM2Pt2Global,2);
    fixedFontHist(hmcRawM2Pt2Global,2.5,2,20);
    hmcRawM2Pt2Global->SetAxisRange(5e-2,5e6,"Y");
    hmcRawM2Pt2Global->Draw();
    hdataRawM2Pt2Global->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.70,0.86,1,24);
    TLegend *leg37 = new TLegend(0.3871854,0.6894175,0.7849639,0.9033743,NULL,"brNDC");
    easyLeg(leg37,"Before p_{T} weighting");
    leg37->AddEntry(hdataRawM2Pt2Global, "Data","pl");
    leg37->AddEntry(hmcRawM2Pt2Global, "MC Reco","pl");
    leg37->Draw();
    c37->cd(2);
    TH1D* hmptRatioPtGlobal = (TH1D*)hdataRawM2Pt2Global->Clone("hmptRatioPtGlobal");
    hmptRatioPtGlobal->Divide(hmcRawM2Pt2Global);
    TH1D* hmptRatioPtsmoothGlobal = (TH1D*)hmptRatioPtGlobal->Clone("hmptRatioPtsmoothGlobal");

    hmptRatioPtsmoothGlobal->Smooth();
    hmptRatioPtGlobal->SetAxisRange(0,70000,"Y");
    hmptRatioPtGlobal->SetXTitle("(m/p_{T}^{Reco})^{2}");
    hmptRatioPtGlobal->SetYTitle("Data/MC");
    hmptRatioPtGlobal->SetNdivisions(505,"X");
    fixedFontHist(hmptRatioPtGlobal,3,2,20);
    hmptRatioPtGlobal->Draw();
    hmptRatioPtsmoothGlobal->Draw("same hist");
    jumSun(-.2,1,0.35,1);
  

    c37->SaveAs(Form("reweightFactors/Plots/MassIntreweighting_kSample%d_icent%d.png",kSample,icent));
    c37->SaveAs(Form("reweightFactors/MassIntreweighting_kSample%d_icent%d.root",kSample,icent));
 
    
    TCanvas* c38=  new TCanvas("c38","",500,500);
    makeEfficiencyCanvas(c38,1, 0.05, 0.01, 0.1, 0.3, 0.01);
    c38->cd(1);

    //    TH1D* hdataRawM2Pt2Global = (TH1D*)hdataRawGlobal->ProjectionY("hdataRawM2Pt2Global");

    handsomeTH1(hmcRawM2Pt2GlobalTruth,1);
    fixedFontHist(hmcRawM2Pt2GlobalTruth,2.5,2,20);
    hmcRawM2Pt2GlobalTruth->SetAxisRange(5e-3,5e6,"Y");
    hmcRawM2Pt2GlobalTruth->Draw();
    hdataRawM2Pt2Global->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.70,0.86,1,24);
    TLegend *leg38 = new TLegend(0.3871854,0.6894175,0.7849639,0.9033843,NULL,"brNDC");
    easyLeg(leg38,"Before p_{T} weighting");
    leg38->AddEntry(hdataRawM2Pt2Global, "Data","pl");
    leg38->AddEntry(hmcRawM2Pt2GlobalTruth, "MC Truth","pl");
    leg38->Draw();
    c38->cd(2);
    TH1D* hmptRatioPtGlobalTruth = (TH1D*)hdataRawM2Pt2Global->Clone("hmptRatioPtGlobalTruth");
    hmptRatioPtGlobalTruth->Divide(hmcRawM2Pt2GlobalTruth);
    TH1D* hmptRatioPtsmoothGlobalTruth = (TH1D*)hmptRatioPtGlobalTruth->Clone("hmptRatioPtsmoothGlobalTruth");

    hmptRatioPtsmoothGlobalTruth->Smooth();
    hmptRatioPtGlobalTruth->SetAxisRange(0,70000,"Y");
    hmptRatioPtGlobalTruth->SetXTitle("(m/p_{T}^{Truth})^{2}");
    hmptRatioPtGlobalTruth->SetYTitle("Data/MC");
    hmptRatioPtGlobalTruth->SetNdivisions(505,"X");
    fixedFontHist(hmptRatioPtGlobalTruth,3,2,20);
    hmptRatioPtGlobalTruth->Draw();
    hmptRatioPtsmoothGlobalTruth->Draw("same hist");
    jumSun(-.2,1,0.35,1);
  

    c38->SaveAs(Form("reweightFactors/Plots/MassIntreweighting_Truth_kSample%d_icent%d.png",kSample,icent));
    c38->SaveAs(Form("reweightFactors/MassIntreweighting_Truth_kSample%d_icent%d.root",kSample,icent));



    TCanvas* c39=  new TCanvas("c39","",500,500);
    makeEfficiencyCanvas(c39,1, 0.05, 0.01, 0.1, 0.3, 0.01);
    c39->cd(1);
    
    handsomeTH1(hmcM2Pt2GlobalCorr,1);
    fixedFontHist(hmcM2Pt2GlobalCorr,2.5,2,20);
    hmcM2Pt2GlobalCorr->SetAxisRange(5e-2,5e6,"Y");
    hmcM2Pt2GlobalCorr->Draw();
    hdataRawM2Pt2Global->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.70,0.86,1,24);
    TLegend *leg39 = new TLegend(0.3871854,0.6894175,0.7849639,0.9033743,NULL,"brNDC");
    easyLeg(leg39,"After p_{T} weighting");
    leg39->AddEntry(hdataRawM2Pt2Global, "Data","pl");
    leg39->AddEntry(hmcM2Pt2GlobalCorr, "MC Reco","pl");
    leg39->Draw();
    c39->cd(2);
    TH1D* hmptRatioPtGlobalCorr = (TH1D*)hdataRawM2Pt2Global->Clone("hmptRatioPtGlobalCorr");
    hmptRatioPtGlobalCorr->Divide(hmcM2Pt2GlobalCorr);
    TH1D* hmptRatioPtsmoothGlobalCorr = (TH1D*)hmptRatioPtGlobalCorr->Clone("hmptRatioPtsmoothGlobalCorr");

    hmptRatioPtsmoothGlobalCorr->Smooth();
    hmptRatioPtGlobalCorr->SetAxisRange(0,5,"Y");
    hmptRatioPtGlobalCorr->SetXTitle("(m/p_{T}^{Reco})^{2}");
    hmptRatioPtGlobalCorr->SetYTitle("Data/MC");
    hmptRatioPtGlobalCorr->SetNdivisions(505,"X");
    fixedFontHist(hmptRatioPtGlobalCorr,3,2,20);
    hmptRatioPtGlobalCorr->Draw();
    hmptRatioPtsmoothGlobalCorr->Draw("same hist");
    jumSun(-.2,1,0.35,1);
  

    c39->SaveAs(Form("reweightFactors/Plots/MassIntreweightingCorr_kSample%d_icent%d.png",kSample,icent));
    c39->SaveAs(Form("reweightFactors/MassIntreweightingCorr_kSample%d_icent%d.root",kSample,icent));
 
    
    TCanvas* c40=  new TCanvas("c40","",500,500);
    makeEfficiencyCanvas(c40,1, 0.05, 0.01, 0.1, 0.3, 0.01);
    c40->cd(1);

    handsomeTH1(hmcM2Pt2GlobalTruthCorr,1);
    
    fixedFontHist(hmcM2Pt2GlobalTruthCorr,2.5,2,20);
    hmcM2Pt2GlobalTruthCorr->SetAxisRange(5e-3,5e6,"Y");
    hmcM2Pt2GlobalTruthCorr->Draw();
    hdataRawM2Pt2Global->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.70,0.86,1,24);
    TLegend *leg40 = new TLegend(0.3871854,0.6894175,0.7849639,0.9033843,NULL,"brNDC");
    easyLeg(leg40,"After p_{T} weighting");
    leg40->AddEntry(hdataRawM2Pt2Global, "Data","pl");
    leg40->AddEntry(hmcM2Pt2GlobalTruthCorr, "MC","pl");
    leg40->Draw();
    c40->cd(2);
    TH1D* hmptRatioPtGlobalTruthCorr = (TH1D*)hdataRawM2Pt2Global->Clone("hmptRatioPtGlobalTruthCorr");
    hmptRatioPtGlobalTruthCorr->Divide(hmcM2Pt2GlobalTruthCorr);
    TH1D* hmptRatioPtsmoothGlobalTruthCorr = (TH1D*)hmptRatioPtGlobalTruthCorr->Clone("hmptRatioPtsmoothGlobalTruthCorr");

    hmptRatioPtsmoothGlobalTruthCorr->Smooth();
    hmptRatioPtGlobalTruthCorr->SetAxisRange(0,5,"Y");
    hmptRatioPtGlobalTruthCorr->SetXTitle("(m/p_{T}^{Reco})^{2}");
    hmptRatioPtGlobalTruthCorr->SetYTitle("Data/MC");
    hmptRatioPtGlobalTruthCorr->SetNdivisions(505,"X");
    fixedFontHist(hmptRatioPtGlobalTruthCorr,3,2,20);
    hmptRatioPtGlobalTruthCorr->Draw();
    hmptRatioPtsmoothGlobalTruthCorr->Draw("same hist");
    jumSun(-.2,1,0.35,1);
  
    c40->SaveAs(Form("reweightFactors/Plots/MassIntreweighting_TruthCorr_kSample%d_icent%d.png",kSample,icent));
    c40->SaveAs(Form("reweightFactors/MassIntreweighting_TruthCorr_kSample%d_icent%d.root",kSample,icent));


    TCanvas* c41=  new TCanvas("c41","",500,500);
    handsomeTH2(hm2pt);
    fixedFontHist(hm2pt,2.5,2,20);
    //    hmcM2Pt2GlobalTruthCorr->SetAxisRange(5e-3,5e6,"Y");
    std::cout << "hm2pt entries " << hm2pt->GetEntries() << std::endl;
    hm2pt->Draw("colz");
    hm2pt->SetContour(99);
    hm2pt->GetXaxis()->SetTitle("m^{2}");
    hm2pt->GetYaxis()->SetTitle("p_{T}");
    hm2pt->GetXaxis()->SetTitle("m^{2}");
    float entries = (float) hm2pt->GetEntries();
    hm2pt->Scale(1/entries);
    drawCentrality(kSample, icent, 0.70,0.86,1,24);
    c41->SaveAs(Form("reweightFactors/Plots/m2ptCorrelation_kSample%d_icent%d.png",kSample,icent));
    c41->SaveAs(Form("reweightFactors/m2ptCorrelation_kSample%d_icent%d.root",kSample,icent));

    TCanvas* c42=  new TCanvas("c42","",500,500);
    handsomeTH2(hm2ptTruth);
    fixedFontHist(hm2ptTruth,2.5,2,20);
    //    hmcM2Pt2GlobalTruthCorr->SetAxisRange(5e-3,5e6,"Y");
    std::cout << "hm2ptTruth entries " << hm2ptTruth->GetEntries() << std::endl;
    hm2ptTruth->Draw("colz");
    hm2ptTruth->SetContour(99);
    hm2ptTruth->GetXaxis()->SetTitle("m^{2}");
    hm2ptTruth->GetYaxis()->SetTitle("p_{T}");
    entries = (float) hm2ptTruth->GetEntries();
    hm2ptTruth->Scale(1/entries);
    drawCentrality(kSample, icent, 0.70,0.86,1,24);
    c42->SaveAs(Form("reweightFactors/Plots/m2ptTruthCorrelation_kSample%d_icent%d.png",kSample,icent));
    c42->SaveAs(Form("reweightFactors/m2ptTruthCorrelation_kSample%d_icent%d.root",kSample,icent));

    TCanvas* c43=  new TCanvas("c43","",500,500);
    makeEfficiencyCanvas(c43,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c43->cd(1); 
    //    ptmcGlobal->Print("ALL");
    //    hmcRawM2Pt2Global->Print("ALL");
    
    handsomeTH1(ptmcGlobal_jzNorm,1);
    handsomeTH1(ptdataGlobal,2);
    ptmcGlobal_jzNorm->SetAxisRange(1e-6,1e8,"Y");
    ptmcGlobal_jzNorm->SetYTitle("dN/dp_{T}^{Reco}");
    ptmcGlobal_jzNorm->Draw();
    ptdataGlobal->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.6,0.86,1,24);

    //    std::cout << "ptmcGlobal_jzNorm entries " << ptmcGlobal_jzNorm->GetEntries() << " ptdataGlobal entries " << ptdataGlobal->GetEntries() << std::endl;
    TLegend *leg43 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

    easyLeg(leg43,"Mass integrated");
    leg43->AddEntry(ptdataGlobal, "Data","pl");
    leg43->AddEntry(ptmcGlobal_jzNorm, "MC Reco jzNorm","pl");
    leg43->Draw();
    ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

    c43->cd(2);
    TH1D* hptRatioGlobal_jzNorm = (TH1D*)ptdataGlobal->Clone("hptRatioGlobal_jzNorm");
    hptRatioGlobal_jzNorm->Divide(ptmcGlobal_jzNorm);
    if ( kSample == 0) hptRatioGlobal_jzNorm->SetAxisRange(0,4e4,"Y");

    hptRatioGlobal_jzNorm->SetXTitle("p_{T} (GeV)");
    hptRatioGlobal_jzNorm->SetYTitle("Data/MC");
    fixedFontHist(hptRatioGlobal_jzNorm,3,2,20);
    hptRatioGlobal_jzNorm->SetNdivisions(505,"X");
    hptRatioGlobal_jzNorm->Draw();
  
    c43->SaveAs(Form("reweightFactors/Plots/pTreweighting_jzNorm_kSample%d_icent%d.png",kSample,icent));
    c43->SaveAs(Form("reweightFactors/pTreweighting_jzNorm_kSample%d_icent%d.root",kSample,icent));

    for (int i = 0; i < nXbins; ++i) {
      std::cout << " bin " << i << " lower bin edge " << ptmcGlobal->GetBinLowEdge(i+1) << " upper bin edge " << ptmcGlobal->GetBinLowEdge(i+1) + ptmcGlobal->GetBinWidth(i+1) << std::endl;
    }
    
    TCanvas* c44=  new TCanvas("c44","",500,500);
    makeEfficiencyCanvas(c44,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c44->cd(1); 
    //    ptmcGlobal_fcalWeight->Print("ALL");
    //    hmcRawM2Pt2Global->Print("ALL");
    
    handsomeTH1(ptmcGlobal_fcalWeight,1);
    handsomeTH1(ptdataGlobal,2);
    ptmcGlobal_fcalWeight->SetAxisRange(1e-6,1e8,"Y");
    ptmcGlobal_fcalWeight->SetYTitle("dN/dp_{T}^{Reco}");
    ptmcGlobal_fcalWeight->Draw();
    ptdataGlobal->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.6,0.86,1,24);

    //    std::cout << "ptmcGlobal_fcalWeight entries " << ptmcGlobal_fcalWeight->GetEntries() << " ptdataGlobal entries " << ptdataGlobal->GetEntries() << std::endl;
    TLegend *leg44 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

    easyLeg(leg44,"Mass integrated");
    leg44->AddEntry(ptdataGlobal, "Data","pl");
    leg44->AddEntry(ptmcGlobal_fcalWeight, "MC Reco fcalWeight","pl");
    leg44->Draw();
    ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

    c44->cd(2);
    TH1D* hptRatioGlobal_fcalWeight = (TH1D*)ptdataGlobal->Clone("hptRatioGlobal_fcalWeight");
    hptRatioGlobal_fcalWeight->Divide(ptmcGlobal_fcalWeight);
    if ( kSample == 0) hptRatioGlobal_fcalWeight->SetAxisRange(0,4e4,"Y");

    hptRatioGlobal_fcalWeight->SetXTitle("p_{T} (GeV)");
    hptRatioGlobal_fcalWeight->SetYTitle("Data/MC");
    fixedFontHist(hptRatioGlobal_fcalWeight,3,2,20);
    hptRatioGlobal_fcalWeight->SetNdivisions(505,"X");
    hptRatioGlobal_fcalWeight->Draw();
  
    c44->SaveAs(Form("reweightFactors/Plots/pTreweighting_fcalWeight_kSample%d_icent%d.png",kSample,icent));
    c44->SaveAs(Form("reweightFactors/pTreweighting_fcalWeight_kSample%d_icent%d.root",kSample,icent));

    TCanvas* c45=  new TCanvas("c45","",500,500);
    makeEfficiencyCanvas(c45,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c45->cd(1); 
    //    ptmcGlobal->Print("ALL");
    //    hmcRawM2Pt2Global->Print("ALL");
    
    handsomeTH1(ptmcGlobal_jetWeight,1);
    handsomeTH1(ptdataGlobal,2);
    ptmcGlobal_jetWeight->SetAxisRange(1e1,1e13,"Y");
    ptmcGlobal_jetWeight->SetYTitle("dN/dp_{T}^{Reco}");
    ptmcGlobal_jetWeight->Draw();
    ptdataGlobal->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.6,0.86,1,24);

    //    std::cout << "ptmcGlobal_jetWeight entries " << ptmcGlobal_jetWeight->GetEntries() << " ptdataGlobal entries " << ptdataGlobal->GetEntries() << std::endl;
    TLegend *leg45 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

    easyLeg(leg45,"Mass integrated");
    leg45->AddEntry(ptdataGlobal, "Data","pl");
    leg45->AddEntry(ptmcGlobal_jetWeight, "MC Reco jet weight","pl");
    leg45->Draw();
    ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

    c45->cd(2);
    TH1D* hptRatioGlobal_jetWeight = (TH1D*)ptdataGlobal->Clone("hptRatioGlobal_jetWeight");
    hptRatioGlobal_jetWeight->Divide(ptmcGlobal_jetWeight);
    if ( kSample == 0) hptRatioGlobal_jetWeight->SetAxisRange(0,4e4,"Y");

    hptRatioGlobal_jetWeight->SetXTitle("p_{T} (GeV)");
    hptRatioGlobal_jetWeight->SetYTitle("Data/MC");
    fixedFontHist(hptRatioGlobal_jetWeight,3,2,20);
    hptRatioGlobal_jetWeight->SetNdivisions(505,"X");
    hptRatioGlobal_jetWeight->Draw();
  
    c45->SaveAs(Form("reweightFactors/Plots/pTreweighting_jetWeight_kSample%d_icent%d.png",kSample,icent));
    c45->SaveAs(Form("reweightFactors/pTreweighting_jetWeight_kSample%d_icent%d.root",kSample,icent));

    TCanvas* c46=  new TCanvas("c46","",500,500);
    makeEfficiencyCanvas(c46,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c46->cd(1); 
    //    ptmcGlobalnoWeight->Print("ALL");
    //    hmcRawM2Pt2Global->Print("ALL");
    
    handsomeTH1(ptmcGlobal_noWeight,1);
    handsomeTH1(ptdataGlobal,2);
    ptmcGlobal_noWeight->SetAxisRange(1e-6,1e8,"Y");
    ptmcGlobal_noWeight->SetYTitle("dN/dp_{T}^{Reco}");
    ptmcGlobal_noWeight->Draw();
    ptdataGlobal->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.6,0.86,1,24);

   std::cout << "ptmcGlobal_noWeight entries " << ptmcGlobal_noWeight->GetEntries() << " ptdataGlobal entries " << ptdataGlobal->GetEntries() << std::endl;
    TLegend *leg46 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

    easyLeg(leg46,"Mass integrated");
    leg46->AddEntry(ptdataGlobal, "Data","pl");
    leg46->AddEntry(ptmcGlobal_noWeight, "MC Reco","pl");
    leg46->Draw();
    ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

    c46->cd(2);
    TH1D* hptRatioGlobal_noWeight = (TH1D*)ptdataGlobal->Clone("hptRatioGlobal_noWeight");
    hptRatioGlobal_noWeight->Divide(ptmcGlobal_noWeight);
    if ( kSample == 0) hptRatioGlobal_noWeight->SetAxisRange(0,4e4,"Y");

    hptRatioGlobal_noWeight->SetXTitle("p_{T} (GeV)");
    hptRatioGlobal_noWeight->SetYTitle("Data/MC");
    fixedFontHist(hptRatioGlobal_noWeight,3,2,20);
    hptRatioGlobal_noWeight->SetNdivisions(505,"X");
    hptRatioGlobal_noWeight->Draw();
  
    c46->SaveAs(Form("reweightFactors/Plots/pTreweighting_noWeight_kSample%d_icent%d.png",kSample,icent));
    c46->SaveAs(Form("reweightFactors/pTreweighting_noWeight_kSample%d_icent%d.root",kSample,icent));

    
    TCanvas* c47=  new TCanvas("c47","",500,500);
    makeEfficiencyCanvas(c47,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c47->cd(1); 
    //    ptmcGlobal->Print("ALL");
    //    hmcRawM2Pt2Global->Print("ALL");
    
    handsomeTH1(ptmcGlobalCorr,1);
    handsomeTH1(ptdataGlobal,2);
    ptmcGlobalCorr->SetAxisRange(1e-6,1e8,"Y");
    ptmcGlobalCorr->SetYTitle("dN/dp_{T}^{Reco}");
    ptmcGlobalCorr->Draw();
    ptdataGlobal->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample, icent, 0.6,0.86,1,24);

    //    std::cout << "ptmcGlobalCorr entries " << ptmcGlobalCorr->GetEntries() << " ptdataGlobal entries " << ptdataGlobal->GetEntries() << std::endl;
    TLegend *leg47 = new TLegend(0.60,0.5583593,1,0.8023161,NULL,"brNDC");

    easyLeg(leg47,"Mass integrated");
    leg47->AddEntry(ptdataGlobal, "Data","pl");
    leg47->AddEntry(ptmcGlobalCorr, "MC Reco Corr","pl");
    leg47->Draw();
    ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

    c47->cd(2);
    TH1D* hptRatioGlobalCorr = (TH1D*)ptdataGlobal->Clone("hptRatioGlobalCorr");
    hptRatioGlobalCorr->Divide(ptmcGlobalCorr);
    if ( kSample == 0) hptRatioGlobalCorr->SetAxisRange(0,4e4,"Y");

    hptRatioGlobalCorr->SetXTitle("p_{T} (GeV)");
    hptRatioGlobalCorr->SetYTitle("Data/MC");
    fixedFontHist(hptRatioGlobalCorr,3,2,20);
    hptRatioGlobalCorr->SetNdivisions(505,"X");
    hptRatioGlobalCorr->Draw();
  
    c47->SaveAs(Form("reweightFactors/Plots/pTreweightingCorr_kSample%d_icent%d.png",kSample,icent));
    c47->SaveAs(Form("reweightFactors/pTreweightingCorr_kSample%d_icent%d.root",kSample,icent));



    
    delete c35; delete c36; delete c37; delete c38; delete c39; delete c40; delete c41; delete c42; delete c43; delete c44; delete c45; delete c46; delete c47;
    
  }

    
  /*
  
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
  
  delete c1; delete c2; delete c3; delete c16; delete c15; delete c1ptCorr; delete c111;
  delete hTemp;
}


void getMCspectra(int kSample, int icent, int etaBin,  TH2D* hmcRaw, TH2D* hmcTruth, TH2D* hm2pt, TH2D* hm2ptTruth, TH2D* hmcRaw_jzNorm, TH2D* hmcRaw_fcalWeight, TH2D* hmcRaw_jetWeight, TH2D* hmcRaw_noWeight, TF1* ptScale, int nSys) {
  
  TRandom3 genRandom;
  genRandom.SetSeed(200);

  TH1::SetDefaultSumw2();
  
  hmcRaw->Reset();
  hmcRaw_jzNorm->Reset();
  hmcRaw_fcalWeight->Reset();
  hmcRaw_jetWeight->Reset();
  hmcRaw_noWeight->Reset();
  
  hmcTruth->Reset();
  hm2pt->Reset();
  hm2ptTruth->Reset();

  TF1* fjmscal[30];

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
      
      genX = tr.jets_genPt;
      genY = tr.jets_genMass2 /  (tr.jets_genPt * tr.jets_genPt);
      recoX = tr.jets_rawPt;
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
      
      hmcRaw->Fill( recoX, recoY, tr.jets_weight * jzNorm * fcalWeight * ptWeight);
      hmcRaw_jzNorm->Fill( recoX, recoY, jzNorm);
      hmcRaw_fcalWeight->Fill( recoX, recoY, fcalWeight);
      hmcRaw_jetWeight->Fill( recoX, recoY, tr.jets_weight);
      hmcRaw_noWeight->Fill( recoX, recoY);
      hmcTruth->Fill( genX, genY, tr.jets_weight * jzNorm * fcalWeight * ptWeight);
      if (ptScale == 0) {
	hm2pt->Fill( tr.jets_massRaw2, tr.jets_rawPt, tr.jets_weight * jzNorm * fcalWeight * ptWeight);
	hm2ptTruth->Fill( tr.jets_genMass2, tr.jets_genPt, tr.jets_weight * jzNorm * fcalWeight * ptWeight);
	//	std::cout << " massRaw2 " << tr.jets_massRaw2 << " rawPt " << tr.jets_rawPt << std::endl;
      }
      ++jetCount;
    }
    cout << "jz " << ijz << " jetCount " << jetCount << endl;
  }
  //  std::cout << " hmcRawGlobal entriesss " << hmcRawGlobal->GetEntries() << std::endl;


  
  delete xBinTemp;  
}

void getDATAspectra(int kSample, int icent, int etaBin, TH2D* hdataRaw, /*TH2D* hdataRawGlobal,*/ int nSys) { 

  TH1::SetDefaultSumw2();
  hdataRaw->Reset();

  TString fname;
  if ( kSample == kPbPb ) {
    fname = pbpbDataString; 
  }
  else if ( kSample == kPP) {
    fname = ppDataString;
  }

  double eventWeight = 1.;
  if (kSample == 0) {
    // 1 / pp lumi
    eventWeight = ppEvtWgt;
  }
  else if (kSample == 1 && icent == 0) {
    // 1/ ( TAA_cent * nEvents_cent)
    eventWeight = hi9EvtWgt_0;
  }
  else if (kSample == 1 && icent == 1) {
    // 1/ ( TAA_cent * nEvents_cent)
    eventWeight = hi9EvtWgt_1;
  }
  else if (kSample == 1 && icent == 2) {
    // 1/ ( TAA_cent * nEvents_cent)
    eventWeight = hi9EvtWgt_2;
  }
  else if (kSample == 1 && icent == 3) {
    // 1/ ( TAA_cent * nEvents_cent)
    eventWeight = hi9EvtWgt_3;
  }
  else if (kSample == 1 && icent == 4) {
    // 1/ ( TAA_cent * nEvents_cent)
    eventWeight = hi9EvtWgt_4;
  }
  else if (kSample == 1 && icent == 5) {
    // 1/ ( TAA_cent * nEvents_cent)
    eventWeight = hi9EvtWgt_5;
  }
  else if (kSample == 1 && icent == 6) {
    // 1/ ( TAA_cent * nEvents_cent)
    eventWeight = hi9EvtWgt_6;
  }
  TF1* fjmscal[30];

  int nXbinsCal;
  double xBinCal[30];
  getXbin(nXbinsCal, xBinCal, 1);
  TH1D* xBinTemp = new TH1D("xBinTemp","", nXbinsCal, xBinCal);
  
  
  TFile* fData = new TFile(Form("../ntuples/%s",fname.Data()));
  TTree* tree = (TTree*)fData->Get("tr");
  tashka tr;
  tr.Init(tree);
  
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
   
    //    cout << " data massRaw2 " << tr.jets_massRaw2 << " rawPt " << tr.jets_rawPt << " recoY " << recoY << endl;
    //    hdataRaw->Fill ( recoX, recoY, eventWeight);
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
  
  getEtaPhiMap();
  
  return 0;
}  // Main program when run stand-alone
#endif
