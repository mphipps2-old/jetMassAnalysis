// Goal: produce ratio plot of Data/MC for dsigma/dpt vs pt. Eta binning integrated from -2.1 < eta < 2.1. Works for pp or PbPb by flipping kSample switch
// KEEP THESE SCRIPTS SMALL AND SINGLE PURPOSED!

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
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
  int nCentIntervals;
  if (kSample == kPP) nCentIntervals = 1;
  else nCentIntervals = 7; 
  
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

  int nXbins;
  double xBin[30];
  double binWidth[30];
  getXbin(nXbins, xBin, 77);
  cout << " nXbins = " << nXbins << endl;
  for (int i = 0; i < nXbins; ++i) {
    cout << " xBin " << i << " = " <<  xBin[i] <<endl;
    binWidth[i] = xBin[i+1]-xBin[i];
  }
  
  TH1D* hCrossSec_data[7];
  TH1D* hCrossSec_mc[7];
  TH1D* hCrossSec_mcTruth[7];
  TH1D* hPtWidth[7];
  
  int cent = 0;
  long long entries;
  TString fname;
  if ( kSample == kPbPb ) {
    fname = pbpbDataString; 
  }
  else if ( kSample == kPP) {
    fname = ppDataString;
  }
  
  double eventWeight = 1.;
  TFile* fData = new TFile(Form("../ntuples/%s",fname.Data()));
  TTree* tree = (TTree*)fData->Get("tr");
  cout << " Setting Tashka for data ..." << endl;
  tashka tr;
  tr.Init(tree);
  entries = tree->GetEntries();
  
  for (int i=0; i < nCentIntervals; ++i) {
    char histName[256];
    sprintf(histName,"hCrossSec_data_cent%d",i);
    hCrossSec_data[i]  = new TH1D(histName,histName,nXbins,xBin);
  }
  
  cout << " data entries " << entries << endl;
  for (Int_t i= 0; i<entries ; i++) {
    tr.GetEntry(i);
    //  way to cut processing off early ie if statFrac == 0.5, process half the file
    if ( i > entries * statFrac) break;
    // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
    if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, false) ) // isMC = false
    // NOTE: CHANGE BACK WITH NEW TTREE!!!
    //if ( ! passEvent(tr.jets_cent, tr.jets_genMass2, tr.jets_ptCalib, tr.jets_centtr.jets_etaCalib, tr.jets_phiCalib, false) ) // isMC = false
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
    hCrossSec_data[cent]->Fill ( tr.jets_ptCalib, eventWeight);
  }
  for (int i=0; i < nCentIntervals; ++i) {
    cout << "cent " << i << " data events that passed cuts " << hCrossSec_data[i]->GetEntries() << endl;
  }
  TH1D* hFcalReweight;
  if ( kSample == kPbPb ) {
    TFile* fcal = new TFile("reweightFactors/FCal_HP_v_MB_weights.root");
    hFcalReweight = (TH1D*)fcal->Get("weight_MBov_to_HP");
  }

  for (int i=0; i < nCentIntervals; ++i) {
    char histName[256];
    sprintf(histName,"hCrossSec_mc_cent%d",i);
    hCrossSec_mc[i] = new TH1D(histName,histName,nXbins,xBin);
    sprintf(histName,"hCrossSec_mcTruth_cent%d",i);
    hCrossSec_mcTruth[i] = new TH1D(histName,histName,nXbins,xBin);
    hPtWidth[i] = (TH1D*) hCrossSec_mc[i]->Clone(Form("hPtWeight_cent%d",i));
  }
  
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
    cout << "Scanning JZ"<<ijz<<" file.  Total events = " << entries << endl;
    int jetCount = 0;
    for (Int_t i= 0; i<entries ; i++) {
      //  way to cut processing off early ie if statFrac == 0.5, process half the file
      if ( i > entries * statFrac ) break;
      
      // why was this here!!!???
      //      if (i%2==0)  continue;
      
      tr.GetEntry(i);
      // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
       if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, true) ) // isMC = true
	 //    if ( ! passEvent(tr.jets_cent, tr.jets_genMass2, tr.jets_ptCalib, tr.jets_cent, tr.jets_etaCalib, tr.jets_phiCalib, true) ) // isMC = true
	continue;
            
      double fcalWeight;
      if ( kSample==kPbPb) {

	//	cout << " jet fcal " << tr.jets_fcal << " rcalreweight " << hFcalReweight->GetXaxis()->FindBin(tr.jets_fcal) << endl;
	//	cout << " jet fcal " << tr.jets_chMassRcSubt << " rcalreweight " << hFcalReweight->GetXaxis()->FindBin(tr.jets_chMassRcSubt) << endl;
	fcalWeight = hFcalReweight->GetBinContent(hFcalReweight->GetXaxis()->FindBin(tr.jets_fcal));
	
      }
      else fcalWeight = 1.0;
      cent = tr.jets_cent;
      cout << " pt " << tr.jets_ptCalib << " jetweight " << tr.jets_weight << " jznorm " << jzNorm << " fcal " << fcalWeight << " total weight " << tr.jets_weight * jzNorm * fcalWeight << endl;
      hCrossSec_mc[cent]->Fill( tr.jets_ptCalib, tr.jets_weight * jzNorm * fcalWeight);
      hCrossSec_mcTruth[cent]->Fill( tr.jets_genPt, tr.jets_weight * jzNorm * fcalWeight);  
      ++jetCount;
    }   
  }
  for (int i=0; i < nCentIntervals; ++i) {
    cout << "cent " << i << " mc events that passed cuts " << hCrossSec_mc[i]->GetEntries() << endl;
    for (int j=0; j < hCrossSec_mc[i]->GetNbinsX(); ++j) {
      hPtWidth[i]->SetBinContent(j+1,binWidth[j]);
      cout << "pt bin " << j << " width " << binWidth[j] << endl;    }
  }
  
  
  for (int i=0; i < nCentIntervals; ++i) {
    hCrossSec_data[i]->Divide(hPtWidth[i]);
    hCrossSec_mc[i]->Divide(hPtWidth[i]);
    hCrossSec_mcTruth[i]->Divide(hPtWidth[i]);
    TCanvas* c1=  new TCanvas("c1","",500,500);
    makeEfficiencyCanvas(c1,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c1->cd(1);
    // defined in commonUtility: 1 is black; 2 is red    
    handsomeTH1(hCrossSec_mc[i],1);
    handsomeTH1(hCrossSec_data[i],2);
    fixedFontHist(hCrossSec_mc[i],2.5,2,20);
    hCrossSec_mc[i]->SetAxisRange(5e-9,5e3,"Y");
    if (kSample == kPP) hCrossSec_mc[i]->SetYTitle("#frac{d#sigma}{dp_{T}} [nb/GeV]");
    else if (kSample == kPbPb) hCrossSec_mc[i]->SetYTitle("#frac{1}{#LT T_{AA} #GT} #frac{1}{N_{evt}} #frac{d#sigma}{dp_{T}} [nb/GeV]");
    hCrossSec_mc[i]->Draw();
    hCrossSec_data[i]->Draw("same");
    cout << "mc " << hCrossSec_mc[i]->GetEntries() << endl;
    cout << "data " << hCrossSec_data[i]->GetEntries() << endl;
    gPad->SetLogy();
    drawCentrality(kSample,i,0.70,0.86,1,24);
    TLegend *leg1 = new TLegend(0.3871854,0.6394175,0.7849639,0.8533743,NULL,"brNDC");
    easyLeg(leg1,"Before p_{T} reweighting");
    leg1->AddEntry(hCrossSec_data[i], "Data","pl");
    leg1->AddEntry(hCrossSec_mc[i], "MC Reco","pl");
    leg1->Draw();
    ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    
    c1->cd(2);
    TH1D* hCrossSec_ratio = (TH1D*)hCrossSec_data[i]->Clone("hCrossSec_ratio");
    hCrossSec_ratio->Divide(hCrossSec_mc[i]);

    if (kSample == kPP) hCrossSec_ratio->SetAxisRange(.4,.85,"Y");
    else if (kSample == kPbPb && i == 0) hCrossSec_ratio->SetAxisRange(.004,.019,"Y");
    else if (kSample == kPbPb && i == 1) hCrossSec_ratio->SetAxisRange(.005,.0275,"Y");
    else if (kSample == kPbPb && i == 2) hCrossSec_ratio->SetAxisRange(.01,.065,"Y");
    else if (kSample == kPbPb && i == 3) hCrossSec_ratio->SetAxisRange(.02,.125,"Y");
    else if (kSample == kPbPb && i == 4) hCrossSec_ratio->SetAxisRange(.0,.355,"Y");
    else if (kSample == kPbPb && i == 5) hCrossSec_ratio->SetAxisRange(.0,.22,"Y");
    else if (kSample == kPbPb && i == 6) hCrossSec_ratio->SetAxisRange(.0,.17,"Y");
    hCrossSec_ratio->SetXTitle("p_{T}^{Reco}");
    hCrossSec_ratio->SetYTitle("Data/MC");
    hCrossSec_ratio->SetNdivisions(505,"X");
    fixedFontHist(hCrossSec_ratio,3,2,20);
    hCrossSec_ratio->Draw();
    //    jumSun(-.2,1,0.35,1);
  
    c1->SaveAs(Form("Plots/CrossSection/CrossSection_kSample%d_icent%d.png",kSample,i));
    c1->SaveAs(Form("Plots/CrossSection/CrossSection_kSample%d_icent%d.root",kSample,i));
    
    TCanvas* c2=  new TCanvas("c2","",500,500);
    makeEfficiencyCanvas(c2,1, 0.0, 0.01, 0.2, 0.25, 0.01);
    c2->cd(1);
    // defined in commonUtility: 1 is black; 2 is red    
    handsomeTH1(hCrossSec_mcTruth[i],1);
    fixedFontHist(hCrossSec_mcTruth[i],2.5,2,20);
    hCrossSec_mcTruth[i]->SetAxisRange(5e-9,5e3,"Y");
    if (kSample == kPP) hCrossSec_mcTruth[i]->SetYTitle("#frac{d#sigma}{dp_{T}} [nb/GeV]");
    else if (kSample == kPbPb) hCrossSec_mcTruth[i]->SetYTitle("#frac{1}{#LT T_{AA} #GT} #frac{1}{N_{evt}} #frac{d#sigma}{dp_{T}} [nb/GeV]");
    
    hCrossSec_mcTruth[i]->Draw();
    hCrossSec_data[i]->Draw("same");
    gPad->SetLogy();
    drawCentrality(kSample,i,0.70,0.86,1,24);
    TLegend *leg2 = new TLegend(0.3871854,0.6394175,0.7849639,0.8533743,NULL,"brNDC");
    easyLeg(leg2,"Before p_{T} reweighting");
    leg2->AddEntry(hCrossSec_data[i], "Data","pl");
    leg2->AddEntry(hCrossSec_mcTruth[i], "MC Truth","pl");
    leg2->Draw();
    ATLASLabel(0.22, 0.88, "Internal");//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    
    c2->cd(2);
    TH1D* hCrossSec_ratioTruth = (TH1D*)hCrossSec_data[i]->Clone("hCrossSec_ratioTruth");
    hCrossSec_ratioTruth->Divide(hCrossSec_mcTruth[i]);
    /*
        if (kSample == kPP) hCrossSec_ratioTruth->SetAxisRange(.07,.25,"Y");
    else if (kSample == kPbPb && i == 0) hCrossSec_ratioTruth->SetAxisRange(.0,.015,"Y");
    else if (kSample == kPbPb && i == 1) hCrossSec_ratioTruth->SetAxisRange(.0,.021,"Y");
    else if (kSample == kPbPb && i == 2) hCrossSec_ratioTruth->SetAxisRange(.0,.015,"Y");
    else if (kSample == kPbPb && i == 3) hCrossSec_ratioTruth->SetAxisRange(.005,.017,"Y");
    else if (kSample == kPbPb && i == 4) hCrossSec_ratioTruth->SetAxisRange(.002,.027,"Y");
    else if (kSample == kPbPb && i == 5) hCrossSec_ratioTruth->SetAxisRange(.0,.07,"Y");
    else if (kSample == kPbPb && i == 6) hCrossSec_ratioTruth->SetAxisRange(.0,.037,"Y");
    */
    
    hCrossSec_ratioTruth->SetXTitle("p_{T}^{Reco}");
    hCrossSec_ratioTruth->SetYTitle("Data/MC");
    hCrossSec_ratioTruth->SetNdivisions(505,"X");
    fixedFontHist(hCrossSec_ratioTruth,3,2,20);
    hCrossSec_ratioTruth->Draw();
    //    jumSun(-.2,1,0.35,1);
  
    c2->SaveAs(Form("Plots/CrossSection/CrossSectionTruth_kSample%d_icent%d.png",kSample,i));
    c2->SaveAs(Form("Plots/CrossSection/CrossSectionTruth_kSample%d_icent%d.root",kSample,i));
    
    delete c1;  delete c2;

    
  }  
  return 0;
}
