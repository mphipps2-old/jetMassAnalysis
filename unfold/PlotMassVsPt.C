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
#include <TProfile.h>
#include <TPaletteAxis.h>
#include "unfoldingUtil.h"
#include "systematicsTool.h"

// fraction of file you want to process ie) 0.5 would process half the file
double statFrac = 1.;
// set this to either kPP or kPbPb before running
int kSample = kPP;

int main(int argc, char *argv[]) {
  

  //   TH2D *h2_jetMass2Pt = new TH2D("jetMass2Pt","jetMass2Pt",250,100,500,250,-100,2500);
  TH2D *h2_jetMass2Pt = new TH2D("jetMass2Pt","jetMass2Pt",250,100,1000,250,-100,22000);
  TH1D *h1_jetMass2OverPt2 = new TH1D("jetMass2OverPt2","jetMass2OverPt2",125,-0.01,.1);
  TH2F *h2_jetMass2OverPt2VsPt2 = new TH2F("jetMass2OverPt2VsPt2","jetMass2Pt2",250,10000,1000000,125,-.005,.1);
  TH2F *h2_jetMass2OverPt2VsPt = new TH2F("jetMass2OverPt2VsPt","jetMass2OverPt2VsPt",250,100,1000,125,-.005,.1);
  TH2F *h2_jetMass2Pt_mc = new TH2F("jetMass2Pt_mc","jetMass2Pt_mc",250,100,1000,250,-100,22000);
  TH2F *h2_jetMass2OverPt2VsPt2_mc = new TH2F("jetMass2OverPt2VsPt2_mc","jetMass2Pt2_mc",250,10000,1000000,125,-.005,.1);
  TH2F *h2_jetMass2OverPt2VsPt_mc = new TH2F("jetMass2OverPt2VsPt_mc","jetMass2OverPt2VsPt_mc",250,100,1000,125,-.005,.1);
  TH2F *h2_jetMass2Pt_truth = new TH2F("jetMass2Pt_truth","jetMass2Pt_truth",250,100,1000,250,-100,22000);
  TH2F *h2_jetMass2OverPt2VsPt2_truth = new TH2F("jetMass2OverPt2VsPt2_truth","jetMass2Pt2_truth",250,10000,1000000,125,-.005,.1);
  TH2F *h2_jetMass2OverPt2VsPt_truth = new TH2F("jetMass2OverPt2VsPt_truth","jetMass2OverPt2VsPt_truth",250,100,1000,125,-.005,.1);

   

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
  
  
  cout << " data entries " << entries << endl;
  for (Int_t i= 0; i<entries ; i++) {
    tr.GetEntry(i);
    //  way to cut processing off early ie if statFrac == 0.5, process half the file
    if ( i > entries * statFrac) break;
    // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
    //    if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_cent, tr.jets_etaCalib, tr.jets_phiCalib, false) ) // isMC = false
    // NOTE: CHANGE BACK WITH NEW TTREE!!!
    if ( ! passEvent(tr.jets_cent, tr.jets_genMass2, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, false) ) // isMC = false
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
    h2_jetMass2OverPt2VsPt2->Fill(tr.jets_ptCalib*tr.jets_ptCalib,tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib), eventWeight);
    //       h2_jetMass2OverPt2VsPt2->Fill(tr.jets_ptCalib*tr.jets_ptCalib,tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib));
    h2_jetMass2OverPt2VsPt->Fill(tr.jets_ptCalib, tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib), eventWeight);
    h1_jetMass2OverPt2->Fill(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib), eventWeight);
    //    h2_jetMass2OverPt2VsPt->Fill(tr.jets_ptCalib, tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib));
       h2_jetMass2Pt->Fill(tr.jets_ptCalib, tr.jets_mass2Calib, eventWeight);
       //     h2_jetMass2Pt->Fill(tr.jets_ptCalib, tr.jets_mass2Calib);
       //       cout << "eventWeight " << eventWeight << endl;
  }
  
  TH1D* hFcalReweight;
  if ( kSample == kPbPb ) {
    TFile* fcal = new TFile("reweightFactors/FCal_HP_v_MB_weights.root");
    hFcalReweight = (TH1D*)fcal->Get("weight_MBov_to_HP");
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
      tr.GetEntry(i);
      // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
      //      if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_cent, tr.jets_etaCalib, tr.jets_phiCalib, true) ) // isMC = true
      if ( ! passEvent(tr.jets_cent, tr.jets_genMass2, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, true) ) // isMC = true
	continue;
            
      double fcalWeight;
      if ( kSample==kPbPb) {
	fcalWeight = hFcalReweight->GetBinContent(hFcalReweight->GetXaxis()->FindBin(tr.jets_chMassRcSubt));
      }
      else fcalWeight = 1.0;

      //      cout << " pt " << tr.jets_ptCalib << " jetweight " << tr.jets_weight << " jznorm " << jzNorm << " fcal " << fcalWeight << " total weight " << tr.jets_weight * jzNorm * fcalWeight << endl;

      h2_jetMass2OverPt2VsPt2_mc->Fill(tr.jets_ptCalib*tr.jets_ptCalib,tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib),tr.jets_weight * jzNorm * fcalWeight);
      h2_jetMass2OverPt2VsPt_mc->Fill(tr.jets_ptCalib, tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib),tr.jets_weight * jzNorm * fcalWeight);
      h2_jetMass2Pt_mc->Fill(tr.jets_ptCalib, tr.jets_mass2Calib,tr.jets_weight * jzNorm * fcalWeight);
      h2_jetMass2OverPt2VsPt2_truth->Fill(tr.jets_genPt*tr.jets_genPt,tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt),tr.jets_weight * jzNorm * fcalWeight);
      h2_jetMass2OverPt2VsPt_truth->Fill(tr.jets_genPt, tr.jets_genMass2/(tr.jets_genPt*tr.jets_genPt),tr.jets_weight * jzNorm * fcalWeight);
      h2_jetMass2Pt_truth->Fill(tr.jets_genPt, tr.jets_genMass2,tr.jets_weight * jzNorm * fcalWeight);      
      ++jetCount;
    }   
  }

  
   TString outFileName;
   gStyle->SetPalette(kBird);
   TString fileName; 

   TCanvas *c1 = new TCanvas("c1","c1");
   cout << " nJets data " << h2_jetMass2OverPt2VsPt2->GetEntries() << endl;
   cout << " nJets mc " << h2_jetMass2OverPt2VsPt2_mc->GetEntries() << endl;
   
   h2_jetMass2OverPt2VsPt2->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");      
   h2_jetMass2OverPt2VsPt2->GetXaxis()->SetTitle("#it{p}_{T}^{2}");
   h2_jetMass2OverPt2VsPt2->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2OverPt2VsPt2->SetContour(99);
   h2_jetMass2OverPt2VsPt2->Draw("colz");
   c1->SetLogz();
   h2_jetMass2OverPt2VsPt2->SetMinimum(1e-5);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppData/mass2OverPt2VsPt2_ppData.png");
   c1->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppData/mass2OverPt2VsPt2_ppData.root");
   c1->SaveAs(outFileName.Data());


   TCanvas *c1b = new TCanvas("c1b","c1b");
   TProfile *p1 = h2_jetMass2OverPt2VsPt2->ProfileX("px1");
   p1->Draw();
   p1->SetMarkerSize(0.25);
   outFileName.Form("../ntuples/plots/ppData/mass2OverPt2VsPt2_ppData_prof.png");
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   c1b->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppData/mass2OverPt2VsPt2_ppData_prof.root");
   c1b->SaveAs(outFileName.Data());

   TCanvas *c1c = new TCanvas("c1c","c1c");
   
   h1_jetMass2OverPt2->GetXaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");      
   h1_jetMass2OverPt2->GetXaxis()->SetTitleOffset(1.55);
   h1_jetMass2OverPt2->Draw();
   //   c1->SetLogz();
   h1_jetMass2OverPt2->SetMinimum(1e-5);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppData/mass2OverPt2_ppData.png");
   c1c->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppData/mass2OverPt2_ppData.root");
   c1c->SaveAs(outFileName.Data());


   
   TCanvas *c2 = new TCanvas("c2","c2");
   h2_jetMass2Pt->GetYaxis()->SetTitle("m^{2}");      
   h2_jetMass2Pt->GetXaxis()->SetTitle("#it{p}_{T}");
   h2_jetMass2Pt->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2Pt->SetContour(99);
   h2_jetMass2Pt->Draw("colz");
   //   h2_jetMass2Pt->GetZaxis()->SetRangeUser(5e-4,2e1);
   c2->SetLogz();
   h2_jetMass2Pt->SetMinimum(1e-5);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppData/mass2Pt_ppData.png");
   c2->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppData/mass2Pt_ppData.root");
   c2->SaveAs(outFileName.Data());

   TCanvas *c2b = new TCanvas("c2b","c2b");
   TProfile *p2 = h2_jetMass2Pt->ProfileX("px2");
   p2->Draw();
   p2->SetMarkerSize(0.25);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppData/mass2Pt_ppData_prof.png");
   c2b->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppData/mass2Pt_ppData_prof.root");
   c2b->SaveAs(outFileName.Data());

   TCanvas *c3 = new TCanvas("c3","c3");
   h2_jetMass2OverPt2VsPt->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");      
   h2_jetMass2OverPt2VsPt->GetXaxis()->SetTitle("#it{p}_{T}");
   h2_jetMass2OverPt2VsPt->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2OverPt2VsPt->SetContour(99);
   h2_jetMass2OverPt2VsPt->Draw("colz");
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   c3->SetLogz();
   h2_jetMass2OverPt2VsPt->SetMinimum(1e-5);
   outFileName.Form("../ntuples/plots/ppData/mass2OverPt2VsPt_ppData.png");
   c3->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppData/mass2OverPt2VsPt_ppData.root");
   c3->SaveAs(outFileName.Data());

   TCanvas *c3b = new TCanvas("c3b","c3b");
   TProfile *p3 = h2_jetMass2OverPt2VsPt->ProfileX("px3");
   p3->Draw();
   p3->SetMarkerSize(0.25);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppData/mass2OverPt2VsPt_ppData_prof.png");
   c3b->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppData/mass2OverPt2VsPt_ppData_prof.root");
   c3b->SaveAs(outFileName.Data());

   TCanvas *c4 = new TCanvas("c4","c4");
   h2_jetMass2OverPt2VsPt2_mc->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");      
   h2_jetMass2OverPt2VsPt2_mc->GetXaxis()->SetTitle("#it{p}_{T}^{2}");
   h2_jetMass2OverPt2VsPt2_mc->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2OverPt2VsPt2_mc->SetContour(99);
   h2_jetMass2OverPt2VsPt2_mc->Draw("colz");
   c4->SetLogz();
   h2_jetMass2OverPt2VsPt2_mc->SetMinimum(1e-9);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppMC/mass2OverPt2VsPt2_ppMC.png");
   c4->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppMC/mass2OverPt2VsPt2_ppMC.root");
   c4->SaveAs(outFileName.Data());

   TCanvas *c4b = new TCanvas("c4b","c4b");
   TProfile *p4 = h2_jetMass2OverPt2VsPt2_mc->ProfileX("px4");
   p4->Draw();
   p4->SetMarkerSize(0.25);
   outFileName.Form("../ntuples/plots/ppMC/mass2OverPt2VsPt2_ppMC_prof.png");
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   c4b->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppMC/mass2OverPt2VsPt2_ppMC_prof.root");
   c4b->SaveAs(outFileName.Data());
   
   TCanvas *c5 = new TCanvas("c5","c5");
   h2_jetMass2Pt_mc->GetYaxis()->SetTitle("m^{2}");      
   h2_jetMass2Pt_mc->GetXaxis()->SetTitle("#it{p}_{T}");
   h2_jetMass2Pt_mc->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2Pt_mc->SetContour(99);
   h2_jetMass2Pt_mc->Draw("colz");
   c5->SetLogz();
   h2_jetMass2Pt_mc->SetMinimum(1e-9);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppMC/mass2Pt_ppMC.png");
   c5->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppMC/mass2Pt_ppMC.root");
   c5->SaveAs(outFileName.Data());

   TCanvas *c5b = new TCanvas("c5b","c5b");
   TProfile *p5 = h2_jetMass2Pt_mc->ProfileX("px5");
   p5->Draw();
   p5->SetMarkerSize(0.25);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppMC/mass2Pt_ppMC_prof.png");
   c5b->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppMC/mass2Pt_ppMC_prof.root");
   c5b->SaveAs(outFileName.Data());

   TCanvas *c6 = new TCanvas("c6","c6");
   h2_jetMass2OverPt2VsPt_mc->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");      
   h2_jetMass2OverPt2VsPt_mc->GetXaxis()->SetTitle("#it{p}_{T}");
   h2_jetMass2OverPt2VsPt_mc->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2OverPt2VsPt_mc->SetContour(99);
   h2_jetMass2OverPt2VsPt_mc->Draw("colz");
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   c6->SetLogz();
   h2_jetMass2OverPt2VsPt_mc->SetMinimum(1e-9);
   outFileName.Form("../ntuples/plots/ppMC/mass2OverPt2VsPt_ppMC.png");
   c6->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppMC/mass2OverPt2VsPt_ppMC.root");
   c6->SaveAs(outFileName.Data());

   TCanvas *c6b = new TCanvas("c6b","c6b");
   TProfile *p6 = h2_jetMass2OverPt2VsPt_mc->ProfileX("px6");
   p6->Draw();
   p6->SetMarkerSize(0.25);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppMC/mass2OverPt2VsPt_ppMC_prof.png");
   c6b->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppMC/mass2OverPt2VsPt_ppMC_prof.root");
   c6b->SaveAs(outFileName.Data());

   
   TCanvas *c7 = new TCanvas("c7","c7");
   h2_jetMass2OverPt2VsPt2_truth->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");
   h2_jetMass2OverPt2VsPt2_truth->GetXaxis()->SetTitle("#it{p}_{T}^{2}");
   h2_jetMass2OverPt2VsPt2_truth->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2OverPt2VsPt2_truth->SetContour(99);
   h2_jetMass2OverPt2VsPt2_truth->Draw("colz");
   c7->SetLogz();
   h2_jetMass2OverPt2VsPt2_truth->SetMinimum(1e-9);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppTruth/mass2OverPt2VsPt2_ppTruth.png");
   c7->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppTruth/mass2OverPt2VsPt2_ppTruth.root");
   c7->SaveAs(outFileName.Data());

   TCanvas *c7b = new TCanvas("c7b","c7b");
   TProfile *p7 = h2_jetMass2OverPt2VsPt2_truth->ProfileX("px4");
   p7->Draw();
   p7->SetMarkerSize(0.25);
   outFileName.Form("../ntuples/plots/ppTruth/mass2OverPt2VsPt2_ppTruth_prof.png");
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   c7b->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppTruth/mass2OverPt2VsPt2_ppTruth_prof.root");
   c7b->SaveAs(outFileName.Data());
   
   TCanvas *c8 = new TCanvas("c8","c8");
   h2_jetMass2Pt_truth->GetYaxis()->SetTitle("m^{2}");      
   h2_jetMass2Pt_truth->GetXaxis()->SetTitle("#it{p}_{T}");
   h2_jetMass2Pt_truth->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2Pt_truth->SetContour(99);
   h2_jetMass2Pt_truth->Draw("colz");
   c8->SetLogz();
   h2_jetMass2Pt_truth->SetMinimum(1e-9);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppTruth/mass2Pt_ppTruth.png");
   c8->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppTruth/mass2Pt_ppTruth.root");
   c8->SaveAs(outFileName.Data());

   TCanvas *c8b = new TCanvas("c8b","c8b");
   TProfile *p8 = h2_jetMass2Pt_truth->ProfileX("px5");
   p8->Draw();
   p8->SetMarkerSize(0.25);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppTruth/mass2Pt_ppTruth_prof.png");
   c8b->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppTruth/mass2Pt_ppTruth_prof.root");
   c8b->SaveAs(outFileName.Data());

   TCanvas *c9 = new TCanvas("c9","c9");
   h2_jetMass2OverPt2VsPt_truth->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");      
   h2_jetMass2OverPt2VsPt_truth->GetXaxis()->SetTitle("#it{p}_{T}");
   h2_jetMass2OverPt2VsPt_truth->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2OverPt2VsPt_truth->SetContour(99);
   h2_jetMass2OverPt2VsPt_truth->Draw("colz");
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   c9->SetLogz();
   h2_jetMass2OverPt2VsPt_truth->SetMinimum(1e-9);
   outFileName.Form("../ntuples/plots/ppTruth/mass2OverPt2VsPt_ppTruth.png");
   c9->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppTruth/mass2OverPt2VsPt_ppTruth.root");
   c9->SaveAs(outFileName.Data());

   TCanvas *c9b = new TCanvas("c9b","c9b");
   TProfile *p9 = h2_jetMass2OverPt2VsPt_truth->ProfileX("px6");
   p9->Draw();
   p9->SetMarkerSize(0.25);
   drawCentrality(kSample,0,0.70,0.86,1,24);
   ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
   outFileName.Form("../ntuples/plots/ppTruth/mass2OverPt2VsPt_ppTruth_prof.png");
   c9b->SaveAs(outFileName.Data());
   outFileName.Form("../ntuples/plots/ppTruth/mass2OverPt2VsPt_ppTruth_prof.root");
   c9b->SaveAs(outFileName.Data());

   delete c1; delete c1b;
   delete c2; delete c2b;
   delete c3; delete c3b;
   delete c4; delete c4b;
   delete c5; delete c5b;
   delete c6; delete c6b;
   delete c7; delete c7b;
   delete c8; delete c8b;
   delete c9; delete c9b;
    
  return 0;
}
