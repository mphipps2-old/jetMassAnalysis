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
  
  int nPtBins;
  double ptBins[30];
  int nMBins;
  double mBins[40];
  double binWidth[30];
  getXbin(nPtBins, ptBins, 77);
  cout << " nXbins = " << nPtBins << endl;
  for (int i = 0; i < nPtBins; ++i) {
    cout << " ptBin " << i << " = " <<  ptBins[i] <<endl;
    binWidth[i] = ptBins[i+1]-ptBins[i];
  }
  getYbin(nMBins, mBins, 772);
  cout << " nMBins = " << nMBins << endl;
  TH1D* hCrossSec_data[6][7];
  TH1D* hCrossSec_mc[6][7];
  TH1D* hPtWidth;
  TH1D* hCrossSec_ratio[6][7];
  TH1D* hCrossSec_m_data[6][7];
  TH1D* hCrossSec_m_mc[6][7];
  TH1D* hCrossSec_m_ratio[6][7];
  int cent = 0;  
  float eta = 0.;
  long long entries;
  TString fname;
  if ( kSample == kPbPb ) {
    fname = pbpbDataString; 
  }
  else if ( kSample == kPP) {
    fname = ppDataString;
  }
  double histoOffset;
  double eventWeight = 1.;
  TFile* fData = new TFile(Form("../ntuples/%s",fname.Data()));
  TTree* tree = (TTree*)fData->Get("tr");
  cout << " Setting Tashka for data ..." << endl;
  tashka tr;
  tr.Init(tree);
  entries = tree->GetEntries();
 
  for (int i=0; i < nEtaIntervals; ++i) { 
    for (int j=0; j < nCentIntervals; ++j) {
      char histName[256];
      sprintf(histName,"hCrossSec_data_eta%dcent%d",i,j);
      hCrossSec_data[i][j]  = new TH1D(histName,histName,nPtBins,ptBins);
      sprintf(histName,"hCrossSec_mc_eta%dcent%d",i,j);
      hCrossSec_mc[i][j]  = new TH1D(histName,histName,nPtBins,ptBins);
      sprintf(histName,"hCrossSec_m2pt2_data_eta%dcent%d",i,j);
      hCrossSec_m_data[i][j]  = new TH1D(histName,histName,nMBins,mBins);
      sprintf(histName,"hCrossSec_m2pt2_mc_eta%dcent%d",i,j);
      hCrossSec_m_mc[i][j]  = new TH1D(histName,histName,nMBins,mBins);
    }
  }

  TH1D* hFcalReweight;
  if ( kSample == kPbPb ) {
    TFile* fcal = new TFile("reweightFactors/FCal_HP_v_MB_weights.root");
    hFcalReweight = (TH1D*)fcal->Get("weight_MBov_to_HP");
  }
  TH2D* hReweightingFactors; 
  TFile *reweightFile = new TFile("ReweightingFactors/ReweightingFactors_kSample0_etaBin5.root");
  hReweightingFactors = (TH2D*) reweightFile->Get("ReweightingHisto");



  for (Int_t i= 0; i<entries ; i++) {
    tr.GetEntry(i);
    //  way to cut processing off early ie if statFrac == 0.5, process half the file
    if ( i > entries * statFrac) break;
    // cuts made here: centrality, pt window, detector hole (0 < eta < 1; pi/4 < phi < 11pi/32)
    cent = tr.jets_cent;
    eta  = TMath::Abs(tr.jets_yCalib);
    if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, false) ) // isMC = false
      //   if ( ! passEvent(tr.jets_cent, tr.jets_genMass2, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, false) ) // isMC = false      
      continue;
    int etaBin = -1;
    if (eta < 0.3) etaBin = 0;
    else if (eta >= 0.3 && eta < 0.8) etaBin = 1;
    else if (eta >= 0.8 && eta < 1.2) etaBin = 2;
    else if (eta >= 1.2 && eta < 1.6) etaBin = 3;
    else if (eta >= 1.6 && eta < 2.1) etaBin = 4;
    else continue;
    //      cout << " passed " << endl;
    if (kSample == kPP && etaBin == 0) {
      // 1 / pp lumi (nb)
      eventWeight = ppEvtWgt * (1e+8);
    }
    else if (kSample == kPP && etaBin == 1) {
      eventWeight = ppEvtWgt*(1e+6);
    }
    else if (kSample == kPP && etaBin == 2) {
      eventWeight = ppEvtWgt*(1e+4);
    }
    else if (kSample == kPP && etaBin == 3) {
      eventWeight = ppEvtWgt*(1e+2);
    }
    else if (kSample == kPP && etaBin == 4) {
      eventWeight = ppEvtWgt*1;
    }
    else if (kSample == kPbPb && cent == 0) {
      // 1/ ( TAA_cent * nEvents_cent) (nb)
      eventWeight = hi9EvtWgt_0;
    }
    else if (kSample == kPbPb && cent == 1) {
      eventWeight = hi9EvtWgt_1*(1e+2);
    }
    else if (kSample == kPbPb && cent == 2) {
      eventWeight = hi9EvtWgt_2*(1e+4);
    }
    else if (kSample == kPbPb && cent == 3) {
      eventWeight = hi9EvtWgt_3*(1e+6);
    }
    else if (kSample == kPbPb && cent == 4) {
      eventWeight = hi9EvtWgt_4*(1e+8);
    }
    else if (kSample == kPbPb && cent == 5) {
      eventWeight = hi9EvtWgt_5*(1e+10);
    }
    else if (kSample == kPbPb && cent == 6) {
      eventWeight = hi9EvtWgt_6*(1e+12);
    }
    hCrossSec_data[etaBin][cent]->Fill ( tr.jets_ptCalib, eventWeight);
    if (kSample == kPP && eta < 2.1) {
      hCrossSec_data[5][0]->Fill( tr.jets_ptCalib, ppEvtWgt*(1e+10));
      hCrossSec_m_data[5][0]->Fill( (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), ppEvtWgt*(1e+10));
    }
    
    hCrossSec_m_data[etaBin][cent]->Fill ( (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), eventWeight);
  }

  
  cout << " Fill MC " << endl;
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
      if ( ! passEvent(tr.jets_cent, tr.jets_genPt, tr.jets_ptCalib, tr.jets_etaCalib, tr.jets_phiCalib, true) ) // isMC = true
      
	continue;
            
      double fcalWeight;
      if ( kSample==kPbPb) {
	//	cout << " jet fcal " << tr.jets_chMassRcSubt << " rcalreweight " << hFcalReweight->GetXaxis()->FindBin(tr.jets_chMassRcSubt) << endl;
	fcalWeight = hFcalReweight->GetBinContent(hFcalReweight->GetXaxis()->FindBin(tr.jets_fcal));
      }
      else fcalWeight = 1.0;
      //      double reweightFactor = hReweightingFactors->GetBinContent(hReweightingFactors->GetXaxis()->FindBin(tr.jets_ptCalib),hReweightingFactors->GetYaxis()->FindBin(tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)));
      double reweightFactor = 1;
      //      if (reweightFactor == 0) continue;
      //      else reweightFactor = 1. / reweightFactor;
      //      cout << " reweightFactor " << reweightFactor << endl;
      cent = tr.jets_cent;
      eta  = TMath::Abs(tr.jets_yCalib);
      int etaBin = -1;
      if (eta < 0.3) etaBin = 0;
      else if (eta >= 0.3 && eta < 0.8) etaBin = 1;
      else if (eta >= 0.8 && eta < 1.2) etaBin = 2;
      else if (eta >= 1.2 && eta < 1.6) etaBin = 3;
      else if (eta >= 1.6 && eta < 2.1) etaBin = 4;   
      else continue;
      if (kSample == kPP && etaBin == 0) {
	// 1 / pp lumi (nb)
	histoOffset = (1e+8);
      }
      else if (kSample == kPP && etaBin == 1) {
	histoOffset = (1e+6);
      }
      else if (kSample == kPP && etaBin == 2) {
	histoOffset = (1e+4);
      }
      else if (kSample == kPP && etaBin == 3) {
	histoOffset = (1e+2);
      }
      else if (kSample == kPP && etaBin == 4) {
	histoOffset = 1;
      }
      else if (kSample == kPbPb && cent == 0) {
	// 1/ ( TAA_cent * nEvents_cent) (nb)
	histoOffset = 1;
      }
      else if (kSample == kPbPb && cent == 1) {
	histoOffset = (1e+2);
      }
      else if (kSample == kPbPb && cent == 2) {
	histoOffset = (1e+4);
      }
      else if (kSample == kPbPb && cent == 3) {
	histoOffset = (1e+6);
      }
      else if (kSample == kPbPb && cent == 4) {
	histoOffset = (1e+8);
      }
      else if (kSample == kPbPb && cent == 5) {
	histoOffset = (1e+10);
      }
     else if (kSample == kPbPb && cent == 6) {
	histoOffset = (1e+12);
      }
      //      cout << " eta " << etaBin << " cent " << cent << " pt " << tr.jets_ptCalib << " weight " << tr.jets_weight << " jznorm " << jzNorm << " fcal " << fcalWeight << " histooffset " << histoOffset << " total weight " << tr.jets_weight * jzNorm * fcalWeight * histoOffset << endl;
      hCrossSec_mc[etaBin][cent]->Fill( tr.jets_ptCalib, tr.jets_weight * jzNorm * fcalWeight * histoOffset * reweightFactor);
      hCrossSec_m_mc[etaBin][cent]->Fill((tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * fcalWeight * histoOffset * reweightFactor);
      if (kSample == kPP && eta < 2.1) {
	histoOffset = (1e+10);
	hCrossSec_mc[5][0]->Fill( tr.jets_ptCalib, tr.jets_weight * jzNorm * fcalWeight * histoOffset * reweightFactor);
	hCrossSec_m_mc[5][0]->Fill( (tr.jets_mass2Calib/(tr.jets_ptCalib*tr.jets_ptCalib)), tr.jets_weight * jzNorm * histoOffset * fcalWeight * reweightFactor);
      }
      
      ++jetCount;
    }   
  }
    

  for (int i=0; i < nEtaIntervals; ++i) {
    for (int j=0; j < nCentIntervals; ++j) {
            cout << "eta " << i << " cent " << j << " data events that passed cuts " << hCrossSec_data[i][j]->GetEntries() << endl;
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
	TH1ScaleByWidth(hCrossSec_data[i][j]);
	TH1ScaleByWidth(hCrossSec_mc[i][j]);
	TH1ScaleByRapidity(hCrossSec_data[i][j],i);
	TH1ScaleByRapidity(hCrossSec_mc[i][j],i);

	// defined in commonUtility: 1 is black; 2 is red    
	if (j == 0) {
	  handsomeTH1(hCrossSec_data[i][j],2);
	  handsomeTH1(hCrossSec_mc[i][j],2);
	  handsomeTH1(hCrossSec_m_data[i][j],2);
	  handsomeTH1(hCrossSec_m_mc[i][j],2);
	}
	else if (j == 1) {
	  handsomeTH1(hCrossSec_data[i][j],3);
	  handsomeTH1(hCrossSec_mc[i][j],3);
	  handsomeTH1(hCrossSec_m_data[i][j],3);
	  handsomeTH1(hCrossSec_m_mc[i][j],3);
	}
	else if (j == 2) {
	  handsomeTH1(hCrossSec_data[i][j],4);
	  handsomeTH1(hCrossSec_mc[i][j],4);
	  handsomeTH1(hCrossSec_m_data[i][j],4);
	  handsomeTH1(hCrossSec_m_mc[i][j],4);
	}
	else if (j == 3) {
	  handsomeTH1(hCrossSec_data[i][j],6);
	  handsomeTH1(hCrossSec_mc[i][j],6);
	  handsomeTH1(hCrossSec_m_data[i][j],6);
	  handsomeTH1(hCrossSec_m_mc[i][j],6);
	}
	else if (j == 4) {
	  handsomeTH1(hCrossSec_data[i][j],7);
	  handsomeTH1(hCrossSec_mc[i][j],7);
	  handsomeTH1(hCrossSec_m_data[i][j],7);
	  handsomeTH1(hCrossSec_m_mc[i][j],7);
	}
	else if (j == 5) {
	  handsomeTH1(hCrossSec_data[i][j],9);
	  handsomeTH1(hCrossSec_mc[i][j],9);
	  handsomeTH1(hCrossSec_m_data[i][j],9);
	  handsomeTH1(hCrossSec_m_mc[i][j],9);	  
	}
	else if (j == 6) {
	  handsomeTH1(hCrossSec_data[i][j],30);
	  handsomeTH1(hCrossSec_mc[i][j],30);
	  handsomeTH1(hCrossSec_m_data[i][j],30);
	  handsomeTH1(hCrossSec_m_mc[i][j],30);
	}
	c1->cd(1);
	TH1ScaleByWidth(hCrossSec_m_data[i][j]);
	TH1ScaleByWidth(hCrossSec_m_mc[i][j]);
	TH1ScaleByRapidity(hCrossSec_m_data[i][j],i);
	TH1ScaleByRapidity(hCrossSec_m_mc[i][j],i);
       	hCrossSec_mc[i][j]->SetMarkerStyle(kOpenCircle);
	fixedFontHist(hCrossSec_data[i][j],2.5,2,20);
	fixedFontHist(hCrossSec_mc[i][j],2.5,2,20);
	hCrossSec_data[i][j]->SetAxisRange(5e-9,5e17,"Y");
	hCrossSec_data[i][j]->SetYTitle("#frac{1}{#LT T_{AA} #GT} #frac{1}{N_{evt}} #frac{d#sigma}{dp_{T}} [nb/GeV]");
	if (j == 0) hCrossSec_data[i][j]->Draw();
	else hCrossSec_data[i][j]->Draw("same");
	hCrossSec_mc[i][j]->Draw("same");
	//	cout << "data eta " << i << hCrossSec_data[i][j]->GetEntries() << endl;
	c1->cd(2);
	char ratioName[256];
	sprintf(ratioName,"hCrossSec_ratio_eta%d_cent%d",i,j);
	hCrossSec_ratio[i][j] = (TH1D*)hCrossSec_data[i][j]->Clone(ratioName);
	hCrossSec_ratio[i][j]->Divide(hCrossSec_mc[i][j]);
	hCrossSec_ratio[i][j]->SetXTitle("#it{p}_{T}^{Reco}");
	hCrossSec_ratio[i][j]->SetYTitle("Data/MC");
	hCrossSec_ratio[i][j]->SetNdivisions(505,"X");
	fixedFontHist(hCrossSec_ratio[i][j],3,2,20);
	if (j == 0) hCrossSec_ratio[i][j]->Draw();
	else hCrossSec_ratio[i][j]->Draw("same");
	
	c2->cd(1);
	hCrossSec_m_mc[i][j]->SetMarkerStyle(kOpenCircle);
	fixedFontHist(hCrossSec_m_data[i][j],2.5,2,20);
	fixedFontHist(hCrossSec_m_mc[i][j],2.5,2,20);
	hCrossSec_m_data[i][j]->SetAxisRange(5e-9,5e17,"Y");
	hCrossSec_m_data[i][j]->SetYTitle("#frac{1}{#LT T_{AA} #GT} #frac{1}{N_{evt}} #frac{d#sigma}{dp_{T}} [nb/GeV]");
	if (j == 0) hCrossSec_m_data[i][j]->Draw();
	else hCrossSec_m_data[i][j]->Draw("same");
	hCrossSec_m_mc[i][j]->Draw("same");
	//	cout << "data eta " << i << hCrossSec_m_data[i][j]->GetEntries() << endl;
	c2->cd(2);

	sprintf(ratioName,"hCrossSec_m_ratio_eta%d_cent%d",i,j);
	hCrossSec_m_ratio[i][j] = (TH1D*)hCrossSec_m_data[i][j]->Clone(ratioName);
	hCrossSec_m_ratio[i][j]->Divide(hCrossSec_m_mc[i][j]);
	hCrossSec_m_ratio[i][j]->SetXTitle("m^{2} / #it{p}_{T}^{2}");
	hCrossSec_m_ratio[i][j]->SetYTitle("Data/MC");
	hCrossSec_m_ratio[i][j]->SetNdivisions(505,"X");
	fixedFontHist(hCrossSec_m_ratio[i][j],3,2,20);
	if (j == 0) hCrossSec_m_ratio[i][j]->Draw();
	else hCrossSec_m_ratio[i][j]->Draw("same");
      }

      //      gPad->SetLogy();
      c1->cd(1);
      gPad->SetLogy();
      drawCentrality(kSample,-1,0.70,0.86,1,24);
      TLegend *leg1 = new TLegend(0.3271854,0.6394175,0.9049639,0.8533743,NULL,"brNDC");      
      easyLeg(leg1,"");
      leg1->SetNColumns(2);
      leg1->AddEntry(hCrossSec_data[i][6], "Data: 60-80% (x10^{12})","plfe");
      leg1->AddEntry(hCrossSec_mc[i][6], "MC: 60-80% (x10^{12})","plfe");
      leg1->AddEntry(hCrossSec_data[i][5], "Data: 50-60% (x10^{10})","plfe");
      leg1->AddEntry(hCrossSec_mc[i][5], "MC: 50-60% (x10^{10})","plfe");
      leg1->AddEntry(hCrossSec_data[i][4], "Data: 40-50% (x10^{8})","plfe");
      leg1->AddEntry(hCrossSec_mc[i][4], "MC: 40-50% (x10^{8})","plfe");
      leg1->AddEntry(hCrossSec_data[i][3], "Data: 30-40% (x10^{6})","plfe");
      leg1->AddEntry(hCrossSec_mc[i][3], "MC: 30-40% (x10^{6})","plfe");
      leg1->AddEntry(hCrossSec_data[i][2], "Data: 20-30% (x10^{4})","plfe");
      leg1->AddEntry(hCrossSec_mc[i][2], "MC: 20-30% (x10^{4})","plfe");
      leg1->AddEntry(hCrossSec_data[i][1], "Data: 10-20% (x10^{2})","plfe");
      leg1->AddEntry(hCrossSec_mc[i][1], "MC: 10-20% (x10^{2})","plfe");      
      leg1->AddEntry(hCrossSec_data[i][0], "Data: 0-10%","plfe");
      leg1->AddEntry(hCrossSec_mc[i][0], "MC: 0-10%","plfe");
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
      leg2->AddEntry(hCrossSec_m_data[i][6], "Data: 60-80% (x10^{12})","plfe");
      leg2->AddEntry(hCrossSec_m_mc[i][6], "MC: 60-80% (x10^{12})","plfe");
      leg2->AddEntry(hCrossSec_m_data[i][5], "Data: 50-60% (x10^{10})","plfe");
      leg2->AddEntry(hCrossSec_m_mc[i][5], "MC: 50-60% (x10^{10})","plfe");
      leg2->AddEntry(hCrossSec_m_data[i][4], "Data: 40-50% (x10^{8})","plfe");
      leg2->AddEntry(hCrossSec_m_mc[i][4], "MC: 40-50% (x10^{8})","plfe");
      leg2->AddEntry(hCrossSec_m_data[i][3], "Data: 30-40% (x10^{6})","plfe");
      leg2->AddEntry(hCrossSec_m_mc[i][3], "MC: 30-40% (x10^{6})","plfe");
      leg2->AddEntry(hCrossSec_m_data[i][2], "Data: 20-30% (x10^{4})","plfe");
      leg2->AddEntry(hCrossSec_m_mc[i][2], "MC: 20-30% (x10^{4})","plfe");
      leg2->AddEntry(hCrossSec_m_data[i][1], "Data: 10-20% (x10^{2})","plfe");
      leg2->AddEntry(hCrossSec_m_mc[i][1], "MC: 10-20% (x10^{2})","plfe");      
      leg2->AddEntry(hCrossSec_m_data[i][0], "Data: 0-10%","plfe");
      leg2->AddEntry(hCrossSec_m_mc[i][0], "MC: 0-10%","plfe");


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

      TH1ScaleByWidth(hCrossSec_data[i][0]);
      TH1ScaleByWidth(hCrossSec_mc[i][0]);
      TH1ScaleByRapidity(hCrossSec_data[i][0], i);
      TH1ScaleByRapidity(hCrossSec_mc[i][0], i);
      TH1ScaleByWidth(hCrossSec_m_data[i][0]);
      TH1ScaleByWidth(hCrossSec_m_mc[i][0]);
      TH1ScaleByRapidity(hCrossSec_m_data[i][0], i);
      TH1ScaleByRapidity(hCrossSec_m_mc[i][0], i);
 
      // defined in commonUtility: 1 is black; 2 is red
      if (i == 5) {
	handsomeTH1(hCrossSec_data[i][0],kMagenta+1);
	handsomeTH1(hCrossSec_mc[i][0],kMagenta+1);
	handsomeTH1(hCrossSec_m_data[i][0],kMagenta+1);
	handsomeTH1(hCrossSec_m_mc[i][0],kMagenta+1);
	hCrossSec_m_mc[i][0]->SetMarkerStyle(kOpenCircle);
	hCrossSec_mc[i][0]->SetMarkerStyle(kOpenCircle);
	hCrossSec_m_data[i][0]->SetMarkerStyle(kFullCircle);
	hCrossSec_data[i][0]->SetMarkerStyle(kFullCircle);
      }
      
      else if (i == 0) {
	handsomeTH1(hCrossSec_data[i][0],kGreen+1);
	handsomeTH1(hCrossSec_mc[i][0],kGreen+1);
	handsomeTH1(hCrossSec_m_data[i][0],kGreen+1);
	handsomeTH1(hCrossSec_m_mc[i][0],kGreen+1);
	hCrossSec_m_mc[i][0]->SetMarkerStyle(kOpenCross);
	hCrossSec_mc[i][0]->SetMarkerStyle(kOpenCross);
	hCrossSec_m_data[i][0]->SetMarkerStyle(kFullCross);
	hCrossSec_data[i][0]->SetMarkerStyle(kFullCross);
      }
      else if (i == 1) {
	handsomeTH1(hCrossSec_data[i][0],kBlue-1);
	handsomeTH1(hCrossSec_mc[i][0],kBlue-1);
	handsomeTH1(hCrossSec_m_data[i][0],kBlue-1);
	handsomeTH1(hCrossSec_m_mc[i][0],kBlue-1);
	hCrossSec_m_mc[i][0]->SetMarkerStyle(kOpenDiamond);
	hCrossSec_mc[i][0]->SetMarkerStyle(kOpenDiamond);
	hCrossSec_m_data[i][0]->SetMarkerStyle(kFullDiamond);
	hCrossSec_data[i][0]->SetMarkerStyle(kFullDiamond);
      }
      else if (i == 2) {
	handsomeTH1(hCrossSec_data[i][0],kRed+1);
	handsomeTH1(hCrossSec_mc[i][0],kRed+1);
	handsomeTH1(hCrossSec_m_data[i][0],kRed+1);
	handsomeTH1(hCrossSec_m_mc[i][0],kRed+1);
	hCrossSec_m_mc[i][0]->SetMarkerStyle(kOpenSquare);
	hCrossSec_mc[i][0]->SetMarkerStyle(kOpenSquare);
	hCrossSec_m_data[i][0]->SetMarkerStyle(kFullSquare);
	hCrossSec_data[i][0]->SetMarkerStyle(kFullSquare);
      }
      else if (i == 3) {
	handsomeTH1(hCrossSec_data[i][0],kOrange+5);
	handsomeTH1(hCrossSec_mc[i][0],kOrange+5);
	handsomeTH1(hCrossSec_m_data[i][0],kOrange+5);
	handsomeTH1(hCrossSec_m_mc[i][0],kOrange+5);
	hCrossSec_m_mc[i][0]->SetMarkerStyle(kOpenStar);
	hCrossSec_mc[i][0]->SetMarkerStyle(kOpenStar);
	hCrossSec_m_data[i][0]->SetMarkerStyle(kFullStar);
	hCrossSec_data[i][0]->SetMarkerStyle(kFullStar);
      }
      else if (i == 4) {
	handsomeTH1(hCrossSec_data[i][0],kPink+1);
	handsomeTH1(hCrossSec_mc[i][0],kPink+1);
	handsomeTH1(hCrossSec_m_data[i][0],kPink+1);
	handsomeTH1(hCrossSec_m_mc[i][0],kPink+1);
	hCrossSec_m_mc[i][0]->SetMarkerStyle(kOpenCrossX);
	hCrossSec_mc[i][0]->SetMarkerStyle(kOpenCrossX);
	hCrossSec_m_data[i][0]->SetMarkerStyle(kFullCrossX);
	hCrossSec_data[i][0]->SetMarkerStyle(kFullCrossX);
      }

      c1->cd(1);
      fixedFontHist(hCrossSec_data[i][0],2.5,2,20);
      fixedFontHist(hCrossSec_mc[i][0],2.5,2,20);

      hCrossSec_data[i][0]->SetAxisRange(9e-9,5e20,"Y");
      hCrossSec_data[i][0]->SetAxisRange(125,1000,"X");
      hCrossSec_data[i][0]->SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T}dy} [nb/GeV]");
      //      hCrossSec_data[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
      hCrossSec_mc[i][0]->SetAxisRange(9e-9,5e20,"Y");
      hCrossSec_mc[i][0]->SetAxisRange(125,1000,"X");
      hCrossSec_mc[i][0]->SetYTitle("#frac{d^{2}#sigma}{d#it{p}_{T}dy} [nb/GeV]");
      //      hCrossSec_mc[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
      if (i == 0) hCrossSec_data[i][0]->Draw();
      else hCrossSec_data[i][0]->Draw("same");
      hCrossSec_mc[i][0]->Draw("same");
      c1->cd(2);
      char ratioName[256];
      sprintf(ratioName,"hCrossSec_ratio_eta%d_cent%d",i,0);
      hCrossSec_ratio[i][0] = (TH1D*)hCrossSec_data[i][0]->Clone(ratioName);
      hCrossSec_ratio[i][0]->Divide(hCrossSec_mc[i][0]);
      hCrossSec_ratio[i][0]->SetXTitle("#it{p}_{T}^{Reco}");
      hCrossSec_ratio[i][0]->SetYTitle("Data/MC");
      //      hCrossSec_ratio[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
      hCrossSec_ratio[i][0]->SetNdivisions(505,"X");
      hCrossSec_ratio[i][0]->SetAxisRange(.75,1.25,"Y");
      fixedFontHist(hCrossSec_ratio[i][0],3,2,20);
      if (i == 0) hCrossSec_ratio[i][0]->Draw();
      else hCrossSec_ratio[i][0]->Draw("same");

      
      c2->cd(1);
      hCrossSec_m_mc[i][0]->SetMarkerStyle(kOpenCircle);
      fixedFontHist(hCrossSec_m_data[i][0],2.5,2,20);
      fixedFontHist(hCrossSec_m_mc[i][0],2.5,2,20);
      //      hCrossSec_m_data[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
       hCrossSec_m_data[i][0]->SetAxisRange(9e-9,5e26,"Y");
      hCrossSec_m_data[i][0]->SetYTitle("#frac{d^{2}#sigma}{d(m/#it{p}_{T})^{2}dy} [nb]");
      //      hCrossSec_m_mc[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
      hCrossSec_m_mc[i][0]->SetAxisRange(9e-9,5e26,"Y");
      hCrossSec_m_mc[i][0]->SetYTitle("#frac{d^{2}#sigma}{d(m/#it{p}_{T})^{2}}_{T}dy} [nb]"); 
      if (i == 0) hCrossSec_m_data[i][0]->Draw();
      else hCrossSec_m_data[i][0]->Draw("same");
      hCrossSec_m_mc[i][0]->Draw("same");
      c2->cd(2);
      sprintf(ratioName,"hCrossSec_m_ratio_eta%d_cent%d",i,0);
      hCrossSec_m_ratio[i][0] = (TH1D*)hCrossSec_m_data[i][0]->Clone(ratioName);
      hCrossSec_m_ratio[i][0]->Divide(hCrossSec_m_mc[i][0]);
      hCrossSec_m_ratio[i][0]->SetXTitle("m^{2} / #it{p}_{T}^{2}");
      hCrossSec_m_ratio[i][0]->SetYTitle("Data/MC");
      hCrossSec_m_ratio[i][0]->SetNdivisions(505,"X");
      //      hCrossSec_m_ratio[i][0]->GetXaxis()->SetMoreLogLabels(kTRUE);
      hCrossSec_m_ratio[i][0]->SetAxisRange(.6,1.4,"Y");
      fixedFontHist(hCrossSec_m_ratio[i][0],3,2,20);
      if (i == 0) hCrossSec_m_ratio[i][0]->Draw();
      else hCrossSec_m_ratio[i][0]->Draw("same");

    }
    //    c1->cd(2);
    //    gPad->SetLogy();
    c1->cd(1);

    drawCentrality(kSample,0,0.77,0.94,1,24); 
    TLegend *leg1 = new TLegend(0.3571854,0.6194175,0.9849639,0.9333743,NULL,"brNDC");
    easyLeg(leg1,"");
    leg1->SetNColumns(2);
    leg1->AddEntry(hCrossSec_data[5][0], "Data: |y| < 2.1 (x10^{10})","plfe");
    leg1->AddEntry(hCrossSec_mc[5][0], "MC: |y| < 2.1 (x10^{10})","plfe");
    leg1->AddEntry(hCrossSec_data[0][0], "Data: |y| < 0.3 (x10^{8})","plfe");
    leg1->AddEntry(hCrossSec_mc[0][0], "MC: |y| < 0.3 (x10^{8})","plfe");    
    leg1->AddEntry(hCrossSec_data[1][0], "Data: 0.3 #leq |y| < 0.8 (x10^{6})","plfe");
    leg1->AddEntry(hCrossSec_mc[1][0], "MC: 0.3 #leq |y| < 0.8 (x10^{6})","plfe");    
    leg1->AddEntry(hCrossSec_data[2][0], "Data: 0.8 #leq |y| < 1.2 (x10^{4})","plfe");
    leg1->AddEntry(hCrossSec_mc[2][0], "MC: 0.8 #leq |y| < 1.2 (x10^{4})","plfe");    
    leg1->AddEntry(hCrossSec_data[3][0], "Data: 1.2 #leq |y| < 1.6 (x10^{2})","plfe");
    leg1->AddEntry(hCrossSec_mc[3][0], "MC: 1.2 #leq |y| < 1.6 (x10^{2})","plfe");    
    leg1->AddEntry(hCrossSec_data[4][0], "Data: 1.6 #leq |y| < 2.1","plfe");
    leg1->AddEntry(hCrossSec_mc[4][0], "MC: 1.6 #leq |y| < 2.1","plfe");
     
    leg1->Draw();
    ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    c1->cd(2);


    c1->SaveAs(Form("Plots/CrossSection/CrossSection_kSample%d_eta_NoReweighting.png",kSample));
    c1->SaveAs(Form("Plots/CrossSection/CrossSection_kSample%d_eta_NoReweighting.root",kSample));
    delete c1;


    c2->cd(1);
    //    c2->SetLogx();
    //    c2->SetLogy();
    drawCentrality(kSample,0,0.77,0.94,1,24); 
    TLegend *leg2 = new TLegend(0.3571854,0.6194175,0.9849639,0.9333743,NULL,"brNDC");
    easyLeg(leg2,"");
    leg2->SetNColumns(2);
    leg2->AddEntry(hCrossSec_m_data[5][0], "Data: |y| < 2.1 (x10^{10})","plfe");
    leg2->AddEntry(hCrossSec_m_mc[5][0], "MC: |y| < 2.1 (x10^{10})","plfe");     
    leg2->AddEntry(hCrossSec_m_data[0][0], "Data: |y| < 0.3 (x10^{8})","plfe");
    leg2->AddEntry(hCrossSec_m_mc[0][0], "MC: |y| < 0.3 (x10^{8})","plfe");    
    leg2->AddEntry(hCrossSec_m_data[1][0], "Data: 0.3 #leq |y| < 0.8 (x10^{6})","plfe");
    leg2->AddEntry(hCrossSec_m_mc[1][0], "MC: 0.3 #leq |y| < 0.8 (x10^{6})","plfe");    
    leg2->AddEntry(hCrossSec_m_data[2][0], "Data: 0.8 #leq |y| < 1.2 (x10^{4})","plfe");
    leg2->AddEntry(hCrossSec_m_mc[2][0], "MC: 0.8 #leq |y| < 1.2 (x10^{4})","plfe");    
    leg2->AddEntry(hCrossSec_m_data[3][0], "Data: 1.2 #leq |y| < 1.6 (x10^{2})","plfe");
    leg2->AddEntry(hCrossSec_m_mc[3][0], "MC: 1.2 #leq |y| < 1.6 (x10^{2})","plfe");    
    leg2->AddEntry(hCrossSec_m_data[4][0], "Data: 1.6 #leq |y| < 2.1","plfe");
    leg2->AddEntry(hCrossSec_m_mc[4][0], "MC: 1.6 #leq |y| < 2.1","plfe");
    leg2->Draw();
    ATLASLabel(0.25, 0.92, "Internal",0.085,0.28);//, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
    c2->SaveAs(Form("Plots/CrossSection/CrossSection_m2pt2_kSample%d_eta_NoReweighting.png",kSample));
    c2->SaveAs(Form("Plots/CrossSection/CrossSection_m2pt2_kSample%d_eta_NoReweighting.root",kSample));
    delete c2;
  }

  return 0;
}
