#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "../getSdHists.C"
#include "../ntupleDefinition_v50.h"
#include "../commonUtility.h"
#include "../jzWeight.h"
#include "unfoldingUtil.h"

#include "../JssUtils.h"
#include <TPaletteAxis.h>

double statFrac = 001;
double fracStstData = 001;
bool doUnfData = true ;

void getFcalData(TH1D* hh=0);
void getFcalMc(TH1D* hh=0);


void getFcalWeights() {
  TH1::SetDefaultSumw2();
  
  TH1D* hdata = new TH1D("hdata",";FCal E_{T} ;",100,-0.5,6);
  TH1D* hmc = (TH1D*)hdata->Clone("hmc");
  
  getFcalMc(hdata);
  getFcalData(hmc);
  TH1D* hWeight = (TH1D*)hdata->Clone("fcalWeight");
  hWeight->Divide(hmc);

  TCanvas* c1 = new TCanvas("c1","",600,600);
  makeEfficiencyCanvas(c1,1, 0.05, 0.01, 0.1, 0.3, 0.01);
  c1->cd(1);
  handsomeTH1(hdata,2);
  handsomeTH1(hmc,1);
  hdata->Draw();
  hmc->Draw("same");
  c1->cd(2);
  hWeight->Draw();
  cleverRange(hWeight);

  TFile* fout = new TFile("reweightFactors/fcalWeight_PbPb_v1.root","recreate");
  hWeight->Write();
  fout->Close();
}



void getFcalMc(TH1D* hh) {
  TH1::SetDefaultSumw2();
  hh->Reset();
  
  TString jz2;
  TString jz3;
  TString jz4;
  
  jz2 = "jetSubstructure_MC_HION9_pbpb_v50_jz2.root";
  jz3 = "jetSubstructure_MC_HION9_pbpb_v50_jz3.root";
  jz4 = "jetSubstructure_MC_HION9_pbpb_v50_jz4.root";
  
  jetSubStr myJetMc;
  TBranch       *b_myJetSubMc;

  TFile* fjz2 = new TFile(Form("../ntuples/%s",jz2.Data()));
  TTree* tr2 = (TTree*)fjz2->Get("tr");
  tr2->SetBranchAddress("jets", &(myJetMc.cent), &b_myJetSubMc);

  TFile* fjz3 = new TFile(Form("../ntuples/%s",jz3.Data()));
  TTree* tr3 = (TTree*)fjz3->Get("tr");
  tr3->SetBranchAddress("jets", &(myJetMc.cent), &b_myJetSubMc);

  TFile* fjz4 = new TFile(Form("../ntuples/%s",jz4.Data()));
  TTree* tr4 = (TTree*)fjz4->Get("tr");
  tr4->SetBranchAddress("jets", &(myJetMc.cent), &b_myJetSubMc);

  for ( int ijz =2 ; ijz<=4 ; ijz++) {
    TTree* tr;
    double jzNorm=0;
    if ( ijz==2)  {
      tr = tr2;
      jzNorm = hi9EvtWgtJZ2;
    }
    else if ( ijz==3)  {
      tr = tr3;
      jzNorm = hi9EvtWgtJZ3;
    }
    else if ( ijz==4)  {
      tr = tr4;
      jzNorm = hi9EvtWgtJZ4;
    }

    cout << "Scanning JZ"<<ijz<<" file.  Total events = " << tr->GetEntries() << endl;
    for (Int_t i= 0; i<tr->GetEntries() ; i++) {
      if ( i > tr->GetEntries() * statFrac ) continue;
      if (i%2==0)  continue;

      tr->GetEntry(i);

      hh->Fill( myJetMc.fcalet,  myJetMc.weight * jzNorm);
    }
  }
}

void getFcalData(TH1D* hh) { 


  TH1::SetDefaultSumw2();
  hh->Reset();
  TString fname;
  fname = "jetSubstructure_data_HION9_v50_r4_pp.root";
  TFile* fData = new TFile(Form("../ntuples/%s",fname.Data()));
  TTree* tr = (TTree*)fData->Get("tr");
  jetSubStr myJet;
  TBranch       *b_myJet;
  tr->SetBranchAddress("jets", &(myJet.cent), &b_myJet);
  cout << "data entries = " << tr->GetEntries() << endl;
  
  for (Int_t i= 0; i<tr->GetEntries() ; i++) {
    tr->GetEntry(i);
 
    hh->Fill( myJet.fcalet,  myJet.weight);
  }
}

