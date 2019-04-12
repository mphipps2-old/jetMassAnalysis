#define getEtaPhiMap_cxx
#include "getEtaPhiMap.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void getEtaPhiMap::Loop()
{
//   In a ROOT session, you can do:
//      root> .L getEtaPhiMap.C
//      root> getEtaPhiMap t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   gStyle->SetPalette(kBird);
   TH2F* etaPhiAvgPt;
   TH2F* etaPhiStats;
   char name[256];
   sprintf(name,"etaPhiAvgPt_pbpb");
   etaPhiAvgPt = new TH2F(name,name, 56,-2.1,2.1, 62,-3.14,3.14);
   sprintf(name,"etaPhiStats_pbpb");
   etaPhiStats = new TH2F(name,name, 56,-2.1,2.1, 62,-3.14,3.14);
   
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;
     if (jets_pt < 100 || jets_eta < -2.1 || jets_eta > 2.1) continue; 
     etaPhiAvgPt->Fill(jets_eta,jets_phi,jets_pt);
     etaPhiStats->Fill(jets_eta,jets_phi);      
   }
   float norm = etaPhiAvgPt->GetEntries();
   etaPhiAvgPt->Divide(etaPhiStats);
   TCanvas *c1 = new TCanvas("c1","c1");
   etaPhiAvgPt->Draw("colz");
   etaPhiAvgPt->GetXaxis()->SetTitle("#eta");
   etaPhiAvgPt->GetYaxis()->SetTitle("#phi");
   etaPhiAvgPt->SetContour(99);
   sprintf(name,"calibration/etaPhiAvgPt_pbpb.png");
   c1->SaveAs(name);
   sprintf(name,"calibration/etaPhiAvgPt_pbpb.root");
   c1->SaveAs(name);
   TCanvas *c2 = new TCanvas("c2","c2");
   etaPhiStats->Draw("colz");
   etaPhiStats->GetXaxis()->SetTitle("#eta");
   etaPhiStats->GetYaxis()->SetTitle("#phi");
   etaPhiStats->SetContour(99);
   sprintf(name,"calibration/etaPhiStats_pbpb.png");
   c2->SaveAs(name);
   sprintf(name,"calibration/etaPhiStats_pbpb.root");
   c2->SaveAs(name);
   delete c1; delete c2;
}
