#define getEtaPhiMap_MCpbpb_cxx
#include "getEtaPhiMap_MCpbpb.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void getEtaPhiMap_MCpbpb::Loop()
{
//   In a ROOT session, you can do:
//      root> .L getEtaPhiMap_MCpbpb.C
//      root> getEtaPhiMap_MCpbpb t
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
   char name[256];
   sprintf(name,"etaPhiAvgPt_cent0_MCpbpb");
   TH2F* etaPhiAvgPt_cent0 = new TH2F(name,name, 56,-2.1,2.1, 62,-3.14,3.14);
   sprintf(name,"etaPhiStats_cent0_MCpbpb");
   TH2F* etaPhiStats_cent0 = new TH2F(name,name, 56,-2.1,2.1, 62,-3.14,3.14);
   sprintf(name,"etaPhiStats_cent1_MCpbpb");
   TH2F* etaPhiStats_cent1 = new TH2F(name,name, 56,-2.1,2.1, 62,-3.14,3.14);
   sprintf(name,"etaPhiStats_cent2_MCpbpb");
   TH2F* etaPhiStats_cent2 = new TH2F(name,name, 56,-2.1,2.1, 62,-3.14,3.14);
   sprintf(name,"etaPhiStats_cent3_MCpbpb");
   TH2F* etaPhiStats_cent3 = new TH2F(name,name, 56,-2.1,2.1, 62,-3.14,3.14);
   sprintf(name,"etaPhiStats_cent4_MCpbpb");
   TH2F* etaPhiStats_cent4 = new TH2F(name,name, 56,-2.1,2.1, 62,-3.14,3.14);
   sprintf(name,"etaPhiStats_cent5_MCpbpb");
   TH2F* etaPhiStats_cent5 = new TH2F(name,name, 56,-2.1,2.1, 62,-3.14,3.14);
   sprintf(name,"etaPhiAvgPt_cent6_MCpbpb");
   TH2F* etaPhiAvgPt_cent6 = new TH2F(name,name, 56,-2.1,2.1, 62,-3.14,3.14);
   sprintf(name,"etaPhiStats_cent6_MCpbpb");
   TH2F* etaPhiStats_cent6 = new TH2F(name,name, 56,-2.1,2.1, 62,-3.14,3.14);
   
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;
     if (jets_pt < 100 || jets_eta < -2.1 || jets_eta > 2.1) continue;
     if (jets_cent == 0) {
       etaPhiAvgPt_cent0->Fill(jets_eta,jets_phi,jets_pt);
       etaPhiStats_cent0->Fill(jets_eta,jets_phi);
     }
     else if (jets_cent == 1) {
       etaPhiStats_cent1->Fill(jets_eta,jets_phi);
     }
     else if (jets_cent == 2) {
       etaPhiStats_cent2->Fill(jets_eta,jets_phi);
     }
     else if (jets_cent == 3) {
       etaPhiStats_cent3->Fill(jets_eta,jets_phi);
     }
     else if (jets_cent == 4) {
       etaPhiStats_cent4->Fill(jets_eta,jets_phi);
     }
     else if (jets_cent == 5) {
       etaPhiStats_cent5->Fill(jets_eta,jets_phi);
     }
     else if (jets_cent == 6) {
       etaPhiAvgPt_cent6->Fill(jets_eta,jets_phi,jets_pt);
       etaPhiStats_cent6->Fill(jets_eta,jets_phi);
     }
   }
   float norm = etaPhiAvgPt_cent0->GetEntries();
   etaPhiAvgPt_cent0->Divide(etaPhiStats_cent0);
   TCanvas *c1 = new TCanvas("c1","c1");
   etaPhiAvgPt_cent0->Draw("colz");
   etaPhiAvgPt_cent0->GetXaxis()->SetTitle("#eta");
   etaPhiAvgPt_cent0->GetYaxis()->SetTitle("#phi");
   etaPhiAvgPt_cent0->SetContour(99);
   sprintf(name,"calibration/etaPhiAvgPtcent0_MCpbpb.png");
   c1->SaveAs(name);
   sprintf(name,"calibration/etaPhiAvgPtcent0_MCpbpb.root");
   c1->SaveAs(name);
   
   TCanvas *c2 = new TCanvas("c2","c2");
   etaPhiStats_cent0->Draw("colz");
   etaPhiStats_cent0->GetXaxis()->SetTitle("#eta");
   etaPhiStats_cent0->GetYaxis()->SetTitle("#phi");
   etaPhiStats_cent0->SetContour(99);
   sprintf(name,"calibration/etaPhiStatscent0_MCpbpb.png");
   c2->SaveAs(name);
   sprintf(name,"calibration/etaPhiStatscent0_MCpbpb.root");
   c2->SaveAs(name);

   TCanvas *c3 = new TCanvas("c3","c3");
   etaPhiStats_cent1->Draw("colz");
   etaPhiStats_cent1->GetXaxis()->SetTitle("#eta");
   etaPhiStats_cent1->GetYaxis()->SetTitle("#phi");
   etaPhiStats_cent1->SetContour(99);
   sprintf(name,"calibration/etaPhiStatscent1_MCpbpb.png");
   c3->SaveAs(name);
   sprintf(name,"calibration/etaPhiStatscent1_MCpbpb.root");
   c3->SaveAs(name);

   TCanvas *c4 = new TCanvas("c4","c4");
   etaPhiStats_cent2->Draw("colz");
   etaPhiStats_cent2->GetXaxis()->SetTitle("#eta");
   etaPhiStats_cent2->GetYaxis()->SetTitle("#phi");
   etaPhiStats_cent2->SetContour(99);
   sprintf(name,"calibration/etaPhiStatscent2_MCpbpb.png");
   c4->SaveAs(name);
   sprintf(name,"calibration/etaPhiStatscent2_MCpbpb.root");
   c4->SaveAs(name);

   TCanvas *c5 = new TCanvas("c5","c5");
   etaPhiStats_cent3->Draw("colz");
   etaPhiStats_cent3->GetXaxis()->SetTitle("#eta");
   etaPhiStats_cent3->GetYaxis()->SetTitle("#phi");
   etaPhiStats_cent3->SetContour(99);
   sprintf(name,"calibration/etaPhiStatscent3_MCpbpb.png");
   c5->SaveAs(name);
   sprintf(name,"calibration/etaPhiStatscent3_MCpbpb.root");
   c5->SaveAs(name);

   TCanvas *c6 = new TCanvas("c6","c6");
   etaPhiStats_cent4->Draw("colz");
   etaPhiStats_cent4->GetXaxis()->SetTitle("#eta");
   etaPhiStats_cent4->GetYaxis()->SetTitle("#phi");
   etaPhiStats_cent4->SetContour(99);
   sprintf(name,"calibration/etaPhiStatscent4_MCpbpb.png");
   c2->SaveAs(name);
   sprintf(name,"calibration/etaPhiStatscent4_MCpbpb.root");
   c2->SaveAs(name);
   
   TCanvas *c7 = new TCanvas("c7","c7");
   etaPhiStats_cent5->Draw("colz");
   etaPhiStats_cent5->GetXaxis()->SetTitle("#eta");
   etaPhiStats_cent5->GetYaxis()->SetTitle("#phi");
   etaPhiStats_cent5->SetContour(99);
   sprintf(name,"calibration/etaPhiStatscent5_MCpbpb.png");
   c7->SaveAs(name);
   sprintf(name,"calibration/etaPhiStatscent5_MCpbpb.root");
   c7->SaveAs(name);

   norm = etaPhiAvgPt_cent6->GetEntries();
   etaPhiAvgPt_cent6->Divide(etaPhiStats_cent6);
   TCanvas *c8 = new TCanvas("c8","c8");
   etaPhiAvgPt_cent6->Draw("colz");
   etaPhiAvgPt_cent6->GetXaxis()->SetTitle("#eta");
   etaPhiAvgPt_cent6->GetYaxis()->SetTitle("#phi");
   etaPhiAvgPt_cent6->SetContour(99);
   sprintf(name,"calibration/etaPhiAvgPt_cent6_MCpbpb.png");
   c8->SaveAs(name);
   sprintf(name,"calibration/etaPhiAvgPt_cent6_MCpbpb.root");
   c8->SaveAs(name);
   
   TCanvas *c9 = new TCanvas("c9","c9");
   etaPhiStats_cent6->Draw("colz");
   etaPhiStats_cent6->GetXaxis()->SetTitle("#eta");
   etaPhiStats_cent6->GetYaxis()->SetTitle("#phi");
   etaPhiStats_cent6->SetContour(99);
   sprintf(name,"calibration/etaPhiStats_MCpbpb.png");
   c9->SaveAs(name);
   sprintf(name,"calibration/etaPhiStats_MCpbpb.root");
   c9->SaveAs(name);
   delete c1; delete c2; delete c3; delete c4; delete c5; delete c6; delete c7; delete c8; delete c9;
}
