#define PlotJetMass_ppMC_cxx
#include "PlotJetMass_ppMC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void PlotJetMass_ppMC::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PlotJetMass_ppMC.C
//      root> PlotJetMass_ppMC t
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

   TH2F *h2_jetMass2Pt_mc = new TH2F("jetMass2Pt_mc","jetMass2Pt_mc",125,100,700,125,-100,300);
   TH2F *h2_jetMass2OverPt2VsPt2_mc = new TH2F("jetMass2OverPt2VsPt2_mc","jetMass2Pt2_mc",125,10000,500000,125,-.005,.1);
   TH2F *h2_jetMass2OverPt2VsPt_mc = new TH2F("jetMass2OverPt2VsPt_mc","jetMass2OverPt2VsPt_mc",125,100,700,125,-.005,.01);

   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      h2_jetMass2OverPt2VsPt2->Fill(jets_ptCalib*jets_ptCalib,jets_mass2Calib/(jets_ptCalib*jets_ptCalib) );
      h2_jetMass2OverPt2VsPt->Fill(jets_ptCalib, jets_mass2Calib/(jets_ptCalib*jets_ptCalib));
      h2_jetMass2Pt->Fill(jets_ptCalib, jets_mass2Calib);
      //      cout << "mass2 " << jets_mass2Calib << " ptCalib " << jets_ptCalib << endl;
   }
   cout << " h1 entries " << h2_jetMass2OverPt2VsPt2->GetEntries() << endl;
   TString outFileName;
   gStyle->SetPalette(kBird);
   TString fileName; 
   TCanvas *c1 = new TCanvas("c1","c1");
   h2_jetMass2OverPt2VsPt2->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");      
   h2_jetMass2OverPt2VsPt2->GetXaxis()->SetTitle("#it{p}_{T}^{2}");
   h2_jetMass2OverPt2VsPt2->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2OverPt2VsPt2->SetContour(99);
   h2_jetMass2OverPt2VsPt2->Draw("colz");
   c1->SetLogz();
   outFileName.Form("plots/ppData/mass2OverPt2VsPt2_ppData.png");
   c1->SaveAs(outFileName.Data());
   outFileName.Form("plots/ppData/mass2OverPt2VsPt2_ppData.root");
   c1->SaveAs(outFileName.Data());

   TCanvas *c1b = new TCanvas("c1b","c1b");
   TProfile *p1 = h2_jetMass2OverPt2VsPt2->ProfileX("px1");
   p1->Draw();
   p1->SetMarkerSize(0.25);
   outFileName.Form("plots/ppData/mass2OverPt2VsPt2_ppData_prof.png");
   c1b->SaveAs(outFileName.Data());
   outFileName.Form("plots/ppData/mass2OverPt2VsPt2_ppData_prof.root");
   c1b->SaveAs(outFileName.Data());
   
   TCanvas *c2 = new TCanvas("c2","c2");
   h2_jetMass2Pt->GetYaxis()->SetTitle("m^{2}");      
   h2_jetMass2Pt->GetXaxis()->SetTitle("#it{p}_{T}");
   h2_jetMass2Pt->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2Pt->SetContour(99);
   h2_jetMass2Pt->Draw("colz");
   c2->SetLogz();
   outFileName.Form("plots/ppData/mass2Pt_ppData.png");
   c2->SaveAs(outFileName.Data());
   outFileName.Form("plots/ppData/mass2Pt_ppData.root");
   c2->SaveAs(outFileName.Data());

   TCanvas *c2b = new TCanvas("c2b","c2b");
   TProfile *p2 = h2_jetMass2Pt->ProfileX("px2");
   p2->Draw();
   p2->SetMarkerSize(0.25);
   outFileName.Form("plots/ppData/mass2Pt_ppData_prof.png");
   c2b->SaveAs(outFileName.Data());
   outFileName.Form("plots/ppData/mass2Pt_ppData_prof.root");
   c2b->SaveAs(outFileName.Data());

   TCanvas *c3 = new TCanvas("c3","c3");
   h2_jetMass2OverPt2VsPt->GetYaxis()->SetTitle("m^{2}/#it{p}_{T}^{2}");      
   h2_jetMass2OverPt2VsPt->GetXaxis()->SetTitle("#it{p}_{T}");
   h2_jetMass2OverPt2VsPt->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2OverPt2VsPt->SetContour(99);
   h2_jetMass2OverPt2VsPt->Draw("colz");
   c3->SetLogz();
   outFileName.Form("plots/ppData/mass2OverPt2VsPt_ppData.png");
   c3->SaveAs(outFileName.Data());
   outFileName.Form("plots/ppData/mass2OverPt2VsPt_ppData.root");
   c3->SaveAs(outFileName.Data());

   TCanvas *c3b = new TCanvas("c3b","c3b");
   TProfile *p3 = h2_jetMass2OverPt2VsPt->ProfileX("px3");
   p3->Draw();
   p3->SetMarkerSize(0.25);
   outFileName.Form("plots/ppData/mass2OverPt2VsPt_ppData_prof.png");
   c3b->SaveAs(outFileName.Data());
   outFileName.Form("plots/ppData/mass2OverPt2VsPt_ppData_prof.root");
   c3b->SaveAs(outFileName.Data());
   
   delete c1; delete c1b;
   delete c2; delete c2b;
   delete c3; delete c3b;
}
