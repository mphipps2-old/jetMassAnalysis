#define MassComparison_cxx
#include "MassComparison.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MassComparison::Loop()
{
//   In a ROOT session, you can do:
//      root> .L MassComparison.C
//      root> MassComparison t
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
   TH1F * unsubM2 = new TH1F("unsubM2","unsubM2",100,-50,10000);
   TH1F * subM2 = new TH1F("subM2","subM2",100,-50,10000);
   TH1F * towerM2 = new TH1F("towerM2","towerM2",100,-50,10000);
   TH1F * calibM2 = new TH1F("calibM2","calibM2",100,-50,10000);
   TH1F * calibM = new TH1F("calibM","calibM",100,-50,100);
   TH1F * calibMDirect = new TH1F("calibMDirect","calibMDirect",100,-50,100);
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      unsubM2->Fill(jets_mass2Unsub);
      subM2->Fill(jets_mass2Sub);
      towerM2->Fill(jets_mass2Constit);
      calibM2->Fill(jets_mass2Calib);
      calibM->Fill(jets_massCalib);
      calibMDirect->Fill(jets_massCalibDirect);
      cout << " Unsubtracted m2 " << jets_mass2Unsub << " Subtracted m2 " << jets_mass2Sub << " Tower m2 " << jets_mass2Constit << " mass2Calib calc " << jets_mass2Calib << " mass2Calib direct " << jets_massCalibDirect * jets_massCalibDirect << " truthM2 " << jets_genMass2 << endl << endl;
   }
   TCanvas *c1 = new TCanvas("c1","c1");
   unsubM2->Draw();
   c1->SaveAs("Plots/massDefinition/UnsubM2_ppData.png");
   delete c1;
   TCanvas *c2 = new TCanvas("c2","c2");
   subM2->Draw();
   c2->SaveAs("Plots/massDefinition/SubM2_ppData.png");
   delete c2;
   TCanvas *c3 = new TCanvas("c3","c3");
   towerM2->Draw();
   c3->SaveAs("Plots/massDefinition/towerM2_ppData.png");
   delete c3;
   TCanvas *c4 = new TCanvas("c4","c4");
   calibM2->Draw();
   c4->SaveAs("Plots/massDefinition/calibM2_ppData.png");
   delete c4;
   TCanvas *c5 = new TCanvas("c5","c5");
   calibM->Draw();
   c5->SaveAs("Plots/massDefinition/calibM_ppData.png");
   delete c5;
   TCanvas *c6 = new TCanvas("c6","c6");
   calibMDirect->Draw();
   c6->SaveAs("Plots/massDefinition/calibMDirect_ppData.png");
   delete c6;
}
