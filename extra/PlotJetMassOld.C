#define PlotJetMass_cxx
#include "PlotJetMass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void PlotJetMass::Loop()
{
  gStyle->SetPalette(kBird);
//   In a ROOT session, you can do:
//      root> .L PlotJetMass.C
//      root> PlotJetMass t
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
   TH1F *h_jetMass = new TH1F("jetMass","jetMass",1000,-50,250);
   TH1F *h_jetMass2 = new TH1F("jetMass2","jetMass2",5000,-2000,30000);
   TH1F *h_jetGenMass = new TH1F("jetGenMass","jetGenMass",1000,-5,50);
   TH1F *h_jetGenMass2 = new TH1F("jetGenMass2","jetGenMass2",5000,-50,20000);
   TH2F *h2_jetMassPt = new TH2F("jetMassPt","jetMassPt",200,0,600,1000,-50,200);
   TH2F *h2_jetMass2Pt = new TH2F("jetMass2Pt","jetMass2Pt",200,0,800,1000,-300,5000);
   //   TH1F *h_jetPt = new TH1F("jetPt","jetPt",200,0,800);
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      h_jetMass->Fill(jets_mass);
      h_jetMass2->Fill(jets_mass2);
      h_jetMass->Fill(jets_mass);
      h_jetMass2->Fill(jets_mass2);
      h_jetGenMass->Fill(jets_genMass);
      h_jetGenMass2->Fill(jets_genMass2);
      h2_jetMassPt->Fill(jets_rawPt,jets_mass);
      h2_jetMass2Pt->Fill(jets_rawPt,jets_mass2);
      h_jetPt->Fill(jets_rawPt);

   }
   TString outFileName;
      
   /*
   TCanvas *c1 = new TCanvas();
   h_jetMass->Draw();
   h_jetMass->GetXaxis()->SetTitle("m calib");
   c1->SetLogy();

   outFileName.Form("plots/%s_mass.png",fileName.Data());
   c1->SaveAs(outFileName.Data());
   outFileName.Form("plots/%s_mass.root",fileName.Data());
   c1->SaveAs(outFileName.Data());
   delete c1;
   c1 = new TCanvas();
   c1->SetLogy();
   h_jetMass2->Draw();
   h_jetMass2->GetXaxis()->SetTitle("m^{2} calib");
   outFileName.Form("plots/%s_mass2.png",fileName.Data());
   c1->SaveAs(outFileName.Data());
   outFileName.Form("plots/%s_mass2.root",fileName.Data());
   c1->SaveAs(outFileName.Data());
   delete c1;
   */
   /*
   TCanvas *c2 = new TCanvas();
   c2->SetLogy();
   
   h2_jetMassPt->Draw("colz");
   h2_jetMassPt->GetXaxis()->SetTitle("p_{T} raw");      
   h2_jetMassPt->GetYaxis()->SetTitle("m raw");
   h2_jetMassPt->SetContour(99);
   outFileName.Form("plots/%s_massPt.png",fileName.Data());
   c2->SaveAs(outFileName.Data());
   outFileName.Form("plots/%s_massPt.root",fileName.Data());
   c2->SaveAs(outFileName.Data());
   delete c2;
   */
   outFileName.Form("plots/PbPbData/%s_Plots.root",fileName.Data());
   TFile *f = new TFile(outFileName.Data(),"RECREATE");
   TCanvas *c1 = new TCanvas("c1","c1");
   c1->Divide(2,1);
   c1->cd(1);
   h_jetMass2->GetXaxis()->SetTitle("m^{2} raw");
   h_jetMass2->Draw();
   c1->SetLogy();
   //   outFileName.Form("plots/%s_mass2.png",fileName.Data());
   //   c1->SaveAs(outFileName.Data());
   //   outFileName.Form("plots/%s_mass2.root",fileName.Data());
   //   c1->SaveAs(outFileName.Data());
   
   //   c1->SetLogy();
   c1->cd(2);
   h2_jetMass2Pt->GetXaxis()->SetTitle("p_{T} raw");      
   h2_jetMass2Pt->GetYaxis()->SetTitle("m^{2} raw");
   h2_jetMass2Pt->GetYaxis()->SetTitleOffset(1.55);
   h2_jetMass2Pt->Draw("colz");
   h2_jetMass2Pt->SetContour(99);
   outFileName.Form("plots/PbPbData/%s_mass2Pt.png",fileName.Data());
   c1->SaveAs(outFileName.Data());
   //   outFileName.Form("plots/ppMC/%s_mass2Pt.root",fileName.Data());
   //   c1->SaveAs(outFileName.Data());

   TCanvas *c2 = new TCanvas("c2","c2");
   h_jetPt->GetXaxis()->SetTitle("m^{2} raw");
   h_jetPt->Draw();
   outFileName.Form("plots/PbPbData/%s_Pt.png",fileName.Data());
   c2->SaveAs(outFileName.Data());
   f->cd();
   h_jetPt->Write("jetPt");
   c1->Write("c1");
   c2->Write("c2");
   f->Close();
   
   delete c1;
   delete c2;

   /*
   c1 = new TCanvas();
   h_jetMass->Draw();
   c1->SetLogy();
   outFileName.Form("plots/%s_mass.png",fileName.Data());
   c1->SaveAs(outFileName.Data());
   outFileName.Form("plots/%s_mass.root",fileName.Data());
   c1->SaveAs(outFileName.Data());
   delete c1;
   */

   /*
   c1 = new TCanvas();
   h_jetGenMass->Draw();
   c1->SetLogy();
   outFileName.Form("plots/%s_genMass.png",fileName.Data());
   c1->SaveAs(outFileName.Data());
   delete c1;
   c1 = new TCanvas();
   h_jetGenMass2->Draw();
   c1->SetLogy();
   outFileName.Form("plots/%s_genMass2.png",fileName.Data());
   c1->SaveAs(outFileName.Data());
   delete c1;
   */
}
