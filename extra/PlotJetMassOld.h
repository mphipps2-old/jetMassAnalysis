//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 19 13:33:06 2018 by ROOT version 6.13/02
// from TTree tr/new tree
// found on file: jetSubstructure_data_HION9_v51_r4_pp_apr24.root
//////////////////////////////////////////////////////////

#ifndef PlotJetMass_h
#define PlotJetMass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class PlotJetMass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           jets_cent;
   Float_t         jets_weight;
   Float_t         jets_rhoCh;
   Float_t         jets_mass;
   Float_t         jets_mass2;
   Float_t         jets_massRaw;
   Float_t         jets_massRaw2;
   Float_t         jets_pt;
   Float_t         jets_rawPt;
   Float_t         jets_eta;
   Float_t         jets_y;
   Float_t         jets_phi;
   Float_t         jets_rcPt;
   Float_t         jets_sdPt;
   Float_t         jets_sdMass;
   Float_t         jets_zg;
   Float_t         jets_theta;
   Float_t         jets_spt1;
   Float_t         jets_sy1;
   Float_t         jets_sphi1;
   Float_t         jets_spt2;
   Float_t         jets_sy2;
   Float_t         jets_sphi2;
   Float_t         jets_sda1;
   Float_t         jets_sda2;
   Float_t         jets_chPtCSubt;
   Float_t         jets_chMassCSubt;
   Float_t         jets_chPtRaw;
   Float_t         jets_chMassRaw;
   Float_t         jets_chMassGm;
   Float_t         jets_chSdPt;
   Float_t         jets_chSdMass;
   Float_t         jets_chZg;
   Float_t         jets_chTheta;
   Float_t         jets_dr;
   Float_t         jets_genMass;
   Float_t         jets_genMass2;
   Float_t         jets_genPt;
   Float_t         jets_genPt2;
   Float_t         jets_genEta;
   Float_t         jets_genRap;
   Float_t         jets_genPhi;
   Float_t         jets_genRcPt;
   Float_t         jets_genSdPt;
   Float_t         jets_genSdMass;
   Float_t         jets_genZg;
   Float_t         jets_genTheta;
   Float_t         jets_genSpt1;
   Float_t         jets_genSy1;
   Float_t         jets_genSphi1;
   Float_t         jets_genSpt2;
   Float_t         jets_genSy2;
   Float_t         jets_genSphi2;
   Float_t         jets_genSda1;
   Float_t         jets_genSda2;
   Float_t         jets_genNch;
   Float_t         jets_genChSdPt;
   Float_t         jets_genChSdMass;
   Float_t         jets_genChZg;
   Float_t         jets_genChTheta;
   Float_t         jets_nTrkRaw;
   Float_t         jets_nTrkBkg;
   Float_t         jets_nTrkBkgNoWgt;
   Float_t         jets_chPtRcSubt;
   Float_t         jets_chMassRcSubt;
   Float_t         jets_drTrkJetBkg;
   Float_t         jets_maxTrkPt;
   Float_t         jets_fcal;
   Float_t         trkJetMass4;
   Float_t         trkJetMass6;
   Float_t         trkJetMass8;
   Float_t         trkJetMass10;
   Float_t         trkJetPt4;
   Float_t         trkJetPt6;
   Float_t         trkJetPt8;
   Float_t         trkJetPt10;
   Float_t         trkJetMassRcSub2;
   Float_t         trkJetMassRcSub4;
   Float_t         trkJetMassRcSubEV;
   Float_t         trkJetMassRcSubMV;
   Int_t           evt_run;
   Int_t           evt_lumi;
   Int_t           evt_event;

   // List of branches
   TBranch        *b_jets;   //!
   TBranch        *b_trkJetMass4;   //!
   TBranch        *b_trkJetMass6;   //!
   TBranch        *b_trkJetMass8;   //!
   TBranch        *b_trkJetMass10;   //!
   TBranch        *b_trkJetPt4;   //!
   TBranch        *b_trkJetPt6;   //!
   TBranch        *b_trkJetPt8;   //!
   TBranch        *b_trkJetPt10;   //!
   TBranch        *b_trkJetMassRcSub2;   //!
   TBranch        *b_trkJetMassRcSub4;   //!
   TBranch        *b_trkJetMassRcSubEV;   //!
   TBranch        *b_trkJetMassRcSubMV;   //!
   TBranch        *b_evt;   //!

   TString        fileName;
   PlotJetMass(TString fName);
   virtual ~PlotJetMass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PlotJetMass_cxx
PlotJetMass::PlotJetMass(TString fName) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  fileName = fName;
  TString name;
  name.Form("%s.root",fileName.Data());
  cout << " file " << name.Data() << endl;
  TFile *f = TFile::Open(name.Data());
    cout << " file success" <<  TFile::Open(name.Data()) << endl;
  TTree *tree = new TTree("tr","tr");
  f->GetObject("tr",tree);
  
  Init(tree);
}

PlotJetMass::~PlotJetMass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PlotJetMass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PlotJetMass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PlotJetMass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("jets", &jets_cent, &b_jets);
   fChain->SetBranchAddress("trkJetMass4", &trkJetMass4, &b_trkJetMass4);
   fChain->SetBranchAddress("trkJetMass6", &trkJetMass6, &b_trkJetMass6);
   fChain->SetBranchAddress("trkJetMass8", &trkJetMass8, &b_trkJetMass8);
   fChain->SetBranchAddress("trkJetMass10", &trkJetMass10, &b_trkJetMass10);
   fChain->SetBranchAddress("trkJetPt4", &trkJetPt4, &b_trkJetPt4);
   fChain->SetBranchAddress("trkJetPt6", &trkJetPt6, &b_trkJetPt6);
   fChain->SetBranchAddress("trkJetPt8", &trkJetPt8, &b_trkJetPt8);
   fChain->SetBranchAddress("trkJetPt10", &trkJetPt10, &b_trkJetPt10);
   fChain->SetBranchAddress("trkJetMassRcSub2", &trkJetMassRcSub2, &b_trkJetMassRcSub2);
   fChain->SetBranchAddress("trkJetMassRcSub4", &trkJetMassRcSub4, &b_trkJetMassRcSub4);
   fChain->SetBranchAddress("trkJetMassRcSubEV", &trkJetMassRcSubEV, &b_trkJetMassRcSubEV);
   fChain->SetBranchAddress("trkJetMassRcSubMV", &trkJetMassRcSubMV, &b_trkJetMassRcSubMV);
   fChain->SetBranchAddress("evt", &evt_run, &b_evt);
   Notify();
}

Bool_t PlotJetMass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PlotJetMass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PlotJetMass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PlotJetMass_cxx
