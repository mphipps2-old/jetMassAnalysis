//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Feb 16 11:40:12 2019 by ROOT version 6.12/06
// from TTree tr/new tree
// found on file: jetSubstructure_data_HION9_v51_r4_pp_apr24.root
//////////////////////////////////////////////////////////

#ifndef tashka_h
#define tashka_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class tashka {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         jets_cent;
   Float_t         jets_weight;
   Float_t         jets_massUnsub;
   Float_t         jets_mass2Unsub;
   Float_t         jets_ptUnsub;
   Float_t         jets_massSub;
   Float_t         jets_mass2Sub;
   Float_t         jets_ptSub;
   Float_t         jets_massCalib;
   Float_t         jets_mass2Calib;
   Float_t         jets_massCalibDirect;
   Float_t         jets_ptCalib;
   Float_t         jets_etaCalib;
   Float_t         jets_yCalib;
   Float_t         jets_phiCalib;
   Float_t         jets_massConstitAlt;
   Float_t         jets_mass2ConstitAlt;
   Float_t         jets_massConstitRaw;
   Float_t         jets_mass2ConstitRaw;
   Float_t         jets_genMass;
   Float_t         jets_genMass2;
   Float_t         jets_genPt;
   Float_t         jets_genPt2;
   Float_t         jets_genEta;
   Float_t         jets_genRap;
   Float_t         jets_genPhi;
   Float_t         jets_fcal;
   Float_t         jets_chPtRcSubt;
   Float_t         jets_chMassRcSubt;
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

   tashka();
   virtual ~tashka();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tashka_cxx
tashka::tashka() : fChain(0) 
{;}

tashka::~tashka()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tashka::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tashka::LoadTree(Long64_t entry)
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

void tashka::Init(TTree *tree)
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

Bool_t tashka::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tashka::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tashka::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tashka_cxx
