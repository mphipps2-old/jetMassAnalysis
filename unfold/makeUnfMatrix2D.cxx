double fracStst=.5;

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
#include <vector>
using std::cout;
using std::endl;

#include <TRandom3.h>
#include "TH1D.h"
#include "TH2D.h"

#include "/home/mike/Desktop/JetAnalysis/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "/home/mike/Desktop/JetAnalysis/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"

#include "../getSdHists.C"
//#include "../ntupleDefinition.h"
#include "../ntupleDefinition.h"   
#include "../commonUtility.h"
#include "../jzWeight.h"
#include "../ntuples/tashka.C"
#endif

//==============================================================================
// Global definitions
//==============================================================================

#include "../JssUtils.h"
#include <TPaletteAxis.h>
#include "unfoldingUtil.h"
#include "systematicsTool.h"



bool useFullMC = false;

RooUnfoldResponse* getResponse( int kSample = kPP, int icent = 0, int optX=77, int optY=772, TH2D* hTruth=0, TH2D* hReco=0, TH1D* hRecoX=0, TH1D* hRecoY=0, TH1D* hRecoXNoWeight=0, TH1D* hRecoXNoReweight = 0, TH2D* respX=0, TH2D* respY=0, int etaBin =1,  bool doReweight = true, int nSys = -1);

// all centrality bins looped through
// doReweight grabs reweighting histogram produced from getMcWeights.C and gets bin specific reweighting factor -- reweighting MC spectra compared to data
void makeUnfMatrix2D(int kSample = kPbPb, int optX =77, int optY=772,int etaBin =0 , bool doReweight=true, int nSys=-1) {
  cout << "sample " << kSample << " optX " << optX << " optY " << optY << " etaBin " << etaBin << " doReweight " << doReweight << endl;  
  if ( nSys < 0 )  
    cout << "Nominal mode " << endl;
  else if ( (nSys >= 0 ) && ( nSys <= 21 ) )
    cout << "pp intrinsic JES sys mode" << endl;
  else if ( (nSys >= 100 ) && ( nSys <= 106 ) )
    cout << "HI JES/JMR sys mode" << endl;
  else if ( (nSys >= 200 ) && ( nSys <= 250 ) ) 
    cout << "HI JMS/JMR sys mode" << endl;
  else if  (nSys == 300 )
    cout << "Jet mass Calibration" << endl;
  else if  (nSys == 401 )
    cout << "Do both 300 and 101" << endl;
  else {
    cout << "Invald nSys option : " << nSys << endl;
    return ;
  }
  
  
  
  
  
  TH1::SetDefaultSumw2();
  int nXbins;
  double xBin[30];
  getXbin(nXbins, xBin, optX);
  cout << " nXbins = " << nXbins << endl; 
  cout << " xBin = " << xBin[0] << ",   " << xBin[1] << ",   " <<xBin[2] << ", ..." <<endl;

  
  int nYbins ;
  double yBin[30] ;
  getYbin(nYbins, yBin, optY);
    cout << " nYbins = " << nYbins << endl; 
  TH2D* hTruthTemp = new TH2D("hTruth","",nXbins,xBin,nYbins,yBin);
  TH2D* hRecoTemp = (TH2D*)hTruthTemp->Clone("hReco");
  RooUnfoldResponse* res[7];
  
  TH2D* hTruth[7];
  TH2D* hReco[7];
  TH1D* hRecoX;
  TH1D* hRecoY;
  TH2D* hResultMc[7];
  
  TH2D* hMatrix[7]; // 4-d matrix unrolled to 2-d
  
  TH2D* hResX[7]; // response matrix for pT ( mass integrated)
  TH2D* hResY[7]; // response matrix for mass ( pT integrated)
  TH2D* h2dtempM = new TH2D("h2dtemp",";Truth (m/p_{T})^{2};Reco (m/p_{T})^{2}",100,0,0.16,100,0,0.16);
  
  hRecoY = new TH1D("hRecoMass2Pt2", "Reco Mass^{2}/p_{T}^{2} (GeV/c)",nXbins,xBin);
  hRecoX = new TH1D("hRecoPt", "p_{T} (GeV/c)",nYbins,yBin);
  hRecoXNoWeight = new TH1D("hRecoPt_No_Weight", "p_{T} (GeV/c): No Weight",nYbins,yBin);
  hRecoXNoReweight = new TH1D("hRecoPt_NoReweight", "p_{T} (GeV/c): No Reweight",nYbins,yBin);
  
  for ( int i=0 ; i<=6; i++) {
    int icent = i;
    if ( !selectedCent(icent))       continue;
    if ( (kSample == kPP) && ( icent != 0 ) )      continue;

    hTruth[i] = (TH2D*)hTruthTemp->Clone(Form("%s_icent%d",hTruthTemp->GetName(),i));
    hReco[i] = (TH2D*)hRecoTemp->Clone(Form("%s_icent%d",hRecoTemp->GetName(),i));

    (TH2D*)hRecoTemp->Clone(Form("%s_icent%d",hRecoTemp->GetName(),i));
    hTruth[i]->Reset();
    hReco[i]->Reset();
    
    hResX[i] = new TH2D(Form("hResPt_icent%d",icent), ";Truth p_{T} (GeV/c);Reco p_{T} (GeV/c)",nXbins,xBin,nXbins,xBin);
    hResY[i] = new TH2D(Form("hResM_icent%d",icent), ";Truth (m/p_{T})^{2};Reco (m/p_{T})^{2}",nYbins,yBin,nYbins,yBin);

    res[i] = getResponse(kSample, i, optX, optY, hTruth[i], hReco[i], hRecoX, hRecoY, hRecoXNoWeight, hRecoXNoReweight, hResX[i], hResY[i], etaBin,doReweight, nSys);
    
    TCanvas* c01 = new TCanvas("c01", "",600,500);
    hMatrix[i] = (TH2D*)res[i]->Hresponse();
    hMatrix[i]->SetXTitle("Bin # of reco (p_{T}, (m/p_{T})^{2})");
    hMatrix[i]->SetYTitle("Bin # of truth (p_{T}, (m/p_{T})^{2})");
    hMatrix[i]->SetTitle("MC Response matrix");
    c01->SetRightMargin(0.2);
    hMatrix[i]->Draw("colz");
    ATLASLabel(0.2,0.9,"Internal",0.05,0.14);

    c01->SaveAs(Form("pdfs/correlation_2dUnf_coll%d_cent%d_etaBin%d_doReweight%d.pdf",kSample,i,etaBin,doReweight));
    TCanvas* c02 = new TCanvas("c02","",1000,500);
    c02->Divide(2,1);
    c02->cd(1);
    hResX[i]->SetNdivisions(505,"X");
    hResX[i]->Draw("colz");
    c02->cd(1)->SetRightMargin(0.2);
    gPad->SetLogz();
    ATLASLabel(0.18,0.9,"Internal",0.05,0.17);

    c02->cd(2);
    h2dtempM->SetNdivisions(505,"X");
    h2dtempM->Draw();
    hResY[i]->Draw("colz same");
    c02->cd(2)->SetRightMargin(0.2);
    gPad->SetLogz();
    ATLASLabel(0.18,0.9,"Internal",0.05,0.17);

    c02->SaveAs(Form("pdfs/PtMassResp_coll%d_cent%d_etaBin%d_doReweight%d.pdf",kSample,i,etaBin,doReweight));
    
  }
  
  TString foutname ;
  if ( nSys < 0 )   
    foutname = Form("spectraFiles/unfoldingMatrix2D_coll%d_optX%d_optY%d_etaBin%d_doReweight%d.root",kSample,optX,optY,etaBin,(int)doReweight);
  else 
    foutname = Form("spectraFiles/sys/unfoldingMatrix2D_coll%d_optX%d_optY%d_etaBin%d_doReweight%d_sys%d.root",kSample,optX,optY,etaBin,(int)doReweight,nSys);

  TFile* fout = new TFile(foutname.Data(),"update");
  for ( int i=0 ; i<=6; i++) {
    int icent = i;
    if ( !selectedCent(icent))      continue;
    if ( (kSample == kPP) && ( icent !=0 ) )      continue;
    res[i]->Write("",TObject::kOverwrite);
    hMatrix[i]->Write("",TObject::kOverwrite);
    hReco[i]->Write("",TObject::kOverwrite);
    hTruth[i]->Write("",TObject::kOverwrite);
  }
  hRecoX->Write("",TObject::kOverwrite);
  hRecoY->Write("",TObject::kOverwrite);
  fout->Close();
}

RooUnfoldResponse* getResponse(int kSample,  int icent,  int optX, int optY, TH2D* hTruth, TH2D* hReco, TH1D* hRecoX, TH1D* hRecoY, TH1D* hRecoXNoWeight, TH1D* hRecoXNoReweight, TH2D* respX, TH2D* respY, int etaBin, bool doReweight, int nSys)
{

  std::cout << " eta binnnnn " << etaBin <<  std::endl;
  RtrkProvider rtrkProv; 
  rtrkProv.Setup(kSample, icent);

  TRandom3 genRandom;
  genRandom.SetSeed(200);

  
  TH1::SetDefaultSumw2();
  TString jz2;
  TString jz3;
  TString jz4;
  if ( kSample == kPbPb ) {
    if ( (nSys>=0) && (nSys<200) ) {
      jz2 = jz2PbPbStringSys;
      jz3 = jz3PbPbStringSys;
      jz4 = jz4PbPbStringSys;
    }
    else { 
      jz2 = jz2PbPbString;
      jz3 = jz3PbPbString;
      jz4 = jz4PbPbString;
    }
  }
  else if ( kSample == kPP ) {
    if ( (nSys>=0) && (nSys<200) ) {
      jz2 = jz2PPStringSys;
      jz3 = jz3PPStringSys;
      jz4 = jz4PPStringSys;
    }
    else {
      jz2 = jz2PPString;
      jz3 = jz3PPString;
      jz4 = jz4PPString;
    }
  }
  
  TH1D* hFcalReweight;
  if ( kSample == kPbPb ) {
    TFile* fcal = new TFile("reweightFactors/FCal_HP_v_MB_weights.root");
    hFcalReweight = (TH1D*)fcal->Get("weight_MBov_to_HP");
  }

  
  TH2D* hReweight;
  if ( doReweight ) {
    hReweight = getRewTable(kSample, icent, etaBin);
  }

  TF1* fjmscal[30];
  if( (nSys == 300)||( nSys==401)) {
    TFile* fin = new TFile(Form("fJMScalibration_kSample%d_icent%d_num.root",kSample,icent));
    for ( int ix = 5 ; ix<=11 ; ix++) {
      fjmscal[ix] = (TF1*)fin->Get(Form("f1_kSample%d_icent%d_ix%d",kSample,icent,ix));
    }
  }
  int nXbinsCal;
  double xBinCal[30];
  getXbin(nXbinsCal, xBinCal, optX);
  TH1D* xBinTempCal = new TH1D("xBinTempCal","", nXbinsCal, xBinCal);



  
  //  TFile* checkEntries = new TFile(Form("checkEntry/entries_kSample%d_icent%d_optX%d_optY%d.root",kSample,icent,optX,optY));
  //  TH2D* recoEntries_jz2 = (TH2D*)checkEntries->Get("reco_jz2");
  //  TH2D* recoEntries_jz3 = (TH2D*)checkEntries->Get("reco_jz3");
  //  TH2D* recoEntries_jz4 = (TH2D*)checkEntries->Get("reco_jz4");
  
  //  jetSubStr  myJetMc;
  //  TBranch  *b_myJetSubMc;

  float ptSys;
  TBranch *b_ptSys;
  TString jetSysName = getPtSysName(nSys);
  if ( nSys == 401) jetSysName = getPtSysName(101);


  cout << " Setting tree branch address..." << endl;
  TFile* fjz2 = new TFile(Form("../ntuples/%s",jz2.Data()));
  TTree* tr2 = (TTree*)fjz2->Get("tr");
  tashka t_jz2;
  //  tr2->SetBranchAddress("jets", &(tr.jets_cent), &b_myJetSubMc);
  if ( (nSys>=0) && (nSys<200) )
    tr2->SetBranchAddress(jetSysName.Data(), &ptSys, &b_ptSys);
  else if (nSys==401)
    tr2->SetBranchAddress(jetSysName.Data(), &ptSys, &b_ptSys);


  TFile* fjz3 = new TFile(Form("../ntuples/%s",jz3.Data()));
  TTree* tr3 = (TTree*)fjz3->Get("tr");
  tashka t_jz3;
  //  tr3->SetBranchAddress("jets", &(tr.jets_cent), &b_myJetSubMc);
  if ( (nSys>=0) && (nSys<200) )
    tr3->SetBranchAddress(jetSysName.Data(), &ptSys, &b_ptSys);
  else if (nSys==401)
    tr3->SetBranchAddress(jetSysName.Data(), &ptSys, &b_ptSys);

  TFile* fjz4 = new TFile(Form("../ntuples/%s",jz4.Data()));
  TTree* tr4 = (TTree*)fjz4->Get("tr");
  tashka t_jz4;
  //  tr4->SetBranchAddress("jets", &(tr.jets_cent), &b_myJetSubMc);
  if ( (nSys>=0) && (nSys<200) )  
    tr4->SetBranchAddress(jetSysName.Data(), &ptSys, &b_ptSys);
  else if (nSys==401)
    tr4->SetBranchAddress(jetSysName.Data(), &ptSys, &b_ptSys);


  int nXbins;
  double xBin[30];
  getXbin(nXbins, xBin, optX);

  int nYbins ;
  double yBin[30] ;
  getYbin(nYbins, yBin, optY);

  TH1D* xBinTemp = new TH1D("xBinTemp","", nXbins, xBin);
  TH1D* yBinTemp = new TH1D("yBinTemp","", nYbins, yBin);

  std::vector<std::vector<TH1D*> > hJES1d; // centrailty // gen pT // mass bin

  TH2D* hJES; // centrailty // gen pT // 
  TH2D* hJER; // centrailty // gen pT // 
  hJES = new TH2D("hJES",";Truth p_{T} (GeV/c); Truth (m/p_{T})^{2}", nXbins,xBin,nYbins,yBin);
  hJER = new TH2D("hJER",";Truth p_{T} (GeV/c); Truth (m/p_{T})^{2}", nXbins,xBin,nYbins,yBin);

  for ( int ix = 0 ; ix< nXbins; ix++) {
    std::vector<TH1D*> tempV;
    for ( int iy = 0 ; iy< nYbins; iy++) { 
      if ( (kSample == kPP) && ( icent != 0 ) )      continue;
      TH1D* tempH = new TH1D(Form("hjes1d_%d_%d",ix,iy),Form("hjes1d_%d_%d",ix,iy),100,0,2);
      tempV.push_back(tempH);
      //      hJES1d[ix][iy]->push_back(new TH1D(Form("hjes1d_%d_%d",ix,iy),Form("hjes1d_%d_%d",ix,iy),100,0,2));
    }
    hJES1d.push_back(tempV);
  }
  std::cout <<" vector1 size " << hJES1d[0].size() << " vector2 size " << hJES1d.size() << std::endl;   

  RooUnfoldResponse* res;
  res = new RooUnfoldResponse( hReco, hTruth );
  res->SetName(Form("responseMatrix_icent%d",icent));

  long long entries;
  tashka *tr = new tashka();  
  for ( int ijz =2 ; ijz<=4 ; ijz++) {
  //  int ijz =4 ;
  //  for ( int ijz =4 ; ijz<=4 ; ijz++) {

    //    TH2D* hRecoEntries;
    double jzNorm=0;
    if ( ijz==2)  {
      *tr = t_jz2;
      tr->Init(tr2);
      entries = tr2->GetEntries();
      jzNorm = hi9EvtWgtJZ2; 
      //      hRecoEntries = recoEntries_jz2;
    }
    else if ( ijz==3)  {
      *tr = t_jz3;
      tr->Init(tr3);
      entries = tr3->GetEntries();
      jzNorm = hi9EvtWgtJZ3; 
      //      hRecoEntries = recoEntries_jz3;
    }
    else if ( ijz==4)  {
      *tr = t_jz4;
      tr->Init(tr4);
      entries = tr4->GetEntries();
      jzNorm = hi9EvtWgtJZ4; 
      //      hRecoEntries = recoEntries_jz4;
    }
    cout << "Scanning JZ"<<ijz<<" file.  Total events = " << entries << endl;
    
    for (Int_t i= 0; i<entries ; i++) {
      if ( i > entries * fracStst ) continue;
      
      tr->GetEntry(i);
      
      //      if ( (!useFullMC) && (i%2 == 0) )
      //continue;
      double recoVarX=0.; double truthVarX=0.;
      double recoVarY=0.; double truthVarY=0.;
      
      getXvalues( recoVarX, truthVarX, tr->jets_genPt, tr->jets_pt, optX); 
      getYvalues( recoVarY, truthVarY, tr->jets_genMass2, tr->jets_massRaw2, tr->jets_genPt, tr->jets_pt, tr->jets_chMassRcSubt, tr->jets_chPtRcSubt, optY);

      
      if ( (nSys>=0) && (nSys<200) )   {
	double extraPtScale = ptSys / tr->jets_pt ; 
	recoVarX = recoVarX * extraPtScale ; //pt 
	tr->jets_pt = ptSys;  // New pT!!! 
	tr->jets_massRaw2 = tr->jets_massRaw2 * extraPtScale ; // new mass so that m/pT is invariant.  This step is necessary for reweighitng (pT, m/pT)
      }
      if ( nSys==401) {
	double extraPtScale = ptSys / tr->jets_pt ;
        recoVarX = recoVarX * extraPtScale ; //pt
        tr->jets_pt = ptSys;  // New pT!!!
        tr->jets_massRaw2 = tr->jets_massRaw2 * extraPtScale ;
      }
      else if (nSys==200) { // JMR 
	// smear by 20% the recoY 
	double theCenter = truthVarY * getJMSscale( kSample, icent, tr->jets_pt);
        double recoDev = recoVarY - theCenter;
	double jmrUnc = getJmrUnc( kSample, icent, tr->jets_pt); 
	double theVariation = sqrt ( (1 +jmrUnc)*(1+jmrUnc) - 1 );  
        double theResol = getJMRsigma( kSample, icent, tr->jets_pt);
	recoVarY = theCenter + recoDev * genRandom.Gaus(1, theVariation * theResol);  //20 percent
	//	recoVarY = theCenter + recoDev * genRandom.Gaus(1, 0.66 * theResol);  //20 percent
      }	
      else if (nSys==201) { // JMR  HI
	double theCenter = truthVarY * getJMSscale( kSample, icent, tr->jets_pt);
        double recoDev = recoVarY - theCenter;
	double jmrUnc = getJmrUncHI( kSample, icent, tr->jets_pt); 
	double theVariation = sqrt ( (1 +jmrUnc)*(1+jmrUnc) - 1 );  
        double theResol = getJMRsigma( kSample, icent, tr->jets_pt);
	recoVarY = theCenter + recoDev * genRandom.Gaus(1, theVariation * theResol);  //20 percent
	//	recoVarY = theCenter + recoDev * genRandom.Gaus(1, 0.66 * theResol);  //20 percent
      }	

      else if (nSys==211) { // JMS
	recoVarY = recoVarY * 1.014;
      }
      else if (nSys==217) { // JMS by Herwig
        recoVarY = recoVarY * 1.05;
      }

      if ( passEvent(tr->jets_cent, tr->jets_genPt, tr->jets_pt,icent, true) == false ) // true = isMC
	continue;
      if ( !passEtaCut(tr->jets_eta, etaBin) )
	continue;
      
      double fcalWeight = 1.0; 
      if ( kSample==kPbPb) {
	fcalWeight = hFcalReweight->GetBinContent(hFcalReweight->GetXaxis()->FindBin(tr->jets_fcal));
      }
      
      // Data/MC reweighting factors 
      double rewFact = 1;
      if ( doReweight) { 
	int rewBin = hReweight->FindBin(tr->jets_pt, (tr->jets_massRaw2) / (tr->jets_pt * tr->jets_pt));
	rewFact = hReweight->GetBinContent(rewBin);
      }
      if (  (nSys==300) || (nSys == 401)) { // JMS calibration
	int xx = xBinTempCal->FindBin( tr->jets_pt);
	if ( xx > 11 )  xx = 11;
	if ( xx < 5 )  xx = 5;
	float mptVal = recoVarY;
	if (mptVal <0 )  mptVal = 0;
	double theFac = 1;
        if ( recoVarY < 0 )
          theFac = 1;
        else if ( recoVarY > 0.25 )
          theFac = fjmscal[xx]->Eval(0.25);
	else
          theFac = fjmscal[xx]->Eval(recoVarY);

	recoVarY = recoVarY / theFac;
      }
      //      cout << " respX truthVarX " << truthVarX << " recoVarX " << recoVarX << " jet weight " << tr->jets_weight << " rewFact " << rewFact << " jzNorm " <<  jzNorm << " fcalWEight " << fcalWeight << endl;

      hTruth->Fill(truthVarX, truthVarY, tr->jets_weight * rewFact * jzNorm * fcalWeight);
      hReco->Fill(recoVarX, recoVarY, tr->jets_weight * rewFact * jzNorm * fcalWeight);
      // pt
      hRecoX->Fill(recoVarX, tr->jets_weight * rewFact * jzNorm * fcalWeight);
            // pt
      hRecoXNoWeight->Fill(recoVarX);
      hRecoXNoReweight->Fill(recoVarX, tr->jets_weight * rewFact * jzNorm * fcalWeight);
      // m2/pt2
      hRecoY->Fill(recoVarY, tr->jets_weight * rewFact * jzNorm * fcalWeight);
      
      //      res->Fill(  recoVarX, recoVarY, truthVarX, truthVarY, tr->jets_weight * rewFact * jzNorm * fcalWeight);

      respX->Fill( truthVarX, recoVarX,  tr->jets_weight * rewFact * jzNorm* fcalWeight);
      respY->Fill( truthVarY, recoVarY,  tr->jets_weight * rewFact * jzNorm* fcalWeight);
    
      //        cout << " respX truthVarY " << truthVarY << " recoVarY " << recoVarY << " val " << tr->jets_weight * rewFact * jzNorm* fcalWeight << " fcalweight " << fcalWeight << endl;
      
      int ix = xBinTemp->FindBin(truthVarX);
      int iy = yBinTemp->FindBin(truthVarY);

      if ( (ix >=0) && (ix<nXbins) && (iy >=0) && (iy<nYbins) )   {
	 (hJES1d[ix][iy])->Fill(recoVarX/truthVarX);	
      }
    }
    //    hJES1d[1][1]->Print("ALL");
  
  }
  /*
  std::cout << " about to printt " << hJES1d[1][1] << std::endl;
  hJES1d[1][1]->Print("ALL");
  TCanvas *c1 = new TCanvas();
  std::cout << " about to printttt555555 " << hJES1d[1][1] << std::endl;
  hJES1d[1][1]->Draw();
  std::cout << " here "  << std::endl;
  c1->SaveAs("tempx5y5.root");
  TCanvas *c2 = new TCanvas();
  std::cout << " here eeeeeeeee"  << std::endl;
  hJES1d[0][0]->Draw();
  c2->SaveAs("tempx10y10.root");
  */

      
  for ( int ix = 1 ; ix<= nXbins; ix++) {
    for ( int iy = 1 ; iy<= nYbins; iy++) {
      if( (kSample == kPP) && ( icent != 0 ) )      continue;
      //      cout << " kSample " << kSample << " icent " << icent << endl;
      //      cout << " ix " << ix << " iy "<< iy << " mean " << " nxbins " << nXbins << " nybins " << nYbins << endl;
      //      cout << "hist " << hJES1d[ix-1][iy-1] << endl;
      //      TCanvas *c1 = new TCanvas();
      //      hJES1d[ix-1][iy-1]->Draw();
      //      c1->SaveAs("temp.root");
      //      cout << " mean " << hJES1d[ix-1][iy-1]->GetMean() << endl;
      //      cout << " rms " << hJES1d[ix-1][iy-1]->GetRMS() << endl;
      hJES->SetBinContent(ix, iy, hJES1d[ix-1][iy-1]->GetMean());
      hJER->SetBinContent(ix, iy, hJES1d[ix-1][iy-1]->GetRMS());
    }
  }
  //  cJes2->cd(2)->Update();
  
  TCanvas* cJESJER = new TCanvas("cjesjer","",1000,500);
  cJESJER->Divide(2,1);
  cJESJER->cd(1);
  TH2D* htemp2 = new TH2D("htemp2",";Truth p_{T};Truth (m/p_{T})^{2};",100,126,400,100,0,0.16);
  //  handsomeTH1(htemp2,1);
  htemp2->GetXaxis()->SetMoreLogLabels();
  htemp2->SetAxisRange(0.7,1.3,"Z");
  htemp2->DrawCopy("colz");
  hJES->SetAxisRange(0.7,1.3,"Z");
  hJES->Draw("colz same");
  drawCentrality(kSample,icent,0.55,0.87,1,20);
  drawText("Jet energy scale",0.5,0.8,1,20);
  ATLASLabel(0.2,0.87,"Internal",0.05,0.17);
  cJESJER->cd(2);
  htemp2->SetAxisRange(0,0.3,"Z");
  htemp2->DrawCopy("colz");
  hJER->SetAxisRange(0,0.3,"Z");
  hJER->Draw("colz same");
  drawCentrality(kSample,icent,0.55,0.87,1,20);
  drawText("Jet energy resolution",0.4,0.8,1,20);
  ATLASLabel(0.2,0.87,"Internal",0.05,0.17);

  cJESJER->cd(1)->Update();
  TPaletteAxis *palette1 = (TPaletteAxis*)hJES->GetListOfFunctions()->FindObject("palette");
  palette1->SetX1NDC(0.87);
  palette1->SetX2NDC(0.92);
  gPad->SetLogx();
  cJESJER->cd(2)->Update();
  TPaletteAxis *palette2 = (TPaletteAxis*)hJER->GetListOfFunctions()->FindObject("palette");
  palette2->SetX1NDC(0.87);
  palette2->SetX2NDC(0.92);
  gPad->SetLogx();
  
  cJESJER->SaveAs(Form("pdfs/JESJER_kSample%d_icent%d_etaBin%d_optX%d_optY%d_applyMDJ0.pdf",kSample,icent,etaBin,optX,optY));
  delete tr;
  return res;
}


#ifndef __CINT__
int main (int argc, char** argv) {
  int kSample    = atoi(argv[1]);
  int optX       = atoi(argv[2]);
  int optY       = atoi(argv[3]);
  int etaBin     = atoi(argv[4]);
  int doReweight = atoi(argv[5]);
  int nSys       = atoi(argv[6]);
  std::cout << "kSample " << kSample << " optX " << optX << " optY " << optY << " etaBin " << etaBin << " doReweight " << doReweight << " nSys " << nSys << std::endl;  
  // do reweight? no centrality? optx and opty?
  makeUnfMatrix2D(kSample, optX, optY, etaBin, doReweight, nSys);
  return 0;
}  // Main program when run stand-alone
#endif


