#ifndef UNFOLDINGUTIL_H
#define UNFOLDINGUTIL_H

bool selectedCent(int icent=0) {
  if ( icent ==0 )  return true;
  if ( icent ==1 )  return true;
  if ( icent ==2 )  return true;
  if ( icent ==3 )  return true;
  if ( icent ==4 )  return true;
  if ( icent ==5 )  return true;
  if ( icent ==6 )  return true;
  if ( icent ==7 )  return true;
  else return false; 
  return 0;
}



int findMbin(double mpt=0) {
  if ( mpt < 0.06 ) return 1;
  else if ( mpt < 0.1 ) return 2;
  else if ( mpt < 0.15 ) return 3;
  else 
    return 4;
}



TH2D* getRewTable(int kSample, int icent, int etaBin)  { 
  TFile* fReweight = new TFile(Form("reweightFactors/reweightingFactor_kSample%d_etaBin%d_factorized_v60.root",kSample,etaBin)); // reweightingFactor_etaBin1_factorized_v60.root
  TH2D* hTemp = (TH2D*)fReweight->Get(Form("factorizedRatio2_kSample%d_icent%d",kSample,icent));  // 00 default
  return hTemp;
}
/*
  TString jz2PbPbString = "jetSubstructure_MC_HION9_jz2_v4.7_v4_Jan23_ptCut90Eta2.1.root";
  TString jz3PbPbString = "jetSubstructure_MC_HION9_jz3_v4.7_v4_Jan23_ptCut90Eta2.1.root";
  TString jz4PbPbString = "jetSubstructure_MC_HION9_jz4_v4.7_v4_Jan23_ptCut90Eta2.1.root";
  TString   jz2PPString = "jetSubstructure_MC_HION9_jz2_v4.7_r4_pp_Jan23_ptCut90Eta2.1.root";
  TString   jz3PPString = "jetSubstructure_MC_HION9_jz3_v4.7_r4_pp_Jan23_ptCut90Eta2.1.root";
  TString   jz4PPString = "jetSubstructure_MC_HION9_jz4_v4.7_r4_pp_Jan23_ptCut90Eta2.1.root";
*/
// systematics

TString jz2PbPbStringSys = "jetSubstructure_pbpbMC_HION9_jz2_v52_april24.root";
TString jz3PbPbStringSys = "jetSubstructure_pbpbMC_HION9_jz3_v52_april24.root";
TString jz4PbPbStringSys = "jetSubstructure_pbpbMC_HION9_jz4_v52_april24.root";

TString jz2PPStringSys = "jetSubstructure_ppMC_HION9_jz2_v60_april24.root";     // _april14 is added only for pp as ptSysHI is fixed
TString jz3PPStringSys = "jetSubstructure_ppMC_HION9_jz3_v60_april24.root";     // _april14 is added only for pp as ptSysHI is fixed
TString jz4PPStringSys = "jetSubstructure_ppMC_HION9_jz4_v60_april24.root";     // _april14 is added only for pp as ptSysHI is fixed

TString jz2PbPbString = "jetSubstructure_pbpbMC_HION9_jz2_v52_april24.root";
TString jz3PbPbString = "jetSubstructure_pbpbMC_HION9_jz3_v52_april24.root";
TString jz4PbPbString = "jetSubstructure_pbpbMC_HION9_jz4_v52_april24.root";

TString jz2PPString = "jetSubstructure_ppMC_HION9_jz2_v60_april24.root";
TString jz3PPString = "jetSubstructure_ppMC_HION9_jz3_v60_april24.root";
TString jz4PPString = "jetSubstructure_ppMC_HION9_jz4_v60_april24.root";

TString ppDataString = "jetSubstructure_data_HION9_v51_r4_pp_apr24.root";
TString pbpbDataString = "jetSubstructure_data_HION9_v51_r4_pbpb_apr24.root";


TString getPtSysName( int nSys) {
  if (nSys < 0) 
    return "";
  else if ( (nSys >=0) && (nSys <= 21) )    // 0 -21 : pp intrinsic
    return Form("ptSysPP%d",nSys);
  else if ( (nSys >=100) && (nSys <= 106) )    // 100 -106 : HI intrinsic
    return Form("ptSysHI%d",nSys-100);
  else 
    return "";
}

int getRefIter( int kSample=0, int icent=0) {

  return 6;

  if ( kSample == 0 ) return 6;
  if ( kSample == 1 ) {
    if ( icent == 0 ) return 6;
    if ( icent == 1 ) return 4;
    if ( icent == 2 ) return 4;
    if ( icent == 3 ) return 4;
    if ( icent == 4 ) return 3;
    if ( icent == 5 ) return 3;
    if ( icent == 6 ) return 3;
    if ( icent == 7 ) return 3;
  }
  
  return -1;
}



void getXbin(int &nBins, double* xBin, int optX) {
  if ( optX == 1 ) {
    //    nBins = 13;  // default
    //    double ptBin[14]={20,40,63.096, 82., 100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  500.,  630.944, 999.970};
    nBins = 12;  // default
    double ptBin[13]={20,40,63.096, 82., 100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  500.,  630.944};
    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
    }
  }
  else if ( optX == 2) {
    nBins = 8;  // default
    double ptBin[9]={100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  500.,  1000};
    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
    }
  }

  else if ( optX == 3) {
    nBins = 8;  // default
    double ptBin[9]={100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  500.,  1000};
    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
    }
  }

  else if ( optX == 21) { //  21 is reserved for calibration 
    nBins = 6;  // default
    double ptBin[7]={100.000, 126, 158, 200,  251,  316.224,  500};
    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
    }
  }
  else if ( optX == 210) { //  21, etabin =2 
    nBins = 6;  // default
    double ptBin[7]={100.000, 126, 158, 200,  251,  316.224, 500};
    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
    }
  }
  else if ( optX == 211) { //  21, etabin =2 
    nBins = 6;  // default
    double ptBin[7]={100.000, 126, 158, 200,  251,  316.224, 500};
    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
    }
  }
  else if ( optX == 212) { //  21, etabin =2 
    nBins = 5;  // default
    double ptBin[6]={100.000, 126, 158, 200,  251,  316.224};
    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
    }
  }
  else if ( optX == 77) {// reweighting default (Mike)
    nBins = 20;  // default
    //nBins = 19;  // default
    //nBins = 24;  // default

    // change bin spacing to logarithmic
    // what I was originally using
    //    double ptBin[14]={100.000, 112.202, 125.893, 141.254, 158.489, 177.830, 199.526, 223.872, 251.189, 281.838, 316.228,  398.107,  501.187, 1000.};
    // new fully logarithmic binning with 14 bins from 100-1000 GeV
    //    double ptBin[15]={100.000,117.877,138.9495,163.7894,193.0698,227.5847,268.2696,316.2278,372.75937,439.3969,517.9478,610.5342,719.6857,848.3401, 1000.};
    //    double ptBin[16]={100.000,114.8698355,131.9507911,151.5716567,174.1101127,200.,229.739671,263.9015822,303.1433133,348.2202253,400.,459.479342,527.8031643,606.2866266,693.4404506,800.};
    double ptBin[21]={100.000, 112.202, 125.893, 141.254, 158.489, 177.828, 199.526, 223.872, 251.189, 281.838, 316.228, 354.813, 398.107, 446.684, 501.187,562.341, 630.957, 707.946, 794.328, 891.251, 1000.};

    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
            std::cout << " xBin " << xBin[i] << std::endl;
    }
  }
   else if ( optX == 80) {// binning for reweighting (since truth values go down to 60s/70s
    nBins = 24;  // default

    double ptBin[25]={63.0957,70.7946,79.4328,89.1251,100.000, 112.202, 125.893, 141.254, 158.489, 177.828, 199.526, 223.872, 251.189, 281.838, 316.228, 354.813, 398.107, 446.684, 501.187,562.341, 630.957, 707.946, 794.328, 891.251, 1000.};

    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
            std::cout << " xBin " << xBin[i] << std::endl;
    }
  }
   else if ( optX == 81) {// binning for reweighting (since truth values go down to 60s/70s
    nBins = 12;  // default

    double ptBin[13]={63.0957,79.4328,100.000, 125.893, 158.489, 199.526, 251.189, 316.228, 398.107, 501.187,630.957, 794.328, 1000.};
    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
           std::cout << " xBin " << xBin[i] << std::endl;
    }
  }
 
  else if ( optX == 78) {// reweighting 
    nBins = 9;  // default
    double ptBin[10]={100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  500.,  630.944, 999.970};
    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
    }
  }
  else if ( (optX == 79) || ( optX==80) )  {
    nBins = 4;
    //    double ptBin[10]={100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  500.,  630.944, 999.970};
    double ptBin[5]={100.000, 125.892, 199.525,  500.,  999.970};
    for ( int i=0 ; i<= nBins ; i++) {
      xBin[i] = ptBin[i] ;
    }
  }
  else { 
    cout << " Error in unfoldingUtil::getXbin! " << endl << " optX  = " << optX << " option is not defined in this function" << endl;
  }  
}

void getXvalues( double &recoVarX, double &truthVarX, float jets_genPt, float jets_recoPt, int optX) {
  if ( (optX == 1) || ( optX==2) || ( optX==3) ) {
    truthVarX = jets_genPt;
    recoVarX = jets_recoPt;
  }
  else if ( optX == 77 ) {
    truthVarX = jets_genPt;
    recoVarX = jets_recoPt;
  }
  else if ( optX == 79 )  {
    recoVarX =  jets_genPt;
    truthVarX = jets_recoPt;
  }
  else if ( optX == 80 )  {
    recoVarX =  jets_genPt;
    truthVarX = jets_recoPt;
  }
  else { 
    cout << " Error in unfoldingUtil::getXvalues! " << endl << " optX  = " << optX << " option is not defined in this function" << endl;
  }  
}

void getYbin(int &nBins, double* yBin, int optY) {

  if ( optY == 1) {
    //    nBins = 18;
    //    double massBin[19] = { -35,-19,-17,-15, -13,-10,0,10,13,15,17,19,21,24,28,35,50,100,200};
    nBins = 13;
    double massBin[14] = { -35,-15,0,10,15,20,25,30,40,50,70,100,150,300};
    for ( int i=0 ; i<= nBins ; i++) {
      yBin[i] = massBin[i];
    }
  }
  else if ( (optY==2) || (optY==8) ) { // m/pT 
    nBins = 9;
    double massBinLog[10] = { -10000, -2.8, -2.6, -2.4, -2.2, -2, -1.8, -1.6, -1.2, -.8};
    yBin[0] = -.5;
    for ( int i=1 ; i<= nBins ; i++) {
      yBin[i] = pow(10, massBinLog[i]/2);
    }
  }
  else if ( optY == 21) { //  21 is reserved for calibration.  Used only for Truth m/pT 
    nBins = 9;
    double massBinLog[20] = { -10000, -2.8, -2.6, -2.4, -2.2, -2, -1.8, -1.6, -1.2, -.8}; 
    yBin[0] = 0;  
    for ( int i=1 ; i<= nBins ; i++) {
      yBin[i] = pow(10, massBinLog[i]/2);
    }
  }
  
  
  else if ( optY == 3) {
    nBins = 10;
    double massBin[11] = { -0.5,-0.05,0,0.05,0.1,0.13,0.16,0.2,0.24,0.3,0.5};
    for ( int i=0 ; i<= nBins ; i++) {
      yBin[i] = massBin[i];
    }
    
  }
  else if ( optY == 4) {
    nBins = 4;
    double massBin[5] = { -0.5, 0, 0.16, 0.2, 0.5};
    
    for ( int i=0 ; i<= nBins ; i++) {
      yBin[i] = massBin[i];
    }
  }
  else if ( optY == 7)   {
    nBins = 12 ;
    for ( int i=0 ; i<= nBins ; i++) {
      yBin[i] = -2000 + (2000+2000) * float(i)/nBins ;
    }
  }
  
  else if ( optY == 77) {
    nBins = 7;
    double massBin[8] = { -50, -25,0,25,50,75,100,300};
    for ( int i=0 ; i<= nBins ; i++) {
      yBin[i] = massBin[i];
    }
  }
  else if ( optY == 771) {
    nBins = 12;
    double massBin[13] = { -50, -30,-10,0,10,20,30,40,50,75,100,200,300};
    for ( int i=0 ; i<= nBins ; i++) {
      yBin[i] = massBin[i];
    }
  }
  else if ( optY == 772) {
    // current default
    nBins = 30;
    double highY = 0.1;
    double lowY = -0.015;
    for ( int i=0 ; i<= nBins ; i++) {
      yBin[i] = lowY + (highY-lowY)*i/nBins;
      //      cout << " i " << i << " masssBin " << massBin[i] << " yBin " << yBin[i] << endl;
    }
  }
  else if ( optY == 773) {
    // m/pt default
    nBins = 7;
    double highY = 0.24;
    double lowY = 0.;
    for ( int i=0 ; i<= nBins ; i++) {
      yBin[i] = lowY + (highY-lowY)*i/nBins;
      //      cout << " i " << i << " masssBin " << massBin[i] << " yBin " << yBin[i] << endl;
    }
  }
  
  else if ( optY == 78) {
    nBins = 8;
    double massBin[9] = { 0.,0.03,0.06,0.09,0.12,0.15,0.18,0.24,0.3};
    for ( int i=0 ; i<= nBins ; i++) {
      yBin[i] = massBin[i];
    }
  }

  else if ( optY == 79) {
    nBins = 10;
    //    double massBin[11] = { -0.2,0,0.015,0.03,0.06,0.09,0.12,0.15,0.18,0.24,0.3};
    double massBin[11] = { -0.2,0,0.01, 0.02, 0.04, 0.06, 0.08, 0.12, 0.16, 0.24, 0.32};
    for ( int i=0 ; i<= nBins ; i++) {
      if ( massBin[i] > 0 )
	yBin[i] = massBin[i] * massBin[i];
      else
	yBin[i] = - massBin[i] * massBin[i];

      yBin[i] = yBin[i] + 1 ;
    }
  }
  else if ( optY == 80) {
    nBins = 8;
    //    double massBin[11] = { -0.2,0,0.015,0.03,0.06,0.09,0.12,0.15,0.18,0.24,0.3};
    //    double massBin[10] = { -0.1, 0 , 0.03, 0.06, 0.09, 0.12, 0.16, 0.2, 0.25, 0.3};
    //    double massBin[8] = { -0.1, 0 , 0.05, 0.10, 0.15, 0.2, 0.25, 0.3};
    double massBin[9] = { -0.05, 0.05, 0.08, 0.12, 0.15, 0.18, 0.24, 0.3,0.5};
    for ( int i=0 ; i<= nBins ; i++) {
      if ( massBin[i] > 0 )
	yBin[i] = massBin[i] * massBin[i];
      else
	yBin[i] = - massBin[i] * massBin[i];

      yBin[i] = sqrt(yBin[i] + 1) ;
    }
  }
  else { 
    cout << " Error in unfoldingUtil::getYbin! " << endl << " optY  = " << optY << " option is not defined in this function" << endl;
  }    
}


void getRebinMpt(int &nYbins, double* yBin, int icent, int ptBin) {

  if ( (icent ==0) && ( ptBin <= 9 ) ) {
    getYbin(nYbins, yBin, 2);
  }
  else if ( (icent ==0) && ( ptBin >= 10 ) ) {
    cout << "rebinning" << endl;
    nYbins = 8;
    double ytemp[9] = { -0.5,-0.05,0,0.05,0.1,0.13,0.2,0.3,0.5};
    for ( int i =0 ; i<=8 ; i++) { 
      yBin[i] = ytemp[i];
    }
  }
  else  {
    getYbin(nYbins, yBin, 2);
  }
  
  
}


const Double_t cutdummy= -99999.0;

void getYvalues( double &recoVarY, double &truthVarY, float jets_genMass2, float jets_recoMass2, float jets_genPt, float jets_recoPt, float jets_recoChMassRcSubt, float jets_recoChPtRcSubt, int optY) {

  double genM = jets_genMass2;
  double genM2 = jets_genMass2;
  if ( genM < 0 ) genM = 0.000001;
  double recoM = jets_recoMass2;
  double recoM2 = jets_recoMass2;
  double genPt = jets_genPt;
  double recoPt = jets_recoPt;

  double genMoverPt = genM / genPt;
  double recoMoverPt = recoM / recoPt;

  double genMoverPtSquared = genM2 / (genPt*genPt);
  double recoMoverPtSquared = recoM2 / (recoPt*recoPt);

  if (optY==1)  {
    recoVarY = recoM;
    truthVarY = genM;
  }
  else if ( (optY==2) || (optY ==3))  {
    truthVarY = genMoverPt;
    recoVarY = recoMoverPt;
  }
  else if (optY == 772)  {
    truthVarY = genMoverPtSquared;
    recoVarY = recoMoverPtSquared;
  }
  else if ( optY == 7) { // charge assisted mass
    truthVarY = genM;
    recoVarY = jets_recoChMassRcSubt * jets_recoPt / jets_recoChPtRcSubt ;
  }
  else if ( optY == 8) { // charge assisted mass
    truthVarY = genMoverPt;
    recoVarY = jets_recoChMassRcSubt / jets_recoChPtRcSubt ;
  }
  else if ( optY == 79) { 
    truthVarY = genMoverPt * genMoverPt;
    if ( recoMoverPt > 0 )
      recoVarY = recoMoverPt * recoMoverPt;
    else 
      recoVarY = - recoVarY;
    
    recoVarY = recoVarY + 1;
    truthVarY = truthVarY + 1;
  }
  else if ( optY == 80) { 
    truthVarY = genMoverPt * genMoverPt;
    if ( recoMoverPt > 0 )
      recoVarY = recoMoverPt * recoMoverPt;
    else 
      recoVarY = - recoVarY;
    
    recoVarY = recoVarY + 1;
    truthVarY = truthVarY + 1;
  }
  else { 
    cout << " Error in unfoldingUtil::getYvalues! " << endl << " optY  = " << optY << " option is not defined in this function" << endl;
  }
  
}


bool passGenEvent( int jets_cent, float jets_genPt)  {

  double ptCutGen = 60;
  double ptCutUpGen = 1000;
  /*
  if ( jets_cent != icent )
    return false;
  */
  if ( jets_genPt < ptCutGen ) 
    return false;
  
  if ( jets_genPt > ptCutUpGen ) 
    return false;

  return true;
  
}


bool passRecoEvent( int jets_cent, float jets_recoPt, float jets_eta, float jets_phi, bool isMC)  {

  double ptCut = 100;
  //  double ptCut = 50;
  double ptCutUp = 1000;
  /*
  if ( jets_cent != icent ) {
    return false;
  }
  */
  if ( jets_recoPt < ptCut ) {
    return false;
  }
  if ( jets_recoPt > ptCutUp ) {
    return false;
  }
  // cut for detector hole. Yields of final result should be scaled by ratio of total acceptance we've eliminated
  if ((jets_eta > 0 && jets_eta < 1) && (jets_phi > TMath::Pi()/4 && jets_phi < (11*TMath::Pi()/32))) {
    //    if (!isMC) cout << " returning false detector hole" << endl; 
    return false;
  }
  return true;

}

bool passReco2dEvent( int jets_cent, float jets_recoPt, float jets_eta, float jets_phi, bool isMC)  {

  //  double ptCut = 0;
  double ptCut = 50;
  double ptCutUp = 1000;
  /*
  if ( jets_cent != icent ) {
    return false;
  }
  */
  if ( jets_recoPt < ptCut ) {
    return false;
  }
  if ( jets_recoPt > ptCutUp ) {
    return false;
  }
  // cut for detector hole. Yields of final result should be scaled by ratio of total acceptance we've eliminated
  if ((jets_eta > 0 && jets_eta < 1) && (jets_phi > TMath::Pi()/4 && jets_phi < (11*TMath::Pi()/32))) {
    //    if (!isMC) cout << " returning false detector hole" << endl; 
    return false;
  }
  return true;

}


bool pass2dEvent( int jets_cent, float jets_genPt, float jets_recoPt, float jets_eta, float jets_phi, bool isMC)  {

  if ( !passReco2dEvent(jets_cent, jets_recoPt, jets_eta, jets_phi, isMC))
    return false;
  
  return true;
}

bool passEvent( int jets_cent, float jets_genPt, float jets_recoPt, float jets_eta, float jets_phi, bool isMC)  {

  if ( !passRecoEvent(jets_cent, jets_recoPt, jets_eta, jets_phi, isMC))
    return false;
  if ( (isMC) && !passGenEvent(jets_cent, jets_genPt) )
    return false;
  
  return true;
  
} 

bool passEtaCut( float jets_recoEta, int etaBin) {
  int absEta = fabs(jets_recoEta);
  if (  (absEta < 0.3) && (etaBin == 0) )
    return true;
  if (  (absEta >= 0.3) &&  (absEta < 1.2) && (etaBin == 1) )
    return true;
  if (  (absEta >= 1.2) &&  (absEta < 2.1) && (etaBin == 2) )
    return true;
  if (  (absEta < 2.1) && (etaBin == 3) )
    return true;
  
  return false;

}

bool passJesEvent( int jets_cent, float jets_genPt, float jets_recoPt, float jets_eta, float jets_phi)  {

  if ( !passGenEvent(jets_cent, jets_genPt) )
    return false;
  
  double ptCut = 100;
  double ptCutUp = 1000;
  /*
  if ( jets_cent != icent )
    return false;
  */
  if ( jets_recoPt < ptCut )
    return false;

  if ( jets_recoPt > ptCutUp )
    return false;

    // cut for detector hole. Yields of final result should be scaled by ratio of total acceptance we've eliminated
  if ((jets_eta > 0 && jets_eta < 1) && (jets_phi > TMath::Pi()/4 && jets_phi < (11*TMath::Pi()/32))) {
    //    if (!isMC) cout << " returning false detector hole" << endl; 
    return false;
  }
  
  return true;

}

#endif
