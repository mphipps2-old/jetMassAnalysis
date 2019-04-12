#ifndef JSSUTIL_Yongsun_H
#define JSSUTIL_Yongsun_H
#include "commonUtility.h"

int kPP   = 0 ;
int kPbPb = 1 ;


TH1D* getHistoRweighted( int kSample = kPP, int ipt=0, int icent=0, int nIter=0, int optY=2) {
  TFile* f = new TFile(Form("spectraFiles_fullForResponse_fullForReweighting/unfoldingResult_coll%d_optX1_optY%d_radius0.4_nIter%d.root",kSample,optY, nIter));
  TH1D* ret = (TH1D*)f->Get(Form("hmassFinal_ix%d_icent%d",ipt,icent));
  return ret;
}
TH1D* getHistoMassDepJES( int kSample = kPP, int ipt=0, int icent=0, int nIter=0, int optY=2) {
  TFile* f = new TFile(Form("spectraFiles_massDependentJESDone/unfoldingResult_coll%d_optX1_optY%d_radius0.4_nIter%d.root",kSample,optY, nIter));
  TH1D* ret = (TH1D*)f->Get(Form("hmassFinal_ix%d_icent%d",ipt,icent));
  return ret;
}

TH1D* getHisto( int kSample = kPP, int ipt=0, int icent=0, int nIter=0, int optY=2) {
  TFile* f = new TFile(Form("spectraFiles/unfoldingResult_coll%d_optX1_optY%d_radius0.4_nIter%d.root",kSample,optY, nIter));
  TH1D* ret = (TH1D*)f->Get(Form("hmassFinal_ix%d_icent%d",ipt,icent));
  return ret;
}

TH1D* getMCRaw(  int kSample = kPP, int ipt=0, int icent=0, int nIter=0, int optY=2) {
  TFile* f = new TFile(Form("spectraFiles/unfoldingResult_coll%d_optX1_optY%d_radius0.4_nIter%d.root",kSample,optY,nIter));
  TH1D* ret = (TH1D*)f->Get(Form("hmassRawMC_ix%d_icent%d",ipt,icent));
  return ret;
}
TH1D* getDataRaw( int kSample = kPP, int ipt=0, int icent=0, int nIter=0, int optY=2) {
  TFile* f = new TFile(Form("spectraFiles/unfoldingResult_coll%d_optX1_optY%d_radius0.4_nIter%d.root",kSample,optY,    nIter));
  TH1D* ret = (TH1D*)f->Get(Form("hmassRawData_ix%d_icent%d",ipt,icent));
  return ret;
}


TH1D* getFinal_test4( int kSample = kPP, int ipt=0, int icent=0, int nIter=0, int optY=2) {
  TFile* f = new TFile(Form("spectraFiles/unfoldingResult_coll%d_optX1_optY%d_radius0.4_nIter%d_test4.root",kSample,optY, nIter));
  TH1D* ret = (TH1D*)f->Get(Form("hmassFinal_ix%d_icent%d",ipt,icent));
  return ret;
}

TH1D* getInitial_test4( int kSample = kPP, int ipt=0, int icent=0, int nIter=0, int optY=2) {
  TFile* f = new TFile(Form("spectraFiles/unfoldingResult_coll%d_optX1_optY%d_radius0.4_nIter%d_test4.root",kSample,optY, nIter));
  TH1D* ret = (TH1D*)f->Get(Form("hmassInitial_ix%d_icent%d",ipt,icent));
  return ret;
}

TH1D* getMCRaw_test4(  int kSample = kPP, int ipt=0, int icent=0, int nIter=0, int optY=2) {
  TFile* f = new TFile(Form("spectraFiles/unfoldingResult_coll%d_optX1_optY%d_radius0.4_nIter%d_test4.root",kSample,optY,nIter));
  TH1D* ret = (TH1D*)f->Get(Form("hmassRawMC_ix%d_icent%d",ipt,icent));
  return ret;
}
TH1D* getDataRaw_test4( int kSample = kPP, int ipt=0, int icent=0, int nIter=0, int optY=2) {
  TFile* f = new TFile(Form("spectraFiles/unfoldingResult_coll%d_optX1_optY%d_radius0.4_nIter%d_test4.root",kSample,optY,    nIter));
  TH1D* ret = (TH1D*)f->Get(Form("hmassRawData_ix%d_icent%d",ipt,icent));
  return ret;
}





void drawCentrality( int kSample = kPP, int icent = 0, float xp=0.2, float yp=0.8, int textColor=kBlack, int textSize=18) {
  if ( kSample == kPP)  drawText("#it{pp}",xp,yp,textColor,textSize) ;
  else {
    if ( icent==-1 )  drawText( "Pb+Pb",xp,yp,textColor,textSize) ;
    if ( icent==0 )  drawText( "Pb+Pb 0-10%",xp,yp,textColor,textSize) ;
    if ( icent==1 )  drawText( "Pb+Pb 10-20%",xp,yp,textColor,textSize) ;
    if ( icent==2 )  drawText( "Pb+Pb 20-30%",xp,yp,textColor,textSize) ;
    if ( icent==3 )  drawText( "Pb+Pb 30-40%",xp,yp,textColor,textSize) ;
    if ( icent==4 )  drawText( "Pb+Pb 40-50%",xp,yp,textColor,textSize) ;
    if ( icent==5 )  drawText( "Pb+Pb 50-60%",xp,yp,textColor,textSize) ;
    if ( icent==6 )  drawText( "Pb+Pb 60-80%",xp,yp,textColor,textSize) ;

  }
}


void drawCentralityRCP(int icent = 0, float xp=0.2, float yp=0.8, int textColor=kBlack, int textSize=18){
  if ( icent==0 )  drawText( "0-10% / 60-80%",xp,yp,textColor,textSize) ;
  if ( icent==1 )  drawText( "10-20% / 60-80%",xp,yp,textColor,textSize) ;
  if ( icent==2 )  drawText( "20-30% / 60-80%",xp,yp,textColor,textSize) ;
  if ( icent==3 )  drawText( "30-40% / 60-80%",xp,yp,textColor,textSize) ;
  if ( icent==4 )  drawText( "40-50% / 60-80%",xp,yp,textColor,textSize) ;
  if ( icent==5 )  drawText( "50-60% / 60-80%",xp,yp,textColor,textSize) ;
}


void drawBin(double *xBin, int ix, TString demen = "GeV", float xp=0.2, float yp=0.8, int textColor=kBlack, int textSize=18){
  if ( xBin[ix-1] >= 1  ) drawText( Form("%d - %d %s", (int)(xBin[ix-1]),  (int)(xBin[ix]), demen.Data()) , xp,yp,textColor,textSize) ;
  else if ( xBin[ix-1] >= 0.1 ) drawText( Form("%.2f - %.2f %s", (float)(xBin[ix-1]),  (float)(xBin[ix]), demen.Data()) , xp,yp,textColor,textSize) ;
  else if ( xBin[ix-1] >= 0.01 ) drawText( Form("%.3f - %.3f %s", (float)(xBin[ix-1]),  (float)(xBin[ix]), demen.Data()) , xp,yp,textColor,textSize) ;
  else if ( xBin[ix-1] >= 0.001 ) drawText( Form("%.4f - %.4f %s", (float)(xBin[ix-1]),  (float)(xBin[ix]), demen.Data()) , xp,yp,textColor,textSize) ;
  else if ( xBin[ix-1] >= 0.0001 ) drawText( Form("%.5f - %.5f %s", (float)(xBin[ix-1]),  (float)(xBin[ix]), demen.Data()) , xp,yp,textColor,textSize) ;
  else  drawText( Form("%.5f - %.5f %s", (float)(xBin[ix-1]),  (float)(xBin[ix]), demen.Data()) , xp,yp,textColor,textSize) ;
}
void drawBinPt(double *xBin, int ix, TString demen = "GeV", float xp=0.2, float yp=0.8, int textColor=kBlack, int textSize=18){
  int lowPtBin =  (int)(xBin[ix-1]); 
  if (lowPtBin==125) lowPtBin = 126;
  drawText( Form("%d < p_{T} < %d %s", lowPtBin,  (int)(xBin[ix]), demen.Data()) , xp,yp,textColor,textSize) ;
}
void drawBinMpt(double *xBin, int ix, TString demen = "", float xp=0.2, float yp=0.8, int textColor=kBlack, int textSize=18){
  if ( (xBin[ix-1]) == 0 )   
    drawText( Form("       m/p_{T} < %.3f %s",  (float)(xBin[ix]), demen.Data()) , xp,yp,textColor,textSize) ;
  else 
    drawText( Form("%.2f < m/p_{T} < %.3f %s", (float)(xBin[ix-1]),  (float)(xBin[ix]), demen.Data()) , xp,yp,textColor,textSize) ;
}

void drawBinMpt2(double *xBin, int ix, TString demen = "", float xp=0.2, float yp=0.8, int textColor=kBlack, int textSize=18){
  if ( (xBin[ix-1]) == 0 )   
    drawText( Form("       m/p_{T} < %.2f %s",  (float)(sqrt(xBin[ix])), demen.Data()) , xp,yp,textColor,textSize) ;
  else 
    drawText( Form("%.2f < m/p_{T} < %.2f %s", (float)(sqrt(xBin[ix-1])),  (float)(sqrt(xBin[ix])), demen.Data()) , xp,yp,textColor,textSize) ;
}

void drawBinPt2(double *xBin, int ix, TString demen = "GeV", float xp=0.2, float yp=0.8, int textColor=kBlack, int textSize=18){
  int lowPtBin =  (int)(xBin[ix-1]);
  if (lowPtBin==125) lowPtBin =126;
  drawText( Form("%d < p_{T} < %d %s", lowPtBin,  (int)(xBin[ix+1]), demen.Data()) , xp,yp,textColor,textSize) ;
}

TString textBin(double *xBin, int ix, TString demen = "GeV") {
  if ( xBin[ix-1] >= 1  ) return Form("%d - %d %s", (int)(xBin[ix-1]),  (int)(xBin[ix]), demen.Data()) ;
  else if ( xBin[ix-1] >= 0.1 ) return Form("%.2f - %.2f %s", (float)(xBin[ix-1]),  (float)(xBin[ix]), demen.Data());
  else if ( xBin[ix-1] >= 0.01 ) return Form("%.3f - %.3f %s", (float)(xBin[ix-1]),  (float)(xBin[ix]), demen.Data());
  else if ( xBin[ix-1] >= 0.001 ) return  Form("%.4f - %.4f %s", (float)(xBin[ix-1]),  (float)(xBin[ix]), demen.Data());
  else if ( xBin[ix-1] >= 0.0001 ) return  Form("%.5f - %.5f %s", (float)(xBin[ix-1]),  (float)(xBin[ix]), demen.Data());
  else  return Form("%.5f - %.5f %s", (float)(xBin[ix-1]),  (float)(xBin[ix]), demen.Data());
}

TString textBin2(double *xBin, int ix, TString demen = "GeV") {
  int lowPtBin =  (int)(xBin[ix-1]);
  if (lowPtBin==125) lowPtBin =126;
  return Form("%d < p_{T} < %d %s", (int)(xBin[ix-1]),  (int)(xBin[ix]), demen.Data()) ;

}


#endif
