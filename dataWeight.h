#ifndef DATAWEIGHT_H
#define DATAWEIGHT_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TROOT.h>

// General 2015 PbPb recommendations: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/HeavyIonAnalysis2015
// Updated centrality TAA: https://twiki.cern.ch/twiki/pub/AtlasProtected/HeavyIonAnalysis2015/centrality_values_UPDATE_30Mar18.txt
// 0-10%: 1/mb
const double TAA_0 = 233.523e-5;
// 10-20%: 1/mb
const double TAA_1 = 143.345e-5;
// 20-30%: 1/mb
const double TAA_2 = 86.3844e-5;
// 30-40%: 1/mb
const double TAA_3 = 49.4611e-5;
// 40-50%: 1/mb
const double TAA_4 = 26.3435e-5;
// 50-60%: 1/mb
const double TAA_5 = 12.811e-5;

// 60-80%: 1/mb
// in 1/mb: 3.9e-1 or .39e0 
// in 1/b: 390
// in 1/nb: .39e-6 or 3.9e-5
//const double TAA_6 = 0.393994e3;
const double TAA_6 = 3.9e-5; 


const double nEvents_0 = 3.74e8;
const double nEvents_1 = 3.74e8;
const double nEvents_2 = 3.74e8;
const double nEvents_3 = 3.74e8;
const double nEvents_4 = 3.74e8;
const double nEvents_5 = 3.74e8;
const double nEvents_6 = 7.48e8;


// in 1/pb: 2.5e1
// in 1/b: 25e-12
//const double lumi_pp = 25e-12;
//const double lumi_pp = 25000e-9;
// in 1/nb
//const double lumi_pp = 25000;
const double lumi_pp = 2.5e4;


// order e-5 for most central
double hi9EvtWgt_0 = 1. / (TAA_0 * nEvents_0) ;
double hi9EvtWgt_1 = 1. / (TAA_1 * nEvents_1) ;
double hi9EvtWgt_2 = 1. / (TAA_2 * nEvents_2) ;
double hi9EvtWgt_3 = 1. / (TAA_3 * nEvents_3) ;
double hi9EvtWgt_4 = 1. / (TAA_4 * nEvents_4) ;
double hi9EvtWgt_5 = 1. / (TAA_5 * nEvents_5) ;
double hi9EvtWgt_6 = 1. / (TAA_6 * nEvents_6) ;
// order e-4
double ppEvtWgt = 1. / lumi_pp ;

#endif 
