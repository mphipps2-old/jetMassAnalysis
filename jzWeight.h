#ifndef JZWEIGHT_H
#define JZWEIGHT_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TROOT.h>

// 19.648 mb in units of nb to correspond to data
const double JZ2_cs     = 1.9648E+07;
//576.130 ub
const double JZ3_cs     = 5.7613E+05;
// 41.522 ub
const double JZ4_cs     = 4.1522E+04;

// these filter efficiencies came from Aaron's email ("Reweighting instructions") and the TWiki (https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HIJetMCSamples#Powheg_Pythia8_dijet_samples) and appear to be right (results in a smoother Data/MC spectrum)
const double JZ2_fltEff = 1.2948E-04;
const double JZ3_fltEff = 4.2129E-05;
const double JZ4_fltEff = 2.8563E-06;

// these filter efficiencies came from AMI and appears to be wrong
//const double JZ2_fltEff = 1.4065E-04 ; 
//const double JZ3_fltEff = 4.2577E-05;
//const double JZ4_fltEff = 2.8815E-06;

double sumJetWeight_jz2 = 4.656080e+06; 
double sumJetWeight_jz3 = 5.805497e+04; 
double sumJetWeight_jz4 = 7.924295e+02; 

const double nEvents_jz2  = 5856854;
const double nEvents_jz3  = 5940974;
const double nEvents_jz4  = 5436963;

double hi9EvtWgtJZ2 = (JZ2_cs * JZ2_fltEff ) / (sumJetWeight_jz2 * nEvents_jz2) ;
double hi9EvtWgtJZ3 = (JZ3_cs * JZ3_fltEff ) / (sumJetWeight_jz3 * nEvents_jz3) ;
double hi9EvtWgtJZ4 = (JZ4_cs * JZ4_fltEff ) / (sumJetWeight_jz4 * nEvents_jz4) ;

// Event weight defined as: event weight = (JZ_xsection * JZ_filterEff * Powheg_jetWeight) / (num_events*sum_jet_weights)
// note: jet weight is applied per jet. Corrections are stored in the ntuple and included in the reweighting section of the analysis scripts
//


// dN/dpT (Nb and then divide by bin widths to get nb/GeV)


double perEvtWgtJZ2 = (JZ2_cs * JZ2_fltEff ) / sumJetWeight_jz2 ; 
double perEvtWgtJZ3 = (JZ3_cs * JZ3_fltEff ) / sumJetWeight_jz3 ; 
double perEvtWgtJZ4 = (JZ4_cs * JZ4_fltEff ) / sumJetWeight_jz4 ; 

#endif 
