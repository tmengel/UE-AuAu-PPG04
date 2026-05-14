#include <iostream>
#include <vector>
#include <utility> 

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TStyle.h>

#include <sPhenixStyle.C>

const std::string sPHENIX_Tag = "#it{#bf{sPHENIX}} Internal";
const std::string DataType_Tag = "Au+Au 200 GeV";
bool OVERRIDE = true;

const float V2_VALUES[] = {2.32, 3.39, 4.76, 6.18, 7.03, 7.4, 7.44, 7.23, 6.96};
const float V3_VALUES[] = {1.45, 1.62, 1.76, 1.9, 1.99, 2.05, 1.92, 1.75, 1.57};
const float X_CENT_BINS[]= {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
const int N_X_CENT_BINS = sizeof(X_CENT_BINS)/sizeof(X_CENT_BINS[0]) - 1;
float MAX_X_CENT = 80;

const float AREA_CONE = TMath::Pi()*0.4*0.4;
const float AREA_TOWER_CEMC = (2.0*TMath::Pi()/256.0)*(2.2/96.0);
const float AREA_HCAL_TOWER = (2.0*TMath::Pi()/64.0)*(2.2/24.0);
const float N_CEMC_TOWERS = 256*94;
const float N_HCALIN_TOWERS = 24*64;
const float N_HCALOUT_TOWERS = 24*64;
int NEVENTS = -1;

const int COLORS[] = {kBlack, kRed, kBlue, kCyan, kGreen+2, kViolet, kOrange+2,kAzure-1,  kMagenta+2, kRed+2};
const int MARKERS[] = { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kOpenCircle, kOpenSquare, kOpenTriangleUp};
const float MARKER_SIZE = 1.2;
const float LINE_WIDTH = 2.0;

std::string output_dir = "plots/conecrosscheck";


