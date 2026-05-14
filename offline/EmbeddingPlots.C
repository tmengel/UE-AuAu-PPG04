#include <iostream>
#include <vector>
#include <utility> 

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h> 

#include <sPhenixStyle.C>

const std::string sPHENIX_Tag = "#it{#bf{sPHENIX}} Internal";

bool DO_OVERRIDE = true;

const int N_SUM_ET_BINS = 10;
float SUM_ET_BINS[N_SUM_ET_BINS+1];
float MAX_SUM_ET = 1200;

const int N_SUM_Q_BINS = 50;
float SUM_Q_BINS[N_SUM_Q_BINS+1];
float MAX_SUM_Q = 2200;

const float V2_VALUES[] = {2.32, 3.39, 4.76, 6.18, 7.03, 7.4, 7.44, 7.23, 6.96};
const float V3_VALUES[] = {1.45, 1.62, 1.76, 1.9, 1.99, 2.05, 1.92, 1.75, 1.57};
const float X_CENT_BINS[]= {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
const int N_X_CENT_BINS = sizeof(X_CENT_BINS)/sizeof(X_CENT_BINS[0]) - 1;
float MAX_X_CENT = 80;

const int N_CONECOMP_BINS = 210;
float CONECOMP_MAX = 210;
float CONECOMP_BINS[N_CONECOMP_BINS+1];

const int N_CONECOMP_SUB1_BINS = 120;
float CONECOMP_SUB1_MAX = 120;
float CONECOMP_SUB1_BINS[N_CONECOMP_SUB1_BINS+1];

const int N_CONE_DET_BINS = 250;
float MAX_CONE_DET = 60;
float CONE_DET_BINS[N_CONE_DET_BINS+1];

const int N_ET_BINS = 15;
float MAX_ET = 150;
float ET_BINS[N_ET_BINS+1];

// const int N_PROBE_DET_BINS = 120;
// float MAX_PROBE_DET = 59.5;
// float MIN_PROBE_DET = -60.5; // to include negative values for backward compatibility
// float PROBE_DET_BINS[N_PROBE_DET_BINS+1];


const int N_PROBE_DET_BINS = 250;
float MAX_PROBE_DET = 100;
float PROBE_DET_BINS[N_PROBE_DET_BINS+1];


const int N_RHO_BINS = 200;
float RHO_M_BINS[N_RHO_BINS+1];
float RHO_A_BINS[N_RHO_BINS+1];
float TOWER_BACKGROUND_BINS[N_RHO_BINS+1];
float BACKGROUND_BINS[N_RHO_BINS+1];
float MAX_RHO_M = 0.08;
float MIN_RHO_M = 0;
float MAX_RHO_A = 120;
float MIN_RHO_A = 0;



void SetBins(){
    for ( int i = 0; i < N_SUM_ET_BINS+1; ++i ) { SUM_ET_BINS[i] = i*MAX_SUM_ET/N_SUM_ET_BINS; }
    for ( int i = 0; i < N_SUM_Q_BINS+1; ++i ) { SUM_Q_BINS[i] = i*MAX_SUM_Q/N_SUM_Q_BINS; }
    for ( int i = 0; i < N_CONECOMP_BINS+1; ++i ) { CONECOMP_BINS[i] =  i*CONECOMP_MAX/N_CONECOMP_BINS; }
    for ( int i = 0; i < N_CONECOMP_SUB1_BINS+1; ++i ) { CONECOMP_SUB1_BINS[i] = i*CONECOMP_SUB1_MAX/N_CONECOMP_SUB1_BINS; }
    for ( int i = 0; i < N_CONE_DET_BINS+1; ++i ) { CONE_DET_BINS[i] = -MAX_CONE_DET + i*2*MAX_CONE_DET/N_CONE_DET_BINS; }
    // for ( int i = 0; i < N_PROBE_DET_BINS+1; ++i ) { PROBE_DET_BINS[i] = MIN_PROBE_DET + i*(MAX_PROBE_DET-MIN_PROBE_DET)/N_PROBE_DET_BINS; }
    for ( int i = 0; i < N_PROBE_DET_BINS+1; ++i ) { PROBE_DET_BINS[i] = -100 + i*2*MAX_PROBE_DET/N_PROBE_DET_BINS; }
    for ( int i = 0; i < N_ET_BINS+1; ++i ) { ET_BINS[i] =  -1 + i*(MAX_ET)/N_ET_BINS; }
    for ( int i = 0; i < N_RHO_BINS+1; ++i ) { 
        RHO_M_BINS[i] = MIN_RHO_M + i*(MAX_RHO_M-MIN_RHO_M)/N_RHO_BINS; 
        RHO_A_BINS[i] = MIN_RHO_A + i*(MAX_RHO_A-MIN_RHO_A)/N_RHO_BINS; 
    }
}

const float AREA_CONE = TMath::Pi()*0.4*0.4;
const float AREA_TOWER_CEMC = (2.0*TMath::Pi()/256.0)*(2.2/96.0);
const float AREA_HCAL_TOWER = (2.0*TMath::Pi()/64.0)*(2.2/24.0);

bool IS_DATA = false;
int NEVENTS = 0;
std::string DataType_Tag;


const int COLORS[] = {kBlack, kRed , kAzure-2, kGreen+2, kViolet, kCyan, kOrange+2, kMagenta+2, kAzure-2};
const int MARKERS[] = { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kOpenCircle, kOpenSquare, kOpenTriangleUp};
const float MARKER_SIZE = 1.2;
const float LINE_WIDTH = 2.0;

std::string probe_plots;

void ConfigureOutputDirs(std::string input_file_base, std::string plotting_dir = "plots/");
std::string MakeGetDir( const std::string & dir ){
    if ( gSystem->AccessPathName(dir.c_str()) ) {
        gSystem->Exec(Form("mkdir -p %s", dir.c_str()));
    }
    return dir;
}

std::string ProcessEmbedTree(const std::string & input_file_10, const std::string & input_file_30, const std::string & prefix);

void EmbeddingPlots() 
{

    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    SetsPhenixStyle();
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kRainBow);

   
    DataType_Tag = "Au+Au 200 GeV";
    const std::string & input_file_30 = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/FINAL_EMBED/jet30.root";
    const std::string & input_file_10 = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/FINAL_EMBED/jet10.root";

    // const std::string & input_file_30 = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/EMBEDDING_OMIT2/jet30.root";
    // const std::string & input_file_10 = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/EMBEDDING_OMIT2/jet10.root";
    // get base name of input file
    std::string input_file_base = input_file_30;
    size_t found = input_file_30.find_last_of("/");
    if (found != std::string::npos) {
        input_file_base = input_file_30.substr(found+1);
    }
    found = input_file_base.find_last_of(".");
    if (found != std::string::npos) {
        input_file_base = input_file_base.substr(0, found);
    }

  
    
    ConfigureOutputDirs("EmbedJets_Test", "prelim/");
    std::cout << "Embed plots: " << probe_plots << std::endl;

    SetBins();
    
    std::string probe_hist = ProcessEmbedTree(input_file_10, input_file_30, "embed");
    
    gSystem->Exit(0);
  
   
}

std::string ProcessEmbedTree(const std::string & input_file_10, const std::string & input_file_30, const std::string & prefix  )
{

    std::string outdir = probe_plots;
    outdir += "/" + prefix;
    if( !gSystem->OpenDirectory(outdir.c_str()) ) {
        gSystem->mkdir(outdir.c_str(), true);
    }

    TFile * f_10 = new TFile(input_file_10.c_str(), "READ");
    if( !f_10->IsOpen() || f_10->IsZombie() ) { std::cout << "File " << input_file_10 << " is zombie" << std::endl;  exit(1); }
    TTree * t_10 = (TTree*)f_10->Get("T");
    
    // tree branches 
      int centrality_10 = 0;
      t_10->SetBranchAddress("centrality", &centrality_10);
  
      float rho_val_TowerRho_AREA_10 = 0;
      float std_rho_val_TowerRho_AREA_10 = 0;
      t_10->SetBranchAddress("rho_val_TowerRho_AREA", &rho_val_TowerRho_AREA_10);
      t_10->SetBranchAddress("std_rho_val_TowerRho_AREA", &std_rho_val_TowerRho_AREA_10);
      
      float rho_val_TowerRho_MULT_10 = 0;
      float std_rho_val_TowerRho_MULT_10 = 0;
      t_10->SetBranchAddress("rho_val_TowerRho_MULT", &rho_val_TowerRho_MULT_10);
      t_10->SetBranchAddress("std_rho_val_TowerRho_MULT", &std_rho_val_TowerRho_MULT_10);
  
      float rho_val_TowerRho_AREA_CEMC_10 = 0;
      float std_rho_val_TowerRho_AREA_CEMC_10 = 0;
      t_10->SetBranchAddress("rho_val_TowerRho_AREA_CEMC", &rho_val_TowerRho_AREA_CEMC_10);
      t_10->SetBranchAddress("std_rho_val_TowerRho_AREA_CEMC", &std_rho_val_TowerRho_AREA_CEMC_10);
  
      float rho_val_TowerRho_MULT_CEMC_10 = 0;
      float std_rho_val_TowerRho_MULT_CEMC_10 = 0;
      t_10->SetBranchAddress("rho_val_TowerRho_MULT_CEMC", &rho_val_TowerRho_MULT_CEMC_10);
      t_10->SetBranchAddress("std_rho_val_TowerRho_MULT_CEMC", &std_rho_val_TowerRho_MULT_CEMC_10);
  
      float rho_val_TowerRho_AREA_HCALIN_10 = 0;
      float std_rho_val_TowerRho_AREA_HCALIN_10 = 0;
      t_10->SetBranchAddress("rho_val_TowerRho_AREA_HCALIN", &rho_val_TowerRho_AREA_HCALIN_10);
      t_10->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALIN", &std_rho_val_TowerRho_AREA_HCALIN_10);
  
      float rho_val_TowerRho_MULT_HCALIN_10 = 0;
      float std_rho_val_TowerRho_MULT_HCALIN_10 = 0;
      t_10->SetBranchAddress("rho_val_TowerRho_MULT_HCALIN", &rho_val_TowerRho_MULT_HCALIN_10);
      t_10->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALIN", &std_rho_val_TowerRho_MULT_HCALIN_10);
  
      float rho_val_TowerRho_AREA_HCALOUT_10 = 0;
      float std_rho_val_TowerRho_AREA_HCALOUT_10 = 0;
      t_10->SetBranchAddress("rho_val_TowerRho_AREA_HCALOUT", &rho_val_TowerRho_AREA_HCALOUT_10);
      t_10->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALOUT", &std_rho_val_TowerRho_AREA_HCALOUT_10);
  
      float rho_val_TowerRho_MULT_HCALOUT_10 = 0;
      float std_rho_val_TowerRho_MULT_HCALOUT_10 = 0;
      t_10->SetBranchAddress("rho_val_TowerRho_MULT_HCALOUT", &rho_val_TowerRho_MULT_HCALOUT_10);
      t_10->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALOUT", &std_rho_val_TowerRho_MULT_HCALOUT_10);
  
      std::vector< float > * emb_jet_eta_10 = 0;
      std::vector< float > * emb_jet_phi_10 = 0;
      std::vector< float > * emb_jet_energy_10 = 0;
      std::vector< float > * emb_jet_area_10 = 0;
      std::vector< float > * emb_jet_energy_cemc_10 = 0;
      std::vector< float > * emb_jet_energy_hcalin_10 = 0;
      std::vector< float > * emb_jet_energy_hcalout_10 = 0;
      std::vector< int > * emb_jet_num_towers_10 = 0;
      std::vector< int > * emb_jet_num_towers_cemc_10 = 0;
      std::vector< int > * emb_jet_num_towers_hcalin_10 = 0;
      std::vector< int > * emb_jet_num_towers_hcalout_10 = 0;
      
      std::vector< float > * emb_jet_sub1_eta_10 = 0;
      std::vector< float > * emb_jet_sub1_phi_10 = 0;
      std::vector< float > * emb_jet_sub1_energy_10 = 0;
      std::vector< float > * emb_jet_sub1_area_10 = 0;
      std::vector< float > * emb_jet_sub1_energy_cemc_10 = 0;
      std::vector< float > * emb_jet_sub1_energy_hcalin_10 = 0;
      std::vector< float > * emb_jet_sub1_energy_hcalout_10 = 0;
      std::vector< int > * emb_jet_sub1_num_towers_10 = 0;
      std::vector< int > * emb_jet_sub1_num_towers_cemc_10 = 0;
      std::vector< int > * emb_jet_sub1_num_towers_hcalin_10 = 0;
      std::vector< int > * emb_jet_sub1_num_towers_hcalout_10 = 0;
      
      std::vector< float > * sim_jet_eta_10 = 0;
      std::vector< float > * sim_jet_phi_10 = 0;
      std::vector< float > * sim_jet_energy_10 = 0;
      std::vector< int > * sim_jet_num_towers_10 = 0;
      std::vector< int > * sim_jet_num_towers_cemc_10 = 0;
      std::vector< int > * sim_jet_num_towers_hcalin_10 = 0;
      std::vector< int > * sim_jet_num_towers_hcalout_10 = 0;
  
      std::vector< float > * sim_retower_jet_eta_10 = 0;
      std::vector< float > * sim_retower_jet_phi_10 = 0;
      std::vector< float > * sim_retower_jet_energy_10 = 0;
      std::vector< int > * sim_retower_jet_num_towers_10 = 0;
      std::vector< int > * sim_retower_jet_num_towers_cemc_10 = 0;
      std::vector< int > * sim_retower_jet_num_towers_hcalin_10 = 0;
      std::vector< int > * sim_retower_jet_num_towers_hcalout_10 = 0;
  
      std::vector< float > * truth_jet_eta_10 = 0;
      std::vector< float > * truth_jet_phi_10 = 0;
      std::vector< float > * truth_jet_energy_10 = 0;
      std::vector< int > * truth_jet_ncomp_10 = 0;
  
      t_10->SetBranchAddress("emb_jet_eta", &emb_jet_eta_10);
      t_10->SetBranchAddress("emb_jet_phi", &emb_jet_phi_10);
      t_10->SetBranchAddress("emb_jet_energy", &emb_jet_energy_10);
      t_10->SetBranchAddress("emb_jet_area", &emb_jet_area_10);
      t_10->SetBranchAddress("emb_jet_energy_cemc", &emb_jet_energy_cemc_10);
      t_10->SetBranchAddress("emb_jet_energy_hcalin", &emb_jet_energy_hcalin_10);
      t_10->SetBranchAddress("emb_jet_energy_hcalout", &emb_jet_energy_hcalout_10);
      t_10->SetBranchAddress("emb_jet_num_towers", &emb_jet_num_towers_10);
      t_10->SetBranchAddress("emb_jet_num_towers_cemc", &emb_jet_num_towers_cemc_10);
      t_10->SetBranchAddress("emb_jet_num_towers_hcalin", &emb_jet_num_towers_hcalin_10);
      t_10->SetBranchAddress("emb_jet_num_towers_hcalout", &emb_jet_num_towers_hcalout_10);
  
      t_10->SetBranchAddress("emb_jet_sub1_eta", &emb_jet_sub1_eta_10);
      t_10->SetBranchAddress("emb_jet_sub1_phi", &emb_jet_sub1_phi_10);
      t_10->SetBranchAddress("emb_jet_sub1_energy", &emb_jet_sub1_energy_10);
      t_10->SetBranchAddress("emb_jet_sub1_area", &emb_jet_sub1_area_10);
      t_10->SetBranchAddress("emb_jet_sub1_energy_cemc", &emb_jet_sub1_energy_cemc_10);
      t_10->SetBranchAddress("emb_jet_sub1_energy_hcalin", &emb_jet_sub1_energy_hcalin_10);
      t_10->SetBranchAddress("emb_jet_sub1_energy_hcalout", &emb_jet_sub1_energy_hcalout_10);
      t_10->SetBranchAddress("emb_jet_sub1_num_towers", &emb_jet_sub1_num_towers_10);
      t_10->SetBranchAddress("emb_jet_sub1_num_towers_cemc", &emb_jet_sub1_num_towers_cemc_10);
      t_10->SetBranchAddress("emb_jet_sub1_num_towers_hcalin", &emb_jet_sub1_num_towers_hcalin_10);
      t_10->SetBranchAddress("emb_jet_sub1_num_towers_hcalout", &emb_jet_sub1_num_towers_hcalout_10); 
      
      t_10->SetBranchAddress("sim_jet_eta", &sim_jet_eta_10);
      t_10->SetBranchAddress("sim_jet_phi", &sim_jet_phi_10);
      t_10->SetBranchAddress("sim_jet_energy", &sim_jet_energy_10);
      t_10->SetBranchAddress("sim_jet_num_towers", &sim_jet_num_towers_10);
      t_10->SetBranchAddress("sim_jet_num_towers_cemc", &sim_jet_num_towers_cemc_10);
      t_10->SetBranchAddress("sim_jet_num_towers_hcalin", &sim_jet_num_towers_hcalin_10);
      t_10->SetBranchAddress("sim_jet_num_towers_hcalout", &sim_jet_num_towers_hcalout_10);
  
      t_10->SetBranchAddress("sim_retower_jet_eta", &sim_retower_jet_eta_10);
      t_10->SetBranchAddress("sim_retower_jet_phi", &sim_retower_jet_phi_10);
      t_10->SetBranchAddress("sim_retower_jet_energy", &sim_retower_jet_energy_10);
      t_10->SetBranchAddress("sim_retower_jet_num_towers", &sim_retower_jet_num_towers_10);
      t_10->SetBranchAddress("sim_retower_jet_num_towers_cemc", &sim_retower_jet_num_towers_cemc_10);
      t_10->SetBranchAddress("sim_retower_jet_num_towers_hcalin", &sim_retower_jet_num_towers_hcalin_10);
      t_10->SetBranchAddress("sim_retower_jet_num_towers_hcalout", &sim_retower_jet_num_towers_hcalout_10);
  
      t_10->SetBranchAddress("truth_jet_eta", &truth_jet_eta_10);
      t_10->SetBranchAddress("truth_jet_phi", &truth_jet_phi_10);
      t_10->SetBranchAddress("truth_jet_energy", &truth_jet_energy_10);
      t_10->SetBranchAddress("truth_jet_ncomp", &truth_jet_ncomp_10);
      
    int nentries_10 = t_10->GetEntries();
    std::cout << "Jet 10: " << nentries_10 << " events" << std::endl;
  
    TFile * f_30 = new TFile(input_file_30.c_str(), "READ");
    if( !f_30->IsOpen() || f_30->IsZombie() ) { std::cout << "File " << input_file_30 << " is zombie" << std::endl;  exit(1); }
    TTree * t_30 = (TTree*)f_30->Get("T");

    // tree branches 
        int centrality_30 = 0;
        t_30->SetBranchAddress("centrality", &centrality_30);

        float rho_val_TowerRho_AREA_30 = 0;
        float std_rho_val_TowerRho_AREA_30 = 0;
        t_30->SetBranchAddress("rho_val_TowerRho_AREA", &rho_val_TowerRho_AREA_30);
        t_30->SetBranchAddress("std_rho_val_TowerRho_AREA", &std_rho_val_TowerRho_AREA_30);
        
        float rho_val_TowerRho_MULT_30 = 0;
        float std_rho_val_TowerRho_MULT_30 = 0;
        t_30->SetBranchAddress("rho_val_TowerRho_MULT", &rho_val_TowerRho_MULT_30);
        t_30->SetBranchAddress("std_rho_val_TowerRho_MULT", &std_rho_val_TowerRho_MULT_30);

        float rho_val_TowerRho_AREA_CEMC_30 = 0;
        float std_rho_val_TowerRho_AREA_CEMC_30 = 0;
        t_30->SetBranchAddress("rho_val_TowerRho_AREA_CEMC", &rho_val_TowerRho_AREA_CEMC_30);
        t_30->SetBranchAddress("std_rho_val_TowerRho_AREA_CEMC", &std_rho_val_TowerRho_AREA_CEMC_30);

        float rho_val_TowerRho_MULT_CEMC_30 = 0;
        float std_rho_val_TowerRho_MULT_CEMC_30 = 0;
        t_30->SetBranchAddress("rho_val_TowerRho_MULT_CEMC", &rho_val_TowerRho_MULT_CEMC_30);
        t_30->SetBranchAddress("std_rho_val_TowerRho_MULT_CEMC", &std_rho_val_TowerRho_MULT_CEMC_30);

        float rho_val_TowerRho_AREA_HCALIN_30 = 0;
        float std_rho_val_TowerRho_AREA_HCALIN_30 = 0;
        t_30->SetBranchAddress("rho_val_TowerRho_AREA_HCALIN", &rho_val_TowerRho_AREA_HCALIN_30);
        t_30->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALIN", &std_rho_val_TowerRho_AREA_HCALIN_30);

        float rho_val_TowerRho_MULT_HCALIN_30 = 0;
        float std_rho_val_TowerRho_MULT_HCALIN_30 = 0;
        t_30->SetBranchAddress("rho_val_TowerRho_MULT_HCALIN", &rho_val_TowerRho_MULT_HCALIN_30);
        t_30->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALIN", &std_rho_val_TowerRho_MULT_HCALIN_30);

        float rho_val_TowerRho_AREA_HCALOUT_30 = 0;
        float std_rho_val_TowerRho_AREA_HCALOUT_30 = 0;
        t_30->SetBranchAddress("rho_val_TowerRho_AREA_HCALOUT", &rho_val_TowerRho_AREA_HCALOUT_30);
        t_30->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALOUT", &std_rho_val_TowerRho_AREA_HCALOUT_30);

        float rho_val_TowerRho_MULT_HCALOUT_30 = 0;
        float std_rho_val_TowerRho_MULT_HCALOUT_30 = 0;
        t_30->SetBranchAddress("rho_val_TowerRho_MULT_HCALOUT", &rho_val_TowerRho_MULT_HCALOUT_30);
        t_30->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALOUT", &std_rho_val_TowerRho_MULT_HCALOUT_30);

        std::vector< float > * emb_jet_eta_30 = 0;
        std::vector< float > * emb_jet_phi_30 = 0;
        std::vector< float > * emb_jet_energy_30 = 0;
        std::vector< float > * emb_jet_area_30 = 0;
        std::vector< float > * emb_jet_energy_cemc_30 = 0;
        std::vector< float > * emb_jet_energy_hcalin_30 = 0;
        std::vector< float > * emb_jet_energy_hcalout_30 = 0;
        std::vector< int > * emb_jet_num_towers_30 = 0;
        std::vector< int > * emb_jet_num_towers_cemc_30 = 0;
        std::vector< int > * emb_jet_num_towers_hcalin_30 = 0;
        std::vector< int > * emb_jet_num_towers_hcalout_30 = 0;
        
        std::vector< float > * emb_jet_sub1_eta_30 = 0;
        std::vector< float > * emb_jet_sub1_phi_30 = 0;
        std::vector< float > * emb_jet_sub1_energy_30 = 0;
        std::vector< float > * emb_jet_sub1_area_30 = 0;
        std::vector< float > * emb_jet_sub1_energy_cemc_30 = 0;
        std::vector< float > * emb_jet_sub1_energy_hcalin_30 = 0;
        std::vector< float > * emb_jet_sub1_energy_hcalout_30 = 0;
        std::vector< int > * emb_jet_sub1_num_towers_30 = 0;
        std::vector< int > * emb_jet_sub1_num_towers_cemc_30 = 0;
        std::vector< int > * emb_jet_sub1_num_towers_hcalin_30 = 0;
        std::vector< int > * emb_jet_sub1_num_towers_hcalout_30 = 0;
        
        std::vector< float > * sim_jet_eta_30 = 0;
        std::vector< float > * sim_jet_phi_30 = 0;
        std::vector< float > * sim_jet_energy_30 = 0;
        std::vector< int > * sim_jet_num_towers_30 = 0;
        std::vector< int > * sim_jet_num_towers_cemc_30 = 0;
        std::vector< int > * sim_jet_num_towers_hcalin_30 = 0;
        std::vector< int > * sim_jet_num_towers_hcalout_30 = 0;

        std::vector< float > * sim_retower_jet_eta_30 = 0;
        std::vector< float > * sim_retower_jet_phi_30 = 0;
        std::vector< float > * sim_retower_jet_energy_30 = 0;
        std::vector< int > * sim_retower_jet_num_towers_30 = 0;
        std::vector< int > * sim_retower_jet_num_towers_cemc_30 = 0;
        std::vector< int > * sim_retower_jet_num_towers_hcalin_30 = 0;
        std::vector< int > * sim_retower_jet_num_towers_hcalout_30 = 0;

        std::vector< float > * truth_jet_eta_30 = 0;
        std::vector< float > * truth_jet_phi_30 = 0;
        std::vector< float > * truth_jet_energy_30 = 0;
        std::vector< int > * truth_jet_ncomp_30 = 0;

        t_30->SetBranchAddress("emb_jet_eta", &emb_jet_eta_30);
        t_30->SetBranchAddress("emb_jet_phi", &emb_jet_phi_30);
        t_30->SetBranchAddress("emb_jet_energy", &emb_jet_energy_30);
        t_30->SetBranchAddress("emb_jet_area", &emb_jet_area_30);
        t_30->SetBranchAddress("emb_jet_energy_cemc", &emb_jet_energy_cemc_30);
        t_30->SetBranchAddress("emb_jet_energy_hcalin", &emb_jet_energy_hcalin_30);
        t_30->SetBranchAddress("emb_jet_energy_hcalout", &emb_jet_energy_hcalout_30);
        t_30->SetBranchAddress("emb_jet_num_towers", &emb_jet_num_towers_30);
        t_30->SetBranchAddress("emb_jet_num_towers_cemc", &emb_jet_num_towers_cemc_30);
        t_30->SetBranchAddress("emb_jet_num_towers_hcalin", &emb_jet_num_towers_hcalin_30);
        t_30->SetBranchAddress("emb_jet_num_towers_hcalout", &emb_jet_num_towers_hcalout_30);

        t_30->SetBranchAddress("emb_jet_sub1_eta", &emb_jet_sub1_eta_30);
        t_30->SetBranchAddress("emb_jet_sub1_phi", &emb_jet_sub1_phi_30);
        t_30->SetBranchAddress("emb_jet_sub1_energy", &emb_jet_sub1_energy_30);
        t_30->SetBranchAddress("emb_jet_sub1_area", &emb_jet_sub1_area_30);
        t_30->SetBranchAddress("emb_jet_sub1_energy_cemc", &emb_jet_sub1_energy_cemc_30);
        t_30->SetBranchAddress("emb_jet_sub1_energy_hcalin", &emb_jet_sub1_energy_hcalin_30);
        t_30->SetBranchAddress("emb_jet_sub1_energy_hcalout", &emb_jet_sub1_energy_hcalout_30);
        t_30->SetBranchAddress("emb_jet_sub1_num_towers", &emb_jet_sub1_num_towers_30);
        t_30->SetBranchAddress("emb_jet_sub1_num_towers_cemc", &emb_jet_sub1_num_towers_cemc_30);
        t_30->SetBranchAddress("emb_jet_sub1_num_towers_hcalin", &emb_jet_sub1_num_towers_hcalin_30);
        t_30->SetBranchAddress("emb_jet_sub1_num_towers_hcalout", &emb_jet_sub1_num_towers_hcalout_30); 
        
        t_30->SetBranchAddress("sim_jet_eta", &sim_jet_eta_30);
        t_30->SetBranchAddress("sim_jet_phi", &sim_jet_phi_30);
        t_30->SetBranchAddress("sim_jet_energy", &sim_jet_energy_30);
        t_30->SetBranchAddress("sim_jet_num_towers", &sim_jet_num_towers_30);
        t_30->SetBranchAddress("sim_jet_num_towers_cemc", &sim_jet_num_towers_cemc_30);
        t_30->SetBranchAddress("sim_jet_num_towers_hcalin", &sim_jet_num_towers_hcalin_30);
        t_30->SetBranchAddress("sim_jet_num_towers_hcalout", &sim_jet_num_towers_hcalout_30);

        t_30->SetBranchAddress("sim_retower_jet_eta", &sim_retower_jet_eta_30);
        t_30->SetBranchAddress("sim_retower_jet_phi", &sim_retower_jet_phi_30);
        t_30->SetBranchAddress("sim_retower_jet_energy", &sim_retower_jet_energy_30);
        t_30->SetBranchAddress("sim_retower_jet_num_towers", &sim_retower_jet_num_towers_30);
        t_30->SetBranchAddress("sim_retower_jet_num_towers_cemc", &sim_retower_jet_num_towers_cemc_30);
        t_30->SetBranchAddress("sim_retower_jet_num_towers_hcalin", &sim_retower_jet_num_towers_hcalin_30);
        t_30->SetBranchAddress("sim_retower_jet_num_towers_hcalout", &sim_retower_jet_num_towers_hcalout_30);

        t_30->SetBranchAddress("truth_jet_eta", &truth_jet_eta_30);
        t_30->SetBranchAddress("truth_jet_phi", &truth_jet_phi_30);
        t_30->SetBranchAddress("truth_jet_energy", &truth_jet_energy_30);
        t_30->SetBranchAddress("truth_jet_ncomp", &truth_jet_ncomp_30);




    int nentries_30 = t_30->GetEntries();
    std::cout << "Jet 30:  " << nentries_30 << " events" << std::endl;

    // get the number of events to process
    int nentries = std::min(nentries_10, nentries_30);
    std::cout << "Processing " << nentries << " events" << std::endl;

    TH2F * h2_area_res_vs_x = new TH2F("h2_area_res", "h2_area_res", N_X_CENT_BINS, X_CENT_BINS, N_PROBE_DET_BINS, PROBE_DET_BINS);
    TH2F * h2_mult_res_vs_x = new TH2F("h2_mult_res", "h2_mult_res", N_X_CENT_BINS, X_CENT_BINS, N_PROBE_DET_BINS, PROBE_DET_BINS);
    TH2F * h2_sub1_res_vs_x = new TH2F("h2_sub1_res", "h2_sub1_res", N_X_CENT_BINS, X_CENT_BINS, N_PROBE_DET_BINS, PROBE_DET_BINS);
    TH3F * h2_et_ntruth_cent = new TH3F("h2_et_ntruth_cent", "h2_et_ntruth_cent", N_X_CENT_BINS, X_CENT_BINS, N_ET_BINS, ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
    TH3F * h2_et_ntruth_cemc_cent = new TH3F("h2_et_ntruth_cemc_cent", "h2_et_ntruth_cemc_cent", N_X_CENT_BINS, X_CENT_BINS, N_ET_BINS, ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
    TH3F * h2_et_ntruth_hcalin_cent = new TH3F("h2_et_ntruth_hcalin_cent", "h2_et_ntruth_hcalin_cent", N_X_CENT_BINS, X_CENT_BINS, N_ET_BINS, ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
    TH3F * h2_et_ntruth_hcalout_cent = new TH3F("h2_et_ntruth_hcalout_cent", "h2_et_ntruth_hcalout_cent", N_X_CENT_BINS, X_CENT_BINS, N_ET_BINS, ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
   
    // initialize tprofiles
    const float X_COURSE_CENT_BINS[]= {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
    const int N_X_COURSE_CENT_BINS = sizeof(X_COURSE_CENT_BINS)/sizeof(X_COURSE_CENT_BINS[0]) - 1;
    TProfile * p_et_ntruth[N_X_COURSE_CENT_BINS];
    TProfile * p_et_ntruth_cemc[N_X_COURSE_CENT_BINS];
    TProfile * p_et_ntruth_hcalin[N_X_COURSE_CENT_BINS];
    TProfile * p_et_ntruth_hcalout[N_X_COURSE_CENT_BINS];

    for (  unsigned ibin = 0; ibin < N_X_COURSE_CENT_BINS; ibin++ ){
        p_et_ntruth[ibin] = new TProfile(Form("p_et_ntruth_%d", ibin), Form("p_et_ntruth_%d", ibin), N_ET_BINS, ET_BINS);
        p_et_ntruth_cemc[ibin] = new TProfile(Form("p_et_ntruth_cemc_%d", ibin), Form("p_et_ntruth_cemc_%d", ibin), N_ET_BINS, ET_BINS);
        p_et_ntruth_hcalin[ibin] = new TProfile(Form("p_et_ntruth_hcalin_%d", ibin), Form("p_et_ntruth_hcalin_%d", ibin), N_ET_BINS, ET_BINS);
        p_et_ntruth_hcalout[ibin] = new TProfile(Form("p_et_ntruth_hcalout_%d", ibin), Form("p_et_ntruth_hcalout_%d", ibin), N_ET_BINS, ET_BINS);
        p_et_ntruth[ibin]->GetXaxis()->SetTitle("E_{T}^{Raw} [GeV]");
        p_et_ntruth[ibin]->GetYaxis()->SetTitle("N_{towers}^{Sim.}");
        p_et_ntruth_cemc[ibin]->GetXaxis()->SetTitle("E_{T}^{Raw} [GeV]");
        p_et_ntruth_cemc[ibin]->GetYaxis()->SetTitle("N_{towers}^{Sim.}");
        p_et_ntruth_hcalin[ibin]->GetXaxis()->SetTitle("E_{T}^{Raw} [GeV]");
        p_et_ntruth_hcalin[ibin]->GetYaxis()->SetTitle("N_{towers}^{Sim.}");
        p_et_ntruth_hcalout[ibin]->GetXaxis()->SetTitle("E_{T}^{Raw} [GeV]");
        p_et_ntruth_hcalout[ibin]->GetYaxis()->SetTitle("N_{towers}^{Sim.}");
    }

    TH2F * h2_rhoM_vs_x = new TH2F("h2_rhoM_vs_x", "h2_rhoM_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_BINS, RHO_M_BINS);
    TH2F * h2_rhoA_vs_x = new TH2F("h2_rhoA_vs_x", "h2_rhoA_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_BINS, RHO_A_BINS);

    std::vector<TH3F*> h3s = {h2_et_ntruth_cent, h2_et_ntruth_cemc_cent, h2_et_ntruth_hcalin_cent, h2_et_ntruth_hcalout_cent};
    for ( auto h3 : h3s ) {
        h3->GetXaxis()->SetTitle("Centrality [%]");
        h3->GetYaxis()->SetTitle("E_{T}^{Raw} [GeV]");
        h3->GetZaxis()->SetTitle("N_{truth}");
    }
   
    std::vector<TH2F*> h2s = {h2_area_res_vs_x, h2_mult_res_vs_x, h2_sub1_res_vs_x};
    for ( auto h2 : h2s ) {
        h2->GetXaxis()->SetTitle("Centrality [%]");
        h2->GetYaxis()->SetTitle("#delta E_{T}^{Jet} [GeV]");
    }
    
    std::vector<std::string> labs = {"Area", "Multiplicity", "Iterative"};
   
    int ntruth_total = 0;
    int ntruth_matched = 0;
    int nemb_total = 0;

    float weight_10 = 3.210e-6; //3.210μb
    float weight_30 = 2.178e-9;
    // float weight_10 = 1.0; //3.210μb
    // float weight_30 = 1.0;

    // jet 10
    std::cout << "Looping over jet 10" << std::endl;
    for ( int i = 0; i < nentries; ++i ) {
       
        t_10->GetEntry(i);
       
        float xaxis_var =1.0*centrality_10;
        if ( xaxis_var <0 || xaxis_var > MAX_X_CENT ) { continue; }

        int icourse_bin = -1;
        for ( unsigned ibin = 0; ibin < N_X_COURSE_CENT_BINS; ibin++ ){
            if ( xaxis_var > X_COURSE_CENT_BINS[ibin] && xaxis_var <= X_COURSE_CENT_BINS[ibin+1] ){
                icourse_bin = ibin;
                break;
            }
        }
    
        unsigned int ntruth = sim_jet_eta_10->size();
        unsigned int nemb = emb_jet_eta_10->size();
        ntruth_total += ntruth;
        nemb_total += nemb;
        // match truth to embedded
        for ( unsigned itruth = 0; itruth < ntruth; itruth++) {
    
            float teta = sim_jet_eta_10->at(itruth);
            float tphi = sim_jet_phi_10->at(itruth);
            float tet = sim_jet_energy_10->at(itruth);
            if ( std::fabs(teta)> 0.6) { continue; }
            if ( tet < 5.0 ) { continue; }
            int ncomp = sim_jet_num_towers_10->at(itruth);
            int ncomp_cemc = sim_jet_num_towers_cemc_10->at(itruth);
            int ncomp_hcalin = sim_jet_num_towers_hcalin_10->at(itruth);
            int ncomp_hcalout = sim_jet_num_towers_hcalout_10->at(itruth);
            bool is_matched = false;   
            float dr_current = 1000.0; 
            int iemb_match = -1;
            for ( unsigned iemb = 0; iemb < nemb; iemb++ ) {
                float reta = emb_jet_eta_10->at(iemb);
                float rphi = emb_jet_phi_10->at(iemb);
                float ret  = emb_jet_energy_10->at(iemb);
                if ( ret < 10.0 ) { continue; }
                if ( std::fabs(reta)> 0.6) { continue; }
                float deta = teta - reta;
                float dphi = tphi - rphi;
                if ( dphi > TMath::Pi() ) { dphi -= 2*TMath::Pi(); }
                if ( dphi < -TMath::Pi() ) { dphi += 2*TMath::Pi(); }
                float dr = TMath::Sqrt(deta*deta + dphi*dphi);
                if ( dr < 0.3 && dr < dr_current ){
                    is_matched = true;
                    iemb_match = iemb;
                    dr_current = dr;
                }
    
            }
            if (!is_matched) { continue; }
    
            float et = emb_jet_energy_10->at(iemb_match);
            float et_cemc = emb_jet_energy_cemc_10->at(iemb_match);
            float et_hcalin = emb_jet_energy_hcalin_10->at(iemb_match);
            float et_hcalout = emb_jet_energy_hcalout_10->at(iemb_match);
            float area = emb_jet_area_10->at(iemb_match);
            int ntowers = emb_jet_num_towers_10->at(iemb_match);
            int ntowers_cemc = emb_jet_num_towers_cemc_10->at(iemb_match);
            int ntowers_hcalin = emb_jet_num_towers_hcalin_10->at(iemb_match);
            int ntowers_hcalout = emb_jet_num_towers_hcalout_10->at(iemb_match);
    
            h2_et_ntruth_cent->Fill(xaxis_var, et, ncomp, weight_10);
            h2_et_ntruth_cemc_cent->Fill(xaxis_var, et_cemc, ncomp_cemc, weight_10);
            h2_et_ntruth_hcalin_cent->Fill(xaxis_var, et_hcalin, ncomp_hcalin, weight_10);
            h2_et_ntruth_hcalout_cent->Fill(xaxis_var, et_hcalout, ncomp_hcalout, weight_10);
            ntruth_matched++;

            p_et_ntruth[icourse_bin]->Fill(et, ncomp, weight_10);
            p_et_ntruth_cemc[icourse_bin]->Fill(et_cemc, ncomp_cemc, weight_10);
            p_et_ntruth_hcalin[icourse_bin]->Fill(et_hcalin, ncomp_hcalin, weight_10);
            p_et_ntruth_hcalout[icourse_bin]->Fill(et_hcalout, ncomp_hcalout, weight_10);

            
    
        }
        
    }
    
    // jet 30
    std::cout << "Looping over jet 30" << std::endl;
    for ( int i = 0; i < nentries; ++i ) {
       
        t_30->GetEntry(i);
       
        float xaxis_var =1.0*centrality_30;
        if ( xaxis_var <0 || xaxis_var > MAX_X_CENT ) { continue; }

        int icourse_bin = -1;
        for ( unsigned ibin = 0; ibin < N_X_COURSE_CENT_BINS; ibin++ ){
            if ( xaxis_var > X_COURSE_CENT_BINS[ibin] && xaxis_var <= X_COURSE_CENT_BINS[ibin+1] ){
                icourse_bin = ibin;
                break;
            }
        }
    
    
    
        unsigned int ntruth = sim_jet_eta_30->size();
        unsigned int nemb = emb_jet_eta_30->size();
        ntruth_total += ntruth;
        nemb_total += nemb;
        // match truth to embedded
        for ( unsigned itruth = 0; itruth < ntruth; itruth++) {
    
            float teta = sim_jet_eta_30->at(itruth);
            float tphi = sim_jet_phi_30->at(itruth);
            float tet = sim_jet_energy_30->at(itruth);
            if ( std::fabs(teta)> 0.6) { continue; }
            if ( tet < 5.0 ) { continue; }
            int ncomp = sim_jet_num_towers_30->at(itruth);
            int ncomp_cemc = sim_jet_num_towers_cemc_30->at(itruth);
            int ncomp_hcalin = sim_jet_num_towers_hcalin_30->at(itruth);
            int ncomp_hcalout = sim_jet_num_towers_hcalout_30->at(itruth);
            bool is_matched = false;   
            float dr_current = 1000.0; 
            int iemb_match = -1;
            for ( unsigned iemb = 0; iemb < nemb; iemb++ ) {
                float reta = emb_jet_eta_30->at(iemb);
                float rphi = emb_jet_phi_30->at(iemb);
                float ret  = emb_jet_energy_30->at(iemb);
                if ( ret < 10.0 ) { continue; }
                if ( std::fabs(reta)> 0.6) { continue; }
                float deta = teta - reta;
                float dphi = tphi - rphi;
                if ( dphi > TMath::Pi() ) { dphi -= 2*TMath::Pi(); }
                if ( dphi < -TMath::Pi() ) { dphi += 2*TMath::Pi(); }
                float dr = TMath::Sqrt(deta*deta + dphi*dphi);
                if ( dr < 0.3 && dr < dr_current ){
                    is_matched = true;
                    iemb_match = iemb;
                    dr_current = dr;
                }
    
            }
            if (!is_matched) { continue; }
    
            float et = emb_jet_energy_30->at(iemb_match);
            float et_cemc = emb_jet_energy_cemc_30->at(iemb_match);
            float et_hcalin = emb_jet_energy_hcalin_30->at(iemb_match);
            float et_hcalout = emb_jet_energy_hcalout_30->at(iemb_match);
            float area = emb_jet_area_30->at(iemb_match);
            int ntowers = emb_jet_num_towers_30->at(iemb_match);
            int ntowers_cemc = emb_jet_num_towers_cemc_30->at(iemb_match);
            int ntowers_hcalin = emb_jet_num_towers_hcalin_30->at(iemb_match);
            int ntowers_hcalout = emb_jet_num_towers_hcalout_30->at(iemb_match);
    
            h2_et_ntruth_cent->Fill(xaxis_var, et, ncomp, weight_30);
            h2_et_ntruth_cemc_cent->Fill(xaxis_var, et_cemc, ncomp_cemc, weight_30);
            h2_et_ntruth_hcalin_cent->Fill(xaxis_var, et_hcalin, ncomp_hcalin, weight_30);
            h2_et_ntruth_hcalout_cent->Fill(xaxis_var, et_hcalout, ncomp_hcalout, weight_30);
            ntruth_matched++;

            p_et_ntruth[icourse_bin]->Fill(et, ncomp, weight_30);
            p_et_ntruth_cemc[icourse_bin]->Fill(et_cemc, ncomp_cemc, weight_30);
            p_et_ntruth_hcalin[icourse_bin]->Fill(et_hcalin, ncomp_hcalin, weight_30);
            p_et_ntruth_hcalout[icourse_bin]->Fill(et_hcalout, ncomp_hcalout, weight_30);
    
        }
        
    }

    std::cout << "Total truth jets: " << ntruth_total << " matched: " << ntruth_matched << " total emb jets: " << nemb_total << std::endl;
    
 
    TH1F * h1_et_vs_ntruth[N_X_CENT_BINS];
    TH1F * h1_et_vs_ntruth_cemc[N_X_CENT_BINS];
    TH1F * h1_et_vs_ntruth_hcalin[N_X_CENT_BINS];
    TH1F * h1_et_vs_ntruth_hcalout[N_X_CENT_BINS];
    TH2F * h2s_cent[N_X_CENT_BINS];
    TH2F * h2s_cemc[N_X_CENT_BINS];
    TH2F * h2s_hcalin[N_X_CENT_BINS];
    TH2F * h2s_hcalout[N_X_CENT_BINS];
    for ( int ibin = 0; ibin < N_X_CENT_BINS; ibin++ ) {
        h2_et_ntruth_cent->GetXaxis()->SetRange(ibin+1, ibin+1);
        TH2F * h2 = (TH2F*)h2_et_ntruth_cent->Project3D("yz");
        h2->SetName(Form("h2_et_vs_ntruth_%d", ibin));
        h2s_cent[ibin] = h2;
        TH1F * h1 = (TH1F*)h2->ProfileY(Form("h1_et_vs_ntruth_%d", ibin));
        h1->SetName(Form("h1_et_vs_ntruth_%d", ibin));
        h1_et_vs_ntruth[ibin] = h1;

        h2_et_ntruth_cemc_cent->GetXaxis()->SetRange(ibin+1, ibin+1);
        TH2F * h2_cemc = (TH2F*)h2_et_ntruth_cemc_cent->Project3D("yz");
        h2_cemc->SetName(Form("h2_et_ntruth_cemc_cent_%d", ibin));
        h2s_cemc[ibin] = h2_cemc;
        TH1F * h1_cemc = (TH1F*)h2_cemc->ProfileY(Form("h1_et_vs_ntruth_cemc_%d", ibin));
        h1_cemc->SetName(Form("h1_et_vs_ntruth_cemc_%d", ibin));
        h1_et_vs_ntruth_cemc[ibin] = h1_cemc;

        h2_et_ntruth_hcalin_cent->GetXaxis()->SetRange(ibin+1, ibin+1);
        TH2F * h2_hcalin = (TH2F*)h2_et_ntruth_hcalin_cent->Project3D("yz");
        h2_hcalin->SetName(Form("h2_et_ntruth_hcalin_cent_%d", ibin));
        h2s_hcalin[ibin] = h2_hcalin;
        TH1F * h1_hcalin = (TH1F*)h2_hcalin->ProfileY(Form("h1_et_vs_ntruth_hcalin_%d", ibin));
        h1_hcalin->SetName(Form("h1_et_vs_ntruth_hcalin_%d", ibin));
        h1_et_vs_ntruth_hcalin[ibin] = h1_hcalin;
        
        h2_et_ntruth_hcalout_cent->GetXaxis()->SetRange(ibin+1, ibin+1);
        TH2F * h2_hcalout = (TH2F*)h2_et_ntruth_hcalout_cent->Project3D("yz");
        h2_hcalout->SetName(Form("h2_et_ntruth_hcalout_cent_%d", ibin));
        h2s_hcalout[ibin] = h2_hcalout;
        TH1F * h1_hcalout = (TH1F*)h2_hcalout->ProfileY(Form("h1_et_vs_ntruth_hcalout_%d", ibin));
        h1_hcalout->SetName(Form("h1_et_vs_ntruth_hcalout_%d", ibin));
        h1_et_vs_ntruth_hcalout[ibin] = h1_hcalout;


    }
    std::cout << "Done with Nsignal" << std::endl;

    // jet 10
    std::cout << "Looping over jet 10" << std::endl;
    for ( int i = 0; i < nentries; ++i ) {
        t_10->GetEntry(i);
           
        float xaxis_var =1.0*centrality_10;
        if ( xaxis_var <0 || xaxis_var > MAX_X_CENT ) { continue; }
    
    
        int THIS_BIN = h2_area_res_vs_x->GetXaxis()->FindBin(xaxis_var);
        THIS_BIN--;// binning starts at 1
        if ( THIS_BIN < 0 || THIS_BIN >= N_X_CENT_BINS ) { continue; }
        // std::cout << xaxis_var << " is in bin " << THIS_BIN << std::endl;
    
    
        TH1F * h1_nsig = h1_et_vs_ntruth[THIS_BIN];
        if (!h1_nsig) { std::cout << "h1_nsig is null" << std::endl; continue; }
        TH1F * h1_nsig_cemc = h1_et_vs_ntruth_cemc[THIS_BIN];
        if (!h1_nsig_cemc) { std::cout << "h1_nsig_cemc is null" << std::endl; continue; }
        TH1F * h1_nsig_hcalin = h1_et_vs_ntruth_hcalin[THIS_BIN];
        if (!h1_nsig_hcalin) { std::cout << "h1_nsig_hcalin is null" << std::endl; continue; }
        TH1F * h1_nsig_hcalout = h1_et_vs_ntruth_hcalout[THIS_BIN];
        if (!h1_nsig_hcalout) { std::cout << "h1_nsig_hcalout is null" << std::endl; continue; }
    
        unsigned int ntruth = sim_jet_eta_10->size();
        unsigned int nemb = emb_jet_eta_10->size();

    
        // match truth to embedded
        for ( unsigned itruth = 0; itruth < ntruth; itruth++) {

            float teta = sim_jet_eta_10->at(itruth);
            float tphi = sim_jet_phi_10->at(itruth);
            float tet = sim_jet_energy_10->at(itruth);
            if ( std::fabs(teta)> 0.6) { continue; }
            if ( tet < 5.0 ) { continue; }
            int ncomp = sim_jet_num_towers_10->at(itruth);
            bool is_matched = false;   
            float dr_current = 1000.0; 
            int iemb_match = -1;
            for ( unsigned iemb = 0; iemb < nemb; iemb++ ) {
                float reta = emb_jet_eta_10->at(iemb);
                float rphi = emb_jet_phi_10->at(iemb);
                float ret  = emb_jet_energy_10->at(iemb);
                if ( ret < 10.0 ) { continue; }
                if ( std::fabs(reta)> 0.6) { continue; }
                float deta = teta - reta;
                float dphi = tphi - rphi;
                if ( dphi > TMath::Pi() ) { dphi -= 2*TMath::Pi(); }
                if ( dphi < -TMath::Pi() ) { dphi += 2*TMath::Pi(); }
                float dr = TMath::Sqrt(deta*deta + dphi*dphi);
                if ( dr < 0.3 && dr < dr_current ){
                    is_matched = true;
                    iemb_match = iemb;
                    dr_current = dr;
                }
    
            }
            if (!is_matched) { continue; }
    
            float et = emb_jet_energy_10->at(iemb_match);
            float et_cemc = emb_jet_energy_cemc_10->at(iemb_match);
            float et_hcalin = emb_jet_energy_hcalin_10->at(iemb_match);
            float et_hcalout = emb_jet_energy_hcalout_10->at(iemb_match);
            float area = emb_jet_area_10->at(iemb_match);
            int ntowers = emb_jet_num_towers_10->at(iemb_match);
            int ntowers_cemc = emb_jet_num_towers_cemc_10->at(iemb_match);
            int ntowers_hcalin = emb_jet_num_towers_hcalin_10->at(iemb_match);
            int ntowers_hcalout = emb_jet_num_towers_hcalout_10->at(iemb_match);

            int et_bin = h1_nsig->GetXaxis()->FindBin(et);
            if (et_bin < 1 || et_bin > N_ET_BINS ) { std::cout << "et_bin: (et) " << et_bin << " (" << et << ")" << std::endl; continue; }
            float ncorr = h1_nsig->GetBinContent(et_bin);

            int et_bin_cemc = h1_nsig_cemc->GetXaxis()->FindBin(et_cemc);
            float ncorr_cemc = 0;
            if (et_bin_cemc < 1 || et_bin_cemc > N_ET_BINS ) { 
                std::cout << "et_bin_cemc: (et_cemc) " << et_bin_cemc << " (" << et_cemc << ")" << std::endl;             
            } else {
                ncorr_cemc = h1_nsig_cemc->GetBinContent(et_bin_cemc);
            }
            
            int et_bin_hcalin = h1_nsig_hcalin->GetXaxis()->FindBin(et_hcalin);
            float ncorr_hcalin = 0;
            if (et_bin_hcalin < 1 || et_bin_hcalin > N_ET_BINS ) { 
                std::cout << "et_bin_hcalin: (et_hcalin) " << et_bin_hcalin << " (" << et_hcalin << ")" << std::endl; 
            } else {
                ncorr_hcalin = h1_nsig_hcalin->GetBinContent(et_bin_hcalin);
            }

            int et_bin_hcalout = h1_nsig_hcalout->GetXaxis()->FindBin(et_hcalout);
            float ncorr_hcalout = 0;
            if (et_bin_hcalout < 1 || et_bin_hcalout > N_ET_BINS ) { 
                std::cout << "et_bin_hcalout: (et_hcalout) " << et_bin_hcalout << " (" << et_hcalout << ")" << std::endl; 
            } else {
                ncorr_hcalout = h1_nsig_hcalout->GetBinContent(et_bin_hcalout);
            }
            // if(ncorr_cemc == 0 || ncorr_hcalin == 0 || ncorr_hcalout == 0) {
            //     // if any of the corrections are zero, we can't use them
            //     std::cout << "Warning: Zero correction for et_cemc/hcalin/hcalout. Skipping." << std::endl;
            //     continue;
            // }

            float rhoM = rho_val_TowerRho_MULT_10;
            float rhoA = rho_val_TowerRho_AREA_10;
            float rhoM_cemc = rho_val_TowerRho_MULT_CEMC_10;
            float rhoA_cemc = rho_val_TowerRho_AREA_CEMC_10;
            float rhoM_hcalin = rho_val_TowerRho_MULT_HCALIN_10;
            float rhoA_hcalin = rho_val_TowerRho_AREA_HCALIN_10;
            float rhoM_hcalout = rho_val_TowerRho_MULT_HCALOUT_10;
            float rhoA_hcalout = rho_val_TowerRho_AREA_HCALOUT_10;

            float mult_et = et - (rhoM*(ntowers - ncorr));
            // float mult_et = et - (rhoM_hcalout * (ntowers_hcalout - ncorr_hcalout) 
            //                         + rhoM_hcalin * (ntowers_hcalin - ncorr_hcalin)
            //                         + rhoM_cemc * (ntowers_cemc - ncorr_cemc));
         
            float area_et = et - rhoA*area;

            float cone_res_area = area_et - tet;
            float cone_res_mult = mult_et - tet;

            if ( area_et > 10 ) {
                h2_area_res_vs_x->Fill(xaxis_var, cone_res_area, weight_10);
                h2_rhoA_vs_x->Fill(xaxis_var, rhoA, weight_10);
            }
            
            if ( mult_et > 10 ) {
                h2_mult_res_vs_x->Fill(xaxis_var, cone_res_mult, weight_10);
                h2_rhoM_vs_x->Fill(xaxis_var, rhoM, weight_10);
            }

        } // end loop over truth jets

        // match truth to embedded sub1
        unsigned int ntruth_retower = sim_retower_jet_eta_10->size();
        unsigned int nsub1 = emb_jet_sub1_eta_10->size();
        for ( unsigned itruth = 0; itruth < ntruth_retower; itruth++) {
            float teta = sim_retower_jet_eta_10->at(itruth);
            float tphi = sim_retower_jet_phi_10->at(itruth);
            float tet = sim_retower_jet_energy_10->at(itruth);
            if ( std::fabs(teta)> 0.6) { continue; }
            if ( tet < 5.0 ) { continue; }
            int ncomp = sim_retower_jet_num_towers_10->at(itruth);
            bool is_matched = false;   
            float dr_current = 1000.0; 
            int iemb_match = -1;
            for ( unsigned iemb = 0; iemb < nsub1; iemb++ ) {
                float reta = emb_jet_sub1_eta_10->at(iemb);
                float rphi = emb_jet_sub1_phi_10->at(iemb);
                float ret  = emb_jet_sub1_energy_10->at(iemb);
                if ( ret < 10.0 ) { continue; }
                if ( std::fabs(reta)> 0.6) { continue; }
                float deta = teta - reta;
                float dphi = tphi - rphi;
                if ( dphi > TMath::Pi() ) { dphi -= 2*TMath::Pi(); }
                if ( dphi < -TMath::Pi() ) { dphi += 2*TMath::Pi(); }
                float dr = TMath::Sqrt(deta*deta + dphi*dphi);
                if ( dr < 0.3 && dr < dr_current ){
                    is_matched = true;
                    iemb_match = iemb;
                    dr_current = dr;
                }

            }
            if (!is_matched) { continue; }

            float et = emb_jet_sub1_energy_10->at(iemb_match);
            float res_sub1 = et - tet;
            if (et > 10 ) {
                h2_sub1_res_vs_x->Fill(xaxis_var, res_sub1, weight_10);
            }
        } // end loop over truth jets (retower)

    
    } // end loop over events
    
    // jet 30
    std::cout << "Looping over jet 30" << std::endl;
    for ( int i = 0; i < nentries; ++i ) {
       
        t_30->GetEntry(i);
           
        float xaxis_var = 1.0*centrality_30;
        if ( xaxis_var <0 || xaxis_var > MAX_X_CENT ) { continue; }
    
    
        int THIS_BIN = h2_area_res_vs_x->GetXaxis()->FindBin(xaxis_var);
        THIS_BIN--;// binning starts at 1
        if ( THIS_BIN < 0 || THIS_BIN >= N_X_CENT_BINS ) { continue; }
    
    
    
        TH1F * h1_nsig = h1_et_vs_ntruth[THIS_BIN];
        if (!h1_nsig) { std::cout << "h1_nsig is null" << std::endl; continue; }
        TH1F * h1_nsig_cemc = h1_et_vs_ntruth_cemc[THIS_BIN];
        if (!h1_nsig_cemc) { std::cout << "h1_nsig_cemc is null" << std::endl; continue; }
        TH1F * h1_nsig_hcalin = h1_et_vs_ntruth_hcalin[THIS_BIN];
        if (!h1_nsig_hcalin) { std::cout << "h1_nsig_hcalin is null" << std::endl; continue; }
        TH1F * h1_nsig_hcalout = h1_et_vs_ntruth_hcalout[THIS_BIN];
        if (!h1_nsig_hcalout) { std::cout << "h1_nsig_hcalout is null" << std::endl; continue; }
    
        unsigned int ntruth = sim_jet_eta_30->size();
        unsigned int nemb = emb_jet_eta_30->size();
        unsigned int ntruth_retower = sim_retower_jet_eta_30->size();
        unsigned int nsub1 = emb_jet_sub1_eta_30->size();
    
        // match truth to embedded
        for ( unsigned itruth = 0; itruth < ntruth; itruth++) {
    
            float teta = sim_jet_eta_30->at(itruth);
            float tphi = sim_jet_phi_30->at(itruth);
            float tet = sim_jet_energy_30->at(itruth);
            if ( std::fabs(teta)> 0.6) { continue; }
            if ( tet < 5.0 ) { continue; }
            int ncomp = sim_jet_num_towers_30->at(itruth);
            bool is_matched = false;   
            float dr_current = 1000.0; 
            int iemb_match = -1;
            for ( unsigned iemb = 0; iemb < nemb; iemb++ ) {
                float reta = emb_jet_eta_30->at(iemb);
                float rphi = emb_jet_phi_30->at(iemb);
                float ret  = emb_jet_energy_30->at(iemb);
                if ( ret < 10.0 ) { continue; }
                if ( std::fabs(reta)> 0.6) { continue; }
                float deta = teta - reta;
                float dphi = tphi - rphi;
                if ( dphi > TMath::Pi() ) { dphi -= 2*TMath::Pi(); }
                if ( dphi < -TMath::Pi() ) { dphi += 2*TMath::Pi(); }
                float dr = TMath::Sqrt(deta*deta + dphi*dphi);
                if ( dr < 0.3 && dr < dr_current ){
                    is_matched = true;
                    iemb_match = iemb;
                    dr_current = dr;
                }
    
            }
            if (!is_matched) { continue; }
    
            float et = emb_jet_energy_30->at(iemb_match);
            float et_cemc = emb_jet_energy_cemc_30->at(iemb_match);
            float et_hcalin = emb_jet_energy_hcalin_30->at(iemb_match);
            float et_hcalout = emb_jet_energy_hcalout_30->at(iemb_match);
            float area = emb_jet_area_30->at(iemb_match);
            int ntowers = emb_jet_num_towers_30->at(iemb_match);
            int ntowers_cemc = emb_jet_num_towers_cemc_30->at(iemb_match);
            int ntowers_hcalin = emb_jet_num_towers_hcalin_30->at(iemb_match);
            int ntowers_hcalout = emb_jet_num_towers_hcalout_30->at(iemb_match);
    
            int et_bin = h1_nsig->GetXaxis()->FindBin(et);
            if (et_bin < 1 || et_bin > N_ET_BINS ) { std::cout << "et_bin: (et) " << et_bin << " (" << et << ")" << std::endl; continue; }
            float ncorr = h1_nsig->GetBinContent(et_bin);
    
            int et_bin_cemc = h1_nsig_cemc->GetXaxis()->FindBin(et_cemc);
            float ncorr_cemc = 0;
            if (et_bin_cemc < 1 || et_bin_cemc > N_ET_BINS ) { 
                std::cout << "et_bin_cemc: (et_cemc) " << et_bin_cemc << " (" << et_cemc << ")" << std::endl;             
            } else {
                ncorr_cemc = h1_nsig_cemc->GetBinContent(et_bin_cemc);
            }
            
            int et_bin_hcalin = h1_nsig_hcalin->GetXaxis()->FindBin(et_hcalin);
            float ncorr_hcalin = 0;
            if (et_bin_hcalin < 1 || et_bin_hcalin > N_ET_BINS ) { 
                std::cout << "et_bin_hcalin: (et_hcalin) " << et_bin_hcalin << " (" << et_hcalin << ")" << std::endl; 
            } else {
                ncorr_hcalin = h1_nsig_hcalin->GetBinContent(et_bin_hcalin);
            }

            int et_bin_hcalout = h1_nsig_hcalout->GetXaxis()->FindBin(et_hcalout);
            float ncorr_hcalout = 0;
            if (et_bin_hcalout < 1 || et_bin_hcalout > N_ET_BINS ) { 
                std::cout << "et_bin_hcalout: (et_hcalout) " << et_bin_hcalout << " (" << et_hcalout << ")" << std::endl; 
            } else {
                ncorr_hcalout = h1_nsig_hcalout->GetBinContent(et_bin_hcalout);
            }

            // if(ncorr_cemc == 0 || ncorr_hcalin == 0 || ncorr_hcalout == 0) {
            //     // if any of the corrections are zero, we can't use them
            //     std::cout << "Warning: Zero correction for et_cemc/hcalin/hcalout. Skipping." << std::endl;
            //     continue;
            // }
    
            float rhoM = rho_val_TowerRho_MULT_30;
            float rhoA = rho_val_TowerRho_AREA_30;
            float rhoM_cemc = rho_val_TowerRho_MULT_CEMC_30;
            float rhoA_cemc = rho_val_TowerRho_AREA_CEMC_30;
            float rhoM_hcalin = rho_val_TowerRho_MULT_HCALIN_30;
            float rhoA_hcalin = rho_val_TowerRho_AREA_HCALIN_30;
            float rhoM_hcalout = rho_val_TowerRho_MULT_HCALOUT_30;
            float rhoA_hcalout = rho_val_TowerRho_AREA_HCALOUT_30;
    
 
            // float mult_et = et - (rhoM*(ntowers - ncorr));
            float mult_et = et - (rhoM*(ntowers - ncorr));
            // float mult_et = et - (rhoM_hcalout * (ntowers_hcalout - ncorr_hcalout) 
            //                         + rhoM_hcalin * (ntowers_hcalin - ncorr_hcalin)
            //                         + rhoM_cemc * (ntowers_cemc - ncorr_cemc));
         
           
            float area_et = et - rhoA*area;

            float cone_res_area = area_et - tet;
            float cone_res_mult = mult_et - tet;

            if ( area_et > 10 ) {
                h2_area_res_vs_x->Fill(xaxis_var, cone_res_area, weight_30);
                h2_rhoA_vs_x->Fill(xaxis_var, rhoA, weight_30);
            }

            if ( mult_et > 10 ) {
                h2_mult_res_vs_x->Fill(xaxis_var, cone_res_mult, weight_30);
                h2_rhoM_vs_x->Fill(xaxis_var, rhoM, weight_30);
            }


        
        } // end loop over truth jets
    
        // match truth to embedded sub1
        for ( unsigned itruth = 0; itruth < ntruth_retower; itruth++) {
            float teta = sim_retower_jet_eta_30->at(itruth);
            float tphi = sim_retower_jet_phi_30->at(itruth);
            float tet = sim_retower_jet_energy_30->at(itruth);
            if ( std::fabs(teta)> 0.6) { continue; }
            if ( tet < 5.0 ) { continue; }
            int ncomp = sim_retower_jet_num_towers_30->at(itruth);
            bool is_matched = false;   
            float dr_current = 1000.0; 
            int iemb_match = -1;
            for ( unsigned iemb = 0; iemb < nsub1; iemb++ ) {
                float reta = emb_jet_sub1_eta_30->at(iemb);
                float rphi = emb_jet_sub1_phi_30->at(iemb);
                float ret  = emb_jet_sub1_energy_30->at(iemb);
                if ( ret < 10.0 ) { continue; }
                if ( std::fabs(reta)> 0.6) { continue; }
                float deta = teta - reta;
                float dphi = tphi - rphi;
                if ( dphi > TMath::Pi() ) { dphi -= 2*TMath::Pi(); }
                if ( dphi < -TMath::Pi() ) { dphi += 2*TMath::Pi(); }
                float dr = TMath::Sqrt(deta*deta + dphi*dphi);
                if ( dr < 0.3 && dr < dr_current ){
                    is_matched = true;
                    iemb_match = iemb;
                    dr_current = dr;
                }
    
            }
            if (!is_matched) { continue; }
    
            float et = emb_jet_sub1_energy_30->at(iemb_match);
            float res_sub1 = et - tet;
            if (et > 10 ) { 
                h2_sub1_res_vs_x->Fill(xaxis_var, res_sub1, weight_30);
            }
        } // end loop over truth jets (retower)
    
    
    } // end loop over events

    
    TCanvas * c;
    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    // start with 3x3 of course x bins
    c = new TCanvas("c", "c", 400*3, 400*3);
    c->Divide(3,3);
    double tx=0.18;
    double ty_start=0.85;
    
    TLegend * leg = new TLegend(0.6,0.7,0.8,0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for ( unsigned islice = 0; islice < N_X_CENT_BINS; islice++ ){
        
        c->cd(islice+1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        

        double miny = 1e-3;
        int ihist = 0;
        std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
        std::string leg_title = Form("%0.0f-%0.0f %%", X_CENT_BINS[islice], X_CENT_BINS[islice+1]); 
        tags.push_back(leg_title);
        for ( auto h2 : h2s ) {
            h2->GetXaxis()->SetRange(islice+1, islice+1);
            TH1F * h1 = (TH1F*)h2->ProjectionY(Form("h1_%s_%d", h2->GetTitle(), islice));
            h1->SetLineColor(COLORS[ihist]);
            h1->SetMarkerColor(COLORS[ihist]);
            h1->SetMarkerStyle(MARKERS[ihist]);
            h1->Scale(1.0/h1->Integral());
            h1->GetYaxis()->SetRangeUser(miny, 1e0);
            h1->GetXaxis()->SetRangeUser(-60,50);
            // h1->GetXaxis()->SetTitle("#delta E_{T}^{Jet} [GeV]");
            h1->GetYaxis()->SetTitle("Probability Density [A.U.]");
            // yaxis offset
            h1->GetYaxis()->SetTitleOffset(1.5);
            if ( ihist == 0 ) { h1->Draw("p"); }
            else { h1->Draw("p same"); }
            if ( islice == 0 ) { leg->AddEntry(h1, labs[ihist].c_str(), "lp"); }
            ihist++;
        }

        double ty = ty_start;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.05;
        }

        leg->Draw("same");
    } 
    c->SaveAs((outdir+"/embedjet_res_cent_slices.png").c_str());

    delete c;
    delete leg;

    leg = new TLegend(0.18,0.5,0.35,0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    TGraphErrors * g_std_devs[h2s.size()];
    for ( unsigned ihist = 0; ihist < h2s.size(); ihist++ ) {
        g_std_devs[ihist] = new TGraphErrors(N_X_CENT_BINS);
    }
    for ( unsigned ibin = 0; ibin < N_X_CENT_BINS; ibin++ ) {
        c = new TCanvas("c", "c", 3*400, 400);
        c->Divide(3,1);
        double miny = 1e-5;
        int ihist = 0;
    
        for ( auto h2 : h2s ) {

            c->cd(ihist+1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.1);
            gPad->SetBottomMargin(0.15);
            gPad->SetTopMargin(0.05);

            std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
            // tags.push_back(labs[ihist]);
            std::string leg_title = Form("%s %0.0f-%0.0f %%", labs[ihist].c_str(),h2->GetXaxis()->GetBinLowEdge(ibin+1), h2->GetXaxis()->GetBinUpEdge(ibin+1)); 
            tags.push_back(leg_title);

            h2->GetXaxis()->SetRange(ibin+1, ibin+1);
            TH1F * h1 = (TH1F*)h2->ProjectionY(Form("h1_%s_%d", h2->GetTitle(), ibin));
            h1->SetLineColor(COLORS[ihist]);
            h1->SetMarkerColor(COLORS[ihist]);
            h1->SetMarkerStyle(MARKERS[ihist]);
            h1->Scale(1.0/h1->Integral());
            h1->GetYaxis()->SetRangeUser(miny, 1e0);
            h1->GetXaxis()->SetRangeUser(-60,50);

            float avg = h1->GetMean();
            float std = h1->GetRMS();

            // int mean_bin = h1->FindBin(avg);
            // h1->GetXaxis()->SetRange(1, mean_bin);
            // TF1 * f1 = new TF1("f1", "gaus", -abs_max_x, abs_max_x);
            // h1->Fit(f1, "RQ", "", -abs_max_x, avg); // left side of peak
            // float mean_left = f1->GetParameter(1);
            // float sigma_left = f1->GetParameter(2);
            // h1->Fit(f1, "RQ", "", mean_left - 1.5*sigma_left, avg); // left side of peak
            // mean_left = f1->GetParameter(1);
            // sigma_left = f1->GetParameter(2);
            // float mean_left_err = f1->GetParError(1);
            // float sigma_left_err = f1->GetParError(2);
            // TF1 * fitfunc = h1->GetFunction("f1");
            // fitfunc->SetLineColor(kRed);
            // fitfunc->SetLineStyle(2);
            // fitfunc->SetLineWidth(2);
            // float chi2 = fitfunc->GetChisquare();

            g_std_devs[ihist]->SetPoint(ibin, h2->GetXaxis()->GetBinCenter(ibin+1), std);
            g_std_devs[ihist]->SetPointError(ibin, h2->GetXaxis()->GetBinWidth(ibin+1)/2.0, 0);

            // h1->GetXaxis()->SetRange(1, h1->GetNbinsX());
            // TF1 * f2_gamma = new TF1("f2_gamma", "[0]*([1]/TMath::Gamma([2]))*TMath::Power([1]*x + [2], [2]-1)*TMath::Exp(-[1]*x - [2])",-1.2*h1->GetBinCenter(lastbin_above_threshold), 1.2*h1->GetBinCenter(lastbin_above_threshold));
            // h1->Fit(f2_gamma, "RQ", "", -abs_max_x, abs_max_x);
            // float gamma_ab = f2_gamma->GetParameter(1);
            // float gamma_ap = f2_gamma->GetParameter(2);
            // float gamma_norm = f2_gamma->GetParameter(0);
            // float gamma_mean = gamma_ap/gamma_ab;
            // float gamma_sigma = TMath::Sqrt(gamma_ap)/gamma_ab;  
            // h1->Fit(f2_gamma, "RQ", "", gamma_mean - 1.5*gamma_sigma, gamma_mean + 1.5*gamma_sigma);
            // gamma_ab = f2_gamma->GetParameter(1);
            // gamma_ap = f2_gamma->GetParameter(2);
            // gamma_norm = f2_gamma->GetParameter(0);
            // gamma_mean = gamma_ap/gamma_ab;
            // gamma_sigma = TMath::Sqrt(gamma_ap)/gamma_ab;
            // float gamma_ab_err = f2_gamma->GetParError(1);
            // float gamma_ap_err = f2_gamma->GetParError(2);
            // float gamma_mean_err = std::sqrt((gamma_ap_err/gamma_ab)*(gamma_ap_err/gamma_ab) + (gamma_ab_err/gamma_ab)*(gamma_ab_err/gamma_ab))*gamma_mean;
            // float gamma_sigma_err = std::sqrt((0.5*(gamma_ap_err/gamma_ap))*(0.5*(gamma_ap_err/gamma_ap)) + (gamma_ab_err/gamma_ab)*(gamma_ab_err/gamma_ab))*gamma_sigma;
            // TF1 * fitfunc_gamma = h1->GetFunction("f2_gamma");
            // fitfunc_gamma->SetLineColor(kAzure);
            // fitfunc_gamma->SetLineStyle(2);
            // fitfunc_gamma->SetLineWidth(2);
            // float chi2_gamma = fitfunc_gamma->GetChisquare();
            h1->Draw("p");
            leg->AddEntry(h1, Form("#mu = %0.2f, #sigma = %0.2f", avg, std), "p");
            // leg->AddEntry(fitfunc, Form("#mu_{LHS} = %0.2f #pm %0.2f, #sigma_{LHS} = %0.2f #pm %0.2f", mean_left, mean_left_err, sigma_left, sigma_left_err), "l");
            // leg->AddEntry(fitfunc_gamma, Form("a_{b} = %0.2f #pm %0.2f, a_{p} = %0.2f #pm %0.2f", gamma_ab, gamma_ab_err, gamma_ap, gamma_ap_err), "l");

            double ty = 0.87;
            for ( auto tag : tags ) {
                tex->DrawLatex(tx, ty, tag.c_str());
                ty -= 0.05;
            }

            leg->Draw("same");
            ihist++;
            c->Update();
            leg->Clear();
        }
        std::cout << "done with slice " << ibin << std::endl;

        c->SaveAs((outdir+"/embedjet_res_cent_slice_"+std::to_string(ibin)+".png").c_str());
        leg->Clear();
        delete c;
    }

    TLegend * leg2 = new TLegend(0.18,0.5,0.35,0.7);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    int mycolors[3] = {kRed, kBlue, kCyan};
    for ( unsigned ihist = 0; ihist < h2s.size(); ihist++ ) {
       
        g_std_devs[ihist]->SetMarkerStyle(MARKERS[ihist+1]);
        g_std_devs[ihist]->SetMarkerColor(mycolors[ihist]);
        g_std_devs[ihist]->SetMarkerSize(1.5);
        g_std_devs[ihist]->GetXaxis()->SetNdivisions(505);
        g_std_devs[ihist]->GetYaxis()->SetTitle("#sigma_{LHS} [GeV]");
        g_std_devs[ihist]->GetXaxis()->SetTitle("Centrality [%]"); 
        g_std_devs[ihist]->GetYaxis()->SetRangeUser(0.1, 12);
        if ( ihist == 0 ) { g_std_devs[ihist]->Draw("ap"); }
        else { g_std_devs[ihist]->Draw("p same"); }
        g_std_devs[ihist]->Draw("p same");
        leg2->AddEntry(g_std_devs[ihist], labs[ihist].c_str(), "p");
    }
    leg2->Draw("same");
    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
    double tx2=0.18;
    double ty2=0.89;
    for ( auto tag : tags ) {
        tex->DrawLatex(tx2, ty2, tag.c_str());
        ty2 -= 0.05;
    }
    c->SaveAs((outdir+"/res_vs_x.png").c_str());
    delete c;
    delete leg2;
    

    // mult curves
    c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);

    leg2 = new TLegend(0.18,0.5,0.35,0.7);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    for (  unsigned ibin = 0; ibin < N_X_COURSE_CENT_BINS; ibin++ ){
        std::string label = Form("%0.0f-%0.0f %%", X_COURSE_CENT_BINS[ibin], X_COURSE_CENT_BINS[ibin+1]);
        p_et_ntruth[ibin]->SetMarkerStyle(MARKERS[ibin]);
        p_et_ntruth[ibin]->SetMarkerColor(COLORS[ibin]);
        p_et_ntruth[ibin]->SetMarkerSize(1.5);
        p_et_ntruth[ibin]->GetXaxis()->SetNdivisions(505);
        p_et_ntruth[ibin]->SetLineColor(COLORS[ibin]);
        p_et_ntruth[ibin]->GetYaxis()->SetTitle("N_{Towers}^{Sim.}");
        p_et_ntruth[ibin]->GetXaxis()->SetTitle("E_{T}^{Uncorr.} [GeV]");
        if ( ibin == 0 ) { p_et_ntruth[ibin]->Draw("p"); }
        else { p_et_ntruth[ibin]->Draw("p same"); }
        leg2->AddEntry(p_et_ntruth[ibin], label.c_str(), "pe");
    }
    leg2->Draw("same");
    tags = {sPHENIX_Tag, DataType_Tag};
    tx2=0.18;
    ty2=0.89;
    for ( auto tag : tags ) {
        tex->DrawLatex(tx2, ty2, tag.c_str());
        ty2 -= 0.05;
    }
    c->SaveAs((outdir+"/et_vs_ntruth.png").c_str());
    delete c;
    delete leg2;







    TFile * fout = new TFile((outdir+"/embed.root").c_str(), "RECREATE");
    for ( auto h2 : h2s ) {
        h2->Write();
    }
    for ( auto h3 : h3s ) {
        h3->Write();
    }
    for ( unsigned ihist = 0; ihist < h2s.size(); ihist++ ) {
        g_std_devs[ihist]->SetName(labs[ihist].c_str());
        g_std_devs[ihist]->Write();
    }
    for ( int ibin = 0; ibin < N_X_CENT_BINS; ibin++ ) {
        h1_et_vs_ntruth[ibin]->Write();
        h2s_cent[ibin]->Write();
    }
    h2_rhoM_vs_x->Write();
    h2_rhoA_vs_x->Write();
    for ( int ibin = 0; ibin < N_X_COURSE_CENT_BINS; ibin++ ) {
        p_et_ntruth[ibin]->Write();
        p_et_ntruth_cemc[ibin]->Write();
        p_et_ntruth_hcalin[ibin]->Write();
        p_et_ntruth_hcalout[ibin]->Write();
    }

    fout->Close();
    
    f_10->Close();
    f_30->Close();
    return outdir+"/embed.root";
    
}

void ConfigureOutputDirs(std::string input_file_base, std::string plotting_dir)
{
    
    if ( !gSystem->OpenDirectory(plotting_dir.c_str()) ) {
        gSystem->mkdir(plotting_dir.c_str(), true);
    }

    plotting_dir += input_file_base + "/";
    if ( !gSystem->OpenDirectory(plotting_dir.c_str()) ) {
        gSystem->mkdir(plotting_dir.c_str(), true);
    }
    
    probe_plots = plotting_dir + "probes";
    std::vector<std::string> plot_directories = {probe_plots};
    for ( auto const& d : plot_directories ) {
        if ( !gSystem->OpenDirectory(d.c_str()) ) {
            gSystem->mkdir(d.c_str(), true);
        }
    }

    return;
}



