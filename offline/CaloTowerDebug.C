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

const float CENT_BINS[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 100};
const int N_CENT_BINS = sizeof(CENT_BINS)/sizeof(CENT_BINS[0]) - 1;

const int N_SUM_ET_BINS = 10;
float SUM_ET_BINS[N_SUM_ET_BINS+1];
float MAX_SUM_ET = 1200;

const int N_SUM_Q_BINS = 50;
float SUM_Q_BINS[N_SUM_Q_BINS+1];
float MAX_SUM_Q = 2200;

// const int N_X_CENT_BINS = 18;
// float X_CENT_BINS[N_X_CENT_BINS+1];
// float MAX_X_CENT = 80;

const float V2_VALUES[] = {2.32, 3.39, 4.76, 6.18, 7.03, 7.4, 7.44, 7.23, 6.96};
const float V3_VALUES[] = {1.45, 1.62, 1.76, 1.9, 1.99, 2.05, 1.92, 1.75, 1.57};
const float X_CENT_BINS[]= {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
const int N_X_CENT_BINS = sizeof(X_CENT_BINS)/sizeof(X_CENT_BINS[0]) - 1;
float MAX_X_CENT = 80;

const float X_COURSE_CENT_BINS[]= {0, 10, 20, 30, 40, 50, 60, 70, 80};
const int N_COURSE_X_BINS = sizeof(X_COURSE_CENT_BINS)/sizeof(X_COURSE_CENT_BINS[0]) - 1;

const int N_ZVTX_BINS = 100;
float ZVTX_BINS[N_ZVTX_BINS+1];
float MAX_ZVTX = 25;

const int N_RHO_M_BINS = 200;
float RHO_M_MAX = 0.1;
float RHO_M_BINS[N_RHO_M_BINS+1];

const int N_RHO_A_BINS = 200;
float RHO_A_MAX = 150;
float RHO_A_BINS[N_RHO_A_BINS+1];

const int N_BKGD_BINS = 200;
float BKGD_MAX = 2;
float BKGD_BINS[N_BKGD_BINS+1];

const int N_CONECOMP_BINS = 500;
float CONECOMP_MAX = 1100;
float CONECOMP_BINS[N_CONECOMP_BINS+1];

const int N_CONECOMP_SUB1_BINS = 100;
float CONECOMP_SUB1_MAX = 100;
float CONECOMP_SUB1_BINS[N_CONECOMP_SUB1_BINS+1];

const int N_CONE_DET_BINS = 250;
float MAX_CONE_DET = 60;
float CONE_DET_BINS[N_CONE_DET_BINS+1];

 
void SetBins(){
    for ( int i = 0; i < N_SUM_ET_BINS+1; ++i ) { SUM_ET_BINS[i] = i*MAX_SUM_ET/N_SUM_ET_BINS; }
    // for ( int i = 0; i < N_X_CENT_BINS+1; ++i ) { X_CENT_BINS[i] = i*MAX_X_CENT/N_X_CENT_BINS; }
    for ( int i = 0; i < N_SUM_Q_BINS+1; ++i ) { SUM_Q_BINS[i] = i*MAX_SUM_Q/N_SUM_Q_BINS; }
    for ( int i = 0; i < N_ZVTX_BINS+1; ++i ) { ZVTX_BINS[i] = -MAX_ZVTX + i*2*MAX_ZVTX/N_ZVTX_BINS; }
    for ( int i = 0; i < N_RHO_M_BINS+1; ++i ) { RHO_M_BINS[i] = i*RHO_M_MAX/N_RHO_M_BINS; }
    for ( int i = 0; i < N_RHO_A_BINS+1; ++i ) { RHO_A_BINS[i] = i*RHO_A_MAX/N_RHO_A_BINS; }
    for ( int i = 0; i < N_BKGD_BINS+1; ++i ) { BKGD_BINS[i] = i*BKGD_MAX/N_BKGD_BINS; }
    for ( int i = 0; i < N_CONECOMP_BINS+1; ++i ) { CONECOMP_BINS[i] = 500 + i*CONECOMP_MAX/N_CONECOMP_BINS; }
    for ( int i = 0; i < N_CONECOMP_SUB1_BINS+1; ++i ) { CONECOMP_SUB1_BINS[i] = i*CONECOMP_SUB1_MAX/N_CONECOMP_SUB1_BINS; }
    for ( int i = 0; i < N_CONE_DET_BINS+1; ++i ) { CONE_DET_BINS[i] = -MAX_CONE_DET + i*2*MAX_CONE_DET/N_CONE_DET_BINS; }
}

const float AREA_CONE = TMath::Pi()*0.4*0.4;
const float AREA_TOWER_CEMC = (2.0*TMath::Pi()/256.0)*(2.2/96.0);
const float AREA_HCAL_TOWER = (2.0*TMath::Pi()/64.0)*(2.2/24.0);

bool IS_DATA = false;
int NEVENTS = 0;
std::string DataType_Tag;


const int COLORS[] = {kBlack, kRed , kBlue, kGreen+2, kViolet, kCyan, kOrange+2, kMagenta+2, kAzure-2};
const int MARKERS[] = { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kOpenCircle, kOpenSquare, kOpenTriangleUp};
const float MARKER_SIZE = 1.2;
const float LINE_WIDTH = 2.0;

std::string global_plots;
std::string calo_plots;
std::string calo_window_plots;
std::string bkgd_plots;
std::string random_cone_plots;

void ConfigureOutputDirs(std::string input_file_base, std::string plotting_dir = "plots/");
std::string MakeGetDir( const std::string & dir ){
    if ( gSystem->AccessPathName(dir.c_str()) ) {
        gSystem->Exec(Form("mkdir -p %s", dir.c_str()));
    }
    return dir;
}

const unsigned int k_window_array_size = 11; 
const std::vector < std::pair < unsigned int, unsigned int > > k_calo_window_dims_cemc_geom = {
        {1,1}, {2,2}, {5,5}, {10,10}, {15,15}, {22,23}, {30,29}, {37,37}, {45,45}, {52,52}, {60,60}
};
const std::vector < std::pair < unsigned int, unsigned int > > k_calo_window_dims_hcal_geom = {
    {1,1}, {2,2}, {3,4}, {5,6}, {7,8}, {9,10}, {11,12}, {13,13}, {15,15}, {17,17}, {20,20}
}; 

std::pair<int,int> GetWindowDimFromString(std::string hist_name);
void CaloWindowMultiPanel(std::vector<TH2F*> h2_vec, const std::vector<std::string> h2_titles,
                        const std::string output_location, const std::string pngbase,
                        const std::string xaxis_title, const std::string yaxis_title, 
                        const double tagx, const double tagy, const double dty, 
                        const double lx1, const double lx2, const double ly1, const double ly2,
                        float minx, float maxx, const float miny, const float maxy,
                        bool sym_x = false );
void CaloWindowMultiPanel3D(std::vector<TH3F*> h3s, const std::vector<std::string> h3_titles,
                        const std::string output_location, const std::string pngbase,
                        const std::string xaxis_title, const std::string yaxis_title, // after projection
                        const double tagx, const double tagy, const double dty, 
                        const double lx1, const double lx2, const double ly1, const double ly2,
                        float minx, float maxx, const float miny, const float maxy,
                        bool sym_x = false, bool logy = false );

void ProcessCaloWindowHistograms(const std::string & input_file);
void ProcessCaloWindowTree(const std::string & input_file, bool x_axis_cent = true);
void ProcessBackgroundTree(const std::string & input_file, bool x_axis_cent = true);
void ProcessRandomConeTree(const std::string & input_file, bool x_axis_cent = true);
void ProcessGlobal(const std::string & input_file);
float CalcPoisson(const float sum2, const float sum, const float n, const float ncomp, const float ncones){
    if ( n == 0 ) { return 0;}
    if ( ncones == 0 ) { return 0;}
    float Na = ncomp/ncones;
    float mu = sum/n;
    float sigma2 = sum2/n - (mu*mu);
    float sigma = TMath::Sqrt(Na*sigma2  + Na*(mu*mu));
    return sigma;
}

float CalcPoissonHarm(const float sum2, const float sum, const float n, const float ncomp, const float ncones, const float v2 , const float v3 = 0){
    if ( n == 0 ) { return 0;}
    if ( ncones == 0 ) { return 0;}
    float Na = ncomp/ncones;
    float mu = sum/n;
    float sigma2 = sum2/n - (mu*mu);
    float vn2 = v2*v2 + v3*v3;
    float harmcont = Na*Na*vn2*2.0*mu*mu;
    float sigma = TMath::Sqrt(Na*sigma2  + Na*(mu*mu) + harmcont);
    return sigma;
}

void RunAll(const std::string & input_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/dsts/DATA/feb7_basic.root") {
   
    IS_DATA= true;
    DataType_Tag = IS_DATA ? "2024 Au+Au 200 GeV" : "HIJING MDC2";

    // get base name of input file
    std::string input_file_base = input_file;
    size_t found = input_file.find_last_of("/");
    if (found != std::string::npos) {
        input_file_base = input_file.substr(found+1);
    }
    found = input_file_base.find_last_of(".");
    if (found != std::string::npos) {
        input_file_base = input_file_base.substr(0, found);
    }

    
    ConfigureOutputDirs(input_file_base);
    std::cout << "Input file: " << input_file << std::endl;
    std::cout << "Global plots: " << global_plots << std::endl;
    std::cout << "Calo plots: " << calo_plots << std::endl;
    std::cout << "Calo window plots: " << calo_window_plots << std::endl;
    std::cout << "Background plots: " << bkgd_plots << std::endl;
    std::cout << "Random cone plots: " << random_cone_plots << std::endl;
    SetBins();
    
    ProcessGlobal(input_file);


    bool do_calo_window = false;
    bool do_background = true;
    bool do_random_cone = true;

    if ( do_background ) {
        ProcessBackgroundTree(input_file, true); // x-axis = cent
        ProcessBackgroundTree(input_file, false); // x-axis = SUM_Et
    }
    

    // if ( do_calo_window ) {
        ProcessCaloWindowTree(input_file, true); // x-axis = cent
        // ProcessCaloWindowTree(input_file, false); // x-axis = SUM_Et
        ProcessCaloWindowHistograms(input_file);
    // }
    
    return ;
   
}

void CaloWindowHistos() {

    


    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    SetsPhenixStyle();
    
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kRainBow);

    RunAll("/sphenix/user/tmengel/UE-AuAu-PPG04/dsts/DATA/feb7_basic.root");
    RunAll("/sphenix/user/tmengel/UE-AuAu-PPG04/dsts/DATA/feb7_random.root");
    gSystem->Exit(0);

}


void ProcessGlobal( const std::string & input_file )
{
    auto outdir_base = MakeGetDir(global_plots);
    auto outdir_png = MakeGetDir(outdir_base + "/png");
    auto outdir_C = MakeGetDir(outdir_base + "/C"); 
    auto outdir_rootifles = MakeGetDir(outdir_base + "/rootfiles");
    
    std::cout <<"ProcessGlobal: Begin" << std::endl;

    TFile * f = new TFile(input_file.c_str(), "READ");
    if( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl;  exit(1); }

    TH1F * h1_num_events = (TH1F*)f->Get("h1_num_events");
    if( !h1_num_events ){ std::cout << "h1_num_events not found!" << std::endl; exit(1); }
    NEVENTS = (int)h1_num_events->GetBinContent(1);
    std::cout << "Number of events: " << NEVENTS << std::endl;

    TTree * t = (TTree*)f->Get("T");
    
    // tree branches 
        float mbd_q_N = 0;
        float mbd_q_S = 0;
        float mbd_time_N = 0;
        float mbd_time_S = 0;
        t->SetBranchAddress("mbd_q_N", &mbd_q_N);
        t->SetBranchAddress("mbd_q_S", &mbd_q_S);
        t->SetBranchAddress("mbd_time_N", &mbd_time_N);
        t->SetBranchAddress("mbd_time_S", &mbd_time_S);

        float zvtx = 0;
        t->SetBranchAddress("zvtx", &zvtx);

        int centrality = 0;
        t->SetBranchAddress("centrality", &centrality);

        std::vector< float > * tower_background_energy_recemc = 0;
        std::vector< float > * tower_background_energy_hcalin = 0;
        std::vector< float > * tower_background_energy_hcalout = 0;
        t->SetBranchAddress("tower_background_energy_recemc", &tower_background_energy_recemc);
        t->SetBranchAddress("tower_background_energy_hcalin", &tower_background_energy_hcalin);
        t->SetBranchAddress("tower_background_energy_hcalout", &tower_background_energy_hcalout);


        float rho_val_TowerRho_AREA = 0;
        float std_rho_val_TowerRho_AREA = 0;
        t->SetBranchAddress("rho_val_TowerRho_AREA", &rho_val_TowerRho_AREA);
        t->SetBranchAddress("std_rho_val_TowerRho_AREA", &std_rho_val_TowerRho_AREA);
        
        float rho_val_TowerRho_MULT = 0;
        float std_rho_val_TowerRho_MULT = 0;
        t->SetBranchAddress("rho_val_TowerRho_MULT", &rho_val_TowerRho_MULT);
        t->SetBranchAddress("std_rho_val_TowerRho_MULT", &std_rho_val_TowerRho_MULT);

        float rho_val_TowerRho_AREA_CEMC = 0;
        float std_rho_val_TowerRho_AREA_CEMC = 0;
        t->SetBranchAddress("rho_val_TowerRho_AREA_CEMC", &rho_val_TowerRho_AREA_CEMC);
        t->SetBranchAddress("std_rho_val_TowerRho_AREA_CEMC", &std_rho_val_TowerRho_AREA_CEMC);

        float rho_val_TowerRho_MULT_CEMC = 0;
        float std_rho_val_TowerRho_MULT_CEMC = 0;
        t->SetBranchAddress("rho_val_TowerRho_MULT_CEMC", &rho_val_TowerRho_MULT_CEMC);
        t->SetBranchAddress("std_rho_val_TowerRho_MULT_CEMC", &std_rho_val_TowerRho_MULT_CEMC);

        float rho_val_TowerRho_AREA_HCALIN = 0;
        float std_rho_val_TowerRho_AREA_HCALIN = 0;
        t->SetBranchAddress("rho_val_TowerRho_AREA_HCALIN", &rho_val_TowerRho_AREA_HCALIN);
        t->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALIN", &std_rho_val_TowerRho_AREA_HCALIN);

        float rho_val_TowerRho_MULT_HCALIN = 0;
        float std_rho_val_TowerRho_MULT_HCALIN = 0;
        t->SetBranchAddress("rho_val_TowerRho_MULT_HCALIN", &rho_val_TowerRho_MULT_HCALIN);
        t->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALIN", &std_rho_val_TowerRho_MULT_HCALIN);

        float rho_val_TowerRho_AREA_HCALOUT = 0;
        float std_rho_val_TowerRho_AREA_HCALOUT = 0;
        t->SetBranchAddress("rho_val_TowerRho_AREA_HCALOUT", &rho_val_TowerRho_AREA_HCALOUT);
        t->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALOUT", &std_rho_val_TowerRho_AREA_HCALOUT);

        float rho_val_TowerRho_MULT_HCALOUT = 0;
        float std_rho_val_TowerRho_MULT_HCALOUT = 0;
        t->SetBranchAddress("rho_val_TowerRho_MULT_HCALOUT", &rho_val_TowerRho_MULT_HCALOUT);
        t->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALOUT", &std_rho_val_TowerRho_MULT_HCALOUT);

        float tower_frac_fired_TOWERINFO_CALIB_CEMC = 0;
        float tower_frac_dead_TOWERINFO_CALIB_CEMC = 0;
        float tower_avg_energy_TOWERINFO_CALIB_CEMC = 0;
        float tower_std_energy_TOWERINFO_CALIB_CEMC = 0;
        float tower_sum_energy_TOWERINFO_CALIB_CEMC = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC", &tower_frac_fired_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC", &tower_frac_dead_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC", &tower_avg_energy_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC", &tower_std_energy_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC", &tower_sum_energy_TOWERINFO_CALIB_CEMC);

        float tower_frac_fired_TOWERINFO_CALIB_HCALIN = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALIN = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALIN = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALIN = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALIN = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN", &tower_frac_fired_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN", &tower_frac_dead_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN", &tower_avg_energy_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN", &tower_std_energy_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN", &tower_sum_energy_TOWERINFO_CALIB_HCALIN);

        float tower_frac_fired_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALOUT = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT", &tower_frac_fired_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT", &tower_frac_dead_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT", &tower_avg_energy_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT", &tower_std_energy_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT", &tower_sum_energy_TOWERINFO_CALIB_HCALOUT);

        float tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER", &tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER", &tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER);

        float tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);

        float tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1", &tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1", &tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1);

        float tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1);

    int nentries = t->GetEntries();

 
    TH2F * h2_SUM_Et_vs_cent = new TH2F("h2_SUM_Et_vs_cent", "SUM_Et_vs_cent", N_X_CENT_BINS, X_CENT_BINS, N_SUM_ET_BINS, SUM_ET_BINS);
    h2_SUM_Et_vs_cent->GetXaxis()->SetTitle("Centrality [%]");
    h2_SUM_Et_vs_cent->GetYaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
    TH2F * h2_sumMBDq_vs_cent = new TH2F("h2_sumMBDq_vs_cent", "sumMBDq_vs_cent", N_X_CENT_BINS, X_CENT_BINS, N_SUM_ET_BINS, SUM_ET_BINS);
    h2_sumMBDq_vs_cent->GetXaxis()->SetTitle("Centrality [%]");
    h2_sumMBDq_vs_cent->GetYaxis()->SetTitle("#Sigma Q_{MBD}");
    TH2F * h2_sumMBDq_vs_SUM_Et = new TH2F("h2_sumMBDq_vs_SUM_Et", "sumMBDq_vs_SUM_Et", N_SUM_ET_BINS, SUM_ET_BINS, N_SUM_ET_BINS, SUM_ET_BINS);
    h2_sumMBDq_vs_SUM_Et->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
    h2_sumMBDq_vs_SUM_Et->GetYaxis()->SetTitle("#Sigma Q_{MBD}");

    TH1F * h1_SUM_Et = new TH1F("h1_SUM_Et", "SUM_Et", N_SUM_ET_BINS, SUM_ET_BINS);
    h1_SUM_Et->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
    h1_SUM_Et->GetYaxis()->SetTitle("1/N_{events} dN/d#Sigma E_{T}^{Raw} [GeV^{-1}]");
    TH1F * h1_sumMBDq = new TH1F("h1_sumMBDq", "sumMBDq", N_SUM_ET_BINS, SUM_ET_BINS);
    h1_sumMBDq->GetXaxis()->SetTitle("#Sigma Q_{MBD}");
    h1_sumMBDq->GetYaxis()->SetTitle("1/N_{events} dN/d#Sigma Q_{MBD}");
    TH1F * h1_cent = new TH1F("h1_cent", "cent", N_X_CENT_BINS, X_CENT_BINS);
    h1_cent->GetXaxis()->SetTitle("Centrality [%]");
    h1_cent->GetYaxis()->SetTitle("Probability Density [a.u.]");
    TH1F * h1_zvtx = new TH1F("h1_zvtx", "zvtx", N_ZVTX_BINS, ZVTX_BINS);
    h1_zvtx->GetXaxis()->SetTitle("z_{vtx} [cm]");
    h1_zvtx->GetYaxis()->SetTitle("Probability Density [a.u.]");


    for ( int i = 0; i < nentries; ++i ) {
        t->GetEntry(i);
        float sum_et = tower_sum_energy_TOWERINFO_CALIB_CEMC + tower_sum_energy_TOWERINFO_CALIB_HCALIN + tower_sum_energy_TOWERINFO_CALIB_HCALOUT;
        float sum_mbdq = mbd_q_N + mbd_q_S;
        h2_SUM_Et_vs_cent->Fill(centrality, sum_et);
        h2_sumMBDq_vs_cent->Fill(centrality, sum_mbdq);
        h2_sumMBDq_vs_SUM_Et->Fill(sum_et, sum_mbdq);
        h1_SUM_Et->Fill(sum_et);
        h1_sumMBDq->Fill(sum_mbdq);
        h1_cent->Fill(centrality);
        h1_zvtx->Fill(zvtx);
    }

    TCanvas * c;
    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    double tx=0.45;
    double ty_start=0.85;
    std::vector<TH2F*> h2s = {h2_SUM_Et_vs_cent, h2_sumMBDq_vs_cent};
    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
    std::string nevents = Form("%0.0e Events", 1.0*NEVENTS);
    tags.push_back(nevents);
    for ( auto h2 : h2s ) {
        double ty = ty_start;
        c = new TCanvas("c", "c", 800, 600);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        h2->GetXaxis()->SetNdivisions(505);
        h2->GetYaxis()->SetNdivisions(505);
       
        h2->Draw("colz");
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.05;
        }
        c->SaveAs((outdir_png+"/"+h2->GetTitle()+".png").c_str());
        c->SaveAs((outdir_C+"/"+h2->GetTitle()+".C").c_str());
        delete c;
    }

    tx=0.18;
    h2s= {h2_sumMBDq_vs_SUM_Et};
    for ( auto h2 : h2s ) {
        double ty = ty_start;
        c = new TCanvas("c", "c", 800, 600);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        h2->GetXaxis()->SetNdivisions(505);
        h2->GetYaxis()->SetNdivisions(505);
       
        h2->Draw("colz");
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.05;
        }

        c->SaveAs((outdir_C+"/"+h2->GetTitle()+".C").c_str());
        c->SaveAs((outdir_png+"/"+h2->GetTitle()+".png").c_str());
        delete c;
    }

    tx = 0.45;
    std::vector<TH1F*> h1s = {h1_SUM_Et, h1_sumMBDq, h1_cent, h1_zvtx};
    std::vector<std::string> h1_titles = {"h1_SUM_Et", "h1_sumMBDq", "h1_cent", "h1_zvtx"};
    for ( unsigned int i = 0; i < h1s.size(); ++i ) {
        c = new TCanvas("c", "c", 800, 600);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogy();
        h1s[i]->GetXaxis()->SetNdivisions(505);
        h1s[i]->GetYaxis()->SetNdivisions(505);
        
        h1s[i]->Scale(1.0/NEVENTS);
        h1s[i]->GetYaxis()->SetRangeUser(1e-3,1e0);
        h1s[i]->Draw();
        double ty = ty_start;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.05;
        }
        c->SaveAs((outdir_C+"/"+h1_titles[i]+".C").c_str());
        c->SaveAs((outdir_png+"/"+h1_titles[i]+".png").c_str());
        delete c;
    }

    TFile * fout = new TFile((outdir_rootifles+"/GlobalHistos.root").c_str(), "RECREATE");
    fout->cd();
    h2_SUM_Et_vs_cent->Write();
    h2_sumMBDq_vs_cent->Write();
    h2_sumMBDq_vs_SUM_Et->Write();
    h1_SUM_Et->Write();
    h1_sumMBDq->Write();
    h1_cent->Write();
    h1_zvtx->Write();

    fout->Close();
    f->Close();

    std::cout <<"ProcessGlobal: End" << std::endl;
    return ;
    
}

void ProcessCaloWindowTree( const std::string & input_file, bool x_axis_cent )
{


    bool x_axis_sum_et = !x_axis_cent;
    std::string outdir = calo_window_plots;
    if (x_axis_cent) { outdir += "/xaxis_cent"; }
    if (x_axis_sum_et) { outdir += "/xaxis_sum_et"; }
    if( !gSystem->OpenDirectory(outdir.c_str()) ) {
        gSystem->mkdir(outdir.c_str(), true);
    }


    TFile * f = new TFile(input_file.c_str(), "READ");
    if( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl;  exit(1); }

    TTree * t = (TTree*)f->Get("T");
    
    // tree branches 
        float mbd_q_N = 0;
        float mbd_q_S = 0;
        t->SetBranchAddress("mbd_q_N", &mbd_q_N);
        t->SetBranchAddress("mbd_q_S", &mbd_q_S);


        int centrality = 0;
        t->SetBranchAddress("centrality", &centrality);

        float tower_frac_fired_TOWERINFO_CALIB_CEMC = 0;
        float tower_frac_dead_TOWERINFO_CALIB_CEMC = 0;
        float tower_avg_energy_TOWERINFO_CALIB_CEMC = 0;
        float tower_std_energy_TOWERINFO_CALIB_CEMC = 0;
        float tower_sum_energy_TOWERINFO_CALIB_CEMC = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC", &tower_frac_fired_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC", &tower_frac_dead_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC", &tower_avg_energy_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC", &tower_std_energy_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC", &tower_sum_energy_TOWERINFO_CALIB_CEMC);

        float tower_frac_fired_TOWERINFO_CALIB_HCALIN = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALIN = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALIN = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALIN = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALIN = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN", &tower_frac_fired_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN", &tower_frac_dead_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN", &tower_avg_energy_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN", &tower_std_energy_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN", &tower_sum_energy_TOWERINFO_CALIB_HCALIN);

        float tower_frac_fired_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALOUT = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT", &tower_frac_fired_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT", &tower_frac_dead_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT", &tower_avg_energy_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT", &tower_std_energy_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT", &tower_sum_energy_TOWERINFO_CALIB_HCALOUT);

        float tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER", &tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER", &tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER);

        float tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);

        float tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1", &tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1", &tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1);

        float tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1);
      
        unsigned int max_window_vector_size = 0;  
        float avg_energy_full[11];
        float std_energy_full[11];
        float avg_frac_energy_recemc_full[11];
        float std_frac_energy_recemc_full[11];
        float avg_frac_energy_hcalin_full[11];
        float std_frac_energy_hcalin_full[11];
        float avg_frac_energy_hcalout_full[11];
        float std_frac_energy_hcalout_full[11];
        float avg_energy_recemc[11];
        float std_energy_recemc[11];
        float avg_energy_hcalin[11];
        float std_energy_hcalin[11];
        float avg_energy_hcalout[11];
        float std_energy_hcalout[11];
        int num_windows_full[11];
        int num_windows_recemc[11];
        int num_windows_hcalin[11];
        int num_windows_hcalout[11];
        float avg_energy_cemc[11];
        float std_energy_cemc[11];
        int num_windows_cemc[11];

        t->SetBranchAddress("num_window_dims", &max_window_vector_size);
        t->SetBranchAddress("avg_energy_full", &avg_energy_full);
        t->SetBranchAddress("std_energy_full", &std_energy_full);
        t->SetBranchAddress("avg_frac_energy_recemc_full", &avg_frac_energy_recemc_full);
        t->SetBranchAddress("std_frac_energy_recemc_full", &std_frac_energy_recemc_full);
        t->SetBranchAddress("avg_frac_energy_hcalin_full", &avg_frac_energy_hcalin_full);
        t->SetBranchAddress("std_frac_energy_hcalin_full", &std_frac_energy_hcalin_full);
        t->SetBranchAddress("avg_frac_energy_hcalout_full", &avg_frac_energy_hcalout_full);
        t->SetBranchAddress("std_frac_energy_hcalout_full", &std_frac_energy_hcalout_full);
        t->SetBranchAddress("avg_energy_recemc", &avg_energy_recemc);
        t->SetBranchAddress("std_energy_recemc", &std_energy_recemc);
        t->SetBranchAddress("avg_energy_hcalin", &avg_energy_hcalin);
        t->SetBranchAddress("std_energy_hcalin", &std_energy_hcalin);
        t->SetBranchAddress("avg_energy_hcalout", &avg_energy_hcalout);
        t->SetBranchAddress("std_energy_hcalout", &std_energy_hcalout);
        t->SetBranchAddress("num_windows_full", &num_windows_full);
        t->SetBranchAddress("num_windows_recemc", &num_windows_recemc);
        t->SetBranchAddress("num_windows_hcalin", &num_windows_hcalin);
        t->SetBranchAddress("num_windows_hcalout", &num_windows_hcalout);

        t->SetBranchAddress("avg_energy_cemc", &avg_energy_cemc);
        t->SetBranchAddress("std_energy_cemc", &std_energy_cemc);
        t->SetBranchAddress("num_windows_cemc", &num_windows_cemc);

    int nentries = t->GetEntries();


    const int AVG_ET_BINS = 200;
    TH2F* h2_avg_et_full_window[k_window_array_size];
    TH2F* h2_sigma_et_full_window[k_window_array_size];
    TProfile * p2_avg_cent_full_window[k_window_array_size];
    TProfile* p2_sigma_cent_full_window[k_window_array_size];

    const int N_COURSE_AVG_ET_BINS = 8;
    float COURSE_AVG_ET_BINS[N_COURSE_AVG_ET_BINS+1];
    float MAX_COURSE_AVG_ET = 2200;
    if (x_axis_cent){ MAX_COURSE_AVG_ET = MAX_X_CENT; }
    for ( int i = 0; i < N_COURSE_AVG_ET_BINS+1; ++i ) { COURSE_AVG_ET_BINS[i] = i*MAX_COURSE_AVG_ET/N_COURSE_AVG_ET_BINS; }
    TProfile* p2_sigma_course_full_window[k_window_array_size];
    for ( unsigned int i = 0; i < k_window_array_size; ++i ) {
        
        float avgbins[AVG_ET_BINS+1];
        float MAX_AVG  = (float)k_calo_window_dims_hcal_geom[i].first*k_calo_window_dims_hcal_geom[i].second*2;
        for ( int j = 0; j < AVG_ET_BINS+1; ++j ) { avgbins[j] = j*MAX_AVG/AVG_ET_BINS; }

        float simgabins[AVG_ET_BINS+1];
        float MAX_STD  = (float)std::sqrt(k_calo_window_dims_hcal_geom[i].first*k_calo_window_dims_hcal_geom[i].second)*2;
        for ( int j = 0; j < AVG_ET_BINS+1; ++j ) { simgabins[j] = j*MAX_STD/AVG_ET_BINS; }
        if ( x_axis_sum_et ){
            h2_avg_et_full_window[i] = new TH2F(Form("h2_avg_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                                Form("h2_avg_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                                N_SUM_ET_BINS, SUM_ET_BINS, AVG_ET_BINS, avgbins);
            h2_avg_et_full_window[i]->GetXaxis()->SetTitle("Centrality [%]");
            h2_avg_et_full_window[i]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
            h2_avg_et_full_window[i]->GetYaxis()->SetTitle(Form("#LT E_{T}^{%dx%d} #GT [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));

            h2_sigma_et_full_window[i] = new TH2F(Form("h2_sigma_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            Form("h2_sigma_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            N_SUM_ET_BINS, SUM_ET_BINS, AVG_ET_BINS, simgabins);
            h2_sigma_et_full_window[i]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
            h2_sigma_et_full_window[i]->GetYaxis()->SetTitle(Form("#sigma^{%dx%d} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));

            p2_avg_cent_full_window[i] = new TProfile(Form("p2_avg_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            Form("p2_avg_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            N_SUM_ET_BINS, SUM_ET_BINS);
            p2_avg_cent_full_window[i]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
            p2_avg_cent_full_window[i]->GetYaxis()->SetTitle(Form("#bar{#LT E_{T}^{%dx%d} #GT} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));

            p2_sigma_cent_full_window[i] = new TProfile(Form("p2_sigma_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            Form("p2_sigma_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            N_SUM_ET_BINS, SUM_ET_BINS);
            p2_sigma_cent_full_window[i]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
            p2_sigma_cent_full_window[i]->GetYaxis()->SetTitle(Form("#bar{#sigma}^{%dx%d} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));


            p2_sigma_course_full_window[i] = new TProfile(Form("p2_sigma_course_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            Form("p2_sigma_course_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            N_COURSE_AVG_ET_BINS, COURSE_AVG_ET_BINS);
            p2_sigma_course_full_window[i]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
            p2_sigma_course_full_window[i]->GetYaxis()->SetTitle(Form("#bar{#sigma}^{%dx%d} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));
        } else 
        {
            h2_avg_et_full_window[i] = new TH2F(Form("h2_avg_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                                Form("h2_avg_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                                N_X_CENT_BINS, X_CENT_BINS, AVG_ET_BINS, avgbins);
            h2_avg_et_full_window[i]->GetXaxis()->SetTitle("Centrality [%]");
            h2_avg_et_full_window[i]->GetYaxis()->SetTitle(Form("#LT E_{T}^{%dx%d} #GT [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));

            h2_sigma_et_full_window[i] = new TH2F(Form("h2_sigma_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            Form("h2_sigma_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            N_X_CENT_BINS, X_CENT_BINS, AVG_ET_BINS, simgabins);
            h2_sigma_et_full_window[i]->GetXaxis()->SetTitle("Centrality [%]");
            h2_sigma_et_full_window[i]->GetYaxis()->SetTitle(Form("#sigma^{%dx%d} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));

            p2_avg_cent_full_window[i] = new TProfile(Form("p2_avg_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            Form("p2_avg_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            N_X_CENT_BINS, X_CENT_BINS);
            p2_avg_cent_full_window[i]->GetXaxis()->SetTitle("Centrality [%]");
            p2_avg_cent_full_window[i]->GetYaxis()->SetTitle(Form("#bar{#LT E_{T}^{%dx%d} #GT} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));

            p2_sigma_cent_full_window[i] = new TProfile(Form("p2_sigma_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            Form("p2_sigma_et_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            N_X_CENT_BINS, X_CENT_BINS);    
            p2_sigma_cent_full_window[i]->GetXaxis()->SetTitle("Centrality [%]");
            p2_sigma_cent_full_window[i]->GetYaxis()->SetTitle(Form("#bar{#sigma}^{%dx%d} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));

            p2_sigma_course_full_window[i] = new TProfile(Form("p2_sigma_course_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            Form("p2_sigma_course_full_window_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
                                            N_COURSE_AVG_ET_BINS, COURSE_AVG_ET_BINS);
            p2_sigma_course_full_window[i]->GetXaxis()->SetTitle("Centrality [%]");
            p2_sigma_course_full_window[i]->GetYaxis()->SetTitle(Form("#bar{#sigma}^{%dx%d} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));
        }

    }

    for ( int i = 0; i < nentries; ++i ) {
        t->GetEntry(i);
        float sum_et = tower_sum_energy_TOWERINFO_CALIB_CEMC + tower_sum_energy_TOWERINFO_CALIB_HCALIN + tower_sum_energy_TOWERINFO_CALIB_HCALOUT;
        sum_et = mbd_q_N + mbd_q_S;
        for ( unsigned int iwindow = 0; iwindow < k_window_array_size; ++iwindow ) {
            if ( num_windows_full[iwindow] == 0 ) { continue; }
            if ( x_axis_sum_et ){
                h2_avg_et_full_window[iwindow]->Fill(sum_et, avg_energy_full[iwindow]);
                h2_sigma_et_full_window[iwindow]->Fill(sum_et, std_energy_full[iwindow]);
                p2_avg_cent_full_window[iwindow]->Fill(sum_et, avg_energy_full[iwindow]);
                p2_sigma_cent_full_window[iwindow]->Fill(sum_et, std_energy_full[iwindow]);
                p2_sigma_course_full_window[iwindow]->Fill(sum_et, std_energy_full[iwindow]);
            } else {
                h2_avg_et_full_window[iwindow]->Fill(centrality, avg_energy_full[iwindow]);
                h2_sigma_et_full_window[iwindow]->Fill(centrality, std_energy_full[iwindow]);
                p2_avg_cent_full_window[iwindow]->Fill(centrality, avg_energy_full[iwindow]);
                p2_sigma_cent_full_window[iwindow]->Fill(centrality, std_energy_full[iwindow]);
                p2_sigma_course_full_window[iwindow]->Fill(centrality, std_energy_full[iwindow]);
            }
        }
    }


    TGraphErrors * g_sigma_vs_towerarea[N_COURSE_AVG_ET_BINS];
    std::vector<float> tower_areas {};
    for ( unsigned int i = 0; i < k_window_array_size; ++i ) {
        // if ( num_windows_full[i] == 0 ) { continue; }
        if ( k_calo_window_dims_hcal_geom[i].first > 15 ) { continue; }
        tower_areas.push_back(1.0*k_calo_window_dims_hcal_geom[i].first*k_calo_window_dims_hcal_geom[i].second);
    }
    for ( int j = 0; j < N_COURSE_AVG_ET_BINS; ++j ) {
        std::vector<float> avg_sigmas {};
        std::vector<float> avg_sigmas_err {};
        for ( unsigned int i = 0; i < k_window_array_size; ++i ) {
            if ( k_calo_window_dims_hcal_geom[i].first > 15 ) { continue; }
            avg_sigmas.push_back(p2_sigma_course_full_window[i]->GetBinContent(j+1));
            avg_sigmas_err.push_back(p2_sigma_course_full_window[i]->GetBinError(j+1));
        }
        g_sigma_vs_towerarea[j] = new TGraphErrors(tower_areas.size(), &tower_areas[0], &avg_sigmas[0], 0, &avg_sigmas_err[0]);
    }

    TCanvas * c;
    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    double tx=0.19;
    if ( x_axis_cent ){ tx = 0.45; }
    double ty_start=0.85;
    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
    std::vector<TH2F*> h2s = {h2_avg_et_full_window[0], h2_sigma_et_full_window[0]};
    for ( unsigned int i = 0; i < k_window_array_size; ++i ) {
        tags = {sPHENIX_Tag, DataType_Tag, Form("%dx%d Towers", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second)};
        h2s.clear();
        h2s.push_back(h2_avg_et_full_window[i]);
        h2s.push_back(h2_sigma_et_full_window[i]);
        double tx=0.19;
        if ( x_axis_cent ){ tx = 0.45; }
        double ty_start=0.85;
        for ( auto h2 : h2s ) {
            double ty = ty_start;
            c = new TCanvas("c", "c", 800, 600);
            gPad->SetLogz();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.1);
            gPad->SetBottomMargin(0.15);
            gPad->SetTopMargin(0.05);
            h2->GetXaxis()->SetNdivisions(505);
            h2->GetYaxis()->SetNdivisions(505);
            // h2->GetXaxis()->SetRangeUser(0, 2200);
            h2->GetYaxis()->SetRangeUser(0, h2->GetYaxis()->GetXmax());
            h2->Draw("colz");
            for ( auto tag : tags ) {
                tex->DrawLatex(tx, ty, tag.c_str());
                ty -= 0.05;
            }

            c->SaveAs((outdir+"/"+h2->GetTitle()+".png").c_str());
            delete c;
        }

        std::vector<TProfile*> p2s = {p2_avg_cent_full_window[i], p2_sigma_cent_full_window[i], p2_sigma_course_full_window[i]};
        for ( auto p2 : p2s ) {
            double ty = ty_start;
            c = new TCanvas("c", "c", 800, 600);
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.1);
            gPad->SetBottomMargin(0.15);
            gPad->SetTopMargin(0.05);
            p2->GetXaxis()->SetNdivisions(505);
            p2->GetYaxis()->SetNdivisions(505);
            // p2->GetXaxis()->SetRangeUser(0, 2200);

            p2->Draw();
            for ( auto tag : tags ) {
                tex->DrawLatex(tx, ty, tag.c_str());
                ty -= 0.05;
            }

            c->SaveAs((outdir+"/"+p2->GetTitle()+".png").c_str());
            delete c;
        }
        
    }

    c = new TCanvas("c", "c", 800, 600);
    tags = {sPHENIX_Tag, DataType_Tag};
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    TLegend * leg = new TLegend(0.18,0.5,0.35,0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetNColumns(2);
    tx = 0.45;

    double ty = 0.35;
    for ( unsigned int i = 0; i < k_window_array_size; ++i ) {
        if ( k_calo_window_dims_hcal_geom[i].first > 15 ) { continue; }
        // std::vector<TProfile*> p2s = {p2_avg_cent_full_window[i], p2_sigma_cent_full_window[i], p2_sigma_course_full_window[i]};
        TProfile * p2 = p2_sigma_course_full_window[i];
        // p2->GetXaxis()->SetNdivisions(505);
        // p2->GetYaxis()->SetNdivisions(505);
        // p2->GetXaxis()->SetRangeUser(0, 2200);
        p2->SetLineColor(COLORS[i]);
        p2->SetMarkerColor(COLORS[i]);
        p2->SetMarkerStyle(MARKERS[i]);
        p2->SetMarkerSize(1.5);
        p2->SetLineWidth(2);
        p2->GetYaxis()->SetRangeUser(0,15);
        p2->GetYaxis()->SetTitle("#LT#sigma#GT [GeV]");
        if ( x_axis_sum_et ) p2->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        else p2->GetXaxis()->SetTitle("Centrality [%]");
        // p2->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        c->cd();
        if ( i == 0 ) {
            p2->Draw("APL");
        } else {
            p2->Draw("SAME");
        }
        leg->AddEntry(p2, Form("%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), "lp");
    }
    leg->Draw("SAME");
    tags = {sPHENIX_Tag, DataType_Tag};
    for ( auto tag : tags ) {
        tex->DrawLatex(tx, ty, tag.c_str());
        ty -= 0.07;
    }
    c->SaveAs((outdir+"/p2_sigma_vs_SUM_Et.png").c_str());
    delete c;


    for ( int j = 0; j < N_COURSE_AVG_ET_BINS; ++j ) {

        TLegend * leg2 = new TLegend(0.18,0.8,0.35,0.9);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        TF1 * f1 = new TF1("f1", "[0]*(x^[1])", 0, 250);
        tags = {sPHENIX_Tag, DataType_Tag};
        if (x_axis_sum_et) { tags.push_back(Form("%0.0f < #Sigma E_{T}^{Raw} < %0.0f GeV", COURSE_AVG_ET_BINS[j], COURSE_AVG_ET_BINS[j+1])); }
        else { tags.push_back(Form("%0.0f - %0.0f %% Central", COURSE_AVG_ET_BINS[j], COURSE_AVG_ET_BINS[j+1])); }
        c = new TCanvas("c", "c", 800, 600);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        g_sigma_vs_towerarea[j]->GetXaxis()->SetTitle("A^{window} / A^{1 #times 1} ");
        g_sigma_vs_towerarea[j]->GetYaxis()->SetTitle("#bar{#sigma} [GeV]");
        g_sigma_vs_towerarea[j]->GetYaxis()->SetRangeUser(0, g_sigma_vs_towerarea[j]->GetHistogram()->GetMaximum()*1.1);
        
        g_sigma_vs_towerarea[j]->GetXaxis()->SetNdivisions(505);
        g_sigma_vs_towerarea[j]->Draw("AP");
        // get first point 
        float y0 = g_sigma_vs_towerarea[j]->GetY()[0];
        f1->FixParameter(0, y0);
        g_sigma_vs_towerarea[j]->Fit("f1", "Q", "", 0, 250);
        // get the fit
        TF1 * fit = g_sigma_vs_towerarea[j]->GetFunction("f1");
        fit->SetLineColor(kRed);
        fit->SetLineWidth(2);
        fit->SetLineStyle(2);
        double tx = 0.45;
        double ty = 0.45;
        // tags.push_back(Form("k = %0.2f #pm %0.2f", fit->GetParameter(1), fit->GetParError(1)));
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.07;
        }
        leg2->AddEntry(fit, Form("k = %0.2f #pm %0.2f", fit->GetParameter(1), fit->GetParError(1)), "l");
        leg2->Draw("SAME");
        c->SaveAs((outdir+"/g_sigma_vs_towerarea_"+std::to_string(j)+".png").c_str());
        delete c;

        TCanvas * c2 = new TCanvas("c2", "c2", 800, 600);
        TGraphErrors * g_sigma_vs_towerarea_ratio = (TGraphErrors*)g_sigma_vs_towerarea[j]->Clone("g_sigma_vs_towerarea_ratio");
        for (int i = 0; i < g_sigma_vs_towerarea_ratio->GetN(); ++i) {
            double x, y;
            g_sigma_vs_towerarea_ratio->GetPoint(i, x, y);
            double f0 = fit->Eval(x);
            g_sigma_vs_towerarea_ratio->SetPoint(i, x, y/f0);
            g_sigma_vs_towerarea_ratio->Draw("AP");
            g_sigma_vs_towerarea_ratio->GetYaxis()->SetTitle("#bar{#sigma} / Fit ");
            c2->SaveAs((outdir+"/g_sigma_vs_towerarea_ratio_"+std::to_string(j)+".png").c_str());
        }
    }

    TFile * fout = new TFile((outdir+"/CaloWindowHistos.root").c_str(), "RECREATE");
    fout->cd();
    for ( unsigned int i = 0; i < k_window_array_size; ++i ) {
        h2_avg_et_full_window[i]->Write();
        h2_sigma_et_full_window[i]->Write();
        p2_avg_cent_full_window[i]->Write();
        p2_sigma_cent_full_window[i]->Write();
        p2_sigma_course_full_window[i]->Write();
    }
    for ( int j = 0; j < N_COURSE_AVG_ET_BINS; ++j ) {
        g_sigma_vs_towerarea[j]->Write();
    }
    fout->Close();


    f->Close();

    return ;
    
}

void ProcessBackgroundTree(const std::string & input_file, bool x_axis_cent )
{


    bool x_axis_sum_et = !x_axis_cent;
    std::string outdir = bkgd_plots;
    if (x_axis_cent) { outdir += "/xaxis_cent"; }
    if (x_axis_sum_et) { outdir += "/xaxis_sum_et"; }
    if( !gSystem->OpenDirectory(outdir.c_str()) ) {
        gSystem->mkdir(outdir.c_str(), true);
    }


    TFile * f = new TFile(input_file.c_str(), "READ");
    if( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl;  exit(1); }

    TTree * t = (TTree*)f->Get("T");
    
    // tree branches 
        int centrality = 0;
        t->SetBranchAddress("centrality", &centrality);

        std::vector< float > * tower_background_energy_recemc = 0;
        std::vector< float > * tower_background_energy_hcalin = 0;
        std::vector< float > * tower_background_energy_hcalout = 0;
        t->SetBranchAddress("tower_background_energy_recemc", &tower_background_energy_recemc);
        t->SetBranchAddress("tower_background_energy_hcalin", &tower_background_energy_hcalin);
        t->SetBranchAddress("tower_background_energy_hcalout", &tower_background_energy_hcalout);


        float rho_val_TowerRho_AREA = 0;
        float std_rho_val_TowerRho_AREA = 0;
        t->SetBranchAddress("rho_val_TowerRho_AREA", &rho_val_TowerRho_AREA);
        t->SetBranchAddress("std_rho_val_TowerRho_AREA", &std_rho_val_TowerRho_AREA);
        
        float rho_val_TowerRho_MULT = 0;
        float std_rho_val_TowerRho_MULT = 0;
        t->SetBranchAddress("rho_val_TowerRho_MULT", &rho_val_TowerRho_MULT);
        t->SetBranchAddress("std_rho_val_TowerRho_MULT", &std_rho_val_TowerRho_MULT);

        float rho_val_TowerRho_AREA_CEMC = 0;
        float std_rho_val_TowerRho_AREA_CEMC = 0;
        t->SetBranchAddress("rho_val_TowerRho_AREA_CEMC", &rho_val_TowerRho_AREA_CEMC);
        t->SetBranchAddress("std_rho_val_TowerRho_AREA_CEMC", &std_rho_val_TowerRho_AREA_CEMC);

        float rho_val_TowerRho_MULT_CEMC = 0;
        float std_rho_val_TowerRho_MULT_CEMC = 0;
        t->SetBranchAddress("rho_val_TowerRho_MULT_CEMC", &rho_val_TowerRho_MULT_CEMC);
        t->SetBranchAddress("std_rho_val_TowerRho_MULT_CEMC", &std_rho_val_TowerRho_MULT_CEMC);

        float rho_val_TowerRho_AREA_HCALIN = 0;
        float std_rho_val_TowerRho_AREA_HCALIN = 0;
        t->SetBranchAddress("rho_val_TowerRho_AREA_HCALIN", &rho_val_TowerRho_AREA_HCALIN);
        t->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALIN", &std_rho_val_TowerRho_AREA_HCALIN);

        float rho_val_TowerRho_MULT_HCALIN = 0;
        float std_rho_val_TowerRho_MULT_HCALIN = 0;
        t->SetBranchAddress("rho_val_TowerRho_MULT_HCALIN", &rho_val_TowerRho_MULT_HCALIN);
        t->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALIN", &std_rho_val_TowerRho_MULT_HCALIN);

        float rho_val_TowerRho_AREA_HCALOUT = 0;
        float std_rho_val_TowerRho_AREA_HCALOUT = 0;
        t->SetBranchAddress("rho_val_TowerRho_AREA_HCALOUT", &rho_val_TowerRho_AREA_HCALOUT);
        t->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALOUT", &std_rho_val_TowerRho_AREA_HCALOUT);

        float rho_val_TowerRho_MULT_HCALOUT = 0;
        float std_rho_val_TowerRho_MULT_HCALOUT = 0;
        t->SetBranchAddress("rho_val_TowerRho_MULT_HCALOUT", &rho_val_TowerRho_MULT_HCALOUT);
        t->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALOUT", &std_rho_val_TowerRho_MULT_HCALOUT);

        float tower_frac_fired_TOWERINFO_CALIB_CEMC = 0;
        float tower_frac_dead_TOWERINFO_CALIB_CEMC = 0;
        float tower_avg_energy_TOWERINFO_CALIB_CEMC = 0;
        float tower_std_energy_TOWERINFO_CALIB_CEMC = 0;
        float tower_sum_energy_TOWERINFO_CALIB_CEMC = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC", &tower_frac_fired_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC", &tower_frac_dead_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC", &tower_avg_energy_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC", &tower_std_energy_TOWERINFO_CALIB_CEMC);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC", &tower_sum_energy_TOWERINFO_CALIB_CEMC);

        float tower_frac_fired_TOWERINFO_CALIB_HCALIN = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALIN = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALIN = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALIN = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALIN = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN", &tower_frac_fired_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN", &tower_frac_dead_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN", &tower_avg_energy_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN", &tower_std_energy_TOWERINFO_CALIB_HCALIN);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN", &tower_sum_energy_TOWERINFO_CALIB_HCALIN);

        float tower_frac_fired_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALOUT = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALOUT = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT", &tower_frac_fired_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT", &tower_frac_dead_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT", &tower_avg_energy_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT", &tower_std_energy_TOWERINFO_CALIB_HCALOUT);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT", &tower_sum_energy_TOWERINFO_CALIB_HCALOUT);

        float tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        float tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER", &tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER", &tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER);

        float tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        float tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);

        float tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1", &tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1", &tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1);

        float tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        float tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
        t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1);
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1);


        float random_cone_R_RandomCones_r04 = 0;
        float random_cone_eta_RandomCones_r04 = 0;
        float random_cone_phi_RandomCones_r04 = 0;
        float random_cone_energy_RandomCones_r04 = 0;
        float random_cone_energy_cemc_RandomCones_r04 = 0;
        float random_cone_energy_hcalin_RandomCones_r04 = 0;
        float random_cone_energy_hcalout_RandomCones_r04 = 0;
        int random_cone_num_towers_RandomCones_r04 = 0;
        int random_cone_num_towers_cemc_RandomCones_r04 = 0;
        int random_cone_num_towers_hcalin_RandomCones_r04 = 0;
        int random_cone_num_towers_hcalout_RandomCones_r04 = 0;
        int random_cone_num_masked_towers_RandomCones_r04 = 0;
        int random_cone_num_masked_towers_cemc_RandomCones_r04 = 0;
        int random_cone_num_masked_towers_hcalin_RandomCones_r04 = 0;
        int random_cone_num_masked_towers_hcalout_RandomCones_r04 = 0;
        t->SetBranchAddress("random_cone_R_RandomCones_r04", &random_cone_R_RandomCones_r04);
        t->SetBranchAddress("random_cone_eta_RandomCones_r04", &random_cone_eta_RandomCones_r04);
        t->SetBranchAddress("random_cone_phi_RandomCones_r04", &random_cone_phi_RandomCones_r04);
        t->SetBranchAddress("random_cone_energy_RandomCones_r04", &random_cone_energy_RandomCones_r04);
        t->SetBranchAddress("random_cone_energy_cemc_RandomCones_r04", &random_cone_energy_cemc_RandomCones_r04);
        t->SetBranchAddress("random_cone_energy_hcalin_RandomCones_r04", &random_cone_energy_hcalin_RandomCones_r04);
        t->SetBranchAddress("random_cone_energy_hcalout_RandomCones_r04", &random_cone_energy_hcalout_RandomCones_r04);
        t->SetBranchAddress("random_cone_num_towers_RandomCones_r04", &random_cone_num_towers_RandomCones_r04);
        t->SetBranchAddress("random_cone_num_towers_cemc_RandomCones_r04", &random_cone_num_towers_cemc_RandomCones_r04);
        t->SetBranchAddress("random_cone_num_towers_hcalin_RandomCones_r04", &random_cone_num_towers_hcalin_RandomCones_r04);
        t->SetBranchAddress("random_cone_num_towers_hcalout_RandomCones_r04", &random_cone_num_towers_hcalout_RandomCones_r04);
        t->SetBranchAddress("random_cone_num_masked_towers_RandomCones_r04", &random_cone_num_masked_towers_RandomCones_r04);
        t->SetBranchAddress("random_cone_num_masked_towers_cemc_RandomCones_r04", &random_cone_num_masked_towers_cemc_RandomCones_r04);
        t->SetBranchAddress("random_cone_num_masked_towers_hcalin_RandomCones_r04", &random_cone_num_masked_towers_hcalin_RandomCones_r04);
        t->SetBranchAddress("random_cone_num_masked_towers_hcalout_RandomCones_r04", &random_cone_num_masked_towers_hcalout_RandomCones_r04);

        float random_cone_R_RandomCones_r04_Sub1 = 0;
        float random_cone_eta_RandomCones_r04_Sub1 = 0;
        float random_cone_phi_RandomCones_r04_Sub1 = 0;
        float random_cone_energy_RandomCones_r04_Sub1 = 0;
        float random_cone_energy_cemc_RandomCones_r04_Sub1 = 0;
        float random_cone_energy_hcalin_RandomCones_r04_Sub1 = 0;
        float random_cone_energy_hcalout_RandomCones_r04_Sub1 = 0;
        int random_cone_num_towers_RandomCones_r04_Sub1 = 0;
        int random_cone_num_towers_cemc_RandomCones_r04_Sub1 = 0;
        int random_cone_num_towers_hcalin_RandomCones_r04_Sub1 = 0;
        int random_cone_num_towers_hcalout_RandomCones_r04_Sub1 = 0;
        int random_cone_num_masked_towers_RandomCones_r04_Sub1 = 0;
        int random_cone_num_masked_towers_cemc_RandomCones_r04_Sub1 = 0;
        int random_cone_num_masked_towers_hcalin_RandomCones_r04_Sub1 = 0;
        int random_cone_num_masked_towers_hcalout_RandomCones_r04_Sub1 = 0;
        t->SetBranchAddress("random_cone_R_RandomCones_r04_Sub1", &random_cone_R_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_eta_RandomCones_r04_Sub1", &random_cone_eta_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_phi_RandomCones_r04_Sub1", &random_cone_phi_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_energy_RandomCones_r04_Sub1", &random_cone_energy_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_energy_cemc_RandomCones_r04_Sub1", &random_cone_energy_cemc_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_energy_hcalin_RandomCones_r04_Sub1", &random_cone_energy_hcalin_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_energy_hcalout_RandomCones_r04_Sub1", &random_cone_energy_hcalout_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_num_towers_RandomCones_r04_Sub1", &random_cone_num_towers_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_num_towers_cemc_RandomCones_r04_Sub1", &random_cone_num_towers_cemc_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_num_towers_hcalin_RandomCones_r04_Sub1", &random_cone_num_towers_hcalin_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_num_towers_hcalout_RandomCones_r04_Sub1", &random_cone_num_towers_hcalout_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_num_masked_towers_RandomCones_r04_Sub1", &random_cone_num_masked_towers_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_num_masked_towers_cemc_RandomCones_r04_Sub1", &random_cone_num_masked_towers_cemc_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_num_masked_towers_hcalin_RandomCones_r04_Sub1", &random_cone_num_masked_towers_hcalin_RandomCones_r04_Sub1);
        t->SetBranchAddress("random_cone_num_masked_towers_hcalout_RandomCones_r04_Sub1", &random_cone_num_masked_towers_hcalout_RandomCones_r04_Sub1);

    int nentries = t->GetEntries();

    TH2F * h2_rho_area_vs_x;
    TH2F * h2_rho_mult_vs_x;
    TH2F * h2_towerbackground_vs_x;

    TH2F * h2_rho_area_times_A_vs_x;
    TH2F * h2_rho_mult_times_N_vs_x;
    TH2F * h2_towerbackground_times_avgT_vs_x;

    TH2F * h2_avgNcone_comp_vs_x;
    TH2F * h2_avgNcone_comp_sub1_vs_x[3]; // 3 calo layers
    if ( x_axis_sum_et ) {
        h2_rho_area_vs_x = new TH2F("h2_rho_area_vs_x", "h2_rho_area_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_RHO_A_BINS, RHO_A_BINS);
        h2_rho_mult_vs_x = new TH2F("h2_rho_mult_vs_x", "h2_rho_mult_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_RHO_M_BINS, RHO_M_BINS);
        h2_towerbackground_vs_x = new TH2F("h2_towerbackground_vs_x", "h2_towerbackground_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_BKGD_BINS, BKGD_BINS);

        h2_rho_area_times_A_vs_x = new TH2F("h2_rho_area_times_A_vs_x", "h2_rho_area_times_A_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_RHO_A_BINS, RHO_A_BINS);
        h2_rho_mult_times_N_vs_x = new TH2F("h2_rho_mult_times_N_vs_x", "h2_rho_mult_times_N_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_RHO_A_BINS, RHO_A_BINS);
        h2_towerbackground_times_avgT_vs_x = new TH2F("h2_towerbackground_times_avgT_vs_x", "h2_towerbackground_times_avgT_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_RHO_A_BINS, RHO_A_BINS);

        h2_avgNcone_comp_vs_x = new TH2F("h2_avgNcone_comp_vs_x", "h2_avgNcone_comp_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);

        for ( int i = 0; i < 3; ++i ) {
            h2_avgNcone_comp_sub1_vs_x[i] = new TH2F(Form("h2_avgNcone_comp_sub1_vs_x_layer%d", i), Form("h2_avgNcone_comp_sub1_vs_x_layer%d", i), 
                N_SUM_ET_BINS, SUM_ET_BINS, N_CONECOMP_SUB1_BINS, CONECOMP_SUB1_BINS);
            h2_avgNcone_comp_sub1_vs_x[i]->GetYaxis()->SetTitle(Form("N_{cone,sub1}^{avg}^{%d}", i));
            h2_avgNcone_comp_sub1_vs_x[i]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        }
        std::vector<TH2F*> h2s = {h2_rho_area_vs_x, h2_rho_mult_vs_x, h2_rho_area_times_A_vs_x, h2_rho_mult_times_N_vs_x, 
            h2_towerbackground_vs_x, h2_towerbackground_times_avgT_vs_x, h2_avgNcone_comp_vs_x};
        for ( auto h2 : h2s ) {
            h2->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        }
    } else {
        h2_rho_area_vs_x = new TH2F("h2_rho_area_vs_x", "h2_rho_area_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_A_BINS, RHO_A_BINS);
        h2_rho_mult_vs_x = new TH2F("h2_rho_mult_vs_x", "h2_rho_mult_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_M_BINS, RHO_M_BINS);
        h2_towerbackground_vs_x = new TH2F("h2_towerbackground_vs_x", "h2_towerbackground_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_BKGD_BINS, BKGD_BINS);

        h2_rho_area_times_A_vs_x = new TH2F("h2_rho_area_times_A_vs_x", "h2_rho_area_times_A_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_A_BINS, RHO_A_BINS);
        h2_rho_mult_times_N_vs_x = new TH2F("h2_rho_mult_times_N_vs_x", "h2_rho_mult_times_N_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_A_BINS, RHO_A_BINS);     
        h2_towerbackground_times_avgT_vs_x = new TH2F("h2_towerbackground_time_avgT_vs_x", "h2_towerbackground_time_avgT_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_A_BINS, RHO_A_BINS);     
       
        h2_avgNcone_comp_vs_x = new TH2F("h2_avgNcone_comp_vs_x", "h2_avgNcone_comp_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
        for ( int i = 0; i < 3; ++i ) {
            h2_avgNcone_comp_sub1_vs_x[i] = new TH2F(Form("h2_avgNcone_comp_sub1_vs_x_layer%d", i), Form("h2_avgNcone_comp_sub1_vs_x_layer%d", i), 
                N_X_CENT_BINS, X_CENT_BINS, N_CONECOMP_SUB1_BINS, CONECOMP_SUB1_BINS);
            h2_avgNcone_comp_sub1_vs_x[i]->GetYaxis()->SetTitle(Form("N_{cone,sub1}^{avg}^{%d}", i));
            h2_avgNcone_comp_sub1_vs_x[i]->GetXaxis()->SetTitle("Centrality [%]");
        }
        std::vector<TH2F*> h2s = {h2_rho_area_vs_x, h2_rho_mult_vs_x, h2_rho_area_times_A_vs_x, h2_rho_mult_times_N_vs_x, h2_towerbackground_vs_x, h2_towerbackground_times_avgT_vs_x, h2_avgNcone_comp_vs_x};
        for ( auto h2 : h2s ) {
            h2->GetXaxis()->SetTitle("Centrality [%]");
        }
    }

    h2_rho_area_vs_x->GetYaxis()->SetTitle("#rho_{A} [GeV]");
    h2_rho_mult_vs_x->GetYaxis()->SetTitle("#rho_{M} [GeV]");
    h2_rho_area_times_A_vs_x->GetYaxis()->SetTitle("#rho_{A} #times A [GeV]");
    h2_rho_mult_times_N_vs_x->GetYaxis()->SetTitle("#rho_{M} #times N [GeV]");
    h2_towerbackground_vs_x->GetYaxis()->SetTitle("E_{T}^{Bkgd}(#eta, #phi) [GeV]");
    h2_towerbackground_times_avgT_vs_x->GetYaxis()->SetTitle("E_{T}^{Bkgd}(#eta, #phi) d#eta d#phi [GeV]");
    h2_avgNcone_comp_vs_x->GetYaxis()->SetTitle("N_{cone}^{avg}");



    const int N_COURSE_X_BINS = 8;
    float COURSE_X_BINS[N_COURSE_X_BINS+1];
    float MAX_COURSE_X = 2000;
    if (x_axis_cent){ MAX_COURSE_X = MAX_X_CENT; }
    for ( int i = 0; i < N_COURSE_X_BINS+1; ++i ) { COURSE_X_BINS[i] = i*MAX_COURSE_X/N_COURSE_X_BINS; }
    TH2F * h2_rho_area_times_A_vs_course_x = new TH2F("h2_rho_area_times_A_vs_course_x", "h2_rho_area_times_A_vs_course_x", 
                    N_COURSE_X_BINS, COURSE_X_BINS, N_RHO_A_BINS, RHO_A_BINS);
    TH2F * h2_rho_mult_times_N_vs_course_x = new TH2F("h2_rho_mult_times_N_vs_course_x", "h2_rho_mult_times_N_vs_course_x",
                    N_COURSE_X_BINS, COURSE_X_BINS, N_RHO_A_BINS, RHO_A_BINS);
    TH2F * h2_towerbackground_times_avgT_vs_course_x = new TH2F("h2_towerbackground_times_avgT_vs_course_x", "h2_towerbackground_times_avgT_vs_course_x",
                    N_COURSE_X_BINS, COURSE_X_BINS, N_RHO_A_BINS, RHO_A_BINS);
  
 
    if (x_axis_sum_et) {
        h2_rho_area_times_A_vs_course_x->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        h2_rho_mult_times_N_vs_course_x->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        h2_towerbackground_times_avgT_vs_course_x->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
    } else {
        h2_rho_area_times_A_vs_course_x->GetXaxis()->SetTitle("Centrality [%]");
        h2_rho_mult_times_N_vs_course_x->GetXaxis()->SetTitle("Centrality [%]");
        h2_towerbackground_times_avgT_vs_course_x->GetXaxis()->SetTitle("Centrality [%]");
    }
    h2_rho_area_times_A_vs_course_x->GetYaxis()->SetTitle("#rho_{A} #times A [GeV]");
    h2_rho_mult_times_N_vs_course_x->GetYaxis()->SetTitle("#rho_{M} #times N [GeV]");
    h2_towerbackground_times_avgT_vs_course_x->GetYaxis()->SetTitle("E_{T}^{Bkgd}(#eta, #phi) d#eta d#phi [GeV]");

    
    for ( int i = 0; i < nentries; ++i ) {
        t->GetEntry(i);
        float sum_et = tower_sum_energy_TOWERINFO_CALIB_CEMC + tower_sum_energy_TOWERINFO_CALIB_HCALIN + tower_sum_energy_TOWERINFO_CALIB_HCALOUT;
        float xaxis_var = sum_et;
        if ( x_axis_cent ) { xaxis_var = 1.0*centrality; }
        h2_avgNcone_comp_vs_x->Fill(xaxis_var, random_cone_num_towers_RandomCones_r04);
        h2_avgNcone_comp_sub1_vs_x[0]->Fill(xaxis_var, random_cone_num_towers_cemc_RandomCones_r04_Sub1);
        h2_avgNcone_comp_sub1_vs_x[1]->Fill(xaxis_var, random_cone_num_towers_hcalin_RandomCones_r04_Sub1);
        h2_avgNcone_comp_sub1_vs_x[2]->Fill(xaxis_var, random_cone_num_towers_hcalout_RandomCones_r04_Sub1);
    }
    TProfile * p2_avgNcone_comp_vs_x = h2_avgNcone_comp_vs_x->ProfileX(Form("p2_avgNcone_comp_vs_x"), 1, -1, "s");
    TProfile * p2_avgNcone_comp_sub1_vs_x[3];
    for ( int i = 0; i < 3; ++i ) {
        p2_avgNcone_comp_sub1_vs_x[i] = h2_avgNcone_comp_sub1_vs_x[i]->ProfileX(Form("p2_avgNcone_comp_sub1_vs_x_layer%d", i), 1, -1, "s");
    }

    for ( int i = 0; i < nentries; ++i ) {
        t->GetEntry(i);
        float sum_et = tower_sum_energy_TOWERINFO_CALIB_CEMC + tower_sum_energy_TOWERINFO_CALIB_HCALIN + tower_sum_energy_TOWERINFO_CALIB_HCALOUT;
        float xaxis_var = sum_et;
        if ( x_axis_cent ) { xaxis_var = 1.0*centrality; }

        int xbin = h2_avgNcone_comp_vs_x->GetXaxis()->FindBin(xaxis_var);
        if ( xbin < 0 ) { continue; }

        h2_rho_area_vs_x->Fill(xaxis_var, rho_val_TowerRho_AREA);
        h2_rho_mult_vs_x->Fill(xaxis_var, rho_val_TowerRho_MULT);

        
        h2_rho_area_times_A_vs_x->Fill(xaxis_var, rho_val_TowerRho_AREA*AREA_CONE);
        
        h2_rho_mult_times_N_vs_x->Fill(xaxis_var, rho_val_TowerRho_MULT*p2_avgNcone_comp_vs_x->GetBinContent(xbin));
        float background_avg =0;
        float background_avg_layer[3] = {0};
        for ( int j = 0; j < tower_background_energy_recemc->size(); ++j ) {
            float bkgd = tower_background_energy_recemc->at(j) + tower_background_energy_hcalin->at(j) + tower_background_energy_hcalout->at(j);
            if ( std::isnan(bkgd) || std::isinf(bkgd) ) { continue; }
            background_avg_layer[0] += (tower_background_energy_recemc->at(j)/64.0);
            background_avg_layer[1] += (tower_background_energy_hcalin->at(j)/64.0);
            background_avg_layer[2] += (tower_background_energy_hcalout->at(j)/64.0);
            background_avg += bkgd;
        }
        background_avg_layer[0] /= 24.0;
        background_avg_layer[1] /= 24.0;
        background_avg_layer[2] /= 24.0;
        h2_towerbackground_vs_x->Fill(xaxis_var, 50*background_avg/(24.0*64.0));
        float avg_ncomps[3] = {0};
        float average_bkdg = 0;
        for ( int j = 0; j < 3; ++j ) {
            int bin = p2_avgNcone_comp_sub1_vs_x[j]->GetXaxis()->FindBin(xaxis_var);
            avg_ncomps[j] = p2_avgNcone_comp_sub1_vs_x[j]->GetBinContent(bin);
            average_bkdg+=background_avg_layer[j]*avg_ncomps[j];
        }
        h2_towerbackground_times_avgT_vs_x->Fill(xaxis_var, average_bkdg);

        h2_rho_area_times_A_vs_course_x->Fill(xaxis_var, rho_val_TowerRho_AREA*AREA_CONE);
        h2_rho_mult_times_N_vs_course_x->Fill(xaxis_var, rho_val_TowerRho_MULT*p2_avgNcone_comp_vs_x->GetBinContent(xbin));
        h2_towerbackground_times_avgT_vs_course_x->Fill(xaxis_var, average_bkdg);
    }


   
    TCanvas * c;
    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kRainBow);

    double tx=0.19;
    if ( x_axis_cent ){ tx = 0.45; }
    double ty_start=0.85;
    std::vector<TH2F*> h2s = {h2_rho_area_vs_x, h2_rho_mult_vs_x, h2_rho_area_times_A_vs_x, h2_rho_mult_times_N_vs_x, 
                             h2_towerbackground_vs_x, h2_towerbackground_times_avgT_vs_x, h2_avgNcone_comp_vs_x};
    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
    for ( auto h2 : h2s ) {
        double ty = ty_start;
        c = new TCanvas("c", "c", 800, 600);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        h2->GetXaxis()->SetNdivisions(505);
        h2->GetYaxis()->SetNdivisions(505);

        h2->Draw("colz");
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.05;
        }

        c->SaveAs((outdir+"/"+h2->GetTitle()+".png").c_str());
        delete c;
    }
    for ( int i = 0; i < 3; ++i ) {
        double ty = ty_start;
        c = new TCanvas("c", "c", 800, 600);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        h2_avgNcone_comp_sub1_vs_x[i]->GetXaxis()->SetNdivisions(505);
        h2_avgNcone_comp_sub1_vs_x[i]->GetYaxis()->SetNdivisions(505);
        tags = {sPHENIX_Tag, DataType_Tag};
        if ( i == 0 ) { tags.push_back("CEMC"); }
        else if ( i == 1 ) { tags.push_back("HCALIN"); }
        else if ( i == 2 ) { tags.push_back("HCALOUT"); }
        h2_avgNcone_comp_sub1_vs_x[i]->Draw("colz");
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.05;
        }

        c->SaveAs((outdir+"/"+h2_avgNcone_comp_sub1_vs_x[i]->GetTitle()+".png").c_str());
        delete c;
    }

    TLegend * leg = new TLegend(0.18,0.7,0.44,0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetNColumns(2);

    std::vector<TH2F*> h2s2 = {h2_rho_area_times_A_vs_course_x, h2_rho_mult_times_N_vs_course_x, h2_towerbackground_times_avgT_vs_course_x};
    double miny = 1e-4;

    for ( auto h2 : h2s2 ) {
        double ty = ty_start;
        tags = {sPHENIX_Tag, DataType_Tag};
        c = new TCanvas("c", "c", 800, 600);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        float maxx= 0;
       for ( int j = 0; j < N_COURSE_X_BINS-1; ++j ) {
            h2->GetXaxis()->SetRange(j+1, j+1);
            TH1D * h1 = h2->ProjectionY(Form("h1_%s_%d", h2->GetTitle(), j), j+1, j+1);
            h1->SetLineColor(COLORS[j]);
            h1->SetMarkerColor(COLORS[j]);
            h1->SetMarkerStyle(MARKERS[j]);
            h1->Scale(1.0/h1->Integral());

            std::string leg_title = Form("%0.0f < #Sigma E_{T}^{Raw} < %0.0f GeV", COURSE_X_BINS[j], COURSE_X_BINS[j+1]);
            if ( x_axis_cent ) { leg_title = Form("%0.0f-%0.0f %%", COURSE_X_BINS[j], COURSE_X_BINS[j+1]); }
            leg->AddEntry(h1, leg_title.c_str(), "lp");
            h1->GetYaxis()->SetRangeUser(miny, 1e1);
            int lastbin_above_threshold = 0;
            lastbin_above_threshold = h1->FindLastBinAbove(miny);
            h1->GetXaxis()->SetRangeUser(0, 1.2*h1->GetBinCenter(lastbin_above_threshold));
            h1->GetXaxis()->SetTitle(h2->GetYaxis()->GetTitle());
            h1->GetYaxis()->SetTitle("1/N dN/dE_{T}^{Bkgd} [GeV^{-1}]");
            
            
            if ( j == 0 ) {
                h1->Draw("P");
            } else {
                h1->Draw("SAME");
            }
       }

       leg->Draw("SAME");
       float txx = 0.5;
         float tyy = 0.85;
       for ( auto tag : tags ) {
            tex->DrawLatex(txx, tyy, tag.c_str());
            tyy -= 0.07;
        }

        c->SaveAs((outdir+"/"+h2->GetTitle()+".png").c_str());
        delete c;
        leg->Clear();
    }
    
    if ( x_axis_cent ) {
        h2s = {h2_rho_area_vs_x, h2_rho_mult_vs_x, h2_towerbackground_vs_x};
        for ( unsigned int ibin = 0; ibin  < N_X_CENT_BINS; ibin++){
            float maxx= 0;
            double miny = 1e-4;
            for ( auto h2 : h2s ) {
                tx = 0.19;
                double ty = ty_start;
                tags = {sPHENIX_Tag, DataType_Tag};
                c = new TCanvas("c", "c", 800, 600);
                gPad->SetLogy();
                gPad->SetLeftMargin(0.15);
                gPad->SetRightMargin(0.1);
                gPad->SetBottomMargin(0.15);
                gPad->SetTopMargin(0.05);
                h2->GetXaxis()->SetRange(ibin+1, ibin+1);
                TH1D * h1 = h2->ProjectionY(Form("h1_%s_%d", h2->GetTitle(), ibin), ibin+1, ibin+1);
                h1->Scale(1.0/h1->Integral());
                std::string leg_title = Form("%0.0f < #Sigma E_{T}^{Raw} < %0.0f GeV", COURSE_X_BINS[ibin], COURSE_X_BINS[ibin+1]);
                if ( x_axis_cent ) { leg_title = Form("%0.0f-%0.0f %%", COURSE_X_BINS[ibin], COURSE_X_BINS[ibin+1]); }
                tags.push_back(leg_title);
                h1->GetYaxis()->SetRangeUser(miny, 1e1);
                int lastbin_above_threshold = 0;
                lastbin_above_threshold = h1->FindLastBinAbove(miny);
                h1->GetXaxis()->SetRangeUser(0, h1->GetBinCenter(lastbin_above_threshold));
                h1->GetXaxis()->SetTitle(h2->GetYaxis()->GetTitle());
                h1->GetYaxis()->SetTitle("Probability Density [a.u.]");
                h1->Draw("P");
                for ( auto tag : tags ) {
                    tex->DrawLatex(tx, ty, tag.c_str());
                    ty -= 0.05;
                }
                c->SaveAs((outdir+"/"+h2->GetTitle()+"_cent"+std::to_string(ibin)+".png").c_str());
                delete c;
            }
        }
    }

    TFile * fout = new TFile((outdir+"/BackgroundEstimates.root").c_str(), "RECREATE");
    fout->cd();
    h2_rho_area_vs_x->Write();
    h2_rho_mult_vs_x->Write();
    h2_rho_area_times_A_vs_x->Write();
    h2_rho_mult_times_N_vs_x->Write();
    h2_towerbackground_vs_x->Write();
    h2_towerbackground_times_avgT_vs_x->Write();
    h2_avgNcone_comp_vs_x->Write();
    for ( int i = 0; i < 3; ++i ) {
        h2_avgNcone_comp_sub1_vs_x[i]->Write();
    }
    h2_rho_area_times_A_vs_course_x->Write();
    h2_rho_mult_times_N_vs_course_x->Write();
    h2_towerbackground_times_avgT_vs_course_x->Write();
    p2_avgNcone_comp_vs_x->Write();
    for ( int i = 0; i < 3; ++i ) {
        p2_avgNcone_comp_sub1_vs_x[i]->Write();
    }

    fout->Close();


    f->Close();

    return ;
    
}

// void ProcessRandomConeTree(const std::string & input_file, bool x_axis_cent )
// {



//     bool x_axis_sum_et = !x_axis_cent;
//     std::string outdir = random_cone_plots;
//     if (x_axis_cent) { outdir += "/xaxis_cent"; }
//     if (x_axis_sum_et) { outdir += "/xaxis_sum_et"; }
//     if( !gSystem->OpenDirectory(outdir.c_str()) ) {
//         gSystem->mkdir(outdir.c_str(), true);
//     }

//     TFile * f = new TFile(input_file.c_str(), "READ");
//     if( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl;  exit(1); }

//     TTree * t = (TTree*)f->Get("T");
    
//     // tree branches 
//         float mbd_q_N = 0;
//         float mbd_q_S = 0;
//         t->SetBranchAddress("mbd_q_N", &mbd_q_N);
//         t->SetBranchAddress("mbd_q_S", &mbd_q_S);
        
//         int centrality = 0;
//         t->SetBranchAddress("centrality", &centrality);

//         float rho_val_TowerRho_AREA = 0;
//         float std_rho_val_TowerRho_AREA = 0;
//         t->SetBranchAddress("rho_val_TowerRho_AREA", &rho_val_TowerRho_AREA);
//         t->SetBranchAddress("std_rho_val_TowerRho_AREA", &std_rho_val_TowerRho_AREA);
        
//         float rho_val_TowerRho_MULT = 0;
//         float std_rho_val_TowerRho_MULT = 0;
//         t->SetBranchAddress("rho_val_TowerRho_MULT", &rho_val_TowerRho_MULT);
//         t->SetBranchAddress("std_rho_val_TowerRho_MULT", &std_rho_val_TowerRho_MULT);

//         float rho_val_TowerRho_AREA_CEMC = 0;
//         float std_rho_val_TowerRho_AREA_CEMC = 0;
//         t->SetBranchAddress("rho_val_TowerRho_AREA_CEMC", &rho_val_TowerRho_AREA_CEMC);
//         t->SetBranchAddress("std_rho_val_TowerRho_AREA_CEMC", &std_rho_val_TowerRho_AREA_CEMC);

//         float rho_val_TowerRho_MULT_CEMC = 0;
//         float std_rho_val_TowerRho_MULT_CEMC = 0;
//         t->SetBranchAddress("rho_val_TowerRho_MULT_CEMC", &rho_val_TowerRho_MULT_CEMC);
//         t->SetBranchAddress("std_rho_val_TowerRho_MULT_CEMC", &std_rho_val_TowerRho_MULT_CEMC);

//         float rho_val_TowerRho_AREA_HCALIN = 0;
//         float std_rho_val_TowerRho_AREA_HCALIN = 0;
//         t->SetBranchAddress("rho_val_TowerRho_AREA_HCALIN", &rho_val_TowerRho_AREA_HCALIN);
//         t->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALIN", &std_rho_val_TowerRho_AREA_HCALIN);

//         float rho_val_TowerRho_MULT_HCALIN = 0;
//         float std_rho_val_TowerRho_MULT_HCALIN = 0;
//         t->SetBranchAddress("rho_val_TowerRho_MULT_HCALIN", &rho_val_TowerRho_MULT_HCALIN);
//         t->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALIN", &std_rho_val_TowerRho_MULT_HCALIN);

//         float rho_val_TowerRho_AREA_HCALOUT = 0;
//         float std_rho_val_TowerRho_AREA_HCALOUT = 0;
//         t->SetBranchAddress("rho_val_TowerRho_AREA_HCALOUT", &rho_val_TowerRho_AREA_HCALOUT);
//         t->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALOUT", &std_rho_val_TowerRho_AREA_HCALOUT);

//         float rho_val_TowerRho_MULT_HCALOUT = 0;
//         float std_rho_val_TowerRho_MULT_HCALOUT = 0;
//         t->SetBranchAddress("rho_val_TowerRho_MULT_HCALOUT", &rho_val_TowerRho_MULT_HCALOUT);
//         t->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALOUT", &std_rho_val_TowerRho_MULT_HCALOUT);

//         float tower_frac_fired_TOWERINFO_CALIB_CEMC = 0;
//         float tower_frac_dead_TOWERINFO_CALIB_CEMC = 0;
//         float tower_avg_energy_TOWERINFO_CALIB_CEMC = 0;
//         float tower_std_energy_TOWERINFO_CALIB_CEMC = 0;
//         float tower_sum_energy_TOWERINFO_CALIB_CEMC = 0;
//         t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC", &tower_frac_fired_TOWERINFO_CALIB_CEMC);
//         t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC", &tower_frac_dead_TOWERINFO_CALIB_CEMC);
//         t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC", &tower_avg_energy_TOWERINFO_CALIB_CEMC);
//         t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC", &tower_std_energy_TOWERINFO_CALIB_CEMC);
//         t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC", &tower_sum_energy_TOWERINFO_CALIB_CEMC);

//         float tower_frac_fired_TOWERINFO_CALIB_HCALIN = 0;
//         float tower_frac_dead_TOWERINFO_CALIB_HCALIN = 0;
//         float tower_avg_energy_TOWERINFO_CALIB_HCALIN = 0;
//         float tower_std_energy_TOWERINFO_CALIB_HCALIN = 0;
//         float tower_sum_energy_TOWERINFO_CALIB_HCALIN = 0;
//         t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN", &tower_frac_fired_TOWERINFO_CALIB_HCALIN);
//         t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN", &tower_frac_dead_TOWERINFO_CALIB_HCALIN);
//         t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN", &tower_avg_energy_TOWERINFO_CALIB_HCALIN);
//         t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN", &tower_std_energy_TOWERINFO_CALIB_HCALIN);
//         t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN", &tower_sum_energy_TOWERINFO_CALIB_HCALIN);

//         float tower_frac_fired_TOWERINFO_CALIB_HCALOUT = 0;
//         float tower_frac_dead_TOWERINFO_CALIB_HCALOUT = 0;
//         float tower_avg_energy_TOWERINFO_CALIB_HCALOUT = 0;
//         float tower_std_energy_TOWERINFO_CALIB_HCALOUT = 0;
//         float tower_sum_energy_TOWERINFO_CALIB_HCALOUT = 0;
//         t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT", &tower_frac_fired_TOWERINFO_CALIB_HCALOUT);
//         t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT", &tower_frac_dead_TOWERINFO_CALIB_HCALOUT);
//         t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT", &tower_avg_energy_TOWERINFO_CALIB_HCALOUT);
//         t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT", &tower_std_energy_TOWERINFO_CALIB_HCALOUT);
//         t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT", &tower_sum_energy_TOWERINFO_CALIB_HCALOUT);

//         float tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER = 0;
//         float tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER = 0;
//         float tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
//         float tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
//         float tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
//         t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER", &tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER);
//         t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER", &tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER);
//         t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER);
//         t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER);
//         t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER);

//         float tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
//         float tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
//         float tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
//         float tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
//         float tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
//         t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
//         t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
//         t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
//         t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
//         t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);

//         float tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
//         float tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
//         float tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
//         float tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
//         float tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
//         t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1", &tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1);
//         t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1", &tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1);
//         t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1);
//         t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1);
//         t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1);

//         float tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
//         float tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
//         float tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
//         float tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
//         float tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
//         t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1);
//         t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1);
//         t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1);
//         t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1);
//         t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1);


//         float random_cone_R_RandomCones_r04 = 0;
//         float random_cone_eta_RandomCones_r04 = 0;
//         float random_cone_phi_RandomCones_r04 = 0;
//         float random_cone_energy_RandomCones_r04 = 0;
//         float random_cone_energy_cemc_RandomCones_r04 = 0;
//         float random_cone_energy_hcalin_RandomCones_r04 = 0;
//         float random_cone_energy_hcalout_RandomCones_r04 = 0;
//         int random_cone_num_towers_RandomCones_r04 = 0;
//         int random_cone_num_towers_cemc_RandomCones_r04 = 0;
//         int random_cone_num_towers_hcalin_RandomCones_r04 = 0;
//         int random_cone_num_towers_hcalout_RandomCones_r04 = 0;
//         int random_cone_num_masked_towers_RandomCones_r04 = 0;
//         int random_cone_num_masked_towers_cemc_RandomCones_r04 = 0;
//         int random_cone_num_masked_towers_hcalin_RandomCones_r04 = 0;
//         int random_cone_num_masked_towers_hcalout_RandomCones_r04 = 0;
//         t->SetBranchAddress("random_cone_R_RandomCones_r04", &random_cone_R_RandomCones_r04);
//         t->SetBranchAddress("random_cone_eta_RandomCones_r04", &random_cone_eta_RandomCones_r04);
//         t->SetBranchAddress("random_cone_phi_RandomCones_r04", &random_cone_phi_RandomCones_r04);
//         t->SetBranchAddress("random_cone_energy_RandomCones_r04", &random_cone_energy_RandomCones_r04);
//         t->SetBranchAddress("random_cone_energy_cemc_RandomCones_r04", &random_cone_energy_cemc_RandomCones_r04);
//         t->SetBranchAddress("random_cone_energy_hcalin_RandomCones_r04", &random_cone_energy_hcalin_RandomCones_r04);
//         t->SetBranchAddress("random_cone_energy_hcalout_RandomCones_r04", &random_cone_energy_hcalout_RandomCones_r04);
//         t->SetBranchAddress("random_cone_num_towers_RandomCones_r04", &random_cone_num_towers_RandomCones_r04);
//         t->SetBranchAddress("random_cone_num_towers_cemc_RandomCones_r04", &random_cone_num_towers_cemc_RandomCones_r04);
//         t->SetBranchAddress("random_cone_num_towers_hcalin_RandomCones_r04", &random_cone_num_towers_hcalin_RandomCones_r04);
//         t->SetBranchAddress("random_cone_num_towers_hcalout_RandomCones_r04", &random_cone_num_towers_hcalout_RandomCones_r04);
//         t->SetBranchAddress("random_cone_num_masked_towers_RandomCones_r04", &random_cone_num_masked_towers_RandomCones_r04);
//         t->SetBranchAddress("random_cone_num_masked_towers_cemc_RandomCones_r04", &random_cone_num_masked_towers_cemc_RandomCones_r04);
//         t->SetBranchAddress("random_cone_num_masked_towers_hcalin_RandomCones_r04", &random_cone_num_masked_towers_hcalin_RandomCones_r04);
//         t->SetBranchAddress("random_cone_num_masked_towers_hcalout_RandomCones_r04", &random_cone_num_masked_towers_hcalout_RandomCones_r04);

//         float random_cone_R_RandomCones_r04_Sub1 = 0;
//         float random_cone_eta_RandomCones_r04_Sub1 = 0;
//         float random_cone_phi_RandomCones_r04_Sub1 = 0;
//         float random_cone_energy_RandomCones_r04_Sub1 = 0;
//         float random_cone_energy_cemc_RandomCones_r04_Sub1 = 0;
//         float random_cone_energy_hcalin_RandomCones_r04_Sub1 = 0;
//         float random_cone_energy_hcalout_RandomCones_r04_Sub1 = 0;
//         int random_cone_num_towers_RandomCones_r04_Sub1 = 0;
//         int random_cone_num_towers_cemc_RandomCones_r04_Sub1 = 0;
//         int random_cone_num_towers_hcalin_RandomCones_r04_Sub1 = 0;
//         int random_cone_num_towers_hcalout_RandomCones_r04_Sub1 = 0;
//         int random_cone_num_masked_towers_RandomCones_r04_Sub1 = 0;
//         int random_cone_num_masked_towers_cemc_RandomCones_r04_Sub1 = 0;
//         int random_cone_num_masked_towers_hcalin_RandomCones_r04_Sub1 = 0;
//         int random_cone_num_masked_towers_hcalout_RandomCones_r04_Sub1 = 0;
//         t->SetBranchAddress("random_cone_R_RandomCones_r04_Sub1", &random_cone_R_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_eta_RandomCones_r04_Sub1", &random_cone_eta_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_phi_RandomCones_r04_Sub1", &random_cone_phi_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_energy_RandomCones_r04_Sub1", &random_cone_energy_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_energy_cemc_RandomCones_r04_Sub1", &random_cone_energy_cemc_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_energy_hcalin_RandomCones_r04_Sub1", &random_cone_energy_hcalin_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_energy_hcalout_RandomCones_r04_Sub1", &random_cone_energy_hcalout_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_num_towers_RandomCones_r04_Sub1", &random_cone_num_towers_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_num_towers_cemc_RandomCones_r04_Sub1", &random_cone_num_towers_cemc_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_num_towers_hcalin_RandomCones_r04_Sub1", &random_cone_num_towers_hcalin_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_num_towers_hcalout_RandomCones_r04_Sub1", &random_cone_num_towers_hcalout_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_num_masked_towers_RandomCones_r04_Sub1", &random_cone_num_masked_towers_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_num_masked_towers_cemc_RandomCones_r04_Sub1", &random_cone_num_masked_towers_cemc_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_num_masked_towers_hcalin_RandomCones_r04_Sub1", &random_cone_num_masked_towers_hcalin_RandomCones_r04_Sub1);
//         t->SetBranchAddress("random_cone_num_masked_towers_hcalout_RandomCones_r04_Sub1", &random_cone_num_masked_towers_hcalout_RandomCones_r04_Sub1);



//     int nentries = t->GetEntries();

//     TH2F * h2_area_cone_res_vs_x;
//     TH2F * h2_mult_cone_res_vs_x;
//     TH2F * h2_sub1_cone_res_vs_x;
//     TH2F * h2_area_cone_ncomp_vs_x;
//     TH2F * h2_mult_cone_ncomp_vs_x;
//     TH2F * h2_sub1_cone_ncomp_vs_x;
    
//     if ( x_axis_sum_et ) {

//         h2_area_cone_res_vs_x = new TH2F("h2_area_cone_res_vs_x", "h2_area_cone_res_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
//         h2_mult_cone_res_vs_x = new TH2F("h2_mult_cone_res_vs_x", "h2_mult_cone_res_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
//         h2_sub1_cone_res_vs_x = new TH2F("h2_sub1_cone_res_vs_x", "h2_sub1_cone_res_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
//         h2_area_cone_ncomp_vs_x = new TH2F("h2_area_cone_ncomp_vs_x", "h2_area_cone_ncomp_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
//         h2_mult_cone_ncomp_vs_x = new TH2F("h2_mult_cone_ncomp_vs_x", "h2_mult_cone_ncomp_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
//         h2_sub1_cone_ncomp_vs_x = new TH2F("h2_sub1_cone_ncomp_vs_x", "h2_sub1_cone_ncomp_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONECOMP_SUB1_BINS, CONECOMP_SUB1_BINS);

//         std::vector<TH2F*> h2s = {h2_area_cone_res_vs_x, h2_mult_cone_res_vs_x, h2_sub1_cone_res_vs_x, h2_area_cone_ncomp_vs_x, h2_mult_cone_ncomp_vs_x, h2_sub1_cone_ncomp_vs_x};
        
//         for ( auto h2 : h2s ) {
//             h2->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
//         }
//     } else {
//         h2_area_cone_res_vs_x = new TH2F("h2_area_cone_res_vs_x", "h2_area_cone_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
//         h2_mult_cone_res_vs_x = new TH2F("h2_mult_cone_res_vs_x", "h2_mult_cone_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
//         h2_sub1_cone_res_vs_x = new TH2F("h2_sub1_cone_res_vs_x", "h2_sub1_cone_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
//         h2_area_cone_ncomp_vs_x = new TH2F("h2_area_cone_ncomp_vs_x", "h2_area_cone_ncomp_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
//         h2_mult_cone_ncomp_vs_x = new TH2F("h2_mult_cone_ncomp_vs_x", "h2_mult_cone_ncomp_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
//         h2_sub1_cone_ncomp_vs_x = new TH2F("h2_sub1_cone_ncomp_vs_x", "h2_sub1_cone_ncomp_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONECOMP_SUB1_BINS, CONECOMP_SUB1_BINS);
//         std::vector<TH2F*> h2s = {h2_area_cone_res_vs_x, h2_mult_cone_res_vs_x, h2_sub1_cone_res_vs_x, h2_area_cone_ncomp_vs_x, h2_mult_cone_ncomp_vs_x, h2_sub1_cone_ncomp_vs_x};
//         for ( auto h2 : h2s ) {
//             h2->GetXaxis()->SetTitle("Centrality [%]");
//         }
//     }

//     h2_area_cone_res_vs_x->GetYaxis()->SetTitle("#delta E_{T}^{Cone} [GeV]");
//     h2_mult_cone_res_vs_x->GetYaxis()->SetTitle("#delta E_{T}^{Cone} [GeV]");
//     h2_sub1_cone_res_vs_x->GetYaxis()->SetTitle("#delta E_{T}^{Cone} [GeV]");
//     h2_area_cone_ncomp_vs_x->GetYaxis()->SetTitle("N_{comp}^{Cone}");
//     h2_mult_cone_ncomp_vs_x->GetYaxis()->SetTitle("N_{comp}^{Cone}");
//     h2_sub1_cone_ncomp_vs_x->GetYaxis()->SetTitle("N_{comp}^{Cone}");

//     const int N_COURSE_X_BINS = 8;
//     float COURSE_X_BINS[N_COURSE_X_BINS+1];
//     float MAX_COURSE_X = 1800;
//     if (x_axis_cent){ MAX_COURSE_X = MAX_X_CENT; }
//     for ( int i = 0; i < N_COURSE_X_BINS+1; ++i ) { COURSE_X_BINS[i] = i*MAX_COURSE_X/N_COURSE_X_BINS; }
//     TH2F * h2_area_cone_res_vs_course_x = new TH2F("h2_area_cone_res_vs_course_x", "h2_area_cone_res_vs_course_x", 
//         N_COURSE_X_BINS, COURSE_X_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
//     TH2F * h2_mult_cone_res_vs_course_x = new TH2F("h2_mult_cone_res_vs_course_x", "h2_mult_cone_res_vs_course_x",
//         N_COURSE_X_BINS, COURSE_X_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
//     TH2F * h2_sub1_cone_res_vs_course_x = new TH2F("h2_sub1_cone_res_vs_course_x", "h2_sub1_cone_res_vs_course_x",
//         N_COURSE_X_BINS, COURSE_X_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
    
//     std::vector<TH2F*> h2s_course_x = {h2_area_cone_res_vs_course_x, h2_mult_cone_res_vs_course_x, h2_sub1_cone_res_vs_course_x};
//     for ( auto h2 : h2s_course_x ) {
//         if (x_axis_sum_et) {
//             h2->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
//         } else {
//             h2->GetXaxis()->SetTitle("Centrality [%]");
//         }
//         h2->GetYaxis()->SetTitle("#delta E_{T}^{Cone} [GeV]");
//     }

    
//     // N_SUM_ET_BINS, SUM_ET_BINS
//     float SUM_ET2[N_X_CENT_BINS];
//     float SUM_E[N_X_CENT_BINS];
//     float TOTAL_TOWERS_NOT_MASKED[N_X_CENT_BINS];
//     float TOTAL_TOWERS_FIRED[N_X_CENT_BINS];
//     float N_TOWERS_TOTAL[N_X_CENT_BINS];
//     float SUM_CONE_COMPS[N_X_CENT_BINS];
//     float N_CONES_THIS_BIN[N_X_CENT_BINS];

//     float SUM_ET2_SUB1[N_X_CENT_BINS];
//     float SUM_E_SUB1[N_X_CENT_BINS];
//     float TOTAL_TOWERS_NOT_MASKED_SUB1[N_X_CENT_BINS];
//     float TOTAL_TOWERS_FIRED_SUB1[N_X_CENT_BINS];
//     float N_TOWERS_TOTAL_SUB1[N_X_CENT_BINS];
//     float SUM_CONE_COMPS_SUB1[N_X_CENT_BINS];
//     float N_CONES_THIS_BIN_SUB1[N_X_CENT_BINS];

//     for ( int i = 0; i < N_X_CENT_BINS; ++i ) {

//         SUM_ET2[i] = 0;
//         SUM_E[i] = 0;
//         TOTAL_TOWERS_NOT_MASKED[i] = 0;
//         TOTAL_TOWERS_FIRED[i] = 0;
//         N_TOWERS_TOTAL[i] = 0;
//         SUM_CONE_COMPS[i] = 0;
//         N_CONES_THIS_BIN[i] = 0;

//         SUM_ET2_SUB1[i] = 0;
//         SUM_E_SUB1[i] = 0;
//         TOTAL_TOWERS_NOT_MASKED_SUB1[i] = 0;
//         TOTAL_TOWERS_FIRED_SUB1[i] = 0;
//         N_TOWERS_TOTAL_SUB1[i] = 0;
//         SUM_CONE_COMPS_SUB1[i] = 0;
//         N_CONES_THIS_BIN_SUB1[i] = 0;

//     }

//     float N_CEMC_TOWERS = 256*94;
//     float N_HCALIN_TOWERS = 24*64;
//     float N_HCALOUT_TOWERS = 24*64;
//     for ( int i = 0; i < nentries; ++i ) {
       
//         t->GetEntry(i);

//         float sum_et = tower_sum_energy_TOWERINFO_CALIB_CEMC + tower_sum_energy_TOWERINFO_CALIB_HCALIN + tower_sum_energy_TOWERINFO_CALIB_HCALOUT;
//         float sum_et_sub1 = tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 + tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1 + tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1;
        
//         float xaxis_var = sum_et;
//         if ( x_axis_cent ) { xaxis_var = 1.0*centrality; }

//         int ncomp_cone_corr = random_cone_num_towers_RandomCones_r04 - random_cone_num_masked_towers_RandomCones_r04;
//         int ncomp_cone_sub1_corr = random_cone_num_towers_RandomCones_r04_Sub1 - random_cone_num_masked_towers_RandomCones_r04_Sub1;

//         float area_bkgd = rho_val_TowerRho_AREA*AREA_CONE;
//         float mult_bkgd = ( (rho_val_TowerRho_MULT_HCALIN * random_cone_num_towers_hcalin_RandomCones_r04) 
//                             + (rho_val_TowerRho_MULT_HCALOUT * random_cone_num_towers_hcalout_RandomCones_r04) 
//                             + (rho_val_TowerRho_MULT_CEMC * random_cone_num_towers_cemc_RandomCones_r04) );
//                             // float mult_bkgd = rho_val_TowerRho_MULT_HCALIN
//         float cone_res_area  = random_cone_energy_RandomCones_r04 - area_bkgd;
//         float cone_res_mult  = random_cone_energy_RandomCones_r04 - mult_bkgd;
//         float cone_res_sub1  = random_cone_energy_RandomCones_r04_Sub1;

        
//         // nominal
//         float sum_et2_cemc = N_CEMC_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_CEMC*tower_std_energy_TOWERINFO_CALIB_CEMC) 
//                                 + ( tower_avg_energy_TOWERINFO_CALIB_CEMC*tower_avg_energy_TOWERINFO_CALIB_CEMC) );
//         float sum_et2_hcalin = N_HCALIN_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_HCALIN*tower_std_energy_TOWERINFO_CALIB_HCALIN)
//                                 + ( tower_avg_energy_TOWERINFO_CALIB_HCALIN*tower_avg_energy_TOWERINFO_CALIB_HCALIN) );
//         float sum_et2_hcalout = N_HCALOUT_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_HCALOUT*tower_std_energy_TOWERINFO_CALIB_HCALOUT)
//                                 + ( tower_avg_energy_TOWERINFO_CALIB_HCALOUT*tower_avg_energy_TOWERINFO_CALIB_HCALOUT) );
//         float total_sum_et2 = sum_et2_cemc + sum_et2_hcalin + sum_et2_hcalout;
        
//         float total_fired = ( ( N_CEMC_TOWERS * 1.0 * tower_frac_fired_TOWERINFO_CALIB_CEMC )     
//                             + ( N_HCALIN_TOWERS * 1.0 * tower_frac_fired_TOWERINFO_CALIB_HCALIN )
//                             + ( N_HCALOUT_TOWERS * 1.0 * tower_frac_fired_TOWERINFO_CALIB_HCALOUT ) );

//         float total_unmasked = ( N_CEMC_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_CEMC) 
//                                 + N_HCALIN_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_HCALIN) 
//                                 + N_HCALOUT_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_HCALOUT) );

//         // sub1
//         float sum_et2_cemc_sub1 = N_HCALIN_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1*tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1) 
//                                 + ( tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1*tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1) );
//         float sum_et2_hcalin_sub1 = N_HCALIN_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1*tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1)
//                                 + ( tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1*tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1) );
//         float sum_et2_hcalout_sub1 = N_HCALOUT_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1*tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1)
//                                 + ( tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1*tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1) );
//         float total_sum_et2_sub1 = sum_et2_cemc_sub1 + sum_et2_hcalin_sub1 + sum_et2_hcalout_sub1;

//         float total_fired_sub1 = ( ( N_HCALIN_TOWERS * tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 )     
//                             + ( N_HCALIN_TOWERS * tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1 )
//                             + ( N_HCALOUT_TOWERS * tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1 ) );
//         float total_unmasked_sub1 = ( N_HCALIN_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1) 
//                                 + N_HCALIN_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1) 
//                                 + N_HCALOUT_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1) );


//         // fill arrays
//         int THIS_BIN = h2_area_cone_res_vs_x->GetXaxis()->FindBin(xaxis_var);
//         THIS_BIN--;// binning starts at 1

//         SUM_ET2[THIS_BIN] += total_sum_et2;
//         SUM_E[THIS_BIN] += sum_et;
//         TOTAL_TOWERS_NOT_MASKED[THIS_BIN] += total_unmasked;
//         TOTAL_TOWERS_FIRED[THIS_BIN] += total_fired;
//         N_TOWERS_TOTAL[THIS_BIN] += (N_CEMC_TOWERS + N_HCALIN_TOWERS + N_HCALOUT_TOWERS);
//         SUM_CONE_COMPS[THIS_BIN] += ncomp_cone_corr;
//         N_CONES_THIS_BIN[THIS_BIN]+=1.0;

//         SUM_ET2_SUB1[THIS_BIN] += total_sum_et2_sub1;
//         SUM_E_SUB1[THIS_BIN] += sum_et_sub1;
//         TOTAL_TOWERS_NOT_MASKED_SUB1[THIS_BIN] += total_unmasked_sub1;
//         TOTAL_TOWERS_FIRED_SUB1[THIS_BIN] += total_fired_sub1;
//         N_TOWERS_TOTAL_SUB1[THIS_BIN] += (N_HCALIN_TOWERS + N_HCALIN_TOWERS + N_HCALOUT_TOWERS);
//         SUM_CONE_COMPS_SUB1[THIS_BIN] += ncomp_cone_sub1_corr;
//         N_CONES_THIS_BIN_SUB1[THIS_BIN]+=1.0;

//         // fill histograms
//         h2_area_cone_res_vs_x->Fill(xaxis_var, cone_res_area);
//         h2_mult_cone_res_vs_x->Fill(xaxis_var, cone_res_mult);
//         h2_sub1_cone_res_vs_x->Fill(xaxis_var, cone_res_sub1);
//         // h2_area_cone_ncomp_vs_x->Fill(xaxis_var, ncomp_cone_corr);
//         // h2_mult_cone_ncomp_vs_x->Fill(xaxis_var, cone_res_mult);
//         // h2_sub1_cone_ncomp_vs_x->Fill(xaxis_var, ncomp_cone_sub1_corr);
//         h2_area_cone_res_vs_course_x->Fill(xaxis_var, cone_res_area);
//         h2_mult_cone_res_vs_course_x->Fill(xaxis_var, cone_res_mult);
//         h2_sub1_cone_res_vs_course_x->Fill(xaxis_var, cone_res_sub1);


//     }


//     // float CalcPoisson(const float sum2, const float sum, const float n, const float ncomp, const float ncones){
//     //     if ( n == 0 ) {return 0;}
//     //     if ( ncones == 0 ) {return 0;}
//     //     float Na = ncomp/n;
//     //     float mu = sum/n;
//     //     float sigma2 = sum2/n - (mu*mu);
//     //     float sigma = TMath::Sqrt(Na*sigma2  + Na*(mu*mu));
//     //     return sigma;
//     // }

//     TGraphErrors * g_poission = new TGraphErrors(N_X_CENT_BINS);
//     TGraphErrors * g_poission_alt1 = new TGraphErrors(N_X_CENT_BINS);
//     TGraphErrors * g_poission_alt2 = new TGraphErrors(N_X_CENT_BINS);
//     TGraphErrors * g_poission_sub1 = new TGraphErrors(N_X_CENT_BINS);

//     for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
//         float x = 0, xerr = 0, y = 0, yerr = 0;
//         x = h2_area_cone_res_vs_x->GetXaxis()->GetBinCenter(i+1);
//         // float Na = SUM_CONE_COMPS[i]/N_TOWERS_TOTAL[i];
//         // float AvgEt = SUM_E[i]/N_TOWERS_TOTAL[i];
//         // float Avg2Et = SUM_ET2[i]/N_TOWERS_TOTAL[i];
//         // float Sigma = TMath::Sqrt( (Avg2Et) - (AvgEt*AvgEt) );

//         xerr = h2_area_cone_res_vs_x->GetXaxis()->GetBinWidth(i+1)/2.0;
       
//         y = CalcPoisson(SUM_ET2[i], SUM_E[i], TOTAL_TOWERS_FIRED[i], SUM_CONE_COMPS[i], N_CONES_THIS_BIN[i]);
//         float y_alt = CalcPoissonHarm(SUM_ET2[i], SUM_E[i], TOTAL_TOWERS_FIRED[i], SUM_CONE_COMPS[i], N_CONES_THIS_BIN[i], V2_VALUES[i]/100.0);
//         float y_alt2 = CalcPoissonHarm(SUM_ET2[i], SUM_E[i], TOTAL_TOWERS_FIRED[i], SUM_CONE_COMPS[i], N_CONES_THIS_BIN[i], V2_VALUES[i]/100.0, V3_VALUES[i]/100.0);
//         float y_sub1 = CalcPoisson(SUM_ET2_SUB1[i], SUM_E_SUB1[i], N_TOWERS_TOTAL_SUB1[i], SUM_CONE_COMPS_SUB1[i], N_CONES_THIS_BIN_SUB1[i]);
//         // std::cout << "sum_et2[i]: " << SUM_ET2[i] << std::endl;
//         // std::cout << "sum_et[i]: " << SUM_E[i] << std::endl;
//         // std::cout << "n_towers_total[i]: " << N_TOWERS_TOTAL[i] << std::endl;
//         // std::cout << "ncomp_cone_corr[i]: " << SUM_CONE_COMPS[i] << std::endl;
//         // std::cout << "ncones_this_bin[i]: " << N_CONES_THIS_BIN[i] << std::endl;
//         // std::cout << "avg" << SUM_CONE_COMPS[i]/N_CONES_THIS_BIN[i] << std::endl;
//         // std::cout << "y: " << y << std::endl;
//         g_poission->SetPoint(i, x, y);
//         g_poission->SetPointError(i, xerr, yerr);

//         g_poission_alt1->SetPoint(i, x, y_alt);
//         g_poission_alt1->SetPointError(i, xerr, yerr);

//         g_poission_alt2->SetPoint(i, x, y_alt2);
//         g_poission_alt2->SetPointError(i, xerr, yerr);

//         g_poission_sub1->SetPoint(i, x, y_sub1);
//         g_poission_sub1->SetPointError(i, xerr, yerr);

//     }

    
   
//     TCanvas * c;
//     TLatex * tex = new TLatex();
//     tex->SetNDC();
//     tex->SetTextFont(42);
//     gStyle->SetOptStat(0);
//     gStyle->SetOptFit(0);
//     gStyle->SetOptTitle(0);


//     // start with 3x3 of course x bins
//     c = new TCanvas("c", "c", 400*3, 400*3);
//     c->Divide(3,3);
//     double tx=0.45;
//     double ty_start=0.85;
    
//     TLegend * leg = new TLegend(0.18,0.5,0.35,0.7);
//     leg->SetBorderSize(0);
//     leg->SetFillStyle(0);

//     std::vector<TH2F*> h2s = {h2_area_cone_res_vs_course_x, h2_mult_cone_res_vs_course_x, h2_sub1_cone_res_vs_course_x};
//     std::vector<std::string> labs = {"Area", "Multiplicity", "Iterative"};
//     for ( unsigned islice = 0; islice < N_COURSE_X_BINS; islice++ ){
//         c->cd(islice+1);
//         gPad->SetLogy();
//         gPad->SetLeftMargin(0.15);
//         gPad->SetRightMargin(0.1);
//         gPad->SetBottomMargin(0.15);
//         gPad->SetTopMargin(0.05);

//         double miny = 1e-4;
//         int ihist = 0;
//         std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
//         std::string leg_title = Form("%0.0f < #Sigma E_{T}^{Raw} < %0.0f GeV", COURSE_X_BINS[islice], COURSE_X_BINS[islice+1]);
//         if ( x_axis_cent ) { leg_title = Form("%0.0f-%0.0f %%", COURSE_X_BINS[islice], COURSE_X_BINS[islice+1]); }
//         tags.push_back(leg_title);
//         for ( auto h2 : h2s ) {
//             h2->GetXaxis()->SetRange(islice+1, islice+1);
//             TH1F * h1 = (TH1F*)h2->ProjectionY(Form("h1_%s_%d", h2->GetTitle(), islice), islice+1, islice+1);
//             h1->SetLineColor(COLORS[ihist]);
//             h1->SetMarkerColor(COLORS[ihist]);
//             h1->SetMarkerStyle(MARKERS[ihist]);
//             h1->Scale(1.0/h1->Integral());
//             h1->GetYaxis()->SetRangeUser(miny, 1e1);
//             int lastbin_above_threshold = 0;
//             lastbin_above_threshold = h1->FindLastBinAbove(miny);
//             // h1->GetXaxis()->SetRangeUser(-1.2*h1->GetBinCenter(lastbin_above_threshold), 1.2*h1->GetBinCenter(lastbin_above_threshold));
//             h1->GetXaxis()->SetRangeUser(-40, 40);
//             if ( ihist == 0 ) { h1->Draw("p"); }
//             else { h1->Draw("p same"); }
//             if ( islice == 0 ) { leg->AddEntry(h1, labs[ihist].c_str(), "lp"); }
//             ihist++;
//         }

//         double ty = ty_start;
//         for ( auto tag : tags ) {
//             tex->DrawLatex(tx, ty, tag.c_str());
//             ty -= 0.05;
//         }

//         leg->Draw("same");
//     } 
//     c->SaveAs((outdir+"/cone_res_vs_course_x_slices.png").c_str());

//     delete c;
//     delete leg;

//     const int NX_BINS = h2_area_cone_res_vs_x->GetNbinsX();
//     h2s.clear();
//     leg = new TLegend(0.18,0.5,0.35,0.7);
//     leg->SetBorderSize(0);
//     leg->SetFillStyle(0);
//     h2s = {h2_area_cone_res_vs_x, h2_mult_cone_res_vs_x, h2_sub1_cone_res_vs_x};
//     TGraphErrors * g_std_devs[h2s.size()];
//     for ( unsigned ihist = 0; ihist < h2s.size(); ihist++ ) {
//         g_std_devs[ihist] = new TGraphErrors(NX_BINS);
//     }
//     for ( unsigned ibin = 0; ibin < NX_BINS; ibin++ ) {
//         c = new TCanvas("c", "c", 3*400, 400);
//         c->Divide(3,1);
             
//         double miny = 1e-3;
//         int ihist = 0;
     
        
//         for ( auto h2 : h2s ) {

//             c->cd(ihist+1);
//             gPad->SetLogy();
//             gPad->SetLeftMargin(0.15);
//             gPad->SetRightMargin(0.1);
//             gPad->SetBottomMargin(0.15);
//             gPad->SetTopMargin(0.05);

//             std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
//             tags.push_back(labs[ihist]);
//             std::string leg_title = Form("%0.0f < #Sigma E_{T}^{Raw} < %0.0f GeV", h2->GetXaxis()->GetBinLowEdge(ibin+1), h2->GetXaxis()->GetBinUpEdge(ibin+1)); 
//             if ( x_axis_cent ) { leg_title = Form("%0.0f-%0.0f %%", h2->GetXaxis()->GetBinLowEdge(ibin+1), h2->GetXaxis()->GetBinUpEdge(ibin+1)); }
//             tags.push_back(leg_title);

//             h2->GetXaxis()->SetRange(ibin+1, ibin+1);
//             TH1F * h1 = (TH1F*)h2->ProjectionY(Form("h1_%s_%d", h2->GetTitle(), ibin), ibin+1, ibin+1);
//             h1->SetLineColor(COLORS[ihist]);
//             h1->SetMarkerColor(COLORS[ihist]);
//             h1->SetMarkerStyle(MARKERS[ihist]);
//             h1->Scale(1.0/h1->Integral());
//             h1->GetYaxis()->SetRangeUser(miny, 1e1);
//             int lastbin_above_threshold = 0;
//             lastbin_above_threshold = h1->FindLastBinAbove(miny);
//             // float abs_max_x = 1.2*h1->GetBinCenter(lastbin_above_threshold);
//             float abs_max_x = 40;
//             h1->GetXaxis()->SetRangeUser(-abs_max_x, abs_max_x);

//             // Fit with Gaussian (do not display)
//             float avg = h1->GetMean();
//             float std = h1->GetRMS();

//             int mean_bin = h1->FindBin(avg);
//             h1->GetXaxis()->SetRange(1, mean_bin);
//             TF1 * f1 = new TF1("f1", "gaus", -abs_max_x, abs_max_x);
//             h1->Fit(f1, "RQ", "", -abs_max_x, avg); // left side of peak
//             float mean_left = f1->GetParameter(1);
//             float sigma_left = f1->GetParameter(2);
//             h1->Fit(f1, "RQ", "", mean_left - 1.5*sigma_left, avg); // left side of peak
//             mean_left = f1->GetParameter(1);
//             sigma_left = f1->GetParameter(2);
//             float mean_left_err = f1->GetParError(1);
//             float sigma_left_err = f1->GetParError(2);
//             TF1 * fitfunc = h1->GetFunction("f1");
//             fitfunc->SetLineColor(kRed);
//             fitfunc->SetLineStyle(2);
//             fitfunc->SetLineWidth(2);
//             float chi2 = fitfunc->GetChisquare();

//             g_std_devs[ihist]->SetPoint(ibin, h2->GetXaxis()->GetBinCenter(ibin+1), std);
//             g_std_devs[ihist]->SetPointError(ibin, h2->GetXaxis()->GetBinWidth(ibin+1)/2.0, 0);

//             // h1->GetXaxis()->SetRange(1, h1->GetNbinsX());
//             // TF1 * f2_gamma = new TF1("f2_gamma", "[0]*([1]/TMath::Gamma([2]))*TMath::Power([1]*x + [2], [2]-1)*TMath::Exp(-[1]*x - [2])",-1.2*h1->GetBinCenter(lastbin_above_threshold), 1.2*h1->GetBinCenter(lastbin_above_threshold));
//             // h1->Fit(f2_gamma, "RQ", "", -abs_max_x, abs_max_x);
//             // float gamma_ab = f2_gamma->GetParameter(1);
//             // float gamma_ap = f2_gamma->GetParameter(2);
//             // float gamma_norm = f2_gamma->GetParameter(0);
//             // float gamma_mean = gamma_ap/gamma_ab;
//             // float gamma_sigma = TMath::Sqrt(gamma_ap)/gamma_ab;  
//             // h1->Fit(f2_gamma, "RQ", "", gamma_mean - 1.5*gamma_sigma, gamma_mean + 1.5*gamma_sigma);
//             // gamma_ab = f2_gamma->GetParameter(1);
//             // gamma_ap = f2_gamma->GetParameter(2);
//             // gamma_norm = f2_gamma->GetParameter(0);
//             // gamma_mean = gamma_ap/gamma_ab;
//             // gamma_sigma = TMath::Sqrt(gamma_ap)/gamma_ab;
//             // float gamma_ab_err = f2_gamma->GetParError(1);
//             // float gamma_ap_err = f2_gamma->GetParError(2);
//             // float gamma_mean_err = std::sqrt((gamma_ap_err/gamma_ab)*(gamma_ap_err/gamma_ab) + (gamma_ab_err/gamma_ab)*(gamma_ab_err/gamma_ab))*gamma_mean;
//             // float gamma_sigma_err = std::sqrt((0.5*(gamma_ap_err/gamma_ap))*(0.5*(gamma_ap_err/gamma_ap)) + (gamma_ab_err/gamma_ab)*(gamma_ab_err/gamma_ab))*gamma_sigma;
//             // TF1 * fitfunc_gamma = h1->GetFunction("f2_gamma");
//             // fitfunc_gamma->SetLineColor(kAzure);
//             // fitfunc_gamma->SetLineStyle(2);
//             // fitfunc_gamma->SetLineWidth(2);
//             // float chi2_gamma = fitfunc_gamma->GetChisquare();
//             h1->GetXaxis()->SetRangeUser(-abs_max_x, abs_max_x);
//             h1->Draw("p");
//             leg->AddEntry(h1, Form("#mu = %0.2f, #sigma = %0.2f", avg, std), "p");
//             leg->AddEntry(fitfunc, Form("#mu_{LHS} = %0.2f #pm %0.2f, #sigma_{LHS} = %0.2f #pm %0.2f", mean_left, mean_left_err, sigma_left, sigma_left_err), "l");
//             // leg->AddEntry(fitfunc_gamma, Form("a_{b} = %0.2f #pm %0.2f, a_{p} = %0.2f #pm %0.2f", gamma_ab, gamma_ab_err, gamma_ap, gamma_ap_err), "l");

//             double ty = ty_start;
//             for ( auto tag : tags ) {
//                 tex->DrawLatex(tx, ty, tag.c_str());
//                 ty -= 0.05;
//             }

//             leg->Draw("same");
//             ihist++;
//             c->Update();
//             leg->Clear();
//         }
//         std::cout << "done with slice " << ibin << std::endl;

//         c->SaveAs((outdir+"/cone_res_vs_x_slice_"+std::to_string(ibin)+".png").c_str());
//         leg->Clear();
//         delete c;
//     }

//     c = new TCanvas("c", "c", 800, 600);
//     gPad->SetLeftMargin(0.15);
//     gPad->SetRightMargin(0.1);
//     gPad->SetBottomMargin(0.15);
//     gPad->SetTopMargin(0.05);
//     TLegend * leg1 = new TLegend(0.18,0.5,0.35,0.7);
//     std::vector<TGraphErrors*> g_std_devs_vec = {g_poission, g_poission_alt1, g_poission_alt2, g_poission_sub1};
//     std::vector<std::string> g_std_devs_labels = {"Poisson", "v2", "v3", "Sub1"};
//     for ( unsigned ihist = 0; ihist < g_std_devs_vec.size(); ihist++ ) {
//         g_std_devs_vec[ihist]->SetMarkerStyle(MARKERS[ihist+1]);
//         g_std_devs_vec[ihist]->SetMarkerColor(COLORS[ihist+1]);
//         g_std_devs_vec[ihist]->SetLineColor(COLORS[ihist+1]);
//         g_std_devs_vec[ihist]->GetYaxis()->SetRangeUser(0.0, 12);
//         // g_std_devs[ihist]->SetLineWidth(0);
//         // g_std_devs[ihist]->SetLineStyle(0);
//         leg1->AddEntry(g_std_devs_vec[ihist], g_std_devs_labels[ihist].c_str(), "p");
//         if ( ihist == 0 ) { g_std_devs_vec[ihist]->Draw("ap"); }
//         else { g_std_devs_vec[ihist]->Draw("p same"); }
//     }
//     leg1->Draw("same");
//     c->SaveAs((outdir+"/cone_res_vs_pois.png").c_str());
//     delete c;



//     TLegend * leg2 = new TLegend(0.18,0.5,0.35,0.7);
//     leg2->SetBorderSize(0);
//     leg2->SetFillStyle(0);
//     c = new TCanvas("c", "c", 800, 600);
//     gPad->SetLeftMargin(0.15);
//     gPad->SetRightMargin(0.1);
//     gPad->SetBottomMargin(0.15);
//     gPad->SetTopMargin(0.05);
//     g_poission->SetMarkerStyle(MARKERS[0]);
//     g_poission->SetMarkerColor(COLORS[0]);
//     g_poission->SetLineColor(COLORS[0]);
//     // g_poission->SetLineWidth(0);
//     // g_poission->SetLineStyle(0);
//     // g_pois->GetXaxis()->SetTitle("#Sigma Q_{MBD}");
//     if ( x_axis_cent ) { g_poission->GetXaxis()->SetTitle("Centrality [%]"); }
//     else { g_poission->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]"); }
//     g_poission->GetYaxis()->SetTitle("#sigma(#delta E_{T}) [GeV]");
//     g_poission->GetYaxis()->SetRangeUser(0.1, 12);
//     g_poission->GetXaxis()->SetRangeUser(0, 1800);
//     g_poission->GetXaxis()->SetNdivisions(505);
//     g_poission->Draw("ap");
//     leg2->AddEntry(g_poission, "Poisson Limit", "pl");

//     int mycolors[3] = {kRed, kBlue, kCyan};
//     for ( unsigned ihist = 0; ihist < h2s.size(); ihist++ ) {
       
//         g_std_devs[ihist]->SetMarkerStyle(MARKERS[ihist+1]);
//         g_std_devs[ihist]->SetMarkerColor(mycolors[ihist]);
//         g_std_devs[ihist]->SetMarkerSize(1.5);
//         g_std_devs[ihist]->GetXaxis()->SetNdivisions(505);
//         // g_std_devs[ihist]->GetXaxis()->SetRangeUser(0, 1800);
//         g_std_devs[ihist]->GetYaxis()->SetTitle("#sigma_{LHS} [GeV]");
//         if ( x_axis_cent ) { g_std_devs[ihist]->GetXaxis()->SetTitle("Centrality [%]"); }
//         else { g_std_devs[ihist]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]"); }
//         g_std_devs[ihist]->GetYaxis()->SetRangeUser(0.1, 12);
//         // if ( ihist == 0 ) { g_std_devs[ihist]->Draw("ap"); }
//         // else { g_std_devs[ihist]->Draw("p same"); }
//         g_std_devs[ihist]->Draw("p same");
//         leg2->AddEntry(g_std_devs[ihist], labs[ihist].c_str(), "p");
        

//     }
//     leg2->Draw("same");
//     std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
//     double tx2=0.18;
//     double ty2=0.89;
//     for ( auto tag : tags ) {
//         tex->DrawLatex(tx2, ty2, tag.c_str());
//         ty2 -= 0.05;
//     }
//     c->SaveAs((outdir+"/cone_res_vs_sumq.png").c_str());



//     f->Close();

//     return ;
    
// }

void ConfigureOutputDirs(std::string input_file_base, std::string plotting_dir)
{
    
    if ( !gSystem->OpenDirectory(plotting_dir.c_str()) ) {
        gSystem->mkdir(plotting_dir.c_str(), true);
    }

    plotting_dir += input_file_base + "/";
    if ( !gSystem->OpenDirectory(plotting_dir.c_str()) ) {
        gSystem->mkdir(plotting_dir.c_str(), true);
    }
    
    global_plots = plotting_dir + "global";
    calo_plots = plotting_dir + "calo_tower";
    calo_window_plots = plotting_dir + "calo_window";
    bkgd_plots = plotting_dir + "background_density";
    random_cone_plots = plotting_dir + "random_cone";
    std::vector<std::string> plot_directories = {global_plots, calo_plots, calo_window_plots, bkgd_plots, random_cone_plots};
    for ( auto const& d : plot_directories ) {
        if ( !gSystem->OpenDirectory(d.c_str()) ) {
            gSystem->mkdir(d.c_str(), true);
        }
    }

    return;
}

std::pair<int,int> GetWindowDimFromString(std::string hist_name)
{
    std::pair<int,int> window_dim;
    std::string sWindow = hist_name.substr(hist_name.find_last_of("_") + 1);
    std::string sWindowX = sWindow.substr(0, sWindow.find("x"));
    std::string sWindowY = sWindow.substr(sWindow.find("x") + 1);
    window_dim.first = std::stoi(sWindowX);
    window_dim.second = std::stoi(sWindowY);
    return window_dim;
}

void ProcessCaloWindowHistograms(const std::string & input_file)
{

    TFile * f = new TFile(input_file.c_str(), "READ");
    if(!f->IsOpen() || f->IsZombie()){ std::cout << "File " << input_file << " is zombie" << std::endl;  exit(1); }

    std::cout <<"Processing HCAL geo Calowindows" << std::endl;
    float miny = 1e-4, maxy = 1.5e1;
    double tagx = 0.59, tagy = 0.91, dty = 0.055;
    double x1 = 0.25, x2 = 0.45, y1 = 0.75, y2 = 0.95;
    double y1_s = 0.85, y2_s = 0.95;

    for ( int idim = 0; idim < k_calo_window_dims_hcal_geom.size(); idim++ ) {

        TH2F * h2_window_energy_cent_recemc_nxm = (TH2F*)f->Get(Form("h2_window_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_window_energy_cent_recemc_nxm ) { std::cout << "h2_window_energy_cent_recemc_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH2F * h2_window_energy_cent_hcalin_nxm = (TH2F*)f->Get(Form("h2_window_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_window_energy_cent_hcalin_nxm ) { std::cout << "h2_window_energy_cent_hcalin_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH2F * h2_window_energy_cent_hcalout_nxm = (TH2F*)f->Get(Form("h2_window_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_window_energy_cent_hcalout_nxm ) { std::cout << "h2_window_energy_cent_hcalout_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH2F * h2_window_energy_cent_full_nxm = (TH2F*)f->Get(Form("h2_window_energy_cent_full_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_window_energy_cent_full_nxm ) { std::cout << "h2_window_energy_cent_full_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        
        CaloWindowMultiPanel({h2_window_energy_cent_recemc_nxm, h2_window_energy_cent_hcalin_nxm, h2_window_energy_cent_hcalout_nxm}, 
                            {"EMCal", "iHCal", "oHCal"},
                            calo_window_plots+"/energy_cent_slices", "sub_calo_windows_et_all_cent_slices", 
                            Form("E_{T}^{%d #times %d} [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            tagx, tagy, dty, x1, x2, y1, y2, 0., 0, miny, maxy);

        CaloWindowMultiPanel({h2_window_energy_cent_full_nxm},
                            {"Total"},
                            calo_window_plots+"/energy_cent_slices", "calo_window_et_all_cent_slices", 
                            Form("E_{T}^{%d #times %d} [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            tagx, tagy, dty, x1, x2, y1_s, y2_s, 0., 0, miny, maxy);

        // delete h2_window_energy_cent_recemc_nxm;
        delete h2_window_energy_cent_hcalin_nxm;
        delete h2_window_energy_cent_hcalout_nxm;
        delete h2_window_energy_cent_full_nxm;

        TH2F * h2_recemc_energy_minus_avg_energy_cent_nxm = (TH2F*)f->Get(Form("h2_window_energy_minus_avg_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_recemc_energy_minus_avg_energy_cent_nxm ) { std::cout << "h2_recemc_energy_minus_avg_energy_cent_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH2F * h2_hcalin_energy_minus_avg_energy_cent_nxm = (TH2F*)f->Get(Form("h2_window_energy_minus_avg_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_hcalin_energy_minus_avg_energy_cent_nxm ) { std::cout << "h2_hcalin_energy_minus_avg_energy_cent_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH2F * h2_hcalout_energy_minus_avg_energy_cent_nxm = (TH2F*)f->Get(Form("h2_window_energy_minus_avg_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_hcalout_energy_minus_avg_energy_cent_nxm ) { std::cout << "h2_hcalout_energy_minus_avg_energy_cent_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH2F * h2_full_energy_minus_avg_energy_cent_nxm = (TH2F*)f->Get(Form("h2_window_energy_minus_avg_energy_cent_full_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_full_energy_minus_avg_energy_cent_nxm ) { std::cout << "h2_full_energy_minus_avg_energy_cent_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        
        CaloWindowMultiPanel({h2_recemc_energy_minus_avg_energy_cent_nxm, h2_hcalin_energy_minus_avg_energy_cent_nxm, h2_hcalout_energy_minus_avg_energy_cent_nxm}, 
                            {"EMCal", "iHCal", "oHCal"},
                            calo_window_plots+"/energy_minus_avg_cent_slices", "sub_calo_windows_et-avget_all_cent_slices", 
                            Form("E_{T}^{%d #times %d} - #LT E_{T}^{%d #times %d} #GT [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second, k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            tagx, tagy, dty, x1, x2, y1, y2, 0., 0, miny, maxy, true);

        CaloWindowMultiPanel({h2_full_energy_minus_avg_energy_cent_nxm},
                            {"Total"},
                            calo_window_plots+"/energy_minus_avg_cent_slices", "calo_window_et-avget_all_cent_slices", 
                            Form("E_{T}^{%d #times %d}- #LT E_{T}^{%d #times %d} #GT [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second, k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            tagx, tagy, dty, x1, x2, y1_s, y2_s, 0., 0, miny, maxy, true);
        
        // delete h2_recemc_energy_minus_avg_energy_cent_nxm;
        delete h2_hcalin_energy_minus_avg_energy_cent_nxm;
        delete h2_hcalout_energy_minus_avg_energy_cent_nxm;
        delete h2_full_energy_minus_avg_energy_cent_nxm;


        // TH3F * h3_energy_frac_energy_cent_recemc_nxm = (TH3F*)f->Get(Form("h3_window_energy_frac_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        // if ( !h3_energy_frac_energy_cent_recemc_nxm ) { std::cout << "h3_energy_frac_energy_cent_recemc_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        // TH3F * h3_energy_frac_energy_cent_hcalin_nxm = (TH3F*)f->Get(Form("h3_window_energy_frac_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        // if ( !h3_energy_frac_energy_cent_hcalin_nxm ) { std::cout << "h3_energy_frac_energy_cent_hcalin_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        // TH3F * h3_energy_frac_energy_cent_hcalout_nxm = (TH3F*)f->Get(Form("h3_window_energy_frac_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        // if ( !h3_energy_frac_energy_cent_hcalout_nxm ) { std::cout << "h3_energy_frac_energy_cent_hcalout_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        
        // CaloWindowMultiPanel3D({h3_energy_frac_energy_cent_recemc_nxm, h3_energy_frac_energy_cent_hcalin_nxm, h3_energy_frac_energy_cent_hcalout_nxm}, 
        //                     {"EMCal", "iHCal", "oHCal"},
        //                     calo_window_plots+"/energy_frac_energy_slices", "sub_calo_windows_frac_et_all_cent_slices", 
        //                     Form("E_{T}^{%d #times %d} [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
        //                     Form("#LT f(E_{T}^{%d #times %d}) #GT [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
        //                     .45, .85, dty, .2, .4, .6, .9, -1, -1, 0, 1.5, false, false);
        
        // delete h3_energy_frac_energy_cent_recemc_nxm;
        // delete h3_energy_frac_energy_cent_hcalin_nxm;
        // delete h3_energy_frac_energy_cent_hcalout_nxm;
        
        TH2F * h2_window_frac_energy_cent_recemc_nxm = (TH2F*)f->Get(Form("h2_window_frac_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_window_frac_energy_cent_recemc_nxm ) { std::cout << "h2_window_frac_energy_cent_recemc_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH2F * h2_window_frac_energy_cent_hcalin_nxm = (TH2F*)f->Get(Form("h2_window_frac_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_window_frac_energy_cent_hcalin_nxm ) { std::cout << "h2_window_frac_energy_cent_hcalin_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH2F * h2_window_frac_energy_cent_hcalout_nxm = (TH2F*)f->Get(Form("h2_window_frac_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_window_frac_energy_cent_hcalout_nxm ) { std::cout << "h2_window_frac_energy_cent_hcalout_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        CaloWindowMultiPanel({h2_window_frac_energy_cent_recemc_nxm, h2_window_frac_energy_cent_hcalin_nxm, h2_window_frac_energy_cent_hcalout_nxm}, 
                            {"EMCal", "iHCal", "oHCal"},
                            calo_window_plots+"/frac_energy_cent_slices", "sub_calo_windows_frac_et_all_cent_slices", 
                            Form("f(E_{T}^{%d #times %d}) [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            "Probability Density",
                            tagx, tagy, dty, x1, x2, y1, y2, 0., 0, miny, maxy);
                            
        delete h2_window_frac_energy_cent_recemc_nxm;
        delete h2_window_frac_energy_cent_hcalin_nxm;
        delete h2_window_frac_energy_cent_hcalout_nxm;
    }

    std::cout << "Processing CEMC geo Calowindows" << std::endl;
    // for ( int idim = 0; idim < k_calo_window_dims_cemc_geom.size(); idim++ ) {

    //     TH2F * h2_window_energy_cent_cemc_nxm = (TH2F*)f->Get(Form("h2_window_energy_cent_cemc_%dx%d", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second));
    //     if ( !h2_window_energy_cent_cemc_nxm ) { 
    //         std::cout << "h2_window_energy_cent_cemc_" << k_calo_window_dims_cemc_geom[idim].first << "x" << k_calo_window_dims_cemc_geom[idim].second << " not found!" << std::endl; 
    //         continue;
    //     }
    //     CaloWindowMultiPanel({h2_window_energy_cent_cemc_nxm},
    //                         {"EMCal (0.025#times 0.025)"},
    //                         calo_window_plots+"/energy_cent_slices/cemc_geo", "calo_window_et_all_cent_slices",
    //                         Form("E_{T}^{%d #times %d} [GeV]", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second),
    //                         Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second),
    //                         tagx, tagy, dty, x1, x2, y1_s, y2_s, 0., 0, miny, maxy);
        
    //     delete h2_window_energy_cent_cemc_nxm;

    //     TH2F * h2_cemc_energy_minus_avg_energy_cent = (TH2F*)f->Get(Form("h2_window_energy_minus_avg_energy_cent_cemc_%dx%d", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second));
    //     if ( !h2_cemc_energy_minus_avg_energy_cent ) { 
    //         std::cout << "h2_cemc_energy_minus_avg_energy_cent_" << k_calo_window_dims_cemc_geom[idim].first << "x" << k_calo_window_dims_cemc_geom[idim].second << " not found!" << std::endl; 
    //         continue;
    //     }
    //     CaloWindowMultiPanel({h2_cemc_energy_minus_avg_energy_cent},
    //                         {"EMCal (0.025#times 0.025)"},
    //                         calo_window_plots+"/energy_minus_avg_cent_slices/cemc_geo", "calo_window_et-avget_all_cent_slices",
    //                         Form("E_{T}^{%d #times %d} [GeV] - #LT E_{T}^{%d #times %d} #GT [GeV]", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second, k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second),
    //                         Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second),
    //                         tagx, tagy, dty, x1, x2, y1_s, y2_s, 0., 0, miny, maxy, true);

    // }

    std::cout << "Processing HCAL geo Calowindows" << std::endl;
    f->Close();
    return;
}

void CaloWindowMultiPanel(std::vector<TH2F*> h2_vec, const std::vector<std::string> h2_titles,
                        const std::string output_location, const std::string pngbase,
                        const std::string xaxis_title, const std::string yaxis_title, 
                        const double tagx, const double tagy, const double dty, 
                        const double lx1, const double lx2, const double ly1, const double ly2,
                        float minx, float maxx, const float miny, const float maxy,
                        bool sym_x )
{

    if ( !gSystem->OpenDirectory(output_location.c_str()) ) {
        gSystem->mkdir(output_location.c_str(), true);
    }

    std::string outrootfile = output_location + "/" + pngbase + "_histograms.root";
    TFile * histo_file = new TFile(outrootfile.c_str(), "RECREATE");


    auto window_dim = GetWindowDimFromString(h2_vec.at(0)->GetName());
    std::vector<std::string> sTags = {sPHENIX_Tag};
    sTags.push_back(Form("%d #times %d Towers", window_dim.first, window_dim.second));

    double tagy_start = tagy;

    int NX_PADS= N_CENT_BINS % 3 == 0 ? N_CENT_BINS / 3 : N_CENT_BINS / 3 + 1;
    int NY_PADS = 3;

    TCanvas * c = new TCanvas("c", "c", 400*NX_PADS, 400*NY_PADS);
    c->Divide(NX_PADS, NY_PADS);
    for (int i = 1; i <= N_CENT_BINS; i++) {
        c->cd(i);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.01);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.01);
    }

    TLegend * leg = new TLegend(lx1, ly1, lx2, ly2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    for ( unsigned int i = 0; i < h2_vec.size(); i++ ){
        h2_vec.at(i)->SetTitle(h2_titles.at(i).c_str());
    }

    // get min and max x values for all histograms
    for ( const auto & h2 : h2_vec ) {
        if ( !h2 ) {
            continue;
        }
        for ( int icent = 0; icent < N_CENT_BINS; icent++ ) {
            h2->GetYaxis()->SetRangeUser(CENT_BINS[icent], CENT_BINS[icent+1]);
            TH1F * h1 = (TH1F *) h2->ProjectionX(Form("h1_%s_cent%d", h2->GetName(), icent));
            if ( h1->Integral() == 0 ) { continue; }
            h1->Scale(1./h1->Integral());
            int max_bin = h1->FindLastBinAbove(miny);
            int min_bin = h1->FindFirstBinAbove(miny);
            if ( max_bin < 0 ) { max_bin = h1->GetNbinsX(); }
            if ( min_bin < 0 ) { min_bin = 1; }
            if ( max_bin + 1 < h1->GetNbinsX() ) { max_bin++; }
            if ( min_bin - 1 > 0 ) { min_bin--; }
            if ( h1->GetXaxis()->GetBinCenter(max_bin) > maxx ) {
                maxx = h1->GetXaxis()->GetBinCenter(max_bin);
            }
            if ( h1->GetXaxis()->GetBinCenter(min_bin) < minx ) {
                minx = h1->GetXaxis()->GetBinCenter(min_bin);
            }
            if (sym_x) {
                minx = -maxx;
            }
        }
    }
    bool is_empty_pannel = true;
    for ( int icent = 0; icent < N_CENT_BINS; icent++ ) {
        tagy_start = tagy;
        std::string sCent = Form("%d-%d%% Central", (int)CENT_BINS[icent], (int)CENT_BINS[icent+1]);
        c->cd(icent+1);
        for ( unsigned int ihist = 0; ihist < h2_vec.size(); ihist++ ) {
            if ( !h2_vec.at(ihist) ) {
                continue;
            }

            h2_vec.at(ihist)->GetYaxis()->SetRangeUser(CENT_BINS[icent], CENT_BINS[icent+1]);
            TH1F * h1 = (TH1F *) h2_vec.at(ihist)->ProjectionX(Form("h1_%s_cent%d", h2_vec.at(ihist)->GetName(), icent));
            if ( h1->Integral() == 0 ) { continue; }
            h1->Scale(1./h1->Integral());
            is_empty_pannel= false;
            
            std::string hist_title = h2_vec.at(ihist)->GetTitle();
            h1->SetTitle("");
            h1->GetXaxis()->SetTitle(xaxis_title.c_str());
            h1->GetYaxis()->SetTitle(yaxis_title.c_str());
            
            h1->SetLineColor(COLORS[ihist]);
            h1->SetLineWidth(2);
            
            h1->SetMarkerColor(COLORS[ihist]);
            h1->SetMarkerStyle(MARKERS[ihist]);
            h1->SetMarkerSize(MARKER_SIZE);

            h1->GetXaxis()->SetRangeUser(minx, maxx);
            h1->GetYaxis()->SetRangeUser(miny, maxy);

            if ( ihist == 0 ) {
                h1->Draw("");
            } else {
                h1->Draw("SAME");
            }
            if ( icent == 0 ) {
                leg->AddEntry( h1, hist_title.c_str(), "lep" );
            }
            histo_file->cd();
            h1->Write();
        }
        leg->Draw("SAME");
        for ( auto const& s : sTags ) {
            tex->DrawLatex(tagx, tagy_start, s.c_str());
            tagy_start -= dty;
        }
        tex->DrawLatex(tagx, tagy_start, sCent.c_str());
        c->Update();
    }
    if ( !is_empty_pannel ) {
        c->SaveAs(Form("%s/%s_%dx%d.png", output_location.c_str(), pngbase.c_str(), window_dim.first, window_dim.second));
    }
    delete c;
    histo_file->Close();
    return;
}

void CaloWindowMultiPanel3D(std::vector<TH3F*> h3s, const std::vector<std::string> h3_titles,
                        const std::string output_location, const std::string pngbase,
                        const std::string xaxis_title, const std::string yaxis_title, // after projection
                        const double tagx, const double tagy, const double dty, 
                        const double lx1, const double lx2, const double ly1, const double ly2,
                        float minx, float maxx, const float miny, const float maxy,
                        bool sym_x, bool logy )
{

    if ( !gSystem->OpenDirectory(output_location.c_str()) ) {
        gSystem->mkdir(output_location.c_str(), true);
    }

    auto window_dim = GetWindowDimFromString(h3s.at(0)->GetName());
    std::vector<std::string> sTags = {sPHENIX_Tag};
    sTags.push_back(Form("%d #times %d Towers", window_dim.first, window_dim.second));

    double tagy_start = tagy;
    int NX_PADS= h3s.size(), NY_PADS = 1;
    const float cent_bins[] = {0, 10, 20, 40, 60};
    const int n_cent_bins = sizeof(cent_bins)/sizeof(cent_bins[0]) - 1;

    TCanvas * c = new TCanvas("c", "c", 400*NX_PADS, 400*NY_PADS);
    c->Divide(NX_PADS, NY_PADS);
    for (int i = 1; i <= n_cent_bins; i++) {
        c->cd(i);
        if (logy) {  gPad->SetLogy(logy); }
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
    }

    TLegend * leg = new TLegend(lx1, ly1, lx2, ly2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    
    TLegend * leg2 = new TLegend(lx1, ly1-.45, lx2, ly2-.45);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    float maxyy=0;
    for ( unsigned int i = 0; i < h3s.size(); i++ ){
        h3s.at(i)->SetTitle(h3_titles.at(i).c_str());
        for ( int icent = 0; icent < n_cent_bins; icent++ ) {
            h3s.at(i)->GetZaxis()->SetRangeUser(cent_bins[icent], cent_bins[icent+1]);
            TH2F * h2 = (TH2F *) h3s.at(i)->Project3D("xy");
            h2->SetName(Form("h2_%s_cent%d", h2->GetName(), icent));
            TH1F * h1 = ( TH1F * ) h2->ProfileY(Form("p1_%s_cent%d", h2->GetName(), icent));
            int max_bin = h1->FindLastBinAbove(miny);
            int min_bin = h1->FindFirstBinAbove(miny);
            if ( max_bin < 0 ) { max_bin = h1->GetNbinsX(); }
            if ( min_bin < 0 ) { min_bin = 1; }
            if ( max_bin + 1 < h1->GetNbinsX() ) { max_bin++; }
            if ( min_bin - 1 > 0 ) { min_bin--; }
            if ( h1->GetXaxis()->GetBinCenter(max_bin) > maxx ) {
                maxx = h1->GetXaxis()->GetBinCenter(max_bin);
            }
            if ( h1->GetXaxis()->GetBinCenter(min_bin) < minx ) {
                minx = h1->GetXaxis()->GetBinCenter(min_bin);
            }
            // if (h1->GetMaximum() > maxyy) {
            //     maxyy = h1->GetMaximum();
            // }
        }
    }

    bool is_empty_pannel = true;
    for ( unsigned int ihist = 0; ihist < h3s.size(); ihist++ ) {
        if ( !h3s.at(ihist) ) {
            continue;
        }
        tagy_start = tagy;
        c->cd(ihist+1);
        h3s.at(ihist)->Scale(1./h3s.at(ihist)->Integral());
        for ( int icent = 0; icent < n_cent_bins; icent++ ) {
            
            std::string sCent = Form("%d-%d%%", (int)cent_bins[icent], (int)cent_bins[icent+1]);
            h3s.at(ihist)->GetZaxis()->SetRangeUser(cent_bins[icent], cent_bins[icent+1]);

            TH2F * h2 = (TH2F *) h3s.at(ihist)->Project3D("xy");
            h2->SetName(Form("h2_%s_cent%d", h2->GetName(), icent));
            TH1F * h1 = ( TH1F * ) h2->ProfileY(Form("p1_%s_cent%d", h2->GetName(), icent));
            if ( h1->Integral() == 0 ) { continue; }
            is_empty_pannel= false;
                
            h1->SetTitle("");
            h1->GetXaxis()->SetTitle(xaxis_title.c_str());
            h1->GetYaxis()->SetTitle(yaxis_title.c_str());
            
            h1->SetLineColor(COLORS[icent]);
            h1->SetLineWidth(2);
            
            h1->SetMarkerColor(COLORS[icent]);
            h1->SetMarkerStyle(MARKERS[icent]);
            h1->SetMarkerSize(1.0);

            if (minx != maxx){
                h1->GetXaxis()->SetRangeUser(minx, maxx);
            }
            h1->GetYaxis()->SetRangeUser(miny, maxyy);

            if ( icent == 0 ) {
                h1->Draw("");
            } else {
                h1->Draw("SAME");
            }
            if ( ihist == 0 ) {
                leg2->AddEntry( h1, sCent.c_str(), "lep" );
                leg->AddEntry( h1, sCent.c_str(), "lep" );
            }
        }
        if ( ihist == 0 ) {
            leg2->Draw("SAME");
        } else {
            leg->Draw("SAME");
        }
   
        if ( ihist == 0 ) {
            for ( auto const& s : sTags ) {
                tex->DrawLatex(tagx, tagy_start-0.5, s.c_str());
                tagy_start -= dty;
            }
            tex->DrawLatex(tagx, tagy_start-0.5, h3_titles.at(ihist).c_str());
        } else {
            for ( auto const& s : sTags ) {
                tex->DrawLatex(tagx, tagy_start, s.c_str());
                tagy_start -= dty;
            }
            tex->DrawLatex(tagx, tagy_start, h3_titles.at(ihist).c_str());
        }
        c->Update();
    }
    if ( !is_empty_pannel ) {
        c->SaveAs(Form("%s/%s_%dx%d.png", output_location.c_str(), pngbase.c_str(), window_dim.first, window_dim.second));
    }
    delete c;
    return;
}

