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

//-------------------------------------------------------------------------------------------
// Style
//-------------------------------------------------------------------------------------------

const float CENT_BINS[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 100};
const int N_CENT_BINS = sizeof(CENT_BINS)/sizeof(CENT_BINS[0]) - 1;

bool IS_DATA = false;
int NEVENTS = 0;
const std::string sPHENIX_Tag = "#it{#bf{sPHENIX}} Internal";
std::string DataType_Tag;

const int COLORS[] = {kBlack, kRed , kBlue, kGreen+2, kViolet, kCyan, kOrange+2, kMagenta+2, kAzure-2};
const int MARKERS[] = { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kOpenCircle, kOpenSquare, kOpenTriangleUp};
const float MARKER_SIZE = 1.2;
const float LINE_WIDTH = 2.0;


//-------------------------------------------------------------------------------------------
// IO
//-------------------------------------------------------------------------------------------

std::string global_plots;
std::string calo_plots;
std::string calo_window_plots;
std::string bkgd_plots;
std::string random_cone_plots;

TFile * histo_file {nullptr};
TFile * f {nullptr};

void ConfigureOutputDirs(std::string input_file_base, std::string plotting_dir = "plots/");

//-------------------------------------------------------------------------------------------
// CaloWindow
//-------------------------------------------------------------------------------------------

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

void ProcessCaloWindowHistograms();

//-------------------------------------------------------------------------------------------
// CaloTowers
//-------------------------------------------------------------------------------------------

std::string GetCaloNameFromNode(std::string node_name){
    if ( node_name.find("CEMC") != std::string::npos ) {
      if ( node_name.find("RETOWER") != std::string::npos ) {
        if ( node_name.find("SUB1") != std::string::npos ) {
          return "Sub. EMCal";
        } else {
            return "EMCal";
        }
        } else {
            return "EMCal 0.025 #times 0.025";
        }
    } else if ( node_name.find("HCALIN") != std::string::npos ) {
        if ( node_name.find("SUB1") != std::string::npos ) {
            return "Sub. iHCal";
        } else {
            return "iHCal";
        }
    } else if ( node_name.find("HCALOUT") != std::string::npos ) {
        if ( node_name.find("SUB1") != std::string::npos ) {
            return "Sub. oHCal";
        } else {
            return "oHCal";
        }
    }
    return "Unknown";
}

void ProcessCaloHistograms()
{
    std::cout <<"Processing Calo Tower Hists" << std::endl;
    
    std::vector<std::string> group_subtitles = {"Uncorr", "CEMC Geo", "Subtracted"};
    std::vector<std::vector<std::string>> calo_groups = {{"TOWERINFO_CALIB_CEMC_RETOWER", "TOWERINFO_CALIB_HCALIN", "TOWERINFO_CALIB_HCALOUT"},
                                                         {"TOWERINFO_CALIB_CEMC"},
                                                         {"TOWERINFO_CALIB_CEMC_RETOWER_SUB1", "TOWERINFO_CALIB_HCALIN_SUB1", "TOWERINFO_CALIB_HCALOUT_SUB1"}};
    
    std::vector<std::string> th3s_title_bases = {"h3_tower_eta_energy_cent", "h3_tower_phi_energy_cent"}; // 3x3 pannels with <Et> vs eta, phi
    std::vector<std::string> th2s_title_bases = {"h2_tower_energy_cent","h2_tower_dead_phi_eta","h2_tower_energy_phi_eta"}; // // 3x3 pannels with Et vs cent, 
    // nX1 pannel of eff vs eta and phi,


    /// Todo: make this a function

    return;

    // TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC");
    // if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC not found!" << std::endl; }
    // TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC");
    // if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC not found!" << std::endl; }
    // TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_CEMC = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC");
    // if ( !h2_tower_energy_cent_TOWERINFO_CALIB_CEMC ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_CEMC not found!" << std::endl; }
    // TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC");
    // if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC not found!" << std::endl; }

    // // TOWERINFO_CALIB_HCALIN
    // TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN");
    // if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN not found!" << std::endl; }
    // TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN");
    // if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN not found!" << std::endl; }
    // TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN");
    // if ( !h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN not found!" << std::endl; }
    // TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN");
    // if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN not found!" << std::endl; }

    // // TOWERINFO_CALIB_HCALOUT
    // TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT");
    // if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT not found!" << std::endl; }
    // TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT");
    // if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT not found!" << std::endl; }
    // TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT");
    // if ( !h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT not found!" << std::endl; }
    // TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT");
    // if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT not found!" << std::endl; }


    // // TOWERINFO_CALIB_CEMC_RETOWER
    // TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER");
    // if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER not found!" << std::endl; }
    // TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER");
    // if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER not found!" << std::endl; }
    // TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER");
    // if ( !h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER not found!" << std::endl; }
    // TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER");
    // if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER not found!" << std::endl; }

    // // TOWERINFO_CALIB_CEMC_RETOWER_SUB1
    // TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    // if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 not found!" << std::endl; }
    // TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    // if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 not found!" << std::endl; }
    // TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    // if ( !h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 not found!" << std::endl; }
    // TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    // if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 not found!" << std::endl; }

    // // TOWERINFO_CALIB_HCALIN_SUB1
    // TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1");
    // if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 not found!" << std::endl; }
    // TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1");
    // if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 not found!" << std::endl; }
    // TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1");
    // if ( !h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 not found!" << std::endl; }
    // TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN_SUB1 = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN_SUB1");
    // if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN_SUB1 ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN_SUB1 not found!" << std::endl; }

    // // TOWERINFO_CALIB_HCALOUT_SUB1
    // TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1");
    // if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 not found!" << std::endl; }
    // TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1");
    // if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 not found!" << std::endl; }
    // TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1");
    // if ( !h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 not found!" << std::endl; }
    // TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT_SUB1 = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT_SUB1");
    // if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT_SUB1 ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT_SUB1 not found!" << std::endl; }

}

void QuickBins(const int N, const float min, const float max, float *bins){
    if (min >= max) {
        std::cout << "Error: min >= max" << std::endl;
        return;
    }
   float step = 1.1*(max-min)/N;
    for ( int i = 0; i < N; ++i ) {
        bins[i] = min + i*step;
    }
}

TH1F * CreateHisto1D(std::string name, const int N, const float min, const float max, std::pair<std::string,std::string> titles = std::make_pair("","")){
    float bins[N+1];
    QuickBins(N, min, max, bins);
    for ( int i = 0; i < N+1; ++i ) {
        if ( bins[i+1] <= bins[i] ) {
            std::cout << "Error: bins[i+1] <= bins[i for histo " << name << std::endl;
            return nullptr;
        }
    }

    TH1F *h = new TH1F(name.c_str(), name.c_str(), N, bins);
    h->GetXaxis()->SetTitle(titles.first.c_str());
    h->GetYaxis()->SetTitle(titles.second.c_str());
    return h;
}

std::pair<std::string, std::string> QuickXYTitles(std::string var, bool normed = false, std::string units=""){
    std::string xtitle = var;
    std::string ytitle = "Counts";
    if ( normed ) {
        ytitle = "1/N dN/d"+var;
    }
    if ( units != "" ) {
        xtitle += " ["+units+"]";
        ytitle += " ["+units+"^{-1}]";
    }
    return std::make_pair(xtitle, ytitle);
}

void ProcessGlobal(){

    
    // VALID_RUN_NUMBERS=(23727 23735 23737 23738 23739 23740 23743 23745)
    std::string  input_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/JAN30/DATA/BASIC/DATA-ProdA_2023-BASIC-";
    std::vector<int> run_numbers = {23727, 23735, 23737, 23738, 23739, 23740, 23743, 23745};
    TChain *chain = new TChain("T");
    for ( auto run_number : run_numbers ) {
        chain->Add((input_file+"0"+std::to_string(run_number)+".root").c_str());
    }
 

    TTree *t = chain;

    float area_rho = 0;
    float area_rho_cemc, area_rho_ihcal, area_rho_ohcal;
    float mult_rho_cemc = 0, mult_rho_ihcal = 0, mult_rho_ohcal = 0;
    float mult_rho = 0;
    float cemc_e = 0;
    float ihcal_e = 0;
    float ohcal_e = 0;
    float ntowers_cemc, ntowers_ihcal, ntowers_ohcal;
    float ndead_cemc;
    int rc_m_cemc = 0, rc_m_ihcal = 0, rc_m_ohcal = 0;
    int rc_m = 0;
    std::vector<float> * tbackground_cemc = 0;
    std::vector<float> * tbackground_ihcal = 0;
    std::vector<float> * tbackground_ohcal = 0;
    float mbd_n,mbs_s;
    int cent = 0;
    t->SetBranchAddress("centrality", &cent);
    t->SetBranchAddress("mbd_q_N", &mbd_n);
    t->SetBranchAddress("mbd_q_S", &mbs_s);
    t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC", &cemc_e);
    t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN", &ihcal_e);
    t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT", &ohcal_e);
    t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC", &ntowers_cemc);
    t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN", &ntowers_ihcal);
    t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT", &ntowers_ohcal);
    t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC", &ndead_cemc);
    t->SetBranchAddress("rho_val_TowerRho_AREA", &area_rho);
    t->SetBranchAddress("rho_val_TowerRho_AREA_CEMC", &area_rho_cemc);
    t->SetBranchAddress("rho_val_TowerRho_AREA_HCALIN", &area_rho_ihcal);
    t->SetBranchAddress("rho_val_TowerRho_AREA_HCALOUT", &area_rho_ohcal);
    t->SetBranchAddress("rho_val_TowerRho_MULT_HCALIN", &mult_rho_ihcal);
    t->SetBranchAddress("rho_val_TowerRho_MULT_HCALOUT", &mult_rho_ohcal);
    t->SetBranchAddress("rho_val_TowerRho_MULT_CEMC", &mult_rho_cemc);
    t->SetBranchAddress("rho_val_TowerRho_MULT", &mult_rho);
    t->SetBranchAddress("random_cone_num_towers_RandomCones_r04", &rc_m);
    t->SetBranchAddress("random_cone_num_towers_cemc_RandomCones_r04", &rc_m_cemc);
    t->SetBranchAddress("random_cone_num_towers_hcalin_RandomCones_r04", &rc_m_ihcal);
    t->SetBranchAddress("random_cone_num_towers_hcalout_RandomCones_r04", &rc_m_ohcal);
    t->SetBranchAddress("tower_background_energy_recemc", &tbackground_cemc);
    t->SetBranchAddress("tower_background_energy_hcalin", &tbackground_ihcal);
    t->SetBranchAddress("tower_background_energy_hcalout", &tbackground_ohcal);

    float avg_energy_full[11];
    float std_energy_full[11];
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
    t->SetBranchAddress("avg_energy_full", &avg_energy_full);
    t->SetBranchAddress("std_energy_full", &std_energy_full);
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

    int nentries = t->GetEntries();
    const int NET_BINS =20;
    float et_bins[NET_BINS+1];
    double step = 100/NET_BINS; 
    for ( int i = 0; i < NET_BINS+1; ++i ) {
        et_bins[i] = i*step;
    }


    const int NRHOM_BINS = 500;
    double rho_m_max = 65;
    double step_m = rho_m_max/NRHOM_BINS;
    float rho_binsm[NRHOM_BINS+1];
    for ( int i = 0; i < NRHOM_BINS+1; ++i ) {
        rho_binsm[i] = i*step_m;
    }
    const int NRHOA_BINS = 500;
    float rho_binsa[NRHOA_BINS+1];
    const double rho_a_max = 65;
    for ( int i = 0; i < NRHOA_BINS+1; ++i ) {
        rho_binsa[i] = i*rho_a_max/NRHOA_BINS;
    }

   
    std::cout <<"Processing Global" << std::endl;
    TH2F* h2_et_rho_mult = new TH2F("h2_et_rho_mult", "h2_et_rho_mult", NET_BINS, et_bins, NRHOM_BINS, rho_binsm);
    h2_et_rho_mult->GetXaxis()->SetTitle("#Sigma E_{T}^{raw} [GeV]");
    h2_et_rho_mult->GetYaxis()->SetTitle("#rho_{M}#times N_{cone} [GeV]");
   
// const std::vector < std::pair < unsigned int, unsigned int > > k_calo_window_dims_hcal_geom = {
//     {1,1}, {2,2}, {3,4}, {5,6}, {7,8}, {9,10}, {11,12}, {13,13}, {15,15}, {17,17}, {20,20}
// }; 

    
    // const int AVG_ET_BINS = 20000;
    // float avgbins[AVG_ET_BINS+1];
    // const double MAX_AVG  = 150;
    // step = MAX_AVG/AVG_ET_BINS;
    // for ( int i = 0; i < AVG_ET_BINS+1; ++i ) {
    //     avgbins[i] = i*step;
    // }

    const std::vector < std::pair < unsigned int, unsigned int > > k_calo_window_dims_hcal_geom = {
    {1,1}, {2,2}, {3,4}, {5,6}, {7,8}, {9,10}, {11,12}, {13,13}, {15,15}, {17,17}, {20,20}
}; 
    const std::vector <std::pair < unsigned int, unsigned int > > calo_window_dims_hcal_geom = {
        {1,1}, {2,2}, {3,4}, {5,6}, {7,8}, {9,10}, {11,12}, {13,13}, {15,15},
    };
    const unsigned int window_array_size = calo_window_dims_hcal_geom.size(); 
    std::vector<unsigned int> window_array_idx{};
    for ( unsigned int i = 0; i < window_array_size; ++i ) {
        for ( unsigned int j = 0; j < k_window_array_size; ++j ) {
            if ( calo_window_dims_hcal_geom[i] == k_calo_window_dims_hcal_geom[j] ) {
                window_array_idx.push_back(j);
                break;
            }
        }
    }

    TH2F* h2_avg_et_full_window[window_array_size];
    TH2F* h2_avg_et_recemc_window[window_array_size];
    TH2F* h2_avg_et_hcalin_window[window_array_size];
    TH2F* h2_avg_et_hcalout_window[window_array_size];
    const int AVG_ET_BINS = 20000;
    for ( unsigned int i = 0; i < window_array_size; ++i ) {
        
        float avgbins[AVG_ET_BINS+1];
        const double MAX_AVG  = std::sqrt(calo_window_dims_hcal_geom[i].first*calo_window_dims_hcal_geom[i].second)*1.5;
        step = MAX_AVG/AVG_ET_BINS;
        for ( int i = 0; i < AVG_ET_BINS+1; ++i ) {
            avgbins[i] = i*step;
        }
        h2_avg_et_full_window[i] = new TH2F(Form("h2_avg_et_full_window_%dx%d", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second), 
                                        Form("h2_avg_et_full_window_%dx%d", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second), 
                                        NET_BINS, et_bins, AVG_ET_BINS, avgbins);
        h2_avg_et_full_window[i]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        h2_avg_et_full_window[i]->GetYaxis()->SetTitle(Form("#LT E_{T}^{%dx%d} #GT [GeV]", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second));
        h2_avg_et_recemc_window[i] = new TH2F(Form("h2_avg_et_recemc_window_%dx%d", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second), 
                                        Form("h2_avg_et_recemc_window_%dx%d", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second), 
                                        NET_BINS, et_bins, AVG_ET_BINS, avgbins);
        h2_avg_et_recemc_window[i]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        h2_avg_et_recemc_window[i]->GetYaxis()->SetTitle(Form("#LT E_{T}^{%dx%d} #GT [GeV]", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second));
        h2_avg_et_hcalin_window[i] = new TH2F(Form("h2_avg_et_hcalin_window_%dx%d", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second), 
                                        Form("h2_avg_et_hcalin_window_%dx%d", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second), 
                                        NET_BINS, et_bins, AVG_ET_BINS, avgbins);
        h2_avg_et_hcalin_window[i]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        h2_avg_et_hcalin_window[i]->GetYaxis()->SetTitle(Form("#LT E_{T}^{%dx%d} #GT [GeV]", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second));
        h2_avg_et_hcalout_window[i] = new TH2F(Form("h2_avg_et_hcalout_window_%dx%d", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second), 
                                        Form("h2_avg_et_hcalout_window_%dx%d", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second), 
                                        NET_BINS, et_bins, AVG_ET_BINS, avgbins);
        h2_avg_et_hcalout_window[i]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        h2_avg_et_hcalout_window[i]->GetYaxis()->SetTitle(Form("#LT E_{T}^{%dx%d} #GT [GeV]", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second));
    }


 for ( int i = 0; i < nentries; ++i ) {
        t->GetEntry(i);
        // float sum_e = cemc_e+ihcal_e+ohcal_e;
        // float sum_e = mbs_s+mbd_n;
        float sum_e = 1.0*cent;
        for ( unsigned int j = 0; j < window_array_size; ++j ) {
            if ( num_windows_full[window_array_idx[j]] > 0 ) {
                h2_avg_et_full_window[j]->Fill(sum_e, std_energy_full[window_array_idx[j]]);
            }
            if ( num_windows_recemc[j] > 0 ) {
                h2_avg_et_recemc_window[j]->Fill(sum_e, std_energy_recemc[window_array_idx[j]]);
            }
            if ( num_windows_hcalin[j] > 0 ) {
                h2_avg_et_hcalin_window[j]->Fill(sum_e, std_energy_hcalin[window_array_idx[j]]);
            }
            if ( num_windows_hcalout[j] > 0 ) {
                h2_avg_et_hcalout_window[j]->Fill(sum_e, std_energy_hcalout[window_array_idx[j]]);
            }
        }
    }

    // set color palette
    gStyle->SetPalette(kRainBow);
    TCanvas *c = new TCanvas("c", "c", 400*3, 400*3);
   
    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    // tex->SetTextSize(0.04);
    TLegend * leg = new TLegend(0.18,0.5,0.35,0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    c->Divide(3,3);
    // std::vector<float> maxys = {1.2,15,150};
    double tx=0.25;
    double ty=0.85;


    // std::vector<float> maxys = {1.2,15,150};
    TH1F* h_1x1 {nullptr};
    TGraph * g1_nxm[window_array_size];
    for (int i = 0; i < window_array_size; ++i) {
        double tx=0.45;
        double ty=0.25;
        c->cd(i+1);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        // TF1 * f1 = new TF1("f1", "[0]+[1]*x", 200, 1700);

        h2_avg_et_full_window[i]->GetXaxis()->SetNdivisions(505);
        // h2_avg_et_full_window[i]->Fit("f1", "Q", "", 200, 1700);
        // TF1 *fitFunc = h2_avg_et_full_window[i]->GetFunction("f1");
        // fitFunc->SetLineColor(kBlack);
        // fitFunc->SetLineWidth(3);
        // fitFunc->SetLineStyle(2);
       
        h2_avg_et_full_window[i]->GetYaxis()->SetRangeUser(0, h2_avg_et_full_window[i]->GetMaximum()*1.1);
        TH1F * hist_fill = (TH1F*)h2_avg_et_full_window[i]->ProfileX("full", 1, -1);
        TGraphErrors * gr_fill = new TGraphErrors(hist_fill);
        // remove points below 200 and above 1700
        for ( int k = 0; k < gr_fill->GetN(); ++k ) {
            if ( gr_fill->GetX()[k] > 1700 ) {
                gr_fill->RemovePoint(k);
                k--;
            }
        }
        if ( i == 0 ) {
            h_1x1 = (TH1F*)hist_fill->Clone("h_1x1");
        }
        gr_fill->SetLineColor(kBlack);
        gr_fill->SetLineWidth(3);
        gr_fill->SetLineStyle(2);
        gr_fill->SetMarkerColor(kBlack);
        gr_fill->SetMarkerStyle(MARKERS[0]);
        gr_fill->SetMarkerSize(1);
        gr_fill->GetYaxis()->SetTitle(Form("#LT #sigma^{%dx%d}#GT [GeV]", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second));
        gr_fill->GetXaxis()->SetTitle("#Sigma E_{T}^{raw} [GeV]");
        gr_fill->GetXaxis()->SetNdivisions(505);
        gr_fill->Draw("APL");
        TGraphErrors * gr_fill2 = new TGraphErrors(h_1x1);
        float Aa = std::sqrt(calo_window_dims_hcal_geom[i].first*calo_window_dims_hcal_geom[i].second);
        gr_fill2->Scale(Aa);
        gr_fill2->SetLineColor(kRed);
        gr_fill2->SetLineWidth(3);
        gr_fill2->SetLineStyle(2);
        gr_fill2->Draw("PLsame");


        hist_fill->Divide(h_1x1);
        for ( int k = 0; k < hist_fill->GetNbinsX(); ++k ) {
            if ( hist_fill->GetBinContent(k+1) == 0 ) {
                hist_fill->SetBinError(k+1, 0);
            }
        }
        g1_nxm[i] = new TGraph(hist_fill);
        // remove points below 200 and above 1700
        for ( int k = 0; k < g1_nxm[i]->GetN(); ++k ) {
            // if ( g1_nxm[i]->GetX()[k] > 1800 || g1_nxm[i]->GetX()[k] < 200 ) {
                 if ( g1_nxm[i]->GetX()[k] > 80 ) {
                g1_nxm[i]->RemovePoint(k);
                k--;
            }
        }

        g1_nxm[i]->SetLineColor(COLORS[i+1]);
        g1_nxm[i]->SetLineWidth(2);
        g1_nxm[i]->SetLineStyle(1);
        g1_nxm[i]->SetMarkerColor(COLORS[i+1]);
        g1_nxm[i]->SetMarkerStyle(MARKERS[i+1]);
        g1_nxm[i]->SetMarkerSize(1.5);
        // g1_nxm[i]->GetYaxis()->SetTitle("#LT #sigma^{%dx%d}#GT / #LT #sigma^{1x1} #GT", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second);

        // hist_fill[i]->Draw("colz");
       
        // TH1F * hist = (TH1F*)h2_avg_et_recemc_window[i]->ProfileX("projA", 1, -1, "s");
        // TH1F * hist2 = (TH1F*)h2_avg_et_hcalin_window[i]->ProfileX("projB", 1, -1, "s");
        // TH1F * hist3 = (TH1F*)h2_avg_et_hcalout_window[i]->ProfileX("projC", 1, -1, "s");
        // std::vector<TH1F*> histos = {hist, hist2, hist3};
        // std::vector<std::string> tags = { "EMCal", "iHCal", "oHCal"};
        // TGraphErrors * grs[3];
        // for ( unsigned int j = 0; j < histos.size(); ++j ) {

        //     histos[j]->GetXaxis()->SetRange(histos[j]->GetXaxis()->FindBin(200), histos[j]->GetXaxis()->FindBin(1700));
        //     histos[j]->GetXaxis()->SetRangeUser(200, 1700);
        //     grs[j] = new TGraphErrors(histos[j]);
        //     // rm poins below 200 and above 1700    
        //     for ( int k = 0; k < grs[j]->GetN(); ++k ) {
        //         if ( grs[j]->GetX()[k] < 200 || grs[j]->GetX()[k] > 1700 ) {
        //             grs[j]->RemovePoint(k);
        //             k--;
        //         }
        //     }
        //     grs[j]->SetLineColor(COLORS[j+1]);
        //     grs[j]->SetLineWidth(3);
        //     grs[j]->SetLineStyle(2);
        //     grs[j]->SetMarkerColor(COLORS[j+1]);
        //     grs[j]->SetMarkerStyle(MARKERS[j]);
        //     grs[j]->SetMarkerSize(1);
        //     grs[j]->GetXaxis()->SetRangeUser(200, 1700);
        //     if ( i == 0 ) {
        //         leg->AddEntry(grs[j], tags[j].c_str(), "l");
        //     }
        //     grs[j]->Draw("same");
        // }
        
        std::vector<std::string> sTags = {sPHENIX_Tag, Form("%d #times %d Towers", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second)};
        for ( unsigned int j = 0; j < sTags.size(); ++j ) {
            tex->DrawLatex(tx, ty, sTags[j].c_str());
            ty -= 0.05;
        }
        // tex->DrawLatex(0.18, 0.85, Form("Fit: (%.1e #pm %.1e) + (%.1e#pm %.1e) #times E_{T}^{raw}", f1->GetParameter(0), f1->GetParError(0), f1->GetParameter(1), f1->GetParError(1)));

        leg->Draw("same");
    }
    c->SaveAs((global_plots+"/avg_et_window.png").c_str());
    leg->Clear();

    TLegend * leg2 = new TLegend(0.2,0.7,0.9,0.9);
    // 2 col
    leg2->SetNColumns(2);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    TCanvas *c2 = new TCanvas("c2", "c2",1200,800);   
    gPad->SetLeftMargin(0.15);
    // gPad->SetLogy();
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05); 
    for ( int i = 0; i < window_array_size; ++i ) {

        if (i == 0) {
            continue; // skip 1x1
        }
        // g1_nxm[i]->GetYaxis()->SetRangeUser(0, 2);
        g1_nxm[i]->SetLineColor(COLORS[i]);
        g1_nxm[i]->SetLineWidth(2);
        g1_nxm[i]->SetMarkerColor(COLORS[i]);
        g1_nxm[i]->SetMarkerStyle(MARKERS[i]);

        g1_nxm[i]->GetXaxis()->SetTitle("#Sigma E_{T}^{raw} [GeV]");
        g1_nxm[i]->GetYaxis()->SetTitle("#LT #sigma^{NxM}#GT / #LT #sigma^{1x1} #GT");
        g1_nxm[i]->GetXaxis()->SetNdivisions(505);
        g1_nxm[i]->GetYaxis()->SetRangeUser(0.0, 100);
        // g1_nxm[i]->Scale(1/std::sqrt(calo_window_dims_hcal_geom[i].first*calo_window_dims_hcal_geom[i].second));
        float A = 1.0*calo_window_dims_hcal_geom[i].first*calo_window_dims_hcal_geom[i].second;
        for ( int k = 0; k < g1_nxm[i]->GetN(); ++k ) {
            float val = g1_nxm[i]->GetY()[k];
            // should be approx A^1/2
            float valn = TMath::Log(val)/TMath::Log(A) - 0.5;
            g1_nxm[i]->GetY()[k] = valn;
        }
        // g1_nxm[i]->Draw("ALP");
        if ( i == 1 ) {
            g1_nxm[i]->Draw("AP");
        } else {
            g1_nxm[i]->Draw("P same");
        }
        // TF1 * f1 = new TF1("f1", "[0]^[1]", 200, 1900);
        // f1->FixParameter(0, calo_window_dims_hcal_geom[i].first*calo_window_dims_hcal_geom[i].second);
        
        // make a constant
        // g1_nxm[i]->Fit("f1", "Q", "", 200, 1900);
        // TF1 *fitFunc = g1_nxm[i]->GetFunction("f1");
        // fitFunc->SetLineColor(0);
        // fitFunc->SetLineWidth(0);
        // TLine *line = new TLine(200, std::sqrt(calo_window_dims_hcal_geom[i].first*calo_window_dims_hcal_geom[i].second), 1900, std::sqrt(calo_window_dims_hcal_geom[i].first*calo_window_dims_hcal_geom[i].second));
        // line->SetLineColor(kBlack);
        // line->SetLineWidth(2);
        // line->SetLineStyle(2);
        // line->Draw("same");

        // leg2->AddEntry(g1_nxm[i], Form("%d #times %d, k = %.1f #pm %.1e", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second, f1->GetParameter(1), f1->GetParError(1)), "p");
        leg2->AddEntry(g1_nxm[i], Form("%d #times %d", calo_window_dims_hcal_geom[i].first, calo_window_dims_hcal_geom[i].second), "p");
    }
    leg2->Draw("same");

    tx=0.45;
    ty=0.25;
    std::vector<std::string> sTags = {sPHENIX_Tag};
        for ( unsigned int j = 0; j < sTags.size(); ++j ) {
            tex->DrawLatex(tx, ty, sTags[j].c_str());
            ty -= 0.05;
        }
    c2->SaveAs((global_plots+"/avg_et_window_ratio.png").c_str());
    
    // TLegend * leg = new TLegend(0.6,0.7,0.9,0.9);
//     // leg->SetBorderSize(0);
//     // leg->SetFillStyle(0);

//     TLatex * tex = new TLatex();
//     tex->SetNDC();
//     tex->SetTextFont(42);
//     tex->SetTextSize(0.04);

//     // TLatex * tex2 = new TLatex();
//     // tex2->SetNDC();
//     // tex2->SetTextFont(42);
//     // tex2->SetTextSize(0.035);

//     gStyle->SetOptStat(0);
//     gStyle->SetOptFit(0);
//     gStyle->SetOptTitle(0);

//     double tx=0.45;
//     double ty=0.25;
//     std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
//     c->cd(1);
//     h2_et_rho->GetXaxis()->SetNdivisions(505);
//     TF1 * f1 = new TF1("f1", "[0]+[1]*x", 200, 1700);
//     h2_et_rho->Fit("f1", "Q", "", 200, 1700);
//     // change fit line color
//    // get the fit function from the histogram
    
//     TF1 *fitFunc = h2_et_rho->GetFunction("f1");
//     fitFunc->SetLineColor(kBlack);
//     fitFunc->SetLineWidth(3);
//     fitFunc->SetLineStyle(2);
//     h2_et_rho->Draw("colz");
//     // f1->Draw("same");
    
//     for ( unsigned int i = 0; i < tags.size(); ++i ) {
//         tex->DrawLatex(tx, ty, tags[i].c_str());
//         ty -= 0.05;
//     }
//     // use 1 siffig
//     tex2->DrawLatex(0.18, 0.85, Form("Fit: (%.1f #pm %.2f) + (%.3f #pm %.1e) #times E_{T}^{raw}", f1->GetParameter(0), f1->GetParError(0), f1->GetParameter(1), f1->GetParError(1)));
    
    
//     c->cd(2);
//     h2_et_rho_mult->GetXaxis()->SetNdivisions(505);
//     h2_et_rho_mult->Draw("colz");
//     h2_et_rho_mult->Fit("f1", "Q", "", 200, 1700);
//     fitFunc = h2_et_rho_mult->GetFunction("f1");
//     fitFunc->SetLineColor(kBlack);
//     fitFunc->SetLineWidth(3);
//     fitFunc->SetLineStyle(2);

//     ty = 0.25;
//     for ( unsigned int i = 0; i < tags.size(); ++i ) {
//         tex->DrawLatex(tx, ty, tags[i].c_str());
//         ty -= 0.05;
//     }
//     tex2->DrawLatex(0.18, 0.85, Form("Fit: (%.1f #pm %.2f) + (%.3f #pm %.1e) #times E_{T}^{raw}", f1->GetParameter(0), f1->GetParError(0), f1->GetParameter(1), f1->GetParError(1)));
//     c->SaveAs((global_plots+"/et_rho.png").c_str());    

//     // TH2F* h2_et_rho = new TH2F("h2_et_rho", "h2_et_rho", NET_BINS, et_bins, NRHOA_BINS, rho_binsa);
//     // h2_et_rho->GetXaxis()->SetTitle("#Sigma E_{T}^{raw} [GeV]");
//     // h2_et_rho->GetYaxis()->SetTitle("#rho_{A}#times A_{cone} [GeV]");
//     // TH2F* h2_avg_et_window = new TH2F("h2_avg_et_window", "h2_avg_et_window",  NET_BINS, et_bins, AVG_ET_BINS, avgbins);
//     // h2_et_rho_mult->GetXaxis()->SetTitle("#Sigma E_{T}^{raw} [GeV]");
//     // h2_et_rho_mult->GetYaxis()->SetTitle("#LT E_{T} #GT [GeV]");

//     // // auto window_dim = GetWindowDimFromString(h3s.at(0)->GetName());
//     // // std::vector<std::string> sTags = {sPHENIX_Tag};
//     // // sTags.push_back(Form("%d #times %d Towers", window_dim.first, window_dim.second));


//     // float ntot_cemc= 256.0*96.0;
//     // float not_hcal = 64.0*24.0;
//     // for ( int i = 0; i < nentries; ++i ) {
//     //     t->GetEntry(i);
    //     // float et = (ntowers_cemc-ndead_cemc)*ntot_cemc + ntowers_ihcal*not_hcal + ntowers_ohcal*not_hcal;
    //     float et = cemc_e + ihcal_e + ohcal_e;
    //     float rho_a = area_rho*TMath::Pi()*0.4*0.4;
    //     float rho_mm = rho_a/1005.3;
    //     float rho_m = rho_mm*rc_m;
    //     h2_et_rho->Fill(et, rho_a);
    //     h2_et_rho_mult->Fill(et, rho_m);

    // }
    // // h2_et_rho->Scale(1.0/h2_et_rho->Integral());
    // // h2_et_rho_mult->Scale(1.0/h2_et_rho_mult->Integral());

//     TCanvas *c = new TCanvas("c", "c", 1200,600);
//     c->Divide(2,1);
//     for (int i = 1; i <= 2; i++) {
//         c->cd(i);
//         // gPad->SetLogy();
//         // gPad->SetLogx();
//         gPad->SetLogz();
//         gPad->SetLeftMargin(0.15);
//         gPad->SetRightMargin(0.15);
//         gPad->SetBottomMargin(0.15);
//         gPad->SetTopMargin(0.1);
//     }

//     // TLegend * leg = new TLegend(0.6,0.7,0.9,0.9);
//     // leg->SetBorderSize(0);
//     // leg->SetFillStyle(0);

//     TLatex * tex = new TLatex();
//     tex->SetNDC();
//     tex->SetTextFont(42);
//     tex->SetTextSize(0.04);

//     TLatex * tex2 = new TLatex();
//     tex2->SetNDC();
//     tex2->SetTextFont(42);
//     tex2->SetTextSize(0.035);

//     gStyle->SetOptStat(0);
//     gStyle->SetOptFit(0);
//     gStyle->SetOptTitle(0);

//     double tx=0.45;
//     double ty=0.25;
//     std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
//     c->cd(1);
//     h2_et_rho->GetXaxis()->SetNdivisions(505);
//     TF1 * f1 = new TF1("f1", "[0]+[1]*x", 200, 1700);
//     h2_et_rho->Fit("f1", "Q", "", 200, 1700);
//     // change fit line color
//    // get the fit function from the histogram
    
//     TF1 *fitFunc = h2_et_rho->GetFunction("f1");
//     fitFunc->SetLineColor(kBlack);
//     fitFunc->SetLineWidth(3);
//     fitFunc->SetLineStyle(2);
//     h2_et_rho->Draw("colz");
//     // f1->Draw("same");
    
//     for ( unsigned int i = 0; i < tags.size(); ++i ) {
//         tex->DrawLatex(tx, ty, tags[i].c_str());
//         ty -= 0.05;
//     }
//     // use 1 siffig
//     tex2->DrawLatex(0.18, 0.85, Form("Fit: (%.1f #pm %.2f) + (%.3f #pm %.1e) #times E_{T}^{raw}", f1->GetParameter(0), f1->GetParError(0), f1->GetParameter(1), f1->GetParError(1)));
    
    
//     c->cd(2);
//     h2_et_rho_mult->GetXaxis()->SetNdivisions(505);
//     h2_et_rho_mult->Draw("colz");
//     h2_et_rho_mult->Fit("f1", "Q", "", 200, 1700);
//     fitFunc = h2_et_rho_mult->GetFunction("f1");
//     fitFunc->SetLineColor(kBlack);
//     fitFunc->SetLineWidth(3);
//     fitFunc->SetLineStyle(2);

//     ty = 0.25;
//     for ( unsigned int i = 0; i < tags.size(); ++i ) {
//         tex->DrawLatex(tx, ty, tags[i].c_str());
//         ty -= 0.05;
//     }
//     tex2->DrawLatex(0.18, 0.85, Form("Fit: (%.1f #pm %.2f) + (%.3f #pm %.1e) #times E_{T}^{raw}", f1->GetParameter(0), f1->GetParError(0), f1->GetParameter(1), f1->GetParError(1)));
//     c->SaveAs((global_plots+"/et_rho.png").c_str());    

    // float mbd_q_N = 0;
    // float mbd_q_S = 0;
    // float mbd_time_N = 0;
    // float mbd_time_S = 0;
    // t->SetBranchAddress("mbd_q_N", &mbd_q_N);
    // t->SetBranchAddress("mbd_q_S", &mbd_q_S);
    // t->SetBranchAddress("mbd_time_N", &mbd_time_N);
    // t->SetBranchAddress("mbd_time_S", &mbd_time_S);
   
    // float frac_fired[calo_nodes.size()];
    // float frac_dead[calo_nodes.size()];
    // float avg_et[calo_nodes.size()];
    // float std_et[calo_nodes.size()];
    // float sum_et[calo_nodes.size()];
    // for ( unsigned int i = 0; i < calo_nodes.size(); ++i ) {
    //     t->SetBranchAddress(Form("tower_frac_fired_%s", calo_nodes[i].c_str()), &frac_fired[i]);
    //     t->SetBranchAddress(Form("tower_frac_dead_%s", calo_nodes[i].c_str()), &frac_dead[i]);
    //     t->SetBranchAddress(Form("tower_avg_energy_%s", calo_nodes[i].c_str()), &avg_et[i]);
    //     t->SetBranchAddress(Form("tower_std_energy_%s", calo_nodes[i].c_str()), &std_et[i]);
    //     t->SetBranchAddress(Form("tower_sum_energy_%s", calo_nodes[i].c_str()), &sum_et[i]);
    // }
   
}


// void ProcessCaloSums(std::vector<std::string> calo_nodes , std::string calo_group )
// {

//     std::cout <<"Processing Calo Sums" << std::endl;
//     const int n_calo_nodes = calo_nodes.size();
//     std::vector<float> n_towers{};
//     std::vector<float> deta_towers{};
//     std::vector<float> dphi_towers{};
//     for ( unsigned int i = 0; i < n_calo_nodes; ++i ) {
//         bool is_cemc = calo_nodes[i].find("CEMC") != std::string::npos;
//         bool is_retower = calo_nodes[i].find("RETOWER") != std::string::npos;
//         bool cemc_geo = is_cemc && !is_retower;
//         if ( cemc_geo ){
//             n_towers.push_back(256.0*96.0);
//             deta_towers.push_back(2.2/96.0);
//             dphi_towers.push_back(2.0*TMath::Pi()/256.0);
//         } else {
//             n_towers.push_back(24.0*64.0);
//             deta_towers.push_back(2.2/24.0);
//             dphi_towers.push_back(2.0*TMath::Pi()/64.0);
//         }
//     }

//     for ( unsigned int i = 0; i < n_calo_nodes; ++i ) {
//         std::cout << "Calo Node: " << calo_nodes[i] << " n_towers: " << n_towers[i] << " deta: " << deta_towers[i] << " dphi: " << dphi_towers[i] << std::endl;
//     }

//     TTree *t = (TTree*)f->Get("T");
//     if(!t){ std::cout << "Tree T not found in file" << std::endl; return -1; }

//     // Centrality
//     int centrality = 0;
//     t->SetBranchAddress("centrality", &centrality);
    
//     // MBD
//     float mbd_q_N = 0;
//     float mbd_q_S = 0;
//     float mbd_time_N = 0;
//     float mbd_time_S = 0;
//     t->SetBranchAddress("mbd_q_N", &mbd_q_N);
//     t->SetBranchAddress("mbd_q_S", &mbd_q_S);
//     t->SetBranchAddress("mbd_time_N", &mbd_time_N);
//     t->SetBranchAddress("mbd_time_S", &mbd_time_S);
   
//     float frac_fired[calo_nodes.size()];
//     float frac_dead[calo_nodes.size()];
//     float avg_et[calo_nodes.size()];
//     float std_et[calo_nodes.size()];
//     float sum_et[calo_nodes.size()];
//     for ( unsigned int i = 0; i < calo_nodes.size(); ++i ) {
//         t->SetBranchAddress(Form("tower_frac_fired_%s", calo_nodes[i].c_str()), &frac_fired[i]);
//         t->SetBranchAddress(Form("tower_frac_dead_%s", calo_nodes[i].c_str()), &frac_dead[i]);
//         t->SetBranchAddress(Form("tower_avg_energy_%s", calo_nodes[i].c_str()), &avg_et[i]);
//         t->SetBranchAddress(Form("tower_std_energy_%s", calo_nodes[i].c_str()), &std_et[i]);
//         t->SetBranchAddress(Form("tower_sum_energy_%s", calo_nodes[i].c_str()), &sum_et[i]);
//     }

//     int nentries = t->GetEntries();
//     std::cout << "Number of entries = " << nentries << std::endl;

//     float mbd_q_min = 1e9, mbd_q_max = -1e9;
//     float mbd_Qmin = 1e9, mbd_Qmax = -1e9;
//     float mbd_dt_min = 1e9, mbd_dt_max = -1e9;
    
//     float sumfired_min = 1e9, sumfired_max = -1e9;
//     float fired_min = 1e9, fired_max = -1e9;
    
//     float dead_min = 1e9, dead_max = -1e9;
    
//     float avg_et_min = 1e9, avg_et_max = -1e9;
//     float tot_avg_et_min = 1e9, tot_avg_et_max = -1e9;

//     float std_et_min = 1e9, std_et_max = -1e9;
//     float tot_std_et_min = 1e9, tot_std_et_max = -1e9;

//     float sum_et_min = 1e9, sum_et_max = -1e9;
//     float tot_sum_et_min = 1e9, tot_sum_et_max = -1e9;
//     for ( int i = 0; i < nentries; ++i ) {
//         t->GetEntry(i);
//         float Qsum = mbd_q_N + mbd_q_S;
//         float dT = std::fabs(mbd_time_N - mbd_time_S);
//         mbd_q_min = std::min(mbd_q_min, std::min(mbd_q_N, mbd_q_S));
//         mbd_q_max = std::max(mbd_q_max, std::max(mbd_q_N, mbd_q_S));
//         mbd_Qmin = std::min(mbd_Qmin, Qsum);
//         mbd_Qmax = std::max(mbd_Qmax, Qsum);
//         mbd_dt_min = std::min(mbd_dt_min, dT);
//         mbd_dt_max = std::max(mbd_dt_max, dT);

//         float sumfired = 0;
//         float sumdead = 0;
//         float sumet = 0;
//         float sumet2 = 0;
    
//         for ( unsigned int j = 0; j < calo_nodes.size(); ++j ) {
//             sumfired += frac_fired[j] * n_towers[j];
//             fired_min = std::min(fired_min, frac_fired[j] * n_towers[j]);
//             fired_max = std::max(fired_max, frac_fired[j] * n_towers[j]);
//             sumdead += frac_dead[j] * n_towers[j];
//             dead_min = std::min(dead_min, frac_dead[j] * n_towers[j]);
//             dead_max = std::max(dead_max, frac_dead[j] * n_towers[j]);
//             sumet += sum_et[j];
//             sumet2 += ((std_et[j]*std_et[j]) + (avg_et[j]*avg_et[j])) * n_towers[j];
//             avg_et_min = std::min(avg_et_min, avg_et[j]);
//             avg_et_max = std::max(avg_et_max, avg_et[j]);
//             std_et_min = std::min(std_et_min, std_et[j]);
//             std_et_max = std::max(std_et_max, std_et[j]);
//             sum_et_min = std::min(sum_et_min, sum_et[j]);
//             sum_et_max = std::max(sum_et_max, sum_et[j]);
//         }
//         sumfired_min = std::min(sumfired_min, sumfired);
//         sumfired_max = std::max(sumfired_max, sumfired);
//         float tot_avg_et = sumet / sumfired;
//         tot_avg_et_min = std::min(tot_avg_et_min, tot_avg_et);
//         tot_avg_et_max = std::max(tot_avg_et_max, tot_avg_et);
//         float tot_std_et = std::sqrt(sumet2/sumfired - (tot_avg_et*tot_avg_et));
//         tot_std_et_min = std::min(tot_std_et_min, tot_std_et);
//         tot_std_et_max = std::max(tot_std_et_max, tot_std_et);
//         float tot_sum_et = sumet;
//         tot_sum_et_min = std::min(tot_sum_et_min, tot_sum_et);
//         tot_sum_et_max = std::max(tot_sum_et_max, tot_sum_et);

//     }

// //  float mbd_q_min = 1e9, mbd_q_max = -1e9;
// //     float mbd_Qmin = 1e9, mbd_Qmax = -1e9;
// //     float mbd_dt_min = 1e9, mbd_dt_max = -1e9;
    
// //     float sumfired_min = 1e9, sumfired_max = -1e9;
// //     float fired_min = 1e9, fired_max = -1e9;
    
// //     float dead_min = 1e9, dead_max = -1e9;
    
// //     float avg_et_min = 1e9, avg_et_max = -1e9;
// //     float tot_avg_et_min = 1e9, tot_avg_et_max = -1e9;

// //     float std_et_min = 1e9, std_et_max = -1e9;
// //     float tot_std_et_min = 1e9, tot_std_et_max = -1e9;

//     mbd_q_min = 0;
//     mbd_q_max = 2000;
//     mbd_Qmin = 0;
//     mbd_Qmax = 2000;
//     mbd_dt_min = 0;
//     mbd_dt_max = 64;
//     sumfired_min = 0;
//     sumfired_max = 12000;
//     dead_min = 0;
//     dead_max = 2000;
//     avg_et_min = 0;
//     avg_et_max = 0.3;
//     tot_avg_et_min = 0;
//     tot_avg_et_max = 0.3;
//     std_et_min = 0;
//     std_et_max = 0.3;
//     tot_std_et_min = 0;
//     tot_std_et_max = 0.3;
//     fired_min = 0;
//     fired_max = 12000;
//     sum_et_min = 0;
//     sum_et_max = 2000;
//     tot_sum_et_min = 0;
//     tot_sum_et_max = 2000;
//     // Make 1D histograms
//     const int N_MBD_Q = 500;
//     // mbd_Qmin=0;
//     TH1F * h1_mbd_sumQ = CreateHisto1D("h1_mbd_sumQ", N_MBD_Q, mbd_Qmin, mbd_Qmax, QuickXYTitles("#Sigma Q_{MBD}", true));
//     TH1F * h1_mbd_q_N = CreateHisto1D("h1_mbd_q_N", N_MBD_Q, mbd_q_min, mbd_q_max, QuickXYTitles("Q_{MBD}^{North}", true));
//     TH1F * h1_mbd_q_S = CreateHisto1D("h1_mbd_q_S", N_MBD_Q, mbd_q_min, mbd_q_max, QuickXYTitles("Q_{MBD}^{South}", true));
    
//     const int N_MBD_DT = 100;
//     mbd_dt_min=0;
//     TH1F * h1_mbd_dt = CreateHisto1D("h1_mbd_dt", N_MBD_DT, mbd_dt_min, mbd_dt_max   , QuickXYTitles("#Delta t_{MBD}", true, "ns"));

//     const int N_CENT = 20;
//     const float MIN_CENT = 0.0;
//     const float MAX_CENT = 90.0;
//     TH1F * h1_cent = CreateHisto1D("h1_cent", N_CENT, MIN_CENT, MAX_CENT, QuickXYTitles("Centrality [%]", false));

//     const int N_SUMFIRED = int((sumfired_max-sumfired_min)/10);
//     TH1F * h1_total_fired = CreateHisto1D("h1_total_fired", N_SUMFIRED, sumfired_min, sumfired_max, QuickXYTitles("N_{towers}^{raw}", true));
//     h1_total_fired->SetTitle(calo_group.c_str());
//     const int N_FIRED = int((fired_max-fired_min)/10);
//     TH1F * h1_fired[calo_nodes.size()];
//     for ( unsigned int i = 0; i < calo_nodes.size(); ++i ) {
//         h1_fired[i] = CreateHisto1D(Form("h1_fired_%s", calo_nodes[i].c_str()), N_FIRED, fired_min, fired_max, QuickXYTitles("N_{towers}^{raw}", true));
//         h1_fired[i]->SetTitle(GetCaloNameFromNode(calo_nodes[i]).c_str());
//     }
//     dead_min=0;
//     const int N_DEAD = int((dead_max-dead_min)/10);
//     TH1F * h1_total_dead = CreateHisto1D("h1_total_dead", N_DEAD, dead_min, dead_max, QuickXYTitles("N_{towers}^{dead}", true));

//     const int N_AVG_ET = 300;
//     TH1F * h1_total_avg_et = CreateHisto1D("h1_total_avg_et", N_AVG_ET, tot_avg_et_min, tot_avg_et_max, QuickXYTitles("#LT E_{T}^{raw} #GT", true, "GeV"));
//     h1_total_avg_et->SetTitle(calo_group.c_str());
//     TH1F * h1_avg_et[calo_nodes.size()];
//     for ( unsigned int i = 0; i < calo_nodes.size(); ++i ) {
//         h1_avg_et[i] = CreateHisto1D(Form("h1_avg_et_%s", calo_nodes[i].c_str()), N_AVG_ET, avg_et_min, avg_et_max, QuickXYTitles("#LT E_{T}^{raw} #GT", true, "GeV"));
//         h1_avg_et[i]->SetTitle(GetCaloNameFromNode(calo_nodes[i]).c_str());
//     }

//     const int N_STD_ET = 300;
//     TH1F * h1_total_std_et = CreateHisto1D("h1_total_std_et", N_STD_ET, tot_std_et_min, tot_std_et_max, QuickXYTitles("#sigma E_{T}^{raw}", true, "GeV"));
//     h1_total_std_et->SetTitle(calo_group.c_str());
//     TH1F * h1_std_et[calo_nodes.size()];
//     for ( unsigned int i = 0; i < calo_nodes.size(); ++i ) {
//         h1_std_et[i] = CreateHisto1D(Form("h1_std_et_%s", calo_nodes[i].c_str()), N_STD_ET, std_et_min, std_et_max, QuickXYTitles("#sigma E_{T}^{raw}", true, "GeV"));
//         h1_std_et[i]->SetTitle(GetCaloNameFromNode(calo_nodes[i]).c_str());
//     }

//     const int N_SUM_ET = 500;
//     TH1F * h1_total_sum_et = CreateHisto1D("h1_total_sum_et", N_SUM_ET, tot_sum_et_min, tot_sum_et_max, QuickXYTitles("#Sigma E_{T}^{raw}", true, "GeV"));
//     h1_total_sum_et->SetTitle(calo_group.c_str());
//     TH1F * h1_sum_et[calo_nodes.size()];
//     for ( unsigned int i = 0; i < calo_nodes.size(); ++i ) {
//         h1_sum_et[i] = CreateHisto1D(Form("h1_sum_et_%s", calo_nodes[i].c_str()), N_SUM_ET, sum_et_min, sum_et_max, QuickXYTitles("#Sigma E_{T}^{raw}", true, "GeV"));
//         h1_sum_et[i]->SetTitle(GetCaloNameFromNode(calo_nodes[i]).c_str());
//     }

//     // Make 2D histograms
//     // x-axes: cent, tot_sum_et, mbd_sumQ
//     float cent_bins[N_CENT];
//     QuickBins(N_CENT, MIN_CENT, MAX_CENT, cent_bins);
//     float tot_et_bins[N_SUM_ET];
//     QuickBins(N_SUM_ET, tot_sum_et_min, tot_sum_et_max, tot_et_bins);
//     float mbd_sumq_bins[N_MBD_Q];
//     QuickBins(N_MBD_Q, mbd_Qmin, mbd_Qmax, mbd_sumq_bins);

//     // all histos will have a y for sum, and one for each calo node
//     float tot_fired_bins[N_SUMFIRED];
//     float fired_bins[N_FIRED];
//     QuickBins(N_SUMFIRED, sumfired_min, sumfired_max, tot_fired_bins);
//     QuickBins(N_FIRED, fired_min, fired_max, fired_bins);

//     float tot_avg_et_bins[N_AVG_ET];
//     float avg_et_bins[N_AVG_ET];
//     QuickBins(N_AVG_ET, tot_avg_et_min, tot_avg_et_max, tot_avg_et_bins);
//     QuickBins(N_AVG_ET, avg_et_min, avg_et_max, avg_et_bins);

//     float tot_std_et_bins[N_STD_ET];
//     float std_et_bins[N_STD_ET];
//     QuickBins(N_STD_ET, tot_std_et_min, tot_std_et_max, tot_std_et_bins);
//     QuickBins(N_STD_ET, std_et_min, std_et_max, std_et_bins);

//     // totals 
//     TH2F * h2_tot_fired_vs_cent = new TH2F("h2_tot_fired_vs_cent", 
//                 "h2_tot_fired_vs_cent;Centrality [%];N_{towers}^{raw}",
//                 N_CENT, cent_bins, N_SUMFIRED, tot_fired_bins);
//     TH2F * h2_tot_fired_vs_sumet = new TH2F("h2_tot_fired_vs_sumet",
//                 "h2_tot_fired_vs_sumet;#Sigma E_{T}^{raw} [GeV];N_{towers}^{raw}",
//                 N_SUM_ET, tot_et_bins, N_SUMFIRED, tot_fired_bins);
//     TH2F * h2_tot_fired_vs_mbdQ = new TH2F("h2_tot_fired_vs_mbdQ",
//                 "h2_tot_fired_vs_mbdQ;#Sigma Q_{MBD};N_{towers}^{raw}",
//                 N_MBD_Q, mbd_sumq_bins, N_SUMFIRED, tot_fired_bins);
//     TH2F * h2_tot_avg_et_vs_cent = new TH2F("h2_tot_avg_et_vs_cent",
//                 "h2_tot_avg_et_vs_cent;Centrality [%];#LT E_{T}^{raw} #GT",
//                 N_CENT, cent_bins, N_AVG_ET, tot_avg_et_bins);
//     TH2F * h2_tot_avg_et_vs_sumet = new TH2F("h2_tot_avg_et_vs_sumet", 
//                 "h2_tot_avg_et_vs_sumet;#Sigma E_{T}^{raw} [GeV];#LT E_{T}^{raw} #GT",
//                 N_SUM_ET, tot_et_bins, N_AVG_ET, tot_avg_et_bins);
//     TH2F * h2_tot_avg_et_vs_mbdQ = new TH2F("h2_tot_avg_et_vs_mbdQ",
//                 "h2_tot_avg_et_vs_mbdQ;#Sigma Q_{MBD};#LT E_{T}^{raw} #GT",
//                 N_MBD_Q, mbd_sumq_bins, N_AVG_ET, tot_avg_et_bins);
//     TH2F * h2_tot_std_et_vs_cent = new TH2F("h2_tot_std_et_vs_cent",
//                 "h2_tot_std_et_vs_cent;Centrality [%];#sigma E_{T}^{raw}",
//                 N_CENT, cent_bins, N_STD_ET, tot_std_et_bins);
//     TH2F * h2_tot_std_et_vs_sumet = new TH2F("h2_tot_std_et_vs_sumet",
//                 "h2_tot_std_et_vs_sumet;#Sigma E_{T}^{raw} [GeV];#sigma E_{T}^{raw}",
//                 N_SUM_ET, tot_et_bins, N_STD_ET, tot_std_et_bins);
//     TH2F * h2_tot_std_et_vs_mbdQ = new TH2F("h2_tot_std_et_vs_mbdQ",
//                 "h2_tot_std_et_vs_mbdQ;#Sigma Q_{MBD};#sigma E_{T}^{raw}",
//                 N_MBD_Q, mbd_sumq_bins, N_STD_ET, tot_std_et_bins); 
//     TH2F * h2_tot_sum_et_vs_cent = new TH2F("h2_tot_sum_et_vs_cent", 
//                 "h2_tot_sum_et_vs_cent;Centrality [%];#Sigma E_{T}^{raw}",
//                 N_CENT, cent_bins, N_SUM_ET, tot_et_bins);
//     TH2F * h2_tot_sum_et_vs_mbdQ = new TH2F("h2_tot_sum_et_vs_mbdQ", 
//                 "h2_tot_sum_et_vs_mbdQ;#Sigma Q_{MBD};#Sigma E_{T}^{raw}",
//                 N_MBD_Q, mbd_sumq_bins, N_SUM_ET, tot_et_bins);
//     TH2F * h2_cent_vs_mbdQ = new TH2F("h2_cent_vs_mbdQ",
//                 "h2_cent_vs_mbdQ;#Sigma Q_{MBD};Centrality [%]",
//                 N_MBD_Q, mbd_sumq_bins, N_CENT, cent_bins);
    
//     // individual calorimeters
//     TH2F * h2_fired_vs_cent[calo_nodes.size()];
//     TH2F * h2_fired_vs_sumet[calo_nodes.size()];
//     TH2F * h2_fired_vs_mbdQ[calo_nodes.size()];
//     TH2F * h2_avg_et_vs_cent[calo_nodes.size()];
//     TH2F * h2_avg_et_vs_sumet[calo_nodes.size()];
//     TH2F * h2_avg_et_vs_mbdQ[calo_nodes.size()];
//     TH2F * h2_std_et_vs_cent[calo_nodes.size()];
//     TH2F * h2_std_et_vs_sumet[calo_nodes.size()];
//     TH2F * h2_std_et_vs_mbdQ[calo_nodes.size()];
//     TH2F * h2_sum_et_vs_cent[calo_nodes.size()];
//     TH2F * h2_sum_et_vs_mbdQ[calo_nodes.size()];
//     for ( unsigned int i = 0; i < calo_nodes.size(); ++i ) {
//         h2_fired_vs_cent[i] = new TH2F(Form("h2_fired_vs_cent_%s", calo_nodes[i].c_str()),
//              "Centrality [%];N_{towers}^{raw}", 
//                 N_CENT, cent_bins, N_FIRED, fired_bins);
//         h2_fired_vs_sumet[i] = new TH2F(Form("h2_fired_vs_sumet_%s", calo_nodes[i].c_str()),
//                 "#Sigma E_{T}^{raw} [GeV];N_{towers}^{raw}", 
//                 N_SUM_ET, tot_et_bins, N_FIRED, fired_bins);
//         h2_fired_vs_mbdQ[i] = new TH2F(Form("h2_fired_vs_mbdQ_%s", calo_nodes[i].c_str()),
//                "#Sigma Q_{MBD};N_{towers}^{raw}",
//                 N_MBD_Q, mbd_sumq_bins, N_FIRED, fired_bins);
//         h2_avg_et_vs_cent[i] = new TH2F(Form("h2_avg_et_vs_cent_%s", calo_nodes[i].c_str()),
//                 "Centrality [\%];#LT E_{T}^{raw} #GT", 
//                 N_CENT, cent_bins, N_AVG_ET, avg_et_bins);
//         h2_avg_et_vs_sumet[i] = new TH2F(Form("h2_avg_et_vs_sumet_%s", calo_nodes[i].c_str()),
//               "#Sigma E_{T}^{raw} [GeV];#LT E_{T}^{raw} #GT", 
//                 N_SUM_ET, tot_et_bins, N_AVG_ET, avg_et_bins);
//         h2_avg_et_vs_mbdQ[i] = new TH2F(Form("h2_avg_et_vs_mbdQ_%s", calo_nodes[i].c_str()),
//              "#Sigma Q_{MBD};#LT E_{T}^{raw} #GT", 
//                 N_MBD_Q, mbd_sumq_bins, N_AVG_ET, avg_et_bins);
//         h2_std_et_vs_cent[i] = new TH2F(Form("h2_std_et_vs_cent_%s", calo_nodes[i].c_str()),
//              "Centrality [\%];#sigma E_{T}^{raw}", 
//                 N_CENT, cent_bins, N_STD_ET, std_et_bins);
//         h2_std_et_vs_sumet[i] = new TH2F(Form("h2_std_et_vs_sumet_%s", calo_nodes[i].c_str()),
//                "#Sigma E_{T}^{raw} [GeV];#sigma E_{T}^{raw}", 
//                 N_SUM_ET, tot_et_bins, N_STD_ET, std_et_bins);
//         h2_std_et_vs_mbdQ[i] = new TH2F(Form("h2_std_et_vs_mbdQ_%s", calo_nodes[i].c_str()),
//                 "#Sigma Q_{MBD};#sigma E_{T}^{raw}", 
//                 N_MBD_Q, mbd_sumq_bins, N_STD_ET, std_et_bins); 
//         h2_sum_et_vs_cent[i] = new TH2F(Form("h2_sum_et_vs_cent_%s", calo_nodes[i].c_str()),
//              "Centrality [%];#Sigma E_{T}^{raw}", 
//                 N_CENT, cent_bins, N_SUM_ET, tot_et_bins);
//         h2_sum_et_vs_mbdQ[i] = new TH2F(Form("h2_sum_et_vs_mbdQ_%s", calo_nodes[i].c_str()),
//           "#Sigma Q_{MBD};#Sigma E_{T}^{raw}", 
//                 N_MBD_Q, mbd_sumq_bins, N_SUM_ET, tot_et_bins);
//     }
//     // Make 1D histograms
//     for ( int i = 0; i < nentries; ++i ) {
//         t->GetEntry(i);

//         float Qsum = mbd_q_N + mbd_q_S;
//         float dT = std::fabs(mbd_time_N - mbd_time_S);
//         float tot_fired = 0;
//         float tot_et = 0;
//         float tot_dead = 0;
//         float tot_et2 = 0;
//         float tot_std_et = 0;
//         float tot_avg_et = 0;
//         for ( unsigned int j = 0; j < calo_nodes.size(); ++j ) {
//             tot_fired += frac_fired[j] * n_towers[j];
//             tot_et += sum_et[j];
//             tot_et2 += ((std_et[j]*std_et[j]) + (avg_et[j]*avg_et[j])) * n_towers[j];
//             tot_dead += frac_dead[j] * n_towers[j];
//             h1_fired[j]->Fill(frac_fired[j] * n_towers[j]);
//             h1_avg_et[j]->Fill(avg_et[j]);
//             h1_std_et[j]->Fill(std_et[j]);
//             h1_sum_et[j]->Fill(sum_et[j]);
//         }
//         tot_avg_et = tot_et / tot_fired;
//         tot_std_et = std::sqrt(tot_et2/tot_fired - (tot_avg_et*tot_avg_et));

//         h1_mbd_sumQ->Fill(Qsum);
//         h1_mbd_q_N->Fill(mbd_q_N);
//         h1_mbd_q_S->Fill(mbd_q_S);
//         h1_mbd_dt->Fill(dT);
//         h1_cent->Fill(centrality);
//         h1_total_fired->Fill(tot_fired);
//         h1_total_avg_et->Fill(tot_avg_et);
//         h1_total_std_et->Fill(tot_std_et);
//         h1_total_sum_et->Fill(tot_et);
//         for ( unsigned int j = 0; j < calo_nodes.size(); ++j ) {
//             h2_fired_vs_cent[j]->Fill(centrality, frac_fired[j] * n_towers[j]);
//             h2_fired_vs_sumet[j]->Fill(tot_et, frac_fired[j] * n_towers[j]);
//             h2_fired_vs_mbdQ[j]->Fill(Qsum, frac_fired[j] * n_towers[j]);
//             h2_avg_et_vs_cent[j]->Fill(centrality, avg_et[j]);
//             h2_avg_et_vs_sumet[j]->Fill(tot_et, avg_et[j]);
//             h2_avg_et_vs_mbdQ[j]->Fill(Qsum, avg_et[j]);
//             h2_std_et_vs_cent[j]->Fill(centrality, std_et[j]);
//             h2_std_et_vs_sumet[j]->Fill(tot_et, std_et[j]);
//             h2_std_et_vs_mbdQ[j]->Fill(Qsum, std_et[j]);
//             h2_sum_et_vs_cent[j]->Fill(centrality, sum_et[j]);
//             h2_sum_et_vs_mbdQ[j]->Fill(Qsum, sum_et[j]);
//         }
//         h2_tot_fired_vs_cent->Fill(centrality, tot_fired);
//         h2_tot_fired_vs_sumet->Fill(tot_et, tot_fired);
//         h2_tot_fired_vs_mbdQ->Fill(Qsum, tot_fired);
//         h2_tot_avg_et_vs_cent->Fill(centrality, tot_avg_et);
//         h2_tot_avg_et_vs_sumet->Fill(tot_et, tot_avg_et);
//         h2_tot_avg_et_vs_mbdQ->Fill(Qsum, tot_avg_et);
//         h2_tot_std_et_vs_cent->Fill(centrality, tot_std_et);
//         h2_tot_std_et_vs_sumet->Fill(tot_et, tot_std_et);
//         h2_tot_std_et_vs_mbdQ->Fill(Qsum, tot_std_et);
//         h2_tot_sum_et_vs_cent->Fill(centrality, tot_et);
//         h2_tot_sum_et_vs_mbdQ->Fill(Qsum, tot_et);
//         h2_cent_vs_mbdQ->Fill(Qsum, centrality);

//     }

//     // Write histograms to file

//     std::cout << "Writing histograms to file" << std::endl;
//     std::string output_file = global_plots + "/CaloSumsGlobal-"+calo_group+"-HISTOGRAMS.root";
//     TFile *fout = new TFile(output_file.c_str(), "RECREATE");
//     fout->cd();
//     std::vector<TH1F*> h1s = {h1_mbd_sumQ, h1_mbd_q_N, h1_mbd_q_S, h1_mbd_dt, h1_cent, h1_total_fired, h1_total_avg_et, 
//         h1_total_std_et, h1_total_sum_et};
//     for ( unsigned int i = 0; i < calo_nodes.size(); ++i ) {
//         h1s.push_back(h1_fired[i]);
//         h1s.push_back(h1_avg_et[i]);
//         h1s.push_back(h1_std_et[i]);
//         h1s.push_back(h1_sum_et[i]);
//     }

//     for ( unsigned int i = 0; i < h1s.size(); ++i ) {
//         h1s[i]->Scale(1.0/1.0*nentries);
//         h1s[i]->Write();
//     }   

//     std::vector<TH2F*> h2s = {h2_tot_fired_vs_cent, h2_tot_fired_vs_sumet, h2_tot_fired_vs_mbdQ, h2_tot_avg_et_vs_cent,
//         h2_tot_avg_et_vs_sumet, h2_tot_avg_et_vs_mbdQ, h2_tot_std_et_vs_cent, h2_tot_std_et_vs_sumet, h2_tot_std_et_vs_mbdQ,
//         h2_tot_sum_et_vs_cent, h2_tot_sum_et_vs_mbdQ, h2_cent_vs_mbdQ};
//     for ( unsigned int i = 0; i < calo_nodes.size(); ++i ) {
//         h2s.push_back(h2_fired_vs_cent[i]);
//         h2s.push_back(h2_fired_vs_sumet[i]);
//         h2s.push_back(h2_fired_vs_mbdQ[i]);
//         h2s.push_back(h2_avg_et_vs_cent[i]);
//         h2s.push_back(h2_avg_et_vs_sumet[i]);
//         h2s.push_back(h2_avg_et_vs_mbdQ[i]);
//         h2s.push_back(h2_std_et_vs_cent[i]);
//         h2s.push_back(h2_std_et_vs_sumet[i]);
//         h2s.push_back(h2_std_et_vs_mbdQ[i]);
//         h2s.push_back(h2_sum_et_vs_cent[i]);
//         h2s.push_back(h2_sum_et_vs_mbdQ[i]);
//     }
    
//     for ( unsigned int i = 0; i < h2s.size(); ++i ) {
//         // h2s[i]->Scale(1.0/1.0*nentries);
//         h2s[i]->Write();
//     }

//     fout->Close();






    // return;
// }

void PPG04_Histograms(const std::string & input_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/JAN30/DATA/BASIC/DATA-ProdA_2023-BASIC-023745.root") 
{

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    SetsPhenixStyle();
    gErrorIgnoreLevel = kWarning;

    IS_DATA= input_file.find("DATA") != std::string::npos;
    DataType_Tag = IS_DATA ? "Au+Au #sqrt{s_{NN}} = 200 GeV" : "HIJING MDC2";

    // get base name of input file
    std::string input_file_base = input_file;
    size_t found = input_file.find_last_of("/");
    if (found != std::string::npos) {
        input_file_base = input_file.substr(found+1);
    }
    // replace .root with -HISTOGRAMS.root
    found = input_file_base.find_last_of(".");
    if (found != std::string::npos) {
        input_file_base = input_file_base.substr(0, found);
    }

    std::cout << "Input file: " << input_file << std::endl;

    
    ConfigureOutputDirs(input_file_base);
    std::cout << "Global plots: " << global_plots << std::endl;
    std::cout << "Calo plots: " << calo_plots << std::endl;
    std::cout << "Calo window plots: " << calo_window_plots << std::endl;
    std::cout << "Background plots: " << bkgd_plots << std::endl;
    std::cout << "Random cone plots: " << random_cone_plots << std::endl;
   
    // open input file
    f = new TFile(input_file.c_str(), "READ");
    if(!f->IsOpen() || f->IsZombie()){ std::cout << "File " << input_file << " is zombie" << std::endl;  return -1; }

    // bool process_calo_sums = false;
    // if (process_calo_sums) {
    //     // Calo Sums
    //     // ProcessGlobal({"TOWERINFO_CALIB_CEMC", "TOWERINFO_CALIB_HCALIN", "TOWERINFO_CALIB_HCALOUT"}, "FullCalos");
    //     ProcessGlobal(input_file);
    // }
    ProcessGlobal();

    bool process_calo_window_histos = true;
    // if ( process_calo_window_histos) {  ProcessCaloWindowHistograms(); }

    f->Close();

    return 0;


   
}


// void PPG04_Histograms(const std::string & input_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/JAN28/DATA/BASIC/DATA-ProdA_2023-BASIC-023745.root") 
// {

//     TH1::SetDefaultSumw2();
//     TH2::SetDefaultSumw2();
//     TH3::SetDefaultSumw2();
//     SetsPhenixStyle();
//     gErrorIgnoreLevel = kWarning;

//     IS_DATA= input_file.find("DATA") != std::string::npos;
//     DataType_Tag = IS_DATA ? "Au+Au #sqrt{s_{NN}} = 200 GeV" : "HIJING MDC2";

//     // get base name of input file
//     std::string input_file_base = input_file;
//     size_t found = input_file.find_last_of("/");
//     if (found != std::string::npos) {
//         input_file_base = input_file.substr(found+1);
//     }
//     // replace .root with -HISTOGRAMS.root
//     found = input_file_base.find_last_of(".");
//     if (found != std::string::npos) {
//         input_file_base = input_file_base.substr(0, found);
//     }

//     std::cout << "Input file: " << input_file << std::endl;

    
//     ConfigureOutputDirs(input_file_base);
//     std::cout << "Global plots: " << global_plots << std::endl;
//     std::cout << "Calo plots: " << calo_plots << std::endl;
//     std::cout << "Calo window plots: " << calo_window_plots << std::endl;
//     std::cout << "Background plots: " << bkgd_plots << std::endl;
//     std::cout << "Random cone plots: " << random_cone_plots << std::endl;
   
//     // open input file
//     f = new TFile(input_file.c_str(), "READ");
//     if(!f->IsOpen() || f->IsZombie()){ std::cout << "File " << input_file << " is zombie" << std::endl;  return -1; }

//     // get tree
//     TTree *t = (TTree*)f->Get("T");
//     if(!t){ std::cout << "Tree T not found in " << input_file << std::endl; return -1; }

//     // -----------------------------------------------------
//     // Global variables
//     // -----------------------------------------------------

//     int event_id = 0;
//     unsigned int random_seed = 0;
//     t->SetBranchAddress("event_id", &event_id);
//     t->SetBranchAddress("random_seed", &random_seed);

//     // nevents histogram
//     TH1F * h1_num_events = (TH1F*)f->Get("h1_num_events");
//     if( !h1_num_events ){ std::cout << "h1_num_events not found!" << std::endl; return -1; }
//     NEVENTS = (int)h1_num_events->GetBinContent(1);
//     std::cout << "Number of events: " << NEVENTS << std::endl;

//     // MBD
//     float mbd_q_N = 0;
//     float mbd_q_S = 0;
//     float mbd_time_N = 0;
//     float mbd_time_S = 0;
//     float mbd_npmt_N = 0;
//     float mbd_npmt_S = 0;
//     t->SetBranchAddress("mbd_q_N", &mbd_q_N);
//     t->SetBranchAddress("mbd_q_S", &mbd_q_S);
//     t->SetBranchAddress("mbd_time_N", &mbd_time_N);
//     t->SetBranchAddress("mbd_time_S", &mbd_time_S);
//     t->SetBranchAddress("mbd_npmt_N", &mbd_npmt_N);
//     t->SetBranchAddress("mbd_npmt_S", &mbd_npmt_S);

//     // Vertex
//     float zvtx = 0;
//     t->SetBranchAddress("zvtx", &zvtx);
//     TH2F * h2_vtrx_cent = new TH2F("h2_vtx_cent", ";z_{vtx} (cm);Centrality (%)", 100, -50, 50, N_CENT_BINS, CENT_BINS);

  
//     // Centrality
//     int centrality = 0;
//     int impact_parameter = 0;
//     t->SetBranchAddress("centrality", &centrality);
//     t->SetBranchAddress("impact_parameter", &impact_parameter);

//     // -----------------------------------------------------
//     // Tower background variables
//     // -----------------------------------------------------
//     float tower_background_v2 = 0;
//     float tower_background_psi2 = 0;
//     std::vector< float > * tower_background_energy_recemc = 0;
//     std::vector< float > * tower_background_energy_hcalin = 0;
//     std::vector< float > * tower_background_energy_hcalout = 0;
//     t->SetBranchAddress("tower_background_v2", &tower_background_v2);
//     t->SetBranchAddress("tower_background_psi2", &tower_background_psi2);
//     t->SetBranchAddress("tower_background_energy_recemc", &tower_background_energy_recemc);
//     t->SetBranchAddress("tower_background_energy_hcalin", &tower_background_energy_hcalin);
//     t->SetBranchAddress("tower_background_energy_hcalout", &tower_background_energy_hcalout);


//     // -----------------------------------------------------
//     // Rho variables
//     // -----------------------------------------------------
//     // TowerRho_AREA
//     float rho_val_TowerRho_AREA = 0;
//     float std_rho_val_TowerRho_AREA = 0;
//     t->SetBranchAddress("rho_val_TowerRho_AREA", &rho_val_TowerRho_AREA);
//     t->SetBranchAddress("std_rho_val_TowerRho_AREA", &std_rho_val_TowerRho_AREA);
    
//     // TH2F * h2_rho_vs_sumEt

//     // TowerRho_MULT
//     float rho_val_TowerRho_MULT = 0;
//     float std_rho_val_TowerRho_MULT = 0;
//     t->SetBranchAddress("rho_val_TowerRho_MULT", &rho_val_TowerRho_MULT);
//     t->SetBranchAddress("std_rho_val_TowerRho_MULT", &std_rho_val_TowerRho_MULT);

//     // TowerRho_AREA_CEMC
//     float rho_val_TowerRho_AREA_CEMC = 0;
//     float std_rho_val_TowerRho_AREA_CEMC = 0;
//     t->SetBranchAddress("rho_val_TowerRho_AREA_CEMC", &rho_val_TowerRho_AREA_CEMC);
//     t->SetBranchAddress("std_rho_val_TowerRho_AREA_CEMC", &std_rho_val_TowerRho_AREA_CEMC);

//     // TowerRho_MULT_CEMC
//     float rho_val_TowerRho_MULT_CEMC = 0;
//     float std_rho_val_TowerRho_MULT_CEMC = 0;
//     t->SetBranchAddress("rho_val_TowerRho_MULT_CEMC", &rho_val_TowerRho_MULT_CEMC);
//     t->SetBranchAddress("std_rho_val_TowerRho_MULT_CEMC", &std_rho_val_TowerRho_MULT_CEMC);

//     // TowerRho_AREA_HCALIN
//     float rho_val_TowerRho_AREA_HCALIN = 0;
//     float std_rho_val_TowerRho_AREA_HCALIN = 0;
//     t->SetBranchAddress("rho_val_TowerRho_AREA_HCALIN", &rho_val_TowerRho_AREA_HCALIN);
//     t->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALIN", &std_rho_val_TowerRho_AREA_HCALIN);

//     // TowerRho_MULT_HCALIN
//     float rho_val_TowerRho_MULT_HCALIN = 0;
//     float std_rho_val_TowerRho_MULT_HCALIN = 0;
//     t->SetBranchAddress("rho_val_TowerRho_MULT_HCALIN", &rho_val_TowerRho_MULT_HCALIN);
//     t->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALIN", &std_rho_val_TowerRho_MULT_HCALIN);

//     // TowerRho_AREA_HCALOUT
//     float rho_val_TowerRho_AREA_HCALOUT = 0;
//     float std_rho_val_TowerRho_AREA_HCALOUT = 0;
//     t->SetBranchAddress("rho_val_TowerRho_AREA_HCALOUT", &rho_val_TowerRho_AREA_HCALOUT);
//     t->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALOUT", &std_rho_val_TowerRho_AREA_HCALOUT);

//     // TowerRho_MULT_HCALOUT
//     float rho_val_TowerRho_MULT_HCALOUT = 0;
//     float std_rho_val_TowerRho_MULT_HCALOUT = 0;
//     t->SetBranchAddress("rho_val_TowerRho_MULT_HCALOUT", &rho_val_TowerRho_MULT_HCALOUT);
//     t->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALOUT", &std_rho_val_TowerRho_MULT_HCALOUT);

//     // -----------------------------------------------------
//     // Calo tower variables
//     // -----------------------------------------------------
//     // TOWERINFO_CALIB_CEMC
//     float tower_frac_fired_TOWERINFO_CALIB_CEMC = 0;
//     float tower_frac_dead_TOWERINFO_CALIB_CEMC = 0;
//     float tower_avg_energy_TOWERINFO_CALIB_CEMC = 0;
//     float tower_std_energy_TOWERINFO_CALIB_CEMC = 0;
//     float tower_sum_energy_TOWERINFO_CALIB_CEMC = 0;
//     t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC", &tower_frac_fired_TOWERINFO_CALIB_CEMC);
//     t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC", &tower_frac_dead_TOWERINFO_CALIB_CEMC);
//     t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC", &tower_avg_energy_TOWERINFO_CALIB_CEMC);
//     t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC", &tower_std_energy_TOWERINFO_CALIB_CEMC);
//     t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC", &tower_sum_energy_TOWERINFO_CALIB_CEMC);

//     TH2F * h2_calotower_frac_fired_FULL_vs_cent = new TH2F("h2_calotower_frac_fired_FULL_vs_cent", ";Centrality (%);Fraction of fired towers", N_CENT, MIN_CENT ,MAX_CENT 100, 0, 1);

//     // TOWERINFO_CALIB_HCALIN
//     float tower_frac_fired_TOWERINFO_CALIB_HCALIN = 0;
//     float tower_frac_dead_TOWERINFO_CALIB_HCALIN = 0;
//     float tower_avg_energy_TOWERINFO_CALIB_HCALIN = 0;
//     float tower_std_energy_TOWERINFO_CALIB_HCALIN = 0;
//     float tower_sum_energy_TOWERINFO_CALIB_HCALIN = 0;
//     t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN", &tower_frac_fired_TOWERINFO_CALIB_HCALIN);
//     t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN", &tower_frac_dead_TOWERINFO_CALIB_HCALIN);
//     t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN", &tower_avg_energy_TOWERINFO_CALIB_HCALIN);
//     t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN", &tower_std_energy_TOWERINFO_CALIB_HCALIN);
//     t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN", &tower_sum_energy_TOWERINFO_CALIB_HCALIN);

//     // TOWERINFO_CALIB_HCALOUT
//     float tower_frac_fired_TOWERINFO_CALIB_HCALOUT = 0;
//     float tower_frac_dead_TOWERINFO_CALIB_HCALOUT = 0;
//     float tower_avg_energy_TOWERINFO_CALIB_HCALOUT = 0;
//     float tower_std_energy_TOWERINFO_CALIB_HCALOUT = 0;
//     float tower_sum_energy_TOWERINFO_CALIB_HCALOUT = 0;
//     t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT", &tower_frac_fired_TOWERINFO_CALIB_HCALOUT);
//     t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT", &tower_frac_dead_TOWERINFO_CALIB_HCALOUT);
//     t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT", &tower_avg_energy_TOWERINFO_CALIB_HCALOUT);
//     t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT", &tower_std_energy_TOWERINFO_CALIB_HCALOUT);
//     t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT", &tower_sum_energy_TOWERINFO_CALIB_HCALOUT);

//     // TOWERINFO_CALIB_CEMC_RETOWER
//     float tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER = 0;
//     float tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER = 0;
//     float tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
//     float tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
//     float tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER = 0;
//     t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER", &tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER);
//     t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER", &tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER);
//     t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER);
//     t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER);
//     t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER", &tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER);

//     // TOWERINFO_CALIB_CEMC_RETOWER_SUB1
//     float tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
//     float tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
//     float tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
//     float tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
//     float tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = 0;
//     t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
//     t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
//     t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
//     t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);
//     t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1);

//     // TOWERINFO_CALIB_HCALIN_SUB1
//     float tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
//     float tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
//     float tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
//     float tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
//     float tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1 = 0;
//     t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1", &tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1);
//     t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1", &tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1);
//     t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1);
//     t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1);
//     t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1", &tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1);

//     // TOWERINFO_CALIB_HCALOUT_SUB1
//     float tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
//     float tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
//     float tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
//     float tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
//     float tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1 = 0;
//     t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1);
//     t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1);
//     t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1);
//     t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1);
//     t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1);

//     // Calo Tower Histograms


//     // -----------------------------------------------------
//     // Random cones
//     // -----------------------------------------------------
//     // RandomCones_r04
//     float random_cone_R_RandomCones_r04 = 0;
//     float random_cone_eta_RandomCones_r04 = 0;
//     float random_cone_phi_RandomCones_r04 = 0;
//     float random_cone_energy_RandomCones_r04 = 0;
//     float random_cone_energy_cemc_RandomCones_r04 = 0;
//     float random_cone_energy_hcalin_RandomCones_r04 = 0;
//     float random_cone_energy_hcalout_RandomCones_r04 = 0;
//     int random_cone_num_towers_RandomCones_r04 = 0;
//     int random_cone_num_towers_cemc_RandomCones_r04 = 0;
//     int random_cone_num_towers_hcalin_RandomCones_r04 = 0;
//     int random_cone_num_towers_hcalout_RandomCones_r04 = 0;
//     int random_cone_num_masked_towers_RandomCones_r04 = 0;
//     int random_cone_num_masked_towers_cemc_RandomCones_r04 = 0;
//     int random_cone_num_masked_towers_hcalin_RandomCones_r04 = 0;
//     int random_cone_num_masked_towers_hcalout_RandomCones_r04 = 0;
//     t->SetBranchAddress("random_cone_R_RandomCones_r04", &random_cone_R_RandomCones_r04);
//     t->SetBranchAddress("random_cone_eta_RandomCones_r04", &random_cone_eta_RandomCones_r04);
//     t->SetBranchAddress("random_cone_phi_RandomCones_r04", &random_cone_phi_RandomCones_r04);
//     t->SetBranchAddress("random_cone_energy_RandomCones_r04", &random_cone_energy_RandomCones_r04);
//     t->SetBranchAddress("random_cone_energy_cemc_RandomCones_r04", &random_cone_energy_cemc_RandomCones_r04);
//     t->SetBranchAddress("random_cone_energy_hcalin_RandomCones_r04", &random_cone_energy_hcalin_RandomCones_r04);
//     t->SetBranchAddress("random_cone_energy_hcalout_RandomCones_r04", &random_cone_energy_hcalout_RandomCones_r04);
//     t->SetBranchAddress("random_cone_num_towers_RandomCones_r04", &random_cone_num_towers_RandomCones_r04);
//     t->SetBranchAddress("random_cone_num_towers_cemc_RandomCones_r04", &random_cone_num_towers_cemc_RandomCones_r04);
//     t->SetBranchAddress("random_cone_num_towers_hcalin_RandomCones_r04", &random_cone_num_towers_hcalin_RandomCones_r04);
//     t->SetBranchAddress("random_cone_num_towers_hcalout_RandomCones_r04", &random_cone_num_towers_hcalout_RandomCones_r04);
//     t->SetBranchAddress("random_cone_num_masked_towers_RandomCones_r04", &random_cone_num_masked_towers_RandomCones_r04);
//     t->SetBranchAddress("random_cone_num_masked_towers_cemc_RandomCones_r04", &random_cone_num_masked_towers_cemc_RandomCones_r04);
//     t->SetBranchAddress("random_cone_num_masked_towers_hcalin_RandomCones_r04", &random_cone_num_masked_towers_hcalin_RandomCones_r04);
//     t->SetBranchAddress("random_cone_num_masked_towers_hcalout_RandomCones_r04", &random_cone_num_masked_towers_hcalout_RandomCones_r04);

//     // RandomCones_r04_Sub1
//     float random_cone_R_RandomCones_r04_Sub1 = 0;
//     float random_cone_eta_RandomCones_r04_Sub1 = 0;
//     float random_cone_phi_RandomCones_r04_Sub1 = 0;
//     float random_cone_energy_RandomCones_r04_Sub1 = 0;
//     float random_cone_energy_cemc_RandomCones_r04_Sub1 = 0;
//     float random_cone_energy_hcalin_RandomCones_r04_Sub1 = 0;
//     float random_cone_energy_hcalout_RandomCones_r04_Sub1 = 0;
//     int random_cone_num_towers_RandomCones_r04_Sub1 = 0;
//     int random_cone_num_towers_cemc_RandomCones_r04_Sub1 = 0;
//     int random_cone_num_towers_hcalin_RandomCones_r04_Sub1 = 0;
//     int random_cone_num_towers_hcalout_RandomCones_r04_Sub1 = 0;
//     int random_cone_num_masked_towers_RandomCones_r04_Sub1 = 0;
//     int random_cone_num_masked_towers_cemc_RandomCones_r04_Sub1 = 0;
//     int random_cone_num_masked_towers_hcalin_RandomCones_r04_Sub1 = 0;
//     int random_cone_num_masked_towers_hcalout_RandomCones_r04_Sub1 = 0;
//     t->SetBranchAddress("random_cone_R_RandomCones_r04_Sub1", &random_cone_R_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_eta_RandomCones_r04_Sub1", &random_cone_eta_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_phi_RandomCones_r04_Sub1", &random_cone_phi_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_energy_RandomCones_r04_Sub1", &random_cone_energy_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_energy_cemc_RandomCones_r04_Sub1", &random_cone_energy_cemc_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_energy_hcalin_RandomCones_r04_Sub1", &random_cone_energy_hcalin_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_energy_hcalout_RandomCones_r04_Sub1", &random_cone_energy_hcalout_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_num_towers_RandomCones_r04_Sub1", &random_cone_num_towers_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_num_towers_cemc_RandomCones_r04_Sub1", &random_cone_num_towers_cemc_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_num_towers_hcalin_RandomCones_r04_Sub1", &random_cone_num_towers_hcalin_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_num_towers_hcalout_RandomCones_r04_Sub1", &random_cone_num_towers_hcalout_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_num_masked_towers_RandomCones_r04_Sub1", &random_cone_num_masked_towers_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_num_masked_towers_cemc_RandomCones_r04_Sub1", &random_cone_num_masked_towers_cemc_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_num_masked_towers_hcalin_RandomCones_r04_Sub1", &random_cone_num_masked_towers_hcalin_RandomCones_r04_Sub1);
//     t->SetBranchAddress("random_cone_num_masked_towers_hcalout_RandomCones_r04_Sub1", &random_cone_num_masked_towers_hcalout_RandomCones_r04_Sub1);

//     // -----------------------------------------------------
//     // Calo window
//     // -----------------------------------------------------
//     unsigned int max_window_vector_size = 0;  
//     float avg_energy_full[11];
//     float std_energy_full[11];
//     float avg_frac_energy_recemc_full[11];
//     float std_frac_energy_recemc_full[11];
//     float avg_frac_energy_hcalin_full[11];
//     float std_frac_energy_hcalin_full[11];
//     float avg_frac_energy_hcalout_full[11];
//     float std_frac_energy_hcalout_full[11];
//     float avg_energy_recemc[11];
//     float std_energy_recemc[11];
//     float avg_energy_hcalin[11];
//     float std_energy_hcalin[11];
//     float avg_energy_hcalout[11];
//     float std_energy_hcalout[11];
//     int num_windows_full[11];
//     int num_windows_recemc[11];
//     int num_windows_hcalin[11];
//     int num_windows_hcalout[11];
//     t->SetBranchAddress("num_window_dims", &max_window_vector_size);
//     t->SetBranchAddress("avg_energy_full", &avg_energy_full);
//     t->SetBranchAddress("std_energy_full", &std_energy_full);
//     t->SetBranchAddress("avg_frac_energy_recemc_full", &avg_frac_energy_recemc_full);
//     t->SetBranchAddress("std_frac_energy_recemc_full", &std_frac_energy_recemc_full);
//     t->SetBranchAddress("avg_frac_energy_hcalin_full", &avg_frac_energy_hcalin_full);
//     t->SetBranchAddress("std_frac_energy_hcalin_full", &std_frac_energy_hcalin_full);
//     t->SetBranchAddress("avg_frac_energy_hcalout_full", &avg_frac_energy_hcalout_full);
//     t->SetBranchAddress("std_frac_energy_hcalout_full", &std_frac_energy_hcalout_full);
//     t->SetBranchAddress("avg_energy_recemc", &avg_energy_recemc);
//     t->SetBranchAddress("std_energy_recemc", &std_energy_recemc);
//     t->SetBranchAddress("avg_energy_hcalin", &avg_energy_hcalin);
//     t->SetBranchAddress("std_energy_hcalin", &std_energy_hcalin);
//     t->SetBranchAddress("avg_energy_hcalout", &avg_energy_hcalout);
//     t->SetBranchAddress("std_energy_hcalout", &std_energy_hcalout);
//     t->SetBranchAddress("num_windows_full", &num_windows_full);
//     t->SetBranchAddress("num_windows_recemc", &num_windows_recemc);
//     t->SetBranchAddress("num_windows_hcalin", &num_windows_hcalin);
//     t->SetBranchAddress("num_windows_hcalout", &num_windows_hcalout);


//     // cemc windows
//     float avg_energy_cemc[11];
//     float std_energy_cemc[11];
//     int num_windows_cemc[11];
//     t->SetBranchAddress("avg_energy_cemc", &avg_energy_cemc);
//     t->SetBranchAddress("std_energy_cemc", &std_energy_cemc);
//     t->SetBranchAddress("num_windows_cemc", &num_windows_cemc);

//     ProcessCaloWindowHistograms();

//     int nentries = t->GetEntries();
//     std::cout << "Number of entries = " << nentries << std::endl;
//     // get first entry
//     t->GetEntry(0);



//     f->Close();

//     return 0;
// }

//-------------------------------------------------------------------------------------------
// IO
//-------------------------------------------------------------------------------------------

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

//-------------------------------------------------------------------------------------------
// CaloWindow
//-------------------------------------------------------------------------------------------
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

void ProcessCaloWindowHistograms()
{
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


        TH3F * h3_energy_frac_energy_cent_recemc_nxm = (TH3F*)f->Get(Form("h3_window_energy_frac_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h3_energy_frac_energy_cent_recemc_nxm ) { std::cout << "h3_energy_frac_energy_cent_recemc_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH3F * h3_energy_frac_energy_cent_hcalin_nxm = (TH3F*)f->Get(Form("h3_window_energy_frac_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h3_energy_frac_energy_cent_hcalin_nxm ) { std::cout << "h3_energy_frac_energy_cent_hcalin_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH3F * h3_energy_frac_energy_cent_hcalout_nxm = (TH3F*)f->Get(Form("h3_window_energy_frac_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h3_energy_frac_energy_cent_hcalout_nxm ) { std::cout << "h3_energy_frac_energy_cent_hcalout_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        
        CaloWindowMultiPanel3D({h3_energy_frac_energy_cent_recemc_nxm, h3_energy_frac_energy_cent_hcalin_nxm, h3_energy_frac_energy_cent_hcalout_nxm}, 
                            {"EMCal", "iHCal", "oHCal"},
                            calo_window_plots+"/energy_frac_energy_slices", "sub_calo_windows_frac_et_all_cent_slices", 
                            Form("E_{T}^{%d #times %d} [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            Form("#LT f(E_{T}^{%d #times %d}) #GT [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            .45, .85, dty, .2, .4, .6, .9, -1, -1, 0, 1.5, false, false);
        
        delete h3_energy_frac_energy_cent_recemc_nxm;
        delete h3_energy_frac_energy_cent_hcalin_nxm;
        delete h3_energy_frac_energy_cent_hcalout_nxm;
        
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
    for ( int idim = 0; idim < k_calo_window_dims_cemc_geom.size(); idim++ ) {

        TH2F * h2_window_energy_cent_cemc_nxm = (TH2F*)f->Get(Form("h2_window_energy_cent_cemc_%dx%d", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second));
        if ( !h2_window_energy_cent_cemc_nxm ) { 
            std::cout << "h2_window_energy_cent_cemc_" << k_calo_window_dims_cemc_geom[idim].first << "x" << k_calo_window_dims_cemc_geom[idim].second << " not found!" << std::endl; 
            continue;
        }
        CaloWindowMultiPanel({h2_window_energy_cent_cemc_nxm},
                            {"EMCal (0.025#times 0.025)"},
                            calo_window_plots+"/energy_cent_slices/cemc_geo", "calo_window_et_all_cent_slices",
                            Form("E_{T}^{%d #times %d} [GeV]", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second),
                            Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second),
                            tagx, tagy, dty, x1, x2, y1_s, y2_s, 0., 0, miny, maxy);
        
        delete h2_window_energy_cent_cemc_nxm;

        TH2F * h2_cemc_energy_minus_avg_energy_cent = (TH2F*)f->Get(Form("h2_window_energy_minus_avg_energy_cent_cemc_%dx%d", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second));
        if ( !h2_cemc_energy_minus_avg_energy_cent ) { 
            std::cout << "h2_cemc_energy_minus_avg_energy_cent_" << k_calo_window_dims_cemc_geom[idim].first << "x" << k_calo_window_dims_cemc_geom[idim].second << " not found!" << std::endl; 
            continue;
        }
        CaloWindowMultiPanel({h2_cemc_energy_minus_avg_energy_cent},
                            {"EMCal (0.025#times 0.025)"},
                            calo_window_plots+"/energy_minus_avg_cent_slices/cemc_geo", "calo_window_et-avget_all_cent_slices",
                            Form("E_{T}^{%d #times %d} [GeV] - #LT E_{T}^{%d #times %d} #GT [GeV]", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second, k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second),
                            Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second),
                            tagx, tagy, dty, x1, x2, y1_s, y2_s, 0., 0, miny, maxy, true);

    }

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

//-------------------------------------------------------------------------------------------
// CaloTowers
//-------------------------------------------------------------------------------------------

