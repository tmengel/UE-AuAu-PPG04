#include <iostream>
#include <vector>
#include <utility> 
#include <sPhenixStyle.C>

const std::string sPHENIX_Tag = "#it{#bf{sPHENIX}} Preliminary";
const std::string DataType_Tag = "Au+Au#kern[0.05]{#sqrt{s_{NN}} = 200 GeV}";
std::string output_dir = "plots/";

bool OVERRIDE = true;

const float V2_VALUES[] = {2.32, 3.39, 4.76, 6.18, 7.03, 7.4, 7.44, 7.23, 6.96};
const float V3_VALUES[] = {1.43, 1.63, 1.82, 1.94, 2.02, 2.03, 1.96, 1.79, 2};

const float X_CENT_BINS[]= {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
const int N_X_CENT_BINS = sizeof(X_CENT_BINS)/sizeof(X_CENT_BINS[0]) - 1;
const int N_X_CENT_BINS_WINDOW = 40;
float X_CENT_BINS_WINDOW[N_X_CENT_BINS_WINDOW+1];
float MAX_X_CENT = 80;

const std::vector < std::pair < unsigned int, unsigned int > > k_calo_window_dims_hcal_geom = {
    {1,1}, {2,2}, {3,4}, {5,6}, {7,8}, {9,10}, {11,12}, {13,13}
}; 

const unsigned int k_window_array_size = k_calo_window_dims_hcal_geom.size();
const int N_WINDOW_BINS = 200;
float MAX_WINDOW_AVG[11];
float MAX_WINDOW_STD[11];
float MAX_WINDOW_NUM[11];

const int N_SUM_ET_Q_BINS = 30;
float SUM_ET_Q_BINS[N_SUM_ET_Q_BINS+1];
float MAX_SUM_ET_Q = 2200;
float MIN_SUM_ET_Q = 200;

const int N_ZVTX_BINS = 100;
float ZVTX_BINS[N_ZVTX_BINS+1];
float MAX_ZVTX = 25;
float MIN_ZVTX = -25;

const int N_RHO_BINS = 200;
float RHO_M_BINS[N_RHO_BINS+1];
float RHO_A_BINS[N_RHO_BINS+1];
float TOWER_BACKGROUND_BINS[N_RHO_BINS+1];
float BACKGROUND_BINS[N_RHO_BINS+1];

float MAX_RHO_M = 0.08;
float MIN_RHO_M = 0;
float MAX_RHO_A = 120;
float MIN_RHO_A = 0;
float MAX_TOWER_BACKGROUND = 1.5;
float MIN_TOWER_BACKGROUND = 0;
float MAX_BACKGROUND = 65;
float MIN_BACKGROUND = 0;

const int N_CONECOMP_BINS = 210;
float CONECOMP_MAX = 210;
float CONECOMP_BINS[N_CONECOMP_BINS+1];

const int N_CONECOMP_SUB1_BINS = 120;
float CONECOMP_SUB1_MAX = 120;
float CONECOMP_SUB1_BINS[N_CONECOMP_SUB1_BINS+1];

const int N_TOWERCOMP_BINS = 500;
float TOWERCOMP_MAX = 1100;
float TOWERCOMP_MIN = 500;
float TOWERCOMP_BINS[N_TOWERCOMP_BINS+1];

const int N_TOWERCOMP_SUB1_BINS = 120;
float TOWERCOMP_SUB1_MAX = 120;
float TOWERCOMP_SUB1_MIN = 0;
float TOWERCOMP_SUB1_BINS[N_TOWERCOMP_SUB1_BINS+1];

const int N_ET_BINS = 15;
float MAX_ET = 150;
float MIN_ET = 0;
float ET_BINS[N_ET_BINS+1];

const int N_RESO_BINS = 120;
float MAX_RESO = 59.5;
float MIN_RESO = -60.5;
float RESO_BINS[N_RESO_BINS+1];

const float AREA_CONE = TMath::Pi()*0.4*0.4;
const float AREA_TOWER_CEMC = (2.0*TMath::Pi()/256.0)*(2.2/96.0);
const float AREA_HCAL_TOWER = (2.0*TMath::Pi()/64.0)*(2.2/24.0);

const float N_CEMC_TOWERS = 256*94;
const float N_HCALIN_TOWERS = 24*64;
const float N_HCALOUT_TOWERS = 24*64;

const int COLORS[] = {kBlack, kRed, kBlue, kCyan, kGreen+2, kViolet, kOrange+2,kAzure-1,  kMagenta+2, kRed+2};
const int MARKERS[] = { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kOpenCircle, kOpenSquare, kOpenTriangleUp};
const float MARKER_SIZE = 1.2;
const float LINE_WIDTH = 2.0;

const float WEIGHT_JET10 = 3.210e-6; //3.210μb
const float WEIGHT_JET30 = 2.178e-9;

const float ABS_ETA_MAX = 0.6;
const float MIN_ET_CUT_TRUTH = 5.0;
const float MIN_ET_CUT_RECO = 10.0;



void SetBins()
{
    for ( int i = 0; i < N_X_CENT_BINS_WINDOW+1; ++i ) { X_CENT_BINS_WINDOW[i] = i*(MAX_X_CENT)/N_X_CENT_BINS_WINDOW; }
    for ( int i = 0; i < N_SUM_ET_Q_BINS+1; ++i ) { SUM_ET_Q_BINS[i] = i*(MAX_SUM_ET_Q-MIN_SUM_ET_Q)/N_SUM_ET_Q_BINS; }
    for ( int i = 0; i < N_ZVTX_BINS+1; ++i ) { ZVTX_BINS[i] = MIN_ZVTX + i*(MAX_ZVTX-MIN_ZVTX)/N_ZVTX_BINS; }
    for ( int i = 0; i < N_TOWERCOMP_BINS+1; ++i ) { TOWERCOMP_BINS[i] = TOWERCOMP_MIN + i*(TOWERCOMP_MAX-TOWERCOMP_MIN)/N_TOWERCOMP_BINS; }
    for ( int i = 0; i < N_TOWERCOMP_SUB1_BINS+1; ++i ) { TOWERCOMP_SUB1_BINS[i] = TOWERCOMP_SUB1_MIN + i*(TOWERCOMP_SUB1_MAX-TOWERCOMP_SUB1_MIN)/N_TOWERCOMP_SUB1_BINS; }
    for ( int i = 0; i < N_ET_BINS+1; ++i ) { ET_BINS[i] = MIN_ET + i*(MAX_ET-MIN_ET)/N_ET_BINS; }
    for ( int i = 0; i < N_RESO_BINS+1; ++i ) { RESO_BINS[i] = MIN_RESO + i*(MAX_RESO-MIN_RESO)/N_RESO_BINS; }
    for ( int i = 0; i < N_RHO_BINS+1; ++i ) { 
        RHO_M_BINS[i] = MIN_RHO_M + i*(MAX_RHO_M-MIN_RHO_M)/N_RHO_BINS; 
        RHO_A_BINS[i] = MIN_RHO_A + i*(MAX_RHO_A-MIN_RHO_A)/N_RHO_BINS; 
        TOWER_BACKGROUND_BINS[i] = MIN_TOWER_BACKGROUND + i*(MAX_TOWER_BACKGROUND-MIN_TOWER_BACKGROUND)/N_RHO_BINS; 
        BACKGROUND_BINS[i] = MIN_BACKGROUND + i*(MAX_BACKGROUND-MIN_BACKGROUND)/N_RHO_BINS;
    }
    for ( unsigned int i = 0; i < k_window_array_size; ++i ) {
        MAX_WINDOW_AVG[i] = 2.0*(k_calo_window_dims_hcal_geom[i].first*k_calo_window_dims_hcal_geom[i].second);
        MAX_WINDOW_STD[i] = 2.0*std::sqrt(1.0*k_calo_window_dims_hcal_geom[i].first*k_calo_window_dims_hcal_geom[i].second);
        MAX_WINDOW_NUM[i] = 2.0*N_HCALIN_TOWERS;
    }
}

float CalcPoisson(const float sum2, const float sum, const float n_tow, const float ncomp, const float ncones,  const float v2 = 0 , const float v3 = 0)
{
    if ( n_tow == 0 ) { return 0;}
    if ( ncones == 0 ) { return 0;}
    float Na = ncomp/ncones;
    float mu = sum/n_tow;
    float sigma2 = sum2/n_tow - (mu*mu);
    float vn2 = v2*v2 + v3*v3;
    float harmcont = Na*Na*vn2*2.0*mu*mu;
    float sigma = TMath::Sqrt(Na*sigma2  + Na*mu*mu + harmcont);
    return sigma;
}

std::string MakeGetDir( const std::string & dir )
{
    if ( gSystem->AccessPathName(dir.c_str()) ) {
        gSystem->Exec(Form("mkdir -p %s", dir.c_str()));
    }
    return dir;
}

std::string GetOutfile( const std::string & outfile , const std::string & dir, const std::string & prefix)
{
    std::string out_dir = MakeGetDir(dir);
    if ( out_dir.back() == '/' ) {
        out_dir.pop_back();
    }
    std::string out_file = MakeGetDir(out_dir + "/" + prefix);
    if ( out_file.back() == '/' ) {
        out_file.pop_back();
    }
    out_file += "/" + outfile;
    if ( out_file.find(".root") == std::string::npos ) {
        out_file += ".root";
    }

    return out_file;
}

void SetStyle(){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    SetsPhenixStyle();
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kRainBow);

    SetBins();
    return;
}

std::string ProcessProbeTree(const std::string & input_file,  const std::string & outfile = "probes.root", const std::string & prefix = "rootfiles");
std::string ProcessEmbedTree(const std::string & input_file, const std::string & outfile = "embed.root", const std::string & prefix = "rootfiles");
std::string ProcessEmbedTree(const std::string & input_file_10, const std::string & input_file_30, const std::string & outfile = "embed.root", const std::string & prefix = "rootfiles");
std::string ProcessConeTree(const std::string & input_file, const std::string & outfile = "cones.root", const std::string & prefix = "rootfiles");
std::string ProcessPoisson(const std::string & input_file, const std::string & outfile = "poisson.root", const std::string & prefix = "rootfiles");
std::string ProcessGlobal(const std::string & input_file, const std::string & outfile = "global.root", const std::string & prefix = "rootfiles");
std::string ProcessRhoTree(const std::string & input_file, const std::string & outfile = "rho.root", const std::string & prefix = "rootfiles");
std::string ProcessWindowTree(const std::string & input_file, const std::string & outfile = "window.root", const std::string & prefix = "rootfiles");


void ProbePlots(const std::string input_file, const std::string & prefix);
void ConePlots(const std::string input_file, const std::string & prefix);
void EmbedPlots(const std::string input_file, const std::string & prefix);

void MultCurves(const std::string input_file, const std::string & prefix);

void SigmaPlots(const std::string input_file_basic, const std::string input_file_random, const std::string & input_file_poisson, const std::string & prefix);
void DeltaPlots(const std::string input_file_basic, const std::string input_file_random, const std::string input_file_probe, const std::string & input_file_embed, const std::string & prefix);

void GlobalPlots(const std::string input_file, const std::string & prefix);
void RhoPlots(const std::string input_file, const std::string & prefix);
void WindowFits(const std::string input_file, const std::string & prefix);
void DeltaFits(const std::string input_file,  const std::string & prefix);

void PPG04_PrelimPlots() {

    SetStyle();
    
    const std::string & input_file_basic = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/NEW/feb10_basic.root";
    const std::string & input_file_random = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/NEW/feb10_random.root";
    const std::string & input_file_embed_10 = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/FEB7/feb7_embed.root";

    output_dir = MakeGetDir("preliminary/");

    std::cout << "PPG04 plots starting." << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;
    std::vector<std::string> infiles = {input_file_basic, input_file_random, input_file_embed};
    std::cout << "Input files: " << std::endl;
    for ( auto infile : infiles ) {
        std::cout << infile << std::endl;
    }

    OVERRIDE = false;
    LoadNevents(input_file_basic);
    auto global_file = ProcessGlobal(input_file_basic, "global.root");
    auto rho_file = ProcessRhoTree(input_file_basic, "rho.root");
    // auto rho_file_random = ProcessRhoTree(input_file_random, "rho-random.root");
    auto cone_file_basic = ProcessConeTree(input_file_basic, "basic.root");
    auto cone_file_random = ProcessConeTree(input_file_random, "random.root");
    // auto poisson_file = ProcessPoisson(input_file_basic, "poisson.root");
    // auto window_file = ProcessWindowTree(input_file_basic, "window.root");
    // auto mbd_window_file = ProcessWindowTreeMBD(input_file_basic, "windowMBD.root");
    // auto window_slice_file = ProcessWindowSlices(input_file_basic, "window-slices.root");
    auto probe_file = ProcessProbeTree(input_file_basic, "probe.root");
    auto embed_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/offline/plots-preliminary/Embed/probes/embed/embed.root";

    std::cout << "PPG04 creating plots finished." << std::endl;

    // GlobalPlots(global_file, "global");
    RhoPlots(rho_file, "rho");
    // RhoPlots(rho_file_random, "rho/random");
    MultCurves(embed_file, "mult_method");
    // WindowFits(window_file, "window");
    // WindowFitsMBD(mbd_window_file, "windowMBD");
    // DeltaPlotsXcheck(cone_file_basic, cone_file_random, "deltaX");
    // DO_GAMMA = true;
    // DeltaPlots(cone_file_basic, cone_file_random, probe_file, embed_file, "delta-gamma");
    DO_GAMMA = false;
    DeltaPlots(cone_file_basic, cone_file_random, probe_file, embed_file, "delta");
    // WindowSlicesPlots(input_file_basic, "window-slices");
    
    gSystem->Exit(0);
  
}

std::string ProcessProbeTree(const std::string & input_file, const std::string & outfile, const std::string & prefix )
{

    std::string output_file = GetOutfile(outfile, output_dir, prefix);
    if ( !OVERRIDE ) {
        std::cout << "ProcessProbeTree: OVERRIDE is false, returning " << output_file << std::endl;
        return output_file;
    } else {
        std::cout << "ProcessProbeTree: " << output_file << std::endl;
    }

    TFile * f = new TFile(input_file.c_str(), "READ");
    if( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl;  exit(1); }

    TTree * t = (TTree*)f->Get("T");
    if ( !t ) { std::cout << "Tree not found in file " << input_file << std::endl; exit(1); }
    
    // tree branches 
        int centrality = 0;
        t->SetBranchAddress("centrality", &centrality);


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

        float probe_jet_truth_energy = 0;
        float probe_jet_sub1_truth_energy = 0;
        float probe_jet_energy = 0;
        float probe_jet_sub1_energy = 0;
        float probe_jet_area = 0;
        float probe_jet_sub1_area = 0;
        int probe_jet_num_towers = 0;
        int probe_jet_sub1_num_towers = 0;
        float probe_jet_truth_eta = 0;
        float probe_jet_eta = 0;
        t->SetBranchAddress("probe_jet_truth_energy", &probe_jet_truth_energy);
        t->SetBranchAddress("probe_jet_sub1_truth_energy", &probe_jet_sub1_truth_energy);
        t->SetBranchAddress("probe_jet_area", &probe_jet_area);
        t->SetBranchAddress("probe_jet_sub1_area", &probe_jet_sub1_area);
        t->SetBranchAddress("probe_jet_num_towers", &probe_jet_num_towers);
        t->SetBranchAddress("probe_jet_sub1_num_towers", &probe_jet_sub1_num_towers);
        t->SetBranchAddress("probe_jet_energy", &probe_jet_energy);
        t->SetBranchAddress("probe_jet_sub1_energy", &probe_jet_sub1_energy);
        t->SetBranchAddress("probe_jet_truth_eta", &probe_jet_truth_eta);
        t->SetBranchAddress("probe_jet_eta", &probe_jet_eta);

    int nentries = t->GetEntries();
    std::cout << "ProcessProbeTree: Processing " << nentries << " events" << std::endl;
  
  
    TH2F * h2_area_res_vs_x = new TH2F("h2_area_res", "h2_area_res", N_X_CENT_BINS, X_CENT_BINS, N_RESO_BINS, RESO_BINS);
    TH2F * h2_mult_res_vs_x = new TH2F("h2_mult_res", "h2_area_res", N_X_CENT_BINS, X_CENT_BINS, N_RESO_BINS, RESO_BINS);
    TH2F * h2_sub1_res_vs_x = new TH2F("h2_sub1_res", "h2_area_res", N_X_CENT_BINS, X_CENT_BINS, N_RESO_BINS, RESO_BINS);
    TH1F * h1_rhoA = new TH1F("h1_rhoA", "h1_rhoA", N_ET_BINS, ET_BINS);
    TH1F * h1_rhoM = new TH1F("h1_rhoM", "h1_rhoM", N_ET_BINS, ET_BINS);
    TH1F * h1_probe_sub1_et = new TH1F("h1_probe_sub1_et", "h1_probe_sub1_et", N_ET_BINS, ET_BINS);
    TH1F * h1_probe_truth_et = new TH1F("h1_probe_truth_et", "h1_probe_truth_et", N_ET_BINS, ET_BINS);
    TH1F * h1_probe_truth_sub1_et = new TH1F("h1_probe_truth_sub1_et", "h1_probe_truth_sub1_et", N_ET_BINS, ET_BINS);
    TH1F * h1_probe_et = new TH1F("h1_probe_et", "h1_probe_et", N_ET_BINS, ET_BINS);

    std::vector<TH1F*> h1s = {h1_rhoA, h1_rhoM, h1_probe_sub1_et, h1_probe_truth_et, h1_probe_truth_sub1_et, h1_probe_et};
    for ( auto h1 : h1s ) {
        h1->GetXaxis()->SetTitle("E_{T} [GeV]");
        h1->GetYaxis()->SetTitle("Counts");
    }
   
    std::vector<TH2F*> h2s = {h2_area_res_vs_x, h2_mult_res_vs_x, h2_sub1_res_vs_x};
    for ( auto h2 : h2s ) {
        h2->GetXaxis()->SetTitle("Centrality [%]");
        h2->GetYaxis()->SetTitle("#delta E_{T}^{Raw}^{Probe} [GeV]");
    }
    std::vector<std::string> labs = {"Area", "Multiplicity", "Iterative"};

    for ( int i = 0; i < nentries; ++i ) {
       
        t->GetEntry(i);
        float xaxis_var =1.0*centrality;
        if ( xaxis_var < X_CENT_BINS[0] || xaxis_var > X_CENT_BINS[N_X_CENT_BINS] ) { continue; }
        if ( probe_jet_energy == 0 || probe_jet_sub1_energy == 0 || probe_jet_truth_energy == 0 ) { continue; }
        if ( std::fabs(probe_jet_eta) > ABS_ETA_MAX|| std::fabs(probe_jet_truth_eta) > ABS_ETA_MAX) { continue; }

        float res_area = probe_jet_energy - rho_val_TowerRho_AREA*probe_jet_area - probe_jet_truth_energy;
        float res_mult = probe_jet_energy - rho_val_TowerRho_MULT*probe_jet_num_towers - probe_jet_truth_energy;
        float res_sub1 = probe_jet_sub1_energy - probe_jet_sub1_truth_energy;

        h2_area_res_vs_x->Fill(xaxis_var, res_area);
        h2_mult_res_vs_x->Fill(xaxis_var, res_mult);
        h2_sub1_res_vs_x->Fill(xaxis_var, res_sub1);
        h1_probe_et->Fill(probe_jet_energy);
        h1_rhoA->Fill(rho_val_TowerRho_AREA*probe_jet_area);
        h1_rhoM->Fill(rho_val_TowerRho_MULT*probe_jet_num_towers);
        h1_probe_sub1_et->Fill(probe_jet_sub1_energy);
        h1_probe_truth_et->Fill(probe_jet_truth_energy);
        h1_probe_truth_sub1_et->Fill(probe_jet_sub1_truth_energy);

    }


    TFile * fout = new TFile(output_file.c_str(), "RECREATE");
    for ( auto h2 : h2s ) {
        h2->Write();
    }
    for ( auto h1 : h1s ) {
        h1->Write();
    }
    fout->Close();
    f->Close();
    std::cout << "ProcessProbeTree: Wrote to " << output_file << std::endl;
    return output_file;
}

std::string ProcessEmbedTree(const std::string & input_file_10, const std::string & input_file_30,  const std::string & outfile, const std::string & prefix  )
{


    std::string output_file = GetOutfile(outfile, output_dir, prefix);
    std::string outdir = MakeGetDir(output_dir+"/"+prefix);
    if ( !OVERRIDE ) {
        std::cout << "ProcessEmbedTree: OVERRIDE is false, returning " << output_file << std::endl;
        return output_file;
    } else {
        std::cout << "ProcessEmbedTree: " << output_file << std::endl;
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

    TH2F * h2_area_res_vs_x = new TH2F("h2_area_res", "h2_area_res", N_X_CENT_BINS, X_CENT_BINS, N_RESO_BINS, RESO_BINS);
    TH2F * h2_mult_res_vs_x = new TH2F("h2_mult_res", "h2_mult_res", N_X_CENT_BINS, X_CENT_BINS, N_RESO_BINS, RESO_BINS);
    TH2F * h2_sub1_res_vs_x = new TH2F("h2_sub1_res", "h2_sub1_res", N_X_CENT_BINS, X_CENT_BINS, N_RESO_BINS, RESO_BINS);
    TH3F * h2_et_ntruth_cent = new TH3F("h2_et_ntruth_cent", "h2_et_ntruth_cent", N_X_CENT_BINS, X_CENT_BINS, N_ET_BINS, ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
    TH3F * h2_et_ntruth_cemc_cent = new TH3F("h2_et_ntruth_cemc_cent", "h2_et_ntruth_cemc_cent", N_X_CENT_BINS, X_CENT_BINS, N_ET_BINS, ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
    TH3F * h2_et_ntruth_hcalin_cent = new TH3F("h2_et_ntruth_hcalin_cent", "h2_et_ntruth_hcalin_cent", N_X_CENT_BINS, X_CENT_BINS, N_ET_BINS, ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
    TH3F * h2_et_ntruth_hcalout_cent = new TH3F("h2_et_ntruth_hcalout_cent", "h2_et_ntruth_hcalout_cent", N_X_CENT_BINS, X_CENT_BINS, N_ET_BINS, ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
   
    // initialize tprofiles
    TProfile * p_et_ntruth[N_X_CENT_BINS];
    TProfile * p_et_ntruth_cemc[N_X_CENT_BINS];
    TProfile * p_et_ntruth_hcalin[N_X_CENT_BINS];
    TProfile * p_et_ntruth_hcalout[N_X_CENT_BINS];

    for (  unsigned ibin = 0; ibin < N_X_CENT_BINS; ibin++ ){
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

    float WEIGHT_JET10 = 3.210e-6; //3.210μb
    float WEIGHT_JET30 = 2.178e-9;

    // jet 10
    std::cout << "Looping over jet 10" << std::endl;
    for ( int i = 0; i < nentries; ++i ) {
       
        t_10->GetEntry(i);
       
        float xaxis_var =1.0*centrality_10;
        if ( xaxis_var <0 || xaxis_var > MAX_X_CENT ) { continue; }

        int icourse_bin = -1;
        for ( unsigned ibin = 0; ibin < N_X_CENT_BINS; ibin++ ){
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
            if ( std::fabs(teta) > ABS_ETA_MAX) { continue; }
            if ( tet < MIN_ET_CUT_TRUTH) { continue; }
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
                if ( ret < MIN_ET_CUT_RECO) { continue; }
                if ( std::fabs(reta) > ABS_ETA_MAX) { continue; }
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
    
            h2_et_ntruth_cent->Fill(xaxis_var, et, ncomp, WEIGHT_JET10);
            h2_et_ntruth_cemc_cent->Fill(xaxis_var, et_cemc, ncomp_cemc, WEIGHT_JET10);
            h2_et_ntruth_hcalin_cent->Fill(xaxis_var, et_hcalin, ncomp_hcalin, WEIGHT_JET10);
            h2_et_ntruth_hcalout_cent->Fill(xaxis_var, et_hcalout, ncomp_hcalout, WEIGHT_JET10);
            ntruth_matched++;

            p_et_ntruth[icourse_bin]->Fill(et, ncomp, WEIGHT_JET10);
            p_et_ntruth_cemc[icourse_bin]->Fill(et_cemc, ncomp_cemc, WEIGHT_JET10);
            p_et_ntruth_hcalin[icourse_bin]->Fill(et_hcalin, ncomp_hcalin, WEIGHT_JET10);
            p_et_ntruth_hcalout[icourse_bin]->Fill(et_hcalout, ncomp_hcalout, WEIGHT_JET10);

            
    
        }
        
    }
    
    // jet 30
    std::cout << "Looping over jet 30" << std::endl;
    for ( int i = 0; i < nentries; ++i ) {
       
        t_30->GetEntry(i);
       
        float xaxis_var =1.0*centrality_30;
        if ( xaxis_var <0 || xaxis_var > MAX_X_CENT ) { continue; }

        int icourse_bin = -1;
        for ( unsigned ibin = 0; ibin < N_X_CENT_BINS; ibin++ ){
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
            if ( std::fabs(teta) > ABS_ETA_MAX) { continue; }
            if ( tet < MIN_ET_CUT_TRUTH) { continue; }
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
                if ( ret < MIN_ET_CUT_RECO) { continue; }
                if ( std::fabs(reta) > ABS_ETA_MAX) { continue; }
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
    
            h2_et_ntruth_cent->Fill(xaxis_var, et, ncomp, WEIGHT_JET30);
            h2_et_ntruth_cemc_cent->Fill(xaxis_var, et_cemc, ncomp_cemc, WEIGHT_JET30);
            h2_et_ntruth_hcalin_cent->Fill(xaxis_var, et_hcalin, ncomp_hcalin, WEIGHT_JET30);
            h2_et_ntruth_hcalout_cent->Fill(xaxis_var, et_hcalout, ncomp_hcalout, WEIGHT_JET30);
            ntruth_matched++;

            p_et_ntruth[icourse_bin]->Fill(et, ncomp, WEIGHT_JET30);
            p_et_ntruth_cemc[icourse_bin]->Fill(et_cemc, ncomp_cemc, WEIGHT_JET30);
            p_et_ntruth_hcalin[icourse_bin]->Fill(et_hcalin, ncomp_hcalin, WEIGHT_JET30);
            p_et_ntruth_hcalout[icourse_bin]->Fill(et_hcalout, ncomp_hcalout, WEIGHT_JET30);
    
        }
        
    }

    std::cout << "Total truth jets: " << ntruth_total << " matched: " << ntruth_matched << " total emb jets: " << nemb_total << std::endl;
    
    WEIGHT_JET10 = 1.0;
    WEIGHT_JET30 = 1.0;
    TH1F * h1_et_vs_ntruth[N_X_CENT_BINS];
    TH1F * h1_et_vs_ntruth_cemc[N_X_CENT_BINS];
    TH1F * h1_et_vs_ntruth_hcalin[N_X_CENT_BINS];
    TH1F * h1_et_vs_ntruth_hcalout[N_X_CENT_BINS];
    TH2F * h2s_cent[N_X_CENT_BINS];
    TH2F * h2s_cemc[N_X_CENT_BINS];
    TH2F * h2s_hcalin[N_X_CENT_BINS];
    TH2F * h2s_hcalout[N_X_CENT_BINS];
    for ( int ibin = 0; ibin < N_X_CENT_BINS; ibin++ ) 
    {
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
            if ( std::fabs(teta) > ABS_ETA_MAX) { continue; }
            if ( tet < MIN_ET_CUT_TRUTH) { continue; }
            int ncomp = sim_jet_num_towers_10->at(itruth);
            bool is_matched = false;   
            float dr_current = 1000.0; 
            int iemb_match = -1;
            for ( unsigned iemb = 0; iemb < nemb; iemb++ ) {
                float reta = emb_jet_eta_10->at(iemb);
                float rphi = emb_jet_phi_10->at(iemb);
                float ret  = emb_jet_energy_10->at(iemb);
                if ( ret < MIN_ET_CUT_RECO) { continue; }
                if ( std::fabs(reta) > ABS_ETA_MAX) { continue; }
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


            float rhoM = rho_val_TowerRho_MULT_10;
            float rhoA = rho_val_TowerRho_AREA_10;
            float rhoM_cemc = rho_val_TowerRho_MULT_CEMC_10;
            float rhoA_cemc = rho_val_TowerRho_AREA_CEMC_10;
            float rhoM_hcalin = rho_val_TowerRho_MULT_HCALIN_10;
            float rhoA_hcalin = rho_val_TowerRho_AREA_HCALIN_10;
            float rhoM_hcalout = rho_val_TowerRho_MULT_HCALOUT_10;
            float rhoA_hcalout = rho_val_TowerRho_AREA_HCALOUT_10;

            float mult_et = 0;
            if ( ncorr_cemc > 0 && ncorr_hcalin > 0 && ncorr_hcalout > 0 ){
                mult_et = et - (rhoM_cemc*(ntowers_cemc - ncorr_cemc) 
                            + rhoM_hcalin*(ntowers_hcalin - ncorr_hcalin) 
                            + rhoM_hcalout*(ntowers_hcalout - ncorr_hcalout));
            } else {
                mult_et = et - (rhoM*(ntowers - ncorr)); // default to full mult if no subdetectors
            }
            // mult_et = et - (rhoM*(ntowers - ncorr_cemc));
            mult_et = et - (rhoM*(ntowers - ncorr)); // default to full mult if no subdetectors
            // mult_et = et - ( rhoM_cemc*(ntowers_cemc - ncorr_cemc) 
            //                 + rhoM_hcalin*(ntowers_hcalin) 
            //                 + rhoM_hcalout*(ntowers_hcalout) );
            float area_et = et - rhoA*area;

            float cone_res_area = area_et - tet;
            float cone_res_mult = mult_et - tet;

            if ( area_et > 10 ) {
                h2_area_res_vs_x->Fill(xaxis_var, cone_res_area, WEIGHT_JET10);
                h2_rhoA_vs_x->Fill(xaxis_var, rhoA, WEIGHT_JET10);
            }
            
            if ( mult_et > 10 ) {
                h2_mult_res_vs_x->Fill(xaxis_var, cone_res_mult, WEIGHT_JET10);
                h2_rhoM_vs_x->Fill(xaxis_var, rhoM, WEIGHT_JET10);
            }

        } // end loop over truth jets

        // match truth to embedded sub1
        unsigned int ntruth_retower = sim_retower_jet_eta_10->size();
        unsigned int nsub1 = emb_jet_sub1_eta_10->size();
        for ( unsigned itruth = 0; itruth < ntruth_retower; itruth++) {
            float teta = sim_retower_jet_eta_10->at(itruth);
            float tphi = sim_retower_jet_phi_10->at(itruth);
            float tet = sim_retower_jet_energy_10->at(itruth);
            if ( std::fabs(teta) > ABS_ETA_MAX) { continue; }
            if ( tet < MIN_ET_CUT_TRUTH) { continue; }
            int ncomp = sim_retower_jet_num_towers_10->at(itruth);
            bool is_matched = false;   
            float dr_current = 1000.0; 
            int iemb_match = -1;
            for ( unsigned iemb = 0; iemb < nsub1; iemb++ ) {
                float reta = emb_jet_sub1_eta_10->at(iemb);
                float rphi = emb_jet_sub1_phi_10->at(iemb);
                float ret  = emb_jet_sub1_energy_10->at(iemb);
                if ( ret < MIN_ET_CUT_RECO) { continue; }
                if ( std::fabs(reta) > ABS_ETA_MAX) { continue; }
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
                h2_sub1_res_vs_x->Fill(xaxis_var, res_sub1, WEIGHT_JET10);
            }
        } // end loop over truth jets (retower)

    
    } // end loop over events
    
    // jet 30
    std::cout << "Looping over jet 30" << std::endl;
    for ( int i = 0; i < nentries; ++i ) {
       
        t_30->GetEntry(i);
           
        float xaxis_var =1.0*centrality_30;
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
            if ( std::fabs(teta) > ABS_ETA_MAX) { continue; }
            if ( tet < MIN_ET_CUT_TRUTH) { continue; }
            int ncomp = sim_jet_num_towers_30->at(itruth);
            bool is_matched = false;   
            float dr_current = 1000.0; 
            int iemb_match = -1;
            for ( unsigned iemb = 0; iemb < nemb; iemb++ ) {
                float reta = emb_jet_eta_30->at(iemb);
                float rphi = emb_jet_phi_30->at(iemb);
                float ret  = emb_jet_energy_30->at(iemb);
                if ( ret < MIN_ET_CUT_RECO) { continue; }
                if ( std::fabs(reta) > ABS_ETA_MAX) { continue; }
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

    
            float rhoM = rho_val_TowerRho_MULT_30;
            float rhoA = rho_val_TowerRho_AREA_30;
            float rhoM_cemc = rho_val_TowerRho_MULT_CEMC_30;
            float rhoA_cemc = rho_val_TowerRho_AREA_CEMC_30;
            float rhoM_hcalin = rho_val_TowerRho_MULT_HCALIN_30;
            float rhoA_hcalin = rho_val_TowerRho_AREA_HCALIN_30;
            float rhoM_hcalout = rho_val_TowerRho_MULT_HCALOUT_30;
            float rhoA_hcalout = rho_val_TowerRho_AREA_HCALOUT_30;
    
 
            float mult_et = 0;
            if ( ncorr_cemc > 0 && ncorr_hcalin > 0 && ncorr_hcalout > 0 ){
                mult_et = et - (rhoM_cemc*(ntowers_cemc - ncorr_cemc) 
                            + rhoM_hcalin*(ntowers_hcalin - ncorr_hcalin) 
                            + rhoM_hcalout*(ntowers_hcalout - ncorr_hcalout));
            } else {
                mult_et = et - (rhoM*(ntowers - ncorr)); // default to full mult if no subdetectors
            }
            // mult_et = et - (rhoM*(ntowers - ncorr));
            // mult_et = et - ( rhoM_cemc*(ntowers_cemc - ncorr_cemc) 
            //                 + rhoM_hcalin*(ntowers_hcalin) 
            //                 + rhoM_hcalout*(ntowers_hcalout) );
            // mult_et = et - (rhoM*(ntowers));

            // mult_et = et - (rhoM*(ntowers - ncorr_cemc));
            mult_et = et - (rhoM*(ntowers - ncorr)); // default to full mult if no subdetectors
            float area_et = et - rhoA*area;

            float cone_res_area = area_et - tet;
            float cone_res_mult = mult_et - tet;

            if ( area_et > 10 ) {
                h2_area_res_vs_x->Fill(xaxis_var, cone_res_area, WEIGHT_JET30);
                h2_rhoA_vs_x->Fill(xaxis_var, rhoA, WEIGHT_JET30);
            }

            if ( mult_et > 10 ) {
                h2_mult_res_vs_x->Fill(xaxis_var, cone_res_mult, WEIGHT_JET30);
                h2_rhoM_vs_x->Fill(xaxis_var, rhoM, WEIGHT_JET30);
            }


        
        } // end loop over truth jets
    
        // match truth to embedded sub1
        for ( unsigned itruth = 0; itruth < ntruth_retower; itruth++) {
            float teta = sim_retower_jet_eta_30->at(itruth);
            float tphi = sim_retower_jet_phi_30->at(itruth);
            float tet = sim_retower_jet_energy_30->at(itruth);
            if ( std::fabs(teta) > ABS_ETA_MAX) { continue; }
            if ( tet < MIN_ET_CUT_TRUTH) { continue; }
            int ncomp = sim_retower_jet_num_towers_30->at(itruth);
            bool is_matched = false;   
            float dr_current = 1000.0; 
            int iemb_match = -1;
            for ( unsigned iemb = 0; iemb < nsub1; iemb++ ) {
                float reta = emb_jet_sub1_eta_30->at(iemb);
                float rphi = emb_jet_sub1_phi_30->at(iemb);
                float ret  = emb_jet_sub1_energy_30->at(iemb);
                if ( ret < MIN_ET_CUT_RECO) { continue; }
                if ( std::fabs(reta) > ABS_ETA_MAX) { continue; }
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
                h2_sub1_res_vs_x->Fill(xaxis_var, res_sub1, WEIGHT_JET30);
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

           
            g_std_devs[ihist]->SetPoint(ibin, h2->GetXaxis()->GetBinCenter(ibin+1), std);
            g_std_devs[ihist]->SetPointError(ibin, h2->GetXaxis()->GetBinWidth(ibin+1)/2.0, 0);

           
            h1->Draw("p");
            leg->AddEntry(h1, Form("#mu = %0.2f, #sigma = %0.2f", avg, std), "p");

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
    for (  unsigned ibin = 0; ibin < N_X_CENT_BINS; ibin++ ){
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

    TFile * fout = new TFile(output_file.c_str(), "RECREATE");
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
    for ( int ibin = 0; ibin < N_X_CENT_BINS; ibin++ ) {
        p_et_ntruth[ibin]->Write();
        p_et_ntruth_cemc[ibin]->Write();
        p_et_ntruth_hcalin[ibin]->Write();
        p_et_ntruth_hcalout[ibin]->Write();
    }

    fout->Close();
    
    f_10->Close();
    f_30->Close();
    return output_file;
    
}

std::string ProcessConeTree(const std::string & input_file, const std::string & outfile, const std::string & prefix  )
{

    
    std::string output_file = GetOutfile(outfile, output_dir, prefix);
    if ( !OVERRIDE ) {
        std::cout << "ProcessConeTree: OVERRIDE is false, returning " << output_file << std::endl;
        return output_file;
    } else {
        std::cout << "ProcessConeTree: " << output_file << std::endl;
    }

    TFile * f = new TFile(input_file.c_str(), "READ");
    if( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl;  exit(1); }

    TTree * t = (TTree*)f->Get("T");
    if ( !t ) { std::cout << "Tree not found in file " << input_file << std::endl; exit(1); }
    
    // tree branches 
        int centrality = 0;
        t->SetBranchAddress("centrality", &centrality);

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
    std::cout << "ProcessConeTree: Processing " << nentries << " events" << std::endl;
  
  
    TH2F * h2_area_res_vs_x = new TH2F("h2_area_res", "h2_area_res", N_X_CENT_BINS, X_CENT_BINS, N_RESO_BINS, RESO_BINS);
    TH2F * h2_mult_res_vs_x = new TH2F("h2_mult_res", "h2_area_res", N_X_CENT_BINS, X_CENT_BINS, N_RESO_BINS, RESO_BINS);
    TH2F * h2_sub1_res_vs_x = new TH2F("h2_sub1_res", "h2_area_res", N_X_CENT_BINS, X_CENT_BINS, N_RESO_BINS, RESO_BINS);
    TProfile * p2_area_res = new TProfile("p2_area_res", "p2_area_res", N_X_CENT_BINS, X_CENT_BINS);
    TProfile * p2_mult_res = new TProfile("p2_mult_res", "p2_mult_res", N_X_CENT_BINS, X_CENT_BINS);
    TProfile * p2_sub1_res = new TProfile("p2_sub1_res", "p2_sub1_res", N_X_CENT_BINS, X_CENT_BINS);
    std::vector<TH2F*> h2s = {h2_area_res_vs_x, h2_mult_res_vs_x, h2_sub1_res_vs_x};
    for ( auto h2 : h2s ) {
        h2->GetXaxis()->SetTitle("Centrality [%]");
        h2->GetYaxis()->SetTitle("#delta E_{T}^{Raw}^{Cone} [GeV]");
        
    }
    std::vector<TProfile*> p2s = {p2_area_res, p2_mult_res, p2_sub1_res};
    for ( auto p2 : p2s ) {
        p2->GetXaxis()->SetTitle("Centrality [%]");
        p2->GetYaxis()->SetTitle("#delta E_{T}^{Raw}^{Cone} [GeV]");
    }
    std::vector<std::string> labs = {"Area", "Multiplicity", "Iterative"};


    for ( int i = 0; i < nentries; ++i ) {
       
        t->GetEntry(i);
        float xaxis_var =1.0*centrality;
        if ( xaxis_var < X_CENT_BINS[0] || xaxis_var > X_CENT_BINS[N_X_CENT_BINS] ) { continue; }
       
        float cone_eta = random_cone_eta_RandomCones_r04;
        float cone_sub1_eta = random_cone_eta_RandomCones_r04_Sub1;
        // if ( std::fabs(cone_eta) > ABS_ETA_MAX|| std::fabs(cone_sub1_eta) > ABS_ETA_MAX) { continue; }
   
        int ncomp_cone_corr = random_cone_num_towers_RandomCones_r04 - random_cone_num_masked_towers_RandomCones_r04;
        int ncomp_cone_corr_sub1 = random_cone_num_towers_RandomCones_r04_Sub1 - random_cone_num_masked_towers_RandomCones_r04_Sub1;

        float area_bkgd = rho_val_TowerRho_AREA*AREA_CONE;
        // float area_bkgd = ( rho_val_TowerRho_AREA_CEMC*AREA_CONE + rho_val_TowerRho_AREA_HCALIN*AREA_CONE + rho_val_TowerRho_AREA_HCALOUT*AREA_CONE );
        float mult_bkgd = ( rho_val_TowerRho_MULT_CEMC*(random_cone_num_towers_cemc_RandomCones_r04 - random_cone_num_masked_towers_cemc_RandomCones_r04) 
                            + rho_val_TowerRho_MULT_HCALIN*(random_cone_num_towers_hcalin_RandomCones_r04 - random_cone_num_masked_towers_hcalin_RandomCones_r04) 
                            + rho_val_TowerRho_MULT_HCALOUT*(random_cone_num_towers_hcalout_RandomCones_r04-random_cone_num_masked_towers_hcalout_RandomCones_r04) );
        float res_area  = random_cone_energy_RandomCones_r04 - area_bkgd;
        float pt_cone = random_cone_energy_cemc_RandomCones_r04 + random_cone_energy_hcalin_RandomCones_r04 + random_cone_energy_hcalout_RandomCones_r04;
       
        // float res_area  = pt_cone - area_bkgd;  
        float res_mult  = random_cone_energy_cemc_RandomCones_r04 + random_cone_energy_hcalin_RandomCones_r04 + random_cone_energy_hcalout_RandomCones_r04 - mult_bkgd;
        float res_sub1  = random_cone_energy_RandomCones_r04_Sub1;


        h2_area_res_vs_x->Fill(xaxis_var, res_area);
        h2_mult_res_vs_x->Fill(xaxis_var, res_mult);
        h2_sub1_res_vs_x->Fill(xaxis_var, res_sub1);
        p2_area_res->Fill(xaxis_var, res_area);
        p2_mult_res->Fill(xaxis_var, res_mult);
        p2_sub1_res->Fill(xaxis_var, res_sub1);

        
    }

    TH1F * h1_area_res[N_X_CENT_BINS];
    TH1F * h1_mult_res[N_X_CENT_BINS];
    TH1F * h1_sub1_res[N_X_CENT_BINS];
    for ( int ibin = 0; ibin < N_X_CENT_BINS; ibin++ ) {
        // h2_area_res_vs_x->GetXaxis()->SetRange(ibin+1, ibin+1);
        // h2_mult_res_vs_x->GetXaxis()->SetRange(ibin+1, ibin+1);
        // h2_sub1_res_vs_x->GetXaxis()->SetRange(ibin+1, ibin+1);
        h1_area_res[ibin] = ( TH1F * )h2_area_res_vs_x->ProjectionY(Form("h1_area_res_%d", ibin), ibin+1, ibin+1);
        h1_mult_res[ibin] = ( TH1F * )h2_mult_res_vs_x->ProjectionY(Form("h1_mult_res_%d", ibin), ibin+1, ibin+1);
        h1_sub1_res[ibin] = ( TH1F * )h2_sub1_res_vs_x->ProjectionY(Form("h1_sub1_res_%d", ibin), ibin+1, ibin+1);
    }


    TFile * fout = new TFile(output_file.c_str(), "RECREATE");
    for ( auto h2 : h2s ) {
        h2->Write();
    }
    for ( auto p2 : p2s ) {
        p2->Write();
    }
    for ( int ibin = 0; ibin < N_X_CENT_BINS; ibin++ ) {
        h1_area_res[ibin]->Write();
        h1_mult_res[ibin]->Write();
        h1_sub1_res[ibin]->Write();
    }
    fout->Close();
    f->Close();

    std::cout << "ProcessConeTree: Wrote to " << output_file << std::endl;

    return output_file;
}

std::string ProcessPoisson(const std::string & input_file, const std::string & outfile, const std::string & prefix  )
{

   
    std::string output_file = GetOutfile(outfile, output_dir, prefix);
    if ( !OVERRIDE ) {
        std::cout << "ProcessPoisson: OVERRIDE is false, returning " << output_file << std::endl;
        return output_file;
    } else {
        std::cout << "ProcessPoisson: " << output_file << std::endl;
    }

    TFile * f = new TFile(input_file.c_str(), "READ");
    if( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl;  exit(1); }

    TTree * t = (TTree*)f->Get("T");
    if ( !t ) { std::cout << "Tree not found in file " << input_file << std::endl; exit(1); }
    
    // tree branches 
        int centrality = 0;
        t->SetBranchAddress("centrality", &centrality);

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
    std::cout << "ProcessPoisson: Processing " << nentries << " events" << std::endl;
  
  
    float SUM_ET2[N_X_CENT_BINS];
    float SUM_E[N_X_CENT_BINS];
    float TOTAL_TOWERS_NOT_MASKED[N_X_CENT_BINS];
    float TOTAL_TOWERS_FIRED[N_X_CENT_BINS];
    float N_TOWERS_TOTAL[N_X_CENT_BINS];
    float SUM_CONE_COMPS[N_X_CENT_BINS];
    float N_CONES_THIS_BIN[N_X_CENT_BINS];

    TH1F * h1_dummy = new TH1F("h1_dummy", "h1_dummy", N_X_CENT_BINS, X_CENT_BINS);

    for ( int i = 0; i < nentries; ++i ) {
       
        t->GetEntry(i);
        float xaxis_var =1.0*centrality;
        if ( xaxis_var < X_CENT_BINS[0] || xaxis_var > X_CENT_BINS[N_X_CENT_BINS] ) { continue; }
        
        float cone_eta = random_cone_eta_RandomCones_r04;
        float cone_sub1_eta = random_cone_eta_RandomCones_r04_Sub1;
        if ( std::fabs(cone_eta) > ABS_ETA_MAX|| std::fabs(cone_sub1_eta) > ABS_ETA_MAX) { continue; }

        int THIS_BIN = -1;
        for ( int ibin = 0; ibin < N_X_CENT_BINS; ibin++ ) {
            if ( xaxis_var >= X_CENT_BINS[ibin] && xaxis_var <= X_CENT_BINS[ibin+1] ) {
                THIS_BIN = ibin;
                break;
            }
        }
        // if ( THIS_BIN < 0 || THIS_BIN >= N_X_CENT_BINS ) { continue; }
        if ( THIS_BIN < 0 || THIS_BIN > N_X_CENT_BINS ) { std::cout << "THIS_BIN " << THIS_BIN << " xaxis_var " << xaxis_var << std::endl; continue; }
        
        float sum_et = tower_sum_energy_TOWERINFO_CALIB_CEMC + tower_sum_energy_TOWERINFO_CALIB_HCALIN + tower_sum_energy_TOWERINFO_CALIB_HCALOUT; 
        
        float nfired_cemc = N_CEMC_TOWERS*tower_frac_fired_TOWERINFO_CALIB_CEMC;
        float nfired_hcalin = N_HCALIN_TOWERS*tower_frac_fired_TOWERINFO_CALIB_HCALIN;
        float nfired_hcalout = N_HCALOUT_TOWERS*tower_frac_fired_TOWERINFO_CALIB_HCALOUT;
        float total_fired = nfired_cemc + nfired_hcalin + nfired_hcalout;

        float total_unmasked = ( N_CEMC_TOWERS*(1.0 - tower_frac_dead_TOWERINFO_CALIB_CEMC) 
            + N_HCALIN_TOWERS*(1.0 - tower_frac_dead_TOWERINFO_CALIB_HCALIN) 
            + N_HCALOUT_TOWERS*(1.0 - tower_frac_dead_TOWERINFO_CALIB_HCALOUT) );

        float total_towers = N_CEMC_TOWERS + N_HCALIN_TOWERS + N_HCALOUT_TOWERS;

        float avgE_cemc = tower_avg_energy_TOWERINFO_CALIB_CEMC;
        float avgE_hcalin = tower_avg_energy_TOWERINFO_CALIB_HCALIN;
        float avgE_hcalout = tower_avg_energy_TOWERINFO_CALIB_HCALOUT;
        
        float stdE_cemc = tower_std_energy_TOWERINFO_CALIB_CEMC;
        float stdE_hcalin = tower_std_energy_TOWERINFO_CALIB_HCALIN;
        float stdE_hcalout = tower_std_energy_TOWERINFO_CALIB_HCALOUT;
        
        float sum_et2_cemc = nfired_cemc*( (stdE_cemc*stdE_cemc) + (avgE_cemc*avgE_cemc) );
        float sum_et2_hcalin = nfired_hcalin*( (stdE_hcalin*stdE_hcalin) + (avgE_hcalin*avgE_hcalin) );
        float sum_et2_hcalout = nfired_hcalout*( (stdE_hcalout*stdE_hcalout) + (avgE_hcalout*avgE_hcalout) );
        float total_sum_et2 = sum_et2_cemc + sum_et2_hcalin + sum_et2_hcalout;

        int ncomp_cone_corr = random_cone_num_towers_RandomCones_r04 - random_cone_num_masked_towers_RandomCones_r04;
        int ncomp_cone_corr_sub1 = random_cone_num_towers_RandomCones_r04_Sub1 - random_cone_num_masked_towers_RandomCones_r04_Sub1;

        SUM_ET2[THIS_BIN] += total_sum_et2;
        SUM_E[THIS_BIN] += sum_et;
        TOTAL_TOWERS_NOT_MASKED[THIS_BIN] += total_unmasked;
        TOTAL_TOWERS_FIRED[THIS_BIN] += total_fired;
        N_TOWERS_TOTAL[THIS_BIN] += total_towers;
        SUM_CONE_COMPS[THIS_BIN] += ncomp_cone_corr;
        N_CONES_THIS_BIN[THIS_BIN]+=1.0;
        
    }

    TGraphErrors * g_poission = new TGraphErrors(N_X_CENT_BINS);
    TGraphErrors * g_poission_v2 = new TGraphErrors(N_X_CENT_BINS);
    TGraphErrors * g_poission_v3 = new TGraphErrors(N_X_CENT_BINS);

    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
       
        float x = h1_dummy->GetBinCenter(i+1);
        float xerr = h1_dummy->GetBinWidth(i+1)/2.0;

        float y = CalcPoisson(SUM_ET2[i], SUM_E[i], TOTAL_TOWERS_FIRED[i], SUM_CONE_COMPS[i], N_CONES_THIS_BIN[i]);
        float y_v2 = CalcPoisson(SUM_ET2[i], SUM_E[i], TOTAL_TOWERS_FIRED[i], SUM_CONE_COMPS[i], N_CONES_THIS_BIN[i], V2_VALUES[i]/100.0);
        float y_v3 = CalcPoisson(SUM_ET2[i], SUM_E[i], TOTAL_TOWERS_FIRED[i], SUM_CONE_COMPS[i], N_CONES_THIS_BIN[i], V2_VALUES[i]/100.0, V3_VALUES[i]/100.0);
       
        g_poission->SetPoint(i, x, y);
        g_poission->SetPointError(i, 0, 0);

        g_poission_v2->SetPoint(i, x, y_v2);
        g_poission_v2->SetPointError(i, 0, 0);

        g_poission_v3->SetPoint(i, x, y_v3);
        g_poission_v3->SetPointError(i, 0, 0);

    }
    g_poission->SetName("g_poission");
    g_poission_v2->SetName("g_poission_v2");
    g_poission_v3->SetName("g_poission_v3");
    std::vector< TGraphErrors * > g_poissions = { g_poission, g_poission_v2, g_poission_v3 };
    for ( auto g : g_poissions ) {
        g->GetXaxis()->SetTitle("Centrality [%]");
        g->GetYaxis()->SetTitle("#sigma(#delta E_{T}^{Raw})");
    }


    TFile * fout = new TFile(output_file.c_str(), "RECREATE");
    for ( auto g : g_poissions ) {
        g->Write();
    }
    fout->Close();
    f->Close();

    std::cout << "ProcessPoisson: Wrote to " << output_file << std::endl;

    return output_file;
    
}

std::string ProcessRhoTree( const std::string & input_file , const std::string & outfile, const std::string & prefix )
{
    std::string output_file = GetOutfile(outfile, output_dir, prefix);
    if ( !OVERRIDE ) {
        std::cout << "ProcessRhoTree: OVERRIDE is false, returning " << output_file << std::endl;
        return output_file;
    } else {
        std::cout << "ProcessRhoTree: " << output_file << std::endl;
    }

    TFile * f = new TFile(input_file.c_str(), "READ");
    if( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl;  exit(1); }


    TTree * t = (TTree*)f->Get("T");
    if ( !t ) { std::cout << "Tree not found in file " << input_file << std::endl; exit(1); }

    
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
    std::cout << "ProcessRhoTree: Processing " << nentries << " events" << std::endl;


    TH2F * h2_rho_area = new TH2F("h2_rho_area", "h2_rho_area", N_X_CENT_BINS, X_CENT_BINS, N_RHO_BINS, RHO_A_BINS);
    h2_rho_area->GetYaxis()->SetTitle("#rho_{A} [GeV]");
    TH2F * h2_rho_mult = new TH2F("h2_rho_mult", "h2_rho_mult", N_X_CENT_BINS, X_CENT_BINS, N_RHO_BINS, RHO_M_BINS);
    h2_rho_mult->GetYaxis()->SetTitle("#rho_{M} [GeV]");
    TH2F * h2_towerbackground = new TH2F("h2_towerbackground", "h2_towerbackground", N_X_CENT_BINS, X_CENT_BINS, N_RHO_BINS, TOWER_BACKGROUND_BINS);
    h2_towerbackground->GetYaxis()->SetTitle("E_{T,Bkgd}^{Iter}[GeV]");
    TH2F * h2_rho_times_A = new TH2F("h2_rho_times_A", "h2_rho_times_A", N_X_CENT_BINS, X_CENT_BINS, N_RHO_BINS, BACKGROUND_BINS);
    h2_rho_times_A->GetYaxis()->SetTitle("#rho_{A} #times A [GeV]");
    TH2F * h2_rho_times_N = new TH2F("h2_rho_times_N", "h2_rho_times_N", N_X_CENT_BINS, X_CENT_BINS, N_RHO_BINS, BACKGROUND_BINS);
    h2_rho_times_N->GetYaxis()->SetTitle("#rho_{M} #times N [GeV]");
    TH2F * h2_towerbackground_times_N = new TH2F("h2_towerbackground_times_N", "h2_towerbackground_times_N", N_X_CENT_BINS, X_CENT_BINS, N_RHO_BINS, BACKGROUND_BINS);
    h2_towerbackground_times_N->GetYaxis()->SetTitle("E_{T,Bkgd}^{Iter} #times N [GeV]");

    TH2F * h2_cone_avgN= new TH2F("h2_cone_avgN", "h2_cone_avgN", N_X_CENT_BINS, X_CENT_BINS, N_TOWERCOMP_BINS, TOWERCOMP_BINS);
    h2_cone_avgN->GetYaxis()->SetTitle("N_{cone}^{avg}");
    TH2F * h2_cone_avgN_sub1[3]; // 3 calo layers
    for ( int i = 0; i < 3; ++i ) {
        h2_cone_avgN_sub1[i] = new TH2F(Form("h2_cone_avgN_sub1_layer%d", i), Form("h2_cone_avgN_sub1_layer%d", i), N_X_CENT_BINS, X_CENT_BINS, N_TOWERCOMP_SUB1_BINS, TOWERCOMP_SUB1_BINS);
        h2_cone_avgN_sub1[i]->GetYaxis()->SetTitle(Form("N_{cone,sub1}^{avg}^{%d}", i));
    }
   
    std::vector< TH2F* > h2s = { h2_rho_area, h2_rho_mult, h2_towerbackground, 
                                h2_rho_times_A, h2_rho_times_N, h2_towerbackground_times_N, 
                                h2_cone_avgN, h2_cone_avgN_sub1[0], h2_cone_avgN_sub1[1], h2_cone_avgN_sub1[2] };
    for ( auto h2 : h2s ) {
        h2->GetXaxis()->SetTitle("Centrality [%]");
    }

    TH2F * h2_rho_area_vs_sumet = new TH2F("h2_rho_area_vs_sumet", "h2_rho_area_vs_sumet", N_SUM_ET_Q_BINS, SUM_ET_Q_BINS, N_RHO_BINS, RHO_A_BINS);
    h2_rho_area_vs_sumet->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
    h2_rho_area_vs_sumet->GetYaxis()->SetTitle("#rho_{A} [GeV]");
    TH2F * h2_rho_mult_vs_sumet = new TH2F("h2_rho_mult_vs_sumet", "h2_rho_mult_vs_sumet", N_SUM_ET_Q_BINS, SUM_ET_Q_BINS, N_RHO_BINS, RHO_M_BINS);
    h2_rho_mult_vs_sumet->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
    h2_rho_mult_vs_sumet->GetYaxis()->SetTitle("#rho_{M} [GeV]");
    TH2F * h2_towerbackground_vs_sumet = new TH2F("h2_towerbackground_vs_sumet", "h2_towerbackground_vs_sumet", N_SUM_ET_Q_BINS, SUM_ET_Q_BINS, N_RHO_BINS, TOWER_BACKGROUND_BINS);
    h2_towerbackground_vs_sumet->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
    h2_towerbackground_vs_sumet->GetYaxis()->SetTitle("E_{T,Bkgd}^{Iter}[GeV]");

    
    for ( int i = 0; i < nentries; ++i ) {

        t->GetEntry(i);
        float xaxis_var =1.0*centrality;
        if ( xaxis_var < X_CENT_BINS[0] || xaxis_var > X_CENT_BINS[N_X_CENT_BINS] ) { continue; }
        h2_cone_avgN->Fill(xaxis_var, random_cone_num_towers_RandomCones_r04);
        h2_cone_avgN_sub1[0]->Fill(xaxis_var, random_cone_num_towers_cemc_RandomCones_r04_Sub1);
        h2_cone_avgN_sub1[1]->Fill(xaxis_var, random_cone_num_towers_hcalin_RandomCones_r04_Sub1);
        h2_cone_avgN_sub1[2]->Fill(xaxis_var, random_cone_num_towers_hcalout_RandomCones_r04_Sub1);
    }

    TH1F * h1_cone_avgN = (TH1F*)h2_cone_avgN->ProfileX("h1_cone_avgN", 1, -1, "s");
    TH1F * h1_cone_avgN_sub1[3];
    for ( int i = 0; i < 3; ++i ) {
        h1_cone_avgN_sub1[i] = (TH1F*)h2_cone_avgN_sub1[i]->ProfileX(Form("h1_cone_avgN_sub1_layer%d", i), 1, -1, "s");
    }

    for ( int i = 0; i < nentries; ++i ) {

        t->GetEntry(i);
        float xaxis_var =1.0*centrality;
        if ( xaxis_var < X_CENT_BINS[0] || xaxis_var > X_CENT_BINS[N_X_CENT_BINS] ) { continue; }

        int THIS_BIN = -1;
        for ( int ibin = 0; ibin < N_X_CENT_BINS; ibin++ ) {
            if ( xaxis_var >= X_CENT_BINS[ibin] && xaxis_var <= X_CENT_BINS[ibin+1] ) {
                THIS_BIN = ibin;
                break;
            }
        }
        // if ( THIS_BIN < 0 || THIS_BIN >= N_X_CENT_BINS ) { continue; }
        if ( THIS_BIN < 0 || THIS_BIN > N_X_CENT_BINS ) { std::cout << "THIS_BIN " << THIS_BIN << " xaxis_var " << xaxis_var << std::endl; continue; }
        
        float avgN = h1_cone_avgN->GetBinContent(THIS_BIN+1);
        float avgN_cemc = h1_cone_avgN_sub1[0]->GetBinContent(THIS_BIN+1);
        float avgN_hcalin = h1_cone_avgN_sub1[1]->GetBinContent(THIS_BIN+1);
        float avgN_hcalout = h1_cone_avgN_sub1[2]->GetBinContent(THIS_BIN+1);
        // std::cout << "avgN " << avgN << " avgN_cemc " << avgN_cemc << " avgN_hcalin " << avgN_hcalin << " avgN_hcalout " << avgN_hcalout << std::endl;

        float sum_et = tower_sum_energy_TOWERINFO_CALIB_CEMC + tower_sum_energy_TOWERINFO_CALIB_HCALIN + tower_sum_energy_TOWERINFO_CALIB_HCALOUT;

        h2_rho_area->Fill(xaxis_var, rho_val_TowerRho_AREA);
        h2_rho_times_A->Fill(xaxis_var, rho_val_TowerRho_AREA*AREA_CONE);
        h2_rho_area_vs_sumet->Fill(sum_et, rho_val_TowerRho_AREA); 

        h2_rho_mult->Fill(xaxis_var, rho_val_TowerRho_MULT);
        h2_rho_mult_vs_sumet->Fill(sum_et, rho_val_TowerRho_MULT);
        h2_rho_times_N->Fill(xaxis_var, rho_val_TowerRho_MULT*avgN);
                
        
        float bkgd_avg = 0;
        float bkgd_avg_cemc = 0;
        float bkgd_avg_hcalin = 0;
        float bkgd_avg_hcalout = 0;
        for ( int j = 0; j < tower_background_energy_recemc->size(); ++j ) {
            float b_cemc = tower_background_energy_recemc->at(j);
            float b_hcalin = tower_background_energy_hcalin->at(j);
            float b_hcalout = tower_background_energy_hcalout->at(j);
            float b_tot = b_cemc + b_hcalin + b_hcalout;
            if ( std::isnan(b_tot) || std::isinf(b_tot) ) { continue; }
            bkgd_avg += b_tot;
            bkgd_avg_cemc += b_cemc;
            bkgd_avg_hcalin += b_hcalin;
            bkgd_avg_hcalout += b_hcalout;
        }
        bkgd_avg_cemc/=tower_background_energy_recemc->size();
        bkgd_avg_hcalin/=tower_background_energy_hcalin->size();
        bkgd_avg_hcalout/=tower_background_energy_hcalout->size();
        bkgd_avg_cemc*=avgN_cemc;
        bkgd_avg_hcalin*=avgN_hcalin;
        bkgd_avg_hcalout*=avgN_hcalout;
        float bkgd_avg_times_N = bkgd_avg_cemc + bkgd_avg_hcalin + bkgd_avg_hcalout;
        bkgd_avg/=24.0;

        h2_towerbackground->Fill(xaxis_var, bkgd_avg);
        h2_towerbackground_vs_sumet->Fill(sum_et, bkgd_avg);
        h2_towerbackground_times_N->Fill(xaxis_var, bkgd_avg_times_N);
       
    }


    TFile * fout = new TFile(output_file.c_str(), "RECREATE");
    fout->cd();
    for ( auto h2 : h2s ) {
        h2->Write();
    }
    h2_rho_area_vs_sumet->Write();
    h2_rho_mult_vs_sumet->Write();
    h2_towerbackground_vs_sumet->Write();
    h1_cone_avgN->Write();
    for ( int i = 0; i < 3; ++i ) {
        h1_cone_avgN_sub1[i]->Write();
    }
    fout->Close();
    f->Close();

    std::cout << "ProcessRhoTree: Wrote to " << output_file << std::endl;

    return output_file;
    
}

std::string ProcessWindowTree(const std::string & input_file, const std::string & outfile, const std::string & prefix  )
{

    
    std::string output_file = GetOutfile(outfile, output_dir, prefix);
    if ( !OVERRIDE ) {
        std::cout << "ProcessWindowTree: OVERRIDE is false, returning " << output_file << std::endl;
        return output_file;
    } else {
        std::cout << "ProcessWindowTree: " << output_file << std::endl;
    }

    TFile * f = new TFile(input_file.c_str(), "READ");
    if( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl;  exit(1); }

    TTree * t = (TTree*)f->Get("T");
    if ( !t ) { std::cout << "Tree not found in file " << input_file << std::endl; exit(1); }
    
    // tree branches 
        int centrality = 0;
        t->SetBranchAddress("centrality", &centrality);

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

        unsigned int max_window_vector_size = 0;  
        float avg_energy_full[11];
        float std_energy_full[11];
        int num_windows_full[11];
        
        t->SetBranchAddress("num_window_dims", &max_window_vector_size);
        t->SetBranchAddress("avg_energy_full", &avg_energy_full);
        t->SetBranchAddress("std_energy_full", &std_energy_full);
        t->SetBranchAddress("num_windows_full", &num_windows_full);

    int nentries = t->GetEntries();
    std::cout << "ProcessWindowTree: Processing " << nentries << " events" << std::endl;

    TH2F * h2_avg_et[k_window_array_size];
    TH2F * h2_sigma[k_window_array_size];
    TH2F * h2_nwindows[k_window_array_size];

    TProfile * p2_avg_et[k_window_array_size];
    TProfile * p2_sigma[k_window_array_size];
    TProfile * p2_avg_et2[k_window_array_size];
    TProfile * p2_nwindows[k_window_array_size];
    int NEVENTS_CENT[N_X_CENT_BINS_WINDOW];
    for ( int i = 0; i < N_X_CENT_BINS_WINDOW; ++i ) { NEVENTS_CENT[i] = 0; }

    for ( unsigned int i = 0; i < k_window_array_size; ++i ) {
    
        float AVG_BINS[N_WINDOW_BINS+1];
        float STD_BINS[N_WINDOW_BINS+1];
        float NUM_BINS[N_WINDOW_BINS+1];
    
        for ( int j = 0; j < N_WINDOW_BINS+1; ++j ) { 
            AVG_BINS[j] = j*MAX_WINDOW_AVG[i]/N_WINDOW_BINS; 
            STD_BINS[j] = j*MAX_WINDOW_STD[i]/N_WINDOW_BINS;
            NUM_BINS[j] = j*MAX_WINDOW_NUM[i]/N_WINDOW_BINS;
        }
        std::string window_str = Form("%dx%d",  k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second);
        
        h2_avg_et[i] = new TH2F(Form("h2_avg_et_%s", window_str.c_str()), Form("h2_avg_et_%s", window_str.c_str()), N_X_CENT_BINS_WINDOW, X_CENT_BINS_WINDOW, N_WINDOW_BINS, AVG_BINS);
        h2_avg_et[i]->GetYaxis()->SetTitle(Form("#LT E_{T}^{%dx%d} #GT [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));
        h2_avg_et[i]->GetXaxis()->SetTitle("Centrality [%]");
        h2_sigma[i] = new TH2F(Form("h2_sigma_%s", window_str.c_str()), Form("h2_sigma_%s", window_str.c_str()), N_X_CENT_BINS_WINDOW, X_CENT_BINS_WINDOW, N_WINDOW_BINS, STD_BINS);
        h2_sigma[i]->GetYaxis()->SetTitle(Form("#sigma^{%d #times %d} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));
        h2_sigma[i]->GetXaxis()->SetTitle("Centrality [%]");
        h2_nwindows[i] = new TH2F(Form("h2_nwindows_%s", window_str.c_str()), Form("h2_nwindows_%s", window_str.c_str()), N_X_CENT_BINS_WINDOW, X_CENT_BINS_WINDOW, N_WINDOW_BINS, NUM_BINS);
        h2_nwindows[i]->GetYaxis()->SetTitle(Form("N_{windows}^{%d #times %d}", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));
        h2_nwindows[i]->GetXaxis()->SetTitle("Centrality [%]");

        p2_avg_et[i] = new TProfile(Form("p2_avg_et_%s", window_str.c_str()), Form("p2_avg_et_%s", window_str.c_str()), N_X_CENT_BINS_WINDOW, X_CENT_BINS_WINDOW);
        p2_avg_et[i]->GetYaxis()->SetTitle(Form("#bar{#LT E_{T}^{%dx%d} #GT} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));
        p2_avg_et[i]->GetXaxis()->SetTitle("Centrality [%]");
        p2_sigma[i] = new TProfile(Form("p2_sigma_%s", window_str.c_str()), Form("p2_sigma_%s", window_str.c_str()), N_X_CENT_BINS_WINDOW, X_CENT_BINS_WINDOW);
        p2_sigma[i]->GetYaxis()->SetTitle(Form("#bar{#sigma^{%d #times %d}} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));
        p2_sigma[i]->GetXaxis()->SetTitle("Centrality [%]");
        p2_nwindows[i] = new TProfile(Form("p2_nwindows_%s", window_str.c_str()), Form("p2_nwindows_%s", window_str.c_str()), N_X_CENT_BINS_WINDOW, X_CENT_BINS_WINDOW);
        p2_nwindows[i]->GetYaxis()->SetTitle(Form("#bar{N_{windows}^{%d #times %d}}", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));
        p2_nwindows[i]->GetXaxis()->SetTitle("Centrality [%]");
        p2_avg_et2[i] = new TProfile(Form("p2_avg_et2_%s", window_str.c_str()), Form("p2_avg_et2_%s", window_str.c_str()), N_X_CENT_BINS_WINDOW, X_CENT_BINS_WINDOW);
        p2_avg_et2[i]->GetYaxis()->SetTitle(Form("#bar{#LT (E_{T}^{%dx%d})^{2} #GT} [GeV]", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second));
        p2_avg_et2[i]->GetXaxis()->SetTitle("Centrality [%]");

    }

    for ( int i = 0; i < nentries; ++i ) {
        t->GetEntry(i);
        float xaxis_var =1.0*centrality;
        if ( xaxis_var < X_CENT_BINS_WINDOW[0] || xaxis_var > X_CENT_BINS_WINDOW[N_X_CENT_BINS_WINDOW] ) { continue; }
        int THIS_BIN = -1;
        for ( int ibin = 0; ibin < N_X_CENT_BINS_WINDOW; ibin++ ) {
            if ( xaxis_var >= X_CENT_BINS_WINDOW[ibin] && xaxis_var <= X_CENT_BINS_WINDOW[ibin+1] ) {
                THIS_BIN = ibin;
                break;
            }
        }

        if ( THIS_BIN < 0 || THIS_BIN > N_X_CENT_BINS_WINDOW ) { std::cout << "THIS_BIN: " << THIS_BIN << std::endl; continue; }
        NEVENTS_CENT[THIS_BIN]++;

        for ( unsigned int iwindow = 0; iwindow < k_window_array_size; ++iwindow ) {
            float avgE = avg_energy_full[iwindow];
            float stdE = std_energy_full[iwindow];
            float numW = num_windows_full[iwindow];
            float avgE2 = numW*(stdE*stdE + avgE*avgE);
            h2_avg_et[iwindow]->Fill(xaxis_var, avgE);
            h2_sigma[iwindow]->Fill(xaxis_var, stdE);
            h2_nwindows[iwindow]->Fill(xaxis_var, numW);
            p2_avg_et[iwindow]->Fill(xaxis_var, avgE);
            p2_sigma[iwindow]->Fill(xaxis_var, stdE);
            p2_nwindows[iwindow]->Fill(xaxis_var, numW);
            p2_avg_et2[iwindow]->Fill(xaxis_var, avgE2);
        }
    }
    std::cout <<"TEST" << std::endl;

    TGraphErrors * g_sigma[N_X_CENT_BINS_WINDOW];
    TGraphErrors * g_avg_et[N_X_CENT_BINS_WINDOW];
    TGraphErrors * g_nwindows[N_X_CENT_BINS_WINDOW];
    TGraphErrors * g_avg_et2[N_X_CENT_BINS_WINDOW];
    for ( int icent = 0; icent < N_X_CENT_BINS_WINDOW; ++icent ) {
        g_sigma[icent] = new TGraphErrors(k_window_array_size);
        g_avg_et[icent] = new TGraphErrors(k_window_array_size);
        g_nwindows[icent] = new TGraphErrors(k_window_array_size);
        g_avg_et2[icent] = new TGraphErrors(k_window_array_size);
        g_sigma[icent]->SetName(Form("g_sigma_cent%d", icent));
        g_avg_et[icent]->SetName(Form("g_avg_et_cent%d", icent));
        g_nwindows[icent]->SetName(Form("g_nwindows_cent%d", icent));
        g_avg_et2[icent]->SetName(Form("g_avg_et2_cent%d", icent));
        for ( unsigned int jdim = 0; jdim < k_window_array_size; ++jdim ) {
            float area = 1.0*k_calo_window_dims_hcal_geom[jdim].first*k_calo_window_dims_hcal_geom[jdim].second;
            h2_avg_et[jdim]->GetXaxis()->SetRangeUser(X_CENT_BINS_WINDOW[icent], X_CENT_BINS_WINDOW[icent+1]);
            h2_sigma[jdim]->GetXaxis()->SetRangeUser(X_CENT_BINS_WINDOW[icent], X_CENT_BINS_WINDOW[icent+1]);
            h2_nwindows[jdim]->GetXaxis()->SetRangeUser(X_CENT_BINS_WINDOW[icent], X_CENT_BINS_WINDOW[icent+1]);
            // h2_avg_et2[jdim]->GetXaxis()->SetRangeUser(X_CENT_BINS_WINDOW[icent], X_CENT_BINS_WINDOW[icent+1]);
            TH1D * h1_avg_et = (TH1D *) h2_avg_et[jdim]->ProjectionY();
            TH1D * h1_sigma = (TH1D*)h2_sigma[jdim]->ProjectionY();
            TH1D * h1_nwindows = (TH1D*)h2_nwindows[jdim]->ProjectionY();
            // TH1D * h1_avg_et2 =(TH1D*) h2_avg_et2[jdim]->ProjectionY();
            float avg_et = h1_avg_et->GetMean();
            float avg_et_err = h1_avg_et->GetMeanError();
            // float sigma = h1_sigma->GetMean();
            // float sigma_err = h1_sigma->GetMeanError();
            float nwindows = h1_nwindows->GetMean();
            float nwindows_err = h1_nwindows->GetMeanError();
            // float avg_et2 = h1_avg_et2->GetMean();
            // float avg_et2_err = h1_avg_et2->GetMeanError();
            


            // float avg_et = p2_avg_et[jdim]->GetBinContent(icent+1);
            // float avg_et_err = p2_avg_et[jdim]->GetBinError(icent+1);
            float sigma = p2_sigma[jdim]->GetBinContent(icent+1);
            float sigma_err = p2_sigma[jdim]->GetBinError(icent+1);
            // float nwindows = p2_nwindows[jdim]->GetBinContent(icent+1);
            // float nwindows_err = p2_nwindows[jdim]->GetBinError(icent+1);
            float avg_et2 = p2_avg_et2[jdim]->GetBinContent(icent+1);
            float avg_et2_err = p2_avg_et2[jdim]->GetBinError(icent+1);
            g_sigma[icent]->SetPoint(jdim, area, sigma);
            g_sigma[icent]->SetPointError(jdim, 0, sigma_err);
            g_avg_et[icent]->SetPoint(jdim, area, avg_et);
            g_avg_et[icent]->SetPointError(jdim, 0, avg_et_err);
            g_nwindows[icent]->SetPoint(jdim, area, nwindows);
            g_nwindows[icent]->SetPointError(jdim, 0, nwindows_err);
            g_avg_et2[icent]->SetPoint(jdim, area, avg_et2);
            g_avg_et2[icent]->SetPointError(jdim, 0, avg_et2_err);
        }
    }

   
    TFile * fout1 = new TFile(output_file.c_str(), "RECREATE");
    fout1->cd();
    for ( unsigned int i = 0; i < k_window_array_size; ++i ) {
        h2_avg_et[i]->Write();
        h2_sigma[i]->Write();
        h2_nwindows[i]->Write();
        p2_avg_et[i]->Write();
        p2_sigma[i]->Write();
        p2_nwindows[i]->Write();
        p2_avg_et2[i]->Write();
    }
    for ( int icent = 0; icent < N_X_CENT_BINS_WINDOW; ++icent ) {
        g_sigma[icent]->Write();
        g_avg_et[icent]->Write();
        g_nwindows[icent]->Write();
        g_avg_et2[icent]->Write();
    }
    fout1->Close();

    std::cout << "ProcessWindowTree: Wrote to " << output_file << std::endl;

    return output_file;
    
}

void RhoPlots(const std::string input_file, const std::string & prefix)
{

    std::string outdir = MakeGetDir(output_dir+"/"+prefix);
    TFile * f = new TFile(input_file.c_str(), "READ");
    if ( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl; exit(1); }

    TH2F * h2_rho_area = (TH2F*)f->Get("h2_rho_area");
    TH2F * h2_rho_mult = (TH2F*)f->Get("h2_rho_mult");
    TH2F * h2_towerbackground = (TH2F*)f->Get("h2_towerbackground");
    TH2F * h2_rho_times_A = (TH2F*)f->Get("h2_rho_times_A");
    TH2F * h2_rho_times_N = (TH2F*)f->Get("h2_rho_times_N");
    TH2F * h2_towerbackground_times_N = (TH2F*)f->Get("h2_towerbackground_times_N");
    TH2F * h2_rho_area_vs_sumet = (TH2F*)f->Get("h2_rho_area_vs_sumet");
    TH2F * h2_rho_mult_vs_sumet = (TH2F*)f->Get("h2_rho_mult_vs_sumet");
    TH2F * h2_towerbackground_vs_sumet = (TH2F*)f->Get("h2_towerbackground_vs_sumet");

    if ( !h2_rho_area || !h2_rho_mult || !h2_towerbackground 
        || !h2_rho_times_A || !h2_rho_times_N || !h2_towerbackground_times_N 
        || !h2_rho_area_vs_sumet || !h2_rho_mult_vs_sumet || !h2_towerbackground_vs_sumet ) {
             std::cout << "RhoPlots::Histograms not found" << std::endl; exit(1);
    }

    TCanvas * c;

    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    double tx=0.19;
    double ty_start=0.85;
    std::vector<TH2F*> h2s = {h2_rho_area, h2_rho_mult, h2_towerbackground, 
                            h2_rho_times_A, h2_rho_times_N, h2_towerbackground_times_N, 
                            h2_rho_area_vs_sumet, h2_rho_mult_vs_sumet, h2_towerbackground_vs_sumet};
    std::vector<float> txs = {0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.19, 0.19, 0.19};
    std::vector<std::string> names = {"rho_area", "rho_mult", "towerbackground", 
                            "rho_times_A", "rho_times_N", "towerbackground_times_N", 
                            "rho_area_vs_sumet", "rho_mult_vs_sumet", "towerbackground_vs_sumet"};

    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};

    for ( unsigned int i = 0; i < h2s.size(); ++i ) {
        c = new TCanvas("c", "c", 800, 600);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogz();
        h2s[i]->GetXaxis()->SetNdivisions(505);
        h2s[i]->GetYaxis()->SetNdivisions(505);
        h2s[i]->Draw("colz");
        float ty = ty_start;
        tx = txs[i];
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.05;
        } 
        c->SaveAs(Form("%s/%s.png", outdir.c_str(), names[i].c_str()));
        delete c;
    }

    TLegend * leg = new TLegend(0.2,0.72,0.46,0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetNColumns(2);

    h2s = {h2_rho_times_A, h2_rho_times_N};
    names = {"rho_times_A_centslices", "rho_times_N_centslices"};
    double miny = 1e-3;
    int ihist = 0;
    tx = 0.48;
    for ( auto h2 : h2s ) {
        double ty = ty_start;
        c = new TCanvas("c", "c", 800, 600);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        float maxx= 0;
        for ( int j = 0; j < N_X_CENT_BINS; ++j ) {
            h2->GetXaxis()->SetRange(j+1, j+1);
            TH1D * h1 = h2->ProjectionY(Form("h1_%s_%d", h2->GetTitle(), j));
            h1->SetLineColor(COLORS[j]);
            h1->SetMarkerColor(COLORS[j]);
            h1->SetMarkerStyle(MARKERS[j]);
            h1->Scale(1.0/h1->Integral(), "width");

            std::string leg_title = Form("%d-%d %%", int(X_CENT_BINS[j]), int(X_CENT_BINS[j+1]));
            // std::string leg_title = Form("%d-%d %%", int(X_CENT_BINS[j]), int(X_CENT_BINS[j+1]));
            leg->AddEntry(h1, leg_title.c_str(), "lp");
            h1->GetYaxis()->SetRangeUser(miny, 2e1);
            int lastbin_above_threshold = 0;
            lastbin_above_threshold = h1->FindLastBinAbove(miny);
            if(lastbin_above_threshold < 0) { lastbin_above_threshold = h1->GetNbinsX();}
            float maxx = 1.2*h1->GetBinCenter(lastbin_above_threshold);
            if ( maxx > h1->GetXaxis()->GetXmax() ) { maxx = h1->GetXaxis()->GetXmax(); }
            h1->GetXaxis()->SetRangeUser(0, maxx);
            h1->GetXaxis()->SetTitle(h2->GetYaxis()->GetTitle());
            h1->GetYaxis()->SetTitle("Probabilty Density [A.U.]");
            
            if ( j == 0 ) {
                h1->Draw("P");
            } else {
                h1->Draw("SAME");
            }
        }

       leg->Draw("SAME");
       
       for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.07;
        }

        c->SaveAs(Form("%s/%s.png", outdir.c_str(), names[ihist].c_str()));
        ihist++;
        delete c;
        leg->Clear();
    }
    delete leg;


    TLegend * leg_mult = new TLegend(0.2,0.72,0.46,0.92);
    leg_mult->SetBorderSize(0);
    leg_mult->SetFillStyle(0);
    leg_mult->SetTextSize(0.04);
    leg_mult->SetNColumns(1);
    leg_mult->SetHeader("#rho_{M} #times N");

    TLegend * leg_area = new TLegend(0.18,0.72,0.5,0.92);
    leg_area->SetBorderSize(0);
    leg_area->SetFillStyle(0);
    // leg_area->SetTextSize(0.04);
    leg_area->SetNColumns(2);
    // leg_area->SetHeader("#rho_{A} \times A");

    const int MY_COLORS[] = {kBlack, kRed,  kOrange+2, kAzure-2, kGreen+2, kViolet, kOrange+2,kAzure-1,  kMagenta+2, kRed+2};
    const int MY_FULL_MARKERS[] = { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross};
    const int MY_OPEN_MARKERS[] = { kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown, kOpenDiamond, kOpenCross, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown, kOpenDiamond, kOpenCross};
    const float MY_MARKER_SIZE = 1.2;
    const float MY_LINE_WIDTH = 1.5;


    h2s = {h2_rho_times_A, h2_rho_times_N};
    names = {"rho_times_A_centslices", "rho_times_N_centslices"};
    std::vector<std::string> hist_labs = {"#rho_{A} #times A", "#rho_{M} #times N"};
    miny = 1e-4;
    ihist = 0;
    tx = 0.53;
    double ty = ty_start;
    c = new TCanvas("c", "c", 800, 600);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    float maxx= 0;
    for ( int j = 0; j < N_X_CENT_BINS-3; ++j ) {
        for ( unsigned int ihist = 0; ihist < h2s.size(); ++ihist ) {
            TH2F * h2 = h2s[ihist];
            h2->GetXaxis()->SetRange(j+1, j+1);
            TH1D * h1 = h2->ProjectionY(Form("h1_%s_%d", h2->GetTitle(), j));
            h1->Scale(1.0/h1->Integral(), "width");
            
            
            
            // double mean = hd1->GetMean();
            // double mean_err = hd1->GetMeanError();
            // const int N_points = hd1->GetNbinsX();

            // TGraphErrors *h1;
            // std::vector<double> xv {};
            // std::vector<double> yv {};
            // std::vector<double> y_errv {};
            // std::vector<double> y_err_newv {};
            // for ( int ibin = 0; ibin < N_points-1; ++ibin ) {
            //     double x = hd1->GetBinCenter(ibin+1);
            //     double y = hd1->GetBinContent(ibin+1);
            //     double y0 = y;
            //     double y_err = hd1->GetBinError(ibin+1);
            //     // subtract mean
            //     y -= mean;
            //     double b = TMath::Sqrt(mean_err*mean_err + y_err*y_err);

            //     y /= mean;
            //     double y_err_new = TMath::Sqrt(((y_err*y_err)/(y0*y0)) + ((mean_err*mean_err)/(mean*mean)));
            //     // h1->SetPoint(ibin, x, y);
            //     // h1->SetPointError(ibin, 0, y_err_new);
            //     xv.push_back(x);
            //     yv.push_back(y);
            //     y_errv.push_back(y_err);
            //     y_err_newv.push_back(y_err_new);


            // }
            // std::cout << "N_points: " << N_points << std::endl;
            // h1 = new TGraphErrors(N_points-1, &xv[0], &yv[0], 0, &y_err_newv[0]);



            h1->SetName(Form("h1_%s_%d", h2->GetTitle(), j)); 
            
            h1->SetLineColor(MY_COLORS[j]);
            h1->SetMarkerColor(MY_COLORS[j]);
            // h1->SetMarkerSize(MY_MARKER_SIZE);
            h1->SetLineWidth(MY_LINE_WIDTH);
            if (ihist == 0) {
                h1->SetMarkerStyle(MY_FULL_MARKERS[j]);
            } else {
                h1->SetMarkerStyle(MY_OPEN_MARKERS[j]);
            }
            // h1->Scale(1.0/h1->Integral(), "width");

            std::string leg_title = Form("%d-%d%%", int(X_CENT_BINS[j]), int(X_CENT_BINS[j+1]));
            if (ihist == 0) {
                leg_area->AddEntry(h1, leg_title.c_str(), "lp");
            } else {
                leg_mult->AddEntry(h1, leg_title.c_str(), "lp");
            }
            h1->GetYaxis()->SetRangeUser(miny, 2e1);
            int lastbin_above_threshold = 0;
            lastbin_above_threshold = h1->FindLastBinAbove(miny);
            if(lastbin_above_threshold < 0) { lastbin_above_threshold = h1->GetNbinsX();}
            float maxx = 1.2*h1->GetBinCenter(lastbin_above_threshold);
            if ( maxx > h1->GetXaxis()->GetXmax() ) { maxx = h1->GetXaxis()->GetXmax(); }
            h1->GetXaxis()->SetRangeUser(0, maxx);
            h1->GetXaxis()->SetTitle(h2->GetYaxis()->GetTitle());
            h1->GetYaxis()->SetTitle("Probabilty Density [A.U.]");
            h1->GetXaxis()->SetTitle("#rho_{A} #times A [GeV], #rho_{M} #times N [GeV]");
            if ( j == 0 && ihist == 0 ) {
                h1->Draw("PE");
            } else {
                h1->Draw("SAME PE");
            }
        }
    }

    leg_area->Draw("SAME");
    // leg_mult->Draw("SAME");

    for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.07;
        }

    TLatex * tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.038);
    tags = {"#it{Closed Points: #rho_{A} #times A}", "#it{Open Points: #rho_{M} #times N}"};
    for ( auto tag : tags ) {
        tex2->DrawLatex(tx, ty, tag.c_str());
        ty -= 0.07;
    }
    c->SaveAs(Form("%s/rho_comp.png", outdir.c_str()));
    // ihist++;
    delete c;
    // leg_area->Clear();
    // leg_mult->Clear();
    delete leg_area;
    delete leg_mult;
    f->Close();

    return;

}

void MultCurves(const std::string input_file, const std::string & prefix)
{

    std::string outdir = MakeGetDir(output_dir+"/"+prefix);
    TFile * f = new TFile(input_file.c_str(), "READ");
    if ( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl; exit(1); }


    TH1F * h1_et_vs_ntruth[N_X_CENT_BINS];
    TH2F * h2_et_vs_ntruth[N_X_CENT_BINS];
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
        h2_et_vs_ntruth[i] = (TH2F*)f->Get(Form("h2_et_vs_ntruth_%d", i));
        if ( !h2_et_vs_ntruth[i] ) { std::cout << "h2_et_vs_ntruth_" << i << " not found" << std::endl; exit(1); }
        h1_et_vs_ntruth[i] = (TH1F*)f->Get(Form("h1_et_vs_ntruth_%d", i));
        if ( !h1_et_vs_ntruth[i] ) { std::cout << "h1_et_vs_ntruth_" << i << " not found" << std::endl; exit(1); }
        h1_et_vs_ntruth[i]->GetXaxis()->SetNdivisions(505);
        h1_et_vs_ntruth[i]->GetYaxis()->SetNdivisions(505);
        h1_et_vs_ntruth[i]->SetTitle(Form("%d-%d %%", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])));
        h1_et_vs_ntruth[i]->SetLineColor(COLORS[i]);
        h1_et_vs_ntruth[i]->SetMarkerColor(COLORS[i]);
        h1_et_vs_ntruth[i]->SetMarkerStyle(MARKERS[i]);
        h1_et_vs_ntruth[i]->GetXaxis()->SetTitle("E_{T,jet}^{Raw} [GeV]");
        h1_et_vs_ntruth[i]->GetYaxis()->SetTitle("#LT N_{comp}^{truth} #GT");
    }
    const float X_COURSE_CENT_BINS[]= {0, 10, 20, 40, 60, 80};
    const int N_X_CENT_BINS = sizeof(X_COURSE_CENT_BINS)/sizeof(X_COURSE_CENT_BINS[0]) - 1;
    
    // merge h2s depending on centrality
    TH2F * h2_et_vs_ntruth_course[N_X_CENT_BINS];
    h2_et_vs_ntruth_course[0] = (TH2F*)h2_et_vs_ntruth[0]->Clone("h2_et_vs_ntruth_course_0"); // 0-5
    h2_et_vs_ntruth_course[0]->Add(h2_et_vs_ntruth[1]); // 5-10
    h2_et_vs_ntruth_course[1] = (TH2F*)h2_et_vs_ntruth[2]->Clone("h2_et_vs_ntruth_course_1"); //10-20
    h2_et_vs_ntruth_course[2] = (TH2F*)h2_et_vs_ntruth[3]->Clone("h2_et_vs_ntruth_course_2"); //20-30
    h2_et_vs_ntruth_course[2]->Add(h2_et_vs_ntruth[4]); //30-40
    h2_et_vs_ntruth_course[3] = (TH2F*)h2_et_vs_ntruth[5]->Clone("h2_et_vs_ntruth_course_3"); //40-50
    h2_et_vs_ntruth_course[3]->Add(h2_et_vs_ntruth[6]); //50-60
    h2_et_vs_ntruth_course[4] = (TH2F*)h2_et_vs_ntruth[7]->Clone("h2_et_vs_ntruth_course_4"); //60-70
    h2_et_vs_ntruth_course[4]->Add(h2_et_vs_ntruth[8]); //70-80

    // for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
    //     int icourse = -1;
    //     for ( int j = 0; j < N_X_CENT_BINS; ++j ) {
    //         if ( X_CENT_BINS[i] >= X_COURSE_CENT_BINS[j] && X_CENT_BINS[i+1] <= X_COURSE_CENT_BINS[j+1] ) {
    //             icourse = j;
    //             break;
    //         }
    //     }
    //     if ( icourse < 0 || icourse >= N_X_CENT_BINS ) { std::cout << "icourse: " << icourse << std::endl; exit(1); }
    //     if ( !h2_et_vs_ntruth_course[icourse] ) {
    //         h2_et_vs_ntruth_course[icourse] = (TH2F*)h2_et_vs_ntruth[i]->Clone(Form("h2_et_vs_ntruth_%d", icourse));
    //     } else {
    //         h2_et_vs_ntruth_course[icourse]->Add(h2_et_vs_ntruth[i]);
    //     }
    

    // }
    float minx = 0, maxx = -1e9;
    float miny = 0, maxy = -1e9;
    TH1F * h1_et_vs_ntruth_course[N_X_CENT_BINS];
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
        h1_et_vs_ntruth_course[i] = (TH1F*)h2_et_vs_ntruth_course[i]->ProfileY(Form("h1_et_vs_ntruth_%d_couse", i), 1, -1, "e");
        if ( !h1_et_vs_ntruth_course[i] ) { std::cout << "h1_et_vs_ntruth_" << i << " not found" << std::endl; exit(1); }
        h1_et_vs_ntruth_course[i]->GetXaxis()->SetNdivisions(505);
        h1_et_vs_ntruth_course[i]->GetYaxis()->SetNdivisions(505);
        h1_et_vs_ntruth_course[i]->SetTitle(Form("%d-%d%%", int(X_COURSE_CENT_BINS[i]), int(X_COURSE_CENT_BINS[i+1])));
        h1_et_vs_ntruth_course[i]->SetLineColor(COLORS[i]);
        h1_et_vs_ntruth_course[i]->SetMarkerColor(COLORS[i]);
        h1_et_vs_ntruth_course[i]->SetMarkerStyle(MARKERS[i]);
        h1_et_vs_ntruth_course[i]->GetXaxis()->SetTitle("E_{T,jet}^{Uncorr.} [GeV]");
        h1_et_vs_ntruth_course[i]->GetYaxis()->SetTitle("#LT N_{Tower}^{Sim.} #GT");
        int last_nonzero_bin = h1_et_vs_ntruth_course[i]->FindLastBinAbove(0);
        float x = h1_et_vs_ntruth_course[i]->GetBinCenter(last_nonzero_bin);
        if ( x > maxx ) { maxx = x; }
        if ( h1_et_vs_ntruth_course[i]->GetMaximum() > maxy ) { maxy = h1_et_vs_ntruth_course[i]->GetMaximum(); }
    }
    
    // for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
    //     h1_et_vs_ntruth[i] = (TH1F*)f->Get(Form("h1_et_vs_ntruth_%d", i));
    //     if ( !h1_et_vs_ntruth[i] ) { std::cout << "h1_et_vs_ntruth_" << i << " not found" << std::endl; exit(1); }
    //     h1_et_vs_ntruth[i]->GetXaxis()->SetNdivisions(505);
    //     h1_et_vs_ntruth[i]->GetYaxis()->SetNdivisions(505);
    //     h1_et_vs_ntruth[i]->SetTitle(Form("%d-%d %%", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])));
    //     h1_et_vs_ntruth[i]->SetLineColor(COLORS[i]);
    //     h1_et_vs_ntruth[i]->SetMarkerColor(COLORS[i]);
    //     h1_et_vs_ntruth[i]->SetMarkerStyle(MARKERS[i]);
    //     h1_et_vs_ntruth[i]->GetXaxis()->SetTitle("E_{T,jet}^{Raw} [GeV]");
    //     h1_et_vs_ntruth[i]->GetYaxis()->SetTitle("#LT N_{comp}^{truth} #GT");
    //     int last_nonzero_bin = h1_et_vs_ntruth[i]->FindLastBinAbove(0);
    //     float x = h1_et_vs_ntruth[i]->GetBinCenter(last_nonzero_bin);
    //     if ( x > maxx ) { maxx = x; }
    //     if ( h1_et_vs_ntruth[i]->GetMaximum() > maxy ) { maxy = h1_et_vs_ntruth[i]->GetMaximum(); }
    // }

    TCanvas * c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);

    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    double tx=0.19;
    double ty_start=0.3;
    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag, "PYTHIA Dijet Embedded"};

    TLegend * leg = new TLegend(0.2,0.7,0.55,0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(2);

    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
        
        h1_et_vs_ntruth_course[i]->GetXaxis()->SetRangeUser(minx, 1.1*maxx);
        h1_et_vs_ntruth_course[i]->GetYaxis()->SetRangeUser(miny,1.5*maxy);
        if ( i == 0 ) {
            h1_et_vs_ntruth_course[i]->Draw("P");
        } else {
            h1_et_vs_ntruth_course[i]->Draw("P SAME");
        }
        leg->AddEntry(h1_et_vs_ntruth_course[i], h1_et_vs_ntruth_course[i]->GetTitle(), "lp");
    }

    // for ( int i = 0; i < N_X_CENT_BINS-3; ++i ) {
        
    //     h1_et_vs_ntruth[i]->GetXaxis()->SetRangeUser(minx, 1.1*maxx);
    //     h1_et_vs_ntruth[i]->GetYaxis()->SetRangeUser(miny,1.3*maxy);
    //     if ( i == 0 ) {
    //         h1_et_vs_ntruth[i]->Draw("P");
    //     } else {
    //         h1_et_vs_ntruth[i]->Draw("P SAME");
    //     }
    //     leg->AddEntry(h1_et_vs_ntruth[i], h1_et_vs_ntruth[i]->GetTitle(), "lp");
    // }
    // N_X_CENT_BINS
    leg->Draw("SAME");
    for ( auto tag : tags ) {
        tex->DrawLatex(tx, ty_start, tag.c_str());
        ty_start -= 0.05;
    }

    c->SaveAs(Form("%s/nsignal_curves.png", outdir.c_str()));

    delete c;
    delete leg;

    f->Close();

    return;

}

void WindowFits(const std::string input_file, const std::string & prefix)
{

    std::string outdir = MakeGetDir(output_dir+"/"+prefix);
    TFile * f = new TFile(input_file.c_str(), "READ");
    if ( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl; exit(1); }

   
    TGraphErrors * g_sigma[N_X_CENT_BINS_WINDOW];

    float fit_p0_full[N_X_CENT_BINS_WINDOW];
    float fit_p0err_full[N_X_CENT_BINS_WINDOW];
    float fit_chi2_full[N_X_CENT_BINS_WINDOW];

    float fit_p0_tight[N_X_CENT_BINS_WINDOW];
    float fit_p0err_tight[N_X_CENT_BINS_WINDOW];
    float fit_chi2_tight[N_X_CENT_BINS_WINDOW];

    float fit_p0_avg[N_X_CENT_BINS_WINDOW];
    float fit_p0err_avg[N_X_CENT_BINS_WINDOW];

    const float fit_range_full[] = {0, 170};
    const float fit_range_tight[] = {1, 50};

    float fit_p0_ratio[N_X_CENT_BINS_WINDOW];
    float fit_p0_sys[N_X_CENT_BINS_WINDOW];


    TF1 * fit_sigma;
    for ( int i = 0; i < N_X_CENT_BINS_WINDOW; ++i ) {

        g_sigma[i] = (TGraphErrors*)f->Get(Form("g_sigma_cent%d", i));
        if ( !g_sigma[i] ) {
            std::cout << "WindowFits::Graphs not found" << std::endl; exit(1);
        }  
        std::vector<TGraphErrors*> graphs = {g_sigma[i]};
        
        for ( auto g : graphs ) {
            g->GetXaxis()->SetNdivisions(505);
            g->GetYaxis()->SetNdivisions(505);
            g->GetXaxis()->SetTitle("A^{nxm} / A^{1x1}");
            g->SetMarkerStyle(20);
            g->SetMarkerSize(1.5);
            g->SetLineColor(kBlack);
            g->SetMarkerColor(kBlack);

        }

        g_sigma[i]->GetYaxis()->SetTitle("#bar{#sigma} / #bar{#sigma}^{1x1}");
        float x0 = g_sigma[i]->GetY()[0];
        float x0err = g_sigma[i]->GetEY()[0];

        for ( int j = 0; j < k_window_array_size; ++j ) {
            float yerr = g_sigma[i]->GetEY()[j];
            float y = g_sigma[i]->GetY()[j]/x0;
            float yerr2 = y*TMath::Sqrt(TMath::Power(yerr/g_sigma[i]->GetY()[j], 2) + TMath::Power(x0err/x0, 2));
            g_sigma[i]->SetPointError(j, 0, yerr2);
            g_sigma[i]->SetPoint(j, g_sigma[i]->GetX()[j], y);
        }
       
        TGraphErrors * g_sigma_fit_copy = (TGraphErrors*)g_sigma[i]->Clone(Form("g_sigma_fit_copy_%d", i));

        fit_sigma = new  TF1("fit_sigma", "x^[0]", 0, 500);
        g_sigma_fit_copy->Fit(fit_sigma, "MRQ", "", fit_range_full[0], fit_range_full[1]);
        fit_p0_full[i] = fit_sigma->GetParameter(0);
        fit_p0err_full[i] = fit_sigma->GetParError(0);
        fit_chi2_full[i] = fit_sigma->GetChisquare()/fit_sigma->GetNDF();

        g_sigma_fit_copy->Fit(fit_sigma, "MRQ", "", fit_range_tight[0], fit_range_tight[1]);
        fit_p0_tight[i] = fit_sigma->GetParameter(0);
        fit_p0err_tight[i] = fit_sigma->GetParError(0);
        fit_chi2_tight[i] = fit_sigma->GetChisquare()/fit_sigma->GetNDF();

        if( fit_p0_full[i] > 0) {
            fit_p0_ratio[i] = fit_p0_tight[i] / fit_p0_full[i];
        } else {
            std::cout << "Warning: fit_p0_full[" << i << "] is zero or negative, setting ratio to 1.0" << std::endl;
            fit_p0_ratio[i] = 1.0;
        }

        fit_p0_sys[i] = (1 - fit_p0_ratio[i])*fit_p0_full[i];

        fit_p0_avg[i] = (fit_p0_full[i] + fit_p0_tight[i]) / 2.0;
        fit_p0err_avg[i] = TMath::Sqrt(TMath::Power(fit_p0err_full[i], 2) + TMath::Power(fit_p0err_tight[i], 2)) / 2.0;

        std::cout << "Cent: " << i 
                  << ", Full Fit: k_fit = " << fit_p0_full[i] 
                  << " +/- " << fit_p0err_full[i] 
                  << ", Chi2/ndf = " << fit_chi2_full[i] 
                  << ", Tight Fit: k_fit = " << fit_p0_tight[i] 
                  << " +/- " << fit_p0err_tight[i] 
                  << ", Chi2/ndf = " << fit_chi2_tight[i]
                  << ", Ratio = " << fit_p0_ratio[i] 
                  << ", Systematic Error = " << fit_p0_sys[i]
                    << ", Average Fit: k_fit = " << fit_p0_avg[i]
                    << " +/- " << fit_p0err_avg[i]
                  << std::endl;
        
        delete g_sigma_fit_copy;

    }

    TCanvas * c;
    TLegend * leg;
    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    double tx=0.19;
    double ty_start=0.85;
    double tx_f = 0.55;
    double ty_f_start = 0.35;
    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};

    for ( int i = 0; i < N_X_CENT_BINS_WINDOW; ++i ) {
        std::string cent_title = Form("%d-%d%% Central", int(X_CENT_BINS_WINDOW[i]), int(X_CENT_BINS_WINDOW[i+1]));
        c = new TCanvas("c", "c", 800, 600);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
       
        g_sigma[i]->Draw("AP");
        g_sigma[i]->GetXaxis()->SetRangeUser(0, 1.1*g_sigma[i]->GetX()[g_sigma[i]->GetN()-1]);
        g_sigma[i]->GetYaxis()->SetRangeUser(0, 1.3*g_sigma[i]->GetY()[g_sigma[i]->GetN()-1]);

        fit_sigma = new TF1("fit_sigma", "x^[0]", fit_range_full[i], fit_range_full[i]);
        fit_sigma->SetParameter(0, fit_p0_full[i]);
        fit_sigma->SetLineColor(kRed);
        fit_sigma->SetLineStyle(2);
        fit_sigma->SetLineWidth(2);
        fit_sigma->Draw("SAME");

        float ty = ty_start;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.06;
        }
        tex->DrawLatex(tx, ty, cent_title.c_str());
        leg = new TLegend(0.2,0.7,0.55,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        tex->DrawLatex(tx_f, ty_f_start, Form("k_{fit} = %0.2f", fit_p0_full[i]));
        c->SaveAs(Form("%s/sigma_fit_cent%d.png", outdir.c_str(), i));
        delete c;
    }

   
    c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    TGraphErrors * g_sigma_fit = new TGraphErrors(N_X_CENT_BINS_WINDOW);
    TGraphAsymmErrors  * g_sigma_fit_sys = new TGraphAsymmErrors (N_X_CENT_BINS_WINDOW);
    for ( int i = 0; i < N_X_CENT_BINS_WINDOW; ++i ) {
        g_sigma_fit->SetPoint(i, X_CENT_BINS_WINDOW[i+1]-0.5*(X_CENT_BINS_WINDOW[i+1]-X_CENT_BINS_WINDOW[i]), fit_p0_full[i]);
        g_sigma_fit->SetPointError(i, 0, fit_p0err_full[i]);
        g_sigma_fit_sys->SetPoint(i, X_CENT_BINS_WINDOW[i+1]-0.5*(X_CENT_BINS_WINDOW[i+1]-X_CENT_BINS_WINDOW[i]), fit_p0_full[i]);
        g_sigma_fit_sys->SetPointEYlow(i, fit_p0_sys[i]); // 
    }

    g_sigma_fit->GetXaxis()->SetNdivisions(505);
    g_sigma_fit->GetYaxis()->SetNdivisions(505);
    g_sigma_fit->GetXaxis()->SetTitle("Centrality [%]");
    g_sigma_fit->GetYaxis()->SetTitle("k_{fit}");
    g_sigma_fit->GetXaxis()->SetRangeUser(-3, 85);
    g_sigma_fit->GetYaxis()->SetRangeUser(0.5, 0.6);
    g_sigma_fit->SetMarkerStyle(20);
    g_sigma_fit->SetMarkerSize(1.2);
    g_sigma_fit->SetLineColor(kBlack);
    g_sigma_fit->SetMarkerColor(kBlack);
    
    g_sigma_fit_sys->SetFillColorAlpha(kGray, 0.6); // semi-transparent fill
    g_sigma_fit_sys->SetMarkerStyle(20);
    g_sigma_fit_sys->SetMarkerSize(1.2); // hide the markers for the band
    g_sigma_fit_sys->GetXaxis()->SetNdivisions(505);
    g_sigma_fit_sys->GetYaxis()->SetNdivisions(505);
    g_sigma_fit_sys->GetXaxis()->SetTitle("Centrality [%]");
    g_sigma_fit_sys->GetYaxis()->SetTitle("k_{fit}");
    g_sigma_fit_sys->GetXaxis()->SetRangeUser(-3, 85);
    g_sigma_fit_sys->GetYaxis()->SetRangeUser(0.5, 0.6);

    g_sigma_fit_sys->Draw("AP3"); // draw as a filled band
    g_sigma_fit->Draw("PSAME"); // draw as a filled band

    tx = 0.19;
    ty_start = 0.25;
    for ( auto tag : tags ) {
        tex->DrawLatex(0.19, ty_start, tag.c_str());
        ty_start -= 0.06;
    }
    c->SaveAs(Form("%s/sigma_fits_sys.png", outdir.c_str()));

    TFile * f_out = new TFile(Form("%s/sigma_fits.root", outdir.c_str()), "RECREATE");
    g_sigma_fit->Write("g_sigma_fit");
    g_sigma_fit_sys->Write("g_sigma_fit_sys");
    for ( int i = 0; i < N_X_CENT_BINS_WINDOW; ++i ) {
        g_sigma[i]->Write(Form("g_sigma_cent%d", i));
    }
    f_out->Close();

    delete c;
    delete leg;
    delete tex;

    f->Close();

    return;

}

void DeltaPlots(const std::string input_file_basic, const std::string input_file_random, const std::string input_file_probe, const std::string & input_file_embed, const std::string & prefix)
{
    std::string outdir = MakeGetDir(output_dir+"/"+prefix);
    TFile * f_basic = new TFile(input_file_basic.c_str(), "READ");
    TFile * f_random = new TFile(input_file_random.c_str(), "READ");
    TFile * f_probe = new TFile(input_file_probe.c_str(), "READ");
    TFile * f_embed = new TFile(input_file_embed.c_str(), "READ");
    if ( !f_basic->IsOpen() || f_basic->IsZombie() ) { std::cout << "File " << input_file_basic << " is zombie" << std::endl; exit(1); }
    if ( !f_random->IsOpen() || f_random->IsZombie() ) { std::cout << "File " << input_file_random << " is zombie" << std::endl; exit(1); }
    if ( !f_probe->IsOpen() || f_probe->IsZombie() ) { std::cout << "File " << input_file_probe << " is zombie" << std::endl; exit(1); }
    if ( !f_embed->IsOpen() || f_embed->IsZombie() ) { std::cout << "File " << input_file_embed << " is zombie" << std::endl; exit(1); }

    TH2F * h2_area_res_basic = (TH2F*)f_basic->Get("h2_area_res");
    h2_area_res_basic->SetName("h2_area_res_basic");
    if ( !h2_area_res_basic ) { std::cout << "h2_area_res_basic not found" << std::endl; exit(1); }
    TH2F * h2_area_res_random = (TH2F*)f_random->Get("h2_area_res");
    h2_area_res_random->SetName("h2_area_res_random");
    if ( !h2_area_res_random ) { std::cout << "h2_area_res_random not found" << std::endl; exit(1); }
    TH2F * h2_area_res_probe = (TH2F*)f_probe->Get("h2_area_res");
    h2_area_res_probe->SetName("h2_area_res_probe");
    if ( !h2_area_res_probe ) { std::cout << "h2_area_res_probe not found" << std::endl; exit(1); }
    TH2F * h2_area_res_embed = (TH2F*)f_embed->Get("h2_area_res");
    h2_area_res_embed->SetName("h2_area_res_embed");
    if ( !h2_area_res_embed ) { std::cout << "h2_area_res_embed not found" << std::endl; exit(1); }

    TH2F * h2_mult_res_basic = (TH2F*)f_basic->Get("h2_mult_res");
    h2_mult_res_basic->SetName("h2_mult_res_basic");
    if ( !h2_mult_res_basic ) { std::cout << "h2_mult_res_basic not found" << std::endl; exit(1); }
    TH2F * h2_mult_res_random = (TH2F*)f_random->Get("h2_mult_res");
    h2_mult_res_random->SetName("h2_mult_res_random");
    if ( !h2_mult_res_random ) { std::cout << "h2_mult_res_random not found" << std::endl; exit(1); }
    TH2F * h2_mult_res_probe = (TH2F*)f_probe->Get("h2_mult_res");
    h2_mult_res_probe->SetName("h2_mult_res_probe");
    if ( !h2_mult_res_probe ) { std::cout << "h2_mult_res_probe not found" << std::endl; exit(1); }
    TH2F * h2_mult_res_embed = (TH2F*)f_embed->Get("h2_mult_res");
    h2_mult_res_embed->SetName("h2_mult_res_embed");
    if ( !h2_mult_res_embed ) { std::cout << "h2_mult_res_embed not found" << std::endl; exit(1); }

    TH2F * h2_sub1_res_basic = (TH2F*)f_basic->Get("h2_sub1_res");
    h2_sub1_res_basic->SetName("h2_sub1_res_basic");
    if ( !h2_sub1_res_basic ) { std::cout << "h2_sub1_res_basic not found" << std::endl; exit(1); }
    TH2F * h2_sub1_res_random = (TH2F*)f_random->Get("h2_sub1_res");
    h2_sub1_res_random->SetName("h2_sub1_res_random");
    if ( !h2_sub1_res_random ) { std::cout << "h2_sub1_res_random not found" << std::endl; exit(1); }
    TH2F * h2_sub1_res_probe = (TH2F*)f_probe->Get("h2_sub1_res");
    h2_sub1_res_probe->SetName("h2_sub1_res_probe");
    if ( !h2_sub1_res_probe ) { std::cout << "h2_sub1_res_probe not found" << std::endl; exit(1); }
    TH2F * h2_sub1_res_embed = (TH2F*)f_embed->Get("h2_sub1_res");
    h2_sub1_res_embed->SetName("h2_sub1_res_embed");
    if ( !h2_sub1_res_embed ) { std::cout << "h2_sub1_res_embed not found" << std::endl; exit(1); }

    TH1F * h1_area_basic[N_X_CENT_BINS];
    TH1F * h1_area_random[N_X_CENT_BINS];
    TH1F * h1_area_probe[N_X_CENT_BINS];
    TH1F * h1_area_embed[N_X_CENT_BINS];
    TH1F * h1_mult_basic[N_X_CENT_BINS];
    TH1F * h1_mult_random[N_X_CENT_BINS];
    TH1F * h1_mult_probe[N_X_CENT_BINS];
    TH1F * h1_mult_embed[N_X_CENT_BINS];
    TH1F * h1_sub1_basic[N_X_CENT_BINS];
    TH1F * h1_sub1_random[N_X_CENT_BINS];
    TH1F * h1_sub1_probe[N_X_CENT_BINS];
    TH1F * h1_sub1_embed[N_X_CENT_BINS];

    std::vector<TH2F*> h2s = {h2_area_res_basic, h2_area_res_random, h2_area_res_probe, h2_area_res_embed, 
                            h2_mult_res_basic, h2_mult_res_random, h2_mult_res_probe, h2_mult_res_embed, 
                            h2_sub1_res_basic, h2_sub1_res_random, h2_sub1_res_probe, h2_sub1_res_embed};
    int colors[] = {kRed, kAzure-2, kGreen+2};
    int markers[] = {kFullSquare, kOpenCircle, kFullCross, kFullTriangleUp};
    int markers_rc[] = {kFullSquare, kFullSquare, kFullSquare};
    int markers_random[] = {kOpenCircle, kOpenCircle, kOpenCircle};
    int markers_embed[] = {kFullCross, kFullCross, kFullCross};
    int markers_probe[] = {kFullTriangleUp, kFullTriangleUp, kFullTriangleUp};
    int ihist = 0;
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
        for ( auto h2 : h2s ) {
            h2->GetXaxis()->SetRangeUser(X_CENT_BINS[i], X_CENT_BINS[i+1]);
        }

        h1_area_basic[i] = (TH1F*)h2_area_res_basic->ProjectionY(Form("h1_area_basic_%d", i));
        h1_area_random[i] = (TH1F*)h2_area_res_random->ProjectionY(Form("h1_area_random_%d", i));
        h1_area_probe[i] = (TH1F*)h2_area_res_probe->ProjectionY(Form("h1_area_probe_%d", i));
        h1_area_embed[i] = (TH1F*)h2_area_res_embed->ProjectionY(Form("h1_area_embed_%d", i));
        
        h1_mult_basic[i] = (TH1F*)h2_mult_res_basic->ProjectionY(Form("h1_mult_basic_%d", i));
        h1_mult_random[i] = (TH1F*)h2_mult_res_random->ProjectionY(Form("h1_mult_random_%d", i));
        h1_mult_probe[i] = (TH1F*)h2_mult_res_probe->ProjectionY(Form("h1_mult_probe_%d", i));
        h1_mult_embed[i] = (TH1F*)h2_mult_res_embed->ProjectionY(Form("h1_mult_embed_%d", i));
        
        h1_sub1_basic[i] = (TH1F*)h2_sub1_res_basic->ProjectionY(Form("h1_sub1_basic_%d", i));
        h1_sub1_random[i] = (TH1F*)h2_sub1_res_random->ProjectionY(Form("h1_sub1_random_%d", i));
        h1_sub1_probe[i] = (TH1F*)h2_sub1_res_probe->ProjectionY(Form("h1_sub1_probe_%d", i));
        h1_sub1_embed[i] = (TH1F*)h2_sub1_res_embed->ProjectionY(Form("h1_sub1_embed_%d", i));

        std::vector<TH1F*> h1s = {h1_area_basic[i], h1_area_random[i], h1_area_probe[i], h1_area_embed[i],
                                h1_mult_basic[i], h1_mult_random[i], h1_mult_probe[i], h1_mult_embed[i],
                                h1_sub1_basic[i], h1_sub1_random[i], h1_sub1_probe[i], h1_sub1_embed[i]};
        

        for ( auto h1 : h1s ) {
            h1->GetXaxis()->SetNdivisions(505);
            h1->GetYaxis()->SetNdivisions(505);
            h1->GetXaxis()->SetTitle("#delta E_{T}^{Raw} [GeV]");
            h1->GetYaxis()->SetTitle("Probability Density [A.U.]");
            h1->Scale(1./h1->Integral());
        }
        std::vector<TH1F*> h1s_basic = {h1_area_basic[i], h1_mult_basic[i], h1_sub1_basic[i]};
        std::vector<TH1F*> h1s_random = {h1_area_random[i], h1_mult_random[i], h1_sub1_random[i]};
        std::vector<TH1F*> h1s_probe = {h1_area_probe[i], h1_mult_probe[i], h1_sub1_probe[i]};
        std::vector<TH1F*> h1s_embed = {h1_area_embed[i], h1_mult_embed[i], h1_sub1_embed[i]};
        for ( auto h1 : h1s_basic ) {
            h1->SetMarkerStyle(markers[0]);
            h1->SetMarkerSize(1.5);
        }
        for ( auto h1 : h1s_random ) {
            h1->SetMarkerStyle(markers[1]);
            h1->SetMarkerSize(1.5);
        }
        for ( auto h1 : h1s_probe ) {
            h1->SetMarkerStyle(markers[2]);
            h1->SetMarkerSize(1.5);
        }
        for ( auto h1 : h1s_embed ) {
            h1->SetMarkerStyle(markers[3]);
            h1->SetMarkerSize(1.5);
        }

        std::vector<TH1F*> h1s_area = {h1_area_basic[i], h1_area_random[i], h1_area_probe[i], h1_area_embed[i]};
        std::vector<TH1F*> h1s_mult = {h1_mult_basic[i], h1_mult_random[i], h1_mult_probe[i], h1_mult_embed[i]};
        std::vector<TH1F*> h1s_sub1 = {h1_sub1_basic[i], h1_sub1_random[i], h1_sub1_probe[i], h1_sub1_embed[i]};
        for ( auto h1 : h1s_area ) {
            h1->SetLineColor(colors[0]);
            h1->SetMarkerColor(colors[0]);
        }
        for ( auto h1 : h1s_mult ) {
            h1->SetLineColor(colors[1]);
            h1->SetMarkerColor(colors[1]);
        }
        for ( auto h1 : h1s_sub1 ) {
            h1->SetLineColor(colors[2]);
            h1->SetMarkerColor(colors[2]);
        }


    }

    std::vector<float> mu_lhs_area_basic {};
    std::vector<float> mu_lhs_area_random {};
    std::vector<float> sigma_lhs_area_basic {}; 
    std::vector<float> sigma_lhs_area_random {};
    std::vector<float> mu_area_basic {};
    std::vector<float> mu_area_random {};
    std::vector<float> mu_area_probe {};
    std::vector<float> mu_area_embed {};
    std::vector<float> sigma_area_basic {};
    std::vector<float> sigma_area_random {};
    std::vector<float> sigma_area_probe {};
    std::vector<float> sigma_area_embed {};

    std::vector<float> mu_lhs_area_basic_err {};
    std::vector<float> mu_lhs_area_random_err {};
    std::vector<float> sigma_lhs_area_basic_err {}; 
    std::vector<float> sigma_lhs_area_random_err {};
    std::vector<float> mu_area_basic_err {};
    std::vector<float> mu_area_random_err {};
    std::vector<float> mu_area_probe_err {};
    std::vector<float> mu_area_embed_err {};
    std::vector<float> sigma_area_basic_err {};
    std::vector<float> sigma_area_random_err {};
    std::vector<float> sigma_area_probe_err {};
    std::vector<float> sigma_area_embed_err {};

    std::vector<float> mu_lhs_mult_basic {};
    std::vector<float> mu_lhs_mult_random {};
    std::vector<float> sigma_lhs_mult_basic {}; 
    std::vector<float> sigma_lhs_mult_random {};
    std::vector<float> mu_mult_basic {};
    std::vector<float> mu_mult_random {};
    std::vector<float> mu_mult_probe {};
    std::vector<float> mu_mult_embed {};
    std::vector<float> sigma_mult_basic {};
    std::vector<float> sigma_mult_random {};
    std::vector<float> sigma_mult_probe {};
    std::vector<float> sigma_mult_embed {};

    std::vector<float> mu_lhs_mult_basic_err {};
    std::vector<float> mu_lhs_mult_random_err {};
    std::vector<float> sigma_lhs_mult_basic_err {}; 
    std::vector<float> sigma_lhs_mult_random_err {};
    std::vector<float> mu_mult_basic_err {};
    std::vector<float> mu_mult_random_err {};
    std::vector<float> mu_mult_probe_err {};
    std::vector<float> mu_mult_embed_err {};
    std::vector<float> sigma_mult_basic_err {};
    std::vector<float> sigma_mult_random_err {};
    std::vector<float> sigma_mult_probe_err {};
    std::vector<float> sigma_mult_embed_err {};

    std::vector<float> mu_lhs_sub1_basic {};
    std::vector<float> mu_lhs_sub1_random {};
    std::vector<float> sigma_lhs_sub1_basic {}; 
    std::vector<float> sigma_lhs_sub1_random {};
    std::vector<float> mu_sub1_basic {};
    std::vector<float> mu_sub1_random {};
    std::vector<float> mu_sub1_probe {};
    std::vector<float> mu_sub1_embed {};
    std::vector<float> sigma_sub1_basic {};
    std::vector<float> sigma_sub1_random {};
    std::vector<float> sigma_sub1_probe {};
    std::vector<float> sigma_sub1_embed {};

    std::vector<float> mu_lhs_sub1_basic_err {};
    std::vector<float> mu_lhs_sub1_random_err {};
    std::vector<float> sigma_lhs_sub1_basic_err {}; 
    std::vector<float> sigma_lhs_sub1_random_err {};
    std::vector<float> mu_sub1_basic_err {};
    std::vector<float> mu_sub1_random_err {};
    std::vector<float> mu_sub1_probe_err {};
    std::vector<float> mu_sub1_embed_err {};
    std::vector<float> sigma_sub1_basic_err {};
    std::vector<float> sigma_sub1_random_err {};
    std::vector<float> sigma_sub1_probe_err {};
    std::vector<float> sigma_sub1_embed_err {};

    TCanvas * c;
    TLegend * leg;
    bool do_gamma = DO_GAMMA;
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {

        TH1F * h1_a_basic = (TH1F*)h1_area_basic[i]->Clone(Form("h1_a_basic_%d", i));
        TH1F * h1_a_random = (TH1F*)h1_area_random[i]->Clone(Form("h1_a_random_%d", i));
        TH1F * h1_a_probe = (TH1F*)h1_area_probe[i]->Clone(Form("h1_a_probe_%d", i));
        TH1F * h1_a_embed = (TH1F*)h1_area_embed[i]->Clone(Form("h1_a_embed_%d", i));

        float mean_a_basic = h1_a_basic->GetMean();
        float mean_a_random = h1_a_random->GetMean();
        float mean_a_probe = h1_a_probe->GetMean();
        float mean_a_embed = h1_a_embed->GetMean();
        float rms_a_basic = h1_a_basic->GetRMS();
        float rms_a_random = h1_a_random->GetRMS();
        float rms_a_probe = h1_a_probe->GetRMS();
        float rms_a_embed = h1_a_embed->GetRMS();
        float mean_a_basic_err = h1_a_basic->GetMeanError();
        float mean_a_random_err = h1_a_random->GetMeanError();
        float mean_a_probe_err = h1_a_probe->GetMeanError();
        float mean_a_embed_err = h1_a_embed->GetMeanError();
        float rms_a_basic_err = h1_a_basic->GetRMSError();
        float rms_a_random_err = h1_a_random->GetRMSError();
        float rms_a_probe_err = h1_a_probe->GetRMSError();
        float rms_a_embed_err = h1_a_embed->GetRMSError();

       
        
        

        TF1 * fit_gaus = new TF1("fit_gaus", "gaus", -40, 40);
        h1_a_basic->Fit(fit_gaus, "RQ", "", -40, 0);
        h1_a_basic->Fit(fit_gaus, "RQ", "", fit_gaus->GetParameter(1)-3*fit_gaus->GetParameter(2), fit_gaus->GetParameter(1)+0.5*fit_gaus->GetParameter(2));
        float mu_lhs_a_basic = fit_gaus->GetParameter(1);
        float sigma_lhs_a_basic = fit_gaus->GetParameter(2);
        float mu_lhs_a_basic_err = fit_gaus->GetParError(1);
        float sigma_lhs_a_basic_err = fit_gaus->GetParError(2);

        TF1 * fit_gaus2 = new TF1("fit_gaus2", "gaus", -40, 40);
        h1_a_random->Fit(fit_gaus2, "RQ", "", -40, 0);
        h1_a_random->Fit(fit_gaus2, "RQ", "", fit_gaus2->GetParameter(1)-3*fit_gaus2->GetParameter(2), fit_gaus2->GetParameter(1)+0.5*fit_gaus2->GetParameter(2));
        float mu_lhs_a_random = fit_gaus2->GetParameter(1);
        float sigma_lhs_a_random = fit_gaus2->GetParameter(2);
        float mu_lhs_a_random_err = fit_gaus2->GetParError(1);
        float sigma_lhs_a_random_err = fit_gaus2->GetParError(2);


        // fill vectors
        mu_lhs_area_basic.push_back(mu_lhs_a_basic);
        mu_lhs_area_random.push_back(mu_lhs_a_random);
        sigma_lhs_area_basic.push_back(sigma_lhs_a_basic);
        sigma_lhs_area_random.push_back(sigma_lhs_a_random);
        mu_lhs_area_basic_err.push_back(mu_lhs_a_basic_err);
        mu_lhs_area_random_err.push_back(mu_lhs_a_random_err);
        sigma_lhs_area_basic_err.push_back(sigma_lhs_a_basic_err);
        sigma_lhs_area_random_err.push_back(sigma_lhs_a_random_err);

        mu_area_basic.push_back(mean_a_basic);
        mu_area_random.push_back(mean_a_random);
        mu_area_probe.push_back(mean_a_probe);
        mu_area_embed.push_back(mean_a_embed);
        sigma_area_basic.push_back(rms_a_basic);
        sigma_area_random.push_back(rms_a_random);
        sigma_area_probe.push_back(rms_a_probe);
        sigma_area_embed.push_back(rms_a_embed);
        mu_area_basic_err.push_back(mean_a_basic_err);
        mu_area_random_err.push_back(mean_a_random_err);
        mu_area_probe_err.push_back(mean_a_probe_err);
        mu_area_embed_err.push_back(mean_a_embed_err);
        sigma_area_basic_err.push_back(rms_a_basic_err);
        sigma_area_random_err.push_back(rms_a_random_err);
        sigma_area_probe_err.push_back(rms_a_probe_err);
        sigma_area_embed_err.push_back(rms_a_embed_err);

   
        // float ap0 =( mean_a_basic*mean_a_basic);
        // float ab0 = ap0/mean_a_basic;

        TF1 * fit_gamma = new TF1("fit_gamma", myGammaFunction, -30, 50, 3);
        fit_gamma->SetParNames("A","a_b","a_p");
        fit_gamma->SetParameters(0.1, 1.18, 120);
        // fit_gamma->SetParLimits(0, 0.2, 1.2);
        // fit_gamma->SetParLimits(1, 1, 2.5);
        // fit_gamma->SetParLimits(2, 50, 130);
        h1_a_random->Fit("fit_gamma", "LM+QR", "", 5, 50);
        float a_a_basic = fit_gamma->GetParameter(0);
        float ab_a_basic = fit_gamma->GetParameter(1);
        float ap_a_basic = fit_gamma->GetParameter(2);
        float m = ap_a_basic/ab_a_basic;
        float s = sqrt(ap_a_basic)/ab_a_basic;
    
        TF1 * fit_gamma2 = new TF1("fit_gamma2", myGammaFunction, -30, 30, 3);
        fit_gamma2->SetParameter(0, a_a_basic);
        fit_gamma2->SetParameter(1, ab_a_basic);
        fit_gamma2->SetParameter(2, ap_a_basic);
        // std::cout << "m = " << m << ", s = " << s << std::endl;
        // h1_a_random->Fit(fit_gamma2, "LQR", "", m-3*s, m+3*s);
        a_a_basic = fit_gamma2->GetParameter(0);
        ab_a_basic = fit_gamma2->GetParameter(1);
        ap_a_basic = fit_gamma2->GetParameter(2);
        float chi2overndf = fit_gamma2->GetChisquare()/fit_gamma2->GetNDF();
        //  std::cout << "m = " << m << ", s = " << s << std::endl;

        c = new TCanvas("c", "c", 800, 800);

        if(i == 0) {
            // c->SetLogy();
            TFile * f = new TFile(Form("%s/fit_area_cent0_results.root", outdir.c_str()), "RECREATE");
            h1_a_basic->Write();
            h1_a_random->Write();
            h1_a_probe->Write();
            h1_a_embed->Write();
            fit_gaus->Write();
            fit_gaus2->Write();
            fit_gamma->Write();
            f->Close();
            delete f;
        }

        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogy();
        // leg = new TLegend(0.15,0.8,0.4,0.92);
        leg = new TLegend(0.5,0.7,0.8,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(1);
        leg->SetTextSize(0.04);

        TLegend * leg2 = new TLegend(0.17,0.66,0.42,0.8);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetNColumns(1);
        leg2->SetTextSize(0.04);

        h1_area_basic[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_area_basic[i]->GetXaxis()->SetNdivisions(510);
        h1_area_basic[i]->GetYaxis()->SetRangeUser(1e-4, 2e0);
        h1_area_basic[i]->SetMarkerSize(1.5);
        h1_area_basic[i]->SetMarkerColor(kAzure-2);
        h1_area_basic[i]->SetLineColor(kAzure-2);
        h1_area_basic[i]->GetXaxis()->SetTitle("#delta E_{T}^{Raw} [GeV]");
        h1_area_basic[i]->Draw("P");
        // leg->AddEntry(h1_area_basic[i], Form("Basic: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_basic, rms_a_basic), "pe");
        leg->AddEntry(h1_area_basic[i], "Basic Cones", "pe");

        fit_gaus->SetLineColor(kAzure-2);
        fit_gaus->SetLineStyle(2);
        fit_gaus->SetLineWidth(2);
        fit_gaus->SetParameter(1, mu_lhs_a_basic);
        fit_gaus->SetParameter(2, sigma_lhs_a_basic);
        // fit_gaus->Draw("SAME");
      
        h1_area_random[i]->SetMarkerSize(1.5);
        h1_area_random[i]->SetMarkerColor(kRed);
        h1_area_random[i]->SetLineColor(kRed);
        h1_area_random[i]->Draw("P SAME");
        // leg->AddEntry(h1_area_random[i], Form("Randomized#kern[0.05]{#eta#phi}: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_random, rms_a_random), "pe");
        leg->AddEntry(h1_area_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        leg->AddEntry(fit_gaus2, "Gaussian Fit", "l");

        fit_gaus2->SetParameter(1, mu_lhs_a_random);
        fit_gaus2->SetParameter(2, sigma_lhs_a_random);
        fit_gaus2->SetLineColor(kRed);
        fit_gaus2->SetLineStyle(3);
        fit_gaus2->SetLineWidth(2);
        fit_gaus2->Draw("SAME");
        
        h1_a_probe->SetMarkerSize(1.5);
        h1_a_probe->SetMarkerColor(kGreen+2);
        h1_a_probe->SetLineColor(kGreen+2);
        h1_a_probe->Draw("P SAME");
        // leg->AddEntry(h1_a_probe, Form("High E_{T} Probe: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_probe, rms_a_probe), "pe");
        leg->AddEntry(h1_a_probe, "High E_{T} Probe", "pe");

        h1_a_embed->SetMarkerSize(1.5);
        h1_a_embed->SetMarkerColor(kBlack);
        h1_a_embed->SetLineColor(kBlack);
        h1_a_embed->Draw("SAME");
        // leg->AddEntry(h1_a_embed, Form("Embed: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_embed, rms_a_embed), "pe");
        leg->AddEntry(h1_a_embed, "Embed", "pe");

        
        // leg->AddEntry(fit_gaus2, Form("#mu^{lhs} = %0.2f GeV, #sigma^{lhs} = %0.2f GeV", mu_lhs_a_random, sigma_lhs_a_random), "l");

        if(do_gamma) {
            fit_gamma2->SetLineColor(kBlack);
            fit_gamma2->SetParameter(0, a_a_basic);
            fit_gamma2->SetParameter(1, ab_a_basic);
            fit_gamma2->SetParameter(2, ap_a_basic);
            fit_gamma2->SetLineStyle(1);
            fit_gamma2->SetLineWidth(2);
            fit_gamma2->Draw("SAME");
            leg->AddEntry(fit_gamma2, Form("f_{#Gamma}: b = %0.2f GeV, p = %0.1f", ab_a_basic, ap_a_basic), "l");
        }
   
        std::vector<std::string> tags = {sPHENIX_Tag,"Au+Au#kern[0.05]{#sqrt{s_{NN}} = 200 GeV}" ,Form("%d-%d%% Central", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])), "#it{Area Method}"};
        TLatex * tex = new TLatex();
        tex->SetNDC();
        tex->SetTextAlign(12); // center alignment
        tex->SetTextSize(0.04);
        float tx = 0.19;
        float ty_start = 0.87;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty_start, tag.c_str());
            ty_start -= 0.05;
        }
        leg->Draw("SAME");
        // leg2->Draw("SAME");
        c->SaveAs(Form("%s/area_basic_%d.png", outdir.c_str(), i));

        delete c;
        delete leg;
        delete fit_gaus;
        delete fit_gamma;
    }

    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {

        TH1F * h1_a_basic = (TH1F*)h1_mult_basic[i]->Clone(Form("h1_a_basic_%d", i));
        TH1F * h1_a_random = (TH1F*)h1_mult_random[i]->Clone(Form("h1_a_random_%d", i));
        TH1F * h1_a_probe = (TH1F*)h1_mult_probe[i]->Clone(Form("h1_a_probe_%d", i));
        TH1F * h1_a_embed = (TH1F*)h1_mult_embed[i]->Clone(Form("h1_a_embed_%d", i));

        float mean_a_basic = h1_a_basic->GetMean();
        float mean_a_random = h1_a_random->GetMean();
        float mean_a_probe = h1_a_probe->GetMean();
        float mean_a_embed = h1_a_embed->GetMean();
        float rms_a_basic = h1_a_basic->GetRMS();
        float rms_a_random = h1_a_random->GetRMS();
        float rms_a_probe = h1_a_probe->GetRMS();
        float rms_a_embed = h1_a_embed->GetRMS();
        float mean_a_basic_err = h1_a_basic->GetMeanError();
        float mean_a_random_err = h1_a_random->GetMeanError();
        float mean_a_probe_err = h1_a_probe->GetMeanError();
        float mean_a_embed_err = h1_a_embed->GetMeanError();
        float rms_a_basic_err = h1_a_basic->GetRMSError();
        float rms_a_random_err = h1_a_random->GetRMSError();
        float rms_a_probe_err = h1_a_probe->GetRMSError();
        float rms_a_embed_err = h1_a_embed->GetRMSError();

       
        
        

        TF1 * fit_gaus = new TF1("fit_gaus", "gaus", -40, 40);
        h1_a_basic->Fit(fit_gaus, "RQ", "", -40, 0);
        h1_a_basic->Fit(fit_gaus, "RQ", "", fit_gaus->GetParameter(1)-3*fit_gaus->GetParameter(2), fit_gaus->GetParameter(1)+0.5*fit_gaus->GetParameter(2));
        float mu_lhs_a_basic = fit_gaus->GetParameter(1);
        float sigma_lhs_a_basic = fit_gaus->GetParameter(2);
        float mu_lhs_a_basic_err = fit_gaus->GetParError(1);
        float sigma_lhs_a_basic_err = fit_gaus->GetParError(2);

        TF1 * fit_gaus2 = new TF1("fit_gaus2", "gaus", -40, 40);
        h1_a_random->Fit(fit_gaus2, "RQ", "", -40, 0);
        h1_a_random->Fit(fit_gaus2, "RQ", "", fit_gaus2->GetParameter(1)-3*fit_gaus2->GetParameter(2), fit_gaus2->GetParameter(1)+0.5*fit_gaus2->GetParameter(2));
        float mu_lhs_a_random = fit_gaus2->GetParameter(1);
        float sigma_lhs_a_random = fit_gaus2->GetParameter(2);
        float mu_lhs_a_random_err = fit_gaus2->GetParError(1);
        float sigma_lhs_a_random_err = fit_gaus2->GetParError(2);


        // fill vectors
        mu_lhs_mult_basic.push_back(mu_lhs_a_basic);
        mu_lhs_mult_random.push_back(mu_lhs_a_random);
        sigma_lhs_mult_basic.push_back(sigma_lhs_a_basic);
        sigma_lhs_mult_random.push_back(sigma_lhs_a_random);
        mu_lhs_mult_basic_err.push_back(mu_lhs_a_basic_err);
        mu_lhs_mult_random_err.push_back(mu_lhs_a_random_err);
        sigma_lhs_mult_basic_err.push_back(sigma_lhs_a_basic_err);
        sigma_lhs_mult_random_err.push_back(sigma_lhs_a_random_err);

        mu_mult_basic.push_back(mean_a_basic);
        mu_mult_random.push_back(mean_a_random);
        mu_mult_probe.push_back(mean_a_probe);
        mu_mult_embed.push_back(mean_a_embed);
        sigma_mult_basic.push_back(rms_a_basic);
        sigma_mult_random.push_back(rms_a_random);
        sigma_mult_probe.push_back(rms_a_probe);
        sigma_mult_embed.push_back(rms_a_embed);
        mu_mult_basic_err.push_back(mean_a_basic_err);
        mu_mult_random_err.push_back(mean_a_random_err);
        mu_mult_probe_err.push_back(mean_a_probe_err);
        mu_mult_embed_err.push_back(mean_a_embed_err);
        sigma_mult_basic_err.push_back(rms_a_basic_err);
        sigma_mult_random_err.push_back(rms_a_random_err);
        sigma_mult_probe_err.push_back(rms_a_probe_err);
        sigma_mult_embed_err.push_back(rms_a_embed_err);

   
        // float ap0 =( mean_a_basic*mean_a_basic);
        // float ab0 = ap0/mean_a_basic;

        TF1 * fit_gamma = new TF1("fit_gamma", myGammaFunction, -30, 50, 3);
        fit_gamma->SetParNames("A","a_b","a_p");
        fit_gamma->SetParameters(0.1, 1.18, 120);
        // fit_gamma->SetParLimits(0, 0.2, 1.2);
        // fit_gamma->SetParLimits(1, 1, 2.5);
        // fit_gamma->SetParLimits(2, 50, 130);
        h1_a_random->Fit("fit_gamma", "LM+QR", "", 5, 50);
        float a_a_basic = fit_gamma->GetParameter(0);
        float ab_a_basic = fit_gamma->GetParameter(1);
        float ap_a_basic = fit_gamma->GetParameter(2);
        float m = ap_a_basic/ab_a_basic;
        float s = sqrt(ap_a_basic)/ab_a_basic;
    
        TF1 * fit_gamma2 = new TF1("fit_gamma2", myGammaFunction, -30, 30, 3);
        fit_gamma2->SetParameter(0, a_a_basic);
        fit_gamma2->SetParameter(1, ab_a_basic);
        fit_gamma2->SetParameter(2, ap_a_basic);
        // std::cout << "m = " << m << ", s = " << s << std::endl;
        // h1_a_random->Fit(fit_gamma2, "LQR", "", m-3*s, m+3*s);
        a_a_basic = fit_gamma2->GetParameter(0);
        ab_a_basic = fit_gamma2->GetParameter(1);
        ap_a_basic = fit_gamma2->GetParameter(2);
        float chi2overndf = fit_gamma2->GetChisquare()/fit_gamma2->GetNDF();
        //  std::cout << "m = " << m << ", s = " << s << std::endl;

        c = new TCanvas("c", "c", 800, 800);

        if(i == 0) {
            // c->SetLogy();
            TFile * f = new TFile(Form("%s/fit_mult_cent0_results.root", outdir.c_str()), "RECREATE");
            h1_a_basic->Write();
            h1_a_random->Write();
            h1_a_probe->Write();
            h1_a_embed->Write();
            fit_gaus->Write();
            fit_gaus2->Write();
            fit_gamma->Write();
            f->Close();
            delete f;
        }

        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogy();
        // leg = new TLegend(0.15,0.8,0.4,0.92);
        leg = new TLegend(0.6,0.7,0.85,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(1);
        leg->SetTextSize(0.04);

        TLegend * leg2 = new TLegend(0.17,0.66,0.42,0.8);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetNColumns(1);
        leg2->SetTextSize(0.04);

        h1_mult_basic[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_mult_basic[i]->GetXaxis()->SetNdivisions(510);
        h1_mult_basic[i]->GetYaxis()->SetRangeUser(1e-4, 2e0);
        h1_mult_basic[i]->SetMarkerSize(1.5);
        h1_mult_basic[i]->SetMarkerColor(kAzure-2);
        h1_mult_basic[i]->SetLineColor(kAzure-2);
        h1_mult_basic[i]->GetXaxis()->SetTitle("#delta E_{T}^{Raw} [GeV]");
        h1_mult_basic[i]->Draw("P");
        // leg->AddEntry(h1_mult_basic[i], Form("Basic: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_basic, rms_a_basic), "pe");
        leg->AddEntry(h1_mult_basic[i], "Basic Cones", "pe");

        fit_gaus->SetLineColor(kAzure-2);
        fit_gaus->SetLineStyle(2);
        fit_gaus->SetLineWidth(2);
        fit_gaus->SetParameter(1, mu_lhs_a_basic);
        fit_gaus->SetParameter(2, sigma_lhs_a_basic);
        // fit_gaus->Draw("SAME");
      
        h1_mult_random[i]->SetMarkerSize(1.5);
        h1_mult_random[i]->SetMarkerColor(kRed);
        h1_mult_random[i]->SetLineColor(kRed);
        h1_mult_random[i]->Draw(" P SAME");
        // leg->AddEntry(h1_mult_random[i], Form("Randomized#kern[0.05]{#eta#phi}: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_random, rms_a_random), "pe");
        leg->AddEntry(h1_mult_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        leg->AddEntry(fit_gaus2, "Gaussian Fit", "l");

        fit_gaus2->SetParameter(1, mu_lhs_a_random);
        fit_gaus2->SetParameter(2, sigma_lhs_a_random);
        fit_gaus2->SetLineColor(kRed);
        fit_gaus2->SetLineStyle(3);
        fit_gaus2->SetLineWidth(2);
        fit_gaus2->Draw("SAME");
        
        h1_a_probe->SetMarkerSize(1.5);
        h1_a_probe->SetMarkerColor(kGreen+2);
        h1_a_probe->SetLineColor(kGreen+2);
        h1_a_probe->Draw("P SAME");
        // leg->AddEntry(h1_a_probe, Form("High E_{T} Probe: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_probe, rms_a_probe), "pe");
        leg->AddEntry(h1_a_probe, "High E_{T} Probe", "pe");

        h1_a_embed->SetMarkerSize(1.5);
        h1_a_embed->SetMarkerColor(kBlack);
        h1_a_embed->SetLineColor(kBlack);
        h1_a_embed->Draw("SAME");
        // leg->AddEntry(h1_a_embed, Form("Embed: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_embed, rms_a_embed), "pe");
        leg->AddEntry(h1_a_embed, "Embed", "pe");

        
        // leg->AddEntry(fit_gaus2, Form("#mu^{lhs} = %0.2f GeV, #sigma^{lhs} = %0.2f GeV", mu_lhs_a_random, sigma_lhs_a_random), "l");

        if(do_gamma) {
            fit_gamma2->SetLineColor(kBlack);
            fit_gamma2->SetParameter(0, a_a_basic);
            fit_gamma2->SetParameter(1, ab_a_basic);
            fit_gamma2->SetParameter(2, ap_a_basic);
            fit_gamma2->SetLineStyle(1);
            fit_gamma2->SetLineWidth(2);
            fit_gamma2->Draw("SAME");
            leg->AddEntry(fit_gamma2, Form("f_{#Gamma}: b = %0.2f GeV, p = %0.1f", ab_a_basic, ap_a_basic), "l");
        }
   
        std::vector<std::string> tags = {sPHENIX_Tag,"Au+Au#kern[0.05]{#sqrt{s_{NN}} = 200 GeV}" ,Form("%d-%d%% Central", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])), "#it{Multiplicity Method}"};
        TLatex * tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);
        float tx = 0.19;
        float ty_start = 0.87;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty_start, tag.c_str());
            ty_start -= 0.05;
        }
        leg->Draw("SAME");
        // leg2->Draw("SAME");
        c->SaveAs(Form("%s/mult_basic_%d.png", outdir.c_str(), i));

        delete c;
        delete leg;
        delete fit_gaus;
        delete fit_gamma;
    }

    
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {

        TH1F * h1_a_basic = (TH1F*)h1_sub1_basic[i]->Clone(Form("h1_a_basic_%d", i));
        TH1F * h1_a_random = (TH1F*)h1_sub1_random[i]->Clone(Form("h1_a_random_%d", i));
        TH1F * h1_a_probe = (TH1F*)h1_sub1_probe[i]->Clone(Form("h1_a_probe_%d", i));
        TH1F * h1_a_embed = (TH1F*)h1_sub1_embed[i]->Clone(Form("h1_a_embed_%d", i));

        float mean_a_basic = h1_a_basic->GetMean();
        float mean_a_random = h1_a_random->GetMean();
        float mean_a_probe = h1_a_probe->GetMean();
        float mean_a_embed = h1_a_embed->GetMean();
        float rms_a_basic = h1_a_basic->GetRMS();
        float rms_a_random = h1_a_random->GetRMS();
        float rms_a_probe = h1_a_probe->GetRMS();
        float rms_a_embed = h1_a_embed->GetRMS();
        float mean_a_basic_err = h1_a_basic->GetMeanError();
        float mean_a_random_err = h1_a_random->GetMeanError();
        float mean_a_probe_err = h1_a_probe->GetMeanError();
        float mean_a_embed_err = h1_a_embed->GetMeanError();
        float rms_a_basic_err = h1_a_basic->GetRMSError();
        float rms_a_random_err = h1_a_random->GetRMSError();
        float rms_a_probe_err = h1_a_probe->GetRMSError();
        float rms_a_embed_err = h1_a_embed->GetRMSError();

       
        
        

        TF1 * fit_gaus = new TF1("fit_gaus", "gaus", -40, 40);
        h1_a_basic->Fit(fit_gaus, "RQ", "", -40, 0);
        h1_a_basic->Fit(fit_gaus, "RQ", "", fit_gaus->GetParameter(1)-3*fit_gaus->GetParameter(2), fit_gaus->GetParameter(1)+0.5*fit_gaus->GetParameter(2));
        float mu_lhs_a_basic = fit_gaus->GetParameter(1);
        float sigma_lhs_a_basic = fit_gaus->GetParameter(2);
        float mu_lhs_a_basic_err = fit_gaus->GetParError(1);
        float sigma_lhs_a_basic_err = fit_gaus->GetParError(2);

        TF1 * fit_gaus2 = new TF1("fit_gaus2", "gaus", -40, 40);
        h1_a_random->Fit(fit_gaus2, "RQ", "", -40, 0);
        h1_a_random->Fit(fit_gaus2, "RQ", "", fit_gaus2->GetParameter(1)-3*fit_gaus2->GetParameter(2), fit_gaus2->GetParameter(1)+0.5*fit_gaus2->GetParameter(2));
        float mu_lhs_a_random = fit_gaus2->GetParameter(1);
        float sigma_lhs_a_random = fit_gaus2->GetParameter(2);
        float mu_lhs_a_random_err = fit_gaus2->GetParError(1);
        float sigma_lhs_a_random_err = fit_gaus2->GetParError(2);


        // fill vectors
        mu_lhs_sub1_basic.push_back(mu_lhs_a_basic);
        mu_lhs_sub1_random.push_back(mu_lhs_a_random);
        sigma_lhs_sub1_basic.push_back(sigma_lhs_a_basic);
        sigma_lhs_sub1_random.push_back(sigma_lhs_a_random);
        mu_lhs_sub1_basic_err.push_back(mu_lhs_a_basic_err);
        mu_lhs_sub1_random_err.push_back(mu_lhs_a_random_err);
        sigma_lhs_sub1_basic_err.push_back(sigma_lhs_a_basic_err);
        sigma_lhs_sub1_random_err.push_back(sigma_lhs_a_random_err);

        mu_sub1_basic.push_back(mean_a_basic);
        mu_sub1_random.push_back(mean_a_random);
        mu_sub1_probe.push_back(mean_a_probe);
        mu_sub1_embed.push_back(mean_a_embed);
        sigma_sub1_basic.push_back(rms_a_basic);
        sigma_sub1_random.push_back(rms_a_random);
        sigma_sub1_probe.push_back(rms_a_probe);
        sigma_sub1_embed.push_back(rms_a_embed);
        mu_sub1_basic_err.push_back(mean_a_basic_err);
        mu_sub1_random_err.push_back(mean_a_random_err);
        mu_sub1_probe_err.push_back(mean_a_probe_err);
        mu_sub1_embed_err.push_back(mean_a_embed_err);
        sigma_sub1_basic_err.push_back(rms_a_basic_err);
        sigma_sub1_random_err.push_back(rms_a_random_err);
        sigma_sub1_probe_err.push_back(rms_a_probe_err);
        sigma_sub1_embed_err.push_back(rms_a_embed_err);

   
        // float ap0 =( mean_a_basic*mean_a_basic);
        // float ab0 = ap0/mean_a_basic;

        TF1 * fit_gamma = new TF1("fit_gamma", myGammaFunction, -30, 50, 3);
        fit_gamma->SetParNames("A","a_b","a_p");
        fit_gamma->SetParameters(0.1, 1.18, 120);
        // fit_gamma->SetParLimits(0, 0.2, 1.2);
        // fit_gamma->SetParLimits(1, 1, 2.5);
        // fit_gamma->SetParLimits(2, 50, 130);
        h1_a_random->Fit("fit_gamma", "LM+QR", "", 5, 50);
        float a_a_basic = fit_gamma->GetParameter(0);
        float ab_a_basic = fit_gamma->GetParameter(1);
        float ap_a_basic = fit_gamma->GetParameter(2);
        float m = ap_a_basic/ab_a_basic;
        float s = sqrt(ap_a_basic)/ab_a_basic;
    
        TF1 * fit_gamma2 = new TF1("fit_gamma2", myGammaFunction, -30, 30, 3);
        fit_gamma2->SetParameter(0, a_a_basic);
        fit_gamma2->SetParameter(1, ab_a_basic);
        fit_gamma2->SetParameter(2, ap_a_basic);
        // std::cout << "m = " << m << ", s = " << s << std::endl;
        // h1_a_random->Fit(fit_gamma2, "LQR", "", m-3*s, m+3*s);
        a_a_basic = fit_gamma2->GetParameter(0);
        ab_a_basic = fit_gamma2->GetParameter(1);
        ap_a_basic = fit_gamma2->GetParameter(2);
        float chi2overndf = fit_gamma2->GetChisquare()/fit_gamma2->GetNDF();
        //  std::cout << "m = " << m << ", s = " << s << std::endl;

        c = new TCanvas("c", "c", 800, 800);

        if(i == 0) {
            // c->SetLogy();
            TFile * f = new TFile(Form("%s/fit_sub1_cent0_results.root", outdir.c_str()), "RECREATE");
            h1_a_basic->Write();
            h1_a_random->Write();
            h1_a_probe->Write();
            h1_a_embed->Write();
            fit_gaus->Write();
            fit_gaus2->Write();
            fit_gamma->Write();
            f->Close();
            delete f;
        }

        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogy();
        // leg = new TLegend(0.15,0.8,0.4,0.92);
        leg = new TLegend(0.6,0.7,0.85,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(1);
        leg->SetTextSize(0.04);

        TLegend * leg2 = new TLegend(0.17,0.66,0.42,0.8);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetNColumns(1);
        leg2->SetTextSize(0.04);

        h1_sub1_basic[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_sub1_basic[i]->GetXaxis()->SetNdivisions(510);
        h1_sub1_basic[i]->GetYaxis()->SetRangeUser(1e-4, 2e0);
        h1_sub1_basic[i]->SetMarkerSize(1.5);
        h1_sub1_basic[i]->SetMarkerColor(kAzure-2);
        h1_sub1_basic[i]->SetLineColor(kAzure-2);
        h1_sub1_basic[i]->GetXaxis()->SetTitle("#delta E_{T}^{Raw} [GeV]");
        h1_sub1_basic[i]->Draw("P");
        // leg->AddEntry(h1_sub1_basic[i], Form("Basic: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_basic, rms_a_basic), "pe");
        leg->AddEntry(h1_sub1_basic[i], "Basic Cones", "pe");

        fit_gaus->SetLineColor(kAzure-2);
        fit_gaus->SetLineStyle(2);
        fit_gaus->SetLineWidth(2);
        fit_gaus->SetParameter(1, mu_lhs_a_basic);
        fit_gaus->SetParameter(2, sigma_lhs_a_basic);
        // fit_gaus->Draw("SAME");
      
        h1_sub1_random[i]->SetMarkerSize(1.5);
        h1_sub1_random[i]->SetMarkerColor(kRed);
        h1_sub1_random[i]->SetLineColor(kRed);
        h1_sub1_random[i]->Draw(" P SAME");
        // leg->AddEntry(h1_sub1_random[i], Form("Randomized#kern[0.05]{#eta#phi}: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_random, rms_a_random), "pe");
        leg->AddEntry(h1_sub1_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        leg->AddEntry(fit_gaus2, "Gaussian Fit", "l");

        fit_gaus2->SetParameter(1, mu_lhs_a_random);
        fit_gaus2->SetParameter(2, sigma_lhs_a_random);
        fit_gaus2->SetLineColor(kRed);
        fit_gaus2->SetLineStyle(3);
        fit_gaus2->SetLineWidth(2);
        fit_gaus2->Draw("SAME");
        
        h1_a_probe->SetMarkerSize(1.5);
        h1_a_probe->SetMarkerColor(kGreen+2);
        h1_a_probe->SetLineColor(kGreen+2);
        h1_a_probe->Draw("P SAME");
        // leg->AddEntry(h1_a_probe, Form("High E_{T} Probe: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_probe, rms_a_probe), "pe");
        leg->AddEntry(h1_a_probe, "High E_{T} Probe", "pe");

        h1_a_embed->SetMarkerSize(1.5);
        h1_a_embed->SetMarkerColor(kBlack);
        h1_a_embed->SetLineColor(kBlack);
        h1_a_embed->Draw("SAME");
        // leg->AddEntry(h1_a_embed, Form("Embed: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_embed, rms_a_embed), "pe");
        leg->AddEntry(h1_a_embed, "Embed", "pe");

        
        // leg->AddEntry(fit_gaus2, Form("#mu^{lhs} = %0.2f GeV, #sigma^{lhs} = %0.2f GeV", mu_lhs_a_random, sigma_lhs_a_random), "l");

        if(do_gamma) {
            fit_gamma2->SetLineColor(kBlack);
            fit_gamma2->SetParameter(0, a_a_basic);
            fit_gamma2->SetParameter(1, ab_a_basic);
            fit_gamma2->SetParameter(2, ap_a_basic);
            fit_gamma2->SetLineStyle(1);
            fit_gamma2->SetLineWidth(2);
            fit_gamma2->Draw("SAME");
            leg->AddEntry(fit_gamma2, Form("f_{#Gamma}: b = %0.2f GeV, p = %0.1f", ab_a_basic, ap_a_basic), "l");
        }
   
        std::vector<std::string> tags = {sPHENIX_Tag,"Au+Au#kern[0.05]{#sqrt{s_{NN}} = 200 GeV}" ,Form("%d-%d%% Central", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])), "#it{Iterative Method}"};
        TLatex * tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);
        float tx = 0.19;
        float ty_start = 0.87;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty_start, tag.c_str());
            ty_start -= 0.05;
        }
        leg->Draw("SAME");
        // leg2->Draw("SAME");
        c->SaveAs(Form("%s/sub1_basic_%d.png", outdir.c_str(), i));

        delete c;
        delete leg;
        delete fit_gaus;
        delete fit_gamma;
    }


    for ( int i = 0; i < 1; ++i ) {

        TH1F * h1_a_basic = (TH1F*)h1_sub1_basic[i]->Clone(Form("h1_a_basic_%d", i));
        TH1F * h1_a_random = (TH1F*)h1_sub1_random[i]->Clone(Form("h1_a_random_%d", i));
        TH1F * h1_a_probe = (TH1F*)h1_sub1_probe[i]->Clone(Form("h1_a_probe_%d", i));
        TH1F * h1_a_embed = (TH1F*)h1_sub1_embed[i]->Clone(Form("h1_a_embed_%d", i));

        float mean_a_basic = h1_a_basic->GetMean();
        float mean_a_random = h1_a_random->GetMean();
        float mean_a_probe = h1_a_probe->GetMean();
        float mean_a_embed = h1_a_embed->GetMean();
        float rms_a_basic = h1_a_basic->GetRMS();
        float rms_a_random = h1_a_random->GetRMS();
        float rms_a_probe = h1_a_probe->GetRMS();
        float rms_a_embed = h1_a_embed->GetRMS();
        float mean_a_basic_err = h1_a_basic->GetMeanError();
        float mean_a_random_err = h1_a_random->GetMeanError();
        float mean_a_probe_err = h1_a_probe->GetMeanError();
        float mean_a_embed_err = h1_a_embed->GetMeanError();
        float rms_a_basic_err = h1_a_basic->GetRMSError();
        float rms_a_random_err = h1_a_random->GetRMSError();
        float rms_a_probe_err = h1_a_probe->GetRMSError();
        float rms_a_embed_err = h1_a_embed->GetRMSError();

       
        
        

        TF1 * fit_gaus = new TF1("fit_gaus", "gaus", -40, 40);
        h1_a_basic->Fit(fit_gaus, "RQ", "", -40, 0);
        h1_a_basic->Fit(fit_gaus, "RQ", "", fit_gaus->GetParameter(1)-3*fit_gaus->GetParameter(2), fit_gaus->GetParameter(1)+0.5*fit_gaus->GetParameter(2));
        float mu_lhs_a_basic = fit_gaus->GetParameter(1);
        float sigma_lhs_a_basic = fit_gaus->GetParameter(2);
        float mu_lhs_a_basic_err = fit_gaus->GetParError(1);
        float sigma_lhs_a_basic_err = fit_gaus->GetParError(2);

        TF1 * fit_gaus2 = new TF1("fit_gaus2", "gaus", -40, 40);
        h1_a_random->Fit(fit_gaus2, "RQ", "", -40, 0);
        h1_a_random->Fit(fit_gaus2, "RQ", "", fit_gaus2->GetParameter(1)-3*fit_gaus2->GetParameter(2), fit_gaus2->GetParameter(1)+0.5*fit_gaus2->GetParameter(2));
        float mu_lhs_a_random = fit_gaus2->GetParameter(1);
        float sigma_lhs_a_random = fit_gaus2->GetParameter(2);
        float mu_lhs_a_random_err = fit_gaus2->GetParError(1);
        float sigma_lhs_a_random_err = fit_gaus2->GetParError(2);


        // fill vectors
        mu_lhs_sub1_basic.push_back(mu_lhs_a_basic);
        mu_lhs_sub1_random.push_back(mu_lhs_a_random);
        sigma_lhs_sub1_basic.push_back(sigma_lhs_a_basic);
        sigma_lhs_sub1_random.push_back(sigma_lhs_a_random);
        mu_lhs_sub1_basic_err.push_back(mu_lhs_a_basic_err);
        mu_lhs_sub1_random_err.push_back(mu_lhs_a_random_err);
        sigma_lhs_sub1_basic_err.push_back(sigma_lhs_a_basic_err);
        sigma_lhs_sub1_random_err.push_back(sigma_lhs_a_random_err);

        mu_sub1_basic.push_back(mean_a_basic);
        mu_sub1_random.push_back(mean_a_random);
        mu_sub1_probe.push_back(mean_a_probe);
        mu_sub1_embed.push_back(mean_a_embed);
        sigma_sub1_basic.push_back(rms_a_basic);
        sigma_sub1_random.push_back(rms_a_random);
        sigma_sub1_probe.push_back(rms_a_probe);
        sigma_sub1_embed.push_back(rms_a_embed);
        mu_sub1_basic_err.push_back(mean_a_basic_err);
        mu_sub1_random_err.push_back(mean_a_random_err);
        mu_sub1_probe_err.push_back(mean_a_probe_err);
        mu_sub1_embed_err.push_back(mean_a_embed_err);
        sigma_sub1_basic_err.push_back(rms_a_basic_err);
        sigma_sub1_random_err.push_back(rms_a_random_err);
        sigma_sub1_probe_err.push_back(rms_a_probe_err);
        sigma_sub1_embed_err.push_back(rms_a_embed_err);

   
        // float ap0 =( mean_a_basic*mean_a_basic);
        // float ab0 = ap0/mean_a_basic;
        TH1F * h1_dummy = (TH1F*)h1_a_random->Clone("h1_dummy");
        h1_dummy->SetMarkerStyle(20);
        h1_dummy->SetMarkerSize(0);
        h1_dummy->SetLineColor(0);
        h1_dummy->SetLineWidth(0);
        h1_dummy->SetFillStyle(0);
        h1_dummy->SetFillColor(0);

        TF1 * fit_gamma = new TF1("fit_gamma", myGammaFunction, -30, 50, 3);
        fit_gamma->SetParNames("A","a_b","a_p");
        fit_gamma->SetParameters(0.1, 1.18, 120);
        // fit_gamma->SetParLimits(0, 0.2, 1.2);
        // fit_gamma->SetParLimits(1, 1, 2.5);
        // fit_gamma->SetParLimits(2, 50, 130);
        h1_a_random->Fit("fit_gamma", "LM+QR", "", 5, 50);
        float a_a_basic = fit_gamma->GetParameter(0);
        float ab_a_basic = fit_gamma->GetParameter(1);
        float ap_a_basic = fit_gamma->GetParameter(2);
        float m = ap_a_basic/ab_a_basic;
        float s = sqrt(ap_a_basic)/ab_a_basic;
    
        TF1 * fit_gamma2 = new TF1("fit_gamma2", myGammaFunction, -30, 30, 3);
        fit_gamma2->SetParameter(0, a_a_basic);
        fit_gamma2->SetParameter(1, ab_a_basic);
        fit_gamma2->SetParameter(2, ap_a_basic);
        // std::cout << "m = " << m << ", s = " << s << std::endl;
        // h1_a_random->Fit(fit_gamma2, "LQR", "", m-3*s, m+3*s);
        a_a_basic = fit_gamma2->GetParameter(0);
        ab_a_basic = fit_gamma2->GetParameter(1);
        ap_a_basic = fit_gamma2->GetParameter(2);
        float chi2overndf = fit_gamma2->GetChisquare()/fit_gamma2->GetNDF();
        //  std::cout << "m = " << m << ", s = " << s << std::endl;

        c = new TCanvas("c", "c", 800, 800);

        if(i == 0) {
            // c->SetLogy();
            TFile * f = new TFile(Form("%s/fit_sub1_cent0_results.root", outdir.c_str()), "RECREATE");
            h1_a_basic->Write();
            h1_a_random->Write();
            h1_a_probe->Write();
            h1_a_embed->Write();
            fit_gaus->Write();
            fit_gaus2->Write();
            fit_gamma->Write();
            f->Close();
            delete f;
        }

        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogy();
        // leg = new TLegend(0.15,0.8,0.4,0.92);
        leg = new TLegend(0.6,0.7,0.85,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(1);
        leg->SetTextSize(0.04);

        TLegend * leg2 = new TLegend(0.17,0.66,0.42,0.8);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetNColumns(1);
        leg2->SetTextSize(0.04);

        h1_sub1_basic[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_sub1_basic[i]->GetXaxis()->SetNdivisions(510);
        h1_sub1_basic[i]->GetYaxis()->SetRangeUser(1e-4, 2e0);
        h1_sub1_basic[i]->SetMarkerSize(1.5);
        h1_sub1_basic[i]->SetMarkerColor(kAzure-2);
        h1_sub1_basic[i]->SetLineColor(kAzure-2);
        h1_sub1_basic[i]->GetXaxis()->SetTitle("#delta E_{T}^{Raw} [GeV]");
        h1_sub1_basic[i]->Draw("P");
        // leg->AddEntry(h1_sub1_basic[i], Form("Basic: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_basic, rms_a_basic), "pe");
        leg->AddEntry(h1_sub1_basic[i], "Basic Cones", "pe");
        for (int j = 0; j < 3; j++)
        {
            leg->AddEntry(h1_dummy, " ", ""); // Add empty entries to keep the legend size consistent
        }

        std::vector<std::string> tags = {sPHENIX_Tag,"Au+Au#kern[0.05]{#sqrt{s_{NN}} = 200 GeV}" ,Form("%d-%d%% Central", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])), "#it{Iterative Method}"};
        TLatex * tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);
        float tx = 0.19;
        float ty_start = 0.87;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty_start, tag.c_str());
            ty_start -= 0.05;
        }

        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/sub1_basic_%d_prelim.png", outdir.c_str(), i));
        leg->Clear();


        fit_gaus->SetLineColor(kAzure-2);
        fit_gaus->SetLineStyle(2);
        fit_gaus->SetLineWidth(2);
        fit_gaus->SetParameter(1, mu_lhs_a_basic);
        fit_gaus->SetParameter(2, sigma_lhs_a_basic);
        // fit_gaus->Draw("SAME");
      
        h1_sub1_random[i]->SetMarkerSize(1.5);
        h1_sub1_random[i]->SetMarkerColor(kRed);
        h1_sub1_random[i]->SetLineColor(kRed);
        h1_sub1_random[i]->Draw(" P SAME");
        // leg->AddEntry(h1_sub1_random[i], Form("Randomized#kern[0.05]{#eta#phi}: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_random, rms_a_random), "pe");
        leg->AddEntry(h1_sub1_basic[i], "Basic Cones", "pe");
        leg->AddEntry(h1_sub1_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        for (int j = 0; j < 2; j++)
        {
            leg->AddEntry(h1_dummy, " ", ""); // Add empty entries to keep the legend size consistent
        }

        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/sub1_basic_%d_prelim_randomized.png", outdir.c_str(), i));
        leg->Clear();
        // leg->AddEntry(fit_gaus2, "Gaussian Fit", "l");

        fit_gaus2->SetParameter(1, mu_lhs_a_random);
        fit_gaus2->SetParameter(2, sigma_lhs_a_random);
        fit_gaus2->SetLineColor(kRed);
        fit_gaus2->SetLineStyle(3);
        fit_gaus2->SetLineWidth(2);
        // fit_gaus2->Draw("SAME");
        
        h1_a_probe->SetMarkerSize(1.5);
        h1_a_probe->SetMarkerColor(kGreen+2);
        h1_a_probe->SetLineColor(kGreen+2);
        h1_a_probe->Draw("P SAME");
        // leg->AddEntry(h1_a_probe, Form("High E_{T} Probe: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_probe, rms_a_probe), "pe");
        leg->AddEntry(h1_sub1_basic[i], "Basic Cones", "pe");
        leg->AddEntry(h1_sub1_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        leg->AddEntry(h1_a_probe, "High E_{T} Probe", "pe");
        leg->AddEntry(h1_dummy, " ", ""); // Add empty entry to keep the legend size consistent
        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/sub1_basic_%d_prelim_probe.png", outdir.c_str(), i));
        leg->Clear();
        h1_a_embed->SetMarkerSize(1.5);
        h1_a_embed->SetMarkerColor(kBlack);
        h1_a_embed->SetLineColor(kBlack);
        h1_a_embed->Draw("SAME");
        // leg->AddEntry(h1_a_embed, Form("Embed: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_embed, rms_a_embed), "pe");
        leg->AddEntry(h1_sub1_basic[i], "Basic Cones", "pe");
        leg->AddEntry(h1_sub1_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        leg->AddEntry(h1_a_probe, "High E_{T} Probe", "pe");
        leg->AddEntry(h1_a_embed, "Embed", "pe");

        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/sub1_basic_%d_prelim_embed.png", outdir.c_str(), i));

        
        // leg->AddEntry(fit_gaus2, Form("#mu^{lhs} = %0.2f GeV, #sigma^{lhs} = %0.2f GeV", mu_lhs_a_random, sigma_lhs_a_random), "l");

        if(do_gamma) {
            fit_gamma2->SetLineColor(kBlack);
            fit_gamma2->SetParameter(0, a_a_basic);
            fit_gamma2->SetParameter(1, ab_a_basic);
            fit_gamma2->SetParameter(2, ap_a_basic);
            fit_gamma2->SetLineStyle(1);
            fit_gamma2->SetLineWidth(2);
            fit_gamma2->Draw("SAME");
            leg->AddEntry(fit_gamma2, Form("f_{#Gamma}: b = %0.2f GeV, p = %0.1f", ab_a_basic, ap_a_basic), "l");
        }
   
       
        leg->Draw("SAME");
        // leg2->Draw("SAME");
        c->SaveAs(Form("%s/sub1_basic_%d_all.png", outdir.c_str(), i));


        delete c;
        delete leg;
        delete fit_gaus;
        delete fit_gamma;
    }


    for ( int i = 0; i < 1; ++i ) {

        TH1F * h1_a_basic = (TH1F*)h1_area_basic[i]->Clone(Form("h1_a_basic_%d", i));
        TH1F * h1_a_random = (TH1F*)h1_area_random[i]->Clone(Form("h1_a_random_%d", i));
        TH1F * h1_a_probe = (TH1F*)h1_area_probe[i]->Clone(Form("h1_a_probe_%d", i));
        TH1F * h1_a_embed = (TH1F*)h1_area_embed[i]->Clone(Form("h1_a_embed_%d", i));
    
        float mean_a_basic = h1_a_basic->GetMean();
        float mean_a_random = h1_a_random->GetMean();
        float mean_a_probe = h1_a_probe->GetMean();
        float mean_a_embed = h1_a_embed->GetMean();
        float rms_a_basic = h1_a_basic->GetRMS();
        float rms_a_random = h1_a_random->GetRMS();
        float rms_a_probe = h1_a_probe->GetRMS();
        float rms_a_embed = h1_a_embed->GetRMS();
        float mean_a_basic_err = h1_a_basic->GetMeanError();
        float mean_a_random_err = h1_a_random->GetMeanError();
        float mean_a_probe_err = h1_a_probe->GetMeanError();
        float mean_a_embed_err = h1_a_embed->GetMeanError();
        float rms_a_basic_err = h1_a_basic->GetRMSError();
        float rms_a_random_err = h1_a_random->GetRMSError();
        float rms_a_probe_err = h1_a_probe->GetRMSError();
        float rms_a_embed_err = h1_a_embed->GetRMSError();
    
       
        
        
    
        TF1 * fit_gaus = new TF1("fit_gaus", "gaus", -40, 40);
        h1_a_basic->Fit(fit_gaus, "RQ", "", -40, 0);
        h1_a_basic->Fit(fit_gaus, "RQ", "", fit_gaus->GetParameter(1)-3*fit_gaus->GetParameter(2), fit_gaus->GetParameter(1)+0.5*fit_gaus->GetParameter(2));
        float mu_lhs_a_basic = fit_gaus->GetParameter(1);
        float sigma_lhs_a_basic = fit_gaus->GetParameter(2);
        float mu_lhs_a_basic_err = fit_gaus->GetParError(1);
        float sigma_lhs_a_basic_err = fit_gaus->GetParError(2);
    
        TF1 * fit_gaus2 = new TF1("fit_gaus2", "gaus", -40, 40);
        h1_a_random->Fit(fit_gaus2, "RQ", "", -40, 0);
        h1_a_random->Fit(fit_gaus2, "RQ", "", fit_gaus2->GetParameter(1)-3*fit_gaus2->GetParameter(2), fit_gaus2->GetParameter(1)+0.5*fit_gaus2->GetParameter(2));
        float mu_lhs_a_random = fit_gaus2->GetParameter(1);
        float sigma_lhs_a_random = fit_gaus2->GetParameter(2);
        float mu_lhs_a_random_err = fit_gaus2->GetParError(1);
        float sigma_lhs_a_random_err = fit_gaus2->GetParError(2);
    
    
        // fill vectors
        mu_lhs_area_basic.push_back(mu_lhs_a_basic);
        mu_lhs_area_random.push_back(mu_lhs_a_random);
        sigma_lhs_area_basic.push_back(sigma_lhs_a_basic);
        sigma_lhs_area_random.push_back(sigma_lhs_a_random);
        mu_lhs_area_basic_err.push_back(mu_lhs_a_basic_err);
        mu_lhs_area_random_err.push_back(mu_lhs_a_random_err);
        sigma_lhs_area_basic_err.push_back(sigma_lhs_a_basic_err);
        sigma_lhs_area_random_err.push_back(sigma_lhs_a_random_err);
    
        mu_area_basic.push_back(mean_a_basic);
        mu_area_random.push_back(mean_a_random);
        mu_area_probe.push_back(mean_a_probe);
        mu_area_embed.push_back(mean_a_embed);
        sigma_area_basic.push_back(rms_a_basic);
        sigma_area_random.push_back(rms_a_random);
        sigma_area_probe.push_back(rms_a_probe);
        sigma_area_embed.push_back(rms_a_embed);
        mu_area_basic_err.push_back(mean_a_basic_err);
        mu_area_random_err.push_back(mean_a_random_err);
        mu_area_probe_err.push_back(mean_a_probe_err);
        mu_area_embed_err.push_back(mean_a_embed_err);
        sigma_area_basic_err.push_back(rms_a_basic_err);
        sigma_area_random_err.push_back(rms_a_random_err);
        sigma_area_probe_err.push_back(rms_a_probe_err);
        sigma_area_embed_err.push_back(rms_a_embed_err);
    
    
        // float ap0 =( mean_a_basic*mean_a_basic);
        // float ab0 = ap0/mean_a_basic;
        TH1F * h1_dummy = (TH1F*)h1_a_random->Clone("h1_dummy");
        h1_dummy->SetMarkerStyle(20);
        h1_dummy->SetMarkerSize(0);
        h1_dummy->SetLineColor(0);
        h1_dummy->SetLineWidth(0);
        h1_dummy->SetFillStyle(0);
        h1_dummy->SetFillColor(0);
    
        TF1 * fit_gamma = new TF1("fit_gamma", myGammaFunction, -30, 50, 3);
        fit_gamma->SetParNames("A","a_b","a_p");
        fit_gamma->SetParameters(0.1, 1.18, 120);
        // fit_gamma->SetParLimits(0, 0.2, 1.2);
        // fit_gamma->SetParLimits(1, 1, 2.5);
        // fit_gamma->SetParLimits(2, 50, 130);
        h1_a_random->Fit("fit_gamma", "LM+QR", "", 5, 50);
        float a_a_basic = fit_gamma->GetParameter(0);
        float ab_a_basic = fit_gamma->GetParameter(1);
        float ap_a_basic = fit_gamma->GetParameter(2);
        float m = ap_a_basic/ab_a_basic;
        float s = sqrt(ap_a_basic)/ab_a_basic;
    
        TF1 * fit_gamma2 = new TF1("fit_gamma2", myGammaFunction, -30, 30, 3);
        fit_gamma2->SetParameter(0, a_a_basic);
        fit_gamma2->SetParameter(1, ab_a_basic);
        fit_gamma2->SetParameter(2, ap_a_basic);
        // std::cout << "m = " << m << ", s = " << s << std::endl;
        // h1_a_random->Fit(fit_gamma2, "LQR", "", m-3*s, m+3*s);
        a_a_basic = fit_gamma2->GetParameter(0);
        ab_a_basic = fit_gamma2->GetParameter(1);
        ap_a_basic = fit_gamma2->GetParameter(2);
        float chi2overndf = fit_gamma2->GetChisquare()/fit_gamma2->GetNDF();
        //  std::cout << "m = " << m << ", s = " << s << std::endl;
    
        c = new TCanvas("c", "c", 800, 800);
    
        if(i == 0) {
            // c->SetLogy();
            TFile * f = new TFile(Form("%s/fit_area_cent0_results.root", outdir.c_str()), "RECREATE");
            h1_a_basic->Write();
            h1_a_random->Write();
            h1_a_probe->Write();
            h1_a_embed->Write();
            fit_gaus->Write();
            fit_gaus2->Write();
            fit_gamma->Write();
            f->Close();
            delete f;
        }
    
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogy();
        // leg = new TLegend(0.15,0.8,0.4,0.92);
        leg = new TLegend(0.6,0.7,0.85,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(1);
        leg->SetTextSize(0.04);
    
        TLegend * leg2 = new TLegend(0.17,0.66,0.42,0.8);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetNColumns(1);
        leg2->SetTextSize(0.04);
    
        h1_area_basic[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_area_basic[i]->GetXaxis()->SetNdivisions(510);
        h1_area_basic[i]->GetYaxis()->SetRangeUser(1e-4, 2e0);
        h1_area_basic[i]->SetMarkerSize(1.5);
        h1_area_basic[i]->SetMarkerColor(kAzure-2);
        h1_area_basic[i]->SetLineColor(kAzure-2);
        h1_area_basic[i]->GetXaxis()->SetTitle("#delta E_{T}^{Raw} [GeV]");
        h1_area_basic[i]->Draw("P");
        // leg->AddEntry(h1_area_basic[i], Form("Basic: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_basic, rms_a_basic), "pe");
        leg->AddEntry(h1_area_basic[i], "Basic Cones", "pe");
        for (int j = 0; j < 3; j++)
        {
            leg->AddEntry(h1_dummy, " ", ""); // Add empty entries to keep the legend size consistent
        }
    
        std::vector<std::string> tags = {sPHENIX_Tag,"Au+Au#kern[0.05]{#sqrt{s_{NN}} = 200 GeV}" ,Form("%d-%d%% Central", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])), "#it{Area Method}"};
        TLatex * tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);
        float tx = 0.19;
        float ty_start = 0.87;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty_start, tag.c_str());
            ty_start -= 0.05;
        }
    
        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/area_basic_%d_prelim.png", outdir.c_str(), i));
        leg->Clear();
    
    
        fit_gaus->SetLineColor(kAzure-2);
        fit_gaus->SetLineStyle(2);
        fit_gaus->SetLineWidth(2);
        fit_gaus->SetParameter(1, mu_lhs_a_basic);
        fit_gaus->SetParameter(2, sigma_lhs_a_basic);
        // fit_gaus->Draw("SAME");
      
        h1_area_random[i]->SetMarkerSize(1.5);
        h1_area_random[i]->SetMarkerColor(kRed);
        h1_area_random[i]->SetLineColor(kRed);
        h1_area_random[i]->Draw(" P SAME");
        // leg->AddEntry(h1_area_random[i], Form("Randomized#kern[0.05]{#eta#phi}: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_random, rms_a_random), "pe");
        leg->AddEntry(h1_area_basic[i], "Basic Cones", "pe");
        leg->AddEntry(h1_area_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        for (int j = 0; j < 2; j++)
        {
            leg->AddEntry(h1_dummy, " ", ""); // Add empty entries to keep the legend size consistent
        }
    
        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/area_basic_%d_prelim_randomized.png", outdir.c_str(), i));
        leg->Clear();
        // leg->AddEntry(fit_gaus2, "Gaussian Fit", "l");
    
        fit_gaus2->SetParameter(1, mu_lhs_a_random);
        fit_gaus2->SetParameter(2, sigma_lhs_a_random);
        fit_gaus2->SetLineColor(kRed);
        fit_gaus2->SetLineStyle(3);
        fit_gaus2->SetLineWidth(2);
        // fit_gaus2->Draw("SAME");
        
        h1_a_probe->SetMarkerSize(1.5);
        h1_a_probe->SetMarkerColor(kGreen+2);
        h1_a_probe->SetLineColor(kGreen+2);
        h1_a_probe->Draw("P SAME");
        // leg->AddEntry(h1_a_probe, Form("High E_{T} Probe: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_probe, rms_a_probe), "pe");
        leg->AddEntry(h1_area_basic[i], "Basic Cones", "pe");
        leg->AddEntry(h1_area_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        leg->AddEntry(h1_a_probe, "High E_{T} Probe", "pe");
        leg->AddEntry(h1_dummy, " ", ""); // Add empty entry to keep the legend size consistent
        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/area_basic_%d_prelim_probe.png", outdir.c_str(), i));
        leg->Clear();
        h1_a_embed->SetMarkerSize(1.5);
        h1_a_embed->SetMarkerColor(kBlack);
        h1_a_embed->SetLineColor(kBlack);
        h1_a_embed->Draw("SAME");
        // leg->AddEntry(h1_a_embed, Form("Embed: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_embed, rms_a_embed), "pe");
        leg->AddEntry(h1_area_basic[i], "Basic Cones", "pe");
        leg->AddEntry(h1_area_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        leg->AddEntry(h1_a_probe, "High E_{T} Probe", "pe");
        leg->AddEntry(h1_a_embed, "Embed", "pe");
    
        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/area_basic_%d_prelim_embed.png", outdir.c_str(), i));
    
        
        // leg->AddEntry(fit_gaus2, Form("#mu^{lhs} = %0.2f GeV, #sigma^{lhs} = %0.2f GeV", mu_lhs_a_random, sigma_lhs_a_random), "l");
    
        if(do_gamma) {
            fit_gamma2->SetLineColor(kBlack);
            fit_gamma2->SetParameter(0, a_a_basic);
            fit_gamma2->SetParameter(1, ab_a_basic);
            fit_gamma2->SetParameter(2, ap_a_basic);
            fit_gamma2->SetLineStyle(1);
            fit_gamma2->SetLineWidth(2);
            fit_gamma2->Draw("SAME");
            leg->AddEntry(fit_gamma2, Form("f_{#Gamma}: b = %0.2f GeV, p = %0.1f", ab_a_basic, ap_a_basic), "l");
        }
    
       
        leg->Draw("SAME");
        // leg2->Draw("SAME");
        c->SaveAs(Form("%s/area_basic_%d_all.png", outdir.c_str(), i));
    
    
        delete c;
        delete leg;
        delete fit_gaus;
        delete fit_gamma;
    }


    for ( int i = 0; i < 1; ++i ) {

        TH1F * h1_a_basic = (TH1F*)h1_mult_basic[i]->Clone(Form("h1_a_basic_%d", i));
        TH1F * h1_a_random = (TH1F*)h1_mult_random[i]->Clone(Form("h1_a_random_%d", i));
        TH1F * h1_a_probe = (TH1F*)h1_mult_probe[i]->Clone(Form("h1_a_probe_%d", i));
        TH1F * h1_a_embed = (TH1F*)h1_mult_embed[i]->Clone(Form("h1_a_embed_%d", i));
    
        float mean_a_basic = h1_a_basic->GetMean();
        float mean_a_random = h1_a_random->GetMean();
        float mean_a_probe = h1_a_probe->GetMean();
        float mean_a_embed = h1_a_embed->GetMean();
        float rms_a_basic = h1_a_basic->GetRMS();
        float rms_a_random = h1_a_random->GetRMS();
        float rms_a_probe = h1_a_probe->GetRMS();
        float rms_a_embed = h1_a_embed->GetRMS();
        float mean_a_basic_err = h1_a_basic->GetMeanError();
        float mean_a_random_err = h1_a_random->GetMeanError();
        float mean_a_probe_err = h1_a_probe->GetMeanError();
        float mean_a_embed_err = h1_a_embed->GetMeanError();
        float rms_a_basic_err = h1_a_basic->GetRMSError();
        float rms_a_random_err = h1_a_random->GetRMSError();
        float rms_a_probe_err = h1_a_probe->GetRMSError();
        float rms_a_embed_err = h1_a_embed->GetRMSError();
    
       
        
        
    
        TF1 * fit_gaus = new TF1("fit_gaus", "gaus", -40, 40);
        h1_a_basic->Fit(fit_gaus, "RQ", "", -40, 0);
        h1_a_basic->Fit(fit_gaus, "RQ", "", fit_gaus->GetParameter(1)-3*fit_gaus->GetParameter(2), fit_gaus->GetParameter(1)+0.5*fit_gaus->GetParameter(2));
        float mu_lhs_a_basic = fit_gaus->GetParameter(1);
        float sigma_lhs_a_basic = fit_gaus->GetParameter(2);
        float mu_lhs_a_basic_err = fit_gaus->GetParError(1);
        float sigma_lhs_a_basic_err = fit_gaus->GetParError(2);
    
        TF1 * fit_gaus2 = new TF1("fit_gaus2", "gaus", -40, 40);
        h1_a_random->Fit(fit_gaus2, "RQ", "", -40, 0);
        h1_a_random->Fit(fit_gaus2, "RQ", "", fit_gaus2->GetParameter(1)-3*fit_gaus2->GetParameter(2), fit_gaus2->GetParameter(1)+0.5*fit_gaus2->GetParameter(2));
        float mu_lhs_a_random = fit_gaus2->GetParameter(1);
        float sigma_lhs_a_random = fit_gaus2->GetParameter(2);
        float mu_lhs_a_random_err = fit_gaus2->GetParError(1);
        float sigma_lhs_a_random_err = fit_gaus2->GetParError(2);
    
    
        // fill vectors
        mu_lhs_mult_basic.push_back(mu_lhs_a_basic);
        mu_lhs_mult_random.push_back(mu_lhs_a_random);
        sigma_lhs_mult_basic.push_back(sigma_lhs_a_basic);
        sigma_lhs_mult_random.push_back(sigma_lhs_a_random);
        mu_lhs_mult_basic_err.push_back(mu_lhs_a_basic_err);
        mu_lhs_mult_random_err.push_back(mu_lhs_a_random_err);
        sigma_lhs_mult_basic_err.push_back(sigma_lhs_a_basic_err);
        sigma_lhs_mult_random_err.push_back(sigma_lhs_a_random_err);
    
        mu_mult_basic.push_back(mean_a_basic);
        mu_mult_random.push_back(mean_a_random);
        mu_mult_probe.push_back(mean_a_probe);
        mu_mult_embed.push_back(mean_a_embed);
        sigma_mult_basic.push_back(rms_a_basic);
        sigma_mult_random.push_back(rms_a_random);
        sigma_mult_probe.push_back(rms_a_probe);
        sigma_mult_embed.push_back(rms_a_embed);
        mu_mult_basic_err.push_back(mean_a_basic_err);
        mu_mult_random_err.push_back(mean_a_random_err);
        mu_mult_probe_err.push_back(mean_a_probe_err);
        mu_mult_embed_err.push_back(mean_a_embed_err);
        sigma_mult_basic_err.push_back(rms_a_basic_err);
        sigma_mult_random_err.push_back(rms_a_random_err);
        sigma_mult_probe_err.push_back(rms_a_probe_err);
        sigma_mult_embed_err.push_back(rms_a_embed_err);
    
    
        // float ap0 =( mean_a_basic*mean_a_basic);
        // float ab0 = ap0/mean_a_basic;
        TH1F * h1_dummy = (TH1F*)h1_a_random->Clone("h1_dummy");
        h1_dummy->SetMarkerStyle(20);
        h1_dummy->SetMarkerSize(0);
        h1_dummy->SetLineColor(0);
        h1_dummy->SetLineWidth(0);
        h1_dummy->SetFillStyle(0);
        h1_dummy->SetFillColor(0);
    
        TF1 * fit_gamma = new TF1("fit_gamma", myGammaFunction, -30, 50, 3);
        fit_gamma->SetParNames("A","a_b","a_p");
        fit_gamma->SetParameters(0.1, 1.18, 120);
        // fit_gamma->SetParLimits(0, 0.2, 1.2);
        // fit_gamma->SetParLimits(1, 1, 2.5);
        // fit_gamma->SetParLimits(2, 50, 130);
        h1_a_random->Fit("fit_gamma", "LM+QR", "", 5, 50);
        float a_a_basic = fit_gamma->GetParameter(0);
        float ab_a_basic = fit_gamma->GetParameter(1);
        float ap_a_basic = fit_gamma->GetParameter(2);
        float m = ap_a_basic/ab_a_basic;
        float s = sqrt(ap_a_basic)/ab_a_basic;
    
        TF1 * fit_gamma2 = new TF1("fit_gamma2", myGammaFunction, -30, 30, 3);
        fit_gamma2->SetParameter(0, a_a_basic);
        fit_gamma2->SetParameter(1, ab_a_basic);
        fit_gamma2->SetParameter(2, ap_a_basic);
        // std::cout << "m = " << m << ", s = " << s << std::endl;
        // h1_a_random->Fit(fit_gamma2, "LQR", "", m-3*s, m+3*s);
        a_a_basic = fit_gamma2->GetParameter(0);
        ab_a_basic = fit_gamma2->GetParameter(1);
        ap_a_basic = fit_gamma2->GetParameter(2);
        float chi2overndf = fit_gamma2->GetChisquare()/fit_gamma2->GetNDF();
        //  std::cout << "m = " << m << ", s = " << s << std::endl;
    
        c = new TCanvas("c", "c", 800, 800);
    
        if(i == 0) {
            // c->SetLogy();
            TFile * f = new TFile(Form("%s/fit_mult_cent0_results.root", outdir.c_str()), "RECREATE");
            h1_a_basic->Write();
            h1_a_random->Write();
            h1_a_probe->Write();
            h1_a_embed->Write();
            fit_gaus->Write();
            fit_gaus2->Write();
            fit_gamma->Write();
            f->Close();
            delete f;
        }
    
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogy();
        // leg = new TLegend(0.15,0.8,0.4,0.92);
        leg = new TLegend(0.6,0.7,0.85,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(1);
        leg->SetTextSize(0.04);
    
        TLegend * leg2 = new TLegend(0.17,0.66,0.42,0.8);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetNColumns(1);
        leg2->SetTextSize(0.04);
    
        h1_mult_basic[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_mult_basic[i]->GetXaxis()->SetNdivisions(510);
        h1_mult_basic[i]->GetYaxis()->SetRangeUser(1e-4, 2e0);
        h1_mult_basic[i]->SetMarkerSize(1.5);
        h1_mult_basic[i]->SetMarkerColor(kAzure-2);
        h1_mult_basic[i]->SetLineColor(kAzure-2);
        h1_mult_basic[i]->GetXaxis()->SetTitle("#delta E_{T}^{Raw} [GeV]");
        h1_mult_basic[i]->Draw("P");
        // leg->AddEntry(h1_mult_basic[i], Form("Basic: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_basic, rms_a_basic), "pe");
        leg->AddEntry(h1_mult_basic[i], "Basic Cones", "pe");
        for (int j = 0; j < 3; j++)
        {
            leg->AddEntry(h1_dummy, " ", ""); // Add empty entries to keep the legend size consistent
        }
    
        std::vector<std::string> tags = {sPHENIX_Tag,"Au+Au#kern[0.05]{#sqrt{s_{NN}} = 200 GeV}" ,Form("%d-%d%% Central", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])), "#it{Multiplicity Method}"};
        TLatex * tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);
        float tx = 0.19;
        float ty_start = 0.87;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty_start, tag.c_str());
            ty_start -= 0.05;
        }
    
        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/mult_basic_%d_prelim.png", outdir.c_str(), i));
        leg->Clear();
    
    
        fit_gaus->SetLineColor(kAzure-2);
        fit_gaus->SetLineStyle(2);
        fit_gaus->SetLineWidth(2);
        fit_gaus->SetParameter(1, mu_lhs_a_basic);
        fit_gaus->SetParameter(2, sigma_lhs_a_basic);
        // fit_gaus->Draw("SAME");
      
        h1_mult_random[i]->SetMarkerSize(1.5);
        h1_mult_random[i]->SetMarkerColor(kRed);
        h1_mult_random[i]->SetLineColor(kRed);
        h1_mult_random[i]->Draw(" P SAME");
        // leg->AddEntry(h1_mult_random[i], Form("Randomized#kern[0.05]{#eta#phi}: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_random, rms_a_random), "pe");
        leg->AddEntry(h1_mult_basic[i], "Basic Cones", "pe");
        leg->AddEntry(h1_mult_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        for (int j = 0; j < 2; j++)
        {
            leg->AddEntry(h1_dummy, " ", ""); // Add empty entries to keep the legend size consistent
        }
    
        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/mult_basic_%d_prelim_randomized.png", outdir.c_str(), i));
        leg->Clear();
        // leg->AddEntry(fit_gaus2, "Gaussian Fit", "l");
    
        fit_gaus2->SetParameter(1, mu_lhs_a_random);
        fit_gaus2->SetParameter(2, sigma_lhs_a_random);
        fit_gaus2->SetLineColor(kRed);
        fit_gaus2->SetLineStyle(3);
        fit_gaus2->SetLineWidth(2);
        // fit_gaus2->Draw("SAME");
        
        h1_a_probe->SetMarkerSize(1.5);
        h1_a_probe->SetMarkerColor(kGreen+2);
        h1_a_probe->SetLineColor(kGreen+2);
        h1_a_probe->Draw("P SAME");
        // leg->AddEntry(h1_a_probe, Form("High E_{T} Probe: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_probe, rms_a_probe), "pe");
        leg->AddEntry(h1_mult_basic[i], "Basic Cones", "pe");
        leg->AddEntry(h1_mult_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        leg->AddEntry(h1_a_probe, "High E_{T} Probe", "pe");
        leg->AddEntry(h1_dummy, " ", ""); // Add empty entry to keep the legend size consistent
        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/mult_basic_%d_prelim_probe.png", outdir.c_str(), i));
        leg->Clear();
        h1_a_embed->SetMarkerSize(1.5);
        h1_a_embed->SetMarkerColor(kBlack);
        h1_a_embed->SetLineColor(kBlack);
        h1_a_embed->Draw("SAME");
        // leg->AddEntry(h1_a_embed, Form("Embed: #mu = %0.2f GeV, #sigma = %0.2f GeV", mean_a_embed, rms_a_embed), "pe");
        leg->AddEntry(h1_mult_basic[i], "Basic Cones", "pe");
        leg->AddEntry(h1_mult_random[i], "Randomized#kern[0.05]{#eta#phi}", "pe");
        leg->AddEntry(h1_a_probe, "High E_{T} Probe", "pe");
        leg->AddEntry(h1_a_embed, "Embed", "pe");
    
        leg->Draw("SAME");
        c->Update();
        c->SaveAs(Form("%s/mult_basic_%d_prelim_embed.png", outdir.c_str(), i));
    
        
        // leg->AddEntry(fit_gaus2, Form("#mu^{lhs} = %0.2f GeV, #sigma^{lhs} = %0.2f GeV", mu_lhs_a_random, sigma_lhs_a_random), "l");
    
        if(do_gamma) {
            fit_gamma2->SetLineColor(kBlack);
            fit_gamma2->SetParameter(0, a_a_basic);
            fit_gamma2->SetParameter(1, ab_a_basic);
            fit_gamma2->SetParameter(2, ap_a_basic);
            fit_gamma2->SetLineStyle(1);
            fit_gamma2->SetLineWidth(2);
            fit_gamma2->Draw("SAME");
            leg->AddEntry(fit_gamma2, Form("f_{#Gamma}: b = %0.2f GeV, p = %0.1f", ab_a_basic, ap_a_basic), "l");
        }
    
       
        leg->Draw("SAME");
        // leg2->Draw("SAME");
        c->SaveAs(Form("%s/mult_basic_%d_all.png", outdir.c_str(), i));
    
    
        delete c;
        delete leg;
        delete fit_gaus;
        delete fit_gamma;
    }


    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {


        c = new TCanvas("c", "c", 800, 800);
      
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogy();

        leg = new TLegend(0.6,0.77,0.85,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(1);
        leg->SetTextSize(0.04);
       
        h1_area_basic[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_area_basic[i]->GetXaxis()->SetNdivisions(510);
        h1_area_basic[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_area_basic[i]->SetMarkerSize(1.5);
        h1_area_basic[i]->SetMarkerColor(kRed);
        h1_area_basic[i]->SetLineColor(kRed);
        h1_area_basic[i]->SetMarkerStyle(kFullCircle);
        h1_area_basic[i]->Draw("P");
        leg->AddEntry(h1_area_basic[i],"Area Method", "pe");
      
        h1_mult_basic[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_mult_basic[i]->GetXaxis()->SetNdivisions(510);
        h1_mult_basic[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_mult_basic[i]->SetMarkerSize(1.5);
        h1_mult_basic[i]->SetMarkerColor(kAzure-2);
        h1_mult_basic[i]->SetLineColor(kAzure-2);
        h1_mult_basic[i]->SetMarkerStyle(kFullSquare);
        h1_mult_basic[i]->Draw("P SAME");
        leg->AddEntry(h1_mult_basic[i], "Multiplicity Method", "pe");

        h1_sub1_basic[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_sub1_basic[i]->GetXaxis()->SetNdivisions(510);
        h1_sub1_basic[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_sub1_basic[i]->SetMarkerSize(1.5);
        h1_sub1_basic[i]->SetMarkerColor(kGreen+2);
        h1_sub1_basic[i]->SetMarkerStyle(kFullTriangleUp);
        h1_sub1_basic[i]->SetLineColor(kGreen+2);
        h1_sub1_basic[i]->Draw("P SAME");
        leg->AddEntry(h1_sub1_basic[i], "Iterative Method", "pe");

        leg->Draw("SAME");
        // std::vector<std::string> tags = {sPHENIX_Tag, "Basic Random Cones",Form("Au+Au %d-%d%%", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1]))};
        std::vector<std::string> tags = {sPHENIX_Tag,"Au+Au#kern[0.05]{#sqrt{s_{NN}} = 200 GeV}" ,Form("%d-%d%% Central", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])),"#it{Basic Random Cones}"};

        TLatex * tex = new TLatex();

        // draw dashed vertical line at 0
        TLine *line = new TLine(0, 1e-4, 0, 1e0);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(kBlack);
        line->Draw("SAME");
        
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.04);
        float tx = 0.19;
        float ty_start = 0.87;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty_start, tag.c_str());
            ty_start -= 0.04;
        }
        c->SaveAs(Form("%s/area_mult_sub1_basic_%d.png", outdir.c_str(), i));


        delete c;
        delete leg;

    }
    
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {


        c = new TCanvas("c", "c", 800, 800);
      
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogy();
    
        leg = new TLegend(0.6,0.77,0.85,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(1);
        leg->SetTextSize(0.04);
       
        h1_area_random[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_area_random[i]->GetXaxis()->SetNdivisions(510);
        h1_area_random[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_area_random[i]->SetMarkerSize(1.5);
        h1_area_random[i]->SetMarkerColor(kRed);
        h1_area_random[i]->SetLineColor(kRed);
        h1_area_random[i]->SetMarkerStyle(kFullCircle);
        h1_area_random[i]->Draw("P");
        leg->AddEntry(h1_area_random[i],"Area Method", "pe");
      
        h1_mult_random[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_mult_random[i]->GetXaxis()->SetNdivisions(510);
        h1_mult_random[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_mult_random[i]->SetMarkerSize(1.5);
        h1_mult_random[i]->SetMarkerColor(kAzure-2);
        h1_mult_random[i]->SetLineColor(kAzure-2);
        h1_mult_random[i]->SetMarkerStyle(kFullSquare);
        h1_mult_random[i]->Draw("P SAME");
        leg->AddEntry(h1_mult_random[i], "Multiplicity Method", "pe");
    
        h1_sub1_random[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_sub1_random[i]->GetXaxis()->SetNdivisions(510);
        h1_sub1_random[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_sub1_random[i]->SetMarkerSize(1.5);
        h1_sub1_random[i]->SetMarkerColor(kGreen+2);
        h1_sub1_random[i]->SetMarkerStyle(kFullTriangleUp);
        h1_sub1_random[i]->SetLineColor(kGreen+2);
        h1_sub1_random[i]->Draw("P SAME");
        leg->AddEntry(h1_sub1_random[i], "Iterative Method", "pe");
    
        leg->Draw("SAME");
        // std::vector<std::string> tags = {sPHENIX_Tag, "Basic Random Cones",Form("Au+Au %d-%d%%", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1]))};
        std::vector<std::string> tags = {sPHENIX_Tag,"Au+Au#kern[0.05]{#sqrt{s_{NN}} = 200 GeV}" ,Form("%d-%d%% Central", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])),"#it{Randomized#kern[0.05]{#eta#phi}}"};
    
        TLatex * tex = new TLatex();
    
        // draw dashed vertical line at 0
        TLine *line = new TLine(0, 1e-4, 0, 1e0);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(kBlack);
        line->Draw("SAME");
        
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.04);
        float tx = 0.19;
        float ty_start = 0.87;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty_start, tag.c_str());
            ty_start -= 0.04;
        }
        c->SaveAs(Form("%s/area_mult_sub1_random_%d.png", outdir.c_str(), i));
    
    
        delete c;
        delete leg;
    
    }

    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {


        c = new TCanvas("c", "c", 800, 800);
      
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogy();
    
        leg = new TLegend(0.6,0.77,0.85,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(1);
        leg->SetTextSize(0.04);
       
        h1_area_probe[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_area_probe[i]->GetXaxis()->SetNdivisions(510);
        h1_area_probe[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_area_probe[i]->SetMarkerSize(1.5);
        h1_area_probe[i]->SetMarkerColor(kRed);
        h1_area_probe[i]->SetLineColor(kRed);
        h1_area_probe[i]->SetMarkerStyle(kFullCircle);
        h1_area_probe[i]->Draw("P");
        leg->AddEntry(h1_area_probe[i],"Area Method", "pe");
      
        h1_mult_probe[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_mult_probe[i]->GetXaxis()->SetNdivisions(510);
        h1_mult_probe[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_mult_probe[i]->SetMarkerSize(1.5);
        h1_mult_probe[i]->SetMarkerColor(kAzure-2);
        h1_mult_probe[i]->SetLineColor(kAzure-2);
        h1_mult_probe[i]->SetMarkerStyle(kFullSquare);
        h1_mult_probe[i]->Draw("P SAME");
        leg->AddEntry(h1_mult_probe[i], "Multiplicity Method", "pe");
    
        h1_sub1_probe[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_sub1_probe[i]->GetXaxis()->SetNdivisions(510);
        h1_sub1_probe[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_sub1_probe[i]->SetMarkerSize(1.5);
        h1_sub1_probe[i]->SetMarkerColor(kGreen+2);
        h1_sub1_probe[i]->SetMarkerStyle(kFullTriangleUp);
        h1_sub1_probe[i]->SetLineColor(kGreen+2);
        h1_sub1_probe[i]->Draw("P SAME");
        leg->AddEntry(h1_sub1_probe[i], "Iterative Method", "pe");
    
        leg->Draw("SAME");
        // std::vector<std::string> tags = {sPHENIX_Tag, "Basic Random Cones",Form("Au+Au %d-%d%%", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1]))};
        std::vector<std::string> tags = {sPHENIX_Tag,"Au+Au#kern[0.05]{#sqrt{s_{NN}} = 200 GeV}" ,Form("%d-%d%% Central", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])),"#it{High E_{T} Probe}"};
    
        TLatex * tex = new TLatex();
    
        // draw dashed vertical line at 0
        TLine *line = new TLine(0, 1e-4, 0, 1e0);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(kBlack);
        line->Draw("SAME");
        
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.04);
        float tx = 0.19;
        float ty_start = 0.87;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty_start, tag.c_str());
            ty_start -= 0.04;
        }
        c->SaveAs(Form("%s/area_mult_sub1_probe_%d.png", outdir.c_str(), i));
    
    
        delete c;
        delete leg;
    
    }
    // for ( int i = 0; i < N_X_CENT_BINS; ++i ) {


    //     c = new TCanvas("c", "c", 800, 800);
      
    //     gPad->SetLeftMargin(0.15);
    //     gPad->SetRightMargin(0.05);
    //     gPad->SetBottomMargin(0.15);
    //     gPad->SetTopMargin(0.05);
    //     gPad->SetLogy();

    //     leg = new TLegend(0.18,0.8,0.4,0.92);
    //     leg->SetBorderSize(0);
    //     leg->SetFillStyle(0);
    //     leg->SetNColumns(1);
    //     leg->SetTextSize(0.04);
       
    //     h1_area_embed[i]->GetXaxis()->SetRangeUser(-30, 30);
    //     h1_area_embed[i]->GetXaxis()->SetNdivisions(510);
    //     h1_area_embed[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
    //     h1_area_embed[i]->SetMarkerSize(1.5);
    //     h1_area_embed[i]->SetMarkerColor(kRed);
    //     h1_area_embed[i]->SetLineColor(kRed);
    //     h1_area_embed[i]->SetMarkerStyle(kFullCircle);
    //     h1_area_embed[i]->Draw("P");
    //     // leg->AddEntry(h1_area_embed[i], Form("Area: #mu = %0.2f GeV, #sigma = %0.2f GeV", h1_area_embed[i]->GetMean(), h1_area_embed[i]->GetRMS()), "pe");
    //     leg->AddEntry(h1_area_embed[i], "Area", "pe");

    //     h1_mult_embed[i]->GetXaxis()->SetRangeUser(-30, 30);
    //     h1_mult_embed[i]->GetXaxis()->SetNdivisions(510);
    //     h1_mult_embed[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
    //     h1_mult_embed[i]->SetMarkerSize(1.5);
    //     h1_mult_embed[i]->SetMarkerColor(kAzure-2);
    //     h1_mult_embed[i]->SetLineColor(kAzure-2);
    //     h1_mult_embed[i]->SetMarkerStyle(kFullSquare);
    //     h1_mult_embed[i]->Draw("P SAME");
    //     // leg->AddEntry(h1_mult_embed[i], Form("Multiplicity: #mu = 
    //     // %0.2f GeV, #sigma = %0.2f GeV", h1_mult_embed[i]->GetMean(), h1_mult_embed[i]->GetRMS()), "pe");
    //     leg->AddEntry(h1_mult_embed[i], "Multiplicity", "pe");
    //     h1_sub1_embed[i]->GetXaxis()->SetRangeUser(-30, 30);
    //     h1_sub1_embed[i]->GetXaxis()->SetNdivisions(510);
    //     h1_sub1_embed[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
    //     h1_sub1_embed[i]->SetMarkerSize(1.5);
    //     h1_sub1_embed[i]->SetMarkerColor(kGreen+2);
    //     h1_sub1_embed[i]->SetMarkerStyle(kFullTriangleUp);
    //     h1_sub1_embed[i]->SetLineColor(kGreen+2);
    //     h1_sub1_embed[i]->Draw("P SAME");
    //     // leg->AddEntry(h1_sub1_embed[i], Form("Iterative: #mu = %0.2f GeV, #sigma = %0.2f GeV", h1_sub1_embed[i]->GetMean(), h1_sub1_embed[i]->GetRMS()), "pe");
    //     leg->AddEntry(h1_sub1_embed[i], "Iterative", "pe");
    //     leg->Draw("SAME");  

    //      // draw dashed vertical line at 0
    //      TLine *line = new TLine(0, 1e-3, 0, 2e-1);
    //      line->SetLineStyle(2);
    //      line->SetLineWidth(2);
    //      line->SetLineColor(kBlack);
    //      line->Draw("SAME");
    //     std::vector<std::string> tags = {sPHENIX_Tag, "Embedded PYTHIA",Form("Au+Au %d-%d%%", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1]))};
    //     TLatex * tex = new TLatex();
    //     tex->SetNDC();
    //     tex->SetTextFont(42);
    //     tex->SetTextSize(0.04);
    //     float tx = 0.19;
    //     float ty_start = 0.75;
    //     for ( auto tag : tags ) {
    //         tex->DrawLatex(tx, ty_start, tag.c_str());
    //         ty_start -= 0.04;
    //     }
    //     c->SaveAs(Form("%s/area_mult_sub1_embed_%d.png", outdir.c_str(), i));


    //     delete c;
    //     delete leg;

    // }
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {


        c = new TCanvas("c", "c", 800, 800);
      
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
        gPad->SetLogy();
    
        leg = new TLegend(0.6,0.77,0.85,0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(1);
        leg->SetTextSize(0.04);
       
        h1_area_embed[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_area_embed[i]->GetXaxis()->SetNdivisions(510);
        h1_area_embed[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_area_embed[i]->SetMarkerSize(1.5);
        h1_area_embed[i]->SetMarkerColor(kRed);
        h1_area_embed[i]->SetLineColor(kRed);
        h1_area_embed[i]->SetMarkerStyle(kFullCircle);
        h1_area_embed[i]->Draw("P");
        leg->AddEntry(h1_area_embed[i],"Area Method", "pe");
      
        h1_mult_embed[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_mult_embed[i]->GetXaxis()->SetNdivisions(510);
        h1_mult_embed[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_mult_embed[i]->SetMarkerSize(1.5);
        h1_mult_embed[i]->SetMarkerColor(kAzure-2);
        h1_mult_embed[i]->SetLineColor(kAzure-2);
        h1_mult_embed[i]->SetMarkerStyle(kFullSquare);
        h1_mult_embed[i]->Draw("P SAME");
        leg->AddEntry(h1_mult_embed[i], "Multiplicity Method", "pe");
    
        h1_sub1_embed[i]->GetXaxis()->SetRangeUser(-30, 30);
        h1_sub1_embed[i]->GetXaxis()->SetNdivisions(510);
        h1_sub1_embed[i]->GetYaxis()->SetRangeUser(1e-4, 1e0);
        h1_sub1_embed[i]->SetMarkerSize(1.5);
        h1_sub1_embed[i]->SetMarkerColor(kGreen+2);
        h1_sub1_embed[i]->SetMarkerStyle(kFullTriangleUp);
        h1_sub1_embed[i]->SetLineColor(kGreen+2);
        h1_sub1_embed[i]->Draw("P SAME");
        leg->AddEntry(h1_sub1_embed[i], "Iterative Method", "pe");
    
        leg->Draw("SAME");
        // std::vector<std::string> tags = {sPHENIX_Tag, "Basic Random Cones",Form("Au+Au %d-%d%%", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1]))};
        std::vector<std::string> tags = {sPHENIX_Tag,"Au+Au#kern[0.05]{#sqrt{s_{NN}} = 200 GeV}" ,Form("%d-%d%% Central", int(X_CENT_BINS[i]), int(X_CENT_BINS[i+1])),"#it{Embedded Jets}"};
    
        TLatex * tex = new TLatex();
    
        // draw dashed vertical line at 0
        TLine *line = new TLine(0, 1e-4, 0, 1e0);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(kBlack);
        line->Draw("SAME");
        
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.04);
        float tx = 0.19;
        float ty_start = 0.87;
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty_start, tag.c_str());
            ty_start -= 0.04;
        }
        c->SaveAs(Form("%s/area_mult_sub1_embed_%d.png", outdir.c_str(), i));
    
    
        delete c;
        delete leg;
    
    }

   
    std::cout << "\\begin{table}[hbt!]" << std::endl;
    std::cout << " \\begin{center}" << std::endl;
    std::cout << "  \\begin{tabular}{|c|c|c|c|c|}" << std::endl;
    std::cout << "   \\hline" << std::endl;
    std::cout << "   & $\\sigma$ & $\\mu$ & $\\sigma^{l.h.s.}$ & $\\mu^{l.h.s.}$ \\\\ "<< std::endl;
    std::cout << "   \\hline \\hline" << std::endl;
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
        std::cout << "   \\hline" << std::endl;
        std::cout << "   \\multicolumn{5}{|c|}{\\textbf{Centrality: " << int(X_CENT_BINS[i]) << "-" << int(X_CENT_BINS[i+1]) << "\\%}} \\\\" << std::endl;
        std::cout << "   \\hline" << std::endl;
        // std::cout << "   Area & " << std::setprecision(4) << mu_area_basic[i] << " & " << sigma_area_basic[i] << " & " << sigma_lhs_area_basic[i] << " & " << mu_lhs_area_basic[i] << " \\\\" << std::endl;
        // std::cout << "   Randomized & "  << std::setprecision(4) << mu_area_random[i] << " & " << sigma_area_random[i] << " & " << sigma_lhs_area_random[i] << " & " << mu_lhs_area_random[i] << " \\\\" << std::endl;
        // std::cout << "   Probe & "  << std::setprecision(4) << mu_area_probe[i] << " & " << sigma_area_probe[i] << " & &  \\\\" << std::endl;
        // std::cout << "   Embeded & "  << std::setprecision(4) << mu_area_embed[i] << " & " << sigma_area_embed[i] << " & &  \\\\" << std::endl;

        // add errors
        std::cout << "   Area & " << std::setprecision(4) << mu_area_basic[i] << " $\\pm$ " << mu_area_basic_err[i] << " & " << sigma_area_basic[i] << " $\\pm$ " << sigma_area_basic_err[i] << " & " << sigma_lhs_area_basic[i] << " $\\pm$ " << sigma_lhs_area_basic_err[i] << " & " << mu_lhs_area_basic[i] << " $\\pm$ " << mu_lhs_area_basic_err[i] << " \\\\" << std::endl;
        std::cout << "   Randomized & "  << std::setprecision(4) << mu_area_random[i] << " $\\pm$ " << mu_area_random_err[i] << " & " << sigma_area_random[i] << " $\\pm$ " << sigma_area_random_err[i] << " & " << sigma_lhs_area_random[i] << " $\\pm$ " << sigma_lhs_area_random_err[i] << " & " << mu_lhs_area_random[i] << " $\\pm$ " << mu_lhs_area_random_err[i] << " \\\\" << std::endl;
        std::cout << "   Probe & "  << std::setprecision(4) << mu_area_probe[i] << " $\\pm$ " << mu_area_probe_err[i] << " & " << sigma_area_probe[i] << " $\\pm$ " << sigma_area_probe_err[i] << " & &  \\\\" << std::endl;
        std::cout << "   Embeded & "  << std::setprecision(4) << mu_area_embed[i] << " $\\pm$ " << mu_area_embed_err[i] << " & " << sigma_area_embed[i] << " $\\pm$ " << sigma_area_embed_err[i] << " & &  \\\\" << std::endl;
        
    }
    std::cout << "   \\hline" << std::endl;
    std::cout << "  \\end{tabular}" << std::endl;
    std::cout << " \\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    std::cout <<std::endl;
    
    std::cout << "\\begin{table}[hbt!]" << std::endl;
    std::cout << " \\begin{center}" << std::endl;
    std::cout << "  \\begin{tabular}{|c|c|c|c|c|}" << std::endl;
    std::cout << "   \\hline" << std::endl;
    std::cout << "   & $\\sigma$ & $\\mu$ & $\\sigma^{l.h.s.}$ & $\\mu^{l.h.s.}$ \\\\ "<< std::endl;
    std::cout << "   \\hline \\hline" << std::endl;
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
        std::cout << "   \\hline" << std::endl;
        std::cout << "   \\multicolumn{5}{|c|}{\\textbf{Centrality: " << int(X_CENT_BINS[i]) << "-" << int(X_CENT_BINS[i+1]) << "\\%}} \\\\" << std::endl;
        std::cout << "   \\hline" << std::endl;
        // std::cout << "   Multiplicity & "  << std::setprecision(4)  << mu_mult_basic[i] << " & " << sigma_mult_basic[i] << " & " << sigma_lhs_mult_basic[i] << " & " << mu_lhs_mult_basic[i] << " \\\\" << std::endl;
        // std::cout << "   Randomized & "  << std::setprecision(4) << mu_mult_random[i] << " & " << sigma_mult_random[i] << " & " << sigma_lhs_mult_random[i] << " & " << mu_lhs_mult_random[i] << " \\\\" << std::endl;
        // std::cout << "   Probe & "  << std::setprecision(4) << mu_mult_probe[i] << " & " << sigma_mult_probe[i] << " & &  \\\\" << std::endl;
        // std::cout << "   Embeded & "  << std::setprecision(4) << mu_mult_embed[i] << " & " << sigma_mult_embed[i] << " & &  \\\\" << std::endl;

        // add errors
        std::cout << "   Multiplicity & "  << std::setprecision(4)  << mu_mult_basic[i] << " $\\pm$ " << mu_mult_basic_err[i] << " & " << sigma_mult_basic[i] << " $\\pm$ " << sigma_mult_basic_err[i] << " & " << sigma_lhs_mult_basic[i] << " $\\pm$ " << sigma_lhs_mult_basic_err[i] << " & " << mu_lhs_mult_basic[i] << " $\\pm$ " << mu_lhs_mult_basic_err[i] << " \\\\" << std::endl;
        std::cout << "   Randomized & "  << std::setprecision(4) << mu_mult_random[i] << " $\\pm$ " << mu_mult_random_err[i] << " & " << sigma_mult_random[i] << " $\\pm$ " << sigma_mult_random_err[i] << " & " << sigma_lhs_mult_random[i] << " $\\pm$ " << sigma_lhs_mult_random_err[i] << " & " << mu_lhs_mult_random[i] << " $\\pm$ " << mu_lhs_mult_random_err[i] << " \\\\" << std::endl;
        std::cout << "   Probe & "  << std::setprecision(4) << mu_mult_probe[i] << " $\\pm$ " << mu_mult_probe_err[i] << " & " << sigma_mult_probe[i] << " $\\pm$ " << sigma_mult_probe_err[i] << " & &  \\\\" << std::endl;
        std::cout << "   Embeded & "  << std::setprecision(4) << mu_mult_embed[i] << " $\\pm$ " << mu_mult_embed_err[i] << " & " << sigma_mult_embed[i] << " $\\pm$ " << sigma_mult_embed_err[i] << " & &  \\\\" << std::endl;

    }
    std::cout << "   \\hline" << std::endl;
    std::cout << "  \\end{tabular}" << std::endl;
    std::cout << " \\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    std::cout <<std::endl;
    
    std::cout << "\\begin{table}[hbt!]" << std::endl;
    std::cout << " \\begin{center}" << std::endl;
    std::cout << "  \\begin{tabular}{|c|c|c|c|c|}" << std::endl;
    std::cout << "   \\hline" << std::endl;
    std::cout << "   & $\\sigma$ & $\\mu$ & $\\sigma^{l.h.s.}$ & $\\mu^{l.h.s.}$ \\ "<< std::endl;
    std::cout << "   \\hline \\hline" << std::endl;
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
        std::cout << "   \\hline" << std::endl;
        std::cout << "   \\multicolumn{5}{|c|}{\\textbf{Centrality: " << int(X_CENT_BINS[i]) << "-" << int(X_CENT_BINS[i+1]) << "\\%}} \\\\" << std::endl;
        std::cout << "   \\hline" << std::endl;
        // std::cout << "   Sub1 & " << std::setprecision(4) << mu_sub1_basic[i] << " & " << sigma_sub1_basic[i] << " & " << sigma_lhs_sub1_basic[i] << " & " << mu_lhs_sub1_basic[i] << " \\\\" << std::endl;
        // std::cout << "   Randomized & "  << std::setprecision(4) << mu_sub1_random[i] << " & " << sigma_sub1_random[i] << " & " << sigma_lhs_sub1_random[i] << " & " << mu_lhs_sub1_random[i] << " \\\\" << std::endl;
        // std::cout << "   Probe & "  << std::setprecision(4) << mu_sub1_probe[i] << " & " << sigma_sub1_probe[i] << " & &  \\\\" << std::endl;
        // std::cout << "   Embeded & "  << std::setprecision(4) << mu_sub1_embed[i] << " & " << sigma_sub1_embed[i] << " & &  \\\\" << std::endl;

        // add errors
        std::cout << "   Sub1 & " << std::setprecision(4) << mu_sub1_basic[i] << " $\\pm$ " << mu_sub1_basic_err[i] << " & " << sigma_sub1_basic[i] << " $\\pm$ " << sigma_sub1_basic_err[i] << " & " << sigma_lhs_sub1_basic[i] << " $\\pm$ " << sigma_lhs_sub1_basic_err[i] << " & " << mu_lhs_sub1_basic[i] << " $\\pm$ " << mu_lhs_sub1_basic_err[i] << " \\\\" << std::endl;
        std::cout << "   Randomized & "  << std::setprecision(4) << mu_sub1_random[i] << " $\\pm$ " << mu_sub1_random_err[i] << " & " << sigma_sub1_random[i] << " $\\pm$ " << sigma_sub1_random_err[i] << " & " << sigma_lhs_sub1_random[i] << " $\\pm$ " << sigma_lhs_sub1_random_err[i] << " & " << mu_lhs_sub1_random[i] << " $\\pm$ " << mu_lhs_sub1_random_err[i] << " \\\\" << std::endl;
        std::cout << "   Probe & "  << std::setprecision(4) << mu_sub1_probe[i] << " $\\pm$ " << mu_sub1_probe_err[i] << " & " << sigma_sub1_probe[i] << " $\\pm$ " << sigma_sub1_probe_err[i] << " & &  \\\\" << std::endl;
        std::cout << "   Embeded & "  << std::setprecision(4) << mu_sub1_embed[i] << " $\\pm$ " << mu_sub1_embed_err[i] << " & " << sigma_sub1_embed[i] << " $\\pm$ " << sigma_sub1_embed_err[i] << " & &  \\\\" << std::endl;

    }
    std::cout << "   \\hline" << std::endl;
    std::cout << "  \\end{tabular}" << std::endl;
    std::cout << " \\end{center}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    std::cout <<std::endl;


 

    

    // std::cout << "  \end{tabular}" << std::endl;
    // std::cout << " \end{center}" << std::endl;
   
   
    TFile * fdebug = new TFile("debug.root", "RECREATE");
    
    
    fdebug->cd();
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
        h1_area_basic[i]->Write();
        h1_area_random[i]->Write();
        h1_area_probe[i]->Write();
        h1_area_embed[i]->Write();
        h1_mult_basic[i]->Write();
        h1_mult_random[i]->Write();
        h1_mult_probe[i]->Write();
        h1_mult_embed[i]->Write();
        h1_sub1_basic[i]->Write();
        h1_sub1_random[i]->Write();
        h1_sub1_probe[i]->Write();
        h1_sub1_embed[i]->Write();
    }
    fdebug->Close();
    

    f_basic->Close();
    f_random->Close();
    f_probe->Close();
    f_embed->Close();
    return;


}



