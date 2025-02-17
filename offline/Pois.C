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
const float V3_VALUES[] = {1.43, 1.63, 1.82, 1.94, 2.02, 2.03, 1.96, 1.79, 2};
const float X_CENT_BINS[]= {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
const int N_X_CENT_BINS = sizeof(X_CENT_BINS)/sizeof(X_CENT_BINS[0]) - 1;
float MAX_X_CENT = X_CENT_BINS[N_X_CENT_BINS];

const int N_CONECOMP_BINS = 500;
float CONECOMP_MAX = 1100;
float CONECOMP_BINS[N_CONECOMP_BINS+1];

const int N_CONECOMP_SUB1_BINS = 100;
float CONECOMP_SUB1_MAX = 100;
float CONECOMP_SUB1_BINS[N_CONECOMP_SUB1_BINS+1];

const int N_DET_BINS = 250;
float MAX_DET = 40;
float DET_BINS[N_DET_BINS+1];

void SetBins()
{
    for ( int i = 0; i < N_SUM_ET_BINS+1; ++i ) { SUM_ET_BINS[i] = i*MAX_SUM_ET/N_SUM_ET_BINS; }
    for ( int i = 0; i < N_SUM_Q_BINS+1; ++i ) { SUM_Q_BINS[i] = i*MAX_SUM_Q/N_SUM_Q_BINS; }
    for ( int i = 0; i < N_CONECOMP_BINS+1; ++i ) { CONECOMP_BINS[i] = 500 + i*CONECOMP_MAX/N_CONECOMP_BINS; }
    for ( int i = 0; i < N_CONECOMP_SUB1_BINS+1; ++i ) { CONECOMP_SUB1_BINS[i] = i*CONECOMP_SUB1_MAX/N_CONECOMP_SUB1_BINS; }
    for ( int i = 0; i < N_DET_BINS+1; ++i ) { DET_BINS[i] = -MAX_DET + i*2*MAX_DET/N_DET_BINS; }
}

const float AREA_CONE = TMath::Pi()*0.4*0.4;
const float AREA_TOWER_CEMC = (2.0*TMath::Pi()/256.0)*(2.2/96.0);
const float AREA_HCAL_TOWER = (2.0*TMath::Pi()/64.0)*(2.2/24.0);
const float N_CEMC_TOWERS = 256*94;
const float N_HCALIN_TOWERS = 24*64;
const float N_HCALOUT_TOWERS = 24*64;

bool IS_DATA = false;
int NEVENTS = 0;
std::string DataType_Tag;


const int COLORS[] = {kBlack, kRed , kAzure-2, kCyan, kViolet, kCyan, kOrange+2, kMagenta+2, kAzure-2};
const int MARKERS[] = { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kOpenCircle, kOpenSquare, kOpenTriangleUp};
const float MARKER_SIZE = 1.2;
const float LINE_WIDTH = 2.0;

std::string plot_plots;

void ConfigureOutputDirs(std::string input_file_base, std::string plotting_dir = "plots/");



std::string ProcessTTree(const std::string & input_file, const std::string & prefix);

void MakeSigmaPlot(const std::string & basic_hist, const std::string & rand_hist,  const std::string & pois_hist);

float PoissonEq(const float sigma_et, const float avg_et, const float avg_n,  const float v2 = 0 , const float v3 = 0)
{
    float harm_cont = 2.0*avg_n*avg_n*avg_et*avg_et*( (v2*v2) + (v3*v3) );
    return TMath::Sqrt((avg_n*sigma_et*sigma_et) +(avg_et*avg_et*avg_n) + harm_cont);
}

std::string CalcPois(const std::string & basic_hist, const std::string & histfile="/sphenix/user/tmengel/UE-AuAu-PPG04/offline/vn.root")
{

    std::string outdir = plot_plots;
    TFile * fvn = new TFile(histfile.c_str(), "READ");
    if( !fvn->IsOpen() || fvn->IsZombie() ) { std::cout << "File " << histfile << " is zombie" << std::endl;  exit(1); }
    TGraphErrors * gv2 = (TGraphErrors*)fvn->Get("g_v2");
    TGraphErrors * gv3 = (TGraphErrors*)fvn->Get("g_v3");
    if (!gv2 || !gv3 ) {
        std::cout << "Error: One or more graphs not found in basic histogram file." << std::endl;
        exit(1);
    }
    gv2->SetName("g_v2_phenix");
    gv3->SetName("g_v3_phenix");
    TGraphErrors * gv2_star = new TGraphErrors(N_X_CENT_BINS);
    TGraphErrors * gv3_star = new TGraphErrors(N_X_CENT_BINS);
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
        gv2_star->SetPoint(i, (X_CENT_BINS[i]+X_CENT_BINS[i+1])/2, V2_VALUES[i]/100);
        gv3_star->SetPoint(i, (X_CENT_BINS[i]+X_CENT_BINS[i+1])/2, V3_VALUES[i]/100);
    }
    gv2_star->SetName("g_v2_star");
    gv3_star->SetName("g_v3_star");

    std::vector<float> xbins_v2 {};
    std::vector<float> xcenter_v2 {};
    std::vector<float> xwidth_v2 {};
    std::vector<float> yvalues_v2 {};
    for ( int i = 0; i < gv2->GetN(); ++i ) {
        double x = gv2->GetX()[i];
        double y = gv2->GetY()[i];
        double ex = gv2->GetEX()[i];
        double ey = gv2->GetEY()[i];
        xbins_v2.push_back(x-ex);
        if ( i == gv2->GetN()-1 ) { xbins_v2.push_back(x+ex); }
        xcenter_v2.push_back(x);        
        xwidth_v2.push_back(ex);
        yvalues_v2.push_back(y);
    }

    std::vector<float> xbins_v3 {};
    std::vector<float> xcenter_v3 {};
    std::vector<float> xwidth_v3 {};
    std::vector<float> yvalues_v3 {};
    for ( int i = 0; i < gv3->GetN(); ++i ) {
        double x = gv3->GetX()[i];
        double y = gv3->GetY()[i];
        double ex = gv3->GetEX()[i];
        double ey = gv3->GetEY()[i];
        xbins_v3.push_back(x-ex);
        if ( i == gv3->GetN()-1 ) { xbins_v3.push_back(x+ex); }
        xcenter_v3.push_back(x);        
        xwidth_v3.push_back(ex);
        yvalues_v3.push_back(y);
    }

    TFile * f = new TFile(basic_hist.c_str(), "READ");
    if( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << basic_hist << " is zombie" << std::endl;  exit(1); }

    TTree * t = (TTree*)f->Get("T");
    if (!t) {
        std::cout << "Error: TTree not found in basic histogram file." << std::endl;
        exit(1);
    }

    TH2F * h2_cemc = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC");
    TH2F * h2_hcalin = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN");
    TH2F * h2_hcalout = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT");
    if (!h2_cemc || !h2_hcalin || !h2_hcalout ) {
        std::cout << "Error: One or more histograms not found in basic histogram file." << std::endl;
        exit(1);
    }


    TH2F * h2_sum_et = (TH2F*)h2_cemc->Clone("h2_sum_et");
    h2_sum_et->Add(h2_hcalin);
    h2_sum_et->Add(h2_hcalout);
    float AVG_ET[N_X_CENT_BINS];
    float SIGMA_ET[N_X_CENT_BINS];
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
        h2_sum_et->GetYaxis()->SetRangeUser(X_CENT_BINS[i], X_CENT_BINS[i+1]);
        TH1F * h1 = (TH1F*)h2_sum_et->ProjectionX("h1_tmp");
        AVG_ET[i] = h1->GetMean();
        SIGMA_ET[i] = h1->GetRMS();
        delete h1;
    }
    delete h2_sum_et;
    delete h2_cemc;
    delete h2_hcalin;
    delete h2_hcalout;

    // tree branches 
    int centrality = 0;
    t->SetBranchAddress("centrality", &centrality);
       
    float tower_frac_fired_TOWERINFO_CALIB_CEMC = 0;
    float tower_frac_dead_TOWERINFO_CALIB_CEMC =  0;
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

    float random_cone_eta_RandomCones_r04 = 0;
    int random_cone_num_towers_RandomCones_r04 = 0;
    int random_cone_num_masked_towers_RandomCones_r04 = 0;
    t->SetBranchAddress("random_cone_eta_RandomCones_r04", &random_cone_eta_RandomCones_r04);
    t->SetBranchAddress("random_cone_num_towers_RandomCones_r04", &random_cone_num_towers_RandomCones_r04);
    t->SetBranchAddress("random_cone_num_masked_towers_RandomCones_r04", &random_cone_num_masked_towers_RandomCones_r04);

    TProfile * p2_ncone_comp = new TProfile("p2_Ncone_comp", "p2_Ncone_comp", N_X_CENT_BINS, X_CENT_BINS);

    int nentries = t->GetEntries();
    std::cout << "Processing " << nentries << " events" << std::endl;
    
    float SUM_ET2[N_X_CENT_BINS];
    float SUM_E[N_X_CENT_BINS];
    float TOTAL_TOWERS_FIRED[N_X_CENT_BINS];
    float TOTAL_TOWERS_ALIVE[N_X_CENT_BINS];
    float N_TOWERS_TOTAL[N_X_CENT_BINS];
    float SUM_CONE_COMPS[N_X_CENT_BINS];
    float N_CONES_THIS_BIN[N_X_CENT_BINS];

    float SUM_ET2_V2[xcenter_v2.size()];
    float SUM_E_V2[xcenter_v2.size()];
    float TOTAL_TOWERS_FIRED_V2[xcenter_v2.size()];
    float N_TOWERS_TOTAL_V2[xcenter_v2.size()];
    float SUM_CONE_COMPS_V2[xcenter_v2.size()];
    float N_CONES_THIS_BIN_V2[xcenter_v2.size()];

    float SUM_ET2_V3[xcenter_v3.size()];
    float SUM_E_V3[xcenter_v3.size()];
    float TOTAL_TOWERS_FIRED_V3[xcenter_v3.size()];
    float N_TOWERS_TOTAL_V3[xcenter_v3.size()];
    float SUM_CONE_COMPS_V3[xcenter_v3.size()];
    float N_CONES_THIS_BIN_V3[xcenter_v3.size()];

    TH1F * h1_normal_cent = new TH1F("h1_normal_cent", "h1_normal_cent", N_X_CENT_BINS, X_CENT_BINS);
    TH1F * h1_normal_v2 = new TH1F("h1_normal_v2", "h1_normal_v2", xcenter_v2.size(), xbins_v2.data());
    TH1F * h1_normal_v3 = new TH1F("h1_normal_v3", "h1_normal_v3", xcenter_v3.size(), xbins_v3.data());

    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
        SUM_ET2[i] = 0;
        SUM_E[i] = 0;
        TOTAL_TOWERS_FIRED[i] = 0;
        N_TOWERS_TOTAL[i] = 0;
        SUM_CONE_COMPS[i] = 0;
        N_CONES_THIS_BIN[i] = 0;
        TOTAL_TOWERS_ALIVE[i] = 0;
    }

    for ( int i = 0; i < xcenter_v2.size(); ++i ) {
        SUM_ET2_V2[i] = 0;
        SUM_E_V2[i] = 0;
        TOTAL_TOWERS_FIRED_V2[i] = 0;
        N_TOWERS_TOTAL_V2[i] = 0;
        SUM_CONE_COMPS_V2[i] = 0;
        N_CONES_THIS_BIN_V2[i] = 0;
    }

    for ( int i = 0; i < xcenter_v3.size(); ++i ) {
        SUM_ET2_V3[i] = 0;
        SUM_E_V3[i] = 0;
        TOTAL_TOWERS_FIRED_V3[i] = 0;
        N_TOWERS_TOTAL_V3[i] = 0;
        SUM_CONE_COMPS_V3[i] = 0;
        N_CONES_THIS_BIN_V3[i] = 0;
    }

    float MAX_V2_CENT = xbins_v2.at(xbins_v2.size()-1);
    float MAX_V3_CENT = xbins_v3.at(xbins_v3.size()-1);

    for ( int i = 0; i < nentries; ++i ) {
       
        t->GetEntry(i);

        float xaxis_var =1.0*centrality;
        if ( xaxis_var <0 || xaxis_var > MAX_X_CENT ) { continue; }
        
        float cone_eta = random_cone_eta_RandomCones_r04;
        if ( std::fabs(cone_eta) > 0.6 ) { continue; }
        
        int ncomp_cone_corr = random_cone_num_towers_RandomCones_r04 - random_cone_num_masked_towers_RandomCones_r04;
        p2_ncone_comp->Fill(xaxis_var, ncomp_cone_corr);

        float sum_et = tower_sum_energy_TOWERINFO_CALIB_CEMC + tower_sum_energy_TOWERINFO_CALIB_HCALIN + tower_sum_energy_TOWERINFO_CALIB_HCALOUT; 
        
        float nfired_cemc = N_CEMC_TOWERS*tower_frac_fired_TOWERINFO_CALIB_CEMC;
        float nfired_hcalin = N_HCALIN_TOWERS*tower_frac_fired_TOWERINFO_CALIB_HCALIN;
        float nfired_hcalout = N_HCALOUT_TOWERS*tower_frac_fired_TOWERINFO_CALIB_HCALOUT;
        float total_fired = nfired_cemc + nfired_hcalin + nfired_hcalout;

        float ndead_cemc = N_CEMC_TOWERS*tower_frac_dead_TOWERINFO_CALIB_CEMC;
        float ndead_hcalin = N_HCALIN_TOWERS*tower_frac_dead_TOWERINFO_CALIB_HCALIN;
        float ndead_hcalout = N_HCALOUT_TOWERS*tower_frac_dead_TOWERINFO_CALIB_HCALOUT;
        
        float nalive_cemc = N_CEMC_TOWERS - ndead_cemc;
        float nalive_hcalin = N_HCALIN_TOWERS - ndead_hcalin;
        float nalive_hcalout = N_HCALOUT_TOWERS - ndead_hcalout;

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


        // fill arrays
        int THIS_BIN = h1_normal_cent->FindBin(xaxis_var);
        THIS_BIN--;
        SUM_ET2[THIS_BIN] += total_sum_et2;
        SUM_E[THIS_BIN] += sum_et;
        TOTAL_TOWERS_FIRED[THIS_BIN] += total_fired;
        N_TOWERS_TOTAL[THIS_BIN] += (N_CEMC_TOWERS + N_HCALIN_TOWERS + N_HCALOUT_TOWERS);
        SUM_CONE_COMPS[THIS_BIN] += ncomp_cone_corr;
        N_CONES_THIS_BIN[THIS_BIN]+=1.0;
        TOTAL_TOWERS_ALIVE[THIS_BIN] += (nalive_cemc + nalive_hcalin + nalive_hcalout);

        if ( xaxis_var < MAX_V2_CENT ) {
            int THIS_BIN_V2 = h1_normal_v2->FindBin(xaxis_var);
            THIS_BIN_V2--;
            SUM_ET2_V2[THIS_BIN_V2] += total_sum_et2;
            SUM_E_V2[THIS_BIN_V2] += sum_et;
            TOTAL_TOWERS_FIRED_V2[THIS_BIN_V2] += total_fired;
            N_TOWERS_TOTAL_V2[THIS_BIN_V2] += (N_CEMC_TOWERS + N_HCALIN_TOWERS + N_HCALOUT_TOWERS);
            SUM_CONE_COMPS_V2[THIS_BIN_V2] += ncomp_cone_corr;
            N_CONES_THIS_BIN_V2[THIS_BIN_V2]+=1.0;
        }
        if ( xaxis_var < MAX_V3_CENT ) {
            int THIS_BIN_V3 = h1_normal_v3->FindBin(xaxis_var);
            THIS_BIN_V3--;
            SUM_ET2_V3[THIS_BIN_V3] += total_sum_et2;
            SUM_E_V3[THIS_BIN_V3] += sum_et;
            TOTAL_TOWERS_FIRED_V3[THIS_BIN_V3] += total_fired;
            N_TOWERS_TOTAL_V3[THIS_BIN_V3] += (N_CEMC_TOWERS + N_HCALIN_TOWERS + N_HCALOUT_TOWERS);
            SUM_CONE_COMPS_V3[THIS_BIN_V3] += ncomp_cone_corr;
            N_CONES_THIS_BIN_V3[THIS_BIN_V3]+=1.0;
        }
    }

    
    TGraphErrors * g_poission_total = new TGraphErrors(N_X_CENT_BINS);
    TGraphErrors * g_poission_nfired = new TGraphErrors(N_X_CENT_BINS);
    TGraphErrors * g_poission_nalive = new TGraphErrors(N_X_CENT_BINS);
    TGraphErrors * g_poission_calo_hists = new TGraphErrors(N_X_CENT_BINS);
    TGraphErrors * g_poission_v2_phenix = new TGraphErrors(xcenter_v2.size());
    TGraphErrors * g_poission_v3_phenix = new TGraphErrors(xcenter_v3.size());
    TGraphErrors * g_poission_v2_star = new TGraphErrors(N_X_CENT_BINS);
    TGraphErrors * g_poission_v3_star = new TGraphErrors(N_X_CENT_BINS);

    // float PoissonEq(const float sigma_et, const float avg_et, const float avg_n,  const float v2 = 0 , const float v3 = 0)
    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
       
        float x = h1_normal_cent->GetBinCenter(i+1);

        float mu_fired = SUM_E[i]/TOTAL_TOWERS_FIRED[i];
        float mu_total = SUM_E[i]/N_TOWERS_TOTAL[i];
        float mu_alive = SUM_E[i]/TOTAL_TOWERS_ALIVE[i];
        float mu_hists = AVG_ET[i];

        float avgEt2_fired = SUM_ET2[i]/TOTAL_TOWERS_FIRED[i];
        float avgEt2_total = SUM_ET2[i]/N_TOWERS_TOTAL[i];
        float avgEt2_alive = SUM_ET2[i]/TOTAL_TOWERS_ALIVE[i];

        float sigma_fired = TMath::Sqrt(avgEt2_fired - (mu_fired*mu_fired));
        float sigma_total = TMath::Sqrt(avgEt2_total - (mu_total*mu_total));
        float sigma_alive = TMath::Sqrt(avgEt2_alive - (mu_alive*mu_alive));
        float sigma_hists = SIGMA_ET[i];

        float Na = SUM_CONE_COMPS[i]/N_CONES_THIS_BIN[i];
        float Na_hists = p2_ncone_comp->GetBinContent(i+1);

        float v2_star = gv2_star->GetY()[i];
        float v3_star = gv3_star->GetY()[i];

        float y_fired = PoissonEq(sigma_fired, mu_fired, Na);
        float y_total = PoissonEq(sigma_total, mu_total, Na);
        float y_alive = PoissonEq(sigma_alive, mu_alive, Na);
        float y_hists = PoissonEq(sigma_hists, mu_hists, Na_hists);
        float y_v2 = PoissonEq(sigma_fired, mu_fired, Na, v2_star);
        float y_v3 = PoissonEq(sigma_fired, mu_fired, Na, v2_star, v3_star);

        g_poission_total->SetPoint(i, x, y_total);
        g_poission_nfired->SetPoint(i, x, y_fired);
        g_poission_nalive->SetPoint(i, x, y_alive);
        g_poission_calo_hists->SetPoint(i, x, y_hists);
        g_poission_v2_star->SetPoint(i, x, y_v2);
        g_poission_v3_star->SetPoint(i, x, y_v3);

        g_poission_total->SetPointError(i, 0, 0);
        g_poission_nfired->SetPointError(i, 0, 0);
        g_poission_nalive->SetPointError(i, 0, 0);
        g_poission_calo_hists->SetPointError(i, 0, 0);
        g_poission_v2_star->SetPointError(i, 0, 0);
        g_poission_v3_star->SetPointError(i, 0, 0);
    }

    for ( int i = 0; i < xcenter_v2.size(); ++i ) {
        float x = h1_normal_v2->GetBinCenter(i+1);
        float v2 = gv2->GetY()[i];
        float mu = SUM_E_V2[i]/TOTAL_TOWERS_FIRED_V2[i];
        float avgEt2 = SUM_ET2_V2[i]/TOTAL_TOWERS_FIRED_V2[i];
        float sigma = TMath::Sqrt(avgEt2 - (mu*mu));
        float Na = SUM_CONE_COMPS_V2[i]/N_CONES_THIS_BIN_V2[i];

        float y = PoissonEq(sigma, mu, Na, v2);
        if ( std::isnan(y) ) { y = 0; }
        g_poission_v2_phenix->SetPoint(i, x, y);
        g_poission_v2_phenix->SetPointError(i, 0, 0);
        std::cout << "V2: " << x << " " << y << std::endl;
    }

    for ( int i = 0; i < xcenter_v3.size(); ++i ) {
        float x = h1_normal_v3->GetBinCenter(i+1);
        float v3 = gv3->GetY()[i];
        float v2 = 0;
        for ( int j = 0; j < xcenter_v2.size(); ++j ) {
            if ( xcenter_v2.at(j) > x ) {
                v2 = gv2->GetY()[j];
                break;
            }
        }
        float mu = SUM_E_V3[i]/TOTAL_TOWERS_FIRED_V3[i];
        float avgEt2 = SUM_ET2_V3[i]/TOTAL_TOWERS_FIRED_V3[i];
        float sigma = TMath::Sqrt(avgEt2 - (mu*mu));
        float Na = SUM_CONE_COMPS_V3[i]/N_CONES_THIS_BIN_V3[i];

        float y = PoissonEq(sigma, mu, Na, v2, v3);
        g_poission_v3_phenix->SetPoint(i, x, y);
        g_poission_v3_phenix->SetPointError(i, 0, 0);
        
    }

    std::vector<TGraphErrors*> graphs = {g_poission_total, g_poission_nfired, g_poission_nalive, g_poission_calo_hists, g_poission_v2_star, g_poission_v3_star, g_poission_v2_phenix, g_poission_v3_phenix};
    for ( auto g : graphs ) {
        for ( int i = 0; i < g->GetN(); ++i ) {
            double x, y;
            g->GetPoint(i, x, y);
            if ( y == 0 ) {
                g->RemovePoint(i);
                i--;
            }
        }
    }

    g_poission_total->SetName("g_poission_total");
    g_poission_nfired->SetName("g_poission_nfired");
    g_poission_nalive->SetName("g_poission_nalive");
    g_poission_calo_hists->SetName("g_poission_calo_hists");
    g_poission_v2_star->SetName("g_poission_v2_star");
    g_poission_v3_star->SetName("g_poission_v3_star");
    g_poission_v2_phenix->SetName("g_poission_v2_phenix");
    g_poission_v3_phenix->SetName("g_poission_v3_phenix");


    
    TFile * fpois = new TFile((outdir + "/poisson.root").c_str(), "RECREATE");
    fpois->cd();
    for ( auto g : graphs ) { g->Write(); }
    gv2->Write();
    gv3->Write();
    gv2_star->Write();
    gv3_star->Write();
    p2_ncone_comp->Write();

    fpois->Close();
    
    fvn->Close();
    f->Close();
    
    return outdir + "/poisson.root";
}

void Pois() 
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

    const std::string & input_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/NEW/feb10_basic.root";
    const std::string & random_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/NEW/feb10_random.root";

   
    DataType_Tag = "Au+Au 200 GeV";

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

    
    ConfigureOutputDirs("RandomCones");
    std::cout << "RandomCone plots: " << plot_plots << std::endl;

    SetBins();
    
    // std::string basic_hist = ProcessTTree(input_file, "basic");
    // std::string rand_hist = ProcessTTree(random_file, "random");
    // std::string pois_hist = CalcPois(input_file);
    std::string pois_hist = "/sphenix/user/tmengel/UE-AuAu-PPG04/offline/plots/RandomCones/poisson.root";
    std::string basic_hist ="/sphenix/user/tmengel/UE-AuAu-PPG04/offline/plots/RandomCones/basic/cones.root";
    std::string rand_hist = "/sphenix/user/tmengel/UE-AuAu-PPG04/offline/plots/RandomCones/random/cones.root";
    MakeSigmaPlot(basic_hist, rand_hist, pois_hist);
    gSystem->Exit(0);
  
   
}

std::string ProcessTTree(const std::string & input_file, const std::string & prefix  )
{

    std::string outdir = plot_plots;
    outdir += "/" + prefix;
    if( !gSystem->OpenDirectory(outdir.c_str()) ) {
        gSystem->mkdir(outdir.c_str(), true);
    }

    TFile * f = new TFile(input_file.c_str(), "READ");
    if( !f->IsOpen() || f->IsZombie() ) { std::cout << "File " << input_file << " is zombie" << std::endl;  exit(1); }

    TTree * t = (TTree*)f->Get("T");
    
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
    std::cout << "Processing " << nentries << " events" << std::endl;
  
    TH2F * h2_area_res_vs_x = new TH2F("h2_area_res_vs_x", "h2_area_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_DET_BINS, DET_BINS);
    TH2F * h2_mult_res_vs_x = new TH2F("h2_mult_res_vs_x", "h2_mult_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_DET_BINS, DET_BINS);
    TH2F * h2_sub1_res_vs_x = new TH2F("h2_sub1_res_vs_x", "h2_sub1_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_DET_BINS, DET_BINS);
  
    std::vector<TH2F*> h2s = {h2_area_res_vs_x, h2_mult_res_vs_x, h2_sub1_res_vs_x};
    for ( auto h2 : h2s ) {
        h2->GetXaxis()->SetTitle("Centrality [%]");
        h2->GetYaxis()->SetTitle("#delta E_{T}^{Cone} [GeV]");
    }
    std::vector<std::string> labs = {"Area", "Multiplicity", "Iterative"};


    for ( int i = 0; i < nentries; ++i ) {
       
        t->GetEntry(i);
        float xaxis_var =1.0*centrality;
        if ( xaxis_var <0 || xaxis_var > MAX_X_CENT ) { continue; }
        
        int ncomp_cone_corr = random_cone_num_towers_RandomCones_r04 - random_cone_num_masked_towers_RandomCones_r04;
        int ncomp_cone_corr_sub1 = random_cone_num_towers_RandomCones_r04_Sub1 - random_cone_num_masked_towers_RandomCones_r04_Sub1;

        float area_bkgd = rho_val_TowerRho_AREA*AREA_CONE;
        // float mult_bkgd = rho_val_TowerRho_MULT*ncomp_cone_corr;
        float mult_bkgd = ( rho_val_TowerRho_MULT_CEMC*(random_cone_num_towers_cemc_RandomCones_r04 - random_cone_num_masked_towers_cemc_RandomCones_r04) 
                            + rho_val_TowerRho_MULT_HCALIN*(random_cone_num_towers_hcalin_RandomCones_r04 - random_cone_num_masked_towers_hcalin_RandomCones_r04) 
                            + rho_val_TowerRho_MULT_HCALOUT*(random_cone_num_towers_hcalout_RandomCones_r04-random_cone_num_masked_towers_hcalout_RandomCones_r04) );
        float cone_res_area  = random_cone_energy_RandomCones_r04 - area_bkgd;
        float cone_res_mult  = random_cone_energy_cemc_RandomCones_r04 + random_cone_energy_hcalin_RandomCones_r04 + random_cone_energy_hcalout_RandomCones_r04 - mult_bkgd;
        // float cone_res_mult  = random_cone_energy_RandomCones_r04 - mult_bkgd;
        float cone_res_sub1  = random_cone_energy_RandomCones_r04_Sub1;

        float cone_eta = random_cone_eta_RandomCones_r04;
        float cone_sub1_eta = random_cone_eta_RandomCones_r04_Sub1;
        if ( std::fabs(cone_eta) > 0.6 || std::fabs(cone_sub1_eta) > 0.6 ) { continue; }
        

        h2_area_res_vs_x->Fill(xaxis_var, cone_res_area);
        h2_mult_res_vs_x->Fill(xaxis_var, cone_res_mult);
        h2_sub1_res_vs_x->Fill(xaxis_var, cone_res_sub1);

        
    }
    
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

        double miny = 1e-5;
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
            // h1->GetYaxis()->SetRangeUser(miny, 1e0);
            h1->GetYaxis()->SetTitle("Probability Density [A.U.]");
            h1->GetXaxis()->SetRangeUser(-30,30);
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
    c->SaveAs((outdir+"/cone_res_vs_x_slices.png").c_str());

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
            

            float avg = h1->GetMean();
            float std = h1->GetRMS();
            h1->GetYaxis()->SetRangeUser(miny, 1e0);
            h1->GetXaxis()->SetRangeUser(-30,30);
            h1->GetYaxis()->SetTitle("Probability Density [A.U.]");
            // int mean_bin = h1->FindBin(avg);
            // // h1->GetXaxis()->SetRange(1, mean_bin);
            // TF1 * f1 = new TF1("f1", "gaus", -MAX_DET, avg);
            // h1->Fit(f1, "RQ", "", -MAX_DET, avg); // left side of peak
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

        c->SaveAs((outdir+"/cone_res_vs_x_slice_"+std::to_string(ibin)+".png").c_str());
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
        // g_std_devs[ihist]->GetYaxis()->SetRangeUser(0.1, 12);
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
    c->SaveAs((outdir+"/cone_res_vs_x.png").c_str());
    delete c;
    TFile * fout = new TFile((outdir+"/cones.root").c_str(), "RECREATE");
    for ( auto h2 : h2s ) {
        h2->Write();
    }
    for ( unsigned ihist = 0; ihist < h2s.size(); ihist++ ) {
        g_std_devs[ihist]->SetName(labs[ihist].c_str());
        g_std_devs[ihist]->Write();
    }
    
    fout->Close();
    f->Close();

    return outdir+"/cones.root";
    
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
    
    plot_plots = plotting_dir;
    std::vector<std::string> plot_directories = {plot_plots};
    for ( auto const& d : plot_directories ) {
        if ( !gSystem->OpenDirectory(d.c_str()) ) {
            gSystem->mkdir(d.c_str(), true);
        }
    }

    return;
}

void MakeSigmaPlot(const std::string & basic_hist, const std::string & rand_hist, const std::string & pois_hist)
{
    std::string outdir = plot_plots;
    TFile * fbasic = new TFile(basic_hist.c_str(), "READ");
    if( !fbasic->IsOpen() || fbasic->IsZombie() ) { std::cout << "File " << basic_hist << " is zombie" << std::endl;  exit(1); }
    TFile * frand = new TFile(rand_hist.c_str(), "READ");
    if( !frand->IsOpen() || frand->IsZombie() ) { std::cout << "File " << rand_hist << " is zombie" << std::endl;  exit(1); }
    std::cout << "Opened files " << basic_hist << " and " << rand_hist << std::endl;
    TFile * fpois = new TFile(pois_hist.c_str(), "READ");
    if( !fpois->IsOpen() || fpois->IsZombie() ) { std::cout << "File " << pois_hist << " is zombie" << std::endl;  exit(1); }

    TGraphErrors * g_area_basic = (TGraphErrors*)fbasic->Get("Area");
    TGraphErrors * g_mult_basic = (TGraphErrors*)fbasic->Get("Multiplicity");
    TGraphErrors * g_sub1_basic = (TGraphErrors*)fbasic->Get("Iterative");
    g_area_basic->SetName("g_area_basic");
    g_mult_basic->SetName("g_mult_basic");
    g_sub1_basic->SetName("g_sub1_basic");

    if (!g_area_basic || !g_mult_basic || !g_sub1_basic) {
        std::cout << "Error: One or more graphs not found in basic histogram file." << std::endl;
        return;
    }
  

    TGraphErrors * g_poission_basic = (TGraphErrors*)fpois->Get("g_poission_total");
    TGraphErrors * g_poission_v2_basic = (TGraphErrors*)fpois->Get("g_poission_v2_star");
    TGraphErrors * g_poission_v3_basic = (TGraphErrors*)fpois->Get("g_poission_v3_star");
    TGraphErrors * g_poission_fired = (TGraphErrors*)fpois->Get("g_poission_nfired");
    TGraphErrors * g_poission_dead = (TGraphErrors*)fpois->Get("g_poission_nalive");
    TGraphErrors * g_poission_calo_hists = (TGraphErrors*)fpois->Get("g_poission_calo_hists");
    TGraphErrors * g_poission_v2_phenix = (TGraphErrors*)fpois->Get("g_poission_v2_phenix");
    TGraphErrors * g_poission_v3_phenix = (TGraphErrors*)fpois->Get("g_poission_v3_phenix");
    std::vector<TGraphErrors*> graphs_pois = {g_poission_basic, g_poission_v2_basic, g_poission_v3_basic, g_poission_fired, g_poission_dead, g_poission_calo_hists, g_poission_v2_phenix, g_poission_v3_phenix};
    int ig =0;
    for ( auto g : graphs_pois ) {
        if ( !g ) { std::cout << "Error: Graph not found in poission histogram file." << std::endl; return; }
        g->SetLineColor(COLORS[ig]);
        g->SetLineStyle(ig);
        ig++;
    }
    
    if (!g_poission_basic || !g_poission_v2_basic || !g_poission_v3_basic) {
        std::cout << "Error: One or more graphs not found in basic histogram file." << std::endl;
        return;
    }


    TGraphErrors * g_area_rand = (TGraphErrors*)frand->Get("Area");
    TGraphErrors * g_mult_rand = (TGraphErrors*)frand->Get("Multiplicity");
    TGraphErrors * g_sub1_rand = (TGraphErrors*)frand->Get("Iterative");
    if (!g_area_rand || !g_mult_rand || !g_sub1_rand) {
        std::cout << "Error: One or more graphs not found in random histogram file." << std::endl;
        return;
    }

    
    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};

    TCanvas * c;
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    TLegend * leg;
    TLegend * leg2;

    int colors[] = {kRed, kAzure-2, kGreen+2};
    int markers[] = {kFullCircle, kFullSquare, kFullTriangleUp};
    float markersize = 1.7;
    int markers_rand[] = {kOpenCircle, kOpenSquare, kOpenTriangleUp};
    std::vector<std::string> labs = {"Area", "Multiplicity", "Iterative"};

    int linestyles_poission[] = {1, 2, 10};
    int linewidth_poission = 3;
    int colors_poission[] = {kBlack, kBlack, kBlack};

    double tx = 0.18;
    double ty = 0.25;

    std::vector<TGraphErrors*> graphs_basic = {g_area_basic, g_mult_basic, g_sub1_basic};
    std::vector<TGraphErrors*> graphs_rand = {g_area_rand, g_mult_rand, g_sub1_rand};
    std::vector<TGraphErrors*> graphs_poission = {g_poission_basic, g_poission_v2_basic, g_poission_v3_basic};
    float xmin = -1, xmax = 80;
    float ymin = 0, ymax = 7;
    std::string xlabel = "Centrality [%]";
    std::string ylabel = "#sigma(#delta E_{T}) [GeV]";
    for ( int i = 0; i < graphs_basic.size(); ++i ) {
        graphs_basic[i]->SetLineColor(colors[i]);
        graphs_basic[i]->SetMarkerColor(colors[i]);
        graphs_basic[i]->SetLineWidth(2);
        graphs_basic[i]->SetMarkerStyle(markers[i]);
        graphs_basic[i]->SetMarkerSize(markersize);
        // graphs_basic[i]->SetMarkerFillStyle(0);
        // graphs_basic[i]->SetMarkerColorAlpha(colors[i], 0.5);

        for (int j = 0; j < graphs_basic[i]->GetN(); ++j) {
            if ( graphs_basic[i]->GetY()[j] > ymax ) { ymax = graphs_basic[i]->GetY()[j]; }
            // set xerr and yer to 0
            graphs_basic[i]->SetPointError(j, 0, 0);
        }
        graphs_basic[i]->GetXaxis()->SetTitle(xlabel.c_str());
        graphs_basic[i]->GetYaxis()->SetTitle(ylabel.c_str());
        

        graphs_rand[i]->SetLineColor(colors[i]);
        graphs_rand[i]->SetMarkerColor(colors[i]);
        graphs_rand[i]->SetMarkerStyle(markers_rand[i]);
        graphs_rand[i]->SetMarkerSize(markersize);
        for (int j = 0; j < graphs_rand[i]->GetN(); ++j) {

            if ( graphs_rand[i]->GetY()[j] > ymax ) { ymax = graphs_rand[i]->GetY()[j]; }

            // set xerr and yer to 0
            graphs_rand[i]->SetPointError(j, 0, 0);
        }
        graphs_rand[i]->GetXaxis()->SetTitle(xlabel.c_str());
        graphs_rand[i]->GetYaxis()->SetTitle(ylabel.c_str());     

    }

    c = new TCanvas("c", "c", 800, 800);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    leg = new TLegend(0.18,0.8,0.89,0.91);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(2);
    leg->SetTextSize(0.03);
    ig = 0;
    for (auto g : graphs_pois) {
        g->GetYaxis()->SetRangeUser(0.9, 6.1);
        g->GetXaxis()->SetRangeUser(xmin, xmax);
        if ( ig == 0 ) { g->Draw("AL"); }
        else { g->Draw("Lsame"); }
        leg->AddEntry(g, Form("%s",g->GetName()), "l");
        ig++;
    }
    for ( int i = 0; i < graphs_basic.size(); ++i ){
        graphs_basic[i]->GetYaxis()->SetRangeUser(0.9, 6.1);
        graphs_basic[i]->GetXaxis()->SetRangeUser(xmin, xmax);
        graphs_basic[i]->Draw("Psame");
        graphs_rand[i]->Draw("Psame");
        leg->AddEntry(graphs_basic[i], "Random Cones", "p");
        leg->AddEntry(graphs_rand[i], "Randomized #eta,#phi", "p");
    } 
    leg->Draw("same");
    c->SaveAs((outdir+"/pois_debug.png").c_str());
    delete c;
    delete leg;



    for ( int i = 0; i < graphs_poission.size(); ++i ) {
        graphs_poission[i]->SetLineColor(colors_poission[i]);
        graphs_poission[i]->SetLineStyle(linestyles_poission[i]);
        graphs_poission[i]->SetLineWidth(linewidth_poission);
        graphs_poission[i]->GetXaxis()->SetTitle(xlabel.c_str());
        graphs_poission[i]->GetYaxis()->SetTitle(ylabel.c_str());
    }
    
    for ( int i = 0; i < graphs_basic.size(); ++i ) {
        
        tx = 0.18;
        ty = 0.3;

        c = new TCanvas("c", "c", 800, 800);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);

        leg = new TLegend(0.38,0.8,0.89,0.91);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetNColumns(1);

        leg2 = new TLegend(0.43,0.64,0.89,0.79);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetNColumns(1);
        
        graphs_basic[i]->GetYaxis()->SetRangeUser(0.9,6.1);
        graphs_basic[i]->GetXaxis()->SetRangeUser(xmin, xmax);
        graphs_basic[i]->Draw("AP");
        graphs_rand[i]->Draw("Psame");
        leg->AddEntry(graphs_basic[i],"Random Cones", "p");
        leg->AddEntry(graphs_rand[i], "Randomized #eta,#phi", "p");
        g_poission_basic->Draw("Lsame");
        g_poission_v2_basic->Draw("Lsame");
        g_poission_v3_basic->Draw("Lsame");
        leg2->AddEntry(g_poission_basic, "#sigma_{P}", "l");
        leg2->AddEntry(g_poission_v2_basic, "#sigma_{P}#oplus#sigma_{NP}(v_{2})", "l");
        leg2->AddEntry(g_poission_v3_basic, "#sigma_{P}#oplus#sigma_{NP}(v_{2}+v_{3})", "l");

        leg->Draw("same");
        leg2->Draw("same");
        for ( auto tag : tags ) {
            tex->DrawLatex(tx, ty, tag.c_str());
            ty -= 0.05;
        }
        tex->DrawLatex(tx, ty, labs[i].c_str());

        c->SaveAs((outdir+"/sigma_et_vs_centrality_"+labs[i]+".png").c_str());
        delete c;
        delete leg;
        delete leg2;
    }
    
    tx = 0.49;
    ty = 0.88;
    c = new TCanvas("c", "c", 800, 800);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    leg = new TLegend(0.18,0.18,0.4,0.35);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(2);
    leg->SetColumnSeparation(0.01);
    
    
    leg2 = new TLegend(0.5,0.58,0.89,0.72);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetNColumns(1);
    std::vector<std::string> labs_poission = {"#sigma_{P}", "#sigma_{P}+#sigma_{NP}(v_{2})", "#sigma_{P}+#sigma_{NP}(v_{2})+#sigma_{NP}(v_{3})"};
    for ( int i = 0; i < graphs_basic.size(); ++i ) {
        if ( i == 0 ) {
            graphs_basic[i]->GetYaxis()->SetRangeUser(0.9, 6.1);
            graphs_basic[i]->GetYaxis()->SetTitle("#sigma(#delta E_{T})");
            graphs_basic[i]->GetXaxis()->SetRangeUser(xmin, xmax);
            graphs_basic[i]->Draw("AP");
            graphs_rand[i]->Draw("Psame");
            leg->AddEntry(graphs_basic[i], " ", "P");
            leg->AddEntry(graphs_rand[i], Form(" %s", labs[i].c_str()), "P");
        } else {
            graphs_basic[i]->Draw("Psame");
            graphs_rand[i]->Draw("Psame");
            leg->AddEntry(graphs_basic[i], " ", "P");
            leg->AddEntry(graphs_rand[i], Form(" %s", labs[i].c_str()), "P");
        }
    }
    g_poission_basic->Draw("Lsame");
    g_poission_v2_basic->Draw("Lsame");
    g_poission_v3_basic->Draw("Lsame");
    leg2->AddEntry(g_poission_basic, "#sigma_{P}", "l");
    leg2->AddEntry(g_poission_v2_basic, "#sigma_{P}#oplus#sigma_{NP}(v_{2})", "l");
    leg2->AddEntry(g_poission_v3_basic, "#sigma_{P}#oplus#sigma_{NP}(v_{2}+v_{3})", "l");

    leg->Draw("same");
    leg2->Draw("same");

    tex->DrawLatex(tx, ty, sPHENIX_Tag.c_str());
    tex->DrawLatex(tx, ty-0.05, DataType_Tag.c_str());
    tex->SetTextSize(0.035);
    tex->DrawLatex(tx, ty-0.1, "#it{Open Points: Randomized#eta,#phi}");
    tex->DrawLatex(tx, ty-0.15, "#it{Closed Points: Random Cones}");
    c->SaveAs((outdir+"/sigma_et_vs_centrality.png").c_str());
    delete c;
    delete leg2;
    delete leg;
    delete tex;


    TGraphErrors* g_poission_v2_over_basic = (TGraphErrors*)g_poission_basic->Clone("g_poission_v2_over_basic");
    TGraphErrors* g_poission_v3_over_basic = (TGraphErrors*)g_poission_basic->Clone("g_poission_v3_over_basic");
    TGraphErrors * g_poission_basic_over_basic = (TGraphErrors*)g_poission_basic->Clone("g_poission_basic_over_basic");
    TGraphErrors * g_area_basic_div = (TGraphErrors*)g_area_basic->Clone("g_area_basic_div");
    TGraphErrors * g_mult_basic_div = (TGraphErrors*)g_mult_basic->Clone("g_mult_basic_div");
    TGraphErrors * g_sub1_basic_div = (TGraphErrors*)g_sub1_basic->Clone("g_sub1_basic_div");
    TGraphErrors * g_area_rand_div = (TGraphErrors*)g_area_rand->Clone("g_area_rand_div");
    TGraphErrors * g_mult_rand_div = (TGraphErrors*)g_mult_rand->Clone("g_mult_rand_div");
    TGraphErrors * g_sub1_rand_div = (TGraphErrors*)g_sub1_rand->Clone("g_sub1_rand_div");
    for ( int i = 0; i < g_poission_basic->GetN(); ++i ) {
        double x = g_poission_basic->GetX()[i];
        
        double y = g_poission_basic->GetY()[i];
        double y_ab = g_area_basic->GetY()[i];
        double y_mb = g_mult_basic->GetY()[i];
        double y_sb = g_sub1_basic->GetY()[i];
        double y_ar = g_area_rand->GetY()[i];
        double y_mr = g_mult_rand->GetY()[i];
        double y_sr = g_sub1_rand->GetY()[i];
        double y_v2 = g_poission_v2_basic->GetY()[i];
        double y_v3 = g_poission_v3_basic->GetY()[i];
        g_poission_v2_over_basic->SetPoint(i, x, y_v2/y);
        g_poission_v3_over_basic->SetPoint(i, x, y_v3/y);
        g_poission_basic_over_basic->SetPoint(i, x, y/y);
        g_poission_v2_over_basic->SetPointError(i, 0, 0);
        g_poission_v3_over_basic->SetPointError(i, 0, 0);
        g_poission_basic_over_basic->SetPointError(i, 0, 0);
        g_area_basic_div->SetPoint(i, x, y_ab/y);
        g_mult_basic_div->SetPoint(i, x, y_mb/y);
        g_sub1_basic_div->SetPoint(i, x, y_sb/y);
        g_area_rand_div->SetPoint(i, x, y_ar/y);
        g_mult_rand_div->SetPoint(i, x, y_mr/y);
        g_sub1_rand_div->SetPoint(i, x, y_sr/y);
        g_area_basic_div->SetPointError(i, 0, 0);
        g_mult_basic_div->SetPointError(i, 0, 0);
        g_sub1_basic_div->SetPointError(i, 0, 0);
        g_area_rand_div->SetPointError(i, 0, 0);
        g_mult_rand_div->SetPointError(i, 0, 0);
        g_sub1_rand_div->SetPointError(i, 0, 0);
    }

    graphs_poission = {g_poission_basic_over_basic, g_poission_v2_over_basic, g_poission_v3_over_basic};
    graphs_rand = {g_area_rand_div, g_mult_rand_div, g_sub1_rand_div};
    graphs_basic = {g_area_basic_div, g_mult_basic_div, g_sub1_basic_div};

    for ( int i = 0; i < graphs_poission.size(); ++i ) {
        graphs_poission[i]->SetLineColor(colors_poission[i]);
        graphs_poission[i]->SetLineStyle(linestyles_poission[i]);
        graphs_poission[i]->SetLineWidth(linewidth_poission);
        graphs_poission[i]->GetXaxis()->SetTitle(xlabel.c_str());
        graphs_poission[i]->GetYaxis()->SetTitle(ylabel.c_str());
    }

    c = new TCanvas("c", "c", 800, 800);
    leg = new TLegend(0.18,0.78,0.3,0.91);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(2);
    leg->SetColumnSeparation(0.01);
    
    
    leg2 = new TLegend(0.5,0.18,0.8,0.32);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetNColumns(1);
    tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    for ( int i = 0; i < graphs_basic.size(); ++i ) {



        if ( i == 0 ) {
            graphs_basic[i]->GetYaxis()->SetRangeUser(0.7, 1.7);
            graphs_basic[i]->GetYaxis()->SetTitle("#sigma(#delta E_{T}) / #sigma_{P}");
            graphs_basic[i]->GetXaxis()->SetRangeUser(xmin, xmax);
            graphs_basic[i]->Draw("AP");
            graphs_rand[i]->Draw("Psame");
            g_poission_basic_over_basic->Draw("Lsame");
            g_poission_v2_over_basic->Draw("Lsame");
            g_poission_v3_over_basic->Draw("Lsame");
            leg->AddEntry(graphs_basic[i], " ", "P");
            leg->AddEntry(graphs_rand[i], Form(" %s", labs[i].c_str()), "P");


        } else {
            graphs_basic[i]->Draw("Psame");
            graphs_rand[i]->Draw("Psame");
            leg->AddEntry(graphs_basic[i], " ", "P");
            leg->AddEntry(graphs_rand[i], Form(" %s", labs[i].c_str()), "P");
        }
       
    }

    // leg2->AddEntry(g_poission_basic_over_basic, "#sigma_{P}/#sigma_{P}", "l");
    leg2->AddEntry(g_poission_v2_over_basic, "#sigma_{P}#oplus#sigma_{NP}(v_{2})/#sigma_{P}", "l");
    leg2->AddEntry(g_poission_v3_over_basic, "#sigma_{P}#oplus#sigma_{NP}(v_{2}+v_{3})/#sigma_{P}", "l");
    leg->Draw("same");
    leg2->Draw("same");

    tex->DrawLatex(tx, ty, sPHENIX_Tag.c_str());
    tex->DrawLatex(tx, ty-0.05, DataType_Tag.c_str());
    tex->SetTextSize(0.035);
    tex->DrawLatex(tx, ty-0.1, "#it{Open Points: Randomized#eta,#phi}");
    tex->DrawLatex(tx, ty-0.15, "#it{Closed Points: Random Cones}");

    c->SaveAs((outdir+"/sigma_et_vs_centrality_oversigma.png").c_str());
    delete c;
    delete leg;
    delete leg2;

    fbasic->Close();
    frand->Close();
    fpois->Close();
    return;

}


