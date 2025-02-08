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

// const float V2_VALUES[] = {2.32, 3.39, 4.76, 6.18, 7.03, 7.4, 7.44, 7.23, 6.96};
// const float V3_VALUES[] = {1.45, 1.62, 1.76, 1.9, 1.99, 2.05, 1.92, 1.75, 1.57};

const float V2_VALUES[] = { 3.39, 4.76, 6.18, 7.03, 7.4, 7.44, 7.23, 6.96};
const float V3_VALUES[] = { 1.62, 1.76, 1.9, 1.99, 2.05, 1.92, 1.75, 1.57};

const float X_CENT_BINS[]= {0, 10, 20, 30, 40, 50, 60, 70, 80};
const int N_X_CENT_BINS = sizeof(X_CENT_BINS)/sizeof(X_CENT_BINS[0]) - 1;
float MAX_X_CENT = 80;

const float X_COURSE_CENT_BINS[]= {0, 10, 20, 30, 40, 50, 60, 70, 80};
const int N_COURSE_X_BINS = sizeof(X_COURSE_CENT_BINS)/sizeof(X_COURSE_CENT_BINS[0]) - 1;

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
    for ( int i = 0; i < N_SUM_Q_BINS+1; ++i ) { SUM_Q_BINS[i] = i*MAX_SUM_Q/N_SUM_Q_BINS; }
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


const int COLORS[] = {kBlack, kRed , kAzure-2, kGreen+2, kViolet, kCyan, kOrange+2, kMagenta+2, kAzure-2};
const int MARKERS[] = { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kOpenCircle, kOpenSquare, kOpenTriangleUp};
const float MARKER_SIZE = 1.2;
const float LINE_WIDTH = 2.0;

std::string random_cone_plots;

void ConfigureOutputDirs(std::string input_file_base, std::string plotting_dir = "plots/");
std::string MakeGetDir( const std::string & dir ){
    if ( gSystem->AccessPathName(dir.c_str()) ) {
        gSystem->Exec(Form("mkdir -p %s", dir.c_str()));
    }
    return dir;
}

std::string ProcessRandomConeTree(const std::string & input_file, const std::string & prefix, bool x_axis_cent = true);

float CalcPoissonHarm(const float sum2, const float sum, const float n, const float ncomp, const float ncones,  const float v2 = 0 , const float v3 = 0){
    if ( n == 0 ) { return 0;}
    if ( ncones == 0 ) { return 0;}
    float Na = ncomp/ncones;
    // float Na = 1;
    float mu = sum/n;
    float vn2 = v2*v2 + v3*v3;
    float harmcont = Na*Na*vn2*2.0*mu*mu;
    float sigma2 = sum2/n - (mu*mu);
    float sigma = TMath::Sqrt(Na*sigma2  + Na*mu*mu + harmcont);
    return sigma;
}
// g_poission->SetName("g_poission");
// g_poission_v2->SetName("g_poission_v2");
// g_poission_v3->SetName("g_poission_v3");

void MakeSigmaPlot(const std::string & basic_hist, const std::string & rand_hist)
{
    std::string outdir = random_cone_plots;
    TFile * fbasic = new TFile(basic_hist.c_str(), "READ");
    if( !fbasic->IsOpen() || fbasic->IsZombie() ) { std::cout << "File " << basic_hist << " is zombie" << std::endl;  exit(1); }
    TFile * frand = new TFile(rand_hist.c_str(), "READ");
    if( !frand->IsOpen() || frand->IsZombie() ) { std::cout << "File " << rand_hist << " is zombie" << std::endl;  exit(1); }
    

    TGraphErrors * g_area_basic = (TGraphErrors*)fbasic->Get("Area");
    TGraphErrors * g_mult_basic = (TGraphErrors*)fbasic->Get("Multiplicity");
    TGraphErrors * g_sub1_basic = (TGraphErrors*)fbasic->Get("Iterative");
    TGraphErrors * g_poission_basic = (TGraphErrors*)fbasic->Get("g_poission");
    TGraphErrors * g_poission_v2_basic = (TGraphErrors*)fbasic->Get("g_poission_v2");
    TGraphErrors * g_poission_v3_basic = (TGraphErrors*)fbasic->Get("g_poission_v3");
    // TGraphErrors
    g_area_basic->SetName("g_area_basic");
    g_mult_basic->SetName("g_mult_basic");
    g_sub1_basic->SetName("g_sub1_basic");
    g_poission_basic->SetName("g_poission_basic");
    
    g_poission_v2_basic->SetName("g_poission_v2_basic");
    g_poission_v3_basic->SetName("g_poission_v3_basic");

    if (!g_area_basic || !g_mult_basic || !g_sub1_basic || !g_poission_basic || !g_poission_v2_basic || !g_poission_v3_basic) {
        std::cout << "Error: One or more graphs not found in basic histogram file." << std::endl;
        return;
    }
    TGraphErrors * g_area_rand = (TGraphErrors*)frand->Get("Area");
    TGraphErrors * g_mult_rand = (TGraphErrors*)frand->Get("Multiplicity");
    TGraphErrors * g_sub1_rand = (TGraphErrors*)frand->Get("Iterative");
    TGraphErrors * g_poission_rand = (TGraphErrors*)frand->Get("g_poission_sub1");
    TGraphErrors * g_poission_v2_rand = (TGraphErrors*)frand->Get("g_poission_v2");
    TGraphErrors * g_poission_v3_rand = (TGraphErrors*)frand->Get("g_poission_v3");
    if (!g_area_rand || !g_mult_rand || !g_sub1_rand || !g_poission_rand || !g_poission_v2_rand || !g_poission_v3_rand) {
        std::cout << "Error: One or more graphs not found in random histogram file." << std::endl;
        return;
    }

    
    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};

    TCanvas * c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    TLegend * leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(1);
    // leg->SetColumnSeparation(0.001);

    TLegend * leg2 = new TLegend(0.5,0.5,0.9,0.7);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetNColumns(1);

    int colors[] = {kRed, kAzure-2, kGreen+2};
    int markers[] = {kFullCircle, kFullSquare, kFullTriangleUp};
    int markersize = 2;
    int markers_rand[] = {kOpenCircle, kOpenSquare, kOpenTriangleUp};
    std::vector<std::string> labs = {"Area", "Multiplicity", "Iterative"};

    int linestyles_poission[] = {1, 2, 4};
    int linewidth_poission = 3;
    int colors_poission[] = {kBlack, kBlack, kBlack};

    double tx = 0.18;
    double ty = 0.85;

    std::vector<TGraphErrors*> graphs_basic = {g_area_basic, g_mult_basic, g_sub1_basic};
    std::vector<TGraphErrors*> graphs_rand = {g_area_rand, g_mult_rand, g_sub1_rand};
    std::vector<TGraphErrors*> graphs_poission = {g_poission_basic, g_poission_v2_basic, g_poission_v3_basic};
    std::vector<TGraphErrors*> graphs_poission_rand = {g_poission_rand, g_poission_v2_rand, g_poission_v3_rand};
    float xmin = -1, xmax = 80;
    float ymin = 0, ymax = 7;
    std::string xlabel = "Centrality [%]";
    std::string ylabel = "#sigma(#delta E_{T}) [GeV]";
    for ( int i = 0; i < graphs_basic.size(); ++i ) {
        graphs_basic[i]->SetLineColor(colors[i]);
        graphs_basic[i]->SetMarkerColor(colors[i]);
        graphs_basic[i]->SetMarkerStyle(markers[i]);
        graphs_basic[i]->SetMarkerSize(markersize);
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

        graphs_poission[i]->SetLineColor(colors_poission[i]);
        graphs_poission[i]->SetLineStyle(linestyles_poission[i]);
        graphs_poission[i]->SetLineWidth(linewidth_poission);
        graphs_poission[i]->GetXaxis()->SetTitle(xlabel.c_str());
        graphs_poission[i]->GetYaxis()->SetTitle(ylabel.c_str());

        graphs_poission_rand[i]->SetLineColor(colors_poission[i]);
        graphs_poission_rand[i]->SetLineStyle(linestyles_poission[i]);
        graphs_poission_rand[i]->SetLineWidth(linewidth_poission);
        graphs_poission_rand[i]->GetXaxis()->SetTitle(xlabel.c_str());
        graphs_poission_rand[i]->GetYaxis()->SetTitle(ylabel.c_str());

    }

    c->cd();
    std::vector<std::string> labs_poission = {"#sigma_{P}", "#sigma_{P}+#sigma_{NP}(v_{2})", "#sigma_{P}+#sigma_{NP}(v_{2}+v_{3})"};
    int ipick = 1;
    for ( int i = 0; i < graphs_basic.size(); ++i ) {
        if ( i == 0 ) {
            // graphs_basic[i]->GetYaxis()->SetRangeUser(ymin, ymax*1.1);
            // graphs_basic[i]->GetXaxis()->SetRangeUser(xmin, xmax);
            // graphs_basic[i]->Draw("AP");
            graphs_poission[i]->GetYaxis()->SetRangeUser(ymin, ymax*1.1);
            graphs_poission[i]->GetXaxis()->SetRangeUser(xmin, xmax);
            graphs_poission[i]->Draw("AL");
            // graphs_basic[i]->Draw("Psame");
            // graphs_poission[i]->Draw("Lsame");
            // graphs_poission_rand[i]->SetLineColor(kBlue);
            // graphs_poission_rand[i]->Draw("Lsame");
            // leg->AddEntry(graphs_basic[i], " ", "p");
            // leg->AddEntry(graphs_rand[i], labs[i].c_str(), "p");
            leg2->AddEntry(graphs_poission[i], labs_poission[i].c_str(), "l");
        } else {
            // graphs_basic[i]->Draw("Psame");
            // graphs_rand[i]->Draw("Psame");
            // graphs_poission_rand[i]->SetLineColor(kBlue);
            // graphs_poission_rand[i]->Draw("Lsame");
            graphs_poission[i]->Draw("Lsame");
            // leg->AddEntry(graphs_basic[i], " ", "p");
            // leg->AddEntry(graphs_rand[i], labs[i].c_str(), "p");
            leg2->AddEntry(graphs_poission[i], labs_poission[i].c_str(), "l");
        }
        if (i==ipick){
            graphs_basic[i]->Draw("Psame");
            graphs_rand[i]->Draw("Psame");
            leg->AddEntry(graphs_basic[i], "Random Cone", "p");
            leg->AddEntry(graphs_rand[i], "Randomized#eta#phi", "p");
        }
    }

    leg->Draw("same");
    leg2->Draw("same");

    tex->DrawLatex(tx, ty, sPHENIX_Tag.c_str());
    tex->DrawLatex(tx, ty-0.05, DataType_Tag.c_str());
    tex->DrawLatex(tx, ty-0.1, Form("%s",labs[ipick].c_str()));
    tex->DrawLatex(tx, ty-0.15, Form("(WRONG-NOT UPDATED)"));

    c->SaveAs("SigmaPlot.png");


  
    
  



}
void RandomConePlot(const std::string & input_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/dsts/DATA/feb7_basic.root",
                    const std::string & random_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/dsts/DATA/feb7_random.root") 
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

   
    DataType_Tag = "2024 Au+Au 200 GeV";

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
    std::cout << "Random cone plots: " << random_cone_plots << std::endl;
    SetBins();
    
    bool override= true;
    std::string basic_hist;
    std::string rand_hist;
    if (override) {
        basic_hist = ProcessRandomConeTree(input_file, "basic", true);
        rand_hist = ProcessRandomConeTree(random_file, "random", true);
    } else {
        basic_hist = random_cone_plots + "/basic/xaxis_cent/random_cone_plots.root";
        rand_hist = random_cone_plots + "/random/xaxis_cent/random_cone_plots.root";
    }
    MakeSigmaPlot(basic_hist, rand_hist);
    gSystem->Exit(0);
  
   
}

std::string ProcessRandomConeTree(const std::string & input_file, const std::string & prefix ,bool x_axis_cent )
{

    bool x_axis_sum_et = !x_axis_cent;
    std::string outdir = random_cone_plots;
    outdir += "/" + prefix;
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

    TH2F * h2_area_cone_res_vs_x;
    TH2F * h2_mult_cone_res_vs_x;
    TH2F * h2_sub1_cone_res_vs_x;
    TH2F * h2_area_cone_ncomp_vs_x;
    TH2F * h2_mult_cone_ncomp_vs_x;
    TH2F * h2_sub1_cone_ncomp_vs_x;
    // TH1F * h1_raw_cone_et[N_CONE_DET_BINS]; 
    
    if ( x_axis_sum_et ) {

        h2_area_cone_res_vs_x = new TH2F("h2_area_cone_res_vs_x", "h2_area_cone_res_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
        h2_mult_cone_res_vs_x = new TH2F("h2_mult_cone_res_vs_x", "h2_mult_cone_res_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
        h2_sub1_cone_res_vs_x = new TH2F("h2_sub1_cone_res_vs_x", "h2_sub1_cone_res_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
        h2_area_cone_ncomp_vs_x = new TH2F("h2_area_cone_ncomp_vs_x", "h2_area_cone_ncomp_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
        h2_mult_cone_ncomp_vs_x = new TH2F("h2_mult_cone_ncomp_vs_x", "h2_mult_cone_ncomp_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
        h2_sub1_cone_ncomp_vs_x = new TH2F("h2_sub1_cone_ncomp_vs_x", "h2_sub1_cone_ncomp_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_CONECOMP_SUB1_BINS, CONECOMP_SUB1_BINS);

        std::vector<TH2F*> h2s = {h2_area_cone_res_vs_x, h2_mult_cone_res_vs_x, h2_sub1_cone_res_vs_x, h2_area_cone_ncomp_vs_x, h2_mult_cone_ncomp_vs_x, h2_sub1_cone_ncomp_vs_x};
        
        for ( auto h2 : h2s ) {
            h2->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        }
    } else {
        h2_area_cone_res_vs_x = new TH2F("h2_area_cone_res_vs_x", "h2_area_cone_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
        h2_mult_cone_res_vs_x = new TH2F("h2_mult_cone_res_vs_x", "h2_mult_cone_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
        h2_sub1_cone_res_vs_x = new TH2F("h2_sub1_cone_res_vs_x", "h2_sub1_cone_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
        h2_area_cone_ncomp_vs_x = new TH2F("h2_area_cone_ncomp_vs_x", "h2_area_cone_ncomp_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
        h2_mult_cone_ncomp_vs_x = new TH2F("h2_mult_cone_ncomp_vs_x", "h2_mult_cone_ncomp_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
        h2_sub1_cone_ncomp_vs_x = new TH2F("h2_sub1_cone_ncomp_vs_x", "h2_sub1_cone_ncomp_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_CONECOMP_SUB1_BINS, CONECOMP_SUB1_BINS);
        std::vector<TH2F*> h2s = {h2_area_cone_res_vs_x, h2_mult_cone_res_vs_x, h2_sub1_cone_res_vs_x, h2_area_cone_ncomp_vs_x, h2_mult_cone_ncomp_vs_x, h2_sub1_cone_ncomp_vs_x};
        for ( auto h2 : h2s ) {
            h2->GetXaxis()->SetTitle("Centrality [%]");
        }
    }

    h2_area_cone_res_vs_x->GetYaxis()->SetTitle("#delta E_{T}^{Cone} [GeV]");
    h2_mult_cone_res_vs_x->GetYaxis()->SetTitle("#delta E_{T}^{Cone} [GeV]");
    h2_sub1_cone_res_vs_x->GetYaxis()->SetTitle("#delta E_{T}^{Cone} [GeV]");
    h2_area_cone_ncomp_vs_x->GetYaxis()->SetTitle("N_{comp}^{Cone}");
    h2_mult_cone_ncomp_vs_x->GetYaxis()->SetTitle("N_{comp}^{Cone}");
    h2_sub1_cone_ncomp_vs_x->GetYaxis()->SetTitle("N_{comp}^{Cone}");

    const int N_COURSE_X_BINS = 8;
    TH1F * h1_raw_cone_et[N_COURSE_X_BINS+1]; 
    float COURSE_X_BINS[N_COURSE_X_BINS+1];
    float MAX_COURSE_X = 1800;
    if (x_axis_cent){ MAX_COURSE_X = MAX_X_CENT; }
    for ( int i = 0; i < N_COURSE_X_BINS+1; ++i ) { COURSE_X_BINS[i] = i*MAX_COURSE_X/N_COURSE_X_BINS; }
    TH2F * h2_area_cone_res_vs_course_x = new TH2F("h2_area_cone_res_vs_course_x", "h2_area_cone_res_vs_course_x", 
        N_COURSE_X_BINS, COURSE_X_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
    TH2F * h2_mult_cone_res_vs_course_x = new TH2F("h2_mult_cone_res_vs_course_x", "h2_mult_cone_res_vs_course_x",
        N_COURSE_X_BINS, COURSE_X_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
    TH2F * h2_sub1_cone_res_vs_course_x = new TH2F("h2_sub1_cone_res_vs_course_x", "h2_sub1_cone_res_vs_course_x",
        N_COURSE_X_BINS, COURSE_X_BINS, N_CONE_DET_BINS, CONE_DET_BINS);
    
    std::vector<TH2F*> h2s_course_x = {h2_area_cone_res_vs_course_x, h2_mult_cone_res_vs_course_x, h2_sub1_cone_res_vs_course_x};
    for ( auto h2 : h2s_course_x ) {
        if (x_axis_sum_et) {
            h2->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        } else {
            h2->GetXaxis()->SetTitle("Centrality [%]");
        }
        h2->GetYaxis()->SetTitle("#delta E_{T}^{Cone} [GeV]");
    }

    for ( int i = 0; i <= N_COURSE_X_BINS; ++i ) {
        h1_raw_cone_et[i] = new TH1F(Form("h1_raw_cone_et_%d", i), Form("h1_raw_cone_et_%d", i), 120, -10, 110);
        h1_raw_cone_et[i]->GetXaxis()->SetTitle("E_{T}^{Cone} [GeV]");
    }

    
    // N_SUM_ET_BINS, SUM_ET_BINS
    float SUM_ET2[N_X_CENT_BINS];
    float SUM_E[N_X_CENT_BINS];
    float TOTAL_TOWERS_NOT_MASKED[N_X_CENT_BINS];
    float TOTAL_TOWERS_FIRED[N_X_CENT_BINS];
    float N_TOWERS_TOTAL[N_X_CENT_BINS];
    float SUM_CONE_COMPS[N_X_CENT_BINS];
    float N_CONES_THIS_BIN[N_X_CENT_BINS];
    float SUM_ET2_SUB1[N_X_CENT_BINS];
    float SUM_E_SUB1[N_X_CENT_BINS];
    float TOTAL_TOWERS_NOT_MASKED_SUB1[N_X_CENT_BINS];
    float TOTAL_TOWERS_FIRED_SUB1[N_X_CENT_BINS];
    float N_TOWERS_TOTAL_SUB1[N_X_CENT_BINS];
    float SUM_CONE_COMPS_SUB1[N_X_CENT_BINS];
    float N_CONES_THIS_BIN_SUB1[N_X_CENT_BINS];
   

    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {

        SUM_ET2[i] = 0;
        SUM_E[i] = 0;
        TOTAL_TOWERS_NOT_MASKED[i] = 0;
        TOTAL_TOWERS_FIRED[i] = 0;
        N_TOWERS_TOTAL[i] = 0;
        SUM_CONE_COMPS[i] = 0;
        N_CONES_THIS_BIN[i] = 0;

        SUM_ET2_SUB1[i] = 0;
        SUM_E_SUB1[i] = 0;
        TOTAL_TOWERS_NOT_MASKED_SUB1[i] = 0;
        TOTAL_TOWERS_FIRED_SUB1[i] = 0;
        N_TOWERS_TOTAL_SUB1[i] = 0;
        SUM_CONE_COMPS_SUB1[i] = 0;
        N_CONES_THIS_BIN_SUB1[i] = 0;

    }

    float N_CEMC_TOWERS = 256*94;
    float N_HCALIN_TOWERS = 24*64;
    float N_HCALOUT_TOWERS = 24*64;
    for ( int i = 0; i < nentries; ++i ) {
       
        t->GetEntry(i);

        float sum_et = tower_sum_energy_TOWERINFO_CALIB_CEMC + tower_sum_energy_TOWERINFO_CALIB_HCALIN + tower_sum_energy_TOWERINFO_CALIB_HCALOUT;
        float sum_et_sub1 = tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 + tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1 + tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1;
        
        float xaxis_var = sum_et;
        if ( x_axis_cent ) { xaxis_var = 1.0*centrality; }

        int ncomp_cone_corr = random_cone_num_towers_RandomCones_r04 - random_cone_num_masked_towers_RandomCones_r04;
        int ncomp_cone_sub1_corr = random_cone_num_towers_RandomCones_r04_Sub1 - random_cone_num_masked_towers_RandomCones_r04_Sub1;

        float area_bkgd = rho_val_TowerRho_AREA*AREA_CONE;
        float mult_bkgd = rho_val_TowerRho_MULT*random_cone_num_towers_RandomCones_r04;
        // float area_bkgd = ( (rho_val_TowerRho_AREA_HCALIN * AREA_CONE) 
        //                     + (rho_val_TowerRho_AREA_HCALOUT * AREA_CONE) 
        //                     + (rho_val_TowerRho_AREA_CEMC * AREA_CONE) );
        // float mult_bkgd = ( (rho_val_TowerRho_MULT_HCALIN * random_cone_num_towers_hcalin_RandomCones_r04) 
        //                     + (rho_val_TowerRho_MULT_HCALOUT * random_cone_num_towers_hcalout_RandomCones_r04) 
        //                     + (rho_val_TowerRho_MULT_CEMC * random_cone_num_towers_cemc_RandomCones_r04) );

        // float area_bkgd = ( (rho_val_TowerRho_AREA_HCALIN * AREA_CONE) 
        //                     + (rho_val_TowerRho_AREA_HCALOUT * AREA_CONE) );
                            
        // float mult_bkgd = ( (rho_val_TowerRho_MULT_HCALIN * random_cone_num_towers_hcalin_RandomCones_r04) 
        //                     + (rho_val_TowerRho_MULT_HCALOUT * random_cone_num_towers_hcalout_RandomCones_r04) );
                            
                            // float mult_bkgd = rho_val_TowerRho_MULT_HCALIN
        float cone_res_area  = random_cone_energy_RandomCones_r04 - area_bkgd;
        float cone_res_mult  = random_cone_energy_RandomCones_r04 - mult_bkgd;
        // float cone_res_area = random_cone_energy_hcalin_RandomCones_r04 + random_cone_energy_hcalout_RandomCones_r04 - area_bkgd;
        // float cone_res_mult = random_cone_energy_hcalin_RandomCones_r04 + random_cone_energy_hcalout_RandomCones_r04 - mult_bkgd;
        float cone_res_sub1  = random_cone_energy_RandomCones_r04_Sub1;

        
        // nominal
        float sum_et2_cemc = N_CEMC_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_CEMC*tower_std_energy_TOWERINFO_CALIB_CEMC) 
                                + ( tower_avg_energy_TOWERINFO_CALIB_CEMC*tower_avg_energy_TOWERINFO_CALIB_CEMC) );
        float sum_et2_hcalin = N_HCALIN_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_HCALIN*tower_std_energy_TOWERINFO_CALIB_HCALIN)
                                + ( tower_avg_energy_TOWERINFO_CALIB_HCALIN*tower_avg_energy_TOWERINFO_CALIB_HCALIN) );
        float sum_et2_hcalout = N_HCALOUT_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_HCALOUT*tower_std_energy_TOWERINFO_CALIB_HCALOUT)
                                + ( tower_avg_energy_TOWERINFO_CALIB_HCALOUT*tower_avg_energy_TOWERINFO_CALIB_HCALOUT) );
        float total_sum_et2 = sum_et2_cemc + sum_et2_hcalin + sum_et2_hcalout;
        
        float total_fired = ( ( N_CEMC_TOWERS * 1.0 * tower_frac_fired_TOWERINFO_CALIB_CEMC )      
                            + ( N_HCALIN_TOWERS * 1.0 * tower_frac_fired_TOWERINFO_CALIB_HCALIN )
                            + ( N_HCALOUT_TOWERS * 1.0 * tower_frac_fired_TOWERINFO_CALIB_HCALOUT ) );
       
                            // float total_fired = ( N_HCALIN_TOWERS * 1.0 * tower_frac_fired_TOWERINFO_CALIB_HCALIN );
        //                     // + ( N_HCALOUT_TOWERS * 1.0 * tower_frac_fired_TOWERINFO_CALIB_HCALOUT ) );

    
        float total_unmasked = ( N_CEMC_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_CEMC) 
                                + N_HCALIN_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_HCALIN) 
                                + N_HCALOUT_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_HCALOUT) );

        // sub1
        float sum_et2_cemc_sub1 = N_HCALIN_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1*tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1) 
                                + ( tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1*tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1) );
        float sum_et2_hcalin_sub1 = N_HCALIN_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1*tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1)
                                + ( tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1*tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1) );
        float sum_et2_hcalout_sub1 = N_HCALOUT_TOWERS*( (tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1*tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1)
                                + ( tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1*tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1) );
        float total_sum_et2_sub1 = sum_et2_cemc_sub1 + sum_et2_hcalin_sub1 + sum_et2_hcalout_sub1;

        float total_fired_sub1 = ( ( N_HCALIN_TOWERS * tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 )     
                            + ( N_HCALIN_TOWERS * tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1 )
                            + ( N_HCALOUT_TOWERS * tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1 ) );
        float total_unmasked_sub1 = ( N_HCALIN_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1) 
                                + N_HCALIN_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1) 
                                + N_HCALOUT_TOWERS*(1.0 - tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1) );


        // fill arrays
        int THIS_BIN = h2_area_cone_res_vs_x->GetXaxis()->FindBin(xaxis_var);
        THIS_BIN--;// binning starts at 1
        int xbin_course =h2_area_cone_res_vs_course_x->GetXaxis()->FindBin(xaxis_var)-1;
        if (xbin_course < 0) { xbin_course = 0; }
        h1_raw_cone_et[xbin_course]->Fill(random_cone_energy_RandomCones_r04);
        SUM_ET2[THIS_BIN] += total_sum_et2;
        SUM_E[THIS_BIN] += sum_et;
        TOTAL_TOWERS_NOT_MASKED[THIS_BIN] += total_unmasked;
        TOTAL_TOWERS_FIRED[THIS_BIN] += total_fired;
        N_TOWERS_TOTAL[THIS_BIN] += (N_CEMC_TOWERS + N_HCALIN_TOWERS + N_HCALOUT_TOWERS);
        SUM_CONE_COMPS[THIS_BIN] += ncomp_cone_corr;
        N_CONES_THIS_BIN[THIS_BIN]+=1.0;

        
        SUM_ET2_SUB1[THIS_BIN] += total_sum_et2_sub1;
        SUM_E_SUB1[THIS_BIN] += sum_et_sub1;
        TOTAL_TOWERS_NOT_MASKED_SUB1[THIS_BIN] += total_unmasked_sub1;
        TOTAL_TOWERS_FIRED_SUB1[THIS_BIN] += total_fired_sub1;
        N_TOWERS_TOTAL_SUB1[THIS_BIN] += (N_HCALIN_TOWERS + N_HCALIN_TOWERS + N_HCALOUT_TOWERS);
        SUM_CONE_COMPS_SUB1[THIS_BIN] += ncomp_cone_sub1_corr;
        N_CONES_THIS_BIN_SUB1[THIS_BIN]+=1.0;


      
        // fill histograms
        h2_area_cone_res_vs_x->Fill(xaxis_var, cone_res_area);
        h2_mult_cone_res_vs_x->Fill(xaxis_var, cone_res_mult);
        h2_sub1_cone_res_vs_x->Fill(xaxis_var, cone_res_sub1);
        // h2_area_cone_ncomp_vs_x->Fill(xaxis_var, ncomp_cone_corr);
        // h2_mult_cone_ncomp_vs_x->Fill(xaxis_var, cone_res_mult);
        // h2_sub1_cone_ncomp_vs_x->Fill(xaxis_var, ncomp_cone_sub1_corr);
        h2_area_cone_res_vs_course_x->Fill(xaxis_var, cone_res_area);
        h2_mult_cone_res_vs_course_x->Fill(xaxis_var, cone_res_mult);
        h2_sub1_cone_res_vs_course_x->Fill(xaxis_var, cone_res_sub1);


    }


    TGraphErrors * g_poission = new TGraphErrors(N_X_CENT_BINS);
    TGraphErrors * g_poission_v2 = new TGraphErrors(N_X_CENT_BINS);
    TGraphErrors * g_poission_v3 = new TGraphErrors(N_X_CENT_BINS);
    TGraphErrors * g_poission_sub1 = new TGraphErrors(N_X_CENT_BINS);

    for ( int i = 0; i < N_X_CENT_BINS; ++i ) {
       
        float x = h2_area_cone_res_vs_x->GetXaxis()->GetBinCenter(i+1); 
        float xerr = h2_area_cone_res_vs_x->GetXaxis()->GetBinWidth(i+1)/2.0;
        float y = CalcPoissonHarm(SUM_ET2[i], SUM_E[i], TOTAL_TOWERS_FIRED[i], SUM_CONE_COMPS[i], N_CONES_THIS_BIN[i]);
        float y_v2 = CalcPoissonHarm(SUM_ET2[i], SUM_E[i], TOTAL_TOWERS_FIRED[i], SUM_CONE_COMPS[i], N_CONES_THIS_BIN[i], V2_VALUES[i]/100.0);
        float y_v3 = CalcPoissonHarm(SUM_ET2[i], SUM_E[i], TOTAL_TOWERS_FIRED[i], SUM_CONE_COMPS[i], N_CONES_THIS_BIN[i], V2_VALUES[i]/100.0, V3_VALUES[i]/100.0);
        float y_sub1 = CalcPoissonHarm(SUM_ET2_SUB1[i], SUM_E_SUB1[i], TOTAL_TOWERS_FIRED_SUB1[i], SUM_CONE_COMPS_SUB1[i], N_CONES_THIS_BIN_SUB1[i]);
        g_poission->SetPoint(i, x, y);
        g_poission->SetPointError(i, 0, 0);

        g_poission_v2->SetPoint(i, x, y_v2);
        g_poission_v2->SetPointError(i, 0, 0);

        g_poission_v3->SetPoint(i, x, y_v3);
        g_poission_v3->SetPointError(i, 0, 0);

        g_poission_sub1->SetPoint(i, x, y_sub1);
        g_poission_sub1->SetPointError(i, 0, 0);


    }
    g_poission->SetName("g_poission");
    g_poission_v2->SetName("g_poission_v2");
    g_poission_v3->SetName("g_poission_v3");
    g_poission_sub1->SetName("g_poission_sub1");

    
    TCanvas * c;
    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);


    // start with 3x3 of course x bins
    c = new TCanvas("c", "c", 400*3, 400*3);
    c->Divide(3,3);
    double tx=0.45;
    double ty_start=0.85;
    
    TLegend * leg = new TLegend(0.18,0.5,0.35,0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    std::vector<TH2F*> h2s = {h2_area_cone_res_vs_course_x, h2_mult_cone_res_vs_course_x, h2_sub1_cone_res_vs_course_x};
    std::vector<std::string> labs = {"Area", "Multiplicity", "Iterative"};
    for ( unsigned islice = 0; islice < N_COURSE_X_BINS; islice++ ){
        c->cd(islice+1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);

        double miny = 1e-4;
        int ihist = 0;
        std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
        std::string leg_title = Form("%0.0f < #Sigma E_{T}^{Raw} < %0.0f GeV", COURSE_X_BINS[islice], COURSE_X_BINS[islice+1]);
        if ( x_axis_cent ) { leg_title = Form("%0.0f-%0.0f %%", COURSE_X_BINS[islice], COURSE_X_BINS[islice+1]); }
        tags.push_back(leg_title);
        for ( auto h2 : h2s ) {
            h2->GetXaxis()->SetRange(islice+1, islice+1);
            TH1F * h1 = (TH1F*)h2->ProjectionY(Form("h1_%s_%d", h2->GetTitle(), islice), islice+1, islice+1);
            h1->SetLineColor(COLORS[ihist]);
            h1->SetMarkerColor(COLORS[ihist]);
            h1->SetMarkerStyle(MARKERS[ihist]);
            h1->Scale(1.0/h1->Integral());
            h1->GetYaxis()->SetRangeUser(miny, 1e1);
            int lastbin_above_threshold = 0;
            lastbin_above_threshold = h1->FindLastBinAbove(miny);
            // h1->GetXaxis()->SetRangeUser(-1.2*h1->GetBinCenter(lastbin_above_threshold), 1.2*h1->GetBinCenter(lastbin_above_threshold));
            // h1->GetXaxis()->SetRangeUser(-40, 40);
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
    c->SaveAs((outdir+"/cone_res_vs_course_x_slices.png").c_str());

    delete c;
    delete leg;

    const int NX_BINS = h2_area_cone_res_vs_x->GetNbinsX();
    h2s.clear();
    leg = new TLegend(0.18,0.5,0.35,0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    h2s = {h2_area_cone_res_vs_x, h2_mult_cone_res_vs_x, h2_sub1_cone_res_vs_x};
    TGraphErrors * g_std_devs[h2s.size()];
    for ( unsigned ihist = 0; ihist < h2s.size(); ihist++ ) {
        g_std_devs[ihist] = new TGraphErrors(NX_BINS);
    }
    for ( unsigned ibin = 0; ibin < NX_BINS; ibin++ ) {
        c = new TCanvas("c", "c", 3*400, 400);
        c->Divide(3,1);
             
        double miny = 1e-3;
        int ihist = 0;
     
        
        for ( auto h2 : h2s ) {

            c->cd(ihist+1);
            gPad->SetLogy();
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.1);
            gPad->SetBottomMargin(0.15);
            gPad->SetTopMargin(0.05);

            std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
            tags.push_back(labs[ihist]);
            std::string leg_title = Form("%0.0f < #Sigma E_{T}^{Raw} < %0.0f GeV", h2->GetXaxis()->GetBinLowEdge(ibin+1), h2->GetXaxis()->GetBinUpEdge(ibin+1)); 
            if ( x_axis_cent ) { leg_title = Form("%0.0f-%0.0f %%", h2->GetXaxis()->GetBinLowEdge(ibin+1), h2->GetXaxis()->GetBinUpEdge(ibin+1)); }
            tags.push_back(leg_title);

            h2->GetXaxis()->SetRange(ibin+1, ibin+1);
            TH1F * h1 = (TH1F*)h2->ProjectionY(Form("h1_%s_%d", h2->GetTitle(), ibin), ibin+1, ibin+1);
            h1->SetLineColor(COLORS[ihist]);
            h1->SetMarkerColor(COLORS[ihist]);
            h1->SetMarkerStyle(MARKERS[ihist]);
            h1->Scale(1.0/h1->Integral());
            h1->GetYaxis()->SetRangeUser(miny, 1e1);
            int lastbin_above_threshold = 0;
            lastbin_above_threshold = h1->FindLastBinAbove(miny);
            // float abs_max_x = 1.2*h1->GetBinCenter(lastbin_above_threshold);
            float abs_max_x = 40;
            h1->GetXaxis()->SetRangeUser(-abs_max_x, abs_max_x);

            // Fit with Gaussian (do not display)
            float avg = h1->GetMean();
            float std = h1->GetRMS();

            int mean_bin = h1->FindBin(avg);
            h1->GetXaxis()->SetRange(1, mean_bin);
            TF1 * f1 = new TF1("f1", "gaus", -abs_max_x, abs_max_x);
            h1->Fit(f1, "RQ", "", -abs_max_x, avg); // left side of peak
            float mean_left = f1->GetParameter(1);
            float sigma_left = f1->GetParameter(2);
            h1->Fit(f1, "RQ", "", mean_left - 1.5*sigma_left, avg); // left side of peak
            mean_left = f1->GetParameter(1);
            sigma_left = f1->GetParameter(2);
            float mean_left_err = f1->GetParError(1);
            float sigma_left_err = f1->GetParError(2);
            TF1 * fitfunc = h1->GetFunction("f1");
            fitfunc->SetLineColor(kRed);
            fitfunc->SetLineStyle(2);
            fitfunc->SetLineWidth(2);
            float chi2 = fitfunc->GetChisquare();

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
            h1->GetXaxis()->SetRangeUser(-abs_max_x, abs_max_x);
            h1->Draw("p");
            leg->AddEntry(h1, Form("#mu = %0.2f, #sigma = %0.2f", avg, std), "p");
            leg->AddEntry(fitfunc, Form("#mu_{LHS} = %0.2f #pm %0.2f, #sigma_{LHS} = %0.2f #pm %0.2f", mean_left, mean_left_err, sigma_left, sigma_left_err), "l");
            // leg->AddEntry(fitfunc_gamma, Form("a_{b} = %0.2f #pm %0.2f, a_{p} = %0.2f #pm %0.2f", gamma_ab, gamma_ab_err, gamma_ap, gamma_ap_err), "l");

            double ty = ty_start;
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

    c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    TLegend * leg1 = new TLegend(0.18,0.5,0.35,0.7);
    std::vector<TGraphErrors*> g_std_devs_vec = {g_poission, g_poission_v2, g_poission_v3};
    std::vector<std::string> g_std_devs_labels = {"Poisson", "v2", "v3"};
    for ( unsigned ihist = 0; ihist < g_std_devs_vec.size(); ihist++ ) {
        g_std_devs_vec[ihist]->SetMarkerStyle(MARKERS[ihist+1]);
        g_std_devs_vec[ihist]->SetMarkerColor(COLORS[ihist+1]);
        g_std_devs_vec[ihist]->SetLineColor(COLORS[ihist+1]);
        g_std_devs_vec[ihist]->GetYaxis()->SetRangeUser(0.0, 12);
        leg1->AddEntry(g_std_devs_vec[ihist], g_std_devs_labels[ihist].c_str(), "p");
        if ( ihist == 0 ) { g_std_devs_vec[ihist]->Draw("ap"); }
        else { g_std_devs_vec[ihist]->Draw("p same"); }
    }
    leg1->Draw("same");
    c->SaveAs((outdir+"/cone_res_vs_pois.png").c_str());
    delete c;



    TLegend * leg2 = new TLegend(0.18,0.5,0.35,0.7);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    g_poission->SetMarkerStyle(MARKERS[0]);
    g_poission->SetMarkerColor(COLORS[0]);
    g_poission->SetLineColor(COLORS[0]);
    // g_poission->SetLineWidth(0);
    g_poission->SetLineStyle(2);
    // g_pois->GetXaxis()->SetTitle("#Sigma Q_{MBD}");
    if ( x_axis_cent ) { g_poission->GetXaxis()->SetTitle("Centrality [%]"); }
    else { g_poission->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]"); }
    g_poission->GetYaxis()->SetTitle("#sigma(#delta E_{T}) [GeV]");
    g_poission->GetYaxis()->SetRangeUser(0.1, 12);
    // g_poission->GetXaxis()->SetRangeUser(0, 1800);
    g_poission->GetXaxis()->SetNdivisions(505);
    g_poission->Draw("al");
    leg2->AddEntry(g_poission, "Poisson Limit", "pl");

    int mycolors[3] = {kRed, kBlue, kCyan};
    for ( unsigned ihist = 0; ihist < h2s.size(); ihist++ ) {
       
        g_std_devs[ihist]->SetMarkerStyle(MARKERS[ihist+1]);
        g_std_devs[ihist]->SetMarkerColor(mycolors[ihist]);
        g_std_devs[ihist]->SetMarkerSize(1.5);
        g_std_devs[ihist]->GetXaxis()->SetNdivisions(505);
        // g_std_devs[ihist]->GetXaxis()->SetRangeUser(0, 1800);
        g_std_devs[ihist]->GetYaxis()->SetTitle("#sigma_{LHS} [GeV]");
        if ( x_axis_cent ) { g_std_devs[ihist]->GetXaxis()->SetTitle("Centrality [%]"); }
        else { g_std_devs[ihist]->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]"); }
        g_std_devs[ihist]->GetYaxis()->SetRangeUser(0.1, 12);
        // if ( ihist == 0 ) { g_std_devs[ihist]->Draw("ap"); }
        // else { g_std_devs[ihist]->Draw("p same"); }
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
    c->SaveAs((outdir+"/cone_res_vs_sumq.png").c_str());
    delete c;
    TFile * fout = new TFile((outdir+"/random_cone_plots.root").c_str(), "RECREATE");
    for ( auto h2 : h2s ) {
        h2->Write();
    }
    for ( unsigned ihist = 0; ihist < N_COURSE_X_BINS; ihist++ ) {
       h1_raw_cone_et[ihist]->Scale(1.0/h1_raw_cone_et[ihist]->Integral());
         h1_raw_cone_et[ihist]->Write();
    }
    g_poission->Write();
    g_poission_v2->Write();
    g_poission_v3->Write();
    g_poission_sub1->Write();
    for ( unsigned ihist = 0; ihist < h2s.size(); ihist++ ) {
        g_std_devs[ihist]->SetName(labs[ihist].c_str());
        g_std_devs[ihist]->Write();
    }
    fout->Close();
    f->Close();

    

    return outdir+"/random_cone_plots.root";
    
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
    
    random_cone_plots = plotting_dir + "random_cone";
    std::vector<std::string> plot_directories = {random_cone_plots};
    for ( auto const& d : plot_directories ) {
        if ( !gSystem->OpenDirectory(d.c_str()) ) {
            gSystem->mkdir(d.c_str(), true);
        }
    }

    return;
}



