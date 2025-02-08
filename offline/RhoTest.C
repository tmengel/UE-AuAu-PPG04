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


const float CENT_BINS[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 100};
const int N_CENT_BINS = sizeof(CENT_BINS)/sizeof(CENT_BINS[0]) - 1;

const int N_SUM_ET_BINS = 50;
float SUM_ET_BINS[N_SUM_ET_BINS+1];
float MAX_SUM_ET = 2200;

const int N_SUM_Q_BINS = 50;
float SUM_Q_BINS[N_SUM_Q_BINS+1];
float MAX_SUM_Q = 2200;

const int N_X_CENT_BINS = 18;
float X_CENT_BINS[N_X_CENT_BINS+1];
float MAX_X_CENT = 80;

const int N_ZVTX_BINS = 100;
float ZVTX_BINS[N_ZVTX_BINS+1];
float MAX_ZVTX = 25;

const int N_RHO_M_BINS = 200;
float RHO_M_MAX = 0.1;
float RHO_M_BINS[N_RHO_M_BINS+1];

const int N_RHO_A_BINS = 200;
float RHO_A_MAX = 150;
float RHO_A_BINS[N_RHO_A_BINS+1];



void SetBins(){
    for ( int i = 0; i < N_SUM_ET_BINS+1; ++i ) { SUM_ET_BINS[i] = i*MAX_SUM_ET/N_SUM_ET_BINS; }
    for ( int i = 0; i < N_X_CENT_BINS+1; ++i ) { X_CENT_BINS[i] = i*MAX_X_CENT/N_X_CENT_BINS; }
    for ( int i = 0; i < N_SUM_Q_BINS+1; ++i ) { SUM_Q_BINS[i] = i*MAX_SUM_Q/N_SUM_Q_BINS; }
    for ( int i = 0; i < N_ZVTX_BINS+1; ++i ) { ZVTX_BINS[i] = -MAX_ZVTX + i*2*MAX_ZVTX/N_ZVTX_BINS; }
    for ( int i = 0; i < N_RHO_M_BINS+1; ++i ) { RHO_M_BINS[i] = i*RHO_M_MAX/N_RHO_M_BINS; }
    for ( int i = 0; i < N_RHO_A_BINS+1; ++i ) { RHO_A_BINS[i] = i*RHO_A_MAX/N_RHO_A_BINS; }
}

const float AREA_CONE = TMath::Pi()*0.4*0.4;

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


void PlotRho(){
    std::string fin_cent = "/sphenix/user/tmengel/UE-AuAu-PPG04/offline/plots/RHOTEST_DATA/background_density/xaxis_cent/BackgroundEstimates.root";
    std::string fin_sum_et = "/sphenix/user/tmengel/UE-AuAu-PPG04/offline/plots/RHOTEST_DATA/background_density/xaxis_sum_et/BackgroundEstimates.root";

    TFile * f_cent = new TFile(fin_cent.c_str(), "READ");
    TFile * f_sum_et = new TFile(fin_sum_et.c_str(), "READ");

    TH2F * h2_rho_new_cent = (TH2F*)f_cent->Get("h2_rho_area_vs_x");
    TH2F * h2_rho_old_cent= (TH2F*)f_cent->Get("h2_rho_area_alt_vs_x");
    TH2F * h2_rhoA_new_cent = (TH2F*)f_cent->Get("h2_rho_area_times_A_vs_x");
    TH2F * h2_rhoA_old_cent = (TH2F*)f_cent->Get("h2_rho_area_alt_times_A_vs_x");
    h2_rho_new_cent->SetName("h2_rho_new_cent");
    h2_rho_old_cent->SetName("h2_rho_old_cent");
    h2_rhoA_new_cent->SetName("h2_rhoA_new_cent");
    h2_rhoA_old_cent->SetName("h2_rhoA_old_cent");


    TH2F * h2_rho_new_et = (TH2F*)f_sum_et->Get("h2_rho_area_vs_x");
    TH2F * h2_rho_old_et= (TH2F*)f_sum_et->Get("h2_rho_area_alt_vs_x");
    TH2F * h2_rhoA_new_et = (TH2F*)f_sum_et->Get("h2_rho_area_times_A_vs_x");
    TH2F * h2_rhoA_old_et = (TH2F*)f_sum_et->Get("h2_rho_area_alt_times_A_vs_x");

    TCanvas * c;

    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    double tx = 0.19;
    double ty = 0.85;
    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};

    c = new TCanvas("c_rho_new_cent", "c_rho_new_cent", 1000, 400);
    c->Divide(2,1);
    for ( int i = 0; i < 2; ++i ) {
        c->cd(i+1);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.05);
    }
    c->cd(1);
    h2_rho_new_et->GetXaxis()->SetNdivisions(505);
    h2_rho_new_et->Draw("colz");
    tex->DrawLatex(tx, ty, tags[0].c_str());
    tex->DrawLatex(tx, ty-0.05, tags[1].c_str());
    tex->DrawLatex(tx, ty-0.1, "Custom Function");
    
    c->cd(2);

    h2_rho_old_et->GetXaxis()->SetNdivisions(505);
    h2_rho_old_et->Draw("colz");
    tex->DrawLatex(tx, ty, tags[0].c_str());
    tex->DrawLatex(tx, ty-0.05, tags[1].c_str());
    tex->DrawLatex(tx, ty-0.1, "FastJet Function");

    c->SaveAs("rho_new.png");


    f_cent->Close();
    f_sum_et->Close();

    return;



}
void ProcessBackgroundTree(const std::string & input_file, bool x_axis_cent = true);
void ProcessGlobal(const std::string & input_file);
void RhoTest(const std::string & input_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/RHOTEST_DATA.root") 
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
    
    // ProcessGlobal(input_file);
    ProcessBackgroundTree(input_file, true); // x-axis = cent
    // ProcessBackgroundTree(input_file, false); // x-axis = SUM_Et
    // PlotRho();
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

        float tower_sum_energy_TOWERINFO_CALIB_CEMC = 0;
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC", &tower_sum_energy_TOWERINFO_CALIB_CEMC);
        float tower_sum_energy_TOWERINFO_CALIB_HCALIN = 0;
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN", &tower_sum_energy_TOWERINFO_CALIB_HCALIN);
        float tower_sum_energy_TOWERINFO_CALIB_HCALOUT = 0;
        t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT", &tower_sum_energy_TOWERINFO_CALIB_HCALOUT);



    int nentries = t->GetEntries();

    TH2F * h2_rho_area_vs_x;
    TH2F * h2_rho_area_alt_vs_x;
    TH2F * h2_rho_mult_vs_x;
    TH2F * h2_rho_area_times_A_vs_x;
    TH2F * h2_rho_area_alt_times_A_vs_x;

    if ( x_axis_sum_et ) {
        h2_rho_area_vs_x = new TH2F("h2_rho_area_vs_x", "h2_rho_area_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_RHO_A_BINS, RHO_A_BINS);
        h2_rho_area_alt_vs_x = new TH2F("h2_rho_area_alt_vs_x", "h2_rho_area_alt_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_RHO_A_BINS, RHO_A_BINS);
        h2_rho_mult_vs_x = new TH2F("h2_rho_mult_vs_x", "h2_rho_mult_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_RHO_M_BINS, RHO_M_BINS);
        h2_rho_area_times_A_vs_x = new TH2F("h2_rho_area_times_A_vs_x", "h2_rho_area_times_A_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_RHO_A_BINS, RHO_A_BINS);
        h2_rho_area_alt_times_A_vs_x = new TH2F("h2_rho_area_alt_times_A_vs_x", "h2_rho_area_alt_times_A_vs_x", N_SUM_ET_BINS, SUM_ET_BINS, N_RHO_A_BINS, RHO_A_BINS);
        std::vector<TH2F*> h2s = {h2_rho_area_vs_x, h2_rho_mult_vs_x, h2_rho_area_times_A_vs_x, h2_rho_area_alt_vs_x, h2_rho_area_alt_times_A_vs_x};
        for ( auto h2 : h2s ) {
            h2->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        }
    } else {
        h2_rho_area_vs_x = new TH2F("h2_rho_area_vs_x", "h2_rho_area_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_A_BINS, RHO_A_BINS);
        h2_rho_mult_vs_x = new TH2F("h2_rho_mult_vs_x", "h2_rho_mult_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_M_BINS, RHO_M_BINS);
        h2_rho_area_times_A_vs_x = new TH2F("h2_rho_area_times_A_vs_x", "h2_rho_area_times_A_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_A_BINS, RHO_A_BINS);
        h2_rho_area_alt_vs_x = new TH2F("h2_rho_area_alt_vs_x", "h2_rho_area_alt_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_A_BINS, RHO_A_BINS);
        h2_rho_area_alt_times_A_vs_x = new TH2F("h2_rho_area_alt_times_A_vs_x", "h2_rho_area_alt_times_A_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_RHO_A_BINS, RHO_A_BINS);

        std::vector<TH2F*> h2s = {h2_rho_area_vs_x, h2_rho_mult_vs_x, h2_rho_area_times_A_vs_x, h2_rho_area_alt_vs_x, h2_rho_area_alt_times_A_vs_x};
        for ( auto h2 : h2s ) {
            h2->GetXaxis()->SetTitle("Centrality [%]");
        }
    }

    h2_rho_area_vs_x->GetYaxis()->SetTitle("#rho_{A} [GeV]");
    h2_rho_mult_vs_x->GetYaxis()->SetTitle("#rho_{M} [GeV]");
    h2_rho_area_times_A_vs_x->GetYaxis()->SetTitle("#rho_{A} #times A [GeV]");
    h2_rho_area_alt_vs_x->GetYaxis()->SetTitle("#rho_{A} [GeV]");
    h2_rho_area_alt_times_A_vs_x->GetYaxis()->SetTitle("#rho_{A} #times A [GeV]");


    const int N_COURSE_X_BINS = 8;
    float COURSE_X_BINS[N_COURSE_X_BINS+1];
    float MAX_COURSE_X = 2000;
    if (x_axis_cent){ MAX_COURSE_X = MAX_X_CENT; }
    for ( int i = 0; i < N_COURSE_X_BINS+1; ++i ) { COURSE_X_BINS[i] = i*MAX_COURSE_X/N_COURSE_X_BINS; }
    TH2F * h2_rho_area_times_A_vs_course_x = new TH2F("h2_rho_area_times_A_vs_course_x", "h2_rho_area_times_A_vs_course_x", 
                    N_COURSE_X_BINS, COURSE_X_BINS, N_RHO_A_BINS, RHO_A_BINS);
    TH2F * h2_rho_area_alt_times_A_vs_course_x = new TH2F("h2_rho_area_alt_times_A_vs_course_x", "h2_rho_area_alt_times_A_vs_course_x", 
                    N_COURSE_X_BINS, COURSE_X_BINS, N_RHO_A_BINS, RHO_A_BINS);
    if (x_axis_sum_et) {
        h2_rho_area_times_A_vs_course_x->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
        h2_rho_area_alt_times_A_vs_course_x->GetXaxis()->SetTitle("#Sigma E_{T}^{Raw} [GeV]");
    } else {
        h2_rho_area_times_A_vs_course_x->GetXaxis()->SetTitle("Centrality [%]");
        h2_rho_area_alt_times_A_vs_course_x->GetXaxis()->SetTitle("Centrality [%]");
    }
    h2_rho_area_times_A_vs_course_x->GetYaxis()->SetTitle("#rho_{A} #times A [GeV]");
    h2_rho_area_alt_times_A_vs_course_x->GetYaxis()->SetTitle("#rho_{A} #times A [GeV]");
    
    for ( int i = 0; i < nentries; ++i ) {
        t->GetEntry(i);
        float sum_et = tower_sum_energy_TOWERINFO_CALIB_CEMC + tower_sum_energy_TOWERINFO_CALIB_HCALIN + tower_sum_energy_TOWERINFO_CALIB_HCALOUT;
        float xaxis_var = sum_et;
        if ( x_axis_cent ) { xaxis_var = 1.0*centrality; }
        h2_rho_area_vs_x->Fill(xaxis_var, rho_val_TowerRho_AREA);
        h2_rho_mult_vs_x->Fill(xaxis_var, rho_val_TowerRho_MULT);
        h2_rho_area_alt_vs_x->Fill(xaxis_var, std_rho_val_TowerRho_AREA);
        h2_rho_area_times_A_vs_x->Fill(xaxis_var, rho_val_TowerRho_AREA*AREA_CONE);
        h2_rho_area_alt_times_A_vs_x->Fill(xaxis_var, std_rho_val_TowerRho_AREA*AREA_CONE);
        h2_rho_area_times_A_vs_course_x->Fill(xaxis_var, rho_val_TowerRho_AREA*AREA_CONE);
        h2_rho_area_alt_times_A_vs_course_x->Fill(xaxis_var, std_rho_val_TowerRho_AREA*AREA_CONE);
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
    std::vector<TH2F*> h2s = {h2_rho_area_vs_x, h2_rho_area_times_A_vs_x, h2_rho_area_alt_vs_x, h2_rho_area_alt_times_A_vs_x};
    std::vector<std::string> labs = {"Custom Function", "Custom Function", "FastJet::JetMedianBackgroundEstimator", "FastJet::JetMedianBackgroundEstimator"};
    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};
    int ilab = 0;
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
        tex->DrawLatex(tx, ty, labs[ilab].c_str());
        ilab++;
        c->SaveAs((outdir+"/"+h2->GetTitle()+".png").c_str());
        delete c;
    }

    TLegend * leg = new TLegend(0.18,0.7,0.44,0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetNColumns(2);

    std::vector<TH2F*> h2s2 = {h2_rho_area_times_A_vs_course_x, h2_rho_area_alt_times_A_vs_course_x};
    double miny = 1e-4;
    labs = {"Custom", "FastJet"};
    ilab = 0;
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
            h1->GetYaxis()->SetTitle("Probability Density [a.u.]");
            
            
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
        tex->DrawLatex(txx, tyy, labs[ilab].c_str());
        ilab++;

        c->SaveAs((outdir+"/"+h2->GetTitle()+".png").c_str());
        delete c;
        leg->Clear();
    }
    
    // if ( x_axis_cent ) {
    //     h2s = {h2_rho_area_vs_x,  h2_rho_area_alt_vs_x};
    //     ilab = 0;
    //     for ( unsigned int ibin = 0; ibin  < N_X_CENT_BINS; ibin++){
    //         float maxx= 0;
    //         double miny = 1e-4;
    //         for ( auto h2 : h2s ) {
    //             tx = 0.19;
    //             double ty = ty_start;
    //             tags = {sPHENIX_Tag, DataType_Tag};
    //             labs = {"Custom", "FastJet"};
    //             c = new TCanvas("c", "c", 800, 600);
    //             gPad->SetLogy();
    //             gPad->SetLeftMargin(0.15);
    //             gPad->SetRightMargin(0.1);
    //             gPad->SetBottomMargin(0.15);
    //             gPad->SetTopMargin(0.05);
    //             h2->GetXaxis()->SetRange(ibin+1, ibin+1);
    //             TH1D * h1 = h2->ProjectionY(Form("h1_%s_%d", h2->GetTitle(), ibin), ibin+1, ibin+1);
    //             h1->Scale(1.0/h1->Integral());
    //             std::string leg_title = Form("%0.0f < #Sigma E_{T}^{Raw} < %0.0f GeV", COURSE_X_BINS[ibin], COURSE_X_BINS[ibin+1]);
    //             if ( x_axis_cent ) { leg_title = Form("%0.0f-%0.0f %%", COURSE_X_BINS[ibin], COURSE_X_BINS[ibin+1]); }
    //             tags.push_back(leg_title);
    //             h1->GetYaxis()->SetRangeUser(miny, 1e1);
    //             int lastbin_above_threshold = 0;
    //             lastbin_above_threshold = h1->FindLastBinAbove(miny);
    //             h1->GetXaxis()->SetRangeUser(0, h1->GetBinCenter(lastbin_above_threshold));
    //             h1->GetXaxis()->SetTitle(h2->GetYaxis()->GetTitle());
    //             h1->GetYaxis()->SetTitle("Probability Density [a.u.]");
    //             h1->Draw("P");
    //             tags.push_back(labs[ilab]);
    //             for ( auto tag : tags ) {
    //                 tex->DrawLatex(tx, ty, tag.c_str());
    //                 ty -= 0.05;
    //             }
    //             c->SaveAs((outdir+"/"+h2->GetTitle()+"_cent"+std::to_string(ibin)+".png").c_str());
    //             delete c;
    //         ilab++;
    //         }
            
    //     }

    // }

    TFile * fout = new TFile((outdir+"/BackgroundEstimates.root").c_str(), "RECREATE");
    fout->cd();
    h2_rho_area_vs_x->Write();
    h2_rho_mult_vs_x->Write();
    h2_rho_area_times_A_vs_x->Write();
    h2_rho_area_times_A_vs_course_x->Write();
    h2_rho_area_alt_vs_x->Write();
    h2_rho_area_alt_times_A_vs_x->Write();
    fout->Close();


    f->Close();

    return ;
    
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
