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

const int N_CONECOMP_BINS = 500;
float CONECOMP_MAX = 1100;
float CONECOMP_BINS[N_CONECOMP_BINS+1];

const int N_CONECOMP_SUB1_BINS = 100;
float CONECOMP_SUB1_MAX = 100;
float CONECOMP_SUB1_BINS[N_CONECOMP_SUB1_BINS+1];

const int N_CONE_DET_BINS = 250;
float MAX_CONE_DET = 60;
float CONE_DET_BINS[N_CONE_DET_BINS+1];

const int N_PROBE_DET_BINS = 250;
float MAX_PROBE_DET = 100;
float PROBE_DET_BINS[N_PROBE_DET_BINS+1];



void SetBins(){
    for ( int i = 0; i < N_SUM_ET_BINS+1; ++i ) { SUM_ET_BINS[i] = i*MAX_SUM_ET/N_SUM_ET_BINS; }
    for ( int i = 0; i < N_SUM_Q_BINS+1; ++i ) { SUM_Q_BINS[i] = i*MAX_SUM_Q/N_SUM_Q_BINS; }
    for ( int i = 0; i < N_CONECOMP_BINS+1; ++i ) { CONECOMP_BINS[i] = 500 + i*CONECOMP_MAX/N_CONECOMP_BINS; }
    for ( int i = 0; i < N_CONECOMP_SUB1_BINS+1; ++i ) { CONECOMP_SUB1_BINS[i] = i*CONECOMP_SUB1_MAX/N_CONECOMP_SUB1_BINS; }
    for ( int i = 0; i < N_CONE_DET_BINS+1; ++i ) { CONE_DET_BINS[i] = -MAX_CONE_DET + i*2*MAX_CONE_DET/N_CONE_DET_BINS; }
    for ( int i = 0; i < N_PROBE_DET_BINS+1; ++i ) { PROBE_DET_BINS[i] = -100 + i*2*MAX_PROBE_DET/N_PROBE_DET_BINS; }
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

std::string ProcessProbeTree(const std::string & input_file, const std::string & prefix);

void ProbePlot() 
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

    const std::string & input_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/FEB7/feb7_basic.root";
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

    
    ConfigureOutputDirs("Probes");
    std::cout << "Probe plots: " << probe_plots << std::endl;

    SetBins();
    
    std::string probe_hist = ProcessProbeTree(input_file, "probe");
    
    gSystem->Exit(0);
  
   
}

std::string ProcessProbeTree(const std::string & input_file, const std::string & prefix  )
{

    std::string outdir = probe_plots;
    outdir += "/" + prefix;
    outdir += "/xaxis_cent";
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

        float probe_jet_truth_energy = 0;
        float probe_jet_sub1_truth_energy = 0;
        float probe_jet_energy = 0;
        float probe_jet_sub1_energy = 0;
        float probe_jet_area = 0;
        float probe_jet_sub1_area = 0;
        int probe_jet_num_towers = 0;
        int probe_jet_sub1_num_towers = 0;
        t->SetBranchAddress("probe_jet_truth_energy", &probe_jet_truth_energy);
        t->SetBranchAddress("probe_jet_sub1_truth_energy", &probe_jet_sub1_truth_energy);
        t->SetBranchAddress("probe_jet_area", &probe_jet_area);
        t->SetBranchAddress("probe_jet_sub1_area", &probe_jet_sub1_area);
        t->SetBranchAddress("probe_jet_num_towers", &probe_jet_num_towers);
        t->SetBranchAddress("probe_jet_sub1_num_towers", &probe_jet_sub1_num_towers);
        t->SetBranchAddress("probe_jet_energy", &probe_jet_energy);
        t->SetBranchAddress("probe_jet_sub1_energy", &probe_jet_sub1_energy);

    int nentries = t->GetEntries();
    std::cout << "Processing " << nentries << " events" << std::endl;
  
  
    TH2F * h2_area_res_vs_x = new TH2F("h2_area_res_vs_x", "h2_area_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_PROBE_DET_BINS, PROBE_DET_BINS);
    TH2F * h2_mult_res_vs_x = new TH2F("h2_mult_res_vs_x", "h2_mult_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_PROBE_DET_BINS, PROBE_DET_BINS);
    TH2F * h2_sub1_res_vs_x = new TH2F("h2_sub1_res_vs_x", "h2_sub1_res_vs_x", N_X_CENT_BINS, X_CENT_BINS, N_PROBE_DET_BINS, PROBE_DET_BINS);
    TH1F * h1_probe_et = new TH1F("h1_probe_et", "h1_probe_et", N_PROBE_DET_BINS, PROBE_DET_BINS);
    TH1F * h1_rhoA = new TH1F("h1_rhoA", "h1_rhoA", N_PROBE_DET_BINS, PROBE_DET_BINS);
    TH1F * h1_rhoM = new TH1F("h1_rhoM", "h1_rhoM", N_PROBE_DET_BINS, PROBE_DET_BINS);
    TH1F * h1_probe_sub1_et = new TH1F("h1_probe_sub1_et", "h1_probe_sub1_et", N_PROBE_DET_BINS, PROBE_DET_BINS);
    TH1F * h1_probe_truth_et = new TH1F("h1_probe_truth_et", "h1_probe_truth_et", N_PROBE_DET_BINS, PROBE_DET_BINS);
    TH1F * h1_probe_truth_sub1_et = new TH1F("h1_probe_truth_sub1_et", "h1_probe_truth_sub1_et", N_PROBE_DET_BINS, PROBE_DET_BINS);
   
   
    std::vector<TH2F*> h2s = {h2_area_res_vs_x, h2_mult_res_vs_x, h2_sub1_res_vs_x};
    for ( auto h2 : h2s ) {
        h2->GetXaxis()->SetTitle("Centrality [%]");
        h2->GetYaxis()->SetTitle("#delta E_{T}^{Probe} [GeV]");
    }
    std::vector<std::string> labs = {"Area", "Multiplicity", "Iterative"};

    for ( int i = 0; i < nentries; ++i ) {
       
        t->GetEntry(i);
        float xaxis_var =1.0*centrality;
        if ( probe_jet_energy == 0 || probe_jet_sub1_energy == 0 || probe_jet_truth_energy == 0 ) { continue; }

        float area_bkgd = rho_val_TowerRho_AREA*probe_jet_area - probe_jet_truth_energy;
        float mult_bkgd = rho_val_TowerRho_MULT*probe_jet_num_towers -probe_jet_truth_energy;
        float sub1_bkdg = -1.0*(probe_jet_sub1_truth_energy);
        float cone_res_area = probe_jet_energy - rho_val_TowerRho_AREA*probe_jet_area - probe_jet_truth_energy;
        float cone_res_mult = probe_jet_energy - rho_val_TowerRho_MULT*probe_jet_num_towers - probe_jet_truth_energy;
        float cone_res_sub1 = probe_jet_sub1_energy - probe_jet_sub1_truth_energy;
        

        h2_area_res_vs_x->Fill(xaxis_var, cone_res_area);
        h2_mult_res_vs_x->Fill(xaxis_var, cone_res_mult);
        h2_sub1_res_vs_x->Fill(xaxis_var, cone_res_sub1);
        h1_probe_et->Fill(probe_jet_energy);
        h1_rhoA->Fill(rho_val_TowerRho_AREA*probe_jet_area);
        h1_rhoM->Fill(rho_val_TowerRho_MULT*probe_jet_num_towers);
        h1_probe_sub1_et->Fill(probe_jet_sub1_energy);
        h1_probe_truth_et->Fill(probe_jet_truth_energy);
        h1_probe_truth_sub1_et->Fill(probe_jet_sub1_truth_energy);



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
    c->SaveAs((outdir+"/probe_res_cent_slices.png").c_str());

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
            // std::string leg_title = Form("%0.0f-%0.0f %%", h2->GetXaxis()->GetBinLowEdge(ibin+1), h2->GetXaxis()->GetBinUpEdge(ibin+1)); 
            // tags.push_back(leg_title);
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
            h1->GetXaxis()->SetRangeUser(-30,30);


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

        c->SaveAs((outdir+"/probe_res_vs_x_slice_"+std::to_string(ibin)+".png").c_str());
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
    c->SaveAs((outdir+"/cone_res_vs_x.png").c_str());
    delete c;
    TFile * fout = new TFile((outdir+"/probes.root").c_str(), "RECREATE");
    for ( auto h2 : h2s ) {
        h2->Write();
    }
    for ( unsigned ihist = 0; ihist < h2s.size(); ihist++ ) {
        g_std_devs[ihist]->SetName(labs[ihist].c_str());
        g_std_devs[ihist]->Write();
    }
    for ( auto h1 : {h1_probe_et, h1_rhoA, h1_rhoM, h1_probe_sub1_et, h1_probe_truth_et, h1_probe_truth_sub1_et} ) {
        h1->Write();
    }
    fout->Close();
    f->Close();

    

    return outdir+"/probes.root";
    
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



