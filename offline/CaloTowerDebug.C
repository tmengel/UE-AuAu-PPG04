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

TH1::SetDefaultSumw2();
TH2::SetDefaultSumw2();
TH3::SetDefaultSumw2();

SetsPhenixStyle();

gErrorIgnoreLevel = kWarning;
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetOptTitle(0);
gStyle->SetPalette(kRainBow);

const int COLORS[] = {kBlack, kRed , kBlue, kGreen+2, kViolet, kCyan, kOrange+2, kMagenta+2, kAzure-2};
const int MARKERS[] = { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kOpenCircle, kOpenSquare, kOpenTriangleUp};
const float MARKER_SIZE = 1.5;
const float LINE_WIDTH = 2.0;

const std::string DATA_DIR = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/";
const std::string sPHENIX_Tag = "#it{#bf{sPHENIX}} Internal";
const std::string DataType_Tag = "2024 Au+Au 200 GeV";

float MAX_CENT = 80;
const float CENT_BINS[]= {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
const int N_CENT_BINS = sizeof(CENT_BINS)/sizeof(CENT_BINS[0]) - 1;
const float COURSE_CENT_BINS[]= {0, 10, 20, 30, 40, 50, 60, 70, 80};
const int N_COURSE_CENT_BINS = sizeof(COURSE_CENT_BINS)/sizeof(COURSE_CENT_BINS[0]) - 1;

const int N_ET_BINS = 210; // et or Q
float ET_BINS[N_ET_BINS+1];
float MAX_ET = 200;


const int N_CONECOMP_BINS = 500;
float CONECOMP_MAX = 1100;
float CONECOMP_BINS[N_CONECOMP_BINS+1];

const int N_CONECOMP_SUB1_BINS = 100;
float CONECOMP_SUB1_MAX = 100;
float CONECOMP_SUB1_BINS[N_CONECOMP_SUB1_BINS+1];



void SetBins(){
    for ( int i = 0; i < N_ET_BINS+1; ++i ) { ET_BINS[i] = -10 +  i*MAX_ET/N_ET_BINS; }
    for ( int i = 0; i < N_CONECOMP_BINS+1; ++i ) { CONECOMP_BINS[i] = 500 + i*CONECOMP_MAX/N_CONECOMP_BINS; }
    for ( int i = 0; i < N_CONECOMP_SUB1_BINS+1; ++i ) { CONECOMP_SUB1_BINS[i] = i*CONECOMP_SUB1_MAX/N_CONECOMP_SUB1_BINS; }
}


const float N_CEMC_TOWERS = 256*96;
const float N_HCAL_TOWERS = 64*24;
const float AREA_CONE = TMath::Pi()*0.4*0.4;
const float AREA_CEMC_TOWER = (2.0*TMath::Pi()/256.0)*(2.2/96.0);
const float AREA_HCAL_TOWER = (2.0*TMath::Pi()/64.0)*(2.2/24.0);

int NEVENTS = 0;
void MakeContPlots( const std::string & input_file,  const std::string & output_file) {
    
    TFile * f = new TFile(Form("%s%s", DATA_DIR.c_str(), input_file.c_str()));
    TTree * t = (TTree*)f->Get("T");
    if (!t) { std::cout << "Tree T not found" << std::endl; exit(1); }

    int cent;
    t->SetBranchAddress("centrality", &cent);
    float rc_energy, rc_energy_sub1;
    float rc_eCemc, rc_eIhcal, rc_eOhcal;
    float rc_eCemc_sub1, rc_eIhcal_sub1, rc_eOhcal_sub1;
    t->SetBranchAddress("random_cone_energy_RandomCones_r04", &rc_energy);
    t->SetBranchAddress("random_cone_energy_RandomCones_r04_Sub1", &rc_energy_sub1);
    t->SetBranchAddress("random_cone_energy_cemc_RandomCones_r04", &rc_eCemc);
    t->SetBranchAddress("random_cone_energy_hcalin_RandomCones_r04", &rc_eIhcal);
    t->SetBranchAddress("random_cone_energy_hcalout_RandomCones_r04", &rc_eOhcal);
    t->SetBranchAddress("random_cone_energy_cemc_RandomCones_r04_Sub1", &rc_eCemc_sub1);
    t->SetBranchAddress("random_cone_energy_hcalin_RandomCones_r04_Sub1", &rc_eIhcal_sub1);
    t->SetBranchAddress("random_cone_energy_hcalout_RandomCones_r04_Sub1", &rc_eOhcal_sub1);
    TH2F * h2_rc_et = new TH2F("h2_rc_et", ";Centrality;E_{T}^{RC} (GeV)", N_CENT_BINS, CENT_BINS, N_ET_BINS, ET_BINS);
    TH2F * h2_rc_et_cemc = new TH2F("h2_rc_et_cemc", ";Centrality;E_{T}^{RC} (GeV)", N_CENT_BINS, CENT_BINS, N_ET_BINS, ET_BINS);
    TH2F * h2_rc_et_hcalin = new TH2F("h2_rc_et_hcalin", ";Centrality;E_{T}^{RC} (GeV)", N_CENT_BINS, CENT_BINS, N_ET_BINS, ET_BINS);
    TH2F * h2_rc_et_hcalout = new TH2F("h2_rc_et_hcalout", ";Centrality;E_{T}^{RC} (GeV)", N_CENT_BINS, CENT_BINS, N_ET_BINS, ET_BINS);
    TH2F * h2_rc_et_sub1 = new TH2F("h2_rc_et_sub1", ";Centrality;E_{T}^{RC} (GeV)", N_CENT_BINS, CENT_BINS, N_ET_BINS, ET_BINS);
    TH2F * h2_rc_et_cemc_sub1 = new TH2F("h2_rc_et_cemc_sub1", ";Centrality;E_{T}^{RC} (GeV)", N_CENT_BINS, CENT_BINS, N_ET_BINS, ET_BINS);
    TH2F * h2_rc_et_hcalout_sub1 = new TH2F("h2_rc_et_hcalout_sub1", ";Centrality;E_{T}^{RC} (GeV)", N_CENT_BINS, CENT_BINS, N_ET_BINS, ET_BINS);
    TH2F * h2_rc_et_hcalin_sub1 = new TH2F("h2_rc_et_hcalin_sub1", ";Centrality;E_{T}^{RC} (GeV)", N_CENT_BINS, CENT_BINS, N_ET_BINS, ET_BINS);


    int rc_nTowers, rc_nTowers_sub1;
    int rc_nCemc, rc_nIhcal, rc_nOhcal;
    int rc_nCemc_sub1, rc_nIhcal_sub1, rc_nOhcal_sub1;
    t->SetBranchAddress("random_cone_num_towers_RandomCones_r04", &rc_nTowers);
    t->SetBranchAddress("random_cone_num_towers_RandomCones_r04_Sub1", &rc_nTowers_sub1);
    t->SetBranchAddress("random_cone_num_towers_cemc_RandomCones_r04", &rc_nCemc);
    t->SetBranchAddress("random_cone_num_towers_hcalin_RandomCones_r04", &rc_nIhcal);
    t->SetBranchAddress("random_cone_num_towers_hcalout_RandomCones_r04", &rc_nOhcal);
    t->SetBranchAddress("random_cone_num_towers_cemc_RandomCones_r04_Sub1", &rc_nCemc_sub1);
    t->SetBranchAddress("random_cone_num_towers_hcalin_RandomCones_r04_Sub1", &rc_nIhcal_sub1);
    t->SetBranchAddress("random_cone_num_towers_hcalout_RandomCones_r04_Sub1", &rc_nOhcal_sub1);
    TH2F * h2_rc_ncomp = new TH2F("h2_rc_ncomp", ";Centrality;N_{Towers}^{RC}", N_CENT_BINS, CENT_BINS, N_CONECOMP_BINS, CONECOMP_BINS);
    TH2F * h2_rc_ncomp_sub1 = new TH2F("h2_rc_ncomp_sub1", ";Centrality;N_{Towers}^{RC}", N_CENT_BINS, CENT_BINS, N_CONECOMP_SUB1_BINS, CONECOMP_SUB1_BINS);

    int rc_nMaskedTowers, rc_nMaskedTowers_sub1;
    // int rc_nMaskedCemc, rc_nMaskedIhcal, rc_nMaskedOhcal;
    // int rc_nMaskedCemc_sub1, rc_nMaskedIhcal_sub1, rc_nMaskedOhcal_sub1;
    t->SetBranchAddress("random_cone_num_masked_towers_RandomCones_r04", &rc_nMaskedTowers);
    // t->SetBranchAddress("random_cone_num_masked_towers_cemc_RandomCones_r04", &rc_nMaskedCemc);
    // t->SetBranchAddress("random_cone_num_masked_towers_hcalin_RandomCones_r04", &rc_nMaskedIhcal);
    // t->SetBranchAddress("random_cone_num_masked_towers_hcalout_RandomCones_r04", &rc_nMaskedOhcal);
    t->SetBranchAddress("random_cone_num_masked_towers_RandomCones_r04_Sub1", &rc_nMaskedTowers_sub1);  
    // t->SetBranchAddress("random_cone_num_masked_towers_cemc_RandomCones_r04_Sub1", &rc_nMaskedCemc_sub1);
    // t->SetBranchAddress("random_cone_num_masked_towers_hcalin_RandomCones_r04_Sub1", &rc_nMaskedIhcal_sub1);
    // t->SetBranchAddress("random_cone_num_masked_towers_hcalout_RandomCones_r04_Sub1", &rc_nMaskedOhcal_sub1);
    TH2F * h2_rc_nmask = new TH2F("h2_rc_nmask", ";Centrality;N_{Masked Towers}^{RC}", N_CENT_BINS, CENT_BINS, N_CONECOMP_SUB1_BINS, CONECOMP_SUB1_BINS);
    TH2F * h2_rc_nmask_sub1 = new TH2F("h2_rc_nmask_sub1", ";Centrality;N_{Masked Towers}^{RC}", N_CENT_BINS, CENT_BINS, N_CONECOMP_SUB1_BINS, CONECOMP_SUB1_BINS);

    int nentries = t->GetEntries();
    for (int i = 0; i < nentries; ++i) {
        t->GetEntry(i);
        if (cent > MAX_CENT){ continue; }
        h2_rc_et->Fill(cent, rc_energy);
        h2_rc_et_sub1->Fill(cent, rc_energy_sub1);
        h2_rc_et_cemc->Fill(cent, rc_eCemc);
        h2_rc_et_hcalin->Fill(cent, rc_eIhcal);
        h2_rc_et_hcalout->Fill(cent, rc_eOhcal);
        h2_rc_et_cemc_sub1->Fill(cent, rc_eCemc_sub1);
        h2_rc_et_hcalin_sub1->Fill(cent, rc_eIhcal_sub1);
        h2_rc_et_hcalout_sub1->Fill(cent, rc_eOhcal_sub1);
        h2_rc_ncomp->Fill(cent, rc_nTowers);
        h2_rc_ncomp_sub1->Fill(cent, rc_nTowers_sub1);
        h2_rc_nmask->Fill(cent, rc_nMaskedTowers);
        h2_rc_nmask_sub1->Fill(cent, rc_nMaskedTowers_sub1);

    }

    TFile * foutput = new TFile(output_file.c_str(), "RECREATE");
    foutput->cd();
    h2_rc_et->Write();
    h2_rc_et_sub1->Write();
    h2_rc_et_cemc->Write();
    h2_rc_et_hcalin->Write();
    h2_rc_et_hcalout->Write();
    h2_rc_et_cemc_sub1->Write();
    h2_rc_et_hcalin_sub1->Write();
    h2_rc_et_hcalout_sub1->Write();
    h2_rc_ncomp->Write();
    h2_rc_ncomp_sub1->Write();
    h2_rc_nmask->Write();
    h2_rc_nmask_sub1->Write();
    foutput->Close();
    f->Close();
    return;



}
void CompPlots( const std::string & input_file_basic, const std::string & input_file_random ){
    TFile * fbasic = new TFile(input_file_basic.c_str());
    // TFile * fbasic = new TFile("tower_debug_basic.root");
    if(!fbasic->IsOpen() || fbasic->IsZombie()){ std::cout << "File " << input_file_basic << " is zombie" << std::endl;  exit(1); }
    TFile * frandom = new TFile(input_file_random.c_str());
    // TFile * frandom = new TFile("tower_debug_random.root");
    if(!frandom->IsOpen() || frandom->IsZombie()){ std::cout << "File " << input_file_random << " is zombie" << std::endl;  exit(1); }

    std::string output_path = "debugtowers/";

    TH2F * h2_cone_et = (TH2F*)fbasic->Get("h2_rc_et");
    h2_cone_et->SetName("h2_rc_et_basic");
    TH2F * h2_cone_et_r = (TH2F*)frandom->Get("h2_rc_et");

    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    TCanvas * c;
    TLegend * leg = new TLegend(0.55, 0.55, 0.9, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    // std::vector<std::string> h2_labs = {"CEMC", "iHCAL", "oHCAL", "CEMC Retower", "CEMC Sub1", "iHCAL Sub1", "oHCAL Sub1"};
    // std::vector<std::string> h2_labs = {"Sum", "Sum Sub1", "Sum Retower"};
    for (int i = 0; i < N_CENT_BINS; i++) {
        leg->Clear();

        h2_cone_et->GetXaxis()->SetRangeUser(CENT_BINS[i], CENT_BINS[i+1]);
        TH1F * h1 = (TH1F*) h2_cone_et->ProjectionY(Form("h1_%s_cent%d", h2_cone_et->GetName(), i));
        h1->GetXaxis()->SetTitle("E_{T} [GeV]");
        h1->GetYaxis()->SetTitle("P(E_{T})");
        h2_cone_et_r->GetXaxis()->SetRangeUser(CENT_BINS[i], CENT_BINS[i+1]);
        TH1F * h1_r = (TH1F*) h2_cone_et_r->ProjectionY(Form("h1_%s_cent%d", h2_cone_et_r->GetName(), i));
        h1_r->SetLineColor(kRed);

        // make courser bin
        // h1->Rebin(5);
        // h1_r->Rebin(5);
        h1->Scale(1.0/h1->Integral(), "width");
        h1_r->Scale(1.0/h1_r->Integral(), "width");
        
        
        c = new TCanvas("c", "c", 1200, 600);
        c->Divide(2,1);

        c->cd(2);
        TH1F * h1_ratio = (TH1F*)h1->Clone(Form("h1_ratio_%s_cent%d", h2_cone_et->GetName(), i));
        h1_ratio->Divide(h1_r);
        h1_ratio->GetYaxis()->SetTitle("P(E_{T}^{Basic}) / P(E_{T}^{Random})");
        h1_ratio->Draw("HIST");
        
        c->cd(1);
        gPad->SetLogy();    
        h1->Draw("HIST");
        leg->AddEntry(h1, "Basic", "l");
        h1_r->Draw("same HIST");
        leg->AddEntry(h1_r, "Random", "l");
        tex->DrawLatex(0.55, 0.85, Form("RCone"));
        tex->DrawLatex(0.55, 0.8, Form("%d-%d%%", int(CENT_BINS[i]), int(CENT_BINS[i+1])));
        leg->Draw("same");

        std::string savename = "rccone_comp_cent" + std::to_string(i);
        savename.erase(remove(savename.begin(), savename.end(), ' '), savename.end());
        c->SaveAs(Form("%s/%s.png", output_path.c_str(), savename.c_str()));
        delete c;
    }

    fbasic->Close();
    frandom->Close();
    return;

}

void CaloTowerDebug() {
    SetBins();

    const std::string & input_file_basic = "FEB7/feb7_basic.root";
    const std::string & input_file_random = "FEB10/feb10_random.root";

    const std::string & output_path = "debugtowers/";
    if (gSystem->AccessPathName(output_path.c_str())) {
        gSystem->mkdir(output_path.c_str(), true);
    }
    // MakeContPlots(input_file_basic, "debugtowers/cone_basic.root");
    // MakeContPlots(input_file_random, "debugtowers/cone_random.root");
    CompPlots("debugtowers/cone_basic.root", "debugtowers/cone_random.root");
   
    // TFile * fbasic = new TFile(Form("%s%s", DATA_DIR.c_str(), input_file_basic.c_str()));
    // // // TFile * fbasic = new TFile("tower_debug_basic.root");
    // if(!fbasic->IsOpen() || fbasic->IsZombie()){ std::cout << "File " << input_file_basic << " is zombie" << std::endl;  exit(1); }
    // TFile * frandom = new TFile(Form("%s%s", DATA_DIR.c_str(), input_file_random.c_str()));
    // // // TFile * frandom = new TFile("tower_debug_random.root");
    // if(!frandom->IsOpen() || frandom->IsZombie()){ std::cout << "File " << input_file_random << " is zombie" << std::endl;  exit(1); }

   



    


    // std::vector<std::string> h2_names = {
    //     "h2_tower_energy_cent_TOWERINFO_CALIB_CEMC",
    //     "h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN",
    //     "h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT",
    //     "h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER",
    //     "h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1",
    //     "h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1",
    //     "h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1"
    // };

 
    // std::vector<TH2F*> h2s_basic;
    // std::vector<TH2F*> h2s_random;
    // std::vector<std::string> h2_labs = {"CEMC", "iHCAL", "oHCAL", "CEMC Retower", "CEMC Sub1", "iHCAL Sub1", "oHCAL Sub1"};
    // std::vector<TH2F*> h2_sums_basic;
    // std::vector<TH2F*> h2_sums_random;


    // TLatex * tex = new TLatex();
    // tex->SetNDC();
    // tex->SetTextFont(42);
    // TCanvas * c;
    // TLegend * leg = new TLegend(0.55, 0.55, 0.9, 0.7);
    // leg->SetBorderSize(0);
    // leg->SetFillStyle(0);
    // int ihist = 0;
    // for (auto name : h2_names) {
    //     std::cout << "Processing " << name << std::endl;
    //     TH2F * h2 = (TH2F*)fbasic->Get(Form("%s", name.c_str()));
    //     if (!h2) { std::cout << "Basic Histogram " << name << " not found" << std::endl; exit(1); }
    //     h2->SetName(Form("%s_basic", h2->GetName()));
    //     h2s_basic.push_back(h2);
    //     TH2F * h2_r = (TH2F*)frandom->Get(Form("%s", name.c_str()));
    //     if (!h2_r) { std::cout << "Basic Histogram " << name << " not found" << std::endl; exit(1); }
    //     h2_r->SetName(Form("%s_random", h2_r->GetName()));
    //     h2s_random.push_back(h2_r);
    //     int ntowers_fired = h2->Integral();
    //     int ntowers_fired_r = h2_r->Integral();
    //     std::cout << "ntowers_fired = " << ntowers_fired << ", ntowers_fired_r = " << ntowers_fired_r << std::endl;


    //     for (int i = 0; i < N_CENT_BINS; i++) {
    //         leg->Clear();

    //         h2->GetYaxis()->SetRangeUser(CENT_BINS[i], CENT_BINS[i+1]);
    //         TH1F * h1 = (TH1F*) h2->ProjectionX(Form("h1_%s_cent%d", h2->GetName(), i));
    //         h1->GetXaxis()->SetTitle("E_{T} [GeV]");
    //         h1->GetYaxis()->SetTitle("P(E_{T})");
    //         h2_r->GetYaxis()->SetRangeUser(CENT_BINS[i], CENT_BINS[i+1]);
    //         TH1F * h1_r = (TH1F*) h2_r->ProjectionX(Form("h1_%s_cent%d", h2_r->GetName(), i));
    //         h1_r->SetLineColor(kRed);

    //         // make courser bin
    //         h1->Rebin(5);
    //         h1_r->Rebin(5);
    //         h1->Scale(1.0/h1->Integral(), "width");
    //         h1_r->Scale(1.0/h1_r->Integral(), "width");
            
           
    //         c = new TCanvas("c", "c", 1200, 600);
    //         c->Divide(2,1);

    //         c->cd(2);
    //         TH1F * h1_ratio = (TH1F*)h1->Clone(Form("h1_ratio_%s_cent%d", h2->GetName(), i));
    //         h1_ratio->Divide(h1_r);
    //         h1_ratio->GetYaxis()->SetTitle("P(E_{T}^{Basic}) / P(E_{T}^{Random})");
    //         h1_ratio->Draw("HIST");
         
    //         c->cd(1);
    //         gPad->SetLogy();    
    //         h1->Draw("HIST");
    //         leg->AddEntry(h1, "Basic", "l");
    //         h1_r->Draw("same HIST");
    //         leg->AddEntry(h1_r, "Random", "l");
    //         tex->DrawLatex(0.55, 0.85, Form("%s", h2_labs[ihist].c_str()));
    //         tex->DrawLatex(0.55, 0.8, Form("%d-%d%%", int(CENT_BINS[i]), int(CENT_BINS[i+1])));
    //         leg->Draw("same");

    //         std::string savename = h2_labs[ihist] + "_cent" + std::to_string(i);
    //         savename.erase(remove(savename.begin(), savename.end(), ' '), savename.end());
    //         c->SaveAs(Form("%s/%s.png", output_path.c_str(), savename.c_str()));
    //         delete c;
    //     }
    //     ihist++;
    // }
    // fbasic->Close();
    // frandom->Close();
    gSystem->Exit(0);
    // TCanvas * c;
    // TFile * foutput = new TFile("tower_debug.root");
    // TFile * foutput2 = new TFile("tower_debug_random.root");
    // if(!foutput->IsOpen() || foutput->IsZombie()){ std::cout << "File " << "tower_debug.root" << " is zombie" << std::endl;  exit(1); }
    // if(!foutput2->IsOpen() || foutput2->IsZombie()){ std::cout << "File " << "tower_debug_random.root" << " is zombie" << std::endl;  exit(1); }
    


    // TH2F * h2_cemc = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC");
    // TH2F * h2_ihcal = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN");
    // TH2F * h2_ohcal = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT");
    // TH2F * h2_cemc_retower = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER");
    // TH2F * h2_cemc_sub1 = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    // TH2F * h2_ihcal_sub1 = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1");
    // TH2F * h2_ohcal_sub1 = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1");
    // for (auto h2 : {h2_cemc, h2_ihcal, h2_ohcal, h2_cemc_retower, h2_cemc_sub1, h2_ihcal_sub1, h2_ohcal_sub1}) {
    //     if (!h2) { std::cout << "Histogram not found" << std::endl; exit(1); }
    //     h2->SetName(Form("%s_basic", h2->GetName()));
    // }
    // int ncemc = h2_cemc->Integral();
    // int nihcal = h2_ihcal->Integral();
    // int nohcal = h2_ohcal->Integral();
    // std ::cout << "ncemc = " << ncemc << ", nihcal = " << nihcal << ", nohcal = " << nohcal << std::endl;
    // // h2_cemc

    // TFile * foutput = new TFile("tower_debug_random.root", "RECREATE");
    // TH2F* h2_cemc_copy = (TH2F*)h2_cemc->Clone("h2_sum_random");
    // h2_cemc_copy->Add(h2_ihcal);
    // h2_cemc_copy->Add(h2_ohcal);
    
    // TH2F* h2_cemc_sub1_copy = (TH2F*)h2_cemc_sub1->Clone("h2_sum_sub1_random");
    // h2_cemc_sub1_copy->Add(h2_ihcal_sub1);
    // h2_cemc_sub1_copy->Add(h2_ohcal_sub1);

    // TH2F* h2_cemc_retower_copy = (TH2F*)h2_cemc_retower->Clone("h2_sum_retower_random");
    // h2_cemc_retower_copy->Add(h2_ihcal);
    // h2_cemc_retower_copy->Add(h2_ohcal);
    // foutput->cd();
    // h2_cemc_copy->Write();
    // h2_cemc_sub1_copy->Write();
    // h2_cemc_retower_copy->Write();
    // foutput->Close();




    // TH2F * h2_sum = (TH2F*)h2_cemc->Clone("h2_sum");
    // h2_sum->Add(h2_ihcal);
    // h2_sum->Add(h2_ohcal);

    // TH2F * h2_sum_sub1 = (TH2F*)h2_cemc_sub1->Clone("h2_sum_sub1");
    // h2_sum_sub1->Add(h2_ihcal_sub1);
    // h2_sum_sub1->Add(h2_ohcal_sub1);

    // TH2F * h2_sum_retower = (TH2F*)h2_cemc_retower->Clone("h2_sum_retower");
    // h2_sum_retower->Add(h2_ihcal);
    // h2_sum_retower->Add(h2_ohcal);

    // TH2F * h2_cemc_r = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC");
    // TH2F * h2_ihcal_r = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN");
    // TH2F * h2_ohcal_r = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT");
    // TH2F * h2_cemc_retower_r = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER");
    // TH2F * h2_cemc_sub1_r = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    // TH2F * h2_ihcal_sub1_r = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1");
    // TH2F * h2_ohcal_sub1_r = (TH2F*)frandom->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1");
    // for (auto h2 : {h2_cemc_r, h2_ihcal_r, h2_ohcal_r, h2_cemc_retower_r, h2_cemc_sub1_r, h2_ihcal_sub1_r, h2_ohcal_sub1_r}) {
    //     if (!h2) { std::cout << "Histogram not found" << std::endl; exit(1); }
    //     h2->SetName(Form("%s_random", h2->GetName()));
    // }

    // TLatex * tex = new TLatex();
    // tex->SetNDC();
    // tex->SetTextFont(42);

    // TLegend * leg = new TLegend(0.55, 0.55, 0.9, 0.7);
    // leg->SetBorderSize(0);
    // leg->SetFillStyle(0);

    // TH2F * h2_sum_r = (TH2F*)h2_cemc_r->Clone("h2_sum");
    // h2_sum_r->Add(h2_ihcal_r);
    // h2_sum_r->Add(h2_ohcal_r);

    // TH2F * h2_sum_sub1_r = (TH2F*)h2_cemc_sub1_r->Clone("h2_sum_sub1");
    // h2_sum_sub1_r->Add(h2_ihcal_sub1_r);
    // h2_sum_sub1_r->Add(h2_ohcal_sub1_r);

    // TH2F * h2_sum_retower_r = (TH2F*)h2_cemc_retower_r->Clone("h2_sum_retower");
    // h2_sum_retower_r->Add(h2_ihcal_r);
    // h2_sum_retower_r->Add(h2_ohcal_r);

    // std::vector<TH2F*> h2_sums = {h2_sum, h2_sum_sub1, h2_sum_retower};
    // std::vector<TH2F*> h2_sums_r = {h2_sum_r, h2_sum_sub1_r, h2_sum_retower_r};
    // std::vector<std::string> h2_sum_labs = {"Sum", "Sum Sub1", "Sum Retower"};
    // // ihist = 0;
    // for (unsigned int ihist = 0; ihist < h2_sums.size(); ihist++) {
    //     TH2F * h2 = h2_sums[ihist];
    //     int ntowers_fired = h2->Integral();
    //     TH2F * h2_r = h2_sums_r[ihist];
    //     int ntowers_fired_r = h2_r->Integral();
    //     std::cout << "ntowers_fired = " << ntowers_fired << ", ntowers_fired_r = " << ntowers_fired_r << std::endl;

    //     for (int i = 0; i < N_CENT_BINS; i++) {
    //         leg->Clear();

    //         h2->GetYaxis()->SetRangeUser(CENT_BINS[i], CENT_BINS[i+1]);
    //         TH1F * h1 = (TH1F*) h2->ProjectionX(Form("h1_%s_cent%d", h2->GetName(), i));
    //         h1->GetXaxis()->SetTitle("E_{T} [GeV]");
    //         h1->GetYaxis()->SetTitle("P(E_{T})");
    //         h1->Scale(1.0/h1->Integral(), "width");
    //         h1->SetLineColor(kBlack);
           
    //         h2_r->GetYaxis()->SetRangeUser(CENT_BINS[i], CENT_BINS[i+1]);
    //         TH1F * h1_r = (TH1F*) h2_r->ProjectionX(Form("h1_%s_cent%d", h2_r->GetName(), i));
    //         h1_r->SetLineColor(kRed); 
    //         h1_r->Scale(1.0/h1_r->Integral(), "width");

                
    //         c = new TCanvas("c", "c", 1200, 600);
    //         c->Divide(2,1);
    //         // make courser bin
    //         // h1->Rebin(3);
    //         // h1_r->Rebin(3);
    //         c->cd(1);
    //         gPad->SetLogy();    
    //         h1->Draw("HIST");
    //         leg->AddEntry(h1, "Basic", "l");
    //         h1_r->Draw("same HIST");
    //         leg->AddEntry(h1_r, "Random", "l");
    //         tex->DrawLatex(0.55, 0.85, Form("%s", h2_sum_labs[ihist].c_str()));
    //         tex->DrawLatex(0.55, 0.8, Form("%d-%d%%", int(CENT_BINS[i]), int(CENT_BINS[i+1])));
    //         leg->Draw("same");
           
            
       

    //         c->cd(2);
    //         TH1F * h1_ratio = (TH1F*)h1->Clone(Form("h1_ratio_%s_cent%d", h2->GetName(), i));
    //         h1_ratio->Divide(h1_r);
    //         h1_ratio->GetYaxis()->SetTitle("P(E_{T}^{Basic}) / P(E_{T}^{Random})");
    //         h1_ratio->Draw("HIST");
         

    //         std::string savename = h2_sum_labs[ihist] + "_cent" + std::to_string(i);
    //         savename.erase(remove(savename.begin(), savename.end(), ' '), savename.end());
    //         c->SaveAs(Form("%s/%s.png", output_path.c_str(), savename.c_str()));
    //         delete c;
    //     }
    // }
    // std::cout << "Done" << std::endl;
    
    // fbasic->Close();
    // frandom->Close();

    // delete tex;
    // delete leg;
    // gSystem->Exit(0);

    // TH2F * h2_cemc = (TH2F*)fbasic->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC");
    // TH2F * h2_ihcal = (TH2F*)fbasic->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN");
    // TH2F * h2_ohcal = (TH2F*)fbasic->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT");
    // TH2F * h2_cemc_retower = (TH2F*)fbasic->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER");
    // TH2F * h2_cemc_sub1 = (TH2F*)fbasic->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    // TH2F * h2_ihcal_sub1 = (TH2F*)fbasic->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1");
    // TH2F * h2_ohcal_sub1 = (TH2F*)fbasic->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1");
    // if (!h2_cemc || !h2_ihcal || !h2_ohcal || !h2_cemc_retower || !h2_cemc_sub1 || !h2_ihcal_sub1 || !h2_ohcal_sub1) {
    //     std::cout << "Histograms not found" << std::endl;
    //     exit(1);
    // }


    // std::vector<TH2F*> h2s = {h2_cemc, h2_ihcal, h2_ohcal, h2_cemc_retower, h2_cemc_sub1, h2_ihcal_sub1, h2_ohcal_sub1};
    // for (auto h2 : h2s) {
    //     h2->SetName(Form("%s_basic", h2->GetName()));
    // }
    // std::vector<std::string> h2_names = {"CEMC", "iHCAL", "oHCAL", "CEMC Retower", "CEMC Sub1", "iHCAL Sub1", "oHCAL Sub1"};
    
    // TLatex * tex = new TLatex();
    // tex->SetNDC();
    // tex->SetTextFont(42);
    // TCanvas * c;
    // int ihist = 0;
    //     for (auto h2 : h2s) {
    //         for (int i = 0; i < N_CENT_BINS; i++) {
    //             h2->GetXaxis()->SetRange(i+1, i+1);
    //             float max = h2->GetMaximum();
    //             TH1F * h1 = h2->ProjectionY(Form("h1_%s_cent%d", h2->GetName(), i), i+1, i+1);
    //             h1->GetXaxis()->SetTitle("E_{T} [GeV]");
    //             h1->GetYaxis()->SetTitle("P(E_{T})");
    //             h1->Scale(1.0/h1->Integral(), "width");
    //             c = new TCanvas();
    //             h1->Draw();
    //             tex->DrawLatex(0.55, 0.9, Form("%s", h2_names[ihist].c_str()));
    //             tex->DrawLatex(0.55, 0.85, Form("%d-%d%%", int(CENT_BINS[i]), int(CENT_BINS[i+1])));
    //             c->SaveAs(Form("%s/%s_cent%d.png", output_path.c_str(), h2_names[ihist].c_str(), i));
    //         }
    //     }


    // TTree * t = (TTree*)f->Get("T");
    // if (!t) { std::cout << "Tree T not found" << std::endl; exit(1); }
    // // tower_background_energy_recemc
    
    // int cent;
    // t->SetBranchAddress("centrality", &cent);

    // std::vector<float> * bkdgE_cemc = 0;
    // std::vector<float> * bkdgE_ihcal = 0;
    // std::vector<float> * bkdgE_ohcal = 0;
    // t->SetBranchAddress("tower_background_energy_recemc", &bkdgE_cemc);
    // t->SetBranchAddress("tower_background_energy_hcalin", &bkdgE_ihcal);
    // t->SetBranchAddress("tower_background_energy_hcalout", &bkdgE_ohcal);

    // float cemcE, ihcalE, ohcalE;
    // float cemcE_sub1, ihcalE_sub1, ohcalE_sub1;
    // float cemcE_retower;
    // t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC", &cemcE);
    // t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN", &ihcalE);
    // t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT", &ohcalE);
    // t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER", &cemcE_retower);
    // t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &cemcE_sub1);
    // t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALIN_SUB1", &ihcalE_sub1);
    // t->SetBranchAddress("tower_sum_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &ohcalE_sub1);

    // float cemcAvgE, ihcalAvgE, ohcalAvgE;
    // float cemcAvgE_sub1, ihcalAvgE_sub1, ohcalAvgE_sub1;
    // float cemcAvgE_retower;
    // t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC", &cemcAvgE);
    // t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN", &ihcalAvgE);
    // t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT", &ohcalAvgE);
    // t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER", &cemcAvgE_retower);
    // t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &cemcAvgE_sub1);
    // t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALIN_SUB1", &ihcalAvgE_sub1);
    // t->SetBranchAddress("tower_avg_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &ohcalAvgE_sub1);

    // float cemcStdE, ihcalStdE, ohcalStdE;
    // float cemcStdE_sub1, ihcalStdE_sub1, ohcalStdE_sub1;
    // float cemcStdE_retower;
    // t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC", &cemcStdE);
    // t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN", &ihcalStdE);
    // t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT", &ohcalStdE);
    // t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER", &cemcStdE_retower);
    // t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &cemcStdE_sub1);
    // t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALIN_SUB1", &ihcalStdE_sub1);
    // t->SetBranchAddress("tower_std_energy_TOWERINFO_CALIB_HCALOUT_SUB1", &ohcalStdE_sub1);

    // float fCemc, fIhcal, fOhcal;
    // float fCemc_sub1, fIhcal_sub1, fOhcal_sub1;
    // float fCemc_retower;
    // t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC", &fCemc);
    // t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN", &fIhcal);
    // t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT", &fOhcal);
    // t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER", &fCemc_retower);
    // t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &fCemc_sub1);
    // t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALIN_SUB1", &fIhcal_sub1);
    // t->SetBranchAddress("tower_frac_fired_TOWERINFO_CALIB_HCALOUT_SUB1", &fOhcal_sub1);

    // float dCemc, dIhcal, dOhcal;
    // float dCemc_sub1, dIhcal_sub1, dOhcal_sub1;
    // float dCemc_retower;
    // t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC", &dCemc);
    // t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN", &dIhcal);
    // t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT", &dOhcal);
    // t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER", &dCemc_retower);
    // t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_CEMC_RETOWER_SUB1", &dCemc_sub1);
    // t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALIN_SUB1", &dIhcal_sub1);
    // t->SetBranchAddress("tower_frac_dead_TOWERINFO_CALIB_HCALOUT_SUB1", &dOhcal_sub1);




    // int num_windows_full[11];
    // float avg_energy_full[11];
    // float std_energy_full[11];
    // t->SetBranchAddress("num_windows_full", &num_windows_full);
    // t->SetBranchAddress("avg_energy_full", &avg_energy_full);
    // t->SetBranchAddress("std_energy_full", &std_energy_full);
    
    // float rhoA, rhoM;
    // t->SetBranchAddress("rho_val_TowerRho_AREA", &rhoA);
    // t->SetBranchAddress("rho_val_TowerRho_MULT", &rhoM);

    // int nentries = t->GetEntries();
    // NEVENTS = nentries;
    // std::cout << "Entries: " << nentries << std::endl;
    // t->GetEntry(0);
    // f->Close();

   
}








