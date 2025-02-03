#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TMath.h>
#include <TString.h>

#include <iostream>
#include <vector>



void MakeHistograms(const std::string & input_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/JAN30/DATA/BASIC/DATA-ProdA_2023-BASIC-023745.root",
                    const std::string & output_file = "data_all_rc_nop_histos.root")
{

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    bool is_data = input_file.find("data") != std::string::npos;
    std::string suffix = is_data ? "data_" : "_hijing_";
   
   // open input file
    TFile *f = new TFile(input_file.c_str(), "READ");
    if(!f->IsOpen() || f->IsZombie()){ std::cout << "File " << input_file << " is zombie" << std::endl;  return -1; }

    // get tree
    TTree *t = (TTree*)f->Get("T");
    if(!t){ std::cout << "Tree event_info not found in " << input_file << std::endl; return -1; }

    // input variables
    int event_id = 0;
    float zvtx = 0;
    int centrality = 0;
    float mbd_NS = 0;

    // rho variables
    float TowerRho_AREA_rho = 0;
    float TowerRho_MULT_rho = 0;

    // random cone variables
    float RandomCones_Basic_r04_Sub1_eta = 0;
    float RandomCones_Basic_r04_Sub1_phi = 0;
    float RandomCones_Basic_r04_Sub1_pt = 0;
    int RandomCones_Basic_r04_Sub1_nclustered = 0;

    

    float RandomCones_Basic_r04_eta = 0;
    float RandomCones_Basic_r04_phi = 0;
    float RandomCones_Basic_r04_pt = 0;
    int RandomCones_Basic_r04_nclustered = 0;

    
    float RandomCones_RandEtaPhi_r04_Sub1_eta = 0;
    float RandomCones_RandEtaPhi_r04_Sub1_phi = 0;
    float RandomCones_RandEtaPhi_r04_Sub1_pt = 0;
    int RandomCones_RandEtaPhi_r04_Sub1_nclustered = 0;

    float RandomCones_RandEtaPhi_r04_eta = 0;
    float RandomCones_RandEtaPhi_r04_phi = 0;
    float RandomCones_RandEtaPhi_r04_pt = 0;
    int RandomCones_RandEtaPhi_r04_nclustered = 0;

   
    double sum_et_ihcal = 0;
    double sum_et_re_emcal = 0;
    double sum_et_ohcal = 0;
    double sum_et_emcal = 0;
    double sum_et_ihcal_sub1 = 0;
    double sum_et_re_emcal_sub1 = 0;
    double sum_et_ohcal_sub1 = 0;

    // set branch addresses
    t->SetBranchAddress("event_id", &event_id);
    t->SetBranchAddress("zvtx", &zvtx);
    t->SetBranchAddress("centrality", &centrality);
    t->SetBranchAddress("mbd_NS", &mbd_NS);
    t->SetBranchAddress("TowerRho_AREA_rho", &TowerRho_AREA_rho);
    t->SetBranchAddress("TowerRho_MULT_rho", &TowerRho_MULT_rho);
    t->SetBranchAddress("RandomCones_Basic_r04_Sub1_eta", &RandomCones_Basic_r04_Sub1_eta);
    t->SetBranchAddress("RandomCones_Basic_r04_Sub1_phi", &RandomCones_Basic_r04_Sub1_phi);
    t->SetBranchAddress("RandomCones_Basic_r04_Sub1_pt", &RandomCones_Basic_r04_Sub1_pt);
    t->SetBranchAddress("RandomCones_Basic_r04_Sub1_nclustered", &RandomCones_Basic_r04_Sub1_nclustered);
   
    t->SetBranchAddress("RandomCones_Basic_r04_eta", &RandomCones_Basic_r04_eta);
    t->SetBranchAddress("RandomCones_Basic_r04_phi", &RandomCones_Basic_r04_phi);
    t->SetBranchAddress("RandomCones_Basic_r04_pt", &RandomCones_Basic_r04_pt);
    t->SetBranchAddress("RandomCones_Basic_r04_nclustered", &RandomCones_Basic_r04_nclustered);
    t->SetBranchAddress("RandomCones_RandEtaPhi_r04_eta", &RandomCones_RandEtaPhi_r04_eta);
    t->SetBranchAddress("RandomCones_RandEtaPhi_r04_phi", &RandomCones_RandEtaPhi_r04_phi);
    t->SetBranchAddress("RandomCones_RandEtaPhi_r04_pt", &RandomCones_RandEtaPhi_r04_pt);
    t->SetBranchAddress("RandomCones_RandEtaPhi_r04_nclustered", &RandomCones_RandEtaPhi_r04_nclustered);
    t->SetBranchAddress("RandomCones_RandEtaPhi_r04_Sub1_eta", &RandomCones_RandEtaPhi_r04_Sub1_eta);
    t->SetBranchAddress("RandomCones_RandEtaPhi_r04_Sub1_phi", &RandomCones_RandEtaPhi_r04_Sub1_phi);
    t->SetBranchAddress("RandomCones_RandEtaPhi_r04_Sub1_pt", &RandomCones_RandEtaPhi_r04_Sub1_pt);
    t->SetBranchAddress("RandomCones_RandEtaPhi_r04_Sub1_nclustered", &RandomCones_RandEtaPhi_r04_Sub1_nclustered);
   
    t->SetBranchAddress("sum_et_ihcal", &sum_et_ihcal);
    t->SetBranchAddress("sum_et_re_emcal", &sum_et_re_emcal);
    t->SetBranchAddress("sum_et_ohcal", &sum_et_ohcal);
    t->SetBranchAddress("sum_et_emcal", &sum_et_emcal);
    t->SetBranchAddress("sum_et_ihcal_sub1", &sum_et_ihcal_sub1);
    t->SetBranchAddress("sum_et_re_emcal_sub1", &sum_et_re_emcal_sub1);
    t->SetBranchAddress("sum_et_ohcal_sub1", &sum_et_ohcal_sub1);

    TH2D * histo_et_cent = new TH2D("histo_et_cent", "E_{T} vs Centrality; Centrality; E_{T} [GeV]", 100, 0, 100, 2000,0,2000);

    // TH1D * h_et_ihcal = (TH1D*)f->Get("et_ihcal");
    // h_et_ihcal->GetXaxis()->SetTitle("E_{T}^{IHCal} (GeV)");
    // h_et_ihcal->GetYaxis()->SetTitle("Counts");
    // TH1D * h_et_re_emcal = (TH1D*)f->Get("et_re_emcal");
    // h_et_re_emcal->GetXaxis()->SetTitle("E_{T}^{Re-Emcal} (GeV)");
    // h_et_re_emcal->GetYaxis()->SetTitle("Counts");
    // TH1D * h_et_ohcal = (TH1D*)f->Get("et_ohcal");
    // h_et_ohcal->GetXaxis()->SetTitle("E_{T}^{OHCal} (GeV)");
    // h_et_ohcal->GetYaxis()->SetTitle("Counts");
    // TH1D * h_et_emcal = (TH1D*)f->Get("et_emcal");
    // h_et_emcal->GetXaxis()->SetTitle("E_{T}^{Emcal} (GeV)");
    // h_et_emcal->GetYaxis()->SetTitle("Counts");
    // TH1D * h_et_ihcal_sub1 = (TH1D*)f->Get("et_ihcal_sub1");
    // h_et_ihcal_sub1->GetXaxis()->SetTitle("E_{T,sub1}^{IHCal} (GeV)");
    // h_et_ihcal_sub1->GetYaxis()->SetTitle("Counts");
    // TH1D * h_et_re_emcal_sub1 = (TH1D*)f->Get("et_re_emcal_sub1");
    // h_et_re_emcal_sub1->GetXaxis()->SetTitle("E_{T,sub1}^{Re-Emcal} (GeV)");
    // h_et_re_emcal_sub1->GetYaxis()->SetTitle("Counts");
    // TH1D * h_et_ohcal_sub1 = (TH1D*)f->Get("et_ohcal_sub1");
    // h_et_ohcal_sub1->GetXaxis()->SetTitle("E_{T,sub1}^{OHCal} (GeV)");
    // h_et_ohcal_sub1->GetYaxis()->SetTitle("Counts");

    // get number of entries
    int nentries = t->GetEntries();

    // declare histograms
    const int resp_N = 250;
    double max_resp = 60;
    double min_resp = -60;
    double delta_resp = (max_resp - min_resp)/resp_N;
    double resp_bins[resp_N+1];
    for(int i = 0; i < resp_N+1; i++)
    {
        resp_bins[i] = min_resp + i*delta_resp;
    }

    // const int sum_eT_N = 10;
    // double sum_eT_max = 1500;
    // double sum_eT_min = 100;

    // double delta_eT = (sum_eT_max - sum_eT_min)/sum_eT_N;
    // double sum_eT_bins[sum_eT_N+1];
    // double sum_eT_bins[] = {46, 60, 80, 127, 194, 290, 433, 642, 938, 1141, 1712};
    // double sum_eT_bins[] =  {50, 100, 200, 300, 500, 750, 1000, 1250, 1700};
    // 1712, 1141, 938, 782, 642, 532, 433, 357, 290, 240, 194, 156, 127, 101, 80, 67, 60, 50, 46
    double sum_eT_bins[] ={ 46, 50, 60, 67, 80, 101, 127, 156, 194, 240, 290, 357, 433, 532, 642, 782, 938, 1141, 1712};
    const int sum_eT_N = sizeof(sum_eT_bins)/sizeof(double) - 1;
    // for(int i = 0; i < sum_eT_N+1; i++)
    // {
    //     sum_eT_bins[i] = sum_eT_min + i*delta_eT;
    // }

    // const double sum_eT_bins[] = {100, 200, 300, 400, 500, 700, 900, 1200, 1500, 2000, 2500};
    // const int sum_eT_N = sizeof(sum_eT_bins)/sizeof(sum_eT_bins[0]) - 1;

    TH2D * h2_cone_residual_pt_basic_area = new TH2D("h2_cone_residual_pt_basic_area", "#delta p_{T} Basic Area; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);
    TH2D * h2_cone_residual_pt_basic_mult = new TH2D("h2_cone_residual_pt_basic_mult", "#delta p_{T} Basic Mult; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);
    TH2D * h2_cone_residual_pt_basic_sub1 = new TH2D("h2_cone_residual_pt_basic_sub1", "#delta p_{T} Basic Sub1; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);
    TH2D * h2_cone_residual_pt_rand_eta_phi_sub1 = new TH2D("h2_cone_residual_pt_rand_eta_phi_sub1", "#delta p_{T} Rand Eta Phi Sub1; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);
    TH2D * h2_cone_residual_pt_rand_eta_phi_area = new TH2D("h2_cone_residual_pt_rand_eta_phi_area", "#delta p_{T} Rand Eta Phi Area; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);
    TH2D * h2_cone_residual_pt_rand_eta_phi_mult = new TH2D("h2_cone_residual_pt_rand_eta_phi_mult", "#delta p_{T} Rand Eta Phi Mult; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);

    TH1D * h_et_sum_calos = new TH1D("h_et_sum_calos", "Sum E_{T} Calos; #Sigma E_{T} [GeV]; Counts", sum_eT_N, sum_eT_bins);
    TH1D * h_et_sum_sub1_calos = new TH1D("h_et_sum_sub1_calos", "Sum E_{T,sub1} Calos; #Sigma E_{T} [GeV]; Counts", sum_eT_N, sum_eT_bins);

    // loop over entries
    double maxet = 0;
    for (int i = 0; i < nentries; i++)
    {
        t->GetEntry(i);
        double sum_et_all = sum_et_ihcal + sum_et_ohcal + sum_et_emcal;
        double sum_et_all_sub1 = sum_et_ihcal_sub1 + sum_et_ohcal_sub1 + sum_et_re_emcal_sub1;
        float RandomCones_Basic_r04_SubRhoAREA_pt = RandomCones_Basic_r04_pt - TowerRho_AREA_rho*(TMath::Pi()*0.4*0.4);
        float RandomCones_Basic_r04_SubRhoMULT_pt = RandomCones_Basic_r04_pt - TowerRho_MULT_rho*RandomCones_Basic_r04_nclustered;
        h2_cone_residual_pt_basic_area->Fill(RandomCones_Basic_r04_SubRhoAREA_pt , sum_et_all);
        h2_cone_residual_pt_basic_mult->Fill(RandomCones_Basic_r04_SubRhoMULT_pt , sum_et_all);
        h2_cone_residual_pt_basic_sub1->Fill(RandomCones_Basic_r04_Sub1_pt , sum_et_all);

        float RandomCones_RandEtaPhi_r04_SubRhoAREA_pt = RandomCones_RandEtaPhi_r04_pt - TowerRho_AREA_rho*(TMath::Pi()*0.4*0.4);
        float RandomCones_RandEtaPhi_r04_SubRhoMULT_pt = RandomCones_RandEtaPhi_r04_pt - TowerRho_MULT_rho*RandomCones_RandEtaPhi_r04_nclustered;
        h2_cone_residual_pt_rand_eta_phi_area->Fill(RandomCones_RandEtaPhi_r04_SubRhoAREA_pt , sum_et_all);
        h2_cone_residual_pt_rand_eta_phi_mult->Fill(RandomCones_RandEtaPhi_r04_SubRhoMULT_pt , sum_et_all);
        h2_cone_residual_pt_rand_eta_phi_sub1->Fill(RandomCones_RandEtaPhi_r04_Sub1_pt , sum_et_all);

        h_et_sum_calos->Fill(sum_et_all);
        h_et_sum_sub1_calos->Fill(sum_et_all_sub1);

        histo_et_cent->Fill(centrality, sum_et_all);
        if (sum_et_all > maxet)
        {
            maxet = sum_et_all;
        }

    }

        // make 1D sigma histograms
    TH1D * h1_cone_residual_pt_basic_area_sigma = new TH1D("h1_cone_residual_pt_basic_area_sigma", "h1_cone_residual_pt_basic_area_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);
    TH1D * h1_cone_residual_pt_basic_mult_sigma = new TH1D("h1_cone_residual_pt_basic_mult_sigma", "h1_cone_residual_pt_basic_mult_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);
    TH1D * h1_cone_residual_pt_basic_sub1_sigma = new TH1D("h1_cone_residual_pt_basic_sub1_sigma", "h1_cone_residual_pt_basic_sub1_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);
    TH1D * h1_cone_residual_pt_rand_eta_phi_area_sigma = new TH1D("h1_cone_residual_pt_rand_eta_phi_area_sigma", "h1_cone_residual_pt_rand_eta_phi_area_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);
    TH1D * h1_cone_residual_pt_rand_eta_phi_mult_sigma = new TH1D("h1_cone_residual_pt_rand_eta_phi_mult_sigma", "h1_cone_residual_pt_rand_eta_phi_mult_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);
    TH1D * h1_cone_residual_pt_rand_eta_phi_sub1_sigma = new TH1D("h1_cone_residual_pt_rand_eta_phi_sub1_sigma", "h1_cone_residual_pt_rand_eta_phi_sub1_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);

    // make vectors of 2D histograms for easier access
    std::vector<TH2D*> h2_cone_residual_pt_vec = {h2_cone_residual_pt_basic_area, 
                                                  h2_cone_residual_pt_basic_mult, 
                                                  h2_cone_residual_pt_basic_sub1,
                                                    h2_cone_residual_pt_rand_eta_phi_area,
                                                    h2_cone_residual_pt_rand_eta_phi_mult,
                                                    h2_cone_residual_pt_rand_eta_phi_sub1};

    // 1D vector
    std::vector<TH1D*> h1_cone_residual_pt_sigma_vec = {h1_cone_residual_pt_basic_area_sigma, 
                                                        h1_cone_residual_pt_basic_mult_sigma, 
                                                        h1_cone_residual_pt_basic_sub1_sigma,
                                                        h1_cone_residual_pt_rand_eta_phi_area_sigma,
                                                        h1_cone_residual_pt_rand_eta_phi_mult_sigma,
                                                        h1_cone_residual_pt_rand_eta_phi_sub1_sigma};



    // create output file
    TFile *fout = new TFile(output_file.c_str(), "RECREATE");
    // TCanvas *c = new TCanvas("c", "c", 800, 600);

    // Calculate nx (columns) and ny (rows)
    int nx = sum_eT_N % 2 == 0 ? sum_eT_N / 2 : sum_eT_N / 2 + 1;
    int ny = 2;
    std::cout << "nx = " << nx << ", ny = " << ny << std::endl;

    // Define desired pad dimensions (in pixels)
    const int pad_width = 300;  // Width of a single pad
    const int pad_height = 300; // Height of a single pad

    // Calculate canvas size
    int canvas_width = nx * pad_width;
    int canvas_height = ny * pad_height;

    // Create output file and canvas
   
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    // gStyle->SetOptTitle(0);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetRightMargin(0.05);
    // gPad->SetBottomMargin(0.15);
    // gPad->SetTopMargin(0.05);
    
    // double margin = 0.05; // General margin size
    // double pad_width_fraction = 1.0 / nx;
    // double pad_height_fraction = 1.0 / ny;

    // for (int i = 1; i <= nx * ny; ++i) {
    //     TPad *pad = (TPad*)cn->GetPad(i);
    //     if (pad) {
    //         pad->SetLeftMargin(margin / pad_width_fraction);
    //         pad->SetRightMargin(margin / pad_width_fraction);
    //         pad->SetBottomMargin(margin / pad_height_fraction);
    //         pad->SetTopMargin(margin / pad_height_fraction);
    //     }
    // }
    
    std::vector<std::string> sphenix_tags = {};
    sphenix_tags.push_back("#it{#bf{sPHENIX}} Internal");
    if(is_data){
        sphenix_tags.push_back("Au+Au #sqrt{s_{NN}} = 200 GeV");
    } else {
        sphenix_tags.push_back("HIJING MDC2");
    }
    sphenix_tags.push_back("Random Cone #it{R} = 0.4");
    sphenix_tags.push_back("|#eta_{cone}| < 0.7");
    double tagx = 0.6;
    double tagy = 0.9;

    TLegend * legp = new TLegend(0.1, 0.7, 0.3, 0.9);
    legp->SetBorderSize(0);
    legp->SetFillStyle(0);

    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(43);


    // TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    // loop over 2D histograms, split by centrality and get RMS
    for (unsigned int ihist = 0; ihist < h2_cone_residual_pt_vec.size(); ihist++)
    {
        // get 2D histogram
        TCanvas * cn = new TCanvas("cn", "cn", canvas_width, canvas_height);
        cn->Divide(nx, ny);
        TH2D * h2_temp = h2_cone_residual_pt_vec.at(ihist);
        for (int ix = 1; ix <= sum_eT_N; ix++)
        {
            h2_temp->GetYaxis()->SetRange(ix, ix);
            std::string hist_name = h2_temp->GetName();
            TH1D * h1_temp = (TH1D*)h2_temp->ProjectionX(Form("%s_sumEtbin%d", hist_name.c_str(), ix));
            hist_name = h1_temp->GetName();
            h1_temp->GetXaxis()->SetTitle("#delta E_{T}");
            h1_temp->GetYaxis()->SetTitle("Counts");

            TF1 * func = new TF1("func", "gaus", -10, 10);
            
            h1_temp->Fit(func, "Q", "", -10, 10);
            h1_temp->Fit(func, "Q", "", func->GetParameter(1) - 2*func->GetParameter(2), func->GetParameter(1) + 2*func->GetParameter(2));

            double sigma = func->GetParameter(2);
            double sigma_err = func->GetParError(2);

            h1_cone_residual_pt_sigma_vec.at(ihist)->SetBinContent(ix, sigma);
            h1_cone_residual_pt_sigma_vec.at(ihist)->SetBinError(ix, sigma_err);

            // write histograms
            // h1_temp->Draw();
            // func->SetLineColor(kRed);
            // func->Draw("same");
            // leg->AddEntry(h1_temp, Form("#sigma = %.2f #pm %.2f", sigma, sigma_err));
            // leg->Draw("same");

            // c->SaveAs(Form("plots/%s_%s_sumEtbin%d.png", suffix.c_str(), h2_temp->GetName(), ix));
            // leg->Clear();
            
            // update all canvases
            cn->cd( (ix-1)/2 + 1);
            gPad->SetLogy();
            std::string etbin_str = Form("%d < #Sigma E_{T} < %d", (int)sum_eT_bins[ix], (int)sum_eT_bins[ix+1]);
            // draw latex
            h1_temp->GetXaxis()->SetRangeUser(-20, 20);
            h1_temp->GetYaxis()->SetRangeUser(0.001, 10*h1_temp->GetMaximum());
            // set axis labels
            h1_temp->GetXaxis()->SetTitle("#delta E_{T}");
            h1_temp->GetYaxis()->SetTitle("A.U.");
            h1_temp->Draw();
            func->Draw("same");
            legp->AddEntry(h1_temp, etbin_str.c_str());
            legp->AddEntry(func, Form("#sigma = %.1f #pm %.1f", sigma, sigma_err));
            legp->Draw("same");

            if(ix == 1){
                tex->DrawLatex(tagx, tagy-0.05, sphenix_tags.at(0).c_str());
                tex->DrawLatex(tagx, tagy-0.1, sphenix_tags.at(1).c_str());
                tex->DrawLatex(tagx, tagy-0.15, sphenix_tags.at(2).c_str());
                tex->DrawLatex(tagx, tagy-0.2, sphenix_tags.at(3).c_str());
            }

            


        }
        std::string tmpname = h2_temp->GetName();
        // remove h2_  from name
        tmpname.erase(0, 3);
        cn->SaveAs(Form("%s_%s.png", suffix.c_str(), tmpname.c_str()));

    }

    // write 2D histograms
    fout->cd();
    for (unsigned int ihist = 0; ihist < h2_cone_residual_pt_vec.size(); ihist++)
    {
        h2_cone_residual_pt_vec.at(ihist)->Write();
    }
    h_et_sum_calos->Write();
    h_et_sum_sub1_calos->Write();
    for (unsigned int ihist = 0; ihist < h1_cone_residual_pt_sigma_vec.size(); ihist++)
    {
        h1_cone_residual_pt_sigma_vec.at(ihist)->Write();
    }
    

    // write 1D sigma histograms
    fout->cd();
    for (unsigned int ihist = 0; ihist < h1_cone_residual_pt_sigma_vec.size(); ihist++)
    {
        h1_cone_residual_pt_sigma_vec.at(ihist)->Write();
    }

    TProfile * p_et_cent = histo_et_cent->ProfileX();
    TH1D * h_et_cent = (TH1D*)p_et_cent->ProjectionX("h_et_cent");
    std::vector<double> cent_bins = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55,
                                     60, 65, 70, 75, 80, 85, 90, 95, 100};
                                    // 1712, 1141, 938, 642, 433, 290, 194, 127, 80, 60, 46, 0,

    std::vector<double> et_bins = {};
    int ncent_bins = cent_bins.size() - 1;
    for (int i = 0; i < ncent_bins+1; i++)
    {
        if (i == 0){
            et_bins.push_back(maxet);
            continue;
        }
        double cent = cent_bins.at(i);
        // find inhistogram bin
        int bin = h_et_cent->FindBin(cent);
        double et = h_et_cent->GetBinContent(bin);
        et_bins.push_back(et);
    }
    std::cout << "Centrality bins: ";
    for (unsigned int i = 0; i < et_bins.size(); i++)
    {
        std::cout << (int)et_bins.at(i) << ", ";
    }
    std::cout << std::endl;
    h_et_cent->GetXaxis()->SetTitle("Centrality");
    h_et_cent->GetYaxis()->SetTitle("E_{T} [GeV]");
    h_et_cent->Write();
    // h_et_ihcal->Write();
    // h_et_re_emcal->Write();
    // h_et_ohcal->Write();
    // h_et_emcal->Write();
    // h_et_ihcal_sub1->Write();
    // h_et_re_emcal_sub1->Write();
    // h_et_ohcal_sub1->Write();

    fout->Close();

    return;

 

}

// void MakeHistograms(){

//     const std::string & input_file_hijing = "/sphenix/user/tmengel/ppg04/condor/hijing_all_rc_nop.root";
//     const std::string & output_file_hijing = "hijing_all_rc_nop_histos.root";
//     MakeHistos(input_file_hijing, output_file_hijing);

//     const std::string & input_file_data = "/sphenix/user/tmengel/ppg04/condor/data_all_rc_nop.root";
//     const std::string & output_file_data = "data_all_rc_nop_histos.root";
//     MakeHistos(input_file_data, output_file_data);


//     return ;

// 
// }

// void MakeHistograms(const std::string & input_file = "/sphenix/user/tmengel/ppg04/condor/hijing_all_rc.root",
//                     const std::string & output_file = "hijing_all_rc_histos.root")
// {

//     TH1::SetDefaultSumw2();
//     TH2::SetDefaultSumw2();
   
//    // open input file
//     TFile *f = new TFile(input_file.c_str(), "READ");
//     if(!f->IsOpen() || f->IsZombie()){ std::cout << "File " << input_file << " is zombie" << std::endl;  return -1; }

//     // get tree
//     TTree *t = (TTree*)f->Get("event_info");
//     if(!t){ std::cout << "Tree event_info not found in " << input_file << std::endl; return -1; }

//     // input variables
//     int event_id = 0;
//     float zvtx = 0;
//     int centrality = 0;
//     float mbd_NS = 0;

//     // rho variables
//     float TowerRho_AREA_rho = 0;
//     float TowerRho_MULT_rho = 0;

//     // random cone variables
//     float RandomCones_Basic_r04_Sub1_eta = 0;
//     float RandomCones_Basic_r04_Sub1_phi = 0;
//     float RandomCones_Basic_r04_Sub1_pt = 0;
//     int RandomCones_Basic_r04_Sub1_nclustered = 0;

//     float RandomCones_Basic_r04_SubRhoAREA_eta = 0;
//     float RandomCones_Basic_r04_SubRhoAREA_phi = 0;
//     float RandomCones_Basic_r04_SubRhoAREA_pt = 0;
//     int RandomCones_Basic_r04_SubRhoAREA_nclustered = 0;

//     float RandomCones_Basic_r04_SubRhoMULT_eta = 0;
//     float RandomCones_Basic_r04_SubRhoMULT_phi = 0;
//     float RandomCones_Basic_r04_SubRhoMULT_pt = 0;
//     int RandomCones_Basic_r04_SubRhoMULT_nclustered = 0;

//     float RandomCones_RandEtaPhi_r04_Sub1_eta = 0;
//     float RandomCones_RandEtaPhi_r04_Sub1_phi = 0;
//     float RandomCones_RandEtaPhi_r04_Sub1_pt = 0;
//     int RandomCones_RandEtaPhi_r04_Sub1_nclustered = 0;

//     float RandomCones_RandEtaPhi_r04_SubRhoAREA_eta = 0;
//     float RandomCones_RandEtaPhi_r04_SubRhoAREA_phi = 0;
//     float RandomCones_RandEtaPhi_r04_SubRhoAREA_pt = 0;
//     int RandomCones_RandEtaPhi_r04_SubRhoAREA_nclustered = 0;

//     float RandomCones_RandEtaPhi_r04_SubRhoMULT_eta = 0;
//     float RandomCones_RandEtaPhi_r04_SubRhoMULT_phi = 0;
//     float RandomCones_RandEtaPhi_r04_SubRhoMULT_pt = 0;
//     int RandomCones_RandEtaPhi_r04_SubRhoMULT_nclustered = 0;


//     double sum_et_ihcal = 0;
//     double sum_et_re_emcal = 0;
//     double sum_et_ohcal = 0;
//     double sum_et_emcal = 0;
//     double sum_et_ihcal_sub1 = 0;
//     double sum_et_re_emcal_sub1 = 0;
//     double sum_et_ohcal_sub1 = 0;

//     // set branch addresses
//     t->SetBranchAddress("event_id", &event_id);
//     t->SetBranchAddress("zvtx", &zvtx);
//     t->SetBranchAddress("centrality", &centrality);
//     t->SetBranchAddress("mbd_NS", &mbd_NS);
//     t->SetBranchAddress("TowerRho_AREA_rho", &TowerRho_AREA_rho);
//     t->SetBranchAddress("TowerRho_MULT_rho", &TowerRho_MULT_rho);
//     t->SetBranchAddress("RandomCones_Basic_r04_Sub1_eta", &RandomCones_Basic_r04_Sub1_eta);
//     t->SetBranchAddress("RandomCones_Basic_r04_Sub1_phi", &RandomCones_Basic_r04_Sub1_phi);
//     t->SetBranchAddress("RandomCones_Basic_r04_Sub1_pt", &RandomCones_Basic_r04_Sub1_pt);
//     t->SetBranchAddress("RandomCones_Basic_r04_Sub1_nclustered", &RandomCones_Basic_r04_Sub1_nclustered);
//     t->SetBranchAddress("RandomCones_Basic_r04_SubRhoAREA_eta", &RandomCones_Basic_r04_SubRhoAREA_eta);
//     t->SetBranchAddress("RandomCones_Basic_r04_SubRhoAREA_phi", &RandomCones_Basic_r04_SubRhoAREA_phi);
//     t->SetBranchAddress("RandomCones_Basic_r04_SubRhoAREA_pt", &RandomCones_Basic_r04_SubRhoAREA_pt);
//     t->SetBranchAddress("RandomCones_Basic_r04_SubRhoAREA_nclustered", &RandomCones_Basic_r04_SubRhoAREA_nclustered);
//     t->SetBranchAddress("RandomCones_Basic_r04_SubRhoMULT_eta", &RandomCones_Basic_r04_SubRhoMULT_eta);
//     t->SetBranchAddress("RandomCones_Basic_r04_SubRhoMULT_phi", &RandomCones_Basic_r04_SubRhoMULT_phi);
//     t->SetBranchAddress("RandomCones_Basic_r04_SubRhoMULT_pt", &RandomCones_Basic_r04_SubRhoMULT_pt);
//     t->SetBranchAddress("RandomCones_Basic_r04_SubRhoMULT_nclustered", &RandomCones_Basic_r04_SubRhoMULT_nclustered);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_Sub1_eta", &RandomCones_RandEtaPhi_r04_Sub1_eta);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_Sub1_phi", &RandomCones_RandEtaPhi_r04_Sub1_phi);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_Sub1_pt", &RandomCones_RandEtaPhi_r04_Sub1_pt);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_Sub1_nclustered", &RandomCones_RandEtaPhi_r04_Sub1_nclustered);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_SubRhoAREA_eta", &RandomCones_RandEtaPhi_r04_SubRhoAREA_eta);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_SubRhoAREA_phi", &RandomCones_RandEtaPhi_r04_SubRhoAREA_phi);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_SubRhoAREA_pt", &RandomCones_RandEtaPhi_r04_SubRhoAREA_pt);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_SubRhoAREA_nclustered", &RandomCones_RandEtaPhi_r04_SubRhoAREA_nclustered);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_SubRhoMULT_eta", &RandomCones_RandEtaPhi_r04_SubRhoMULT_eta);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_SubRhoMULT_phi", &RandomCones_RandEtaPhi_r04_SubRhoMULT_phi);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_SubRhoMULT_pt", &RandomCones_RandEtaPhi_r04_SubRhoMULT_pt);
//     t->SetBranchAddress("RandomCones_RandEtaPhi_r04_SubRhoMULT_nclustered", &RandomCones_RandEtaPhi_r04_SubRhoMULT_nclustered);
//     t->SetBranchAddress("sum_et_ihcal", &sum_et_ihcal);
//     t->SetBranchAddress("sum_et_re_emcal", &sum_et_re_emcal);
//     t->SetBranchAddress("sum_et_ohcal", &sum_et_ohcal);
//     t->SetBranchAddress("sum_et_emcal", &sum_et_emcal);
//     t->SetBranchAddress("sum_et_ihcal_sub1", &sum_et_ihcal_sub1);
//     t->SetBranchAddress("sum_et_re_emcal_sub1", &sum_et_re_emcal_sub1);
//     t->SetBranchAddress("sum_et_ohcal_sub1", &sum_et_ohcal_sub1);



//     TH1D * h_et_ihcal = (TH1D*)f->Get("et_ihcal");
//     h_et_ihcal->GetXaxis()->SetTitle("E_{T}^{IHCal} (GeV)");
//     h_et_ihcal->GetYaxis()->SetTitle("Counts");
//     TH1D * h_et_re_emcal = (TH1D*)f->Get("et_re_emcal");
//     h_et_re_emcal->GetXaxis()->SetTitle("E_{T}^{Re-Emcal} (GeV)");
//     h_et_re_emcal->GetYaxis()->SetTitle("Counts");
//     TH1D * h_et_ohcal = (TH1D*)f->Get("et_ohcal");
//     h_et_ohcal->GetXaxis()->SetTitle("E_{T}^{OHCal} (GeV)");
//     h_et_ohcal->GetYaxis()->SetTitle("Counts");
//     TH1D * h_et_emcal = (TH1D*)f->Get("et_emcal");
//     h_et_emcal->GetXaxis()->SetTitle("E_{T}^{Emcal} (GeV)");
//     h_et_emcal->GetYaxis()->SetTitle("Counts");
//     TH1D * h_et_ihcal_sub1 = (TH1D*)f->Get("et_ihcal_sub1");
//     h_et_ihcal_sub1->GetXaxis()->SetTitle("E_{T,sub1}^{IHCal} (GeV)");
//     h_et_ihcal_sub1->GetYaxis()->SetTitle("Counts");
//     TH1D * h_et_re_emcal_sub1 = (TH1D*)f->Get("et_re_emcal_sub1");
//     h_et_re_emcal_sub1->GetXaxis()->SetTitle("E_{T,sub1}^{Re-Emcal} (GeV)");
//     h_et_re_emcal_sub1->GetYaxis()->SetTitle("Counts");
//     TH1D * h_et_ohcal_sub1 = (TH1D*)f->Get("et_ohcal_sub1");
//     h_et_ohcal_sub1->GetXaxis()->SetTitle("E_{T,sub1}^{OHCal} (GeV)");
//     h_et_ohcal_sub1->GetYaxis()->SetTitle("Counts");

//     // get number of entries
//     int nentries = t->GetEntries();

//     // declare histograms
//     const int resp_N = 250;
//     double max_resp = 60;
//     double min_resp = -60;
//     double delta_resp = (max_resp - min_resp)/resp_N;
//     double resp_bins[resp_N+1];
//     for(int i = 0; i < resp_N+1; i++)
//     {
//         resp_bins[i] = min_resp + i*delta_resp;
//     }

//     // const int sum_eT_N = 20;
//     // double sum_eT_max = 2700;
//     // double sum_eT_min = 0;

//     // double delta_eT = (sum_eT_max - sum_eT_min)/sum_eT_N;
//     // double sum_eT_bins[sum_eT_N+1];
//     // for(int i = 0; i < sum_eT_N+1; i++)
//     // {
//     //     sum_eT_bins[i] = sum_eT_min + i*delta_eT;
//     // }

//     const double sum_eT_bins[] = {100, 200, 300, 400, 500, 700, 900,1200,1500,2700};
//     const int sum_eT_N = sizeof(sum_eT_bins)/sizeof(sum_eT_bins[0]) - 1;

//     TH2D * h2_cone_residual_pt_basic_area = new TH2D("h2_cone_residual_pt_basic_area", "#delta p_{T} Basic Area; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);
//     TH2D * h2_cone_residual_pt_basic_mult = new TH2D("h2_cone_residual_pt_basic_mult", "#delta p_{T} Basic Mult; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);
//     TH2D * h2_cone_residual_pt_basic_sub1 = new TH2D("h2_cone_residual_pt_basic_sub1", "#delta p_{T} Basic Sub1; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);
//     TH2D * h2_cone_residual_pt_rand_eta_phi_sub1 = new TH2D("h2_cone_residual_pt_rand_eta_phi_sub1", "#delta p_{T} Rand Eta Phi Sub1; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);
//     TH2D * h2_cone_residual_pt_rand_eta_phi_area = new TH2D("h2_cone_residual_pt_rand_eta_phi_area", "#delta p_{T} Rand Eta Phi Area; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);
//     TH2D * h2_cone_residual_pt_rand_eta_phi_mult = new TH2D("h2_cone_residual_pt_rand_eta_phi_mult", "#delta p_{T} Rand Eta Phi Mult; #delta E_{T} [GeV]; #Sigma E_{T} [GeV]", resp_N, resp_bins, sum_eT_N, sum_eT_bins);

//     TH1D * h_et_sum_calos = new TH1D("h_et_sum_calos", "Sum E_{T} Calos; #Sigma E_{T} [GeV]; Counts", sum_eT_N, sum_eT_bins);
//     TH1D * h_et_sum_sub1_calos = new TH1D("h_et_sum_sub1_calos", "Sum E_{T,sub1} Calos; #Sigma E_{T} [GeV]; Counts", sum_eT_N, sum_eT_bins);

//     // loop over entries
//     for (int i = 0; i < nentries; i++)
//     {
//         t->GetEntry(i);
//         double sum_et_all = sum_et_ihcal + sum_et_ohcal + sum_et_emcal;
//         double sum_et_all_sub1 = sum_et_ihcal_sub1 + sum_et_ohcal_sub1 + sum_et_re_emcal_sub1;
//         h2_cone_residual_pt_basic_area->Fill(RandomCones_Basic_r04_SubRhoAREA_pt , sum_et_all);
//         h2_cone_residual_pt_basic_mult->Fill(RandomCones_Basic_r04_SubRhoMULT_pt , sum_et_all);
//         h2_cone_residual_pt_basic_sub1->Fill(RandomCones_Basic_r04_Sub1_pt , sum_et_all);

//         h2_cone_residual_pt_rand_eta_phi_area->Fill(RandomCones_RandEtaPhi_r04_SubRhoAREA_pt , sum_et_all);
//         h2_cone_residual_pt_rand_eta_phi_mult->Fill(RandomCones_RandEtaPhi_r04_SubRhoMULT_pt , sum_et_all);
//         h2_cone_residual_pt_rand_eta_phi_sub1->Fill(RandomCones_RandEtaPhi_r04_Sub1_pt , sum_et_all);

//         h_et_sum_calos->Fill(sum_et_all);
//         h_et_sum_sub1_calos->Fill(sum_et_all_sub1);

//     }

//         // make 1D sigma histograms
//     TH1D * h1_cone_residual_pt_basic_area_sigma = new TH1D("h1_cone_residual_pt_basic_area_sigma", "h1_cone_residual_pt_basic_area_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);
//     TH1D * h1_cone_residual_pt_basic_mult_sigma = new TH1D("h1_cone_residual_pt_basic_mult_sigma", "h1_cone_residual_pt_basic_mult_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);
//     TH1D * h1_cone_residual_pt_basic_sub1_sigma = new TH1D("h1_cone_residual_pt_basic_sub1_sigma", "h1_cone_residual_pt_basic_sub1_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);
//     TH1D * h1_cone_residual_pt_rand_eta_phi_area_sigma = new TH1D("h1_cone_residual_pt_rand_eta_phi_area_sigma", "h1_cone_residual_pt_rand_eta_phi_area_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);
//     TH1D * h1_cone_residual_pt_rand_eta_phi_mult_sigma = new TH1D("h1_cone_residual_pt_rand_eta_phi_mult_sigma", "h1_cone_residual_pt_rand_eta_phi_mult_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);
//     TH1D * h1_cone_residual_pt_rand_eta_phi_sub1_sigma = new TH1D("h1_cone_residual_pt_rand_eta_phi_sub1_sigma", "h1_cone_residual_pt_rand_eta_phi_sub1_sigma; #Sigma E_{T} [GeV]; #sigma(#delta p_{T})", sum_eT_N, sum_eT_bins);

//     // make vectors of 2D histograms for easier access
//     std::vector<TH2D*> h2_cone_residual_pt_vec = {h2_cone_residual_pt_basic_area, 
//                                                   h2_cone_residual_pt_basic_mult, 
//                                                   h2_cone_residual_pt_basic_sub1,
//                                                     h2_cone_residual_pt_rand_eta_phi_area,
//                                                     h2_cone_residual_pt_rand_eta_phi_mult,
//                                                     h2_cone_residual_pt_rand_eta_phi_sub1};

//     // 1D vector
//     std::vector<TH1D*> h1_cone_residual_pt_sigma_vec = {h1_cone_residual_pt_basic_area_sigma, 
//                                                         h1_cone_residual_pt_basic_mult_sigma, 
//                                                         h1_cone_residual_pt_basic_sub1_sigma,
//                                                         h1_cone_residual_pt_rand_eta_phi_area_sigma,
//                                                         h1_cone_residual_pt_rand_eta_phi_mult_sigma,
//                                                         h1_cone_residual_pt_rand_eta_phi_sub1_sigma};



//     // create output file
//     TFile *fout = new TFile(output_file.c_str(), "RECREATE");
//     TCanvas *c = new TCanvas("c", "c", 800, 600);
//     TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
//     // loop over 2D histograms, split by centrality and get RMS
//     for (unsigned int ihist = 0; ihist < h2_cone_residual_pt_vec.size(); ihist++)
//     {
//         // get 2D histogram
//         TH2D * h2_temp = h2_cone_residual_pt_vec.at(ihist);
//         for (int ix = 1; ix <= sum_eT_N; ix++)
//         {
//             h2_temp->GetYaxis()->SetRange(ix, ix);
//             std::string hist_name = h2_temp->GetName();
//             TH1D * h1_temp = (TH1D*)h2_temp->ProjectionX(Form("%s_sumEtbin%d", hist_name.c_str(), ix));
//             hist_name = h1_temp->GetName();
//             h1_temp->GetXaxis()->SetTitle("#delta E_{T}");
//             h1_temp->GetYaxis()->SetTitle("Counts");

//             TF1 * func = new TF1("func", "gaus", -10, 10);
            
//             h1_temp->Fit(func, "Q", "", -10, 10);
//             h1_temp->Fit(func, "Q", "", func->GetParameter(1) - 2*func->GetParameter(2), func->GetParameter(1) + 2*func->GetParameter(2));

//             double sigma = func->GetParameter(2);
//             double sigma_err = func->GetParError(2);

//             h1_cone_residual_pt_sigma_vec.at(ihist)->SetBinContent(ix, sigma);
//             h1_cone_residual_pt_sigma_vec.at(ihist)->SetBinError(ix, sigma_err);

//             // write histograms
//             h1_temp->Draw();
//             func->SetLineColor(kRed);
//             func->Draw("same");
//             leg->AddEntry(h1_temp, Form("#sigma = %.2f #pm %.2f", sigma, sigma_err));
//             leg->Draw("same");

//             c->SaveAs(Form("plots/%s_%s_sumEtbin%d.png", h1_temp->GetName(), h2_temp->GetName(), ix));
//             leg->Clear();
//             delete h1_temp;
//         }

//     }

//     // write 2D histograms
//     fout->cd();
//     for (unsigned int ihist = 0; ihist < h2_cone_residual_pt_vec.size(); ihist++)
//     {
//         h2_cone_residual_pt_vec.at(ihist)->Write();
//     }
//     h_et_sum_calos->Write();
//     h_et_sum_sub1_calos->Write();
//     for (unsigned int ihist = 0; ihist < h1_cone_residual_pt_sigma_vec.size(); ihist++)
//     {
//         h1_cone_residual_pt_sigma_vec.at(ihist)->Write();
//     }
    

//     // write 1D sigma histograms
//     fout->cd();
//     for (unsigned int ihist = 0; ihist < h1_cone_residual_pt_sigma_vec.size(); ihist++)
//     {
//         h1_cone_residual_pt_sigma_vec.at(ihist)->Write();
//     }

//     h_et_ihcal->Write();
//     h_et_re_emcal->Write();
//     h_et_ohcal->Write();
//     h_et_emcal->Write();
//     h_et_ihcal_sub1->Write();
//     h_et_re_emcal_sub1->Write();
//     h_et_ohcal_sub1->Write();

//     fout->Close();

//     return;

 

// }