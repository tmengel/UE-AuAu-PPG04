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


const unsigned int k_window_array_size = 11; 
const std::vector < std::pair < unsigned int, unsigned int > > k_calo_window_dims_cemc_geom = {
        {1,1}, {2,2}, {5,5}, {10,10}, {15,15}, {22,23}, {30,29}, {37,37}, {45,45}, {52,52}, {60,60}
};
const std::vector < std::pair < unsigned int, unsigned int > > k_calo_window_dims_hcal_geom = {
    {1,1}, {2,2}, {3,4}, {5,6}, {7,8}, {9,10}, {11,12}, {13,13}, {15,15}, {17,17}, {20,20}
}; 

const float CENT_BINS[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 100};
const int N_CENT_BINS = sizeof(CENT_BINS)/sizeof(CENT_BINS[0]) - 1;

bool IS_DATA = false;
const std::string sPHENIX_Tag = "#it{#bf{sPHENIX}} Internal";
std::string DataType_Tag;


const int COLORS[] = {kBlack, kRed , kBlue, kGreen+2, kViolet, kCyan, kOrange+2, kMagenta+2, kAzure-2};
const int MARKERS[] = { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kOpenCircle, kOpenSquare, kOpenTriangleUp};

const float MARKER_SIZE = 1.2;
const float LINE_WIDTH = 2.0;


int NEVENTS;
std::string global_plots;
std::string calo_plots;
std::string calo_window_plots;
std::string bkgd_plots;
std::string random_cone_plots;

TFile * histo_file {nullptr};
TFile * f {nullptr};

void ConfigureOutputDirs(std::string input_file_base, std::string plotting_dir = "plots/"){
    
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

std::pair<int,int> GetWindowDimFromString(std::string hist_name){
    std::pair<int,int> window_dim;
    std::string sWindow = hist_name.substr(hist_name.find_last_of("_") + 1);
    std::string sWindowX = sWindow.substr(0, sWindow.find("x"));
    std::string sWindowY = sWindow.substr(sWindow.find("x") + 1);
    window_dim.first = std::stoi(sWindowX);
    window_dim.second = std::stoi(sWindowY);
    return window_dim;
}

void CaloWindowMultiPanel(std::vector<TH2F*> h2_vec, const std::vector<std::string> h2_titles,
                        const std::string output_location, const std::string pngbase,
                        const std::string xaxis_title, const std::string yaxis_title, 
                        const double tagx, const double tagy, const double dty, 
                        const double lx1, const double lx2, const double ly1, const double ly2,
                        float minx, float maxx, const float miny, const float maxy,
                        bool sym_x = false )
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
                        bool sym_x = false, bool logy = false )
{

    if ( !gSystem->OpenDirectory(output_location.c_str()) ) {
        gSystem->mkdir(output_location.c_str(), true);
    }

    auto window_dim = GetWindowDimFromString(h3s.at(0)->GetName());
    std::vector<std::string> sTags = {sPHENIX_Tag};
    sTags.push_back(Form("%d #times %d Towers", window_dim.first, window_dim.second));

    double tagy_start = tagy;
    int NX_PADS= 3, NY_PADS = 1;

    TCanvas * c = new TCanvas("c", "c", 400*NX_PADS, 400*NY_PADS);
    c->Divide(NX_PADS, NY_PADS);
    for (int i = 1; i <= N_CENT_BINS; i++) {
        c->cd(i);
        // gPad->SetLogy(logy);
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

    for ( unsigned int i = 0; i < h3s.size(); i++ ){
        h3s.at(i)->SetTitle(h3_titles.at(i).c_str());
    }


    // // get min and max x values for all histograms
    // if ( minx == 0 & maxx ==0) {
    //     for ( const auto & h3 : h3s ) {
    //         if ( !h3 ) {
    //             continue;
    //         }
    //         for ( int icent = 0; icent < N_CENT_BINS; icent++ ) {
    //             h3->GetZaxis()->SetRangeUser(CENT_BINS[icent], CENT_BINS[icent+1]);
    //             TH2F * h2 = (TH2F *) h3->Project3D("xy");
    //             h2->SetName(Form("h2_%s_cent%d", h3->GetName(), icent));
    //             TProfile * p1 = (TProfile *) h2->ProfileX(Form("p1_%s_cent%d", h2->GetName(), icent));
    //             TH1F * h1 = (TH1F *) p1->ProjectionX( Form("h1_%s_cent%d", p1->GetName(), icent));

    //             if ( h1->Integral() == 0 ) { continue; }
    //             h1->Scale(1./h1->Integral());
    //             int max_bin = h1->FindLastBinAbove(miny);
    //             int min_bin = h1->FindFirstBinAbove(miny);
    //             if ( max_bin < 0 ) { max_bin = h1->GetNbinsX(); }
    //             if ( min_bin < 0 ) { min_bin = 1; }
    //             if ( max_bin + 1 < h1->GetNbinsX() ) { max_bin++; }
    //             if ( min_bin - 1 > 0 ) { min_bin--; }
    //             if ( h1->GetXaxis()->GetBinCenter(max_bin) > maxx ) {
    //                 maxx = h1->GetXaxis()->GetBinCenter(max_bin);
    //             }
    //             if ( h1->GetXaxis()->GetBinCenter(min_bin) < minx ) {
    //                 minx = h1->GetXaxis()->GetBinCenter(min_bin);
    //             }
    //             if (sym_x) {
    //                 minx = -maxx;
    //             }
    //         }
    //     }
    // }
    bool is_empty_pannel = true;
    for ( unsigned int ihist = 0; ihist < h3s.size(); ihist++ ) {
            if ( !h3s.at(ihist) ) {
                continue;
            }
        tagy_start = tagy;
        c->cd(ihist+1);
        for ( int icent = 0; icent < N_CENT_BINS; icent++ ) {
            
            std::string sCent = Form("%d-%d%% Central", (int)CENT_BINS[icent], (int)CENT_BINS[icent+1]);
            h3s.at(ihist)->GetZaxis()->SetRangeUser(CENT_BINS[icent], CENT_BINS[icent+1]);
            TH2F * h2 = (TH2F *) h3s.at(ihist)->Project3D("xy");
            h2->SetName(Form("h2_%s_cent%d", h2->GetName(), icent));
            TProfile * p1 = (TProfile *) h2->ProfileX(Form("p1_%s_cent%d", h2->GetName(), icent));
            TH1F * h1 = (TH1F *) p1->ProjectionX( Form("h1_%s_cent%d", h2->GetName(), icent) );
            if ( h1->Integral() == 0 ) { continue; }
            h1->Scale(1./h1->Integral());
            is_empty_pannel= false;
                
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

            if ( icent == 0 ) {
                h1->Draw("");
            } else {
                h1->Draw("SAME");
            }
            if ( ihist == 0 ) {
                leg->AddEntry( h1, sCent.c_str(), "lep" );
            }
        }
        leg->Draw("SAME");
        for ( auto const& s : sTags ) {
            tex->DrawLatex(tagx, tagy_start, s.c_str());
            tagy_start -= dty;
        }
        tex->DrawLatex(tagx, tagy_start, h3_titles.at(ihist).c_str());
        c->Update();
    }
    if ( !is_empty_pannel ) {
        c->SaveAs(Form("%s/%s_%dx%d.png", output_location.c_str(), pngbase.c_str(), window_dim.first, window_dim.second));
    }
    delete c;
    return;
}





void ProcessCaloWindowHistograms()
{
    std::cout <<"Processing HCAL geo Calowindows" << std::endl;
    float miny = 1e-5, maxy = 1.5e1;
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
        
        // CaloWindowMultiPanel({h2_window_energy_cent_recemc_nxm, h2_window_energy_cent_hcalin_nxm, h2_window_energy_cent_hcalout_nxm}, 
        //                     {"EMCal", "iHCal", "oHCal"},
        //                     calo_window_plots+"/energy_cent_slices", "sub_calo_windows_et_all_cent_slices", 
        //                     Form("E_{T}^{%d #times %d} [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
        //                     Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
        //                     tagx, tagy, dty, x1, x2, y1, y2, 0., 0, miny, maxy);

        // CaloWindowMultiPanel({h2_window_energy_cent_full_nxm},
        //                     {"Total"},
        //                     calo_window_plots+"/energy_cent_slices", "calo_window_et_all_cent_slices", 
        //                     Form("E_{T}^{%d #times %d} [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
        //                     Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
        //                     tagx, tagy, dty, x1, x2, y1_s, y2_s, 0., 0, miny, maxy);

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
        
        // CaloWindowMultiPanel({h2_recemc_energy_minus_avg_energy_cent_nxm, h2_hcalin_energy_minus_avg_energy_cent_nxm, h2_hcalout_energy_minus_avg_energy_cent_nxm}, 
        //                     {"EMCal", "iHCal", "oHCal"},
        //                     calo_window_plots+"/energy_minus_avg_cent_slices", "sub_calo_windows_et-avget_all_cent_slices", 
        //                     Form("E_{T}^{%d #times %d} [GeV] - #LT E_{T}^{%d #times %d} #GT [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second, k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
        //                     Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
        //                     tagx, tagy, dty, x1, x2, y1, y2, 0., 0, miny, maxy, true);

        // CaloWindowMultiPanel({h2_full_energy_minus_avg_energy_cent_nxm},
        //                     {"Total"},
        //                     calo_window_plots+"/energy_minus_avg_cent_slices", "calo_window_et-avget_all_cent_slices", 
        //                     Form("E_{T}^{%d #times %d} [GeV] - #LT E_{T}^{%d #times %d} #GT [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second, k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
        //                     Form("1/N dN/dE_{T}^{%d #times %d} [GeV^{-1}]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
        //                     tagx, tagy, dty, x1, x2, y1_s, y2_s, 0., 0, miny, maxy, true);
        
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
                            Form("f(E_{T}^{%d #times %d}) [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
                            "f(E_{T}^{%d #times %d}) [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second,
                            tagx, tagy, dty, x1, x2, y1, y2, 0., 0, miny, maxy);
        
        delete h3_energy_frac_energy_cent_recemc_nxm;
        delete h3_energy_frac_energy_cent_hcalin_nxm;
        delete h3_energy_frac_energy_cent_hcalout_nxm;
        
        TH2F * h2_window_frac_energy_cent_recemc_nxm = (TH2F*)f->Get(Form("h2_window_frac_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_window_frac_energy_cent_recemc_nxm ) { std::cout << "h2_window_frac_energy_cent_recemc_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH2F * h2_window_frac_energy_cent_hcalin_nxm = (TH2F*)f->Get(Form("h2_window_frac_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_window_frac_energy_cent_hcalin_nxm ) { std::cout << "h2_window_frac_energy_cent_hcalin_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        TH2F * h2_window_frac_energy_cent_hcalout_nxm = (TH2F*)f->Get(Form("h2_window_frac_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
        if ( !h2_window_frac_energy_cent_hcalout_nxm ) { std::cout << "h2_window_frac_energy_cent_hcalout_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
        // CaloWindowMultiPanel({h2_window_frac_energy_cent_recemc_nxm, h2_window_frac_energy_cent_hcalin_nxm, h2_window_frac_energy_cent_hcalout_nxm}, 
        //                     {"EMCal", "iHCal", "oHCal"},
        //                     calo_window_plots+"/frac_energy_cent_slices", "sub_calo_windows_frac_et_all_cent_slices", 
        //                     Form("f(E_{T}^{%d #times %d}) [GeV]", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second),
        //                     "Probability Density",
        //                     tagx, tagy, dty, x1, x2, y1, y2, 0., 0, miny, maxy);
                            
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

    return;
}


void PPG04_Histograms(const std::string & input_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/rootfiles/JAN28/DATA/BASIC/DATA-ProdA_2023-BASIC-023745.root") 
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
    std::string output_file = input_file_base + "-HISTOGRAMS.root";

    std::cout << "Input file: " << input_file << std::endl;
    std::cout << "Output file: " << output_file << std::endl;

    
    ConfigureOutputDirs(input_file_base);
    std::cout << "Global plots: " << global_plots << std::endl;
    std::cout << "Calo plots: " << calo_plots << std::endl;
    std::cout << "Calo window plots: " << calo_window_plots << std::endl;
    std::cout << "Background plots: " << bkgd_plots << std::endl;
    std::cout << "Random cone plots: " << random_cone_plots << std::endl;
   
    // open input file
    f = new TFile(input_file.c_str(), "READ");
    if(!f->IsOpen() || f->IsZombie()){ std::cout << "File " << input_file << " is zombie" << std::endl;  return -1; }

    // get tree
    TTree *t = (TTree*)f->Get("T");
    if(!t){ std::cout << "Tree T not found in " << input_file << std::endl; return -1; }

    // -----------------------------------------------------
    // Global variables
    // -----------------------------------------------------

    int event_id = 0;
    unsigned int random_seed = 0;
    t->SetBranchAddress("event_id", &event_id);
    t->SetBranchAddress("random_seed", &random_seed);

    // nevents histogram
    TH1F * h1_num_events = (TH1F*)f->Get("h1_num_events");
    if( !h1_num_events ){ std::cout << "h1_num_events not found!" << std::endl; return -1; }
    NEVENTS = (int)h1_num_events->GetBinContent(1);
    std::cout << "Number of events: " << NEVENTS << std::endl;

    // MBD
    float mbd_q_N = 0;
    float mbd_q_S = 0;
    float mbd_time_N = 0;
    float mbd_time_S = 0;
    float mbd_npmt_N = 0;
    float mbd_npmt_S = 0;
    t->SetBranchAddress("mbd_q_N", &mbd_q_N);
    t->SetBranchAddress("mbd_q_S", &mbd_q_S);
    t->SetBranchAddress("mbd_time_N", &mbd_time_N);
    t->SetBranchAddress("mbd_time_S", &mbd_time_S);
    t->SetBranchAddress("mbd_npmt_N", &mbd_npmt_N);
    t->SetBranchAddress("mbd_npmt_S", &mbd_npmt_S);

    // Vertex
    float zvtx = 0;
    t->SetBranchAddress("zvtx", &zvtx);

    // Centrality
    int centrality = 0;
    int impact_parameter = 0;
    t->SetBranchAddress("centrality", &centrality);
    t->SetBranchAddress("impact_parameter", &impact_parameter);

    // -----------------------------------------------------
    // Tower background variables
    // -----------------------------------------------------
    float tower_background_v2 = 0;
    float tower_background_psi2 = 0;
    std::vector< float > * tower_background_energy_recemc = 0;
    std::vector< float > * tower_background_energy_hcalin = 0;
    std::vector< float > * tower_background_energy_hcalout = 0;
    t->SetBranchAddress("tower_background_v2", &tower_background_v2);
    t->SetBranchAddress("tower_background_psi2", &tower_background_psi2);
    t->SetBranchAddress("tower_background_energy_recemc", &tower_background_energy_recemc);
    t->SetBranchAddress("tower_background_energy_hcalin", &tower_background_energy_hcalin);
    t->SetBranchAddress("tower_background_energy_hcalout", &tower_background_energy_hcalout);
    // -----------------------------------------------------
    // Rho variables
    // -----------------------------------------------------
    // TowerRho_AREA
    float rho_val_TowerRho_AREA = 0;
    float std_rho_val_TowerRho_AREA = 0;
    t->SetBranchAddress("rho_val_TowerRho_AREA", &rho_val_TowerRho_AREA);
    t->SetBranchAddress("std_rho_val_TowerRho_AREA", &std_rho_val_TowerRho_AREA);

    // TowerRho_MULT
    float rho_val_TowerRho_MULT = 0;
    float std_rho_val_TowerRho_MULT = 0;
    t->SetBranchAddress("rho_val_TowerRho_MULT", &rho_val_TowerRho_MULT);
    t->SetBranchAddress("std_rho_val_TowerRho_MULT", &std_rho_val_TowerRho_MULT);

    // TowerRho_AREA_CEMC
    float rho_val_TowerRho_AREA_CEMC = 0;
    float std_rho_val_TowerRho_AREA_CEMC = 0;
    t->SetBranchAddress("rho_val_TowerRho_AREA_CEMC", &rho_val_TowerRho_AREA_CEMC);
    t->SetBranchAddress("std_rho_val_TowerRho_AREA_CEMC", &std_rho_val_TowerRho_AREA_CEMC);

    // TowerRho_MULT_CEMC
    float rho_val_TowerRho_MULT_CEMC = 0;
    float std_rho_val_TowerRho_MULT_CEMC = 0;
    t->SetBranchAddress("rho_val_TowerRho_MULT_CEMC", &rho_val_TowerRho_MULT_CEMC);
    t->SetBranchAddress("std_rho_val_TowerRho_MULT_CEMC", &std_rho_val_TowerRho_MULT_CEMC);

    // TowerRho_AREA_HCALIN
    float rho_val_TowerRho_AREA_HCALIN = 0;
    float std_rho_val_TowerRho_AREA_HCALIN = 0;
    t->SetBranchAddress("rho_val_TowerRho_AREA_HCALIN", &rho_val_TowerRho_AREA_HCALIN);
    t->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALIN", &std_rho_val_TowerRho_AREA_HCALIN);

    // TowerRho_MULT_HCALIN
    float rho_val_TowerRho_MULT_HCALIN = 0;
    float std_rho_val_TowerRho_MULT_HCALIN = 0;
    t->SetBranchAddress("rho_val_TowerRho_MULT_HCALIN", &rho_val_TowerRho_MULT_HCALIN);
    t->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALIN", &std_rho_val_TowerRho_MULT_HCALIN);

    // TowerRho_AREA_HCALOUT
    float rho_val_TowerRho_AREA_HCALOUT = 0;
    float std_rho_val_TowerRho_AREA_HCALOUT = 0;
    t->SetBranchAddress("rho_val_TowerRho_AREA_HCALOUT", &rho_val_TowerRho_AREA_HCALOUT);
    t->SetBranchAddress("std_rho_val_TowerRho_AREA_HCALOUT", &std_rho_val_TowerRho_AREA_HCALOUT);

    // TowerRho_MULT_HCALOUT
    float rho_val_TowerRho_MULT_HCALOUT = 0;
    float std_rho_val_TowerRho_MULT_HCALOUT = 0;
    t->SetBranchAddress("rho_val_TowerRho_MULT_HCALOUT", &rho_val_TowerRho_MULT_HCALOUT);
    t->SetBranchAddress("std_rho_val_TowerRho_MULT_HCALOUT", &std_rho_val_TowerRho_MULT_HCALOUT);

    // -----------------------------------------------------
    // Calo tower variables
    // -----------------------------------------------------
    // TOWERINFO_CALIB_CEMC
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

    // TOWERINFO_CALIB_HCALIN
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

    // TOWERINFO_CALIB_HCALOUT
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

    // TOWERINFO_CALIB_CEMC_RETOWER
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

    // TOWERINFO_CALIB_CEMC_RETOWER_SUB1
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

    // TOWERINFO_CALIB_HCALIN_SUB1
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

    // TOWERINFO_CALIB_HCALOUT_SUB1
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

    // Calo Tower Histograms
    // TOWERINFO_CALIB_CEMC
    TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC");
    if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC not found!" << std::endl; }
    TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC");
    if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC not found!" << std::endl; }
    TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_CEMC = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC");
    if ( !h2_tower_energy_cent_TOWERINFO_CALIB_CEMC ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_CEMC not found!" << std::endl; }
    TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC");
    if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC not found!" << std::endl; }

    // TOWERINFO_CALIB_HCALIN
    TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN");
    if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN not found!" << std::endl; }
    TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN");
    if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN not found!" << std::endl; }
    TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN");
    if ( !h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN not found!" << std::endl; }
    TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN");
    if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN not found!" << std::endl; }

    // TOWERINFO_CALIB_HCALOUT
    TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT");
    if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT not found!" << std::endl; }
    TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT");
    if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT not found!" << std::endl; }
    TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT");
    if ( !h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT not found!" << std::endl; }
    TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT");
    if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT not found!" << std::endl; }

    // TOWERINFO_CALIB_CEMC_RETOWER
    TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER");
    if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER not found!" << std::endl; }
    TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER");
    if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER not found!" << std::endl; }
    TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER");
    if ( !h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER not found!" << std::endl; }
    TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER");
    if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER not found!" << std::endl; }

    // TOWERINFO_CALIB_CEMC_RETOWER_SUB1
    TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 not found!" << std::endl; }
    TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 not found!" << std::endl; }
    TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    if ( !h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 not found!" << std::endl; }
    TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_CEMC_RETOWER_SUB1 not found!" << std::endl; }

    // TOWERINFO_CALIB_HCALIN_SUB1
    TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1");
    if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 not found!" << std::endl; }
    TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1");
    if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 not found!" << std::endl; }
    TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1");
    if ( !h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_HCALIN_SUB1 not found!" << std::endl; }
    TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN_SUB1 = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN_SUB1");
    if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN_SUB1 ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALIN_SUB1 not found!" << std::endl; }

    // TOWERINFO_CALIB_HCALOUT_SUB1
    TH3F * h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 = (TH3F*)f->Get("h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1");
    if ( !h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 ) { std::cout << " h3_tower_eta_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 not found!" << std::endl; }
    TH3F * h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 = (TH3F*)f->Get("h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1");
    if ( !h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 ) { std::cout << " h3_tower_phi_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 not found!" << std::endl; }
    TH2F * h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 = (TH2F*)f->Get("h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1");
    if ( !h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 ) { std::cout << " h2_tower_energy_cent_TOWERINFO_CALIB_HCALOUT_SUB1 not found!" << std::endl; }
    TH2F * h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT_SUB1 = (TH2F*)f->Get("h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT_SUB1");
    if ( !h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT_SUB1 ) { std::cout << " h2_tower_dead_phi_eta_TOWERINFO_CALIB_HCALOUT_SUB1 not found!" << std::endl; }

    // -----------------------------------------------------
    // Random cones
    // -----------------------------------------------------
    // RandomCones_r04
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

    // RandomCones_r04_Sub1
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

    // -----------------------------------------------------
    // Calo window
    // -----------------------------------------------------
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


    // cemc windows
    float avg_energy_cemc[11];
    float std_energy_cemc[11];
    int num_windows_cemc[11];
    t->SetBranchAddress("avg_energy_cemc", &avg_energy_cemc);
    t->SetBranchAddress("std_energy_cemc", &std_energy_cemc);
    t->SetBranchAddress("num_windows_cemc", &num_windows_cemc);

    ProcessCaloWindowHistograms();

    // std::cout <<"Processing HCAL geo Calowindows" << std::endl;
    // for ( int idim = 0; idim < k_calo_window_dims_hcal_geom.size(); idim++ ) {

    //     TH2F * h2_window_energy_cent_recemc_nxm = (TH2F*)f->Get(Form("h2_window_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
    //     if ( !h2_window_energy_cent_recemc_nxm ) { std::cout << "h2_window_energy_cent_recemc_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
    //     TH2F * h2_window_energy_cent_hcalin_nxm = (TH2F*)f->Get(Form("h2_window_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
    //     if ( !h2_window_energy_cent_hcalin_nxm ) { std::cout << "h2_window_energy_cent_hcalin_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
    //     TH2F * h2_window_energy_cent_hcalout_nxm = (TH2F*)f->Get(Form("h2_window_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
    //     if ( !h2_window_energy_cent_hcalout_nxm ) { std::cout << "h2_window_energy_cent_hcalout_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
    //     TH2F * h2_window_energy_cent_full_nxm = (TH2F*)f->Get(Form("h2_window_energy_cent_full_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
    //     if ( !h2_window_energy_cent_full_nxm ) { std::cout << "h2_window_energy_cent_full_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
    //     CaloWindowEnergySlicesGrid(h2_window_energy_cent_full_nxm, h2_window_energy_cent_recemc_nxm, h2_window_energy_cent_hcalin_nxm, h2_window_energy_cent_hcalout_nxm, calo_window_plots+"/energy_cent_slices");
    //     CaloWindowEnergySlices(h2_window_energy_cent_full_nxm, h2_window_energy_cent_recemc_nxm, h2_window_energy_cent_hcalin_nxm, h2_window_energy_cent_hcalout_nxm, calo_window_plots+"/energy_cent_slices/singles");
    //     delete h2_window_energy_cent_recemc_nxm;
    //     delete h2_window_energy_cent_hcalin_nxm;
    //     delete h2_window_energy_cent_hcalout_nxm;
    //     delete h2_window_energy_cent_full_nxm;
        
        
    //     TH2F * h2_window_frac_energy_cent_recemc_nxm = (TH2F*)f->Get(Form("h2_window_frac_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
    //     if ( !h2_window_frac_energy_cent_recemc_nxm ) { std::cout << "h2_window_frac_energy_cent_recemc_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
    //     TH2F * h2_window_frac_energy_cent_hcalin_nxm = (TH2F*)f->Get(Form("h2_window_frac_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
    //     if ( !h2_window_frac_energy_cent_hcalin_nxm ) { std::cout << "h2_window_frac_energy_cent_hcalin_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
    //     TH2F * h2_window_frac_energy_cent_hcalout_nxm = (TH2F*)f->Get(Form("h2_window_frac_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[idim].first, k_calo_window_dims_hcal_geom[idim].second));
    //     if ( !h2_window_frac_energy_cent_hcalout_nxm ) { std::cout << "h2_window_frac_energy_cent_hcalout_" << k_calo_window_dims_hcal_geom[idim].first << "x" << k_calo_window_dims_hcal_geom[idim].second << " not found!" << std::endl; }
    //     CaloWindowFracSlicesGrid(h2_window_frac_energy_cent_recemc_nxm, h2_window_frac_energy_cent_hcalin_nxm, h2_window_frac_energy_cent_hcalout_nxm, calo_window_plots+"/energy_cent_slices");
    //     CaloWindowFracSlices(h2_window_frac_energy_cent_recemc_nxm, h2_window_frac_energy_cent_hcalin_nxm, h2_window_frac_energy_cent_hcalout_nxm, calo_window_plots+"/energy_cent_slices/singles");
    //     delete h2_window_frac_energy_cent_recemc_nxm;
    //     delete h2_window_frac_energy_cent_hcalin_nxm;
    //     delete h2_window_frac_energy_cent_hcalout_nxm;
    // }

    // std::cout << "Processing CEMC geo Calowindows" << std::endl;
    // for ( int idim = 0; idim < k_calo_window_dims_cemc_geom.size(); idim++ ) {

    //     TH2F * h2_window_energy_cent_cemc_nxm = (TH2F*)f->Get(Form("h2_window_energy_cent_cemc_%dx%d", k_calo_window_dims_cemc_geom[idim].first, k_calo_window_dims_cemc_geom[idim].second));
    //     if ( !h2_window_energy_cent_cemc_nxm ) { 
    //         std::cout << "h2_window_energy_cent_cemc_" << k_calo_window_dims_cemc_geom[idim].first << "x" << k_calo_window_dims_cemc_geom[idim].second << " not found!" << std::endl; 
    //         continue;
    //     }
    //     CaloWindowEnergySlicesGridCEMC(h2_window_energy_cent_cemc_nxm, calo_window_plots+"/energy_cent_slices/cemc_geom");
    //     CaloWindowEnergySlicesCEMC(h2_window_energy_cent_cemc_nxm, calo_window_plots+"/energy_cent_slices/cemc_geom/singles");
    //     delete h2_window_energy_cent_cemc_nxm;
    // }

    // get the number of events
    int nentries = t->GetEntries();
    std::cout << "Number of entries = " << nentries << std::endl;
    // get first entry
    t->GetEntry(0);



    f->Close();

    return 0;
}

