#include <string>
#include <iostream>
#include <sPhenixStyle.C>

int plot_sgimas(){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();
    SetsPhenixStyle();

    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kRainBow);


    TCanvas * c;
    TLegend * leg;
    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    double tx=0.19;
    double ty_start=0.85;
    double tx_f = 0.55;
    double ty_f_start = 0.35;

    const std::string sPHENIX_Tag = "#it{#bf{sPHENIX}} Internal";
    const std::string DataType_Tag = "Au+Au #sqrt{s_{NN}} = 200 GeV";


    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};

    c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    // set reverse x-axis

    leg = new TLegend(0.4,0.2,0.7,0.3);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    TFile * fin0to170 = new TFile("sigma_0to170.root", "READ"); 
    TFile * fin10to100 = new TFile("sigma_10to100.root", "READ");
    if (!fin0to170 || !fin0to170->IsOpen() || !fin10to100 || !fin10to100->IsOpen()) {
        std::cerr << "Error: Could not open input files!" << std::endl;
        return 1;
    }

    TGraphErrors * g_sigma_fit = (TGraphErrors*)fin0to170->Get("g_sigma_fit");
    if (!g_sigma_fit) {
        std::cerr << "Error: g_sigma_fit not found in sigma_0to170.root!" << std::endl;
        fin0to170->Close();
        return 1;
    }
    g_sigma_fit->SetName("g_sigma_fit_0to170");
    TGraphErrors * g_sigma_fit_10to100 = (TGraphErrors*)fin10to100->Get("g_sigma_fit");
    if (!g_sigma_fit_10to100) {
        std::cerr << "Error: g_sigma_fit not found in sigma_10to100.root!" << std::endl;
        fin10to100->Close();
        return 1;
    }
    g_sigma_fit_10to100->SetName("g_sigma_fit_10to100");


    g_sigma_fit->GetXaxis()->SetNdivisions(505);
    g_sigma_fit->GetYaxis()->SetNdivisions(505);
    g_sigma_fit->GetXaxis()->SetTitle("Centrality [%]");
    g_sigma_fit->GetYaxis()->SetTitle("k_{fit}");
    g_sigma_fit->GetXaxis()->SetRangeUser(-3, 85);
    g_sigma_fit->GetYaxis()->SetRangeUser(0.5, 0.6);
    g_sigma_fit->SetMarkerStyle(20);
    g_sigma_fit->SetMarkerSize(1.5);
    g_sigma_fit->SetLineColor(kBlack);
    g_sigma_fit->SetMarkerColor(kBlack);
    g_sigma_fit->Draw("AP");

    leg->AddEntry(g_sigma_fit, "0-170 Area fit range", "pe");

    g_sigma_fit_10to100->GetXaxis()->SetNdivisions(505);
    g_sigma_fit_10to100->GetYaxis()->SetNdivisions(505);
    g_sigma_fit_10to100->GetXaxis()->SetTitle("Centrality [%]");
    g_sigma_fit_10to100->GetYaxis()->SetTitle("k_{fit}");
    g_sigma_fit_10to100->GetXaxis()->SetRangeUser(-3, 85);
    g_sigma_fit_10to100->GetYaxis()->SetRangeUser(0.5, 0.6);
    g_sigma_fit_10to100->SetMarkerStyle(21);
    g_sigma_fit_10to100->SetMarkerSize(1.5);
    g_sigma_fit_10to100->SetLineColor(kRed);
    g_sigma_fit_10to100->SetMarkerColor(kRed);
    g_sigma_fit_10to100->Draw("P SAME");

    leg->AddEntry(g_sigma_fit_10to100, "10-100 Area fit range", "pe");

    tx = 0.19;
    ty_start = 0.85;
    for ( auto tag : tags ) {
        tex->DrawLatex(0.19, ty_start, tag.c_str());
        ty_start -= 0.06;
    }

    leg->Draw();
    c->SaveAs("sigmas_plot_comp.png");
    // delete c;

    // c = new TCanvas("c_ratio", "c_ratio", 800, 200);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetRightMargin(0.15);
    // gPad->SetBottomMargin(0.15);
    // gPad->SetTopMargin(0.05);
    c->Clear();
    
    TGraphErrors * g_sigma_fit_ratio = new TGraphErrors(*g_sigma_fit_10to100);
    g_sigma_fit_ratio->SetName("g_sigma_fit_ratio");
    for (int i = 0; i < g_sigma_fit_ratio->GetN(); ++i) {
        double x, y;
        g_sigma_fit->GetPoint(i, x, y);
        double y_10to100 = g_sigma_fit_10to100->Eval(x);
        if (y_10to100 != 0) {
            double ratio = y / y_10to100;
            g_sigma_fit_ratio->SetPoint(i, x, ratio);
            double ex = g_sigma_fit->GetErrorX(i);
            double ey = g_sigma_fit->GetErrorY(i) / y_10to100;
            g_sigma_fit_ratio->SetPointError(i, ex, ey);
        } else {
            g_sigma_fit_ratio->SetPoint(i, x, 0); // avoid division by zero
            g_sigma_fit_ratio->SetPointError(i, 0, 0);
        }
    }
    g_sigma_fit_ratio->GetXaxis()->SetNdivisions(505);
    g_sigma_fit_ratio->GetYaxis()->SetNdivisions(505);
    g_sigma_fit_ratio->GetXaxis()->SetTitle("Centrality [%]");
    g_sigma_fit_ratio->GetYaxis()->SetTitle("k_{fit} Ratio (0-170 / 10-100)");
    g_sigma_fit_ratio->GetXaxis()->SetRangeUser(-3, 85);
    g_sigma_fit_ratio->GetYaxis()->SetRangeUser(0.98, 1.01);
    g_sigma_fit_ratio->SetMarkerStyle(20);
    g_sigma_fit_ratio->SetMarkerSize(1.5);
    g_sigma_fit_ratio->SetLineColor(kBlue);
    g_sigma_fit_ratio->SetMarkerColor(kBlue);
    g_sigma_fit_ratio->Draw("AP");
    leg->Clear();
    ty_f_start = 0.75; // Reset ty_f_start for the ratio plot
    leg->AddEntry(g_sigma_fit_ratio, "k_{fit} Ratio (0-170 / 10-100)", "pe");
    // tex->DrawLatex(0.19, ty_f_start, sPHENIX_Tag.c_str());
    ty_f_start -= 0.06;
    // tex->DrawLatex(0.19, ty_f_start, DataType_Tag.c_str());
    ty_f_start -= 0.06;
    leg->Draw();
    c->SaveAs("sigmas_plot_ratio.png");


    // Clean up
    fin0to170->Close();
    fin10to100->Close();


    return 0;
}