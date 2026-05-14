#include <string>


#include <sPhenixStyle.C>

const std::string sPHENIX_Tag = "#it{#bf{sPHENIX}} Internal";
const std::string DataType_Tag = "Au+Au 200 GeV";


void WindowFitRatio()
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

    
    std::string old_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/offline/plots/MAR17_OLD/window/sigma_fits.root";
    std::string new_file = "/sphenix/user/tmengel/UE-AuAu-PPG04/offline/plots/MAR17_NEW/window/sigma_fits.root";

    TFile * fold = new TFile(old_file.c_str(), "READ");
    TFile * fnew = new TFile(new_file.c_str(), "READ");
    if ( !fold || !fnew ) {
        std::cout << "Failed to open file" << std::endl;
        return;
    }



    TGraphErrors * g_old = (TGraphErrors*)fold->Get("g_sigma_fit");
    g_old->SetName("g_sigma_fit_old");
    TGraphErrors * g_new = (TGraphErrors*)fnew->Get("g_sigma_fit");
    g_new->SetName("g_sigma_fit_new");

    int N = g_old->GetN();
    TGraphErrors * g_ratio = new TGraphErrors(N);
    g_ratio->SetName("g_sigma_fit_ratio");
    for ( int i = 0; i < N; i++ ) {
        double x, y;
        g_old->GetPoint(i, x, y);
        double yerr = g_old->GetErrorY(i);
        double ynew, ynewerr;
        g_new->GetPoint(i, x, ynew);
        ynewerr = g_new->GetErrorY(i);
        double ratio = ynew / y;
        double ratio_err = ratio * sqrt( pow(yerr/y, 2) + pow(ynewerr/ynew, 2) );
        g_ratio->SetPoint(i, x, ratio);
        g_ratio->SetPointError(i, 0, ratio_err);
    }



    TCanvas * c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
    TLatex * tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    g_ratio->GetXaxis()->SetTitle("Centrality [%]");
    g_ratio->GetYaxis()->SetTitle("New k^{fit}/Old k^{fit}");

    double tx=0.19;
    double ty_start=0.85;
    double tx_f = 0.55;
    double ty_f_start = 0.35;
    std::vector<std::string> tags = {sPHENIX_Tag, DataType_Tag};

    g_ratio->GetXaxis()->SetNdivisions(505);
    g_ratio->GetYaxis()->SetNdivisions(505);
    g_ratio->GetXaxis()->SetRangeUser(0, 81);
    g_ratio->GetYaxis()->SetRangeUser(0.95, 1.05);
    g_ratio->SetMarkerStyle(20);
    g_ratio->SetMarkerSize(1.5);
    g_ratio->SetLineColor(kBlack);
    g_ratio->SetMarkerColor(kBlack);
    g_ratio->Draw("AP");
    tx = 0.19;
    ty_start = 0.85;
    for ( auto tag : tags ) {
        tex->DrawLatex(0.19, ty_start, tag.c_str());
        ty_start -= 0.05;
    }
    c->SaveAs("fit_ratio.png");


    delete c;
    fold->Close();
    fnew->Close();



    return ;

}






