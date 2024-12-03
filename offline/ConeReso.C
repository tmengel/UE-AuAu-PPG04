#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TMath.h>
#include <TString.h>


#include <iostream>
#include <vector>
#include "sPhenixStyle.h"
#include "sPhenixStyle.C"



int ConeReso()
{

  SetsPhenixStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
   
  // in/out file names
  TString data_file = "data_all_rc_nop_histos.root";
  TString hijing_file = "hijing_all_rc_nop_histos.root";

  // open input file
  TFile *f_data = new TFile(data_file, "READ");
  if(!f_data->IsOpen() || f_data->IsZombie()){ std::cout << "File " << data_file << " is zombie" << std::endl;  return -1; }

  // open hijing file
  TFile *f_hijing = new TFile(hijing_file, "READ");
  if(!f_hijing->IsOpen() || f_hijing->IsZombie()){ std::cout << "File " << hijing_file << " is zombie" << std::endl;  return -1; }

  // get histos
  TH1D * h1_cone_residual_pt_basic_area_sigma_data = (TH1D*)f_data->Get("h1_cone_residual_pt_basic_area_sigma");
  TH1D * h1_cone_residual_pt_basic_mult_sigma_data = (TH1D*)f_data->Get("h1_cone_residual_pt_basic_mult_sigma");
  TH1D * h1_cone_residual_pt_basic_sub1_sigma_data = (TH1D*)f_data->Get("h1_cone_residual_pt_basic_sub1_sigma");
  TH1D * h1_cone_residual_pt_rand_eta_phi_area_sigma_data = (TH1D*)f_data->Get("h1_cone_residual_pt_rand_eta_phi_area_sigma");
  TH1D * h1_cone_residual_pt_rand_eta_phi_mult_sigma_data = (TH1D*)f_data->Get("h1_cone_residual_pt_rand_eta_phi_mult_sigma");
  TH1D * h1_cone_residual_pt_rand_eta_phi_sub1_sigma_data = (TH1D*)f_data->Get("h1_cone_residual_pt_rand_eta_phi_sub1_sigma");
  
  h1_cone_residual_pt_basic_area_sigma_data->SetMarkerColor(kBlack);
  h1_cone_residual_pt_basic_area_sigma_data->SetLineColor(kBlack);
  h1_cone_residual_pt_basic_area_sigma_data->SetMarkerStyle(20);
  h1_cone_residual_pt_basic_area_sigma_data->SetMarkerSize(1.5);

  h1_cone_residual_pt_basic_mult_sigma_data->SetMarkerColor(kRed);
  h1_cone_residual_pt_basic_mult_sigma_data->SetLineColor(kRed);
  h1_cone_residual_pt_basic_mult_sigma_data->SetMarkerStyle(21);
  h1_cone_residual_pt_basic_mult_sigma_data->SetMarkerSize(1.5);

  h1_cone_residual_pt_basic_sub1_sigma_data->SetMarkerColor(kBlue-2);
  h1_cone_residual_pt_basic_sub1_sigma_data->SetLineColor(kBlue-2);
  h1_cone_residual_pt_basic_sub1_sigma_data->SetMarkerStyle(22);
  h1_cone_residual_pt_basic_sub1_sigma_data->SetMarkerSize(1.5);

  h1_cone_residual_pt_rand_eta_phi_area_sigma_data->SetMarkerColor(kBlack);
  h1_cone_residual_pt_rand_eta_phi_area_sigma_data->SetLineColor(kBlack);
  h1_cone_residual_pt_rand_eta_phi_area_sigma_data->SetMarkerStyle(24);
  h1_cone_residual_pt_rand_eta_phi_area_sigma_data->SetMarkerSize(1.5);

  h1_cone_residual_pt_rand_eta_phi_mult_sigma_data->SetMarkerColor(kRed);
  h1_cone_residual_pt_rand_eta_phi_mult_sigma_data->SetLineColor(kRed);
  h1_cone_residual_pt_rand_eta_phi_mult_sigma_data->SetMarkerStyle(25);
  h1_cone_residual_pt_rand_eta_phi_mult_sigma_data->SetMarkerSize(1.5);

  h1_cone_residual_pt_rand_eta_phi_sub1_sigma_data->SetMarkerColor(kBlue-2);
  h1_cone_residual_pt_rand_eta_phi_sub1_sigma_data->SetLineColor(kBlue-2);
  h1_cone_residual_pt_rand_eta_phi_sub1_sigma_data->SetMarkerStyle(26);
  h1_cone_residual_pt_rand_eta_phi_sub1_sigma_data->SetMarkerSize(1.5);

  // // get histos
  TH1D * h1_cone_residual_pt_basic_area_sigma_hijing = (TH1D*)f_hijing->Get("h1_cone_residual_pt_basic_area_sigma");
  TH1D * h1_cone_residual_pt_basic_mult_sigma_hijing = (TH1D*)f_hijing->Get("h1_cone_residual_pt_basic_mult_sigma");
  TH1D * h1_cone_residual_pt_basic_sub1_sigma_hijing = (TH1D*)f_hijing->Get("h1_cone_residual_pt_basic_sub1_sigma");
  TH1D * h1_cone_residual_pt_rand_eta_phi_area_sigma_hijing = (TH1D*)f_hijing->Get("h1_cone_residual_pt_rand_eta_phi_area_sigma");
  TH1D * h1_cone_residual_pt_rand_eta_phi_mult_sigma_hijing = (TH1D*)f_hijing->Get("h1_cone_residual_pt_rand_eta_phi_mult_sigma");
  TH1D * h1_cone_residual_pt_rand_eta_phi_sub1_sigma_hijing = (TH1D*)f_hijing->Get("h1_cone_residual_pt_rand_eta_phi_sub1_sigma");

  h1_cone_residual_pt_rand_eta_phi_area_sigma_hijing->SetMarkerColor(kBlack);
  h1_cone_residual_pt_rand_eta_phi_area_sigma_hijing->SetLineColor(kBlack);
  h1_cone_residual_pt_rand_eta_phi_area_sigma_hijing->SetMarkerStyle(24);
  h1_cone_residual_pt_rand_eta_phi_area_sigma_hijing->SetMarkerSize(1.5);

  h1_cone_residual_pt_rand_eta_phi_mult_sigma_hijing->SetMarkerColor(kRed);
  h1_cone_residual_pt_rand_eta_phi_mult_sigma_hijing->SetLineColor(kRed);
  h1_cone_residual_pt_rand_eta_phi_mult_sigma_hijing->SetMarkerStyle(25);
  h1_cone_residual_pt_rand_eta_phi_mult_sigma_hijing->SetMarkerSize(1.5);

  h1_cone_residual_pt_rand_eta_phi_sub1_sigma_hijing->SetMarkerColor(kBlue-2);
  h1_cone_residual_pt_rand_eta_phi_sub1_sigma_hijing->SetLineColor(kBlue-2);
  h1_cone_residual_pt_rand_eta_phi_sub1_sigma_hijing->SetMarkerStyle(26);
  h1_cone_residual_pt_rand_eta_phi_sub1_sigma_hijing->SetMarkerSize(1.5);

  // use open symbols
  h1_cone_residual_pt_basic_area_sigma_hijing->SetMarkerColor(kBlack);
  h1_cone_residual_pt_basic_area_sigma_hijing->SetLineColor(kBlack);
  h1_cone_residual_pt_basic_area_sigma_hijing->SetMarkerStyle(24);
  h1_cone_residual_pt_basic_area_sigma_hijing->SetMarkerSize(1.5);

  h1_cone_residual_pt_basic_mult_sigma_hijing->SetMarkerColor(kRed);
  h1_cone_residual_pt_basic_mult_sigma_hijing->SetLineColor(kRed);
  h1_cone_residual_pt_basic_mult_sigma_hijing->SetMarkerStyle(25);
  h1_cone_residual_pt_basic_mult_sigma_hijing->SetMarkerSize(1.5);

  h1_cone_residual_pt_basic_sub1_sigma_hijing->SetMarkerColor(kBlue-2);
  h1_cone_residual_pt_basic_sub1_sigma_hijing->SetLineColor(kBlue-2);
  h1_cone_residual_pt_basic_sub1_sigma_hijing->SetMarkerStyle(26);
  h1_cone_residual_pt_basic_sub1_sigma_hijing->SetMarkerSize(1.5);

  double maxX = std::max(h1_cone_residual_pt_basic_area_sigma_data->GetMaximum(), h1_cone_residual_pt_basic_area_sigma_hijing->GetMaximum());
  h1_cone_residual_pt_basic_area_sigma_data->GetYaxis()->SetRangeUser(0, 1.1*maxX);
  h1_cone_residual_pt_basic_area_sigma_data->GetXaxis()->SetRangeUser(0, 2700);

  TCanvas *c = new TCanvas("c", "c", 800, 600);
  gStyle->SetOptStat(0);
  gPad->SetLogx();
  // gPad->SetLogy();
  // set x axis range
  TLegend *leg = new TLegend(0.2, 0.7, 0.7, 0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  // two columns
  leg->SetNColumns(2);

  
  std::vector< TH1D* > h1s = {h1_cone_residual_pt_basic_area_sigma_data, 
                              h1_cone_residual_pt_basic_area_sigma_hijing,
                              h1_cone_residual_pt_basic_mult_sigma_data,
                              h1_cone_residual_pt_basic_mult_sigma_hijing,
                              h1_cone_residual_pt_basic_sub1_sigma_data,
                              h1_cone_residual_pt_basic_sub1_sigma_hijing};
  TGraphErrors * gr[h1s.size()];
  std::vector<std::string> titles = {"#rho_{A}*A", "#rho_{M}*N", "Hijing #rho_{A}*A", "Hijing #rho_{M}*N", "Iterative", "Hijing Iterative"};
  int colors[] = {kBlack, kBlack, kRed, kRed, kBlue-2, kBlue-2};
  int markers[] = {20, 24, 21, 25, 22, 26};
  std::pair<double, double> xbounds = {0.5e2, 0.3e4};
  double maxy = 0, miny = 1e9;
  for (auto h1 : h1s)
  {
    maxy = std::max(maxy, h1->GetMaximum());
    miny = std::min(miny, h1->GetMinimum());
  }
  std::pair<double, double> ybounds = {0, 1.1*maxy};
  std::string xlab = "#Sigma E_{T}^{Calo} [GeV]";
  std::string ylab = "#sigma(#delta E_{T})";
  double marker_size = 1.5;

  for ( unsigned int i = 0; i < h1s.size(); i++)
  {
    gr[i] = new TGraphErrors(h1s[i]->GetNbinsX());
    for (int ibin = 1; ibin <= h1s[i]->GetNbinsX(); ibin++)
    {
      gr[i]->SetPoint(ibin-1, h1s[i]->GetBinCenter(ibin), h1s[i]->GetBinContent(ibin));
      gr[i]->SetPointError(ibin-1, 0, h1s[i]->GetBinError(ibin));
    }
    gr[i]->SetMarkerColor(colors[i]);
    gr[i]->SetLineColor(colors[i]);
    gr[i]->SetLineWidth(2);
    gr[i]->SetMarkerStyle(markers[i]);
    gr[i]->SetMarkerSize(marker_size);
    gr[i]->GetXaxis()->SetRangeUser(xbounds.first, xbounds.second);
    gr[i]->GetYaxis()->SetRangeUser(ybounds.first, ybounds.second);
    gr[i]->GetXaxis()->SetTitle(xlab.c_str());
    gr[i]->GetYaxis()->SetTitle(ylab.c_str());
    if (i == 0)
    {
      gr[i]->Draw("AP");
    }
    else
    {
      gr[i]->Draw("P SAME");
    }
    leg->AddEntry(gr[i], titles[i].c_str(), "PL");
  }
  leg->Draw("same");


//   h1_cone_residual_pt_basic_mult_sigma_hijing->GetYaxis()->SetRangeUser(0, 9);
//   h1_cone_residual_pt_basic_mult_sigma_hijing->GetXaxis()->SetRangeUser(1e-2, 2000);
//   h1_cone_residual_pt_basic_mult_sigma_hijing->GetXaxis()->SetLimits(1e-2, 2000);  
//   h1_cone_residual_pt_basic_mult_sigma_hijing->Draw("APL");
//   h1_cone_residual_pt_basic_mult_sigma_hijing->GetXaxis()->SetTitle("#Sigma E_{T}^{Calo} [GeV]");
//   h1_cone_residual_pt_basic_mult_sigma_hijing->GetYaxis()->SetTitle("#sigma(#delta E_{T})");

//   h1_cone_residual_pt_basic_mult_sigma_data->Draw("PL SAME");

//   h1_cone_residual_pt_basic_area_sigma_hijing->Draw("PL SAME");
//     h1_cone_residual_pt_basic_area_sigma_data->Draw("PL SAME");
//   h1_cone_residual_pt_basic_area_sigma_hijing->GetXaxis()->SetRangeUser(1, 2000);
//   h1_cone_residual_pt_basic_sub1_sigma_hijing->Draw("PL SAME");
//   h1_cone_residual_pt_basic_sub1_sigma_hijing->GetXaxis()->SetRangeUser(1, 2000);

  

//   h1_cone_residual_pt_basic_area_sigma_data->GetXaxis()->SetRangeUser(1, 2000);
//   h1_cone_residual_pt_basic_sub1_sigma_data->Draw("PL SAME");
//   h1_cone_residual_pt_basic_sub1_sigma_data->GetXaxis()->SetRangeUser(1, 2000);
//   h1_cone_residual_pt_basic_mult_sigma_data->GetXaxis()->SetRangeUser(1, 2000);



//   leg->AddEntry(h1_cone_residual_pt_basic_area_sigma_data, "#rho_{A}*A", "PL");
//     leg->AddEntry(h1_cone_residual_pt_basic_area_sigma_hijing, "Hijing #rho_{A}*A", "PL");
//   leg->AddEntry(h1_cone_residual_pt_basic_mult_sigma_data, "#rho_{M}*N", "PL");
//     leg->AddEntry(h1_cone_residual_pt_basic_mult_sigma_hijing, "Hijing #rho_{M}*N", "PL");
//   leg->AddEntry(h1_cone_residual_pt_basic_sub1_sigma_data, "Iterative", "PL");
//   // add one open symbol for hijing

// ;
//   leg->AddEntry(h1_cone_residual_pt_basic_sub1_sigma_hijing, "Hijing Iterative", "PL");
//   leg->Draw("same");


  c->SaveAs("cone_residual_pt_basic_sigma.png");

    


    TLegend *leg2 = new TLegend(0.6, 0.6, 0.8, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);



    TH1D * h_et_sum_calos_data = (TH1D*)f_data->Get("h_et_sum_calos");
    TH1D * h_et_sum_calos_hijing = (TH1D*)f_hijing->Get("h_et_sum_calos");
    gPad->SetLogy();
    h_et_sum_calos_data->SetLineColor(kBlack);
    h_et_sum_calos_data->SetMarkerColor(kBlack);
    h_et_sum_calos_data->SetMarkerStyle(20);
    h_et_sum_calos_data->SetMarkerSize(1.5);
    h_et_sum_calos_data->Scale(1./h_et_sum_calos_data->Integral());
    h_et_sum_calos_hijing->SetLineColor(kRed);
    h_et_sum_calos_hijing->SetMarkerColor(kRed);
    h_et_sum_calos_hijing->SetMarkerStyle(21);
    h_et_sum_calos_hijing->SetMarkerSize(1.5);
    h_et_sum_calos_hijing->Scale(1./h_et_sum_calos_hijing->Integral());

    h_et_sum_calos_hijing->Draw("PL");
    h_et_sum_calos_hijing->GetXaxis()->SetTitle("#Sigma E_{T} [GeV]");
    h_et_sum_calos_hijing->GetYaxis()->SetTitle("Prob.");
    h_et_sum_calos_data->Draw("PL SAME");
    leg2->Clear();
    leg2->AddEntry(h_et_sum_calos_data, "Data", "PL");
    leg2->AddEntry(h_et_sum_calos_hijing, "Hijing", "PL");
    leg2->Draw("same");

    c->SaveAs("et_sum_calos.png");




    // h_et_sum_calos_data->GetYaxis()->SetRangeUser(0, 1.1*h_et_sum_calos_data->GetMaximum());



    f_data->Close();
    f_hijing->Close();


    return 0;

}