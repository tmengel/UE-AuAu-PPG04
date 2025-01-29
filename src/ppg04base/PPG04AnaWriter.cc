#include "PPG04AnaWriter.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <mbd/MbdOut.h>
#include <mbd/MbdOutV1.h>
#include <mbd/MbdOutV2.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexMapv1.h>

#include <centrality/CentralityInfo.h>
#include <centrality/CentralityInfov1.h>
#include <centrality/CentralityInfov2.h>

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h> 
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoDefs.h>

#include <jetbackground/TowerRho.h>
#include <jetbackground/TowerRhov1.h>
#include <jetbackground/TowerBackground.h>
#include <jetbackground/TowerBackgroundv1.h>

#include <underlyingevent/RandomCone.h>
#include <underlyingevent/RandomConev1.h>
#include <underlyingevent/CaloWindowMap.h>
#include <underlyingevent/CaloWindowMapv1.h>

#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include <cstdlib>

PPG04AnaWriter::PPG04AnaWriter( const std::string &outputfile )
  : SubsysReco("PPG04AnaWriter")
  , m_output_filename(outputfile)
{}


int PPG04AnaWriter::Init( PHCompositeNode * /*topNode*/ )
{

  // create output file
  PHTFileServer::get().open(m_output_filename, "RECREATE");

  // Histograms
  m_event_id = -1;
  m_h1_num_events = new TH1F("h1_num_events", "Number of events", 1, 0, 1);
  m_h1_num_events->GetXaxis()->SetBinLabel(1, "Number of events");

  if ( Verbosity () > 0 ){
    std::cout << "PPG04AnaWriter::Init - opening file " << m_output_filename << std::endl;
  } 

  // create tree
  m_tree = new TTree("T", "PPG04 Analysis Tree");
  m_tree->Branch("event_id", &m_event_id, "event_id/I");
  m_tree->Branch("random_seed", &m_random_seed, "random_seed/i");

  // MBD
  if ( m_do_mbd ) {
    m_tree->Branch("mbd_q_N", &m_mbd_q_N, "mbd_q_N/F");
    m_tree->Branch("mbd_q_S", &m_mbd_q_S, "mbd_q_S/F");
    m_tree->Branch("mbd_time_N", &m_mbd_time_N, "mbd_time_N/F");
    m_tree->Branch("mbd_time_S", &m_mbd_time_S, "mbd_time_S/F");
    m_tree->Branch("mbd_npmt_N", &m_mbd_npmt_N, "mbd_npmt_N/F");
    m_tree->Branch("mbd_npmt_S", &m_mbd_npmt_S, "mbd_npmt_S/F");

    if ( Verbosity() > 0 ) {
      std::cout << "PPG04AnaWriter::Init - MBD node: " << m_mbd_node << std::endl;
    }
  }

  // Zvtx
  if ( m_do_zvrtx ) {
    m_tree->Branch("zvtx", &m_zvtx, "zvtx/F");
    if ( Verbosity() > 0 ) {
      std::cout << "PPG04AnaWriter::Init - Zvtx node: " << m_zvrtx_node << std::endl;
    }
  }
 
  // Centrality
  if ( m_do_cent ) {
    m_tree->Branch("centrality", &m_centrality, "centrality/I");
    m_tree->Branch("impact_parameter", &m_impact_parameter, "impact_parameter/I");
    if ( Verbosity() > 0 ) {
      std::cout << "PPG04AnaWriter::Init - Centrality node: " << m_cent_node << std::endl;
    }
  }

  // Tower background
  if ( m_do_tower_background ) {
    m_tree->Branch("tower_background_v2", &m_tower_background_v2, "tower_background_v2/F");
    m_tree->Branch("tower_background_psi2", &m_tower_background_psi2, "tower_background_psi2/F");
    m_tree->Branch("tower_background_energy_recemc", &m_tower_background_energy_recemc);
    m_tree->Branch("tower_background_energy_hcalin", &m_tower_background_energy_hcalin);
    m_tree->Branch("tower_background_energy_hcalout", &m_tower_background_energy_hcalout);

    if ( Verbosity() > 0 ) {
      std::cout << "PPG04AnaWriter::Init - Tower background node: " << m_tower_background_node << std::endl;
    }
  }  

  // Rhos
  if ( m_do_rho ) {
    for ( unsigned int i = 0; i < m_rho_nodes.size(); ++i ) {
      m_tree->Branch(Form("rho_val_%s", m_rho_nodes[i].c_str()), &m_rho_val[i], Form("rho_val_%s/F", m_rho_nodes[i].c_str()));
      m_tree->Branch(Form("std_rho_val_%s", m_rho_nodes[i].c_str()), &m_std_rho_val[i], Form("std_rho_val_%s/F", m_rho_nodes[i].c_str()));
    }

    if ( Verbosity() > 0 ) {
      std::cout << "PPG04AnaWriter::Init - Rho nodes: ";
      for ( unsigned int i = 0; i < m_rho_nodes.size(); ++i ) {
        std::cout << m_rho_nodes[i] << " ";
      }
      std::cout << std::endl;
    }
  }

  // Calo nodes
  if ( m_do_calo ) {
    // check that its only central calos
    for ( unsigned int i = 0; i < m_calo_nodes.size(); ++i ) {
      if ( ! ( m_calo_nodes[i].find("CEMC") != std::string::npos || m_calo_nodes[i].find("HCALIN") != std::string::npos || m_calo_nodes[i].find("HCALOUT") != std::string::npos ) ) {
        std::cout << "PPG04AnaWriter::Init - Error - Calo node " << m_calo_nodes[i] << " is not a central calo node. Ignoring" << std::endl;
        m_calo_nodes.erase(m_calo_nodes.begin() + i);
      }
    }
    
  
    for ( unsigned int i = 0; i < m_calo_nodes.size(); ++i ) {
      m_tree->Branch(Form("tower_frac_fired_%s", m_calo_nodes[i].c_str()), &m_tower_frac_fired[i], Form("tower_frac_fired_%s/F", m_calo_nodes[i].c_str()));
      m_tree->Branch(Form("tower_frac_dead_%s", m_calo_nodes[i].c_str()), &m_tower_frac_dead[i], Form("tower_frac_dead_%s/F", m_calo_nodes[i].c_str()));
      m_tree->Branch(Form("tower_avg_energy_%s", m_calo_nodes[i].c_str()), &m_tower_avg_energy[i], Form("tower_avg_energy_%s/F", m_calo_nodes[i].c_str()));
      m_tree->Branch(Form("tower_std_energy_%s", m_calo_nodes[i].c_str()), &m_tower_std_energy[i], Form("tower_std_energy_%s/F", m_calo_nodes[i].c_str()));
      m_tree->Branch(Form("tower_sum_energy_%s", m_calo_nodes[i].c_str()), &m_tower_sum_energy[i], Form("tower_sum_energy_%s/F", m_calo_nodes[i].c_str()));
    }

    m_h3_tower_eta_energy_cent.clear();
    m_h3_tower_phi_energy_cent.clear();
    m_h2_tower_energy_cent.clear();
    m_h2_tower_dead_phi_eta.clear();
    m_h2_tower_energy_phi_eta.clear();
    int N_calo_e_bins = 170;
    float calo_e_min = -2.0;
    float calo_e_max = 15.0;
    for ( unsigned int i = 0; i < m_calo_nodes.size(); ++i ) {
      int n_phi = 64;
      int n_eta = 24;
      // find "CEMC" but not "RETOWER"
      if ( m_calo_nodes[i].find("CEMC") != std::string::npos && m_calo_nodes[i].find("RETOWER") == std::string::npos ) {
        n_phi = 256;
        n_eta = 96;
      }

      TH3F * h3_eta_energy_cent  = new TH3F(Form("h3_tower_eta_energy_cent_%s", m_calo_nodes[i].c_str()), Form("Tower eta vs energy vs centrality %s", m_calo_nodes[i].c_str()), n_eta, -0.5, n_eta-0.5, N_calo_e_bins, calo_e_min, calo_e_max, 100, 0, 100);
      h3_eta_energy_cent->GetXaxis()->SetTitle("ieta");
      h3_eta_energy_cent->GetYaxis()->SetTitle("Energy (GeV)");
      h3_eta_energy_cent->GetZaxis()->SetTitle("Centrality (%)");

      TH3F * h3_phi_energy_cent  = new TH3F(Form("h3_tower_phi_energy_cent_%s", m_calo_nodes[i].c_str()), Form("Tower phi vs energy vs centrality %s", m_calo_nodes[i].c_str()), n_phi, -0.5, n_phi-0.5, N_calo_e_bins, calo_e_min, calo_e_max, 100, 0, 100);
      h3_phi_energy_cent->GetXaxis()->SetTitle("iphi");
      h3_phi_energy_cent->GetYaxis()->SetTitle("Energy (GeV)");
      h3_phi_energy_cent->GetZaxis()->SetTitle("Centrality (%)");

      TH2F * h2_energy_cent = new TH2F(Form("h2_tower_energy_cent_%s", m_calo_nodes[i].c_str()), Form("Tower energy vs centrality %s", m_calo_nodes[i].c_str()), N_calo_e_bins, calo_e_min, calo_e_max, 100, 0, 100);
      h2_energy_cent->GetXaxis()->SetTitle("Energy (GeV)");
      h2_energy_cent->GetYaxis()->SetTitle("Centrality (%)");

      TH2F * h2_dead_phi_eta = new TH2F(Form("h2_tower_dead_phi_eta_%s", m_calo_nodes[i].c_str()), Form("Dead tower phi vs eta %s", m_calo_nodes[i].c_str()), n_phi, -0.5, n_phi-0.5, n_eta, -0.5, n_eta-0.5);
      h2_dead_phi_eta->GetXaxis()->SetTitle("iphi");
      h2_dead_phi_eta->GetYaxis()->SetTitle("ieta");

      TH2F * h2_energy_phi_eta = new TH2F(Form("h2_tower_energy_phi_eta_%s", m_calo_nodes[i].c_str()), Form("Tower energy phi vs eta %s", m_calo_nodes[i].c_str()), n_phi, -0.5, n_phi-0.5, n_eta, -0.5, n_eta-0.5);
      h2_energy_phi_eta->GetXaxis()->SetTitle("iphi");
      h2_energy_phi_eta->GetYaxis()->SetTitle("ieta");


      m_h3_tower_eta_energy_cent.push_back(h3_eta_energy_cent);
      m_h3_tower_phi_energy_cent.push_back(h3_phi_energy_cent);
      m_h2_tower_energy_cent.push_back(h2_energy_cent);
      m_h2_tower_dead_phi_eta.push_back(h2_dead_phi_eta);  
      m_h2_tower_energy_phi_eta.push_back(h2_energy_phi_eta);
    }

    if ( Verbosity() > 0 ) {
      std::cout << "PPG04AnaWriter::Init - Registered Calo nodes: ";
      for ( unsigned int i = 0; i < m_calo_nodes.size(); ++i ) {
        std::cout << m_calo_nodes[i] << " ";
      }
      std::cout << std::endl;
    }

  }

  // Random cone
  if ( m_do_randomcone ) { 
    for ( unsigned int i = 0; i < m_randomcone_nodes.size(); ++i ) {
      m_tree->Branch(Form("random_cone_R_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_R[i], Form("random_cone_R_%s/F", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_eta_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_eta[i], Form("random_cone_eta_%s/F", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_phi_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_phi[i], Form("random_cone_phi_%s/F", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_energy_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_energy[i], Form("random_cone_energy_%s/F", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_energy_cemc_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_energy_cemc[i], Form("random_cone_energy_cemc_%s/F", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_energy_hcalin_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_energy_hcalin[i], Form("random_cone_energy_hcalin_%s/F", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_energy_hcalout_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_energy_hcalout[i], Form("random_cone_energy_hcalout_%s/F", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_num_towers_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_num_towers[i], Form("random_cone_num_towers_%s/I", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_num_towers_cemc_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_num_towers_cemc[i], Form("random_cone_num_towers_cemc_%s/I", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_num_towers_hcalin_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_num_towers_hcalin[i], Form("random_cone_num_towers_hcalin_%s/I", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_num_towers_hcalout_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_num_towers_hcalout[i], Form("random_cone_num_towers_hcalout_%s/I", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_num_masked_towers_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_num_masked_towers[i], Form("random_cone_num_masked_towers_%s/I", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_num_masked_towers_cemc_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_num_masked_towers_cemc[i], Form("random_cone_num_masked_towers_cemc_%s/I", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_num_masked_towers_hcalin_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_num_masked_towers_hcalin[i], Form("random_cone_num_masked_towers_hcalin_%s/I", m_randomcone_nodes[i].c_str()));
      m_tree->Branch(Form("random_cone_num_masked_towers_hcalout_%s", m_randomcone_nodes[i].c_str()), &m_random_cone_num_masked_towers_hcalout[i], Form("random_cone_num_masked_towers_hcalout_%s/I", m_randomcone_nodes[i].c_str()));
    }
    if ( Verbosity() > 0 ) {
      std::cout << "PPG04AnaWriter::Init - Random cone nodes: ";
      for ( unsigned int i = 0; i < m_randomcone_nodes.size(); ++i ) {
        std::cout << m_randomcone_nodes[i] << " ";
      }
      std::cout << std::endl;
    }
  }
  
  // Calo window
  if ( m_do_calo_window ) {

    m_tree->Branch("num_window_dims", &m_max_window_vector_size, "num_window_dims/i");
    m_tree->Branch("num_windows_full", &m_num_windows_full, "num_windows_full[num_window_dims]/I");
    m_tree->Branch("avg_energy_full", &m_avg_energy_full, "avg_energy_full[num_window_dims]/F");
    m_tree->Branch("std_energy_full", &m_std_energy_full, "std_energy_full[num_window_dims]/F");
    m_tree->Branch("avg_frac_energy_recemc_full", &m_avg_frac_energy_recemc_full, "avg_frac_energy_recemc_full[num_window_dims]/F");
    m_tree->Branch("std_frac_energy_recemc_full", &m_std_frac_energy_recemc_full, "std_frac_energy_recemc_full[num_window_dims]/F");
    m_tree->Branch("avg_frac_energy_hcalin_full", &m_avg_frac_energy_hcalin_full, "avg_frac_energy_hcalin_full[num_window_dims]/F");
    m_tree->Branch("std_frac_energy_hcalin_full", &m_std_frac_energy_hcalin_full, "std_frac_energy_hcalin_full[num_window_dims]/F");
    m_tree->Branch("avg_frac_energy_hcalout_full", &m_avg_frac_energy_hcalout_full, "avg_frac_energy_hcalout_full[num_window_dims]/F");
    m_tree->Branch("std_frac_energy_hcalout_full", &m_std_frac_energy_hcalout_full, "std_frac_energy_hcalout_full[num_window_dims]/F");
    m_tree->Branch("num_windows_recemc", &m_num_windows_recemc, "num_windows_recemc[num_window_dims]/I");
    m_tree->Branch("avg_energy_recemc", &m_avg_energy_recemc, "avg_energy_recemc[num_window_dims]/F");
    m_tree->Branch("std_energy_recemc", &m_std_energy_recemc, "std_energy_recemc[num_window_dims]/F");
    m_tree->Branch("num_windows_hcalin", &m_num_windows_hcalin, "num_windows_hcalin[num_window_dims]/I");
    m_tree->Branch("avg_energy_hcalin", &m_avg_energy_hcalin, "avg_energy_hcalin[num_window_dims]/F");
    m_tree->Branch("std_energy_hcalin", &m_std_energy_hcalin, "std_energy_hcalin[num_window_dims]/F");
    m_tree->Branch("num_windows_hcalout", &m_num_windows_hcalout, "num_windows_hcalout[num_window_dims]/I");
    m_tree->Branch("avg_energy_hcalout", &m_avg_energy_hcalout, "avg_energy_hcalout[num_window_dims]/F");
    m_tree->Branch("std_energy_hcalout", &m_std_energy_hcalout, "std_energy_hcalout[num_window_dims]/F");

    int N_window_e_bins = 500;
    float window_e_min =-5.0;
    float window_e_max = 250.0;
    int N_window_e_minus_avg_bins = 500;
    float window_e_minus_avg_min = -100.0;
    float window_e_minus_avg_max = 100.0;


    m_h2_recemc_window_energy_cent.clear();
    m_h2_hcalin_window_energy_cent.clear();
    m_h2_hcalout_window_energy_cent.clear();
    m_h2_full_window_energy_cent.clear();
    m_h2_recemc_window_frac_energy_cent.clear();
    m_h2_hcalin_window_frac_energy_cent.clear();
    m_h2_hcalout_window_frac_energy_cent.clear();
    m_h3_recemc_window_energy_frac_energy_cent.clear();
    m_h3_hcalin_window_energy_frac_energy_cent.clear();
    m_h3_hcalout_window_energy_frac_energy_cent.clear();
    m_h2_recemc_window_energy_minus_avg_energy_cent.clear();
    m_h2_hcalin_window_energy_minus_avg_energy_cent.clear();
    m_h2_hcalout_window_energy_minus_avg_energy_cent.clear();
    m_h2_full_window_energy_minus_avg_energy_cent.clear();
    for ( unsigned int i = 0; i < m_max_window_vector_size; i++) {

     TH2F * h2_cemc = new TH2F(Form("h2_window_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window energy vs centrality CEMC %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          N_window_e_bins, window_e_min, window_e_max, 100, 0, 100);
      h2_cemc->GetXaxis()->SetTitle("Energy (GeV)");
      h2_cemc->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_recemc_window_energy_cent.push_back(h2_cemc);

      TH2F * h2_hcalin = new TH2F(Form("h2_window_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window energy vs centrality HCALIN %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          N_window_e_bins, window_e_min, window_e_max, 100, 0, 100);
      h2_hcalin->GetXaxis()->SetTitle("Energy (GeV)");
      h2_hcalin->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_hcalin_window_energy_cent.push_back(h2_hcalin);

      TH2F * h2_hcalout = new TH2F(Form("h2_window_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window energy vs centrality HCALOUT %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          N_window_e_bins, window_e_min, window_e_max, 100, 0, 100);
      h2_hcalout->GetXaxis()->SetTitle("Energy (GeV)");
      h2_hcalout->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_hcalout_window_energy_cent.push_back(h2_hcalout);

      TH2F * h2_full = new TH2F(Form("h2_window_energy_cent_full_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window energy vs centrality FULL %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          N_window_e_bins, window_e_min, window_e_max, 100, 0, 100);
      h2_full->GetXaxis()->SetTitle("Energy (GeV)");
      h2_full->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_full_window_energy_cent.push_back(h2_full);

      TH2F * h2_frac_energy_recemc = new TH2F(Form("h2_window_frac_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window frac energy vs centrality CEMC %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          100, 0, 1.1, 100, 0, 100);
      h2_frac_energy_recemc->GetXaxis()->SetTitle("Fractional energy");
      h2_frac_energy_recemc->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_recemc_window_frac_energy_cent.push_back(h2_frac_energy_recemc);

      TH2F * h2_frac_energy_hcalin = new TH2F(Form("h2_window_frac_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window frac energy vs centrality HCALIN %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          100, 0, 1.1, 100, 0, 100);
      h2_frac_energy_hcalin->GetXaxis()->SetTitle("Fractional energy");
      h2_frac_energy_hcalin->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_hcalin_window_frac_energy_cent.push_back(h2_frac_energy_hcalin);

      TH2F * h2_frac_energy_hcalout = new TH2F(Form("h2_window_frac_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window frac energy vs centrality HCALOUT %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          100, 0, 1.1, 100, 0, 100);
      h2_frac_energy_hcalout->GetXaxis()->SetTitle("Fractional energy");
      h2_frac_energy_hcalout->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_hcalout_window_frac_energy_cent.push_back(h2_frac_energy_hcalout);

      TH3F * h3_energy_frac_energy_cent_recemc = new TH3F(Form("h3_window_energy_frac_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window energy vs frac energy vs centrality CEMC %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          N_window_e_bins, window_e_min, window_e_max, 100, 0, 1.1, 100, 0, 100);
      h3_energy_frac_energy_cent_recemc->GetXaxis()->SetTitle("Energy (GeV)");
      h3_energy_frac_energy_cent_recemc->GetYaxis()->SetTitle("Fractional energy");
      h3_energy_frac_energy_cent_recemc->GetZaxis()->SetTitle("Centrality (%)");
      m_h3_recemc_window_energy_frac_energy_cent.push_back(h3_energy_frac_energy_cent_recemc);

      TH3F * h3_energy_frac_energy_cent_hcalin = new TH3F(Form("h3_window_energy_frac_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window energy vs frac energy vs centrality HCALIN %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          N_window_e_bins, window_e_min, window_e_max, 100, 0, 1.1, 100, 0, 100);
      h3_energy_frac_energy_cent_hcalin->GetXaxis()->SetTitle("Energy (GeV)");
      h3_energy_frac_energy_cent_hcalin->GetYaxis()->SetTitle("Fractional energy");
      h3_energy_frac_energy_cent_hcalin->GetZaxis()->SetTitle("Centrality (%)");
      m_h3_hcalin_window_energy_frac_energy_cent.push_back(h3_energy_frac_energy_cent_hcalin);

      TH3F * h3_energy_frac_energy_cent_hcalout = new TH3F(Form("h3_window_energy_frac_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window energy vs frac energy vs centrality HCALOUT %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          N_window_e_bins, window_e_min, window_e_max, 100, 0, 1.1, 100, 0, 100);
      h3_energy_frac_energy_cent_hcalout->GetXaxis()->SetTitle("Energy (GeV)");
      h3_energy_frac_energy_cent_hcalout->GetYaxis()->SetTitle("Fractional energy");
      h3_energy_frac_energy_cent_hcalout->GetZaxis()->SetTitle("Centrality (%)");
      m_h3_hcalout_window_energy_frac_energy_cent.push_back(h3_energy_frac_energy_cent_hcalout);

      TH2F * h2_recemc_energy_minus_avg_energy_cent = new TH2F(Form("h2_window_energy_minus_avg_energy_cent_recemc_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window energy minus avg energy vs centrality CEMC %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          N_window_e_minus_avg_bins, window_e_minus_avg_min, window_e_minus_avg_max, 100, 0, 100);
      h2_recemc_energy_minus_avg_energy_cent->GetXaxis()->SetTitle("Energy (GeV)");
      h2_recemc_energy_minus_avg_energy_cent->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_recemc_window_energy_minus_avg_energy_cent.push_back(h2_recemc_energy_minus_avg_energy_cent);

      TH2F * h2_hcalin_energy_minus_avg_energy_cent = new TH2F(Form("h2_window_energy_minus_avg_energy_cent_hcalin_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window energy minus avg energy vs centrality HCALIN %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          N_window_e_minus_avg_bins, window_e_minus_avg_min, window_e_minus_avg_max, 100, 0, 100);
      h2_hcalin_energy_minus_avg_energy_cent->GetXaxis()->SetTitle("Energy (GeV)");
      h2_hcalin_energy_minus_avg_energy_cent->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_hcalin_window_energy_minus_avg_energy_cent.push_back(h2_hcalin_energy_minus_avg_energy_cent);

      TH2F * h2_hcalout_energy_minus_avg_energy_cent = new TH2F(Form("h2_window_energy_minus_avg_energy_cent_hcalout_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window energy minus avg energy vs centrality HCALOUT %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          N_window_e_minus_avg_bins, window_e_minus_avg_min, window_e_minus_avg_max, 100, 0, 100);
      h2_hcalout_energy_minus_avg_energy_cent->GetXaxis()->SetTitle("Energy (GeV)");
      h2_hcalout_energy_minus_avg_energy_cent->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_hcalout_window_energy_minus_avg_energy_cent.push_back(h2_hcalout_energy_minus_avg_energy_cent);

      TH2F * h2_full_energy_minus_avg_energy_cent = new TH2F(Form("h2_window_energy_minus_avg_energy_cent_full_%dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), Form("Window energy minus avg energy vs centrality FULL %dx%d", k_calo_window_dims_hcal_geom[i].first, k_calo_window_dims_hcal_geom[i].second), 
          N_window_e_minus_avg_bins, window_e_minus_avg_min, window_e_minus_avg_max, 100, 0, 100);
      h2_full_energy_minus_avg_energy_cent->GetXaxis()->SetTitle("Energy (GeV)");
      h2_full_energy_minus_avg_energy_cent->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_full_window_energy_minus_avg_energy_cent.push_back(h2_full_energy_minus_avg_energy_cent);

    }

  }

  if ( m_do_calo_cemc_window ) {
    if ( !m_do_calo_window ) {
      m_tree->Branch("num_window_dims", &m_max_window_vector_size, "num_window_dims/i");
    }

    m_tree->Branch("num_windows_cemc", &m_num_windows_cemc, "num_windows_cemc[num_window_dims]/I");
    m_tree->Branch("avg_energy_cemc", &m_avg_energy_cemc, "avg_energy_cemc[num_window_dims]/F");
    m_tree->Branch("std_energy_cemc", &m_std_energy_cemc, "std_energy_cemc[num_window_dims]/F");
    
    int N_window_e_bins = 500;
    float window_e_min = -5.0;
    float window_e_max = 200.0;
    int N_window_e_minus_avg_bins = 500;
    float window_e_minus_avg_min = -100.0;
    float window_e_minus_avg_max = 100.0;
    m_h2_cemc_window_energy_cent.clear();
    m_h2_cemc_window_energy_minus_avg_energy_cent.clear();
    for ( unsigned int i = 0; i < m_max_window_vector_size; i++) {
      TH2F * h2_cemc = new TH2F(Form("h2_window_energy_cent_cemc_%dx%d", k_calo_window_dims_cemc_geom[i].first, k_calo_window_dims_cemc_geom[i].second), Form("Window energy vs centrality CEMC %dx%d", k_calo_window_dims_cemc_geom[i].first, k_calo_window_dims_cemc_geom[i].second), 
          N_window_e_bins, window_e_min, window_e_max, 100, 0, 100);
      h2_cemc->GetXaxis()->SetTitle("Energy (GeV)");
      h2_cemc->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_cemc_window_energy_cent.push_back(h2_cemc);

      TH2F * h2_cemc_energy_minus_avg_energy_cent = new TH2F(Form("h2_window_energy_minus_avg_energy_cent_cemc_%dx%d", k_calo_window_dims_cemc_geom[i].first, k_calo_window_dims_cemc_geom[i].second), Form("Window energy minus avg energy vs centrality CEMC %dx%d", k_calo_window_dims_cemc_geom[i].first, k_calo_window_dims_cemc_geom[i].second), 
          N_window_e_minus_avg_bins, window_e_minus_avg_min, window_e_minus_avg_max, 100, 0, 100);
      h2_cemc_energy_minus_avg_energy_cent->GetXaxis()->SetTitle("Energy (GeV)");
      h2_cemc_energy_minus_avg_energy_cent->GetYaxis()->SetTitle("Centrality (%)");
      m_h2_cemc_window_energy_minus_avg_energy_cent.push_back(h2_cemc_energy_minus_avg_energy_cent);
    }
  }

  if ( Verbosity () > 0 ){
    std::cout << "PPG04AnaWriter::Init - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PPG04AnaWriter::process_event( PHCompositeNode *topNode )
{

  m_event_id++; 

  if(Verbosity() > 1) {
    std::cout << "PPG04AnaWriter::process_event - Process event " << m_event_id << std::endl;
  }

  if ( m_event_id == 0 ) {
    auto rc = recoConsts::instance(); // try to get random seed from flags
    if ( !rc->FlagExist("PPG04RANDOMSEED") ) { // if not set, use time
      std::cout << PHWHERE << "PPG04RANDOMSEED flag not set, not saving random seed." << std::endl;
      m_random_seed = 0;
    }
    else {
      m_random_seed = rc->get_IntFlag("PPG04RANDOMSEED");
    }
  }

  if ( m_do_mbd ) {
    GetMbdInfo(topNode);
  }

  if ( m_do_zvrtx ) {
    if ( GetZvtx(topNode) == Fun4AllReturnCodes::ABORTEVENT ) {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if ( m_do_cent ) {
    if ( GetCentInfo(topNode) == Fun4AllReturnCodes::ABORTEVENT ) {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if ( m_do_tower_background ) {
    GetTowerBackground(topNode);
  }

  if ( m_do_rho ) {
    GetRhoInfo(topNode);
  }

  if ( m_do_calo ) {
    GetCaloInfo(topNode);
  }

  if ( m_do_randomcone ) {
    GetRandomConeInfo(topNode);
  }

  if ( m_do_calo_window ) {
    GetCaloWindowInfo(topNode);
  }

  if ( m_do_calo_cemc_window ) {
    GetCaloCemcWindowInfo(topNode);
  }

  m_tree->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PPG04AnaWriter::End( PHCompositeNode * /*topNode*/ )
{
  
  if(Verbosity() > 0) {
    std::cout << "PPG04AnaWriter::EndRun - End run " << std::endl;
    std::cout << "PPG04AnaWriter::EndRun - Writing to " << m_output_filename << std::endl;
  }

  PHTFileServer::get().cd(m_output_filename);
  
  m_h1_num_events->Fill(0.5, 1.0*m_event_id);
  m_h1_num_events->Write();

  if ( m_do_calo ) {
    for ( unsigned int i = 0; i < m_h3_tower_eta_energy_cent.size(); ++i ) {
      m_h3_tower_eta_energy_cent[i]->Write();
      m_h3_tower_phi_energy_cent[i]->Write();
      m_h2_tower_energy_cent[i]->Write();
      m_h2_tower_dead_phi_eta[i]->Write();
      m_h2_tower_energy_phi_eta[i]->Write();
    }
  }

  if ( m_do_calo_window ) {
    for ( unsigned int i = 0; i < m_h2_recemc_window_energy_cent.size(); ++i ) {
      m_h2_recemc_window_energy_cent[i]->Write();
      m_h2_hcalin_window_energy_cent[i]->Write();
      m_h2_hcalout_window_energy_cent[i]->Write();
      m_h2_full_window_energy_cent[i]->Write();
      m_h2_recemc_window_frac_energy_cent[i]->Write();
      m_h2_hcalin_window_frac_energy_cent[i]->Write();
      m_h2_hcalout_window_frac_energy_cent[i]->Write();
      m_h3_recemc_window_energy_frac_energy_cent[i]->Write();
      m_h3_hcalin_window_energy_frac_energy_cent[i]->Write();
      m_h3_hcalout_window_energy_frac_energy_cent[i]->Write();
      m_h2_recemc_window_energy_minus_avg_energy_cent[i]->Write();
      m_h2_hcalin_window_energy_minus_avg_energy_cent[i]->Write();
      m_h2_hcalout_window_energy_minus_avg_energy_cent[i]->Write();
      m_h2_full_window_energy_minus_avg_energy_cent[i]->Write();
    }
  }

  if ( m_do_calo_cemc_window ) {
    for ( unsigned int i = 0; i < m_h2_cemc_window_energy_cent.size(); ++i ) {
      m_h2_cemc_window_energy_cent[i]->Write();
      m_h2_cemc_window_energy_minus_avg_energy_cent[i]->Write();
    }
  }
  
  m_tree->Write();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PPG04AnaWriter::ResetEvent( PHCompositeNode * /*topNode*/ )
{ 
  // mbd
  m_mbd_q_N = 0;
  m_mbd_q_S = 0;
  m_mbd_time_N = 0;
  m_mbd_time_S = 0;
  m_mbd_npmt_N = 0;
  m_mbd_npmt_S = 0;

  // zvtx
  m_zvtx = 0;

  // centrality
  m_impact_parameter = -1;
  m_centrality = -1;

  // tower background
  m_tower_background_v2 = 0;
  m_tower_background_psi2 = 0;
  m_tower_background_energy_recemc.clear();
  m_tower_background_energy_hcalin.clear();
  m_tower_background_energy_hcalout.clear();

  // rho nodes
  m_rho_val.fill(0);
  m_std_rho_val.fill(0);

  // calo nodes
  m_tower_frac_fired.fill(0);
  m_tower_frac_dead.fill(0);
  m_tower_avg_energy.fill(0);
  m_tower_std_energy.fill(0);
  m_tower_sum_energy.fill(0);

  // random cone
  m_random_cone_R.fill(0);
  m_random_cone_eta.fill(0);
  m_random_cone_phi.fill(0);
  m_random_cone_energy.fill(0);
  m_random_cone_energy_cemc.fill(0);
  m_random_cone_energy_hcalin.fill(0);
  m_random_cone_energy_hcalout.fill(0);
  m_random_cone_num_towers.fill(-1);
  m_random_cone_num_towers_cemc.fill(-1);
  m_random_cone_num_towers_hcalin.fill(-1);
  m_random_cone_num_towers_hcalout.fill(-1);
  m_random_cone_num_masked_towers.fill(-1);
  m_random_cone_num_masked_towers_cemc.fill(-1);
  m_random_cone_num_masked_towers_hcalin.fill(-1);
  m_random_cone_num_masked_towers_hcalout.fill(-1);

  // calo window
  m_num_windows_full.fill(0);
  m_avg_energy_full.fill(0);
  m_std_energy_full.fill(0);
  m_avg_frac_energy_recemc_full.fill(0);
  m_std_frac_energy_recemc_full.fill(0);
  m_avg_frac_energy_hcalin_full.fill(0);
  m_std_frac_energy_hcalin_full.fill(0);
  m_avg_frac_energy_hcalout_full.fill(0);
  m_std_frac_energy_hcalout_full.fill(0);
  m_num_windows_recemc.fill(0);
  m_avg_energy_recemc.fill(0);
  m_std_energy_recemc.fill(0);
  m_num_windows_hcalin.fill(0);
  m_avg_energy_hcalin.fill(0);
  m_std_energy_hcalin.fill(0);
  m_num_windows_hcalout.fill(0);
  m_avg_energy_hcalout.fill(0);
  m_std_energy_hcalout.fill(0);
  m_num_windows_cemc.fill(0);
  m_avg_energy_cemc.fill(0);
  m_std_energy_cemc.fill(0);


  return Fun4AllReturnCodes::EVENT_OK;

}

int PPG04AnaWriter::GetMbdInfo( PHCompositeNode *topNode )
{
  // get MBD info
  auto mbd = findNode::getClass<MbdOutV2>(topNode, m_mbd_node);
  if ( !mbd ) {
    static bool once = true;
    if ( once ) {
      std::cout << PHWHERE << m_mbd_node << " node missing, supressing further warnings." << std::endl;
      once = false;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_mbd_q_N = mbd->get_q(k_mbd_north_code);
  m_mbd_q_S = mbd->get_q(k_mbd_south_code);
  m_mbd_time_N = mbd->get_time(k_mbd_north_code);
  m_mbd_time_S = mbd->get_time(k_mbd_south_code);
  m_mbd_npmt_N = mbd->get_npmt(k_mbd_north_code);
  m_mbd_npmt_S = mbd->get_npmt(k_mbd_south_code);

  if ( Verbosity() > 1 ) {
    std::cout << "PPG04AnaWriter::GetMbdInfo - MBD info N:(q,t,n), S:(q,t,n) = (" << m_mbd_q_N << ", " << m_mbd_time_N << ", " << m_mbd_npmt_N << "), (" << m_mbd_q_S << ", " << m_mbd_time_S << ", " << m_mbd_npmt_S << ")" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PPG04AnaWriter::GetCentInfo( PHCompositeNode *topNode )
{
  // get centrality
  auto cent_node = findNode::getClass< CentralityInfo >(topNode, m_cent_node);
  if ( !cent_node ) {
    std::cout << PHWHERE << m_cent_node << " node missing, doing nothing." << std::endl;
    exit(-1);
  }
  
  if ( cent_node -> has_centrality_bin(CentralityInfo::PROP::mbd_NS) ) {
    m_centrality = (int)(cent_node->get_centrality_bin(CentralityInfo::PROP::mbd_NS));
  } else if ( cent_node -> has_centile(CentralityInfo::PROP::mbd_NS) ){
    m_centrality = (int)(cent_node->get_centile(CentralityInfo::PROP::mbd_NS));
  } 
  if ( cent_node->has_centile(CentralityInfo::PROP::bimp) ){
    m_impact_parameter = (int)(cent_node->get_centile(CentralityInfo::PROP::bimp));
  } 

  if ( m_impact_parameter < 0 && m_centrality < 0 ) {
    std::cerr << PHWHERE << "CentralityInfo Node missing centrality information" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if ( Verbosity() > 1 ) {
    std::cout << "PPG04AnaWriter::GetCentInfo - Centrality = " << m_centrality << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PPG04AnaWriter::GetZvtx( PHCompositeNode *topNode )
{
  // get zvtx
  auto vertexmap = findNode::getClass<GlobalVertexMap>(topNode, m_zvrtx_node);
  if ( !vertexmap  || vertexmap->empty() ) {
    std::cout << PHWHERE << m_zvrtx_node << " node missing, doing nothing." << std::endl;
    exit(-1);
  }

  auto vtx = vertexmap->begin()->second;
  m_zvtx = NAN;
  if ( vtx ) {
    m_zvtx = vtx->get_z();
  } 
  if ( std::isnan(m_zvtx) || m_zvtx > 1e3) {
    static bool once = true;
    if (once) {
      once = false;
      std::cout << PHWHERE << "vertex is " << m_zvtx << ". Drop all tower inputs (further vertex warning will be suppressed)." << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if ( Verbosity() > 1 ) {
    std::cout << "PPG04AnaWriter::GetZvtx - zvtx = " << m_zvtx << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PPG04AnaWriter::GetTowerBackground( PHCompositeNode *topNode )
{
  // get tower background
  auto tower_background = findNode::getClass<TowerBackgroundv1>(topNode, m_tower_background_node);
  if ( !tower_background ) {
    static bool once = true;
    if ( once ) {
      once = false;
      std::cout << PHWHERE << m_tower_background_node << " node missing, doing nothing." << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_tower_background_v2 = tower_background->get_v2();
  m_tower_background_psi2 = tower_background->get_Psi2();
  m_tower_background_energy_recemc = tower_background->get_UE(0);
  m_tower_background_energy_hcalin = tower_background->get_UE(1);
  m_tower_background_energy_hcalout = tower_background->get_UE(2);

  if ( Verbosity() > 1 ) {
    std::cout << "PPG04AnaWriter::GetTowerBackground - v2 = " << m_tower_background_v2 << ", psi2 = " << m_tower_background_psi2 << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

int PPG04AnaWriter::GetRhoInfo( PHCompositeNode *topNode )
{

  // get rho nodes
  for ( unsigned int i = 0; i < m_rho_nodes.size(); ++i ) {

    auto rho = findNode::getClass<TowerRhov1>(topNode, m_rho_nodes[i]);
    if ( !rho ) {
      std::cout << PHWHERE << " Input node " << m_rho_nodes[i] << " Node missing, doing nothing." << std::endl;
      exit(-1); // fatal error
    }

    m_rho_val[i] = rho->get_rho();
    m_std_rho_val[i] = rho->get_sigma();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PPG04AnaWriter::GetCaloInfo( PHCompositeNode *topNode )
{
  for (unsigned int i = 0; i < m_calo_nodes.size(); i++) {
    
    auto towerinfo = findNode::getClass<TowerInfoContainer>( topNode, m_calo_nodes[i] );

    if ( ! towerinfo ) {
      std::cout << PHWHERE << " Input node " << m_calo_nodes[i] << " Node missing, doing nothing." << std::endl;
      exit(-1); // fatal error
    }

    unsigned int ntowers = towerinfo->size();
    int n_fired = 0;
    int n_dead = 0;
    float sum_energy = 0;
    float sum_energy2 = 0;
    for (unsigned int channel = 0; channel < ntowers; channel++) {
      
      auto tower = towerinfo->get_tower_at_channel(channel);
      assert(tower);
      unsigned int key = towerinfo->encode_key(channel);
      float ieta =1.0*towerinfo->getTowerEtaBin(key);
      float iphi = 1.0*towerinfo->getTowerPhiBin(key);
      float energy = tower->get_energy();
      bool is_masked = tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr() || tower->get_isBadChi2() || std::isnan(energy);
      if ( is_masked ) {
        n_dead++;
        m_h2_tower_dead_phi_eta[i]->Fill(iphi, ieta);
      } else if ( energy != 0 ) {        
        n_fired++; 
        sum_energy += energy;
        sum_energy2 += energy*energy;         
        m_h3_tower_eta_energy_cent[i]->Fill(ieta, energy, m_centrality);
        m_h3_tower_phi_energy_cent[i]->Fill(iphi, energy, m_centrality);
        m_h2_tower_energy_phi_eta[i]->Fill(iphi, ieta, energy);
        m_h2_tower_energy_cent[i]->Fill(energy, m_centrality);
      }

    }

    float num_towers = 1.0*ntowers;
    float num_fired = 1.0*n_fired;
    float num_dead = 1.0*n_dead;

    m_tower_frac_fired[i] = num_fired/num_towers;
    m_tower_frac_dead[i] = num_dead/num_towers;
    m_tower_avg_energy[i] = sum_energy/num_fired;
    m_tower_std_energy[i] = sqrt((sum_energy2/num_fired) - (sum_energy/num_fired)*(sum_energy/num_fired));
    m_tower_sum_energy[i] = sum_energy;

    if ( Verbosity() > 1 ) {
      std::cout << "PPG04AnaWriter::GetCaloInfo - Calo " << m_calo_nodes[i] << " fired, dead, avg energy, std energy, sum energy = " << num_fired << ", " << num_dead << ", " << m_tower_avg_energy[i] << ", " << m_tower_std_energy[i] << ", " << m_tower_sum_energy[i] << std::endl;
    }
  
  }  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PPG04AnaWriter::GetRandomConeInfo( PHCompositeNode *topNode )
{
  // get random cone nodes
  for ( unsigned int i = 0; i < m_randomcone_nodes.size(); ++i ) {
    auto rc = findNode::getClass<RandomConev1>(topNode, m_randomcone_nodes[i]);
    if ( !rc ) {
     std::cout << PHWHERE << " Input node " << m_randomcone_nodes[i] << " Node missing, doing nothing." << std::endl;
     exit(-1); // fatal error 
    }

    std::vector<Jet::SRC> srcs = rc->get_src_vec();
    m_random_cone_R[i] = rc->get_R();
    m_random_cone_eta[i] = rc->get_eta();
    m_random_cone_phi[i] = rc->get_phi();
    m_random_cone_energy[i] = rc->get_pt();
    m_random_cone_num_towers[i] = rc->n_clustered();
    m_random_cone_num_masked_towers[i] = rc->n_masked();
    for ( unsigned int j = 0; j < srcs.size(); ++j ) {
      if ( srcs[j] == Jet::SRC::CEMC_TOWERINFO || srcs[j] == Jet::SRC::CEMC_TOWERINFO_RETOWER || srcs[j] == Jet::SRC::CEMC_TOWERINFO_SUB1 || srcs[j] == Jet::SRC::CEMC_TOWERINFO_EMBED || srcs[j] == Jet::SRC::CEMC_TOWERINFO_SIM ) {
        m_random_cone_energy_cemc[i] = rc->get_pt(srcs[j]);
        m_random_cone_num_towers_cemc[i] = rc->n_clustered(srcs[j]);
        m_random_cone_num_masked_towers_cemc[i] = rc->n_masked(srcs[j]);
      } else if ( srcs[j] == Jet::SRC::HCALIN_TOWERINFO || srcs[j] == Jet::SRC::HCALIN_TOWERINFO_SUB1 || srcs[j] == Jet::SRC::HCALIN_TOWERINFO_EMBED || srcs[j] == Jet::SRC::HCALIN_TOWERINFO_SIM ) {
        m_random_cone_energy_hcalin[i] = rc->get_pt(srcs[j]);
        m_random_cone_num_towers_hcalin[i] = rc->n_clustered(srcs[j]);
        m_random_cone_num_masked_towers_hcalin[i] = rc->n_masked(srcs[j]);
      } else if ( srcs[j] == Jet::SRC::HCALOUT_TOWERINFO || srcs[j] == Jet::SRC::HCALOUT_TOWERINFO_SUB1 || srcs[j] == Jet::SRC::HCALOUT_TOWERINFO_EMBED || srcs[j] == Jet::SRC::HCALOUT_TOWERINFO_SIM ) {
        m_random_cone_energy_hcalout[i] = rc->get_pt(srcs[j]);
        m_random_cone_num_towers_hcalout[i] = rc->n_clustered(srcs[j]);
        m_random_cone_num_masked_towers_hcalout[i] = rc->n_masked(srcs[j]);
      }
    }

    if ( Verbosity() > 2 ) {
      std::cout << "PPG04AnaWriter::GetRandomConeInfo - Random cone " << m_randomcone_nodes[i] << " R, eta, phi, energy = " << m_random_cone_R[i] << ", " << m_random_cone_eta[i] << ", " << m_random_cone_phi[i] << ", " << m_random_cone_energy[i] << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PPG04AnaWriter::GetCaloWindowInfo( PHCompositeNode *topNode )
{
  // get calo window nodes
  CaloWindowMapv1 * hcalin_window = findNode::getClass<CaloWindowMapv1>(topNode, m_window_hcalin_node);
  CaloWindowMapv1 * hcalout_window = findNode::getClass<CaloWindowMapv1>(topNode, m_window_hcalout_node);
  CaloWindowMapv1 * recemc_window = findNode::getClass<CaloWindowMapv1>(topNode, m_window_recemc_node);
  if ( !hcalin_window || !hcalout_window || !recemc_window ) {
    std::cout << PHWHERE << " Input node " << m_window_hcalin_node << " or " << m_window_hcalout_node << " or " << m_window_recemc_node << " Node missing, doing nothing." << std::endl;
    exit(-1); // fatal error
  }

  for ( unsigned int i = 0; i < m_max_window_vector_size; ++i ) {

    auto dims = k_calo_window_dims_hcal_geom[i];    
    auto cemc_vec = recemc_window->get_calo_windows(dims.second, dims.first);
    auto hcalin_vec = hcalin_window->get_calo_windows(dims.second, dims.first);
    auto hcalout_vec = hcalout_window->get_calo_windows(dims.second, dims.first);
    if ( dims.second != dims.first ) { // not symmetric
      auto cemc_vec_rot = recemc_window->get_calo_windows(dims.first, dims.second);
      auto hcalin_vec_rot = hcalin_window->get_calo_windows(dims.first, dims.second);
      auto hcalout_vec_rot = hcalout_window->get_calo_windows(dims.first, dims.second);
      cemc_vec.insert(cemc_vec.end(), cemc_vec_rot.begin(), cemc_vec_rot.end());
      hcalin_vec.insert(hcalin_vec.end(), hcalin_vec_rot.begin(), hcalin_vec_rot.end());
      hcalout_vec.insert(hcalout_vec.end(), hcalout_vec_rot.begin(), hcalout_vec_rot.end());
    }

    if ( cemc_vec.size() != hcalin_vec.size() || cemc_vec.size() != hcalout_vec.size() ) {
      std::cout << PHWHERE << " Calo window vectors are not the same size, ignoring" << std::endl;
      continue;
    }
    
    float sum_energy_full2 = 0;
    float sum_energy_recemc2 = 0;
    float sum_energy_hcalout2 = 0;
    float sum_energy_hcalin2 = 0;
      
    float sum_frac_energy_recemc2 = 0;
    float sum_frac_energy_hcalout2 = 0;
    float sum_frac_energy_hcalin2 = 0;

    m_num_windows_full[i] = 0;
    m_avg_energy_full[i] = 0;
    m_std_energy_full[i] = 0;
    m_avg_frac_energy_recemc_full[i] = 0;
    m_std_frac_energy_recemc_full[i] = 0;
    m_avg_frac_energy_hcalin_full[i] = 0;
    m_std_frac_energy_hcalin_full[i] = 0;
    m_avg_frac_energy_hcalout_full[i] = 0;
    m_std_frac_energy_hcalout_full[i] = 0;

    m_num_windows_recemc[i] = 0;
    m_avg_energy_recemc[i] = 0;
    m_std_energy_recemc[i] = 0;

    m_num_windows_hcalin[i] = 0;
    m_avg_energy_hcalin[i] = 0;
    m_std_energy_hcalin[i] = 0;

    m_num_windows_hcalout[i] = 0;
    m_avg_energy_hcalout[i] = 0;
    m_std_energy_hcalout[i] = 0;

    unsigned int n_vec_length = cemc_vec.size();
    for ( unsigned int j = 0; j < n_vec_length; ++j ) {

      int active = 0;
      if ( cemc_vec.at(j) != CaloWindowMap::kMASK_ENERGY ){
        m_num_windows_recemc[i]++;
        m_avg_energy_recemc[i] += cemc_vec.at(j);
        sum_energy_recemc2 += cemc_vec.at(j)*cemc_vec.at(j);
        m_h2_recemc_window_energy_cent[i]->Fill(cemc_vec.at(j), m_centrality);
        active++;
      }
      if ( hcalin_vec.at(j) != CaloWindowMap::kMASK_ENERGY ){
        m_num_windows_hcalin[i]++;
        m_avg_energy_hcalin[i] += hcalin_vec.at(j);
        sum_energy_hcalin2 += hcalin_vec.at(j)*hcalin_vec.at(j);
        m_h2_hcalin_window_energy_cent[i]->Fill(hcalin_vec.at(j), m_centrality);
        active++;
      }
      if ( hcalout_vec.at(j) != CaloWindowMap::kMASK_ENERGY ){
        m_num_windows_hcalout[i]++;
        m_avg_energy_hcalout[i] += hcalout_vec.at(j);
        sum_energy_hcalout2 += hcalout_vec.at(j)*hcalout_vec.at(j);
        m_h2_hcalout_window_energy_cent[i]->Fill(hcalout_vec.at(j), m_centrality);
        active++;
      }
      if ( active == 3 ) {

        m_num_windows_full[i]++;
        float e_full = cemc_vec.at(j) + hcalin_vec.at(j) + hcalout_vec.at(j);
        float cemc_frac = cemc_vec.at(j)/e_full;
        float hcalin_frac = hcalin_vec.at(j)/e_full;
        float hcalout_frac = hcalout_vec.at(j)/e_full;
        m_avg_energy_full[i] += e_full;
        sum_energy_full2 += e_full*e_full;
        m_avg_frac_energy_recemc_full[i] += cemc_frac;
        sum_frac_energy_recemc2 += cemc_frac*cemc_frac;
        m_avg_frac_energy_hcalin_full[i] += hcalin_frac;
        sum_frac_energy_hcalin2 += hcalin_frac*hcalin_frac;
        m_avg_frac_energy_hcalout_full[i] += hcalout_frac;
        sum_frac_energy_hcalout2 += hcalout_frac*hcalout_frac;

        m_h2_full_window_energy_cent[i]->Fill(e_full, m_centrality);
        m_h2_recemc_window_frac_energy_cent[i]->Fill(cemc_frac, m_centrality);
        m_h2_hcalin_window_frac_energy_cent[i]->Fill(hcalin_frac, m_centrality);
        m_h2_hcalout_window_frac_energy_cent[i]->Fill(hcalout_frac, m_centrality);
        m_h3_recemc_window_energy_frac_energy_cent[i]->Fill(e_full, cemc_frac, m_centrality);
        m_h3_hcalin_window_energy_frac_energy_cent[i]->Fill(e_full, hcalin_frac, m_centrality);
        m_h3_hcalout_window_energy_frac_energy_cent[i]->Fill(e_full, hcalout_frac, m_centrality);
      }
      
    }

    m_avg_energy_full[i] /= m_num_windows_full[i];
    m_std_energy_full[i] = sqrt((sum_energy_full2/m_num_windows_full[i]) - (m_avg_energy_full[i]*m_avg_energy_full[i]));
    m_avg_frac_energy_recemc_full[i] /= m_num_windows_full[i];
    m_std_frac_energy_recemc_full[i] = sqrt((sum_frac_energy_recemc2/m_num_windows_full[i]) - (m_avg_frac_energy_recemc_full[i]*m_avg_frac_energy_recemc_full[i]));
    m_avg_frac_energy_hcalin_full[i] /= m_num_windows_full[i];
    m_std_frac_energy_hcalin_full[i] = sqrt((sum_frac_energy_hcalin2/m_num_windows_full[i]) - (m_avg_frac_energy_hcalin_full[i]*m_avg_frac_energy_hcalin_full[i]));
    m_avg_frac_energy_hcalout_full[i] /= m_num_windows_full[i];
    m_std_frac_energy_hcalout_full[i] = sqrt((sum_frac_energy_hcalout2/m_num_windows_full[i]) - (m_avg_frac_energy_hcalout_full[i]*m_avg_frac_energy_hcalout_full[i]));

    m_avg_energy_recemc[i] /= m_num_windows_recemc[i];
    m_std_energy_recemc[i] = sqrt((sum_energy_recemc2/m_num_windows_recemc[i]) - (m_avg_energy_recemc[i]*m_avg_energy_recemc[i]));

    m_avg_energy_hcalin[i] /= m_num_windows_hcalin[i];
    m_std_energy_hcalin[i] = sqrt((sum_energy_hcalin2/m_num_windows_hcalin[i]) - (m_avg_energy_hcalin[i]*m_avg_energy_hcalin[i]));

    m_avg_energy_hcalout[i] /= m_num_windows_hcalout[i];
    m_std_energy_hcalout[i] = sqrt((sum_energy_hcalout2/m_num_windows_hcalout[i]) - (m_avg_energy_hcalout[i]*m_avg_energy_hcalout[i]));

    // reloop to fill minus avg energy
    for ( unsigned int j = 0; j < n_vec_length; ++j ) {
      int active = 0;
      if ( cemc_vec.at(j) != CaloWindowMap::kMASK_ENERGY ){
        m_h2_recemc_window_energy_minus_avg_energy_cent[i]->Fill(cemc_vec.at(j) - m_avg_energy_recemc[i], m_centrality);
        active++;
      }
      if ( hcalin_vec.at(j) != CaloWindowMap::kMASK_ENERGY ){
        m_h2_hcalin_window_energy_minus_avg_energy_cent[i]->Fill(hcalin_vec.at(j) - m_avg_energy_hcalin[i], m_centrality);
        active++;
      }
      if ( hcalout_vec.at(j) != CaloWindowMap::kMASK_ENERGY ){
        m_h2_hcalout_window_energy_minus_avg_energy_cent[i]->Fill(hcalout_vec.at(j) - m_avg_energy_hcalout[i], m_centrality);
        active++;
      }
      if ( active == 3 ) {
        float e_full = cemc_vec.at(j) + hcalin_vec.at(j) + hcalout_vec.at(j);
        m_h2_full_window_energy_minus_avg_energy_cent[i]->Fill(e_full - m_avg_energy_full[i], m_centrality);
      }
    }

  }


  return Fun4AllReturnCodes::EVENT_OK;  
}

int PPG04AnaWriter::GetCaloCemcWindowInfo( PHCompositeNode *topNode )
{
  // get calo window nodes
  CaloWindowMapv1 * cemc_window = findNode::getClass<CaloWindowMapv1>(topNode, m_window_cemc_node);
  if ( !cemc_window ) {
    std::cout << PHWHERE << " Input node " << m_window_cemc_node << " Node missing, doing nothing." << std::endl;
    exit(-1); // fatal error
  }

  for ( unsigned int i = 0; i < m_max_window_vector_size; ++i ) {

    auto dims = k_calo_window_dims_cemc_geom[i];    
    auto cemc_vec = cemc_window->get_calo_windows(dims.second, dims.first);
    if ( dims.second != dims.first ) { // not symmetric
      auto cemc_vec_rot = cemc_window->get_calo_windows(dims.first, dims.second);
      cemc_vec.insert(cemc_vec.end(), cemc_vec_rot.begin(), cemc_vec_rot.end());
    }

    float sum_energy_cemc2 = 0;
    m_num_windows_cemc[i] = 0;
    m_avg_energy_cemc[i] = 0;
    m_std_energy_cemc[i] = 0;
    for ( unsigned int j = 0; j < cemc_vec.size(); ++j ) {
      if ( cemc_vec.at(j) != CaloWindowMap::kMASK_ENERGY ) {
        m_num_windows_cemc[i]++;
        m_avg_energy_cemc[i] += cemc_vec.at(j);
        sum_energy_cemc2 += cemc_vec.at(j)*cemc_vec.at(j);
        m_h2_cemc_window_energy_cent[i]->Fill(cemc_vec.at(j), m_centrality);
      }
    }
    m_avg_energy_cemc[i] /= m_num_windows_cemc[i];
    m_std_energy_cemc[i] = sqrt((sum_energy_cemc2/m_num_windows_cemc[i]) - (m_avg_energy_cemc[i]*m_avg_energy_cemc[i]));
    // reloop to fill minus avg energy
    for ( unsigned int j = 0; j < cemc_vec.size(); ++j ) {
      if ( cemc_vec.at(j) != CaloWindowMap::kMASK_ENERGY ) {
        m_h2_cemc_window_energy_minus_avg_energy_cent[i]->Fill(cemc_vec.at(j) - m_avg_energy_cemc[i], m_centrality);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}