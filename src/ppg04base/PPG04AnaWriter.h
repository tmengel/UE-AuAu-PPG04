#ifndef PPG04BASE_PPG04ANAWRITER_H
#define PPG04BASE_PPG04ANAWRITER_H

#include <fun4all/SubsysReco.h>

#include <centrality/CentralityInfo.h>

#include <string>
#include <utility>  
#include <vector>
#include <array>
#include <cmath>

class PHCompositeNode;
class TTree;
class TH1;
class TH2;
class TH3;

class PPG04AnaWriter : public SubsysReco
{
 public:

  PPG04AnaWriter( const std::string &outputfile = "output.root" );
  ~PPG04AnaWriter() override {};

  int Init(PHCompositeNode * /*topNode*/) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;

  void set_mbd_node(const bool b, const std::string &name = "MbdOut" ) {
    m_do_mbd = b;
    m_mbd_node = name; 
  }

  void set_zvrtx_node(const bool b,  const std::string &name = "GlobalVertexMap" ) {
    m_do_zvrtx = b;
    m_zvrtx_node = name;
  }

  void do_cent_node( const bool b, const std::string &name = "CentralityInfo" ) { 
    m_do_cent = b;
    m_cent_node = name;
  }

  void do_tower_background_node( const bool b, const std::string &name = "TowerInfoBackground_Sub2" ) { 
    m_do_tower_background = b;
    m_tower_background_node = name; 
  }

  void add_calo_node( const std::string &name ) { 
    m_do_calo = true;
    m_calo_nodes.push_back(name); 
  }

  void add_rho_node( const std::string &name ) { 
    m_do_rho = true;
    m_rho_nodes.push_back(name); 
  }

  void add_random_cone_node( const std::string &name ) { 
    m_do_randomcone = true;
    m_randomcone_nodes.push_back(name); 
  }

  void do_calo_window_ana( const std::string &recemc, const std::string &hcalin, const std::string &hcalout ) {
    m_do_calo_window = true;
    m_window_recemc_node = recemc;
    m_window_hcalin_node = hcalin;
    m_window_hcalout_node = hcalout;
  }

  void do_calo_window_cemc_ana( const std::string &name ) { 
    m_do_calo_cemc_window = true;
    m_window_cemc_node = name;
  }


 private:
    
  std::string m_output_filename { "output.root" };

  TTree * m_tree {nullptr};

  // general
  int m_event_id {-1};
  unsigned int m_random_seed { 0 };
  TH1 * m_h1_num_events {nullptr};

  // mbd
  bool m_do_mbd {false};
  std::string m_mbd_node { "MbdOut" };
  const int k_mbd_south_code = 0;
  const int k_mbd_north_code = 1;
  float m_mbd_q_N { 0 };
  float m_mbd_q_S { 0 };
  float m_mbd_time_N { 0 };
  float m_mbd_time_S { 0 };
  float m_mbd_npmt_N { 0 };
  float m_mbd_npmt_S { 0 };

  // zvtx
  bool m_do_zvrtx {true};
  std::string m_zvrtx_node { "GlobalVertexMap" };
  float m_zvtx { 0 };

  // centrality
  bool m_do_cent {true};
  std::string m_cent_node { "CentralityInfo" };
  int m_centrality {-1};
  int m_impact_parameter {-1};

  // towerbackground
  bool m_do_tower_background {false};
  std::string m_tower_background_node { "TowerInfoBackground_Sub2" };
  float m_tower_background_v2 { 0 };
  float m_tower_background_psi2 { 0 };
  std::vector< float > m_tower_background_energy_recemc {};
  std::vector< float > m_tower_background_energy_hcalin {};
  std::vector< float > m_tower_background_energy_hcalout {};
  
  // rho nodes
  bool m_do_rho { false };
  std::vector< std::string > m_rho_nodes {};
  std::array< float, 20 > m_rho_val {};
  std::array< float, 20 > m_std_rho_val {};

  // calo nodes
  bool m_do_calo {false};
  std::vector< std::string > m_calo_nodes {};
  std::vector< TH3 * > m_h3_tower_eta_energy_cent {};
  std::vector< TH3 * > m_h3_tower_phi_energy_cent {};
  std::vector< TH2 * > m_h2_tower_energy_cent {};
  std::vector< TH2 * > m_h2_tower_dead_phi_eta {};
  std::vector< TH2 * > m_h2_tower_energy_phi_eta {};
  std::array < float , 20 > m_tower_frac_fired {};
  std::array < float , 20 > m_tower_frac_dead {};
  std::array < float , 20 > m_tower_avg_energy {};
  std::array < float , 20 > m_tower_std_energy {};
  std::array < float , 20 > m_tower_sum_energy {};
  

  // random cone
  bool m_do_randomcone {false};
  std::vector< std::string > m_randomcone_nodes {};
  std::array< float, 20 > m_random_cone_R {};
  std::array< float, 20 > m_random_cone_eta {};
  std::array< float, 20 > m_random_cone_phi {};
  std::array< float, 20 > m_random_cone_energy {};
  std::array< float, 20 > m_random_cone_energy_cemc {};
  std::array< float, 20 > m_random_cone_energy_hcalin {};
  std::array< float, 20 > m_random_cone_energy_hcalout {};
  std::array< int, 20 >  m_random_cone_num_towers {};
  std::array< int, 20 > m_random_cone_num_towers_cemc {};
  std::array< int, 20 > m_random_cone_num_towers_hcalin {};
  std::array< int, 20 > m_random_cone_num_towers_hcalout {};
  std::array< int, 20 >  m_random_cone_num_masked_towers {};
  std::array< int, 20 > m_random_cone_num_masked_towers_cemc {};
  std::array< int, 20 > m_random_cone_num_masked_towers_hcalin {};
  std::array< int, 20 > m_random_cone_num_masked_towers_hcalout {};

  // calo window
  bool m_do_calo_window { false };
  unsigned int m_max_window_vector_size = 11;  
  const std::vector < std::pair < unsigned int, unsigned int > > k_calo_window_dims_hcal_geom = {
    {1,1}, {2,2}, {3,4}, {5,6}, {7,8}, {9,10}, {11,12}, {13,13}, {15,15}, {17,17}, {20,20}
  }; 
  
  std::string m_window_hcalin_node = "";
  std::string m_window_hcalout_node = "";
  std::string m_window_recemc_node = ""; 

  std::vector < TH2* > m_h2_recemc_window_energy_cent {};
  std::vector < TH2* > m_h2_hcalin_window_energy_cent {};
  std::vector < TH2* > m_h2_hcalout_window_energy_cent {};
  std::vector < TH2* > m_h2_full_window_energy_cent {};
  std::vector < TH2* > m_h2_recemc_window_frac_energy_cent {};
  std::vector < TH2* > m_h2_hcalin_window_frac_energy_cent {};
  std::vector < TH2* > m_h2_hcalout_window_frac_energy_cent {};
  std::vector < TH3* > m_h3_recemc_window_energy_frac_energy_cent {};
  std::vector < TH3* > m_h3_hcalin_window_energy_frac_energy_cent {};
  std::vector < TH3* > m_h3_hcalout_window_energy_frac_energy_cent {};
  std::vector < TH2* > m_h2_recemc_window_energy_minus_avg_energy_cent {};
  std::vector < TH2* > m_h2_hcalin_window_energy_minus_avg_energy_cent {};
  std::vector < TH2* > m_h2_hcalout_window_energy_minus_avg_energy_cent {};
  std::vector < TH2* > m_h2_full_window_energy_minus_avg_energy_cent {};

  std::array < int, 11 > m_num_windows_full {};
  std::array < float, 11 > m_avg_energy_full{};
  std::array < float, 11 > m_std_energy_full{};
  std::array < float, 11 > m_avg_frac_energy_recemc_full{};
  std::array < float, 11 > m_std_frac_energy_recemc_full{};
  std::array < float, 11 > m_avg_frac_energy_hcalin_full{};
  std::array < float, 11 > m_std_frac_energy_hcalin_full{};
  std::array < float, 11 > m_avg_frac_energy_hcalout_full{};
  std::array < float, 11 > m_std_frac_energy_hcalout_full{};

  std::array < int, 11 > m_num_windows_recemc {};
  std::array < float, 11 > m_avg_energy_recemc{};
  std::array < float, 11 > m_std_energy_recemc{};
  std::array < int, 11 > m_num_windows_hcalin {};
  std::array < float, 11 > m_avg_energy_hcalin{};
  std::array < float, 11 > m_std_energy_hcalin{};
  std::array < int, 11 > m_num_windows_hcalout {};
  std::array < float, 11 > m_avg_energy_hcalout{};
  std::array < float, 11 > m_std_energy_hcalout{};

  // calo window
  bool m_do_calo_cemc_window { false };  
  std::string m_window_cemc_node = ""; 
  const std::vector < std::pair < unsigned int, unsigned int > > k_calo_window_dims_cemc_geom = {
    {1,1}, {2,2}, {5,5}, {10,10}, {15,15}, {22,23}, {30,29}, {37,37}, {45,45}, {52,52}, {60,60}
  };
  std::vector < TH2* > m_h2_cemc_window_energy_cent {};
  std::vector < TH2* > m_h2_cemc_window_energy_minus_avg_energy_cent {};

  std::array < int, 11 > m_num_windows_cemc {};
  std::array < float, 11 > m_avg_energy_cemc{};
  std::array < float, 11 > m_std_energy_cemc{};

  int GetMbdInfo( PHCompositeNode *topNode );
  int GetCentInfo( PHCompositeNode *topNode );
  int GetZvtx( PHCompositeNode *topNode );
  int GetTowerBackground( PHCompositeNode *topNode );
  int GetCaloInfo( PHCompositeNode *topNode );
  int GetRhoInfo( PHCompositeNode *topNode );
  int GetRandomConeInfo( PHCompositeNode *topNode );
  int GetCaloWindowInfo( PHCompositeNode *topNode );
  int GetCaloCemcWindowInfo( PHCompositeNode *topNode );


};


#endif
