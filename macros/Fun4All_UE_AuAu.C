// standard includes
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

// macro includes
#include <G4_ActsGeom.C>
#include <G4_Centrality.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <GlobalVariables.C>
#include <Calo_Calib.C> 

// added by me
#include <PPG04.C>

// coresoftware headers

// coresoftware headers
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

// fun4all headerss
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

// centrality headers
#include <centrality/CentralityReco.h>
#include <g4centrality/PHG4CentralityReco.h>

// load libraries
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libcalotrigger.so)


void Fun4All_RandomCones(
  const std::string & dst_calo = "dst_calo_cluster.list",
  const std::string & dst_global = "dst_global.list",
  const std::string & dst_mbd_epd = "dst_mbd_epd.list",
  const std::string & dst_data = "/sphenix/user/tmengel/ppg04/filelists/ana412_2023p015/23745/ppg04/filelists/ana412_2023p015/23745/DST_CALO_RUN1AUAU_run1auau_ana412_2023p015_00023745-0000.list",
  const std::string & outfile = "DST-JET-00042586-0000.root",
  const std::string & mode = "HIJING") 
{
    // set analysis parameters
    Enable::VERBOSITY = 0;
    Enable::DSTOUT = false;
    Enable::IS_DATA = (mode == "DATA");
    Enable::IS_MC = !Enable::IS_DATA;
    Enable::IS_TRUTH_JETS = false;

    PPG04::VERBOSITY = 0;
    PPG04::tower_info_prefix = "TOWERINFO_CALIB";
    PPG04::output_file_name = outfile;
    PPG04::emcal_min_e = 0.05;

    PPG04::do_event_selection = true;
    EventSelection::VERBOSITY = 1;

    EventSelection::do_zvrtx_cut = true;
    EventSelection::SetZertexCutRange(10, -10);
    EventSelection::zvrtx_node = "GlobalVertexMap";

    EventSelection::do_min_bias_cut = true;
    EventSelection::min_bias_node = "MinimumBiasInfo";

    EventSelection::do_truth_jet_cut = false;
    EventSelection::SetPtHardRange(10);
    EventSelection::truth_jet_node = "AntiKt_Truth_r04";

    EventSelection::do_tower_chi2_cut = true;
    EventSelection::tower_chi2_nodes = { "TOWERINFO_CALIB_CEMC", "TOWERINFO_CALIB_HCALIN", "TOWERINFO_CALIB_HCALOUT" };

    PPG04::do_iter_sub = true;
    PPG04::do_area_sub = true;
    PPG04::do_mult_sub = true;
    JetBackground::VERBOSITY = 0;
    JetBackground::did_background = false;
    JetBackground::do_flow = false;
    JetBackground::retower_frac_cut = 0.5;
    JetBackground::seed_jet_d = 3;
    JetBackground::seed_jet_pt = 7;
    JetBackground::area_rho_node = "TowerRho_AREA";
    JetBackground::mult_rho_node = "TowerRho_MULT";
    JetBackground::rho_est_inputs = {
        new TowerJetInput(Jet::CEMC_TOWERINFO, PPG04::tower_info_prefix, PPG04::emcal_min_e),
        new TowerJetInput(Jet::HCALIN_TOWERINFO, PPG04::tower_info_prefix),
        new TowerJetInput(Jet::HCALOUT_TOWERINFO, PPG04::tower_info_prefix)
    };

    PPG04::do_jet_reco = false;
    PPG04::do_truth_jets = false;
    Jets::VERBOSITY = 0;
    Jets::SetJetParams({0.2, 0.3, 0.4, 0.5});

    PPG04::do_random_cone = true;
    RandomCones::VERBOSITY = 1;
    RandomCones::do_basic = true;
    RandomCones::do_rand_eta_phi = true;
    RandomCones::do_avoid_lead_jet = false;
    RandomCones::input_node_type = "CALIB";
    RandomCones::iter_cone_inputs = {
        new TowerJetInput(Jet::CEMC_TOWERINFO_SUB1, PPG04::tower_info_prefix, PPG04::emcal_min_e),
        new TowerJetInput(Jet::HCALIN_TOWERINFO_SUB1, PPG04::tower_info_prefix),
        new TowerJetInput(Jet::HCALOUT_TOWERINFO_SUB1, PPG04::tower_info_prefix)
    };
    RandomCones::rho_cone_inputs = {
        new TowerJetInput(Jet::CEMC_TOWERINFO, PPG04::tower_info_prefix, PPG04::emcal_min_e),
        new TowerJetInput(Jet::HCALIN_TOWERINFO, PPG04::tower_info_prefix),
        new TowerJetInput(Jet::HCALOUT_TOWERINFO, PPG04::tower_info_prefix)
    };
    
  


    // initialize F4A server
    Fun4AllServer * se = Fun4AllServer::instance();
    se -> Verbosity( Enable::VERBOSITY );

    // set up global variables
    recoConsts * rc = recoConsts::instance();
    
    rc -> set_StringFlag( "CDB_GLOBALTAG", "ProdA_2023" );
    rc -> set_uint64Flag( "TIMESTAMP", 23745 );
     CDBInterface::instance()->Verbosity(1);


      SyncReco *sync = new SyncReco();
  se->registerSubsystem(sync);

  HeadReco *head = new HeadReco();
  se->registerSubsystem(head);


//     // connect to conditions database
   
//     // set up flag handler
    FlagHandler* flag = new FlagHandler();
    se -> registerSubsystem(flag);


    // read in filelists
    if ( Enable::IS_DATA ) {
        
        Fun4AllDstInputManager * input = new Fun4AllDstInputManager( "DSTDATA" );
        input -> AddListFile( dst_data );
        se -> registerInputManager( input );
    
    } else {
        
        Fun4AllDstInputManager * input = new Fun4AllDstInputManager( "DSTCALO" );
        input -> AddListFile( dst_calo );
        se -> registerInputManager( input );

        input = new Fun4AllDstInputManager( "DSTGLOBAL" );
        input -> AddListFile( dst_global );
        se -> registerInputManager( input );

        input = new Fun4AllDstInputManager( "DSTMBDEPD" );
        input -> AddListFile( dst_mbd_epd );
        se -> registerInputManager( input );

    }

  
    // do vertex & centrality reconstruction
    Global_Reco();
    // if ( Enable::IS_DATA ) {   Process_Calo_Calib(); }
    // CentralityReco * cent = new CentralityReco();
    // cent -> Verbosity( Enable::VERBOSITY );
    // se -> registerSubsystem( cent );
     PHG4CentralityReco *cent = new PHG4CentralityReco();
    cent->Verbosity(0);
    cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0, string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
    se->registerSubsystem( cent );

    // ppg04
    RunPPG04();

    // run4all
    se -> run( -1 );
    se -> End();
    delete se;

    std::cout << "Done!" << std::endl;
    gSystem -> Exit(0);
  
}
