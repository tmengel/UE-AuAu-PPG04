// standard includes
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <G4_Centrality.C>
#include <G4_Global.C>
#include <GlobalVariables.C>
#include <Calo_Calib.C> 

// coresoftware headers
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

// phool headers
#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>

// fun4all headers
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/InputFileHandler.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

// // event selection headers
#include <eventselection/EventSelector.h>
#include <eventselection/MinBiasCut.h>
#include <eventselection/TowerChi2Cut.h>
#include <eventselection/ZVertexCut.h>
#include <eventselection/LeadTruthJetCut.h>

#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>

#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/DetermineTowerRho.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/TowerRho.h>

#include <ppg04base/CaloTowerManip.h>
#include <ppg04base/CaloSpy.h>
#include <ppg04base/PPG04AnaWriter.h>

#include <underlyingevent/RandomConeTowerReco.h>
#include <underlyingevent/CaloWindowTowerReco.h>

#include <mbd/MbdReco.h>

#include <Sys_Calo.C>



// load libraries
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libppg04base.so)
R__LOAD_LIBRARY(libeventselection.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libunderlyingevent.so)

namespace Enable
{
    bool IS_DATA = false;
    bool IS_MC = false;
    bool IS_TRUTH_JETS = false;
}  // namespace Enable

namespace PPG04 {

    int VERBOSITY = 0;
    bool do_event_selection = false;
    bool do_calo_manip = false;

} // namespace PPG04

namespace EventSelection  
{
    int VERBOSITY = 0;
    bool do_zvrtx_cut = false;
    std::pair<float,float> z_vrtx_cut_range = {20,-20};
    std::string zvrtx_node = "GlobalVertexMap";
    void SetZertexCutRange(float max, float min = -999.0) {
        if (min == -999.0){ min = -max; }
        z_vrtx_cut_range = {max, min};
        return;
    }
    
    bool do_min_bias_cut = false;
    std::string min_bias_node = "MinimumBiasInfo";
    
    bool do_truth_jet_cut = false;
    std::pair<float,float> truth_jet_pT_hat_range =  {10, 30};
    std::string truth_jet_node = "AntiKt_Truth_r04";
    void SetPtHardRange(float pT_hat) {
        float pT_min = 0;
        float pT_max = 0;
        if (pT_hat == 10){
            pT_min = 10;
            pT_max = 30;
        } else if (pT_hat == 30){
            pT_min = 30;
            pT_max = 1000;
        } else {
            std::cout << "Invalid pT range for truth jets" << std::endl;
            gSystem -> Exit(0);
        }
        truth_jet_pT_hat_range = {pT_min, pT_max};
        return;
    }
    
    bool do_tower_chi2_cut = false;
    std::vector<std::string> tower_chi2_nodes = {
        "TOWERINFO_CALIB_CEMC",
        "TOWERINFO_CALIB_HCALIN",
        "TOWERINFO_CALIB_HCALOUT"
    };    

} // namespace EventSelection

namespace CaloManip  {
    
  int VERBOSITY = 0;
  bool do_min_emcal_e_cut = false;
  double emcal_min_e = 0.05; // 50 MeV
  
  bool do_calo_randomizer = false;
}

void GetRunSegment(const std::string & filelist, int & run_number, int & run_segment) 
{
    auto ifm = new InputFileHandler();
    ifm -> AddListFile( filelist );
    auto f = (ifm -> GetFileList()).front();
    std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(f);
    run_number = runseg.first;
    run_segment = runseg.second;
    delete ifm;
    return;
}

std::vector<std::string> GetDSTLists(const std::string & dst_input_lists) 
{
    std::vector<std::string> dst_files {};
    // seperate by space or comma
    std::string delimiter = " ";
    size_t pos = 0;
    std::string token;
    std::string dst_input = dst_input_lists;
    while ((pos = dst_input.find(delimiter)) != std::string::npos) {
        token = dst_input.substr(0, pos);
        dst_files.push_back(token);
        dst_input.erase(0, pos + delimiter.length());
    }
    dst_files.push_back(dst_input);

    for ( auto & f : dst_files ) {
        std::cout << "DST list: " << f << std::endl;
    } 
    return dst_files;
}

void Fun4All_PPG04_Sys( 
    const std::string & mode = "DATA",
    const std::string & prodTag = "ProdA_2023",
    const int timeStamp = 23745,
    const int nEvents = 10,
    const int nsyst = 1,
    const std::string & outdir = ".",
    const std::string & dst_input_list0 = "dst_data.list",
    const std::string & dst_input_list1 = "dst_calo_waveform.list",
    const std::string & dst_input_list2 = "dst_mbd_epd.list"
)
{
    std::vector<std::string> dst_files {};
    if ( mode == "DATA" ) {
        dst_files.push_back(dst_input_list0);
    } else {
        dst_files.push_back(dst_input_list0);
        dst_files.push_back(dst_input_list1);
        dst_files.push_back(dst_input_list2);
    }
    
    // split dst_lists into individual filelists dst1.list, ..., dstN.list
    // std::vector<std::string> dst_files = GetDSTLists(dst_input_lists);
    // if ( dst_files.size() == 0 ) {
    //     std::cout << "No DST files found" << std::endl;
    //     return;
    // }

    int run_number, run_segment;   
    GetRunSegment(dst_files[0], run_number, run_segment);
    CDB::global_tag = prodTag;
    CDB::timestamp = static_cast<uint64_t>(timeStamp);

    std::ostringstream oss;
    oss << outdir << "/" << mode << "-Sys-" << prodTag << "-" << std::setw(6) << std::setfill('0') << timeStamp << "_"
             << std::setw(10) << std::setfill('0') << run_number 
             << "-" << std::setw(6) << std::setfill('0') << run_segment << ".root";
    std::string outfile = oss.str();
    // clear oss
    oss.str("");
    oss << outdir << "/CALO-" << mode << "-Sys-" << prodTag << "-" << std::setw(6) << std::setfill('0') << timeStamp << "_"
             << std::setw(10) << std::setfill('0') << run_number 
             << "-" << std::setw(6) << std::setfill('0') << run_segment << ".root";
    std::string calibfile = oss.str();
    std::cout << "Output file: " << outfile << std::endl;
    
    // --------------------------------------------------------------------------------------------------------
    // set analysis parameters
    // --------------------------------------------------------------------------------------------------------
    Enable::VERBOSITY = 0;
    Enable::CENTRALITY_VERBOSITY = 0;
    Enable::DSTOUT = false;
    PPG04::VERBOSITY = 1;
    Enable::CDB = true;
   
    // set data or mc
    Enable::IS_DATA = (mode == "DATA");
    Enable::IS_MC = !Enable::IS_DATA;
    Enable::IS_TRUTH_JETS = false;

    // event selection
    PPG04::do_event_selection = true;
    EventSelection::VERBOSITY = PPG04::VERBOSITY;
    EventSelection::do_zvrtx_cut = true;
    EventSelection::SetZertexCutRange(20);
    EventSelection::zvrtx_node = "GlobalVertexMap";
    EventSelection::do_min_bias_cut = true;
    EventSelection::min_bias_node = "MinimumBiasInfo";
    EventSelection::do_tower_chi2_cut = true;
    EventSelection::do_truth_jet_cut = false;

    // calo manipulation
    PPG04::do_calo_manip = true;
    CaloManip::VERBOSITY = PPG04::VERBOSITY;
    CaloManip::do_min_emcal_e_cut = true;
    CaloManip::emcal_min_e = 0.05;
    CaloManip::do_calo_randomizer = false;

    // initialize F4A server
    auto se = Fun4AllServer::instance();
    se -> Verbosity( Enable::VERBOSITY );

    // set up global variables
    recoConsts * rc = recoConsts::instance();
    CDBInterface::instance()->Verbosity(1);
    rc -> set_StringFlag( "CDB_GLOBALTAG", CDB::global_tag );
    rc -> set_uint64Flag( "TIMESTAMP", CDB::timestamp );
    rc -> set_IntFlag( "PPG04RANDOMSEED", PHRandomSeed() );

    auto sync = new SyncReco();
    se -> registerSubsystem(sync);
    auto head = new HeadReco();
    se -> registerSubsystem(head);
    auto flag = new FlagHandler();
    se -> registerSubsystem(flag);

    // read in filelists
    int idx = 0;
    for ( auto & dst_file : dst_files ) {
        Fun4AllDstInputManager * input = new Fun4AllDstInputManager( "DSTINPUT_" + std::to_string(idx) );
        input -> AddListFile( dst_file );
        input -> Verbosity( 1 );
        se -> registerInputManager( input );
        idx++;
    }

    // reconstruct
    Global_Reco();
    Process_Calo_Calib();
    Register_Tower_sys();
    if ( !Enable::IS_DATA ) {
        Centrality();
    } 

    if ( PPG04::do_event_selection ) {
      // event selector
      auto es = new EventSelector();
      es -> Verbosity( EventSelection::VERBOSITY );
      if ( Enable::IS_DATA  && EventSelection::do_min_bias_cut ) {
          MinBiasCut * mbc = new MinBiasCut();
          mbc -> SetNodeName( EventSelection::min_bias_node );
          es -> AddCut( mbc );
      }

      if ( EventSelection::do_tower_chi2_cut ) {
        auto tcc = new TowerChi2Cut();
        tcc -> SetNodeNames( {
                        "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "CEMC",
                        "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALIN",
                        "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALOUT"} );
        es -> AddCut( tcc );
      }

      if ( EventSelection::do_zvrtx_cut ) {
        ZVertexCut * zvc =  new ZVertexCut( EventSelection::z_vrtx_cut_range.first, EventSelection::z_vrtx_cut_range.second );
        zvc -> SetNodeName( EventSelection::zvrtx_node );
        es -> AddCut( zvc );
      }
      
      if( Enable::IS_TRUTH_JETS && EventSelection::do_truth_jet_cut ) {
        LeadTruthJetCut * ltjc = new LeadTruthJetCut( EventSelection::truth_jet_pT_hat_range.first, 
                                                      EventSelection::truth_jet_pT_hat_range.second );
        ltjc -> SetNodeName( EventSelection::truth_jet_node );
        es -> AddCut( ltjc );
      }
      
      
      se->registerSubsystem(es);
      es -> PrintCuts();
    }
     
    if ( PPG04::do_calo_manip ) {

        auto ctm = new CaloTowerManip();
        ctm -> Verbosity( CaloManip::VERBOSITY );
        ctm -> SetInputNode( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "CEMC" );
        ctm -> SaveCopyOutputNode( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_CEMC" );
        se -> registerSubsystem( ctm );


        ctm = new CaloTowerManip();
        ctm -> Verbosity( CaloManip::VERBOSITY );
        ctm -> SetInputNode( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_CEMC" );
        ctm -> SaveCopyOutputNode( "TOWERINFO_CALIB_CEMC_ORIGINAL" );
        if ( CaloManip::do_min_emcal_e_cut ) {
            ctm -> SetMinEnergy( CaloManip::emcal_min_e );
        }
        if ( CaloManip::do_calo_randomizer ) {
            ctm -> RandomizeTowers( true );
        }
        se -> registerSubsystem( ctm );

        auto rcemc = new RetowerCEMC(); 
        rcemc -> Verbosity( 0 );
        rcemc -> set_towerinfo( true );
        rcemc -> set_frac_cut( 0.5 );
        rcemc -> set_towerNodePrefix( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) );
        se -> registerSubsystem( rcemc );

        if ( CaloManip::do_calo_randomizer ) {

            ctm = new CaloTowerManip();
            ctm -> Verbosity( CaloManip::VERBOSITY );
            ctm -> SetInputNode( "TOWERINFO_CALIB_CEMC_RETOWER" );
            ctm -> SaveCopyOutputNode( "TOWERINFO_CALIB_CEMC_RETOWER_ORIGINAL" );
            // if ( CaloManip::do_min_emcal_e_cut ) {
            //     ctm -> SetMinEnergy( CaloManip::emcal_min_e );
            // }
            ctm -> RandomizeTowers( true );
            se -> registerSubsystem( ctm );

            ctm = new CaloTowerManip();
            ctm -> Verbosity( CaloManip::VERBOSITY );
            ctm -> SetInputNode( "TOWERINFO_CALIB_HCALIN" );
            ctm -> SaveCopyOutputNode( "TOWERINFO_CALIB_HCALIN_ORIGINAL" );
            ctm -> RandomizeTowers( true );
            se -> registerSubsystem( ctm );

            ctm = new CaloTowerManip();
            ctm -> Verbosity( CaloManip::VERBOSITY );
            ctm -> SetInputNode( "TOWERINFO_CALIB_HCALOUT" );
            ctm -> SaveCopyOutputNode( "TOWERINFO_CALIB_HCALOUT_ORIGINAL" );
            ctm -> RandomizeTowers( true );
            se -> registerSubsystem( ctm );
        }

    } else {

        auto rcemc = new RetowerCEMC(); 
        rcemc -> Verbosity( 0 );
        rcemc -> set_towerinfo( true );
        rcemc -> set_frac_cut( 0.5 );
        rcemc -> set_towerNodePrefix( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) );
        se -> registerSubsystem( rcemc );
    }

    JetReco * ijr = new JetReco();
    ijr -> add_input( new TowerJetInput( Jet::CEMC_TOWERINFO_RETOWER, "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) ) );
    ijr -> add_input( new TowerJetInput( Jet::HCALIN_TOWERINFO, "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) ) );
    ijr -> add_input( new TowerJetInput( Jet::HCALOUT_TOWERINFO, "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) ) );
    ijr -> add_algo( new FastJetAlgoSub( Jet::ANTIKT, 0.2 ), "AntiKt_TowerInfo_HIRecoSeedsRaw_r02" );
    ijr -> set_algo_node("ANTIKT");
    ijr -> set_input_node("TOWER");
    ijr -> Verbosity( Enable::VERBOSITY );
    se -> registerSubsystem(ijr);

    DetermineTowerBackground *dtb = new DetermineTowerBackground();
    dtb->SetBackgroundOutputName("TowerInfoBackground_Sub1");
    dtb->SetFlow(false);
    dtb->SetSeedType(0);
    dtb->SetSeedJetD(3);
    dtb->set_towerinfo(true);
    dtb->Verbosity( Enable::VERBOSITY ); 
    dtb->set_towerNodePrefix("TOWERINFO_CALIB_SYST" + std::to_string(nsyst));
    se->registerSubsystem(dtb);

    CopyAndSubtractJets *casj = new CopyAndSubtractJets();
    casj->SetFlowModulation( false );
    casj->Verbosity(Enable::VERBOSITY); 
    casj->set_towerinfo(true);
    casj->set_towerNodePrefix("TOWERINFO_CALIB_SYST" + std::to_string(nsyst));
    se->registerSubsystem(casj);
  
    DetermineTowerBackground *dtb2 = new DetermineTowerBackground();
    dtb2->SetBackgroundOutputName("TowerInfoBackground_Sub2");
    dtb2->SetFlow( false );
    dtb2->SetSeedType(1);
    dtb2->SetSeedJetPt(7);
    dtb2->Verbosity(Enable::VERBOSITY); 
    dtb2->set_towerinfo(true);
    dtb2->set_towerNodePrefix("TOWERINFO_CALIB_SYST" + std::to_string(nsyst));
    se->registerSubsystem(dtb2);

    SubtractTowers *st = new SubtractTowers();
    st->SetFlowModulation( false );
    st->Verbosity( Enable::VERBOSITY );
    st->set_towerinfo(true);
    st->set_towerNodePrefix("TOWERINFO_CALIB_SYST" + std::to_string(nsyst));
    se->registerSubsystem(st);

    DetermineTowerRho * dtr = new DetermineTowerRho();
    dtr -> add_method( TowerRho::Method::AREA, "TowerRho_AREA" );
    dtr -> add_method( TowerRho::Method::MULT, "TowerRho_MULT" );
    dtr -> Verbosity( Enable::VERBOSITY  );
    dtr -> add_tower_input( new TowerJetInput( Jet::CEMC_TOWERINFO, "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) ) );
    dtr -> add_tower_input( new TowerJetInput( Jet::HCALIN_TOWERINFO, "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) ) );
    dtr -> add_tower_input( new TowerJetInput( Jet::HCALOUT_TOWERINFO, "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) ) );
    se -> registerSubsystem( dtr );

    dtr = new DetermineTowerRho();
    dtr -> add_method( TowerRho::Method::AREA, "TowerRho_AREA_CEMC" );
    dtr -> add_method( TowerRho::Method::MULT, "TowerRho_MULT_CEMC" );
    dtr -> Verbosity( Enable::VERBOSITY );
    dtr -> add_tower_input( new TowerJetInput( Jet::CEMC_TOWERINFO, "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) ) );
    se -> registerSubsystem( dtr );

    dtr = new DetermineTowerRho();
    dtr -> add_method( TowerRho::Method::AREA, "TowerRho_AREA_HCALIN" );
    dtr -> add_method( TowerRho::Method::MULT, "TowerRho_MULT_HCALIN" );
    dtr -> Verbosity( Enable::VERBOSITY );
    dtr -> add_tower_input( new TowerJetInput( Jet::HCALIN_TOWERINFO, "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) ) );
    se -> registerSubsystem( dtr );

    dtr = new DetermineTowerRho();
    dtr -> add_method( TowerRho::Method::AREA, "TowerRho_AREA_HCALOUT" );
    dtr -> add_method( TowerRho::Method::MULT, "TowerRho_MULT_HCALOUT" );
    dtr -> Verbosity( Enable::VERBOSITY  );
    dtr -> add_tower_input( new TowerJetInput( Jet::HCALOUT_TOWERINFO, "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) ) );
    se -> registerSubsystem( dtr );


    CaloSpy * ch = new CaloSpy( calibfile );
    ch -> AddCaloNode( "TOWERINFO_CALIB_SYST" +std::to_string(nsyst)+ "_CEMC", 256, 96 );
    ch -> AddCaloNode( "TOWERINFO_CALIB_SYST" +std::to_string(nsyst)+ "_HCALIN", 64, 24 );
    ch -> AddCaloNode( "TOWERINFO_CALIB_SYST" +std::to_string(nsyst)+ "_HCALOUT", 64, 24 );
    ch -> AddCaloNode( "TOWERINFO_CALIB_SYST" +std::to_string(nsyst)+ "_CEMC_RETOWER", 64, 24 );
    ch -> AddCaloNode( "TOWERINFO_CALIB_SYST" +std::to_string(nsyst)+ "_CEMC_RETOWER_SUB1", 64, 24 );
    ch -> AddCaloNode( "TOWERINFO_CALIB_SYST" +std::to_string(nsyst)+ "_HCALIN_SUB1", 64, 24 );
    ch -> AddCaloNode( "TOWERINFO_CALIB_SYST" +std::to_string(nsyst)+ "_HCALOUT_SUB1", 64, 24 );
    if( PPG04::do_calo_manip ) {
        ch -> AddCaloNode( "TOWERINFO_CALIB_SYST" +std::to_string(nsyst)+ "_CEMC_ORIGINAL", 256, 96 );
        if ( CaloManip::do_calo_randomizer ) {
            ch -> AddCaloNode( "TOWERINFO_CALIB_CEMC_RETOWER_ORIGINAL", 64, 24 );
            ch -> AddCaloNode( "TOWERINFO_CALIB_HCALIN_ORIGINAL", 64, 24 );
            ch -> AddCaloNode( "TOWERINFO_CALIB_HCALOUT_ORIGINAL", 64, 24 );
        }
    }
    ch -> Normalize( false );
    se -> registerSubsystem( ch );

    RandomConeTowerReco * rct = new RandomConeTowerReco();
    rct -> Verbosity( PPG04::VERBOSITY );
    rct -> add_input_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_CEMC" , "TOWERGEOM_CEMC" , Jet::CEMC_TOWERINFO );
    rct -> add_input_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALIN" , "TOWERGEOM_HCALIN" , Jet::HCALIN_TOWERINFO );
    rct -> add_input_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALOUT" , "TOWERGEOM_HCALOUT" , Jet::HCALOUT_TOWERINFO );
    rct -> set_R( 0.4 );
    rct -> set_abs_eta( 1.1 );
    rct -> set_output_node("RandomCones_r04");
    rct -> set_masked_threshold(0.05);
    se -> registerSubsystem( rct );

    rct = new RandomConeTowerReco();
    rct -> Verbosity( PPG04::VERBOSITY );
    rct -> add_input_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_CEMC_RETOWER" , "TOWERGEOM_HCALIN" , Jet::CEMC_TOWERINFO_SUB1 );
    rct -> add_input_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALIN" , "TOWERGEOM_HCALIN" , Jet::HCALIN_TOWERINFO_SUB1 );
    rct -> add_input_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALOUT" , "TOWERGEOM_HCALOUT" , Jet::HCALOUT_TOWERINFO_SUB1 );   
    rct -> set_R( 0.4 );
    rct -> set_abs_eta( 1.1 );
    rct -> set_output_node("RandomCones_r04_Sub1");
    rct -> set_masked_threshold(0.05);
    se -> registerSubsystem( rct );


    CaloWindowTowerReco * cwt = new CaloWindowTowerReco();
    cwt -> Verbosity( PPG04::VERBOSITY );
    cwt -> add_input_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_CEMC" , "TOWERGEOM_CEMC" , Jet::CEMC_TOWERINFO );
    cwt -> add_input_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALIN" , "TOWERGEOM_HCALIN" , Jet::HCALIN_TOWERINFO );
    cwt -> add_input_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALOUT" , "TOWERGEOM_HCALOUT" , Jet::HCALOUT_TOWERINFO );
    cwt -> add_input_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_CEMC_RETOWER" , "TOWERGEOM_HCALIN" , Jet::CEMC_TOWERINFO );
    cwt -> add_window_node( "CaloWindowMap_CEMC" );
    cwt -> add_window_node( "CaloWindowMap_HCALIN" );
    cwt -> add_window_node( "CaloWindowMap_HCALOUT" );
    cwt -> add_window_node( "CaloWindowMap_CEMC_RETOWER" );
    se -> registerSubsystem( cwt );

    PPG04AnaWriter * aw = new PPG04AnaWriter( outfile );
    aw -> Verbosity( PPG04::VERBOSITY );
    aw -> set_mbd_node( true, "MbdOut" );
    aw -> set_zvrtx_node( true, "GlobalVertexMap");
    aw -> do_cent_node( true );
    aw -> do_tower_background_node(true, "TowerInfoBackground_Sub2");
    aw -> add_calo_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_CEMC" );
    aw -> add_calo_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALIN" );
    aw -> add_calo_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALOUT" );
    aw -> add_calo_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_CEMC_RETOWER" );
    aw -> add_calo_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_CEMC_RETOWER_SUB1" );
    aw -> add_calo_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALIN_SUB1" );
    aw -> add_calo_node( "TOWERINFO_CALIB_SYST" + std::to_string(nsyst) + "_HCALOUT_SUB1" );
    // aw -> add_calo_node( "TOWERINFO_CALIB_HCALIN" );
    // aw -> add_calo_node( "TOWERINFO_CALIB_HCALOUT" );
    // aw -> add_calo_node( "TOWERINFO_CALIB_CEMC_RETOWER" );
    // aw -> add_calo_node( "TOWERINFO_CALIB_CEMC_RETOWER_SUB1" );
    // aw -> add_calo_node( "TOWERINFO_CALIB_HCALIN_SUB1" );
    // aw -> add_calo_node( "TOWERINFO_CALIB_HCALOUT_SUB1" );
    aw -> add_rho_node( "TowerRho_AREA" );
    aw -> add_rho_node( "TowerRho_MULT" );
    aw -> add_rho_node( "TowerRho_AREA_CEMC" );
    aw -> add_rho_node( "TowerRho_MULT_CEMC" );
    aw -> add_rho_node( "TowerRho_AREA_HCALIN" );
    aw -> add_rho_node( "TowerRho_MULT_HCALIN" );
    aw -> add_rho_node( "TowerRho_AREA_HCALOUT" );
    aw -> add_rho_node( "TowerRho_MULT_HCALOUT" );
    aw -> add_random_cone_node( "RandomCones_r04" );
    aw -> add_random_cone_node( "RandomCones_r04_Sub1" );
    aw -> do_calo_window_ana( "CaloWindowMap_CEMC_RETOWER", "CaloWindowMap_HCALIN", "CaloWindowMap_HCALOUT" );
    aw -> do_calo_window_cemc_ana( "CaloWindowMap_CEMC" );
    se -> registerSubsystem( aw );




    // run4all
    se -> run( nEvents );
    se -> End();

    std::cout << "Done!" << std::endl;
    gSystem -> Exit(0);
  
}
