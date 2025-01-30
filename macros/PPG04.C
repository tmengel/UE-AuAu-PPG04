#ifndef MACRO_PPG04_C
#define MACRO_PPG04_C

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

// // event selection headers
#include <eventselection/EventSelector.h>
#include <eventselection/MinBiasCut.h>
#include <eventselection/TowerChi2Cut.h>
#include <eventselection/ZVertexCut.h>
#include <eventselection/LeadJetCut.h>

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

R__LOAD_LIBRARY( libfun4all.so )
R__LOAD_LIBRARY( libppg04base.so )
R__LOAD_LIBRARY( libeventselection.so )
R__LOAD_LIBRARY( libjetbackground.so )
R__LOAD_LIBRARY( libjetbase.so )
R__LOAD_LIBRARY( libunderlyingevent.so )

namespace PPG04 {

    int VERBOSITY = 0;

    bool isDATA = false;
    bool isMC = false;
    bool isTRUTHJETS = false;

    int PPG04RandomSeed = 42;

    std::string MbdNode = "MbdOut";
    std::string MinBiasNode = "MinimumBiasInfo";
    std::string ZVrtxNode = "GlobalVertexMap";
    std::string CentNode = "CentralityInfo";
    std::string TowerBackgroundNode = "TowerInfoBackground_Sub2";

    std::string TowerPrefix = "TOWERINFO_CALIB";

    bool doEventSelect = false;
    bool doCaloManip = false;

    bool doIterBackground = false;
    bool doAreaRho = false;
    bool doMultRho = false;
    
    bool doRandomCones = false;
    bool doCaloWindows = false;

    bool doCaloSpy = false;
    bool doAnaWriter = false;

    EventSelector * EventSelectorHandler {nullptr};
    
    std::vector< DetermineTowerRho* > RhoRecos {};
    std::vector< RandomConeTowerReco* > RandomConesRecos {};
    std::vector< CaloWindowTowerReco* > CaloWindowRecos {};

    CaloSpy * CaloSpyHandler {nullptr};
    PPG04AnaWriter * AnaWriterHandler {nullptr};

    
} // namespace PPG04

namespace EventSelect {
    int VERBOSITY = 0;
    bool doZVrtxCut = false;
    std::pair<float,float> ZVrtxCutRange = {20,-20};
    bool doMinBiasCut = false;
    bool doTowerChi2Cut = false;
    std::vector<std::string> TowerChi2Nodes = {
        "TOWERINFO_CALIB_CEMC",
        "TOWERINFO_CALIB_HCALIN",
        "TOWERINFO_CALIB_HCALOUT"
    };    
} // namespace EventSelect

namespace CaloManip {  
    int VERBOSITY = 0;
    bool doMinEMCalEnergy = false;
    double MinEMCalEnergy = 0.05; // 50 MeV
    bool doTowerRandomizer = false;
} // namespace CaloManip

namespace RandomCones {
    int VERBOSITY = 0;
    float ConeRadius = 0.4;
    float ConeAbsEta = 1.1;
    float ConeMaskedThreshold = 0.05;
} // namespace RandomCones

namespace CaloWindows {
    int VERBOSITY = 0;
    std::string WindowPrefix = "CaloWindowMap";    
} // namespace CaloWindows

namespace PPG04CaloSpy {
    bool Normalize = false;
    std::string outfile = "calo-spy-output.root";
    std::vector<std::string> CaloSpyNodes {};   
} // namespace PPG04CaloSpy

namespace PPG04Output {
    std::string outfile = "ppg04-ana-output.root";
    bool writeMBD = false;
    bool writeZVtx = false;
    bool writeCent = false;
    bool writeIterBackground = false;
    std::vector<std::string> CaloNodes {};
    std::vector<std::string> RhoNodes {};
    std::vector<std::string> RandomConesNodes {};
    bool doFullWindow = false;
    std::string FullWindowRetowerNode;
    std::string FullWindowHCALINNode;
    std::string FullWindowHCALOUTNode;
    bool doCemcOnlyWindow = false;
    std::string CemcOnlyNode;
} // namespace PPG04Output

void AddRhoReco( const std::vector< Jet::SRC > srcs, const std::string & output_suffix = "", const bool doArea = PPG04::doAreaRho, const bool doMult = PPG04::doMultRho )
{
    unsigned int idx = PPG04::RhoRecos.size();
    auto dtr = new DetermineTowerRho( "DetermineTowerRho_" + std::to_string( idx ) );
    dtr -> Verbosity( Enable::VERBOSITY );
    if ( doArea ) {
        std::string outnode = "TowerRho_AREA";
        if ( output_suffix != "" ){ outnode += "_" + output_suffix; }
        dtr -> add_method( TowerRho::Method::AREA, outnode );
        PPG04Output::RhoNodes.push_back( outnode );
    }
    if ( doMult ) {
        std::string outnode = "TowerRho_MULT";
        if ( output_suffix != "" ){ outnode += "_" + output_suffix; }
        dtr -> add_method( TowerRho::Method::MULT, outnode );
        PPG04Output::RhoNodes.push_back( outnode );
    }
    for ( auto src : srcs ) {
        dtr -> add_tower_input( new TowerJetInput( src, PPG04::TowerPrefix ) );
    }
    PPG04::RhoRecos.push_back( dtr );
    return;
}

void AddRandomConeReco( const std::vector< Jet::SRC > srcs, const std::string & output_node, float R = RandomCones::ConeRadius, float abs_eta = RandomCones::ConeAbsEta, float masked_threshold = RandomCones::ConeMaskedThreshold )
{
    PPG04::doRandomCones = true;
    unsigned int idx = PPG04::RandomConesRecos.size();
    auto rct = new RandomConeTowerReco("RandomConeTowerReco_" + std::to_string( idx ) );
    rct -> Verbosity( PPG04::VERBOSITY );
    for ( auto src : srcs ) {
        rct -> add_input( src, PPG04::TowerPrefix );
    }
    rct -> set_R( R );
    rct -> set_abs_eta( abs_eta );
    rct -> set_masked_threshold( masked_threshold );   
    PPG04Output::RandomConesNodes.push_back( output_node );
    rct -> set_output_node( output_node );
    PPG04::RandomConesRecos.push_back( rct );
    return;
}

void AddFullCaloWindowReco()
{
    std::vector< Jet::SRC > srcs = { Jet::SRC::CEMC_TOWERINFO_RETOWER, Jet::SRC::HCALIN_TOWERINFO, Jet::SRC::HCALOUT_TOWERINFO };
    unsigned int idx = PPG04::CaloWindowRecos.size();
    auto cwr = new CaloWindowTowerReco( "CaloWindowTowerReco_" + std::to_string( idx ) );
    cwr -> Verbosity( PPG04::VERBOSITY );
    for ( auto src : srcs ) {
        cwr -> add_input( src, PPG04::TowerPrefix );
    }
    // cwr -> set_window_prefix( CaloWindows::WindowPrefix );
    PPG04Output::FullWindowRetowerNode = CaloWindows::WindowPrefix + "_CEMC_RETOWER";
    PPG04Output::FullWindowHCALINNode = CaloWindows::WindowPrefix + "_HCALIN";
    PPG04Output::FullWindowHCALOUTNode = CaloWindows::WindowPrefix + "_HCALOUT";
    PPG04::CaloWindowRecos.push_back( cwr );

    return;
}

void AddCemcCaloWindowReco()
{
    std::vector< Jet::SRC > srcs = { Jet::SRC::CEMC_TOWERINFO };
    unsigned int idx = PPG04::CaloWindowRecos.size();
    auto cwr = new CaloWindowTowerReco( "CaloWindowTowerReco_" + std::to_string( idx ) );
    cwr -> Verbosity( PPG04::VERBOSITY );
    for ( auto src : srcs ) {
        cwr -> add_input( src, PPG04::TowerPrefix );
    }
    cwr -> set_window_prefix( CaloWindows::WindowPrefix );
    PPG04Output::CemcOnlyNode = CaloWindows::WindowPrefix + "_CEMC";
    PPG04::CaloWindowRecos.push_back( cwr );

    return;
}

void InitPPG04()
{
    
    if ( PPG04::doAreaRho || PPG04::doMultRho ) {

        AddRhoReco( { Jet::CEMC_TOWERINFO, Jet::HCALIN_TOWERINFO, Jet::HCALOUT_TOWERINFO });
        AddRhoReco( { Jet::CEMC_TOWERINFO }, "CEMC" );
        AddRhoReco( { Jet::HCALIN_TOWERINFO }, "HCALIN" );
        AddRhoReco( { Jet::HCALOUT_TOWERINFO }, "HCALOUT" );
    
    }

    if ( PPG04::doRandomCones ) {
    
        AddRandomConeReco( { Jet::CEMC_TOWERINFO, Jet::HCALIN_TOWERINFO, Jet::HCALOUT_TOWERINFO }, "RandomCones_r04" );
        
        if ( PPG04::doIterBackground ) {
            AddRandomConeReco( { Jet::CEMC_TOWERINFO_SUB1, Jet::HCALIN_TOWERINFO_SUB1, Jet::HCALOUT_TOWERINFO_SUB1 }, "RandomCones_r04_Sub1" );
        }

    }

    if ( PPG04::doCaloWindows ) {
    
        if ( PPG04Output::doFullWindow ) {
            AddFullCaloWindowReco();
        }

        if ( PPG04Output::doCemcOnlyWindow ) {
            AddCemcCaloWindowReco();
        }

    }
}

void RunPPG04()
{
    if ( PPG04::VERBOSITY > 0 ) { std::cout << "InitPPG04 - Initializing PPG04 Analysis" << std::endl; }
    if ( PPG04::isDATA ) {
        if ( PPG04::VERBOSITY > 0 ) { std::cout << "InitPPG04 - Running on DATA" << std::endl; }
    }
    if ( PPG04::isMC ) {
        if ( PPG04::VERBOSITY > 0 ) { std::cout << "InitPPG04 - Running on MC" << std::endl; }
    }
    if ( PPG04::isTRUTHJETS ) {
        if ( PPG04::VERBOSITY > 0 ) { std::cout << "InitPPG04 - Running on TRUTH JETS" << std::endl; }
    }

    auto se = Fun4AllServer::instance();
    
    if ( PPG04::doEventSelect ) { 
        PPG04::EventSelectorHandler = new EventSelector();
        PPG04::EventSelectorHandler -> Verbosity( EventSelect::VERBOSITY );
        if ( PPG04::isDATA  && EventSelect::doMinBiasCut ) {
            auto mbc = new MinBiasCut( );
            mbc -> SetNodeName( PPG04::MinBiasNode );
            PPG04::EventSelectorHandler -> AddCut( mbc );
        }
        if ( EventSelect::doTowerChi2Cut ) {
            auto tcc = new TowerChi2Cut( );
            tcc -> SetNodeNames( EventSelect::TowerChi2Nodes );
            PPG04::EventSelectorHandler -> AddCut( tcc );
        }
        if ( EventSelect::doZVrtxCut ) {
            auto zvc =  new ZVertexCut( EventSelect::ZVrtxCutRange.first, EventSelect::ZVrtxCutRange.second );
            zvc -> SetNodeName( PPG04::ZVrtxNode );
            PPG04::EventSelectorHandler -> AddCut( zvc );
        }
        
        PPG04::EventSelectorHandler -> PrintCuts();
        se -> registerSubsystem( PPG04::EventSelectorHandler );
    }

    if ( PPG04::doCaloManip ) { 
        int verbosity = std::max( PPG04::VERBOSITY, CaloManip::VERBOSITY );
        auto ctm = new CaloTowerManip( );
        ctm -> Verbosity( verbosity );
        ctm -> SetInputNode( "TOWERINFO_CALIB_CEMC" );
        ctm -> SaveCopyOutputNode( "TOWERINFO_CALIB_CEMC_ORIGINAL" );
        if ( CaloManip::doMinEMCalEnergy ) {
            ctm -> SetMinEnergy( CaloManip::MinEMCalEnergy );
        }
        if ( CaloManip::doTowerRandomizer ) {
            ctm -> RandomizeTowers( true );
        }
        se -> registerSubsystem( ctm );
        PPG04CaloSpy::CaloSpyNodes.push_back( "TOWERINFO_CALIB_CEMC" );
        PPG04CaloSpy::CaloSpyNodes.push_back( "TOWERINFO_CALIB_CEMC_ORIGINAL" );
    }

    auto rcemc = new RetowerCEMC();
    rcemc -> Verbosity( Enable::VERBOSITY );
    rcemc -> set_towerinfo( true );
    rcemc -> set_frac_cut( 0.5 );
    rcemc -> set_towerNodePrefix( PPG04::TowerPrefix );
    se -> registerSubsystem( rcemc );

    if ( PPG04::doCaloManip && CaloManip::doTowerRandomizer ) { 
        
        int verbosity = std::max( PPG04::VERBOSITY, CaloManip::VERBOSITY );

        auto ctm = new CaloTowerManip( );
        ctm -> Verbosity( verbosity );
        ctm -> SetInputNode( "TOWERINFO_CALIB_CEMC_RETOWER" );
        ctm -> SaveCopyOutputNode( "TOWERINFO_CALIB_CEMC_RETOWER_ORIGINAL" );
        ctm -> RandomizeTowers( true );
        se -> registerSubsystem( ctm );
        PPG04CaloSpy::CaloSpyNodes.push_back( "TOWERINFO_CALIB_CEMC_RETOWER" );
        PPG04CaloSpy::CaloSpyNodes.push_back( "TOWERINFO_CALIB_CEMC_RETOWER_ORIGINAL" );

        ctm = new CaloTowerManip( );
        ctm -> Verbosity( verbosity );
        ctm -> SetInputNode( "TOWERINFO_CALIB_HCALIN" );
        ctm -> SaveCopyOutputNode( "TOWERINFO_CALIB_HCALIN_ORIGINAL" );
        ctm -> RandomizeTowers( true );
        se -> registerSubsystem( ctm );
        PPG04CaloSpy::CaloSpyNodes.push_back( "TOWERINFO_CALIB_HCALIN" );
        PPG04CaloSpy::CaloSpyNodes.push_back( "TOWERINFO_CALIB_HCALIN_ORIGINAL" );

        ctm = new CaloTowerManip( );
        ctm -> Verbosity( verbosity );
        ctm -> SetInputNode( "TOWERINFO_CALIB_HCALOUT" );
        ctm -> SaveCopyOutputNode( "TOWERINFO_CALIB_HCALOUT_ORIGINAL" );
        ctm -> RandomizeTowers( true );
        se -> registerSubsystem( ctm );
        PPG04CaloSpy::CaloSpyNodes.push_back( "TOWERINFO_CALIB_HCALOUT" );
        PPG04CaloSpy::CaloSpyNodes.push_back( "TOWERINFO_CALIB_HCALOUT_ORIGINAL" );
    }

    // Add Calo nodes to AnaWriter
    PPG04Output::CaloNodes.push_back( "TOWERINFO_CALIB_CEMC" );
    PPG04Output::CaloNodes.push_back( "TOWERINFO_CALIB_HCALIN" );
    PPG04Output::CaloNodes.push_back( "TOWERINFO_CALIB_HCALOUT" );
    PPG04Output::CaloNodes.push_back( "TOWERINFO_CALIB_CEMC_RETOWER" );
    PPG04Output::CaloNodes.push_back( "TOWERINFO_CALIB_CEMC_RETOWER_SUB1" );
    PPG04Output::CaloNodes.push_back( "TOWERINFO_CALIB_HCALIN_SUB1" );
    PPG04Output::CaloNodes.push_back( "TOWERINFO_CALIB_HCALOUT_SUB1" );

    if ( PPG04::doIterBackground ) {
        
        auto ijr = new JetReco();
        ijr -> add_input( new TowerJetInput( Jet::CEMC_TOWERINFO_RETOWER, PPG04::TowerPrefix ) );
        ijr -> add_input( new TowerJetInput( Jet::HCALIN_TOWERINFO, PPG04::TowerPrefix ) );
        ijr -> add_input( new TowerJetInput( Jet::HCALOUT_TOWERINFO, PPG04::TowerPrefix ) );
        ijr -> add_algo( new FastJetAlgoSub( Jet::ANTIKT, 0.2 ), "AntiKt_TowerInfo_HIRecoSeedsRaw_r02" );
        ijr -> set_algo_node( "ANTIKT" );
        ijr -> set_input_node( "TOWER" );
        ijr -> Verbosity( Enable::VERBOSITY );
        se -> registerSubsystem( ijr );

        auto dtb = new DetermineTowerBackground( );
        dtb -> SetBackgroundOutputName( "TowerInfoBackground_Sub1" );
        dtb -> SetFlow( false );
        dtb -> SetSeedType( 0 );
        dtb -> SetSeedJetD( 3 );
        dtb -> set_towerinfo( true );
        dtb -> Verbosity( Enable::VERBOSITY ); 
        dtb -> set_towerNodePrefix( PPG04::TowerPrefix );
        se -> registerSubsystem( dtb );

        auto casj = new CopyAndSubtractJets( );
        casj -> SetFlowModulation( false );
        casj -> Verbosity( Enable::VERBOSITY ); 
        casj -> set_towerinfo( true );
        casj -> set_towerNodePrefix( PPG04::TowerPrefix );
        se -> registerSubsystem( casj );

        auto dtb2 = new DetermineTowerBackground( );
        dtb2 -> SetBackgroundOutputName( PPG04::TowerBackgroundNode );
        dtb2 -> SetFlow( false );
        dtb2 -> SetSeedType( 1 );
        dtb2 -> SetSeedJetPt( 7 );
        dtb2 -> Verbosity( Enable::VERBOSITY ); 
        dtb2 -> set_towerinfo( true );
        dtb2 -> set_towerNodePrefix( PPG04::TowerPrefix );
        se -> registerSubsystem( dtb2 );

        auto st = new SubtractTowers( );
        st -> SetFlowModulation( false );
        st -> Verbosity( Enable::VERBOSITY );
        st -> set_towerinfo( true );
        st -> set_towerNodePrefix( PPG04::TowerPrefix );
        se -> registerSubsystem( st );
    }

    for ( auto dtr : PPG04::RhoRecos ) {
        se -> registerSubsystem( dtr );
    }
    for ( auto rct : PPG04::RandomConesRecos ) {
        se -> registerSubsystem( rct );
    }
    for ( auto cwr : PPG04::CaloWindowRecos ) {
        se -> registerSubsystem( cwr );
    }

    std::cout << "PPG04::Running output modules" << std::endl;
    if ( PPG04::doCaloSpy ) {
        std::cout << "CaloSpy file: " << PPG04CaloSpy::outfile << std::endl;
        PPG04::CaloSpyHandler = new CaloSpy( PPG04CaloSpy::outfile );
        PPG04::CaloSpyHandler -> Verbosity( PPG04::VERBOSITY );
        for ( auto node : PPG04CaloSpy::CaloSpyNodes ) { PPG04::CaloSpyHandler -> AddCaloNode( node ); }
        PPG04::CaloSpyHandler -> Normalize( PPG04CaloSpy::Normalize );
        se -> registerSubsystem( PPG04::CaloSpyHandler );
    }

    if ( PPG04::doAnaWriter ){
        std::cout << "Output file: " << PPG04Output::outfile << std::endl; 
        PPG04::AnaWriterHandler = new PPG04AnaWriter(PPG04Output::outfile);
        PPG04::AnaWriterHandler -> Verbosity( PPG04::VERBOSITY );
        if ( PPG04Output::writeMBD ) { PPG04::AnaWriterHandler -> set_mbd_node( true, PPG04::MbdNode ); }
        if ( PPG04Output::writeZVtx ) { PPG04::AnaWriterHandler -> set_zvrtx_node( true, PPG04::ZVrtxNode ); }
        if ( PPG04Output::writeCent ) { PPG04::AnaWriterHandler -> do_cent_node( true, PPG04::CentNode ); }
        if ( PPG04Output::writeIterBackground ) { PPG04::AnaWriterHandler -> do_tower_background_node( true, PPG04::TowerBackgroundNode ); }
        for ( auto node : PPG04Output::CaloNodes ) { PPG04::AnaWriterHandler -> add_calo_node( node ); }
        for ( auto node : PPG04Output::RhoNodes ) { PPG04::AnaWriterHandler -> add_rho_node( node ); }
        for ( auto node : PPG04Output::RandomConesNodes ) { PPG04::AnaWriterHandler -> add_random_cone_node( node ); }
        if ( PPG04Output::doFullWindow ) { PPG04::AnaWriterHandler -> do_calo_window_ana( PPG04Output::FullWindowRetowerNode, PPG04Output::FullWindowHCALINNode, PPG04Output::FullWindowHCALOUTNode ); }
        if ( PPG04Output::doCemcOnlyWindow ) { PPG04::AnaWriterHandler -> do_calo_window_cemc_ana( PPG04Output::CemcOnlyNode ); }
        se -> registerSubsystem( PPG04::AnaWriterHandler );
    }  
    return;

}


#endif