#ifndef MACRO_PPG04_C
#define MACRO_PPG04_C

#include <GlobalVariables.C>

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
#include <jetbackground/TowerConeInput.h>
#include <jetbackground/RandomConeAlgoSub.h>
#include <jetbackground/JetRhoSubtractor.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/DetermineTowerRho.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/TowerRho.h>

#include <fun4all/Fun4AllServer.h>

#include <randomconetree/RandomConeTree.h>

#include <algorithm>
#include <vector>

R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libeventselection.so)
R__LOAD_LIBRARY(libRandomConeTree.so)

namespace Enable
{
    bool IS_DATA = false;
    bool IS_MC = false;
    bool IS_TRUTH_JETS = false;
}  // namespace Enable

namespace PPG04 
{
    int VERBOSITY = 0;
    bool do_event_selection = true;
    bool do_jet_reco = true;
    bool do_truth_jets = false;
    bool do_random_cone = true;

    bool do_iter_sub = true;
    bool do_area_sub = true;
    bool do_mult_sub = true;

    std::string tower_info_prefix = "TOWERINFO_CALIB";

    double emcal_min_e = 0.05;  

    std::string output_file_name = "output.root";
} // namespace PPG04

namespace EventSelection 
{
    int VERBOSITY = 0;
    bool do_zvrtx_cut = true;
    std::pair<float,float> z_vrtx_cut_range = {10,-10};
    std::string zvrtx_node = "GlobalVertexMap";
    
    bool do_min_bias_cut = true;
    std::string min_bias_node = "MinimumBiasInfo";
    
    bool do_truth_jet_cut = true;
    std::pair<float,float> truth_jet_pT_hat_range =  {10, 30};
    std::string truth_jet_node = "AntiKt_Truth_r04";
    
    bool do_tower_chi2_cut = true;
    std::vector<std::string> tower_chi2_nodes = {
        "TOWERINFO_CALIB_CEMC",
        "TOWERINFO_CALIB_HCALIN",
        "TOWERINFO_CALIB_HCALOUT"
    };    

    void SetPtHardRange(float pT_hat)
    {
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

    void SetZertexCutRange(float max, float min = -999.0)
    {
        if (min == -999.0){
            min = -max;
        }
        z_vrtx_cut_range = {max, min};
        return;
    }

} // namespace EventSelection

namespace JetBackground
{
    int VERBOSITY = 0;
    bool did_background = false; // flag to check if background subtraction has been done
    bool do_flow = false;
    float retower_frac_cut = 0.5;
    float seed_jet_d = 3;
    float seed_jet_pt = 7;
    std::string area_rho_node = "TowerRho_AREA";
    std::string mult_rho_node = "TowerRho_MULT";

    std::vector<TowerJetInput*> rho_est_inputs = {
        new TowerJetInput(Jet::CEMC_TOWERINFO, PPG04::tower_info_prefix, PPG04::emcal_min_e),
        new TowerJetInput(Jet::HCALIN_TOWERINFO, PPG04::tower_info_prefix),
        new TowerJetInput(Jet::HCALOUT_TOWERINFO, PPG04::tower_info_prefix)
    };
}  // namespace JetBackground

namespace Jets
{
    int VERBOSITY = 0;
    std::vector<float> jet_params = {0.2, 0.3, 0.4, 0.5};
    void SetJetParams(std::vector<float> params) { jet_params = params; }
    void AddJetParam(float param) { if (std::find(jet_params.begin(), jet_params.end(), param) == jet_params.end()) jet_params.push_back(param); }
}  // namespace Jets

namespace RandomCones
{
    int VERBOSITY = 0;
    bool do_basic = true;
    bool do_rand_eta_phi = true;
    bool do_avoid_lead_jet = true;
    std::string input_node_type = "CALIB";
    std::vector<TowerJetInput*> iter_cone_inputs = {
        new TowerJetInput(Jet::CEMC_TOWERINFO_SUB1, PPG04::tower_info_prefix, PPG04::emcal_min_e),
        new TowerJetInput(Jet::HCALIN_TOWERINFO_SUB1, PPG04::tower_info_prefix),
        new TowerJetInput(Jet::HCALOUT_TOWERINFO_SUB1, PPG04::tower_info_prefix)
    };
    std::vector<TowerJetInput*> rho_cone_inputs = {
        new TowerJetInput(Jet::CEMC_TOWERINFO, PPG04::tower_info_prefix, PPG04::emcal_min_e),
        new TowerJetInput(Jet::HCALIN_TOWERINFO, PPG04::tower_info_prefix),
        new TowerJetInput(Jet::HCALOUT_TOWERINFO, PPG04::tower_info_prefix)
    };
}  // namespace RandomCones

void RunEventSelector()
{

    if ( !PPG04::do_event_selection ) {
        std::cout << "EventSelector: Event selection not enabled. Skipping." << std::endl;
        return ;
    }

    int verbosity = std::max( PPG04::VERBOSITY, EventSelection::VERBOSITY );

    Fun4AllServer * se = Fun4AllServer::instance();

    EventSelector * es = new EventSelector();
    es -> Verbosity( verbosity );
    if ( Enable::IS_DATA && EventSelection::do_min_bias_cut ) {
        MinBiasCut * mbc = new MinBiasCut();
        mbc -> SetNodeName( EventSelection::min_bias_node );
        es -> AddCut( mbc );
    }
    if ( EventSelection::do_tower_chi2_cut ) {
        TowerChi2Cut * tcc = new TowerChi2Cut();
        tcc -> SetNodeNames( EventSelection::tower_chi2_nodes );
        es -> AddCut( tcc );
    }
    if ( EventSelection::do_zvrtx_cut ) {
        ZVertexCut * zvc =  new ZVertexCut( EventSelection::z_vrtx_cut_range.first,
                                            EventSelection::z_vrtx_cut_range.second );
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

    if ( verbosity > 0 ) { es -> PrintCuts(); }

    return ;

}

void RunBackgroundReco()
{

    if ( JetBackground::did_background ) { 
        std::cout << "BackgroundReco: Background subtraction already done. Skipping." << std::endl;
        return; 
    }

    int verbosity = std::max( PPG04::VERBOSITY, JetBackground::VERBOSITY );

    Fun4AllServer * se = Fun4AllServer::instance();
    
    RetowerCEMC * rcemc = new RetowerCEMC(); 
    rcemc -> Verbosity( verbosity ); 
    rcemc -> set_towerinfo( true );
    rcemc -> set_frac_cut( JetBackground::retower_frac_cut ); 
    rcemc -> set_towerNodePrefix( PPG04::tower_info_prefix );
    se -> registerSubsystem( rcemc );

    // iterative background subtraction
    if( PPG04::do_iter_sub ) {

        JetReco * ijr = new JetReco();
        ijr -> add_input( new TowerJetInput( Jet::CEMC_TOWERINFO_RETOWER, PPG04::tower_info_prefix , PPG04::emcal_min_e ) );
        ijr -> add_input( new TowerJetInput( Jet::HCALIN_TOWERINFO, PPG04::tower_info_prefix ) );
        ijr -> add_input( new TowerJetInput( Jet::HCALOUT_TOWERINFO, PPG04::tower_info_prefix ) );
        ijr -> add_algo( new FastJetAlgoSub( Jet::ANTIKT, 0.2 ), "AntiKt_TowerInfo_HIRecoSeedsRaw_r02" );
        ijr -> set_algo_node("ANTIKT");
        ijr -> set_input_node("TOWER");
        ijr -> Verbosity(verbosity);
        se -> registerSubsystem(ijr);

        DetermineTowerBackground * dtb = new DetermineTowerBackground();
        dtb -> SetBackgroundOutputName( "TowerInfoBackground_Sub1" );
        dtb -> SetFlow( JetBackground::do_flow );
        dtb -> SetSeedType( 0 );
        dtb -> SetSeedJetD( JetBackground::seed_jet_d );
        dtb -> set_towerinfo( true );
        dtb -> Verbosity( verbosity ); 
        dtb -> set_towerNodePrefix( PPG04::tower_info_prefix );
        se -> registerSubsystem( dtb );

        CopyAndSubtractJets * casj = new CopyAndSubtractJets();
        casj -> SetFlowModulation( JetBackground::do_flow );
        casj -> Verbosity( verbosity ); 
        casj -> set_towerinfo( true );
        casj -> set_towerNodePrefix( PPG04::tower_info_prefix );
        se -> registerSubsystem( casj );

        DetermineTowerBackground * dtb2 = new DetermineTowerBackground();
        dtb2 -> SetBackgroundOutputName( "TowerInfoBackground_Sub2" );
        dtb2 -> SetFlow( JetBackground::do_flow );
        dtb2 -> SetSeedType( 1 );
        dtb2 -> SetSeedJetPt( JetBackground::seed_jet_pt );
        dtb2 -> set_towerinfo( true );
        dtb2 -> Verbosity( verbosity ); 
        dtb2 -> set_towerNodePrefix( PPG04::tower_info_prefix );
        se -> registerSubsystem( dtb2 );

        SubtractTowers * st = new SubtractTowers();
        st -> SetFlowModulation( JetBackground::do_flow );
        st -> Verbosity( verbosity );
        st -> set_towerinfo( true );
        st -> set_towerNodePrefix( PPG04::tower_info_prefix );
        se -> registerSubsystem( st );

    }
    
    // rho background subtraction
    if( PPG04::do_area_sub || PPG04::do_mult_sub ) {

        DetermineTowerRho * dtr = new DetermineTowerRho();
        if( PPG04::do_area_sub ) { dtr -> add_method( TowerRho::Method::AREA, JetBackground::area_rho_node ); }
        if ( PPG04::do_mult_sub ) { dtr -> add_method( TowerRho::Method::MULT, JetBackground::mult_rho_node ); }
        dtr -> add_tower_input( new TowerJetInput( Jet::CEMC_TOWERINFO, PPG04::tower_info_prefix, PPG04::emcal_min_e ) );
        dtr -> add_tower_input( new TowerJetInput( Jet::HCALIN_TOWERINFO, PPG04::tower_info_prefix ) );
        dtr -> add_tower_input( new TowerJetInput( Jet::HCALOUT_TOWERINFO, PPG04::tower_info_prefix ) );
        // for ( auto input : JetBackground::rho_est_inputs ) { dtr -> add_tower_input( input ); }
        se -> registerSubsystem( dtr );
    
    }
        
    JetBackground::did_background = true;

    return ;
}

void RunJetReco()
{

    if( !PPG04::do_jet_reco ) { 
        std::cout << "JetReco: Jet reconstruction not enabled. Skipping." << std::endl;
        return ; 
    }

    if( !JetBackground::did_background ) { 
        std::cout << "JetReco: Background calculation not done. Calculation now." << std::endl;
        RunBackgroundReco();
    }

    int verbosity = std::max( PPG04::VERBOSITY, Jets::VERBOSITY );

    Fun4AllServer * se = Fun4AllServer::instance();

    if ( Enable::IS_MC && Enable::IS_TRUTH_JETS ) 
    {
        JetReco * tjr = new JetReco();
        TruthJetInput * tji = new TruthJetInput( Jet::PARTICLE );
        tji -> add_embedding_flag( 0 );  // changes depending on signal vs. embedded
        tjr->add_input( tji );
        for ( auto param : Jets::jet_params ) { tjr -> add_algo( new FastJetAlgo( Jet::ANTIKT, param), "AntiKt_Truth_r0" + std::to_string(int(param*10) ) ); }
        tjr -> set_algo_node( "ANTIKT" );
        tjr -> set_input_node( "TRUTH" );
        tjr -> Verbosity( verbosity );
        se -> registerSubsystem( tjr );
    }

    if( PPG04::do_iter_sub ) {
        
        JetReco * ijr  = new JetReco();
        ijr -> add_input( new TowerJetInput( Jet::CEMC_TOWERINFO_SUB1, PPG04::tower_info_prefix ) );
        ijr -> add_input( new TowerJetInput( Jet::HCALIN_TOWERINFO_SUB1, PPG04::tower_info_prefix ) );
        ijr -> add_input( new TowerJetInput( Jet::HCALOUT_TOWERINFO_SUB1, PPG04::tower_info_prefix ) );
        for (auto param : Jets::jet_params) { ijr -> add_algo( new FastJetAlgoSub( Jet::ANTIKT, param, verbosity ), "AntiKt_Tower_r0" + std::to_string(int(param*10)) + "_Sub1" ); }
        ijr -> set_algo_node( "ANTIKT" );
        ijr -> set_input_node( "TOWER" );
        ijr -> Verbosity( verbosity );
        se -> registerSubsystem( ijr );

    }

    if( PPG04::do_area_sub || PPG04::do_mult_sub ) {
        
        JetReco * rjr = new JetReco();
        rjr -> add_input( new TowerJetInput( Jet::CEMC_TOWERINFO_RETOWER, PPG04::tower_info_prefix) );
        rjr -> add_input( new TowerJetInput( Jet::HCALIN_TOWERINFO, PPG04::tower_info_prefix) );
        rjr -> add_input( new TowerJetInput( Jet::HCALOUT_TOWERINFO, PPG04::tower_info_prefix) );
        for ( auto param : Jets::jet_params ) { rjr -> add_algo( new FastJetAlgoSub( Jet::ANTIKT, param, verbosity ),  "AntiKt_Tower_r0" + std::to_string(int(param*10) ) ); }
        rjr -> calc_jet_area(  PPG04::do_area_sub );
        rjr -> set_algo_node( "ANTIKT" );
        rjr -> set_input_node( "TOWER" );
        rjr -> Verbosity( verbosity );
        se->registerSubsystem( rjr );

        JetRhoSubtractor * jrs = new JetRhoSubtractor();
        if ( PPG04::do_area_sub ) { jrs -> add_rho_method( TowerRho::Method::AREA, JetBackground::area_rho_node ); }
        if ( PPG04::do_mult_sub ) { jrs -> add_rho_method( TowerRho::Method::MULT, JetBackground::mult_rho_node ); }
        for ( auto param : Jets::jet_params ) { jrs -> add_input_jets( "AntiKt_Tower_r0" + std::to_string(int(param*10)) ); }
        jrs -> Verbosity( verbosity );
        se -> registerSubsystem( jrs );

    }

    return ;

}

void RunRandomConeReco()
{

    if( !PPG04::do_random_cone ) { 
        std::cout << "RandomConeReco: Random cone not enabled. Skipping." << std::endl;
        return ; 
    }
    if( !JetBackground::did_background ) { 
        std::cout << "JetReco: Background calculation not done. Calculation now." << std::endl;
        RunBackgroundReco();
    }

    int verbosity = std::max( PPG04::VERBOSITY, RandomCones::VERBOSITY );

    Fun4AllServer * se = Fun4AllServer::instance();

    if ( PPG04::do_iter_sub ) {
    
        JetReco * rcr = new JetReco();
        rcr -> add_input( new TowerConeInput( Jet::CEMC_TOWERINFO_SUB1, RandomCones::input_node_type, PPG04::emcal_min_e ) );
        rcr -> add_input( new TowerConeInput( Jet::HCALIN_TOWERINFO_SUB1, RandomCones::input_node_type ) );
        rcr -> add_input( new TowerConeInput( Jet::HCALOUT_TOWERINFO_SUB1, RandomCones::input_node_type ) );
        if( RandomCones::do_basic ) { rcr -> add_algo( new RandomConeAlgoSub( Jet::RANDOMCONE_BASIC, 0.4, verbosity ), "RandomCones_Basic_r04_Sub1" ); }
        if( RandomCones::do_rand_eta_phi ) { rcr -> add_algo( new RandomConeAlgoSub( Jet::RANDOMCONE_RANDETAPHI, 0.4, verbosity ), "RandomCones_RandEtaPhi_r04_Sub1" ); }
        if( RandomCones::do_avoid_lead_jet ) { rcr -> add_algo( new RandomConeAlgoSub( Jet::RANDOMCONE_AVOIDLEADJET, 0.4, verbosity ), "RandomCones_AvoidLeadJet_r04_Sub1" ); }
        rcr -> set_algo_node( "RANDOMCONE" );
        rcr -> set_input_node( "CONES" );
        rcr -> Verbosity( verbosity );
        se -> registerSubsystem( rcr );
    }
    
    if ( PPG04::do_area_sub || PPG04::do_mult_sub ) {

        JetReco * rcr = new JetReco();
        rcr -> add_input( new TowerConeInput( Jet::CEMC_TOWERINFO, RandomCones::input_node_type, PPG04::emcal_min_e ) );
        rcr -> add_input( new TowerConeInput( Jet::HCALIN_TOWERINFO, RandomCones::input_node_type ) );
        rcr -> add_input( new TowerConeInput( Jet::HCALOUT_TOWERINFO, RandomCones::input_node_type ) );
        if( RandomCones::do_basic ) { rcr -> add_algo( new RandomConeAlgoSub( Jet::RANDOMCONE_BASIC, 0.4, verbosity ), "RandomCones_Basic_r04" ); }
        if( RandomCones::do_rand_eta_phi ) { rcr -> add_algo( new RandomConeAlgoSub( Jet::RANDOMCONE_RANDETAPHI, 0.4, verbosity ), "RandomCones_RandEtaPhi_r04" ); }
        if( RandomCones::do_avoid_lead_jet ) { rcr -> add_algo( new RandomConeAlgoSub( Jet::RANDOMCONE_AVOIDLEADJET, 0.4, verbosity ), "RandomCones_AvoidLeadJet_r04" ); }
        rcr -> set_algo_node( "RANDOMCONE" );
        rcr -> set_input_node( "CONES" );
        rcr -> Verbosity( verbosity );
        se -> registerSubsystem( rcr );
        
        // JetRhoSubtractor * jrs = new JetRhoSubtractor();
        // if ( PPG04::do_area_sub ) { jrs -> add_rho_method( TowerRho::Method::AREA, JetBackground::area_rho_node ); }
        // if ( PPG04::do_mult_sub ) { jrs -> add_rho_method( TowerRho::Method::MULT, JetBackground::mult_rho_node ); }
        // if( RandomCones::do_basic ) { jrs -> add_input_jets( "RandomCones_Basic_r04" ); }
        // if( RandomCones::do_rand_eta_phi ) { jrs -> add_input_jets( "RandomCones_RandEtaPhi_r04" ); }
        // if( RandomCones::do_avoid_lead_jet ) { jrs -> add_input_jets( "RandomCones_AvoidLeadJet_r04" ); }
        // jrs -> Verbosity( verbosity );
        // se -> registerSubsystem( jrs );

    }

    return ;

}

void RunPPG04()
{
    RunEventSelector();
    RunBackgroundReco();
    RunJetReco();
    RunRandomConeReco();


    Fun4AllServer * se = Fun4AllServer::instance();
    RandomConeTree * rct = new RandomConeTree( PPG04::output_file_name );
    rct -> Verbosity( PPG04::VERBOSITY );
    if ( PPG04::do_area_sub ) { rct -> add_rho_node( JetBackground::area_rho_node ); }
    if ( PPG04::do_mult_sub ) { rct -> add_rho_node( JetBackground::mult_rho_node ); }
    if ( PPG04::do_random_cone ) { 
        if ( RandomCones::do_basic ) {
            if ( PPG04::do_iter_sub ) { rct -> add_random_cone_node( "RandomCones_Basic_r04_Sub1" ); }
            // if ( PPG04::do_area_sub ) { rct -> add_random_cone_node( "RandomCones_Basic_r04_SubRhoAREA" ); }
            // if ( PPG04::do_mult_sub ) { rct -> add_random_cone_node( "RandomCones_Basic_r04_SubRhoMULT" ); }
        }
        if ( RandomCones::do_rand_eta_phi ) {
            if ( PPG04::do_iter_sub ) { rct -> add_random_cone_node( "RandomCones_RandEtaPhi_r04_Sub1" ); }
        //     if ( PPG04::do_area_sub ) { rct -> add_random_cone_node( "RandomCones_RandEtaPhi_r04_SubRhoAREA" ); }
        //     if ( PPG04::do_mult_sub ) { rct -> add_random_cone_node( "RandomCones_RandEtaPhi_r04_SubRhoMULT" ); }
        }
        // if ( RandomCones::do_avoid_lead_jet ) {
        //     if ( PPG04::do_iter_sub ) { rct -> add_random_cone_node( "RandomCones_AvoidLeadJet_r04_Sub1" ); }
        //     if ( PPG04::do_area_sub ) { rct -> add_random_cone_node( "RandomCones_AvoidLeadJet_r04_SubRhoAREA" ); }
        //     if ( PPG04::do_mult_sub ) { rct -> add_random_cone_node( "RandomCones_AvoidLeadJet_r04_SubRhoMULT" ); }
        // }
        // add raw jets
        if ( PPG04::do_iter_sub ) { rct -> add_random_cone_node( "RandomCones_Basic_r04" ); }
        if ( PPG04::do_area_sub ) { rct -> add_random_cone_node( "RandomCones_RandEtaPhi_r04" ); }

    }
    se -> registerSubsystem( rct );

    return ;
}



#endif // MACRO_PPG04_C