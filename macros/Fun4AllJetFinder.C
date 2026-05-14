#include <GlobalVariables.C>
#include <HIJetReco.C>


#include <ffamodules/CDBInterface.h>

#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>

#include <phparameter/PHParameterUtils.h>

#include <fun4all/InputFileHandler.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllSyncManager.h>


#include <eventselection/EventSelector.h>
#include <eventselection/ZVertexCut.h>

#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>

#include <jetbackground/RetowerCEMC.h>

#include <jetwriter/JetWriter.h>

R__LOAD_LIBRARY( libfun4all.so )
R__LOAD_LIBRARY( libffamodules.so )
R__LOAD_LIBRARY( libeventselection.so )
R__LOAD_LIBRARY( libjetbackground.so )
R__LOAD_LIBRARY( libjetbase.so )
R__LOAD_LIBRARY( libjetwriter.so )


void Fun4AllJetFinder( 
    const int nEvents = 10,
    const std::string & data_type = "calo",
    const std::string & outputfile = "dummy.root",
    const std::string & dst_input_list0 = "filelists/dst_calo_run2pp-00047478.list",
    const std::string & dst_input_list1 = "filelists/dst_calo_run2pp-00047479.list"
)
{

    std::cout << "Starting Fun4All_PPG04" << std::endl;
    
    // Enable
    Enable::VERBOSITY = 0;

    // CDB
    const std::string & mode = "DATA";
    const std::string & prodTag = "2024prodA";
    const int timeStamp = 47478;

    CDB::global_tag = prodTag;
    CDB::timestamp = static_cast<uint64_t>( timeStamp );

 
    ///---------------------------------------------------------------------------------------------------------------------
    // Set up F4A
    auto se = Fun4AllServer::instance();
    se -> Verbosity( Enable::VERBOSITY );

    std::vector<std::string> dst_files = { dst_input_list0 };
    if ( data_type == "jetcalo"){
        dst_files.push_back( dst_input_list1 );
    }
   

    // set up recoConsts
    auto rc = recoConsts::instance();
    CDBInterface::instance( ) -> Verbosity( 0 );
    rc -> set_StringFlag( "CDB_GLOBALTAG", CDB::global_tag );
    rc -> set_uint64Flag( "TIMESTAMP", CDB::timestamp );
    
    for ( unsigned int idx = 0; idx < dst_files.size( ); idx++ ) {

        auto input = new Fun4AllDstInputManager( "DSTINPUT_" + std::to_string( idx ) );
        input -> AddListFile( dst_files[ idx ] );
        input -> Verbosity( 0 );
        se -> registerInputManager( input );
    
    }


    // Event selection
    auto es = new EventSelector( "EventSelector" );
    auto zvc =  new ZVertexCut( 30, -30 );
    zvc -> SetNodeName( "GlobalVertexMap" );
    zvc -> Verbosity( 0 );
    es -> AddCut( zvc );
    es -> PrintCuts();
    se -> registerSubsystem( es );

    if ( data_type != "calo"  ){

        RetowerCEMC *rcemc = new RetowerCEMC(); 
        rcemc->Verbosity(0); 
        rcemc->set_towerinfo(true);
        rcemc->set_frac_cut(0.5); //fraction of retower that must be masked to mask the full retower
        rcemc->set_towerNodePrefix("TOWERINFO_CALIB" );
        se->registerSubsystem(rcemc);

        // RunPPG04();
        auto jr = new JetReco();
        TowerJetInput *incemc = new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER,"TOWERINFO_CALIB");
        TowerJetInput *inihcal = new TowerJetInput(Jet::HCALIN_TOWERINFO,"TOWERINFO_CALIB");
        TowerJetInput *inohcal = new TowerJetInput(Jet::HCALOUT_TOWERINFO,"TOWERINFO_CALIB");
        if (HIJETS::do_vertex_type)
        {
            incemc->set_GlobalVertexType(HIJETS::vertex_type);
            inihcal->set_GlobalVertexType(HIJETS::vertex_type);
            inohcal->set_GlobalVertexType(HIJETS::vertex_type);
        }
        jr -> add_algo( new FastJetAlgoSub( Jet::ANTIKT, 0.4 ), "AntiKt_Tower_r04" );
        jr -> set_algo_node( "ANTIKT" );
        jr -> set_input_node( "TOWER" );
        jr -> Verbosity( Enable::VERBOSITY );
        se -> registerSubsystem( jr );

    }

    JetWriter *jw = new JetWriter( outputfile );
    jw -> set_jet_node( "AntiKt_Tower_r04" );
    jw -> set_min_energy_frac( 0.3 );
    jw -> set_min_dphi( 0.75 * M_PI );
    jw -> set_min_eT( 10.0 );
    jw -> Verbosity( Enable::VERBOSITY );
    se -> registerSubsystem( jw );
    
    se -> run( nEvents );
    se -> End();   

    std::cout << "Done!" << std::endl;
    gSystem -> Exit( 0 );
}
