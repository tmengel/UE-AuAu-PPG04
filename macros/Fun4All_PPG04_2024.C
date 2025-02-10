#include <GlobalVariables.C>

#include "Calo_Calib.C"
#include "Calo_Fitting.C"

#include "PPG04.C"

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

#include <mbd/MbdReco.h>

#include <centrality/CentralityReco.h>
#include <g4centrality/PHG4CentralityReco.h>

#include <calotrigger/MinimumBiasClassifier.h>

#include <globalvertex/GlobalVertexReco.h>

#include <zdcinfo/ZdcReco.h>


R__LOAD_LIBRARY( libfun4all.so )
R__LOAD_LIBRARY( libffamodules.so )
R__LOAD_LIBRARY( libmbd.so )
R__LOAD_LIBRARY( libg4mbd.so )
R__LOAD_LIBRARY( libcalotrigger.so )

R__LOAD_LIBRARY( libg4centrality.so )
R__LOAD_LIBRARY( libcentrality.so )
R__LOAD_LIBRARY( libglobalvertex.so )
R__LOAD_LIBRARY( libg4vertex.so )

R__LOAD_LIBRARY( libzdcinfo.so )

void GetRunSegment( const std::string & filelist, int & run_number, int & run_segment ) 
{
    auto ifm = new InputFileHandler( );
    ifm -> AddListFile( filelist );
    auto f = ( ifm -> GetFileList( ) ).front( );
    std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment( f );
    run_number = runseg.first;
    run_segment = runseg.second;
    delete ifm;
    return;
}

std::string GetOutputFile( const std::string & mode,
                    const std::string & prodTag,
                    const int timeStamp,
                    const int doRunMode,
                    const std::string & outdir,
                    const std::string & dst_input_list,
                    const std::string & prefix = "" )
{
    std::ostringstream oss;
    int run_number, run_segment;
    GetRunSegment( dst_input_list, run_number, run_segment );
    oss << outdir;
    if( outdir.back( ) != '/' ) { oss << "/"; }
    if ( prefix != "" ) { 
        oss << prefix; 
        if ( prefix.back( ) != '-' ) { oss << "-"; }
    }
    oss << mode << "-";
    if ( doRunMode == 1 ) { oss << "RandEtaPhi-"; }
    else if ( doRunMode == 2 ) { oss << "Embed-"; }
    oss << prodTag << "-" << std::setw( 6 ) << std::setfill( '0' ) << timeStamp << "_"
        << std::setw( 10 ) << std::setfill( '0' ) << run_number 
        << "-" << std::setw( 6 ) << std::setfill( '0' ) << run_segment << ".root";
    std::string outfile = oss.str( );
    return outfile;
}

void Fun4All_PPG04_2024( 
    const std::string & mode = "DATA",
    const std::string & prodTag = "2024p009",
    const int timeStamp = 54912,
    const int nEvents = 10,
    const int doRunMode = 1,
    const std::string & outdir = "./",
    const std::string & dst_input_list0 = "dst_triggered_event_run2auau-00054912.list",
    const std::string & dst_input_list1 = "dst_calo_waveform.list",
    const std::string & dst_input_list2 = "dst_mbd_epd.list",
    const std::string & dst_input_list_emb0 = "dst_calo_cluster_dijet.list",
    const std::string & dst_input_list_emb1 = "dst_truth_jet_dijet.list",
    const std::string & dst_input_list_emb2 = "dst_mbd_epd_dijet.list",
    const std::string & dst_input_list_emb3 = "dst_global_dijet.list"
)
{

    std::cout << "Starting Fun4All_PPG04" << std::endl;
    
    // Enable
    Enable::VERBOSITY = 0;

    // CDB
    CDB::global_tag = prodTag;
    CDB::timestamp = static_cast<uint64_t>( timeStamp );

    ///---------------------------------------------------------------------------------------------------------------------
    //  PPG04
    PPG04::VERBOSITY = 1;
    PPG04::isDATA = ( mode == "DATA" );
    PPG04::isMC = !PPG04::isDATA;

    PPG04::PPG04RandomSeed = PHRandomSeed();

    // calo manipulation
    PPG04::doCaloManip = ( doRunMode == 1 );
    CaloManip::doTowerRandomizer = PPG04::doCaloManip;
    CaloManip::doMinEMCalEnergy = false;

    PPG04::doEmbedding = ( doRunMode == 2 );
    PPG04::isTRUTHJETS = PPG04::doEmbedding && true;
    Embdedding::doSim = PPG04::doEmbedding && false;
    Embdedding::doTruth = PPG04::doEmbedding;
    Embdedding::SrcTOP = "TOPData";
    Embdedding::TgtTOP = "TOP";
    Embdedding::TruthJetNode = "AntiKt_Truth_r04";

    // calo calib settings
    CALOCALIB::EMBED = PPG04::doEmbedding;
    CALOCALIB::isData = PPG04::isDATA;
    CALOCALIB::EmbdeddingSrcTOP =  Embdedding::SrcTOP;
    CALOCALIB::EmbdeddingTgtTOP = Embdedding::TgtTOP;

    CALOFITTING::cemc_software_zs = 60;
    CALOFITTING::ohcal_software_zs = 30;
    CALOFITTING::ihcal_software_zs = 30;

    // event selection
    PPG04::doEventSelect = true;
    EventSelect::doZVrtxCut = true;
    EventSelect::doMinBiasCut = !PPG04::isMC;
    EventSelect::doTowerChi2Cut = true;
    EventSelect::ZVrtxCutRange = {20,-20};
    
    // background subtraction
    PPG04::doIterBackground = true;
    PPG04::doAreaRho = true;
    PPG04::doMultRho = true;

    // random cones
    PPG04::doRandomCones = !PPG04::doEmbedding;
    RandomCones::ConeRadius = 0.4;
    RandomCones::ConeAbsEta = 1.1;
    RandomCones::ConeMaskedThreshold = 0.05;

    // probes
    PPG04::doJetProbe = !PPG04::doEmbedding;

    // calo windows
    PPG04::doCaloWindows = !PPG04::doEmbedding;
   
    // analysis writer
    PPG04::doAnaWriter = true; 
    PPG04Output::outfile = GetOutputFile( mode, prodTag, timeStamp, doRunMode, outdir, dst_input_list0 );
    PPG04Output::writeMBD = true;
    PPG04Output::writeZVtx = true;
    PPG04Output::writeCent = true;
    PPG04Output::writeIterBackground = true;
    PPG04Output::doFullWindow = PPG04::doCaloWindows && true;
    PPG04Output::doCemcOnlyWindow = PPG04::doCaloWindows && false;

    // calo spy
    PPG04::doCaloSpy = PPG04::doCaloManip && false;
    PPG04CaloSpy::outfile = GetOutputFile( mode, prodTag, timeStamp, doRunMode, outdir, dst_input_list0 , "CALOSPY-");
    PPG04CaloSpy::Normalize = false;
    
    ///---------------------------------------------------------------------------------------------------------------------
    // Set up F4A
    auto se = Fun4AllServer::instance();
    se -> Verbosity( Enable::VERBOSITY );

    std::vector<std::string> dst_files = { dst_input_list0 };
    if ( !PPG04::isDATA ) {
        dst_files.push_back( dst_input_list1 );
        dst_files.push_back( dst_input_list2 );
    }

    // set up recoConsts
    auto rc = recoConsts::instance();
    CDBInterface::instance( ) -> Verbosity( 0 );
    rc -> set_StringFlag( "CDB_GLOBALTAG", CDB::global_tag );
    rc -> set_uint64Flag( "TIMESTAMP", CDB::timestamp );
    rc -> set_IntFlag( "PPG04RANDOMSEED", PPG04::PPG04RandomSeed );

    // read in filelists
    for ( unsigned int idx = 0; idx < dst_files.size( ); idx++ ) {

        auto input = new Fun4AllDstInputManager( "DSTINPUT_" + std::to_string( idx ) );
        input -> AddListFile( dst_files[ idx ] );
        input -> Verbosity( 0 );
        se -> registerInputManager( input );
    
    }
    if ( PPG04::doEmbedding ) {

        auto sync = se->getSyncManager();
        sync->MixRunsOk(true);
        
        auto input_sim = new Fun4AllNoSyncDstInputManager("DSTSimCalo","DST", CALOCALIB::EmbdeddingSrcTOP );
        input_sim->AddListFile( dst_input_list_emb0 );
        input_sim->Verbosity( 0 );
        se->registerInputManager(input_sim);
        
        auto input_truth = new Fun4AllNoSyncDstInputManager("DSTTruthJet","DST", CALOCALIB::EmbdeddingSrcTOP );
        input_truth->AddListFile( dst_input_list_emb1 );
        input_truth->Verbosity( 0 );
        se->registerInputManager(input_truth);
        
        auto input_mbd = new Fun4AllNoSyncDstInputManager("DSTMBD","DST", CALOCALIB::EmbdeddingSrcTOP );
        input_mbd->AddListFile( dst_input_list_emb2 );
        input_mbd->Verbosity( 0 );
        se->registerInputManager(input_mbd);
        
        // auto input_global = new Fun4AllNoSyncDstInputManager("DSTGlobal","DST", CALOCALIB::EmbdeddingSrcTOP );
        // input_global->AddListFile( dst_input_list_emb3 );
        // input_global->Verbosity( 0 );
        // se->registerInputManager(input_global);
    }

    if ( PPG04::isDATA ) {

        Process_Calo_Fitting();
    }

    Process_Calo_Calib();

    if ( PPG04::isDATA ) {
     
        auto mbdreco = new MbdReco();
        se->registerSubsystem( mbdreco );
    
    }

    auto gvertex = new GlobalVertexReco();
    se->registerSubsystem( gvertex );
    
    if ( PPG04::isDATA ) {

        auto zdcreco = new ZdcReco();
        zdcreco->set_zdc1_cut(0.0);
        zdcreco->set_zdc2_cut(0.0);
        se->registerSubsystem( zdcreco );

        auto mb = new MinimumBiasClassifier();
        mb->Verbosity( Enable::VERBOSITY );
        mb->setOverwriteScale("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/scales/cdb_centrality_scale_54912.root"); // will change run by run
        mb->setOverwriteVtx("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/vertexscales/cdb_centrality_vertex_scale_54912.root"); // will change run by run
        se->registerSubsystem( mb );

        auto cent = new CentralityReco();
        cent->setOverwriteScale("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/scales/cdb_centrality_scale_54912.root"); // will change run by run
        cent->setOverwriteVtx("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/vertexscales/cdb_centrality_vertex_scale_54912.root"); // will change run by run
        cent->setOverwriteDivs("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/divs/cdb_centrality_54912.root");
        se->registerSubsystem( cent );
    }

    if ( PPG04::isMC ) {

        
        auto cent = new PHG4CentralityReco();
        cent->Verbosity(Enable::VERBOSITY);
        if ( Enable::CDB ) {
            PHParameterUtils::FillPHParametersFromCDB( cent->GetCalibrationParameters(),"CENTRALITY" );
        } else {
            cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0, string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
        }
        se->registerSubsystem( cent );
    }
      

    InitPPG04();

    RunPPG04();

    se -> run( nEvents );

    se -> End();   
    std::cout << "Done!" << std::endl;
    gSystem -> Exit( 0 );
}
