#include <GlobalVariables.C>

#include "PPG04_Calo_Calib.C"
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
    const int nEvents = 100,
    const int doRunMode = 0,
    const std::string & outdir = "./",
    const std::string & dst_input_list0 = "filelists/dst_triggered_event_run2auau-00054912.list",
    const std::string & dst_input_list1 = "filelists/dst_calo_waveform.list",
    const std::string & dst_input_list2 = "filelists/dst_mbd_epd.list",
    const std::string & dst_input_list_emb0 = "filelists/dst_calo_cluster_dijet.list",
    const std::string & dst_input_list_emb1 = "filelists/dst_truth_jet_dijet.list",
    const std::string & dst_input_list_emb2 = "filelists/dst_mbd_epd_dijet.list",
    const std::string & dst_input_list_emb3 = "filelists/dst_global_dijet.list"
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

    PPG04::doEmbedding =  false;
    PPG04::isTRUTHJETS = PPG04::doEmbedding && true;
    // Embdedding::doSim = PPG04::doEmbedding && false;
    // Embdedding::doTruth = PPG04::doEmbedding;
    // Embdedding::SrcTOP = "TOPData";
    // Embdedding::TgtTOP = "TOP";
    // Embdedding::TruthJetNode = "AntiKt_Truth_r04";

    // calo calib settings
    CALOCALIB::isData = PPG04::isDATA;
    CALOCALIB::is2024 = true;
    CALOCALIB::cemc_software_zs = 60;
    CALOCALIB::ohcal_software_zs = 30;
    CALOCALIB::ihcal_software_zs = 30;
    CALOCALIB::CalibVersion = 0; //


    // event selection
    PPG04::doEventSelect = true;
    EventSelect::doZVrtxCut = true;
    EventSelect::doMinBiasCut = !PPG04::isMC;
    EventSelect::doTowerChi2Cut = false;
    EventSelect::ZVrtxCutRange = {30,-30};
    
    // background subtraction
    PPG04::doIterBackground = true;
    PPG04::doAreaRho = true;
    PPG04::doMultRho = true;

    // random cones
    PPG04::doRandomCones = true;
    RandomCones::ConeRadius = 0.2;
    RandomCones::ConeAbsEta = 0.8;
    RandomCones::ConeMaskedThreshold = 0.00;

    // probes
    PPG04::doJetProbe = false;// !PPG04::doEmbedding;

    // calo windows
    PPG04::doCaloWindows = true;//  !PPG04::doEmbedding;
   
    // analysis writer
    PPG04::doAnaWriter = true; 
    PPG04Output::outfile = GetOutputFile( mode, prodTag, timeStamp, doRunMode, outdir, dst_input_list0 );
    PPG04Output::writeMBD = true;
    PPG04Output::writeZVtx = true;
    PPG04Output::writeCent = true;
    PPG04Output::writeIterBackground = true;
    PPG04Output::doFullWindow = PPG04::doCaloWindows && false;
    PPG04Output::doCemcOnlyWindow = PPG04::doCaloWindows && false;

    // calo spy
    PPG04::doCaloSpy = false;
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

    
    for ( unsigned int idx = 0; idx < dst_files.size( ); idx++ ) {

        auto input = new Fun4AllDstInputManager( "DSTINPUT_" + std::to_string( idx ) );
        input -> AddListFile( dst_files[ idx ] );
        input -> Verbosity( 0 );
        se -> registerInputManager( input );
    
    }

    if ( PPG04::isDATA ) {
        FitTowers();
    }

    CalibTowers();

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
        // mb->setOverwriteScale("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/scales/cdb_centrality_scale_54912.root"); // will change run by run
        // mb->setOverwriteVtx("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/vertexscales/cdb_centrality_vertex_scale_54912.root"); // will change run by run
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
