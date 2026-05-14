#include <GlobalVariables.C>

#include "PPG04_Calo_Calib.C"

#include "PPG04.C"

#include <ffamodules/CDBInterface.h>

#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>

#include <phparameter/PHParameterUtils.h>

// coresoftware headers
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/InputFileHandler.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllSyncManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>

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

    //get ana build tag from evn
    std::string command = "echo $BUILDVER";
    std::string buildTag = "";
    FILE *fp = popen( command.c_str( ), "r" );
    if ( fp == NULL ) {
        std::cerr << "Error: Could not get build tag" << std::endl;
        return "";
    }
    char buffer[128];
    while ( fgets( buffer, sizeof( buffer ), fp ) != NULL ) {
        buildTag += buffer;
    }
    pclose( fp );
    // remove \n from buildTag

    buildTag = buildTag.substr( 0, buildTag.find( "\n" ) );
    // remove . from buildTag
    std::string::iterator end_pos = std::remove( buildTag.begin( ), buildTag.end( ), '.' );
    buildTag.erase( end_pos, buildTag.end( ) );
    //capitalize first letter only
    // std::transform(buildTag.begin( ), buildTag.end( ), buildTag.begin( ), ::toupper );

    std::string prodBuild="ana450";
    oss << outdir;
    if( outdir.back( ) != '/' ) { oss << "/"; }
    if ( prefix != "" ) { 
        oss << prefix; 
        if ( prefix.back( ) != '_' ) { oss << "_"; }
    }
    
    oss << mode << "_";
    if ( doRunMode == 1 ) { oss << "RandEtaPhi_"; }
    oss << "prod." << prodBuild << "_cdb." << prodTag << "_build." << buildTag << "-"
        << std::setw( 10 ) << std::setfill( '0' ) << run_number 
        << "-" << std::setw( 6 ) << std::setfill( '0' ) << run_segment << ".root";
    std::string outfile = oss.str( );
    return outfile;
}

void Fun4All_PPG04_2024Prod( 
  const std::string & mode = "DATA",
    const std::string & prodTag = "2024p009",
    const int timeStamp = 54912,
    const int nEvents = 10,
    const int doRunMode = 0,
    const int zVrtxBin = -1,
    const std::string & outdir = "./",
    const std::string & dst_input_list0 = "dst_triggered_event_run2auau-00054912.list",
    const std::string & dst_input_list1 = "dst_truth_jet_dijet.list",
    const std::string & dst_input_list2 = "dst_mbd_epd_dijet.list",
    const std::string & dst_input_list3 = "dst_global_dijet.list"
)
{

    std::cout << "Starting Fun4All_PPG04" << std::endl;
    
    // Enable
    Enable::VERBOSITY = 0;

    // CDB
    CDB::global_tag = prodTag;
    CDB::timestamp = static_cast<uint64_t>( timeStamp );



    //  PPG04
    PPG04::VERBOSITY = 0;
    PPG04::isDATA = ( mode == "DATA" );
    PPG04::isMC = !PPG04::isDATA;
    PPG04::isTRUTHJETS = ( mode == "PYTHIA" );

    PPG04::PPG04RandomSeed = PHRandomSeed();

    // calo manipulation
    PPG04::doCaloManip = ( doRunMode == 1 );
    CaloManip::doTowerRandomizer = PPG04::doCaloManip;
    CaloManip::doMinEMCalEnergy = false;

    // event selection
    PPG04::doEventSelect = true;

    CALOCALIB::isData = PPG04::isDATA;
    CALOCALIB::is2024 = true;
    CALOCALIB::cemc_software_zs = 60;
    CALOCALIB::ohcal_software_zs = 30;
    CALOCALIB::ihcal_software_zs = 30;
    CALOCALIB::CalibVersion = 0;

    EventSelect::doZVrtxCut = true;
    EventSelect::doMinBiasCut = PPG04::isDATA;
    EventSelect::doTowerChi2Cut = false;

    if (zVrtxBin == 0) { EventSelect::ZVrtxCutRange = {20,15}; } // 1%
    else if (zVrtxBin == 1) { EventSelect::ZVrtxCutRange = {15,10}; } // 7%
    else if (zVrtxBin == 2) { EventSelect::ZVrtxCutRange = {10,5}; } // 22 %
    else if (zVrtxBin == 3) { EventSelect::ZVrtxCutRange = {5,0}; } // 23 %
    else if (zVrtxBin == 4) { EventSelect::ZVrtxCutRange = {0,-5}; } // 29 %
    else if (zVrtxBin == 5) { EventSelect::ZVrtxCutRange = {-5,-10}; } // 16 %
    else if (zVrtxBin == 6) { EventSelect::ZVrtxCutRange = {-10,-15}; } // 16 %
    else if (zVrtxBin == 7) { EventSelect::ZVrtxCutRange = {-15,-20}; }
    else { 
        EventSelect::ZVrtxCutRange = {20,-20}; 
    }
    
    // Set up F4A
    auto se = Fun4AllServer::instance();
    se -> Verbosity( Enable::VERBOSITY );

    std::vector<std::string> dst_files = { dst_input_list0 };
    if ( !PPG04::isDATA ) {
        dst_files.push_back( dst_input_list1 );
        dst_files.push_back( dst_input_list2 );
        dst_files.push_back( dst_input_list3 );
    }

    // set up recoConsts
    auto rc = recoConsts::instance();
    CDBInterface::instance( ) -> Verbosity( 0 );
    rc -> set_StringFlag( "CDB_GLOBALTAG", CDB::global_tag );
    rc -> set_uint64Flag( "TIMESTAMP", CDB::timestamp );
    rc -> set_IntFlag( "PPG04RANDOMSEED", PPG04::PPG04RandomSeed );
    se -> Verbosity( Enable::VERBOSITY );

    // PPG04::doIterBackground = true;
    PPG04::doAreaRho = true;
    PPG04::doMultRho = true;

   
    auto sync = new SyncReco();
    se -> registerSubsystem(sync);
    auto head = new HeadReco();
    se -> registerSubsystem(head);
    auto flag = new FlagHandler();
    se -> registerSubsystem(flag);

    // read in filelists
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
        mb->Verbosity( 0 );
        mb->setOverwriteScale("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/scales/cdb_centrality_scale_54912.root"); // will change run by run
        mb->setOverwriteVtx("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/vertexscales/cdb_centrality_vertex_scale_54912.root"); // will change run by run
        se->registerSubsystem( mb );

        auto cent = new CentralityReco();
        cent->setOverwriteScale("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/scales/cdb_centrality_scale_54912.root"); // will change run by run
        cent->setOverwriteVtx("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/vertexscales/cdb_centrality_vertex_scale_54912.root"); // will change run by run
        cent->setOverwriteDivs("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/divs/cdb_centrality_54912.root");
        se->registerSubsystem( cent );
    }


    InitPPG04();
    RunPPG04();


    std::string FullOutFile =  GetOutputFile("PPG04_"+ mode + "_DSTS", prodTag, timeStamp, doRunMode, outdir, dst_input_list0 );
    Fun4AllOutputManager *out = new Fun4AllDstOutputManager("PPG04OUT_CALOS",FullOutFile);
    se->registerOutputManager(out);
    se -> run( nEvents );
    se -> End();   
    std::cout << "Done!" << std::endl;
    gSystem -> Exit( 0 );
}
