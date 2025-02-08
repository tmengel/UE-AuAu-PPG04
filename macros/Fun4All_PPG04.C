#include <G4_Centrality.C>
#include <G4_Global.C>
#include <GlobalVariables.C>
#include <Calo_Calib.C> 

#include "PPG04.C"

// coresoftware headers
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

// phool headers
#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>

// fun4all headers
#include <fun4all/InputFileHandler.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <mbd/MbdReco.h>

// load libraries
R__LOAD_LIBRARY( libfun4all.so )
R__LOAD_LIBRARY( libffamodules.so )
R__LOAD_LIBRARY( libmbd.so )
R__LOAD_LIBRARY( libcalotrigger.so )

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
                    const int doRandomEtaPhi,
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
    if ( doRandomEtaPhi == 1 ) { oss << "RandEtaPhi-"; }
    oss << prodTag << "-" << std::setw( 6 ) << std::setfill( '0' ) << timeStamp << "_"
        << std::setw( 10 ) << std::setfill( '0' ) << run_number 
        << "-" << std::setw( 6 ) << std::setfill( '0' ) << run_segment << ".root";
    std::string outfile = oss.str( );
    return outfile;
}

void Fun4All_PPG04( 
    const std::string & mode = "HIJING",
    const std::string & prodTag = "ProdA_2023",
    const int timeStamp = 23745,
    const int nEvents = 10,
    const int doRandomEtaPhi = 0,
    const std::string & outdir = "./",
    const std::string & dst_input_list0 = "dst_calo_cluster.list",
    const std::string & dst_input_list1 = "dst_calo_waveform.list",
    const std::string & dst_input_list2 = "dst_mbd_epd.list"
)
{

    std::cout << "Starting Fun4All_PPG04" << std::endl;
    // Enable
    Enable::VERBOSITY = 0;
    Enable::CENTRALITY_VERBOSITY = 0;
    Enable::DSTOUT = false;
    // Enable::CDB = true;

    // CDB
    CDB::global_tag = prodTag;
    CDB::timestamp = static_cast<uint64_t>( timeStamp );

    ///---------------------------------------------------------------------------------------------------------------------
    //  PPG04
    PPG04::VERBOSITY = 1;
    PPG04::isDATA = ( mode == "DATA" );
    PPG04::isMC = !PPG04::isDATA;
    PPG04::isTRUTHJETS = false;
    PPG04::PPG04RandomSeed = PHRandomSeed();
    // analysis writer
    PPG04::doAnaWriter = true; 
    PPG04Output::outfile = GetOutputFile( mode, prodTag, timeStamp, doRandomEtaPhi, outdir, dst_input_list0 );
    PPG04Output::writeMBD = true;
    PPG04Output::writeZVtx = true;
    PPG04Output::writeCent = true;
    PPG04Output::writeIterBackground = true;
    PPG04Output::doFullWindow = true;
    PPG04Output::doCemcOnlyWindow = true;    
      // calo spy
    PPG04::doCaloSpy = true;
    PPG04CaloSpy::outfile = GetOutputFile( mode, prodTag, timeStamp, doRandomEtaPhi, outdir, dst_input_list0 , "CALOSPY-");
    PPG04CaloSpy::Normalize = false;
    // event selection
    PPG04::doEventSelect = true;
    EventSelect::doZVrtxCut = true;
    EventSelect::doMinBiasCut = true;
    EventSelect::doTowerChi2Cut = true;
    EventSelect::ZVrtxCutRange = {20,-20};
    // calo manipulation
    PPG04::doCaloManip = true;
    CaloManip::doMinEMCalEnergy = false;
    // CaloManip::MinEMCalEnergy = 0.05; // nominal
    CaloManip::MinEMCalEnergy = 0.150; // for 0.15 GeV
    CaloManip::doTowerRandomizer = ( doRandomEtaPhi == 1 );
    // background subtraction
    PPG04::doIterBackground = true;
    PPG04::doAreaRho = true;
    PPG04::doMultRho = true;
    // random cones
    PPG04::doRandomCones = true;
    // calo windows
    PPG04::doCaloWindows = true;
    
    ///---------------------------------------------------------------------------------------------------------------------
    // Set up F4A
    auto se = Fun4AllServer::instance( );
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

    // auto sync = new SyncReco( );
    // se -> registerSubsystem( sync );
    // auto head = new HeadReco( );
    // se -> registerSubsystem( head );
    // auto flag = new FlagHandler( );
    // se -> registerSubsystem( flag );

    // read in filelists
    for ( unsigned int idx = 0; idx < dst_files.size( ); idx++ ) {
        auto input = new Fun4AllDstInputManager( "DSTINPUT_" + std::to_string( idx ) );
        input -> AddListFile( dst_files[ idx ] );
        input -> Verbosity( 0 );
        se -> registerInputManager( input );
    }

    ///---------------------------------------------------------------------------------------------------------------------
    // Run4All
    Global_Reco();
    Process_Calo_Calib();
    if ( !PPG04::isDATA  ) { Centrality(); } 
    InitPPG04();
    RunPPG04();
    se -> run( nEvents );

    ///---------------------------------------------------------------------------------------------------------------------
    se -> End();   
    std::cout << "Done!" << std::endl;
    gSystem -> Exit( 0 );
}
