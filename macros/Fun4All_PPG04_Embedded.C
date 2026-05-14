#include <GlobalVariables.C>

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


void GetRunSegment( const std::string & f, int & run_number, int & run_segment ) 
{
    std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment( f );
    run_number = runseg.first;
    run_segment = runseg.second;
    return;
}

void Fun4All_PPG04_Embedded( 
    const int nEvents = 10,
    const std::string & outputfile = "dummy.root",
    const std::string & dst_input_list0 = "/sphenix/tg/tg01/jets/bkimelman/embed_Feb24_2025/DST_TRUTH_JET_PPG04_EMBED_Jet30-2024p009-00054912-000480.root",
    const std::string & dst_input_list1 = "/sphenix/tg/tg01/jets/bkimelman/embed_Feb24_2025/DST_GLOBAL_PPG04_EMBED_Jet30-2024p009-00054912-000480.root",
    const std::string & dst_input_list2 = "/sphenix/tg/tg01/jets/bkimelman/embed_Feb24_2025/DST_CALO_PPG04_EMBED_Jet30-2024p009-00054912-000480.root"
)
{

    std::cout << "Starting Fun4All_PPG04" << std::endl;
    
    // Enable
    Enable::VERBOSITY = 0;

    // CDB
    const int timeStamp = 54912;
    CDB::timestamp = static_cast<uint64_t>( timeStamp );
    CDB::global_tag = "2024p009";

    int run_number, run_segment;
    GetRunSegment( dst_input_list0, run_number, run_segment );
    std::cout << "Run number: " << run_number << std::endl;
    std::cout << "Run segment: " << run_segment << std::endl;

    //  PPG04
    PPG04::VERBOSITY = 1;
    PPG04::isDATA = true;
    PPG04::isMC = !PPG04::isDATA;

    // random cones
    PPG04::doRandomCones = false;

    // probes
    PPG04::doJetProbe = false;

    // calo windows
    PPG04::doCaloWindows = false;

    // embedding
    PPG04::doEmbedding = true;
    Embdedding::doTruth = true;
    Embdedding::doSimRetower = true;
    Embdedding::doSim = true;


    // background subtraction
    PPG04::doIterBackground = true;
    PPG04::doAreaRho = true;
    PPG04::doMultRho = true;
   
    // // analysis writer
    PPG04::doAnaWriter = true; 
    PPG04Output::outfile = outputfile;
    PPG04Output::writeMBD = true;
    PPG04Output::writeZVtx = true;
    PPG04Output::writeCent = true;


  
    // Set up F4A
    auto se = Fun4AllServer::instance();
    se -> Verbosity( Enable::VERBOSITY );
    std::vector<std::string> dst_files = { dst_input_list0 , dst_input_list1, dst_input_list2 };

    // set up recoConsts
    auto rc = recoConsts::instance();
    CDBInterface::instance( ) -> Verbosity( 0 );
    rc -> set_StringFlag( "CDB_GLOBALTAG", CDB::global_tag );
    rc -> set_uint64Flag( "TIMESTAMP", CDB::timestamp );
    rc -> set_IntFlag( "PPG04RANDOMSEED", PPG04::PPG04RandomSeed );

    // read in filelists
    for ( unsigned int idx = 0; idx < dst_files.size( ); idx++ ) {

        auto input = new Fun4AllDstInputManager( "DSTINPUT_" + std::to_string( idx ) );
        input -> AddFile( dst_files[ idx ] );
        input -> Verbosity( 1 );
        se -> registerInputManager( input );

    }

    InitPPG04();
    RunPPG04();

    se -> run( nEvents );

    se -> End();   
    std::cout << "Done!" << std::endl;
    gSystem -> Exit( 0 );
}
