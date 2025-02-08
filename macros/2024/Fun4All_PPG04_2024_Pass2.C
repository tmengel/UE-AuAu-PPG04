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




void Fun4All_PPG04_2024_Pass2( 
    const std::string & mode = "DATA",
    const std::string & prodTag = "2024p009",
    const int timeStamp = 54912,
    const int nEvents = -1,
    const int doRunMode = 0,
    const std::string & outdir = "./",
    const std::string & dst_input_list0 = "/sphenix/user/tmengel/UE-AuAu-PPG04/dsts/DATA/ZVTX_BIN3/BASIC/CALOSDATA_CALO_DSTS-2024p009-054912_0000054912-000001.root"
    // const std::string & dst_input_list0 = "CALOSDATA_CALO_DSTS-2024p009-054912_0000054912-000000.root"
)
{

    std::cout << "Starting Fun4All_PPG04" << std::endl;
    
    // Enable
    Enable::VERBOSITY = 1;

    // CDB
    CDB::global_tag = prodTag;
    CDB::timestamp = static_cast<uint64_t>( timeStamp );

    //  PPG04
    PPG04::VERBOSITY = 1;
    PPG04::isDATA = ( mode == "DATA" );
    PPG04::isMC = !PPG04::isDATA;

    // random cones
    PPG04::doRandomCones = true;

    // probes
    PPG04::doJetProbe = true;

    // calo windows
    PPG04::doCaloWindows = true;
   
    // // analysis writer
    // PPG04::doAnaWriter = true; 
    // PPG04Output::outfile = "dummy.root";
    // PPG04Output::writeMBD = true;
    // PPG04Output::writeZVtx = true;
    // PPG04Output::writeCent = true;
    // PPG04Output::writeIterBackground = true;
    // PPG04Output::doFullWindow = PPG04::doCaloWindows && true;
  
    // Set up F4A
    auto se = Fun4AllServer::instance();
    se -> Verbosity( Enable::VERBOSITY );
    auto sync = new SyncReco();
    se -> registerSubsystem(sync);
    auto head = new HeadReco();
    se -> registerSubsystem(head);
    auto flag = new FlagHandler();
    se -> registerSubsystem(flag);
    std::vector<std::string> dst_files = { dst_input_list0  };

    // set up recoConsts
    auto rc = recoConsts::instance();
    CDBInterface::instance( ) -> Verbosity( 0 );
    rc -> set_StringFlag( "CDB_GLOBALTAG", CDB::global_tag );
    rc -> set_uint64Flag( "TIMESTAMP", CDB::timestamp );
    rc -> set_IntFlag( "PPG04RANDOMSEED", PPG04::PPG04RandomSeed );
    // read in filelists
    for ( unsigned int idx = 0; idx < dst_files.size( ); idx++ ) {

        auto input = new Fun4AllDstInputManager( "DSTINPUT_" + std::to_string( idx ) );
        // input -> AddListFile( dst_files[ idx ] );
        input -> AddFile( dst_files[ idx ] );
        input -> Verbosity( 1 );
        se -> registerInputManager( input );

        // auto sync = input->TotalEvents(  );
        }

    
    auto gvertex = new GlobalVertexReco();
    se->registerSubsystem( gvertex );
    // InitPPG04();
    // RunPPG04();

    se -> run( -1 );

    se -> End();   
    std::cout << "Done!" << std::endl;
    gSystem -> Exit( 0 );
}
