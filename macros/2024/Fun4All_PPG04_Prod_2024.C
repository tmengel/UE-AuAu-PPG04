#include <GlobalVariables.C>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>

#include <phparameter/PHParameterUtils.h>

#include <fun4all/InputFileHandler.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllSyncManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>

#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterDeadHotMask.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/RawClusterPositionCorrection.h>

#include <eventselection/EventSelector.h>
#include <eventselection/MinBiasCut.h>
#include <eventselection/TowerChi2Cut.h>
#include <eventselection/ZVertexCut.h>
#include <eventselection/LeadJetCut.h>

#include <g4jets/TruthJetInput.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbase/JetProbeInput.h>
#include <jetbase/JetProbeMaker.h>

#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/DetermineTowerRho.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/TowerRho.h>

#include <ppg04base/CaloTowerManip.h>

#include <mbd/MbdReco.h>

#include <centrality/CentralityReco.h>
#include <g4centrality/PHG4CentralityReco.h>

#include <calotrigger/MinimumBiasClassifier.h>

#include <globalvertex/GlobalVertexReco.h>

#include <zdcinfo/ZdcReco.h>


R__LOAD_LIBRARY( libfun4all.so )
R__LOAD_LIBRARY( libffamodules.so )
R__LOAD_LIBRARY( libfun4allutils.so )
R__LOAD_LIBRARY( libmbd.so )
R__LOAD_LIBRARY( libg4mbd.so )
R__LOAD_LIBRARY( libcalotrigger.so )
R__LOAD_LIBRARY( libg4centrality.so )
R__LOAD_LIBRARY( libcentrality.so )
R__LOAD_LIBRARY( libglobalvertex.so )
R__LOAD_LIBRARY( libg4vertex.so )
R__LOAD_LIBRARY( libzdcinfo.so )
R__LOAD_LIBRARY( libcalo_reco.so )
R__LOAD_LIBRARY( libeventselection.so )
R__LOAD_LIBRARY( libppg04base.so )
R__LOAD_LIBRARY( libjetbackground.so )
R__LOAD_LIBRARY( libjetbase.so )



namespace PPG04 {

    int VERBOSITY = 0;

    bool isDATA = false;
    bool isMC = false;
    bool isTRUTHJETS = false;

    int PPG04RandomSeed = 42;

    bool doEventSelect = false;
    bool doCaloManip = false;

    bool doIterBackground = true;
    bool doAreaRho = false;
    bool doMultRho = false;

    bool doRandomCones = false;
     
    std::string MbdNode = "MbdOut";
    std::string MinBiasNode = "MinimumBiasInfo";
    std::string ZVrtxNode = "GlobalVertexMap";
    std::string CentNode = "CentralityInfo";
    std::string TowerBackgroundNode = "TowerInfoBackground_Sub2";

    EventSelector * EventSelectorHandler {nullptr};

    std::string TowerPrefix = "TOWERINFO_CALIB";
    
} // namespace PPG04

namespace OUTPUTMANAGER{
  set<string> outfiles;
} // namespace OUTPUTMANAGER

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

namespace CALOFITTING {
  int cemc_software_zs = 60;
  int ohcal_software_zs = 30;
  int ihcal_software_zs = 30;
}// namespace CALOFITTING


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
    oss << "DST_PPG04-" << mode << "-";
    // oss << mode << "-";
    if ( doRunMode == 1 ) { oss << "RandEtaPhi-"; }
    else if ( doRunMode == 2 ) { oss << "Embed-"; }
    oss << prodTag << "-" << std::setw( 6 ) << std::setfill( '0' ) << timeStamp << "_"
        << std::setw( 10 ) << std::setfill( '0' ) << run_number 
        << "-" << std::setw( 6 ) << std::setfill( '0' ) << run_segment << ".root";
    std::string outfile = oss.str( );
    return outfile;
}

void Fun4All_PPG04_Prod_2024( 
    const std::string & mode = "DATA",
    const std::string & prodTag = "2024p009",
    const int timeStamp = 54912,
    const int nEvents = 10,
    const int doRunMode = 0,
    const std::string & outdir = "./",
    const std::string & dst_input_list0 = "dst_triggered_event_run2auau-00054912.list",
    const std::string & dst_input_list1 = "dst_calo_waveform.list",
    const std::string & dst_input_list2 = "dst_mbd_epd.list",
    const std::string & dst_input_list3 = "dst_global_dijet.list"
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
    PPG04::isTRUTHJETS = false;

    PPG04::PPG04RandomSeed = PHRandomSeed();

    // calo manipulation
    PPG04::doCaloManip = ( doRunMode == 1 );
    CaloManip::doTowerRandomizer = PPG04::doCaloManip;
    CaloManip::doMinEMCalEnergy = false;

    // event selection
    PPG04::doEventSelect = true;
    EventSelect::doZVrtxCut = true;
    EventSelect::doMinBiasCut = !PPG04::isMC;
    EventSelect::doTowerChi2Cut = true;
    EventSelect::ZVrtxCutRange = {20,-20};
    
    // background subtraction
    PPG04::doIterBackground = true;
    PPG04::doAreaRho = false;
    PPG04::doMultRho = false;

    // random cones
    PPG04::doRandomCones = false;
    RandomCones::ConeRadius = 0.4;
    RandomCones::ConeAbsEta = 1.1;
    RandomCones::ConeMaskedThreshold = 0.05;

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

    // read in filelists
    for ( unsigned int idx = 0; idx < dst_files.size( ); idx++ ) {

        auto input = new Fun4AllDstInputManager( "DSTINPUT_" + std::to_string( idx ) );
        input -> AddListFile( dst_files[ idx ] );
        input -> Verbosity( 0 );
        se -> registerInputManager( input );
    
    }

    CALOFITTING::cemc_software_zs = 60;
    CALOFITTING::ohcal_software_zs = 30;
    CALOFITTING::ihcal_software_zs = 30;
    if ( PPG04::isDATA ) {

        CaloTowerDefs::BuilderType buildertype = CaloTowerDefs::kPRDFTowerv4;

        // build towers
        CaloTowerBuilder *caZDC = new CaloTowerBuilder("ZDCBUILDER");
        caZDC->set_detector_type(CaloTowerDefs::ZDC);
        caZDC->set_builder_type(buildertype);
        caZDC->set_processing_type(CaloWaveformProcessing::FAST);
        caZDC->set_nsamples(16);
        caZDC->set_offlineflag();
        se->registerSubsystem(caZDC);

        CaloTowerBuilder *ctbEMCal = new CaloTowerBuilder("EMCalBUILDER");
        ctbEMCal->set_detector_type(CaloTowerDefs::CEMC);
        ctbEMCal->set_processing_type(CaloWaveformProcessing::TEMPLATE);
        ctbEMCal->set_builder_type(buildertype);
        ctbEMCal->set_offlineflag(true);
        ctbEMCal->set_nsamples(12);
        ctbEMCal->set_softwarezerosuppression(true, 100);
        ctbEMCal->set_bitFlipRecovery(true);
        se->registerSubsystem(ctbEMCal);

        CaloTowerBuilder *ctbIHCal = new CaloTowerBuilder("HCALINBUILDER");
        ctbIHCal->set_detector_type(CaloTowerDefs::HCALIN);
        ctbIHCal->set_processing_type(CaloWaveformProcessing::TEMPLATE);
        ctbIHCal->set_builder_type(buildertype);
        ctbIHCal->set_offlineflag();
        ctbIHCal->set_nsamples(12);
        ctbIHCal->set_softwarezerosuppression(true, 50);
        ctbIHCal->set_bitFlipRecovery(true);
        se->registerSubsystem(ctbIHCal);

        CaloTowerBuilder *ctbOHCal = new CaloTowerBuilder("HCALOUTBUILDER");
        ctbOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
        ctbOHCal->set_processing_type(CaloWaveformProcessing::TEMPLATE);
        ctbOHCal->set_builder_type(buildertype);
        ctbOHCal->set_offlineflag();
        ctbOHCal->set_nsamples(12);
        ctbOHCal->set_softwarezerosuppression(true, 50);
        ctbOHCal->set_bitFlipRecovery(true);
        se->registerSubsystem(ctbOHCal);

        CaloTowerBuilder *caEPD = new CaloTowerBuilder("SEPDBUILDER");
        caEPD->set_detector_type(CaloTowerDefs::SEPD);
        caEPD->set_builder_type(buildertype);
        caEPD->set_processing_type(CaloWaveformProcessing::FAST);
        caEPD->set_nsamples(12);
        caEPD->set_offlineflag();
        se->registerSubsystem(caEPD);
    }


    bool isSimCalib = true;
    int data_sim_runnumber_thres = 1000;
    if (rc->get_uint64Flag("TIMESTAMP") > data_sim_runnumber_thres) {
        isSimCalib = false;
    }
    std::cout << "Calo Calib uses runnumber " << rc->get_uint64Flag("TIMESTAMP") << std::endl;

    // Input geometry node
    std::cout << "Adding Geometry file" << std::endl;
    Fun4AllInputManager *ingeo = new Fun4AllRunNodeInputManager("DST_GEO");
    std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
    ingeo->AddFile(geoLocation);
    se->registerInputManager(ingeo);

    std::cout << "status setters" << std::endl;
    CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
    statusEMC->set_detector_type(CaloTowerDefs::CEMC);
    statusEMC->set_time_cut(1);
    // MC Towers Status
    if( isSimCalib ) {
        // Uses threshold of 50% for towers be considered frequently bad.
        std::string calibName_hotMap = "CEMC_hotTowers_status";
        /* Systematic options (to be used as needed). */
        /* Uses threshold of 40% for towers be considered frequently bad. */
        // std::string calibName_hotMap = "CEMC_hotTowers_status_40";

        /* Uses threshold of 60% for towers be considered frequently bad. */
        // std::string calibName_hotMap = "CEMC_hotTowers_status_60";

        std::string calibdir = CDBInterface::instance()->getUrl("calibName_hotMap");
        statusEMC->set_directURL_hotMap(calibdir);
    }
    se->registerSubsystem(statusEMC);

    CaloTowerStatus *statusHCalIn = new CaloTowerStatus("HCALINSTATUS");
    statusHCalIn->set_detector_type(CaloTowerDefs::HCALIN);
    statusHCalIn->set_time_cut(2);
    se->registerSubsystem(statusHCalIn);

    CaloTowerStatus *statusHCALOUT = new CaloTowerStatus("HCALOUTSTATUS");
    statusHCALOUT->set_detector_type(CaloTowerDefs::HCALOUT);
    statusHCALOUT->set_time_cut(2);
    se->registerSubsystem(statusHCALOUT);

    if ( !isSimCalib ) {

        // Calibrate towers
        std::cout << "Calibrating EMCal" << std::endl;
        CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
        calibEMC->set_detector_type(CaloTowerDefs::CEMC);
        calibEMC->setFieldName("Femc_datadriven_qm1_correction");
        calibEMC->set_directURL("/sphenix/user/egm2153/calib_study/emcal_calib_year1/54908_54921/local_calib_copy_iter33.root");
        // calibEMC->set_directURL("/sphenix/user/egm2153/calib_study/emcal_calib_year1/ana450_2024p009_54912_54921/local_calib_copy_iter15.root");
        if ( PPG04::isDATA ) { 
        calibEMC->set_directURL_ZScrosscalib("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/CEMC_ZSCrossCalib_ana450_2024p009_54912.root");
        }
        se->registerSubsystem(calibEMC);

        std::cout << "Calibrating OHcal" << std::endl;
        CaloTowerCalib *calibOHCal = new CaloTowerCalib("HCALOUT");
        calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
        calibOHCal->set_directURL("/sphenix/u/bseidlitz/work/macros/calibrations/calo/hcal_towerSlope_y2/tsc_cos_comb/AuAuOutput/ohcal_cdb_tsc_cos_calib.root");
        if ( PPG04::isDATA ) { 
        calibOHCal->set_directURL_ZScrosscalib("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/HCALOUT_ZSCrossCalib_ana450_2024p009_54912.root");
        }
        se->registerSubsystem(calibOHCal);

        std::cout << "Calibrating IHcal" << std::endl;
        CaloTowerCalib *calibIHCal = new CaloTowerCalib("HCALIN");
        calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
        if ( PPG04::isDATA ) { 
        calibIHCal->set_directURL("/sphenix/u/bseidlitz/work/macros/calibrations/calo/hcal_towerSlope_y2/tsc_cos_comb/AuAuOutput/ihcal_cdb_tsc_cos_calib.root");
        }
        calibIHCal->set_directURL_ZScrosscalib("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/HCALIN_ZSCrossCalib_ana450_2024p009_54912.root");
        se->registerSubsystem(calibIHCal);

    } else {
        // Calibrate towers
        std::cout << "Calibrating EMCal" << std::endl;
        CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
        calibEMC->set_detector_type(CaloTowerDefs::CEMC);
        se->registerSubsystem(calibEMC);

        std::cout << "Calibrating OHcal" << std::endl;
        CaloTowerCalib *calibOHCal = new CaloTowerCalib("HCALOUT");
        calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
        se->registerSubsystem(calibOHCal);

        std::cout << "Calibrating IHcal" << std::endl;
        CaloTowerCalib *calibIHCal = new CaloTowerCalib("HCALIN");
        calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
        se->registerSubsystem(calibIHCal);
    } 

    // Clusters
    std::cout << "Building clusters" << std::endl;
    RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
    ClusterBuilder->Detector("CEMC");
    ClusterBuilder->set_threshold_energy(0.070);  // for when using basic calibration
    std::string emc_prof = getenv("CALIBRATIONROOT");
    emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
    ClusterBuilder->LoadProfile(emc_prof);
    ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
    se->registerSubsystem(ClusterBuilder);

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
    } else { 

        auto cent = new PHG4CentralityReco();
        cent->Verbosity(Enable::VERBOSITY);
        if ( Enable::CDB ) {
            PHParameterUtils::FillPHParametersFromCDB( cent->GetCalibrationParameters(),"CENTRALITY" );
        } else {
            cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0, string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
        }
        se->registerSubsystem( cent );
    }
    
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
        if ( CaloManip::doTowerRandomizer ) {
            ctm -> RandomizeTowers( true );
        }
        se -> registerSubsystem( ctm );
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

        ctm = new CaloTowerManip( );
        ctm -> Verbosity( verbosity );
        ctm -> SetInputNode( "TOWERINFO_CALIB_HCALIN" );
        ctm -> SaveCopyOutputNode( "TOWERINFO_CALIB_HCALIN_ORIGINAL" );
        ctm -> RandomizeTowers( true );
        se -> registerSubsystem( ctm );

        ctm = new CaloTowerManip( );
        ctm -> Verbosity( verbosity );
        ctm -> SetInputNode( "TOWERINFO_CALIB_HCALOUT" );
        ctm -> SaveCopyOutputNode( "TOWERINFO_CALIB_HCALOUT_ORIGINAL" );
        ctm -> RandomizeTowers( true );
        se -> registerSubsystem( ctm );
    }

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


    // std::string FullOutFile = PPG04Output::outfile;
    // Fun4AllOutputManager *out = new Fun4AllDstOutputManager("PPG04OUT",FullOutFile);
    // out->AddNode("Sync");
    // out->AddNode("EventHeader");
    // out->AddNode("TOWERINFO_CALIB_CEMC");
    // out->AddNode("TOWERINFO_CALIB_CEMC_RETOWER");
    // out->AddNode("TOWERINFO_CALIB_HCALIN");
    // out->AddNode("TOWERINFO_CALIB_HCALOUT");
    // out->AddNode("GlobalVertexMap");
    // out->AddNode("MbdVertexMap");
    // out->AddNode("MinimumBiasInfo");
    // out->AddNode("Zdcinfo");
    // out->AddNode("SEPD");
    // out->AddNode("CentralityInfo");
    // out->AddNode("MbdOut");
    // se->registerOutputManager(out);
    // OUTPUTMANAGER::outfiles.insert(FullOutFile);
    // se->registerOutputManager(out);


    se -> run( nEvents );

    se -> End();   
    std::cout << "Done!" << std::endl;
    gSystem -> Exit( 0 );
}
