#ifndef CALO_CALIB_H
#define CALO_CALIB_H

#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterDeadHotMask.h>
#include <caloreco/RawClusterPositionCorrection.h>


#include <caloembedding/caloTowerEmbed.h>


#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <phool/recoConsts.h>

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>  // for Fun4AllServer

#include <TSystem.h>  // for gSystem

R__LOAD_LIBRARY( libcalo_reco.so )
R__LOAD_LIBRARY( libffamodules.so )
R__LOAD_LIBRARY( libfun4allutils.so )
R__LOAD_LIBRARY( libCaloEmbedding.so )

namespace CALOCALIB {
  bool EMBED = false;
  bool isData = false;
  std::string EmbdeddingSrcTOP = "TOPData";
  std::string EmbdeddingTgtTOP = "TOP";   
}

void Process_Calo_Calib()
{

  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc = recoConsts::instance();


  // set MC or data
  bool isSim = true;
  int data_sim_runnumber_thres = 1000;
  if (rc->get_uint64Flag("TIMESTAMP") > data_sim_runnumber_thres)
  {
    isSim = false;
  }

  bool is2024 = false;
  int data_run_2024_ths = 50000;
  if (rc->get_uint64Flag("TIMESTAMP") > data_run_2024_ths)
  {
    is2024 = true;
  }
  std::cout << "Calo Calib uses runnumber " << rc->get_uint64Flag("TIMESTAMP") << std::endl;

    // Input geometry node
  std::cout << "Adding Geometry file" << std::endl;
  Fun4AllInputManager *ingeo = new Fun4AllRunNodeInputManager("DST_GEO");
  std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
  ingeo->AddFile(geoLocation);
  se->registerInputManager(ingeo);

  std::cout << "status setters" << std::endl;
  //////////////////////////////
  // set statuses on raw towers
  std::cout << "status setters" << std::endl;
  CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
  statusEMC->set_detector_type(CaloTowerDefs::CEMC);
  statusEMC->set_time_cut(1);
  // MC Towers Status
  if(isSim) {
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

  if ( is2024 ) {
      ////////////////////
  // Calibrate towers
    std::cout << "Calibrating EMCal" << std::endl;
    CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
    calibEMC->set_detector_type(CaloTowerDefs::CEMC);
    calibEMC->setFieldName("Femc_datadriven_qm1_correction");
    calibEMC->set_directURL("/sphenix/user/egm2153/calib_study/emcal_calib_year1/54908_54921/local_calib_copy_iter33.root");
     // calibEMC->set_directURL("/sphenix/user/egm2153/calib_study/emcal_calib_year1/ana450_2024p009_54912_54921/local_calib_copy_iter15.root");
    if ( CALOCALIB::isData ) { 
      calibEMC->set_directURL_ZScrosscalib("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/CEMC_ZSCrossCalib_ana450_2024p009_54912.root");
    }
    se->registerSubsystem(calibEMC);

    std::cout << "Calibrating OHcal" << std::endl;
    CaloTowerCalib *calibOHCal = new CaloTowerCalib("HCALOUT");
    calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
    calibOHCal->set_directURL("/sphenix/u/bseidlitz/work/macros/calibrations/calo/hcal_towerSlope_y2/tsc_cos_comb/AuAuOutput/ohcal_cdb_tsc_cos_calib.root");
    if ( CALOCALIB::isData ) { 
      calibOHCal->set_directURL_ZScrosscalib("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/HCALOUT_ZSCrossCalib_ana450_2024p009_54912.root");
    }
    se->registerSubsystem(calibOHCal);

    std::cout << "Calibrating IHcal" << std::endl;
    CaloTowerCalib *calibIHCal = new CaloTowerCalib("HCALIN");
    calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
    if ( CALOCALIB::isData ) { 
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

  //////////////////
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


  if ( CALOCALIB::EMBED ) {
  
    caloTowerEmbed * embedder_CEMC = new caloTowerEmbed("embedder_CEMC");
    embedder_CEMC->set_detector_type(CaloTowerDefs::CEMC);
    embedder_CEMC->set_inputTop( CALOCALIB::EmbdeddingSrcTOP );
    embedder_CEMC->set_targetTop( CALOCALIB::EmbdeddingTgtTOP );
    embedder_CEMC->set_propgateStatus(true);
    se->registerSubsystem(embedder_CEMC);
    
    caloTowerEmbed * embedder_IHCAL = new caloTowerEmbed("embedder_IHCAL");
    embedder_IHCAL->set_detector_type(CaloTowerDefs::HCALIN);
    embedder_IHCAL->set_inputTop( CALOCALIB::EmbdeddingSrcTOP );
    embedder_IHCAL->set_targetTop( CALOCALIB::EmbdeddingTgtTOP );
    embedder_IHCAL->set_propgateStatus(true);
    se->registerSubsystem(embedder_IHCAL);

    caloTowerEmbed * embedder_OHCAL = new caloTowerEmbed("embedder_OHCal");
    embedder_OHCAL->set_detector_type(CaloTowerDefs::HCALOUT);
    embedder_OHCAL->set_inputTop( CALOCALIB::EmbdeddingSrcTOP );
    embedder_OHCAL->set_targetTop( CALOCALIB::EmbdeddingTgtTOP );
    embedder_OHCAL->set_propgateStatus(true);
    se->registerSubsystem(embedder_OHCAL);
  
  }

  // currently NOT included!
  // std::cout << "Applying Position Dependent Correction" << std::endl;
  // RawClusterPositionCorrection *clusterCorrection = new RawClusterPositionCorrection("CEMC");
  // clusterCorrection->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
  // se->registerSubsystem(clusterCorrection);
}

#endif
