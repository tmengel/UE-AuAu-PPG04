#ifndef PPG04_CALO_CALIB_H
#define PPG04_CALO_CALIB_H

#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterDeadHotMask.h>
#include <caloreco/RawClusterPositionCorrection.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <phool/recoConsts.h>

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>  

R__LOAD_LIBRARY( libcalo_reco.so )
R__LOAD_LIBRARY( libffamodules.so )
R__LOAD_LIBRARY( libfun4allutils.so )
R__LOAD_LIBRARY( libcalo_reco.so )


namespace CALOCALIB {

  bool isData = false;
  bool is2024 = false;
 
  unsigned int CalibVersion = 0;
  const unsigned int MaxCalibVersion = 1;

  const char * emcal_feild_name = "Femc_datadriven_qm1_correction";
  const char * emcal_calib_files[] = { 
    "/sphenix/user/egm2153/calib_study/emcal_calib_year1/ana450_2024p009_54912_54921/local_calib_copy_iter15.root", 
    "/sphenix/user/egm2153/calib_study/emcal_calib_year1/54908_54921/local_calib_copy_iter33.root" };
  // alternate for v0: "/sphenix/user/egm2153/calib_study/emcal_calib_year1/ana450_2024p009_54912_54921/local_calib_copy_iter26.root"
  const char * emcal_zscrosscheck_file = "/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/CEMC_ZSCrossCalib_ana450_2024p009_54912.root";

  const char * ihcal_field_name = "HCALIN_calib_ADC_to_ETower";
  const char * ihcal_calib_files[]  = { 
    "/sphenix/u/bseidlitz/work/macros/calibrations/calo/hcal_towerSlope_y2/tsc_cos_comb/AuAuOutput/ihcal_cdb_tsc_cos_calib.root",
    "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/calib_files/HCALIN_calib_ADC_to_ETower_old_mc_digi_scale_54912.root" };
  const char * ihcal_zscrosscheck_file = "/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/HCALIN_ZSCrossCalib_ana450_2024p009_54912.root";

  const char * ohcal_field_name = "HCALOUT_calib_ADC_to_ETower";    
  const char * ohcal_calib_files[]  = { 
    "/sphenix/u/bseidlitz/work/macros/calibrations/calo/hcal_towerSlope_y2/tsc_cos_comb/AuAuOutput/ohcal_cdb_tsc_cos_calib.root",
    "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/calib_files/HCALOUT_calib_ADC_to_ETower_old_mc_digi_scale_54912.root"};
  const char * ohcal_zscrosscheck_file = "/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/HCALOUT_ZSCrossCalib_ana450_2024p009_54912.root";

  // default is set at default values 
  int cemc_software_zs = 60;
  int ohcal_software_zs = 30;
  int ihcal_software_zs = 30;

} // namespace CALOCALIB

void FitTowers() {

  Fun4AllServer *se = Fun4AllServer::instance();
  CaloTowerDefs::BuilderType buildertype = CaloTowerDefs::kPRDFTowerv4;

  /////////////////
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
  ctbEMCal->set_softwarezerosuppression(true, CALOCALIB::cemc_software_zs);
  ctbEMCal->set_bitFlipRecovery(true);
  se->registerSubsystem(ctbEMCal);

  CaloTowerBuilder *ctbIHCal = new CaloTowerBuilder("HCALINBUILDER");
  ctbIHCal->set_detector_type(CaloTowerDefs::HCALIN);
  ctbIHCal->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  ctbIHCal->set_builder_type(buildertype);
  ctbIHCal->set_offlineflag();
  ctbIHCal->set_nsamples(12);
  ctbIHCal->set_softwarezerosuppression(true, CALOCALIB::ihcal_software_zs);
  ctbIHCal->set_bitFlipRecovery(true);
  se->registerSubsystem(ctbIHCal);

  CaloTowerBuilder *ctbOHCal = new CaloTowerBuilder("HCALOUTBUILDER");
  ctbOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
  ctbOHCal->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  ctbOHCal->set_builder_type(buildertype);
  ctbOHCal->set_offlineflag();
  ctbOHCal->set_nsamples(12);
  ctbOHCal->set_softwarezerosuppression(true, CALOCALIB::ohcal_software_zs);
  ctbOHCal->set_bitFlipRecovery(true);
  se->registerSubsystem(ctbOHCal);

  CaloTowerBuilder *caEPD = new CaloTowerBuilder("SEPDBUILDER");
  caEPD->set_detector_type(CaloTowerDefs::SEPD);
  caEPD->set_builder_type(buildertype);
  caEPD->set_processing_type(CaloWaveformProcessing::FAST);
  caEPD->set_nsamples(12);
  caEPD->set_offlineflag();
  se->registerSubsystem(caEPD);

} // FitTowers


void CalibTowers() {
  
  if ( CALOCALIB::CalibVersion > CALOCALIB::MaxCalibVersion ) {
    std::cout << "CalibVersion " << CALOCALIB::CalibVersion << " is not valid. Exiting." << std::endl;
    return;
  }

  Fun4AllServer *se = Fun4AllServer::instance();

  //////////////////////////////
  // Input geometry node
  std::cout << "Adding Geometry file" << std::endl;
  Fun4AllInputManager *ingeo = new Fun4AllRunNodeInputManager("DST_GEO");
  std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
  ingeo->AddFile(geoLocation);
  se->registerInputManager(ingeo);


  //////////////////////////////
  // set statuses on raw towers
  std::cout << "status setters" << std::endl;
  CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
  statusEMC->set_detector_type(CaloTowerDefs::CEMC);
  statusEMC->set_time_cut(1);
  // MC Towers Status
  if( !CALOCALIB::isData ) {
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

  ////////////////////
  // Calibrate towers
  if ( CALOCALIB::is2024 ) {

    std::cout << "Calibrating EMCal" << std::endl;
    CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
    calibEMC->set_detector_type(CaloTowerDefs::CEMC);
    calibEMC->setFieldName(CALOCALIB::emcal_feild_name);
    calibEMC->set_directURL(CALOCALIB::emcal_calib_files[CALOCALIB::CalibVersion]);
    if ( CALOCALIB::isData ) { calibEMC->set_directURL_ZScrosscalib(CALOCALIB::emcal_zscrosscheck_file); }
    se->registerSubsystem(calibEMC);

    std::cout << "Calibrating OHcal" << std::endl;
    CaloTowerCalib *calibOHCal = new CaloTowerCalib("HCALOUT");
    calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
    if ( CALOCALIB::CalibVersion == 1) { calibOHCal->setFieldName(CALOCALIB::ohcal_field_name); }
    calibOHCal->set_directURL(CALOCALIB::ohcal_calib_files[CALOCALIB::CalibVersion]);
    if ( CALOCALIB::isData ) { calibOHCal->set_directURL_ZScrosscalib(CALOCALIB::ohcal_zscrosscheck_file); }
    se->registerSubsystem(calibOHCal);

    std::cout << "Calibrating IHcal" << std::endl;
    CaloTowerCalib *calibIHCal = new CaloTowerCalib("HCALIN");
    calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
    // if ( CALOCALIB::CalibVersion == 1) { calibIHCal->setFieldName(CALOCALIB::ihcal_field_name); }
    calibIHCal->set_directURL(CALOCALIB::ihcal_calib_files[CALOCALIB::CalibVersion]);
    if ( CALOCALIB::isData ) { calibIHCal->set_directURL_ZScrosscalib(CALOCALIB::ihcal_zscrosscheck_file); }
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

} // CalibTowers

#endif // PPG04_CALO_CALIB_H
