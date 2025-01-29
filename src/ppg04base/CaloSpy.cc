#include "CaloSpy.h"

// fun4all includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

// phool includes
#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// calobase includes
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>  // for TowerInfo
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoDefs.h>

#include <TH2F.h>
#include <TH1F.h>

#include <iostream>

CaloSpy::CaloSpy(const std::string &outputfilename)
  : SubsysReco("CaloSpy")
  , m_output_filename(outputfilename)
{
}


int CaloSpy::Init(PHCompositeNode * /*topNode*/)
{
    // create output file
    PHTFileServer::get().open(m_output_filename, "RECREATE");

    for (unsigned int i = 0; i < m_caloNodes.size(); i++)
    {
      TH2F *h2d = new TH2F(Form("h2d_%s", m_caloNodes[i].c_str()), Form("%s;#eta;#phi", m_caloNodes[i].c_str()),
          m_nphi[i], -0.5, m_nphi[i] - 0.5, m_neta[i], -0.5, m_neta[i] - 0.5);
      h2d->GetYaxis()->SetTitle("#eta");
      h2d->GetXaxis()->SetTitle("#phi");
      h2d->GetZaxis()->SetTitle("Energy [GeV]");
      m_h2d_tower_e_eta_phi.push_back(h2d);

      TH2F *h2dead = new TH2F(Form("h2d_%s_DEAD", m_caloNodes[i].c_str()), Form("%s DEAD;#eta;#phi", m_caloNodes[i].c_str()),
          m_nphi[i], -0.5, m_nphi[i] - 0.5, m_neta[i], -0.5, m_neta[i] - 0.5);
      h2dead->GetYaxis()->SetTitle("#eta");
      h2dead->GetXaxis()->SetTitle("#phi");
      h2dead->GetZaxis()->SetTitle("Energy [GeV]");
      m_h2d_tower_e_eta_phi_dead.push_back(h2dead);

    }

    if (Verbosity() > 0) {
      std::cout << "CaloSpy::Init - Histograms created" << std::endl;
    }

    m_nevents = 0;
    return Fun4AllReturnCodes::EVENT_OK;
}

int CaloSpy::process_event(PHCompositeNode *topNode)
{
  
  if(Verbosity() > 1 ) {
    std::cout << "CaloSpy::process_event - Process event..." << std::endl;
  }
  
  for (unsigned int i = 0; i < m_caloNodes.size(); i++) {
  
    auto towerinfo = findNode::getClass<TowerInfoContainer>( topNode, m_caloNodes[i] );
    if ( ! towerinfo ) {
      std::cout << PHWHERE << " Input node " << m_caloNodes[i] << " Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    unsigned int ntowers = towerinfo->size();
    for ( unsigned int channel = 0; channel < ntowers; channel++ ) {
        auto tower = towerinfo->get_tower_at_channel(channel);
        assert(tower);
        unsigned int key = towerinfo->encode_key(channel);
        float ieta = 1.0*towerinfo->getTowerEtaBin(key);
        float iphi = 1.0*towerinfo->getTowerPhiBin(key);        
        float energy = tower->get_energy();
        if ( tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr() || tower->get_isBadChi2() || std::isnan(energy) ) {
          m_h2d_tower_e_eta_phi_dead[i]->Fill(iphi, ieta);
        } else {
          m_h2d_tower_e_eta_phi[i]->Fill(iphi, ieta, energy);
        }
    }
  }

  m_nevents++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloSpy::End(PHCompositeNode */*topNode*/)
{
  
  if(Verbosity() > 0) {
    std::cout << "CaloSpy::EndRun - End run. Writing to " << m_output_filename << std::endl;
  }

  PHTFileServer::get().cd(m_output_filename);
  if (m_do_norm) {
    float norm = 1.0/m_nevents;
    for (unsigned int i = 0; i < m_caloNodes.size(); i++) {
      m_h2d_tower_e_eta_phi[i]->Scale(norm);
      m_h2d_tower_e_eta_phi_dead[i]->Scale(norm);
    }
  }
  for (unsigned int i = 0; i < m_caloNodes.size(); i++){
    m_h2d_tower_e_eta_phi[i]->Write();
    m_h2d_tower_e_eta_phi_dead[i]->Write();
  }

  if(Verbosity() > 0){
      std::cout << "CaloSpy::EndRun - Done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

