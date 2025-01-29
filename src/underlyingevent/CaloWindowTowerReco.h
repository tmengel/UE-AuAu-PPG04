#ifndef UNDERLYINGEVENT_CALOWINDOWTOWERRECO_H
#define UNDERLYINGEVENT_CALOWINDOWTOWERRECO_H

//===========================================================
/// \file CaloWindowTowerReco.h 
/// \brief SubsysReco module to create a group of calo windows
/// \author Tanner Mengel
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <iostream> 

#include <jetbase/Jet.h>

class PHCompositeNode;

class CaloWindowTowerReco : public SubsysReco
{
  public:
    
    CaloWindowTowerReco(const std::string &name = "CaloWindowTowerReco");
    ~CaloWindowTowerReco() override;

    // standard Fun4All methods
    int Init(PHCompositeNode * topNode) override;
    int process_event(PHCompositeNode * topNode) override;

    void add_input_node( const std::string &name, const std::string &geom_name, Jet::SRC src){ m_inputs.push_back(name); m_geom_names.push_back(geom_name); m_srcs.push_back(src); }
    void add_window_node( const std::string &name){ m_window_names.push_back(name); }

  private:

    std::vector<std::string> m_inputs {}; // input node names
    std::vector<std::string> m_geom_names {}; // geometry node names
    std::vector<Jet::SRC> m_srcs {}; // source of input
    std::vector<std::string> m_window_names {}; // window names

    static const int neta_ihcal = 24;
    static const int neta_emcal = 96;
    static const int nphi_ihcal = 64;
    static const int nphi_emcal = 256;

    int CreateNodes(PHCompositeNode *topNode);

};

#endif // UNDERLYINGEVENT_CALOWINDOWTOWERRECO_H
