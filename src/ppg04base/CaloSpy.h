#ifndef PPG04BASE_CALOSPY_H
#define PPG04BASE_CALOSPY_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <array>
#include <cmath>

class PHCompositeNode;
class TH2;
class TH1;

class CaloSpy : public SubsysReco
{
 public:

   CaloSpy(const std::string &outputfile = "output.root");
   ~CaloSpy() override {};

   void AddCaloNode(const std::string &name, int nphi, int neta) { m_caloNodes.push_back(name); 
      m_nphi.push_back(nphi); 
      m_neta.push_back(neta);
   }

   void Normalize( bool do_norm = true ) { m_do_norm = do_norm; }

   // standard Fun4All functions
   int Init(PHCompositeNode */*topNode*/) override;
   int process_event(PHCompositeNode *topNode) override;
   int End(PHCompositeNode *topNode) override;

 private:
    
   std::string m_output_filename { "output.root" };
   std::vector< std::string > m_caloNodes {};

   bool m_do_norm {false};
   std::vector< TH2 *> m_h2d_tower_e_eta_phi {};
   std::vector< TH2 *> m_h2d_tower_e_eta_phi_dead {};
   std::vector< int > m_nphi {};
   std::vector< int > m_neta {};
   int m_nevents {0};
};


#endif // PPG04BASE_CALOSPY_H
