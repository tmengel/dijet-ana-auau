#ifndef _CALOTREE_H_
#define _CALOTREE_H_

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TTree;

class CaloTree : public SubsysReco
{
 public:

  CaloTree( const std::string & outputfile = "output.root" ) 
    : SubsysReco("CaloTree")
    , m_output_filename( outputfile ) {} 
  ~CaloTree() override {}


  int Init( PHCompositeNode * /*topNode*/) override;
  int process_event( PHCompositeNode * topNode ) override;
  int End( PHCompositeNode * /*topNode*/ ) override;
 
  void add_event_header ( const std::string & name = "EventHeader" ) { m_eventhead_node = name; }
  void add_cemc_node    ( const std::string & name = "TOWERINFO_CALIB_CEMC_RETOWER" ) { m_cemc_node = name; }
  void add_hcalin_node  ( const std::string & name = "TOWERINFO_CALIB_HCALIN" ) { m_hcalin_node = name; }
  void add_hcalout_node ( const std::string & name = "TOWERINFO_CALIB_HCALOUT" ) { m_hcalout_node = name; }

  
 private:
    
  std::string m_output_filename { "" };

  std::string m_eventhead_node { "" };

  std::string m_cemc_node { "" };
  std::string m_hcalin_node { "" };
  std::string m_hcalout_node { "" };

  TTree * m_tree {nullptr};

  int m_event_id {-1};


  int k_ncemc { -1 };
  int k_nhcalin { -1 };
  int k_nhcalout { -1 };
  std::vector< float > m_cemc_tower_energy {};
  std::vector< int > m_cemc_tower_status {};
  std::vector< float > m_hcalin_tower_energy {};
  std::vector< int > m_hcalin_tower_status {};
  std::vector< float > m_hcalout_tower_energy {};
  std::vector< int > m_hcalout_tower_status {};

  float m_b { 0.0 };
  float m_ep_angle { 0.0 };
  float m_ecc { 0.0 };
  float m_ncoll { 0.0 };
  float m_npart { 0.0 };
  float m_psi1 { 0.0 };
  float m_psi2 { 0.0 };
  float m_psi3 { 0.0 };
  int m_runnumber { -1 };
  int m_evtsequence { -1 };

};


#endif
