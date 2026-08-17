#ifndef _ANATREEV1_H_
#define _ANATREEV1_H_

#include <fun4all/SubsysReco.h>

#include <calobase/RawTowerDefs.h>

#include <string>
#include <vector>
#include <cstdint>

class PHCompositeNode;

class TowerInfoContainer;
class RawTowerGeomContainer;

class TTree;
class TH1D;

class AnaTreev1 : public SubsysReco
{
 public:

  AnaTreev1( const std::string & outputfile = "output.root" );
  ~AnaTreev1() override {}

  int Init( PHCompositeNode */*topNode*/  ) override;
  int InitRun( PHCompositeNode * topNode ) override;
  int process_event( PHCompositeNode * topNode ) override;
  int End( PHCompositeNode * /*topNode*/ ) override;

  void add_gl1_node ( const std::string & name = "GL1Packet"){  m_gl1_node = name; }
  void add_zvrtx_node ( const std::string & name = "GlobalVertexMap" ){ m_zvrtx_node = name; }
  void add_cent_node ( const std::string & name = "CentralityInfo" ) { m_cent_node = name;  }
  void add_mbd_node( const std::string & name = "MbdOut" ) { m_mbd_node = name; }
  void add_minbias_node( const std::string & name = "MinimumBiasInfo" ) { m_minbias_node = name; }

  void add_event_header( const std::string & name =  "EventHeader"){ m_eventhead_node = name; }

  void add_phHep_node ( const std::string & name = "PHHepMCGenEventMap" ){  m_phHep_node = name; }
  void add_truth_jet_node ( const std::string & name ) { m_truth_jet_node = name; }

  void add_sub1_jet_node ( 
    const std::string & name ,
    const std::string & towerbkgd_node 
  ) 
  { 
    m_sub1_jet_node = name; 
    m_sub1_jet_towerbkgd_node = towerbkgd_node;
  }  
  void add_towerbkgd_v2_node( const std::string & name ) { m_towerbkgd_v2_node = name; }

  void add_rho_jet_node ( 
    const std::string & name,
    const std::string & cemcnode ,
    const std::string & hcalincode ,
    const std::string & hcaloutcode 
  )
  {
    m_rho_jet_node = name;
    m_rho_jet_cemc_rho_node = cemcnode;
    m_rho_jet_hcalin_rho_node = hcalincode;
    m_rho_jet_hcalout_rho_node = hcaloutcode;
  }
  void add_rho_nodes ( const std::string & node )
  {
    m_rho_nodes.push_back( node );
  }

  void add_jet_node ( const std::string & node ) { m_jet_node = node; }

  void add_cemc_node ( const std::string & node  , const bool b = true )
  {
    m_cemc_node = node;
    m_save_full_cemc = b;
  }

  void add_ihcal_node ( const std::string & node  , const bool b = true )
  {
    m_ihcal_node = node;
    m_save_full_ihcal = b;
  }

  void add_ohcal_node ( const std::string & node  , const bool b = true )
  {
    m_ohcal_node = node;
    m_save_full_ohcal = b;
  }

 private:
    
  std::string m_output_filename { "" };

  int m_event_id {-1};
  TTree * m_tree { nullptr };

  // TH1D * h_centrality { nullptr };
  // TH1D * h_jet_spectra[100];
  // TH1D * h_jet_spectra_etacut[100];

  std::string m_gl1_node {""};
  static const int k_gl1_max = 64;
  int m_scaled_triggervec[k_gl1_max] {};
  int m_live_triggervec[k_gl1_max] {};

  std::string m_minbias_node { "" };
  int m_is_minbias { 0 };

  std::string m_zvrtx_node { "" };
  float m_zvtx { 0.0 };

  std::string m_cent_node { "" };
  int m_cent {-1};

  std::string m_eventhead_node { "" };
  float m_b { 0.0 };
  float m_ep_angle { 0.0 };
  float m_ecc { 0.0 };
  float m_psi1 { 0.0 };
  float m_psi2 { 0.0 };
  float m_psi3 { 0.0 };
  float m_ncoll { 0.0 };
  float m_npart { 0.0 };
  int m_runnumber { 0 };
  int m_evtsequence { 0 };

  std::string m_mbd_node { "" };
  float m_mbd_q_N { -999.0 };
  float m_mbd_q_S { -999.0 };
  float m_mbd_t_N { -999.0 };
  float m_mbd_t_S { -999.0 };

  std::string m_phHep_node { "" };
  std::string m_truth_jet_node { "" };
  float m_truth_zvtx { 0.0 };
  float m_truth_jet_R { 0.0 };
  float m_truth_jet_maxpT_r04 { 0.0 };
  std::vector < float > m_truth_jet_E {};
  std::vector < float > m_truth_jet_phi {};
  std::vector < float > m_truth_jet_eta {};
  std::vector < float > m_truth_jet_pT {};
  std::vector < int >   m_truth_jet_flavor {};

  static const int k_ieta = 24;
  static const int k_iphi = 64;
  const double m_calo_abs_z[3] = {130.23, 170.299, 301.683};
  const double m_calo_r[3] = {93.5, 127.503, 225.87};
 
  std::string  m_cemc_node {""};
  bool m_save_full_cemc { false };
  float m_sumeT_cemc {0.0};
  float m_cemc_tower_E[k_ieta][k_iphi];
  int m_cemc_tower_isgood[k_ieta][k_iphi];
  
  std::string  m_ihcal_node {""};
  bool m_save_full_ihcal { false };
  float m_sumeT_ihcal {0.0};
  float m_ihcal_tower_E[k_ieta][k_iphi];
  int m_ihcal_tower_isgood[k_ieta][k_iphi];

  std::string  m_ohcal_node {""};
  bool m_save_full_ohcal { false };
  float m_sumeT_ohcal {0.0};
  float m_ohcal_tower_E[k_ieta][k_iphi];
  int m_ohcal_tower_isgood[k_ieta][k_iphi];

  // jet info
  std::string m_jet_node { "" };
  float m_jet_R { 0.0 };
  std::vector < float > m_jet_E {};
  std::vector < float > m_jet_phi {};
  std::vector < float > m_jet_eta {};
  std::vector < float > m_jet_pT {};

  std::string m_sub1_jet_node { "" };
  // particular towerbkgd node for sub1 jet, if not specified, will use m_towerbkgd_node
  std::string m_sub1_jet_towerbkgd_node { "" };
  float m_sub1_jet_towerbkgd_v2 { 0.0 };
  int m_sub1_jet_towerbkgd_flowfaliure { 0 };
  float m_sub1_jet_towerbkgd_psi2 { 0.0 };
  std::vector< std::vector < float > > m_sub1_jet_towerbkgd_ue {};
  float m_sub1_jet_R { 0.0 };
  std::vector < float > m_sub1_jet_E {};
  std::vector < float > m_sub1_jet_phi {};
  std::vector < float > m_sub1_jet_eta {};
  std::vector < float > m_sub1_jet_pT {};
  std::vector < float > m_sub1_jet_unsub_pT {};
  std::vector < float > m_sub1_jet_unsub_E {};
  std::vector < std::vector < float > > m_sub1_jet_constituent_E {};
  std::vector < std::vector < float > > m_sub1_jet_constituent_phi {};
  std::vector < std::vector < float > > m_sub1_jet_constituent_eta {};
  std::vector < std::vector < float > > m_sub1_jet_constituent_pT {};
  std::vector < std::vector < int > > m_sub1_jet_constituent_srcID {};

  // general towerbkgd node for sub2, if not specified, will use m_towerbkgd_node
  std::string m_towerbkgd_v2_node { "" };
  float m_sub2_v2 { 0.0 };
  int m_sub2_flowfaliure { 0 };
  float m_sub2_psi2 { 0.0 };
  std::vector< std::vector < float > > m_sub2_towerbkgd_ue {};

  std::string m_rho_jet_node { "" };
  // particular towerbkgd node for rho jet, if not specified, will use m_towerbkgd_node
  std::string m_rho_jet_cemc_rho_node { "" };
  std::string m_rho_jet_hcalin_rho_node { "" };
  std::string m_rho_jet_hcalout_rho_node { "" };
  float m_rho_jet_cemc_rho { 0.0 };
  float m_rho_jet_hcalin_rho { 0.0 };
  float m_rho_jet_hcalout_rho { 0.0 };
  float m_rho_jet_R { 0.0 };
  std::vector < float > m_rho_jet_E {};
  std::vector < float > m_rho_jet_phi {};
  std::vector < float > m_rho_jet_eta {};
  std::vector < float > m_rho_jet_pT {};
  std::vector < float > m_rho_jet_unsub_pT {};
  std::vector < float > m_rho_jet_unsub_E {};
  std::vector < std::vector < float > > m_rho_jet_constituent_E {};
  std::vector < std::vector < float > > m_rho_jet_constituent_phi {};
  std::vector < std::vector < float > > m_rho_jet_constituent_eta {};
  std::vector < std::vector < float > > m_rho_jet_constituent_pT {};
  std::vector < std::vector < int > > m_rho_jet_constituent_srcID {};

  // general rho
  std::vector< std::string > m_rho_nodes {};
  std::vector< float > m_rho_vals {};
  std::vector< float > m_rho_sigmas {};

  RawTowerDefs::CalorimeterId m_caloid = RawTowerDefs::CalorimeterId::NONE;
  TowerInfoContainer    * m_towerinfos = nullptr;
  RawTowerGeomContainer * m_towergeom  = nullptr;
 
  TowerInfoContainer * LoadTowerInfoContainer( PHCompositeNode *topNode, const std::string & name );
  RawTowerGeomContainer * LoadTowerGeomContainer( PHCompositeNode *topNode, const std::string & name );
};


#endif // _ANATREEV1_H_
