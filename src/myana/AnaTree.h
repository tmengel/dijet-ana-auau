#ifndef _ANATREE_H_
#define _ANATREE_H_

#include <fun4all/SubsysReco.h>

#include <calobase/RawTowerDefs.h>

#include <string>
#include <vector>
#include <cstdint>

class PHCompositeNode;

class TowerInfoContainer;
class RawTowerGeomContainer;

class TTree;

class AnaTree : public SubsysReco
{
 public:

  AnaTree( const std::string & outputfile = "output.root" );
  ~AnaTree() override {}

  int Init( PHCompositeNode */*topNode*/  ) override;
  int InitRun( PHCompositeNode *topNode ) override;
  int process_event( PHCompositeNode * topNode ) override;
  int End( PHCompositeNode * /*topNode*/ ) override;

  void add_gl1_node ( const std::string & name = "GL1Packet"){  m_gl1_node = name; }
  void add_zvrtx_node ( const std::string & name = "GlobalVertexMap" ){ m_zvrtx_node = name; }
  void add_cent_node ( const std::string & name = "CentralityInfo" ) { m_cent_node = name;  }
  void add_mbd_node( const std::string & name = "MbdOut" ) { m_mbd_node = name; }
  void add_minbias_node( const std::string & name = "MinimumBiasInfo" ) { m_minbias_node = name; }

  void add_sub1jet_node ( const std::string & name ) { m_sub1jet_node = name; }  
  void add_towerbkgd_v2_node( const std::string & name ) { m_towerbkgd_v2_node = name; }

  void add_rhojet_node ( const std::string & name ) { m_rhojet_node = name; }
  void add_rho_nodes ( const std::string & cemcnode, const std::string & hcalincode, const std::string & hcaloutcode , const bool b = true )
  {
    m_isrho_subbed = b;
    m_rho_nodes = { cemcnode, hcalincode, hcaloutcode };
  }
  void add_rho_node( const std::string & name, const std::string & nick_name = "" ) 
  {
    m_single_rho_nodes.push_back( name );
    if ( nick_name != "" ) { m_rho_nicknames.push_back( nick_name ); }
    else { m_rho_nicknames.push_back( name ); }
  }

  void add_arearho_jet_node( const std::string & name ) { m_arearhojet_node = name; }

  void save_full_calo( const std::string & cemcnode, const std::string & hcalincode, const std::string & hcaloutcode, const bool b = true ) 
  { 
    m_save_full_calo = b; 
    m_fullcalo_nodes = { cemcnode, hcalincode, hcaloutcode };
  }
  void save_sum_eT( const std::string & cemcsumeTnode, const std::string & hcalinsumeTnode, const std::string & hcaloutsumeTnode, const bool b = true ) 
  { 
    m_save_sumeT = b; 
    m_sumeT_nodes = { cemcsumeTnode, hcalinsumeTnode, hcaloutsumeTnode };
  }

  void add_truthjet_node ( const std::string & name ) { m_truthjet_node = name; }
  void add_event_header( const std::string & name =  "EventHeader"){ m_eventhead_node = name; }


 private:
    
  // output file name
  std::string m_output_filename { "" };

  TTree * m_tree {nullptr};

  int m_event_id {-1};

  // gl1
  std::string m_gl1_node {""};
  int m_scaled_triggervec[64] {};
  int m_live_triggervec[64] {};

  std::string m_minbias_node { "" };
  int m_is_minbias { 0 };

  // z vertex info
  std::string m_zvrtx_node { "" };
  float m_zvtx { 0.0 };

  // centrality info
  std::string m_cent_node { "" };
  int m_cent {-1};
  
  // calo info
  static const int k_ieta = 24;
  static const int k_iphi = 64;
  bool m_save_full_calo { false };
  std::vector< std::string > m_fullcalo_nodes {};
  float m_tower_E[3][k_ieta][k_iphi] {};
  int   m_tower_isgood[3][k_ieta][k_iphi] {};
  
  bool m_save_sumeT { false };
  std::vector< std::string > m_sumeT_nodes {};
  float m_sumeT[3] { 0.0 }; // 0: cemc, 1: hcalin, 2: hcalout
  std::vector<float> SumCaloE( PHCompositeNode *topNode, const std::vector< std::string > &towerinfo_nodes );
  
  // sub1 jet info
  std::string m_sub1jet_node { "" };
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

  std::string m_towerbkgd_v2_node { "" };
  float m_sub2_v2 { 0.0 };
  int m_sub2_flowfaliure { 0 };
  float m_sub2_psi2 { 0.0 };
  float m_sub2_towerbkgd_ue[3][k_ieta] {}; // 0: cemc, 1: hcalin, 2: hcalout
  
  std::string m_rhojet_node { "" };
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

  std::vector< std::string > m_rho_nodes {};
  bool m_isrho_subbed { false };
  float m_rho_vals[3] {}; // 0: cemc, 1: hcalin, 2: hcalout
  float m_rho_sigmas[3] {}; // 0: cemc, 1: hcalin, 2: hcalout

  std::vector< std::string > m_single_rho_nodes {};
  std::vector< std::string > m_rho_nicknames {};
  float m_single_rho_vals[8] {}; // max of 8
  float m_single_rho_sigmas[8] {}; // max of 8

  std::string m_arearhojet_node { "" };
  std::vector < float > m_arearho_jet_E {};
  std::vector < float > m_arearho_jet_phi {};
  std::vector < float > m_arearho_jet_eta {};
  std::vector < float > m_arearho_jet_pT {};
  std::vector < float > m_arearho_jet_area {};
  float m_arearho_jet_rho;
  std::vector < float > m_arearho_jet_unsub_pT {};
  std::vector < float > m_arearho_jet_unsub_E {};
  std::vector < std::vector < float > > m_arearho_jet_constituent_E {};
  std::vector < std::vector < float > > m_arearho_jet_constituent_phi {};
  std::vector < std::vector < float > > m_arearho_jet_constituent_eta {};
  std::vector < std::vector < float > > m_arearho_jet_constituent_pT {};
  std::vector < std::vector < int > > m_arearho_jet_constituent_srcID {};

  std::string m_truthjet_node { "" };
  std::vector < float > m_truth_jet_E {};
  std::vector < float > m_truth_jet_phi {};
  std::vector < float > m_truth_jet_eta {};
  std::vector < float > m_truth_jet_pT {};

  // event header info
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

  // mbd info
  std::string m_mbd_node { "" };
  float m_mbd_q_N { -999.0 };
  float m_mbd_q_S { -999.0 };
  float m_mbd_t_N { -999.0 };
  float m_mbd_t_S { -999.0 };

  RawTowerDefs::CalorimeterId m_caloid = RawTowerDefs::CalorimeterId::NONE;
  TowerInfoContainer    * m_towerinfos = nullptr;
  RawTowerGeomContainer * m_towergeom  = nullptr;
  double m_calo_r[3] {}; // 0: cemc, 1: hcalin, 2: hcalout
  TowerInfoContainer * LoadTowerInfoContainer( PHCompositeNode *topNode, const std::string &name );
  RawTowerGeomContainer * LoadTowerGeomContainer( PHCompositeNode *topNode, const std::string &name );
};


#endif // _ANATREEV1_H_
