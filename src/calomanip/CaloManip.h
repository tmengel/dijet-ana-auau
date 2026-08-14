#ifndef CaloManip_CaloManip_H
#define CaloManip_CaloManip_H

#include <fun4all/SubsysReco.h>

#include <calobase/RawTowerDefs.h>

#include <string>
#include <vector>

class TFile;
class TTree;

class PHCompositeNode;
class RawTowerGeomContainer;

class CaloManip : public SubsysReco {

 public:

  CaloManip( const std::string & infile )
    : SubsysReco("CaloManip")
    , m_input_file( infile )
  {}
  ~CaloManip() override {}

  int Init( PHCompositeNode * /*topNode*/ ) override;
  int InitRun( PHCompositeNode * topNode ) override;
  int process_event( PHCompositeNode * topNode ) override;
  int End( PHCompositeNode * /*topNode*/ ) override;

  void add_event_header ( const std::string & name = "EventHeader" ) { m_eventhead_node = name; }
  void add_cemc_node    ( const std::string & name = "TOWERINFO_CALIB_CEMC_RETOWER" ) { m_cemc_node = name; }
  void add_hcalin_node  ( const std::string & name = "TOWERINFO_CALIB_HCALIN" ) { m_hcalin_node = name; }
  void add_hcalout_node ( const std::string & name = "TOWERINFO_CALIB_HCALOUT" ) { m_hcalout_node = name; }

  void set_scale_factor ( const float factor ) { m_scale_factor = factor; }

  // Single/multi-event tower-level diagnostic mode. When enabled, for every
  // event_id in m_debug_events this module writes a tower-level TTree
  // ("DebugTowers") and an event-level summary TTree ("DebugEvent") to
  // m_debug_outfile, and prints a numerical sanity-check summary to stdout.
  void set_debug_mode ( const bool b = true ) { m_debug_mode = b; }
  void add_debug_event ( const int event_id ) { m_debug_events.push_back( event_id ); }
  void set_debug_outfile ( const std::string & name ) { m_debug_outfile = name; }

 private:

  bool is_debug_event() const;
  void debug_record_tower( const std::string & calo, int ieta, int iphi, float eta,
                            float e_before, float e_hijing, float e_after );
  void debug_flush_event();

  std::string m_input_file { "" };

  std::string m_eventhead_node { "EventHeader" };

  std::string m_cemc_node { "" };
  std::string m_hcalin_node { "" };
  std::string m_hcalout_node { "" };

  TFile * m_tfile {nullptr};
  TTree * m_ttree {nullptr};

  int m_event_id { -1 };
  int m_nentries_in_tree { -1 };
  int m_ttree_runnumber { -1 };
  int m_ttree_evtsequence { -1 };
  float m_ttree_b { 0.0 };
  float m_ttree_ep_angle { 0.0 };
  float m_ttree_ecc { 0.0 };
  float m_ttree_ncoll { 0.0 };
  float m_ttree_psi1 { 0.0 };
  float m_ttree_psi2 { 0.0 };
  float m_ttree_psi3 { 0.0 };

  // Per-channel (not per ieta,iphi) storage, matching CaloTree's layout --
  // see the note in CaloTree.h. Matching Pass1's reference tower to Pass2's
  // live tower is done purely by channel index (ich), never by re-deriving
  // (ieta,iphi) on one side and hoping it lines up with the other.
  int k_ncemc { -1 };
  int k_nhcalin { -1 };
  int k_nhcalout { -1 };
  std::vector< float > * m_cemc_tower_energy {};
  std::vector< int >   * m_cemc_tower_status {};
  std::vector< float > * m_hcalin_tower_energy {};
  std::vector< int >  * m_hcalin_tower_status {};
  std::vector< float > * m_hcalout_tower_energy {};
  std::vector< int > * m_hcalout_tower_status {};

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

  float m_scale_factor { 1.0 };

  // geometry, used only for the ET diagnostics below (eta comes from the
  // framework geometry node, never hard-coded)
  RawTowerGeomContainer * m_towergeom_cemc {nullptr};
  RawTowerGeomContainer * m_towergeom_hcalin {nullptr};
  RawTowerGeomContainer * m_towergeom_hcalout {nullptr};
  // Calorimeter ID to use when encoding the lookup key into m_towergeom_cemc.
  // When CEMC is retowered onto the IHCal grid, m_towergeom_cemc actually
  // points at TOWERGEOM_HCALIN, whose entries were themselves keyed with
  // RawTowerDefs::HCALIN -- so the *geometry* calo id (HCALIN) differs from
  // the *tower-container* calo id (CEMC) in that case.
  RawTowerDefs::CalorimeterId m_cemc_geom_caloid { RawTowerDefs::CEMC };

  // debug/diagnostic mode
  bool m_debug_mode { false };
  std::vector< int > m_debug_events {};
  std::string m_debug_outfile { "CaloManip_debug.root" };
  TFile * m_debug_tfile {nullptr};
  TTree * m_debug_tower_tree {nullptr};
  TTree * m_debug_event_tree {nullptr};

  // per-tower debug branch variables
  int d_event_id {-1};
  char d_calo[16] {};
  int d_ieta {-1};
  int d_iphi {-1};
  float d_eta {0.0};
  float d_e_before {0.0};
  float d_e_hijing {0.0};
  float d_e_after {0.0};
  float d_scale {0.0};
  float d_et_before {0.0};
  float d_et_after {0.0};

  // per-event summary branch variables ( index 0=CEMC, 1=HCALIN, 2=HCALOUT, 3=TOTAL )
  int e_event_id {-1};
  float e_sumET_hijing[4] {0.0};
  float e_sumET_before[4] {0.0};
  float e_sumET_after[4] {0.0};
  int e_ntowers[4] {0};
  int e_ntowers_modified[4] {0};
  float e_maxAbsDeltaET[4] {0.0};
  int e_maxAbsDeltaET_ieta[4] {-1};
  int e_maxAbsDeltaET_iphi[4] {-1};

  struct DebugRow
  {
    std::string calo;
    int ieta {-1};
    int iphi {-1};
    float eta {0.0};
    float e_before {0.0};
    float e_hijing {0.0};
    float e_after {0.0};
    float et_before {0.0};
    float et_after {0.0};
  };
  std::vector< DebugRow > m_event_debug_rows {};

};

#endif  // CaloManip_CaloManip_H
