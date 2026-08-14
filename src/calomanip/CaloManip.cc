#include "CaloManip.h"

#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <phool/getClass.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerDefs.h>

#include <ffaobjects/EventHeaderv1.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

int CaloManip::Init( PHCompositeNode * /*topNode*/ )
{
  
  m_tfile = TFile::Open( m_input_file.c_str(), "READ" );
  if ( !m_tfile || m_tfile->IsZombie() ) 
  {
    std::cout << "CaloManip::Init - Fatal Error - cannot open input tree file " << m_input_file << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_ttree = dynamic_cast<TTree*>( m_tfile->Get("T") );
  if ( !m_ttree ) 
  {
    std::cout << "CaloManip::Init - Fatal Error - cannot find tree T in input file " << m_input_file << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (Verbosity() > 2 )
  {
    m_ttree -> Show(0);
  }

  m_ttree -> SetBranchStatus( "*", false );
  m_ttree -> SetBranchStatus( "runnumber", true );
  m_ttree -> SetBranchStatus( "evtsequence", true );
  m_ttree -> SetBranchStatus( "b", true );
  m_ttree -> SetBranchStatus( "ep_angle", true );
  m_ttree -> SetBranchStatus( "ecc", true );
  m_ttree -> SetBranchStatus( "ncoll", true );
  m_ttree -> SetBranchStatus( "psi1", true );
  m_ttree -> SetBranchStatus( "psi2", true );
  m_ttree -> SetBranchStatus( "psi3", true );
  m_ttree -> SetBranchStatus( "ncemc", true );
  m_ttree -> SetBranchStatus( "nhcalin", true );
  m_ttree -> SetBranchStatus( "nhcalout", true );
  m_ttree -> SetBranchStatus( "cemc_tower_energy", true );
  m_ttree -> SetBranchStatus( "cemc_tower_status", true );
  m_ttree -> SetBranchStatus( "hcalin_tower_energy", true );
  m_ttree -> SetBranchStatus( "hcalin_tower_status", true );
  m_ttree -> SetBranchStatus( "hcalout_tower_energy", true );
  m_ttree -> SetBranchStatus( "hcalout_tower_status", true );
  m_ttree -> SetBranchAddress( "ncemc", &k_ncemc );
  m_ttree -> SetBranchAddress( "nhcalin", &k_nhcalin );
  m_ttree -> SetBranchAddress( "nhcalout", &k_nhcalout );
  m_ttree -> SetBranchAddress( "runnumber", &m_ttree_runnumber );
  m_ttree -> SetBranchAddress( "evtsequence", &m_ttree_evtsequence );
  m_ttree -> SetBranchAddress( "b", &m_ttree_b );
  m_ttree -> SetBranchAddress( "ep_angle", &m_ttree_ep_angle );
  m_ttree -> SetBranchAddress( "ecc", &m_ttree_ecc );
  m_ttree -> SetBranchAddress( "ncoll", &m_ttree_ncoll );
  m_ttree -> SetBranchAddress( "psi1", &m_ttree_psi1 );
  m_ttree -> SetBranchAddress( "psi2", &m_ttree_psi2 );
  m_ttree -> SetBranchAddress( "psi3", &m_ttree_psi3 );
  m_ttree -> SetBranchAddress( "cemc_tower_energy", &m_cemc_tower_energy );
  m_ttree -> SetBranchAddress( "cemc_tower_status", &m_cemc_tower_status );
  m_ttree -> SetBranchAddress( "hcalin_tower_energy", &m_hcalin_tower_energy );
  m_ttree -> SetBranchAddress( "hcalin_tower_status", &m_hcalin_tower_status );
  m_ttree -> SetBranchAddress( "hcalout_tower_energy", &m_hcalout_tower_energy );
  m_ttree -> SetBranchAddress( "hcalout_tower_status", &m_hcalout_tower_status );

  m_event_id = -1;
  m_nentries_in_tree =  m_ttree -> GetEntries();
  if ( Verbosity() > 0 )
  {
    std::cout << "CaloManip: opened input tree file " << m_input_file << " with " << m_nentries_in_tree << " entries" << std::endl;
  }

  if ( m_debug_mode )
  {
    m_debug_tfile = new TFile( m_debug_outfile.c_str(), "RECREATE" );

    m_debug_tower_tree = new TTree( "DebugTowers", "per-tower scaling diagnostics" );
    m_debug_tower_tree -> Branch( "event_id", &d_event_id, "event_id/I" );
    m_debug_tower_tree -> Branch( "calo", d_calo, "calo/C" );
    m_debug_tower_tree -> Branch( "ieta", &d_ieta, "ieta/I" );
    m_debug_tower_tree -> Branch( "iphi", &d_iphi, "iphi/I" );
    m_debug_tower_tree -> Branch( "eta", &d_eta, "eta/F" );
    m_debug_tower_tree -> Branch( "e_before", &d_e_before, "e_before/F" );
    m_debug_tower_tree -> Branch( "e_hijing", &d_e_hijing, "e_hijing/F" );
    m_debug_tower_tree -> Branch( "e_after", &d_e_after, "e_after/F" );
    m_debug_tower_tree -> Branch( "scale_factor", &d_scale, "scale_factor/F" );
    m_debug_tower_tree -> Branch( "et_before", &d_et_before, "et_before/F" );
    m_debug_tower_tree -> Branch( "et_after", &d_et_after, "et_after/F" );

    m_debug_event_tree = new TTree( "DebugEvent", "event-level SigmaET diagnostics" );
    m_debug_event_tree -> Branch( "event_id", &e_event_id, "event_id/I" );
    m_debug_event_tree -> Branch( "sumET_hijing", e_sumET_hijing, "sumET_hijing[4]/F" );
    m_debug_event_tree -> Branch( "sumET_before", e_sumET_before, "sumET_before[4]/F" );
    m_debug_event_tree -> Branch( "sumET_after", e_sumET_after, "sumET_after[4]/F" );
    m_debug_event_tree -> Branch( "ntowers", e_ntowers, "ntowers[4]/I" );
    m_debug_event_tree -> Branch( "ntowers_modified", e_ntowers_modified, "ntowers_modified[4]/I" );
    m_debug_event_tree -> Branch( "maxAbsDeltaET", e_maxAbsDeltaET, "maxAbsDeltaET[4]/F" );
    m_debug_event_tree -> Branch( "maxAbsDeltaET_ieta", e_maxAbsDeltaET_ieta, "maxAbsDeltaET_ieta[4]/I" );
    m_debug_event_tree -> Branch( "maxAbsDeltaET_iphi", e_maxAbsDeltaET_iphi, "maxAbsDeltaET_iphi[4]/I" );

    std::cout << "CaloManip::Init - debug mode enabled, writing diagnostics for "
              << m_debug_events.size() << " event(s) to " << m_debug_outfile << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

int CaloManip::InitRun( PHCompositeNode * topNode )
{
  if ( m_debug_mode )
  {
    // Geometry for the ET diagnostics, read from the RUN node (populated by
    // the DST_GEO Fun4AllRunNodeInputManager), never hard-coded. The
    // retowered CEMC container shares the IHCal (HCALIN) 24x64 projective
    // grid (see jetbackground/RetowerCEMC), so its eta must be looked up
    // from TOWERGEOM_HCALIN rather than the native TOWERGEOM_CEMC.
    const bool is_retowered = ( std::string::npos != m_cemc_node.find("RETOWER") );
    const std::string cemc_geom_node = is_retowered ? "TOWERGEOM_HCALIN" : "TOWERGEOM_CEMC";
    m_cemc_geom_caloid = is_retowered ? RawTowerDefs::HCALIN : RawTowerDefs::CEMC;
    m_towergeom_cemc    = findNode::getClass<RawTowerGeomContainer>( topNode, cemc_geom_node );
    m_towergeom_hcalin  = findNode::getClass<RawTowerGeomContainer>( topNode, "TOWERGEOM_HCALIN" );
    m_towergeom_hcalout = findNode::getClass<RawTowerGeomContainer>( topNode, "TOWERGEOM_HCALOUT" );

    if ( !m_towergeom_cemc || !m_towergeom_hcalin || !m_towergeom_hcalout )
    {
      std::cout << PHWHERE << " CaloManip::InitRun - debug mode requested but tower geometry node(s) missing "
                << "( " << cemc_geom_node << "=" << m_towergeom_cemc
                << ", TOWERGEOM_HCALIN=" << m_towergeom_hcalin
                << ", TOWERGEOM_HCALOUT=" << m_towergeom_hcalout << " ), aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloManip::process_event( PHCompositeNode * topNode )
{

  m_event_id++; 

  if( Verbosity() > 1 ) 
  {
    std::cout << "CaloManip::process_event - Process event " << m_event_id << std::endl;
  }

  auto * eventhead = findNode::getClass<EventHeader>( topNode, m_eventhead_node );
  if ( !eventhead )
  {
    std::cout << PHWHERE << " Input node " << m_eventhead_node << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_runnumber = eventhead->get_RunNumber();
  m_evtsequence = eventhead->get_EvtSequence();
  m_b = eventhead->get_ImpactParameter();
  m_ep_angle = eventhead->get_EventPlaneAngle();
  m_ecc = eventhead->get_eccentricity();
  m_ncoll = eventhead->get_ncoll();
  m_npart = eventhead->get_npart();
  m_psi1 = eventhead->get_FlowPsiN(1);
  m_psi2 = eventhead->get_FlowPsiN(2);
  m_psi3 = eventhead->get_FlowPsiN(3);
  if ( Verbosity() > 1 ) 
  {
    std::cout << PHWHERE << " - runnumber = " << m_runnumber << ", b = " << m_b << ", ep_angle = " << m_ep_angle << ", ecc = " << m_ecc << ", psi2 = " << m_psi2 << ", ncoll = " << m_ncoll << ", npart = " << m_npart << std::endl;
  }

  auto * cemc_towers = findNode::getClass<TowerInfoContainer>( topNode, m_cemc_node );
  if ( !cemc_towers )
  {
    std::cout << PHWHERE << " Input node " << m_cemc_node << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  auto * hcalin_towers = findNode::getClass<TowerInfoContainer>( topNode, m_hcalin_node );
  if ( !hcalin_towers )
  {
    std::cout << PHWHERE << " Input node " << m_hcalin_node << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  auto * hcalout_towers = findNode::getClass<TowerInfoContainer>( topNode, m_hcalout_node );
  if ( !hcalout_towers )
  {
    std::cout << PHWHERE << " Input node " << m_hcalout_node << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if ( m_event_id >= m_nentries_in_tree )
  {
    std::cout << PHWHERE << " Event ID " << m_event_id << " exceeds number of entries in tree " << m_nentries_in_tree << ", skipping. " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_ttree -> GetEntry( m_event_id );
  // Event matching: the HIJING-only reference for event m_event_id must
  // belong to the *same* underlying HIJING event as the HIJING+PYTHIA8
  // event currently on topNode. Entry-number correspondence alone is not
  // proof of that (a dropped/reordered event in either pass would silently
  // desynchronize the two streams), so we additionally require an exact
  // match of evtsequence (the framework's per-event sequence identifier)
  // and, as a robustness cross-check, the HIJING generator-level event
  // observables (b, ep_angle, ecc, ncoll, psi1-3) carried in EventHeader.
  if ( m_runnumber != m_ttree_runnumber
    || m_evtsequence != m_ttree_evtsequence
    || m_b != m_ttree_b
    || m_ep_angle != m_ttree_ep_angle
    || m_ecc != m_ttree_ecc
    || m_ncoll != m_ttree_ncoll
    || m_psi1 != m_ttree_psi1
    || m_psi2 != m_ttree_psi2
    || m_psi3 != m_ttree_psi3
  )
  {
    std::cout << PHWHERE << " FATAL: Event ID " << m_event_id
              << " HIJING+PYTHIA8 event does not match the HIJING-only reference tree entry -- "
              << "event correspondence between the two passes cannot be demonstrated, skipping. " << std::endl;
    std::cout << "  runnumber: " << m_runnumber << " vs " << m_ttree_runnumber << std::endl;
    std::cout << "  evtsequence: " << m_evtsequence << " vs " << m_ttree_evtsequence << std::endl;
    std::cout << "  b: " << m_b << " vs " << m_ttree_b << std::endl;
    std::cout << "  ep_angle: " << m_ep_angle << " vs " << m_ttree_ep_angle << std::endl;
    std::cout << "  ecc: " << m_ecc << " vs " << m_ttree_ecc << std::endl;
    std::cout << "  ncoll: " << m_ncoll << " vs " << m_ttree_ncoll << std::endl;
    std::cout << "  psi1: " << m_psi1 << " vs " << m_ttree_psi1 << std::endl;
    std::cout << "  psi2: " << m_psi2 << " vs " << m_ttree_psi2 << std::endl;
    std::cout << "  psi3: " << m_psi3 << " vs " << m_ttree_psi3 << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  const bool debug_event = is_debug_event();
  if ( debug_event )
  {
    m_event_debug_rows.clear();
  }

  // Channel-count consistency check: Pass1's reference tree and the live
  // embedded container must have been built from the SAME node type (same
  // retowered-vs-native grid) to be matched channel-for-channel below. A
  // mismatch here (e.g. one side native 96x256, the other retowered 24x64)
  // would silently misassociate towers if we matched by (ieta,iphi) instead
  // -- so we assert on channel count and refuse to guess.
  if ( (int) cemc_towers->size() != k_ncemc
    || (int) hcalin_towers->size() != k_nhcalin
    || (int) hcalout_towers->size() != k_nhcalout )
  {
    std::cout << PHWHERE << " FATAL: channel count mismatch between the live embedded towers and the "
              << "HIJING-only reference tree at event " << m_event_id
              << " -- CEMC: " << cemc_towers->size() << " vs " << k_ncemc
              << ", HCALIN: " << hcalin_towers->size() << " vs " << k_nhcalin
              << ", HCALOUT: " << hcalout_towers->size() << " vs " << k_nhcalout
              << ". This means the two passes used different tower grids (e.g. retowered vs native CEMC) "
              << "-- refusing to scale towers that cannot be reliably matched." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // cemc -- matched to the reference tree by raw channel index (see the note
  // in CaloManip.h / CaloTree.h), not by re-deriving (ieta,iphi) on each side.
  for ( auto ich = 0; ich < k_ncemc; ++ich )
  {
    auto * tower   = cemc_towers -> get_tower_at_channel( ich );
    if ( !tower )
    {
      continue;
    }

    float tree_energy = m_cemc_tower_energy->at(ich);
    int tree_isgood = m_cemc_tower_status->at(ich);
    if ( tree_isgood != 1  || !tower->get_isGood() )
    {
      continue;
    }

    float this_E = tower -> get_energy(); // embedded HIJING+PYTHIA8, calibrated
    float new_E = this_E +  (tree_energy * ( m_scale_factor - 1.0 )); // apply energy scale to tree energy (UE) and add to embedded
    tower -> set_energy( new_E ); // BUGFIX: this write-back was previously missing, making the scaling a no-op
    if ( Verbosity() > 3 )
    {
      std::cout << PHWHERE << " CEMC Tower channel " << ich
        << " new energy: " << new_E
        << " ( original energy: " << this_E << " ) "
        << std::endl;
    }

    if ( debug_event )
    {
      // ieta/iphi/eta here are for the debug display ONLY -- they play no
      // role in the physics matching above.
      const auto key = cemc_towers -> encode_key( ich );
      auto ieta      = cemc_towers -> getTowerEtaBin( key );
      auto iphi      = cemc_towers -> getTowerPhiBin( key );
      float eta = m_towergeom_cemc -> get_tower_geometry(
        RawTowerDefs::encode_towerid( m_cemc_geom_caloid, ieta, iphi ) ) -> get_eta();
      debug_record_tower( "CEMC", ieta, iphi, eta, this_E, tree_energy, new_E );
    }
  }

  for ( auto ich = 0; ich < k_nhcalin; ++ich )
  {
    auto * tower   = hcalin_towers -> get_tower_at_channel( ich );
    if ( !tower )
    {
      continue;
    }

    float tree_energy = m_hcalin_tower_energy->at(ich);
    int tree_isgood = m_hcalin_tower_status->at(ich);
    if ( tree_isgood != 1  || !tower->get_isGood() )
    {
      continue;
    }

    float this_E = tower -> get_energy(); // embedded HIJING+PYTHIA8, calibrated
    float new_E = this_E +  (tree_energy * ( m_scale_factor - 1.0 )); // apply energy scale to tree energy (UE) and add to embedded
    tower -> set_energy( new_E ); // BUGFIX: this write-back was previously missing, making the scaling a no-op
    if ( Verbosity() > 3 )
    {
      std::cout << PHWHERE << " HCALIN Tower channel " << ich
        << " new energy: " << new_E
        << " ( original energy: " << this_E << " ) "
        << std::endl;
    }

    if ( debug_event )
    {
      const auto key = hcalin_towers -> encode_key( ich );
      auto ieta      = hcalin_towers -> getTowerEtaBin( key );
      auto iphi      = hcalin_towers -> getTowerPhiBin( key );
      float eta = m_towergeom_hcalin -> get_tower_geometry(
        RawTowerDefs::encode_towerid( RawTowerDefs::HCALIN, ieta, iphi ) ) -> get_eta();
      debug_record_tower( "HCALIN", ieta, iphi, eta, this_E, tree_energy, new_E );
    }
  }

  for ( auto ich = 0; ich < k_nhcalout; ++ich )
  {
    auto * tower   = hcalout_towers -> get_tower_at_channel( ich );
    if ( !tower )
    {
      continue;
    }

    float tree_energy = m_hcalout_tower_energy->at(ich);
    int tree_isgood = m_hcalout_tower_status->at(ich);
    if ( tree_isgood != 1  || !tower->get_isGood() )
    {
      continue;
    }

    float this_E = tower -> get_energy(); // embedded HIJING+PYTHIA8, calibrated
    float new_E = this_E +  (tree_energy * ( m_scale_factor - 1.0 )); // apply energy scale to tree energy (UE) and add to embedded
    tower -> set_energy( new_E ); // BUGFIX: this write-back was previously missing, making the scaling a no-op
    if ( Verbosity() > 3 )
    {
      std::cout << PHWHERE << " HCALOUT Tower channel " << ich
        << " new energy: " << new_E
        << " ( original energy: " << this_E << " ) "
        << std::endl;
    }

    if ( debug_event )
    {
      const auto key = hcalout_towers -> encode_key( ich );
      auto ieta      = hcalout_towers -> getTowerEtaBin( key );
      auto iphi      = hcalout_towers -> getTowerPhiBin( key );
      float eta = m_towergeom_hcalout -> get_tower_geometry(
        RawTowerDefs::encode_towerid( RawTowerDefs::HCALOUT, ieta, iphi ) ) -> get_eta();
      debug_record_tower( "HCALOUT", ieta, iphi, eta, this_E, tree_energy, new_E );
    }
  }

  if ( debug_event )
  {
    debug_flush_event();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool CaloManip::is_debug_event() const
{
  if ( !m_debug_mode )
  {
    return false;
  }
  return std::find( m_debug_events.begin(), m_debug_events.end(), m_event_id ) != m_debug_events.end();
}

void CaloManip::debug_record_tower( const std::string & calo, int ieta, int iphi, float eta,
                                     float e_before, float e_hijing, float e_after )
{
  const float coshEta = std::cosh( eta );
  const float et_before = e_before / coshEta;
  const float et_after  = e_after  / coshEta;

  DebugRow row;
  row.calo      = calo;
  row.ieta      = ieta;
  row.iphi      = iphi;
  row.eta       = eta;
  row.e_before  = e_before;
  row.e_hijing  = e_hijing;
  row.e_after   = e_after;
  row.et_before = et_before;
  row.et_after  = et_after;
  m_event_debug_rows.push_back( row );

  d_event_id = m_event_id;
  std::snprintf( d_calo, sizeof(d_calo), "%s", calo.c_str() );
  d_ieta = ieta;
  d_iphi = iphi;
  d_eta = eta;
  d_e_before = e_before;
  d_e_hijing = e_hijing;
  d_e_after = e_after;
  d_scale = m_scale_factor;
  d_et_before = et_before;
  d_et_after = et_after;
  m_debug_tower_tree -> Fill();
}

void CaloManip::debug_flush_event()
{
  // index 0=CEMC, 1=HCALIN, 2=HCALOUT, 3=TOTAL
  e_event_id = m_event_id;
  for ( int i = 0; i < 4; ++i )
  {
    e_sumET_hijing[i] = 0.0;
    e_sumET_before[i] = 0.0;
    e_sumET_after[i] = 0.0;
    e_ntowers[i] = 0;
    e_ntowers_modified[i] = 0;
    e_maxAbsDeltaET[i] = 0.0;
    e_maxAbsDeltaET_ieta[i] = -1;
    e_maxAbsDeltaET_iphi[i] = -1;
  }

  auto caloIndex = [] ( const std::string & calo ) -> int
  {
    if ( calo == "CEMC" )    return 0;
    if ( calo == "HCALIN" )  return 1;
    if ( calo == "HCALOUT" ) return 2;
    return -1;
  };

  for ( const auto & row : m_event_debug_rows )
  {
    const int i = caloIndex( row.calo );
    if ( i < 0 ) continue;

    const float et_hijing = row.e_hijing / std::cosh( row.eta );
    const float deltaET = row.et_after - row.et_before;

    e_sumET_hijing[i] += et_hijing;
    e_sumET_before[i] += row.et_before;
    e_sumET_after[i]  += row.et_after;
    e_ntowers[i]++;

    e_sumET_hijing[3] += et_hijing;
    e_sumET_before[3] += row.et_before;
    e_sumET_after[3]  += row.et_after;
    e_ntowers[3]++;

    if ( row.e_after != row.e_before )
    {
      e_ntowers_modified[i]++;
      e_ntowers_modified[3]++;
    }

    if ( std::fabs( deltaET ) > e_maxAbsDeltaET[i] )
    {
      e_maxAbsDeltaET[i] = std::fabs( deltaET );
      e_maxAbsDeltaET_ieta[i] = row.ieta;
      e_maxAbsDeltaET_iphi[i] = row.iphi;
    }
    if ( std::fabs( deltaET ) > e_maxAbsDeltaET[3] )
    {
      e_maxAbsDeltaET[3] = std::fabs( deltaET );
      e_maxAbsDeltaET_ieta[3] = row.ieta;
      e_maxAbsDeltaET_iphi[3] = row.iphi;
    }
  }

  m_debug_event_tree -> Fill();

  // Flush to disk immediately rather than waiting for End(): some inputs
  // make the framework read ahead into the next event's I/O record purely
  // to check for end-of-file, which can throw before End() is ever
  // reached. Writing here guarantees the diagnostics for this event survive
  // regardless of what happens afterward.
  m_debug_tfile -> cd();
  m_debug_tower_tree -> Write( nullptr, TObject::kOverwrite );
  m_debug_event_tree -> Write( nullptr, TObject::kOverwrite );

  // ------------------------------------------------------------------
  // console summary
  // ------------------------------------------------------------------
  const char * names[4] = { "CEMC", "HCALIN", "HCALOUT", "TOTAL" };
  std::cout << "\n================ CaloManip debug summary : Event " << m_event_id
            << " ================\n";
  for ( int i = 0; i < 4; ++i )
  {
    const float ratio = ( e_sumET_before[i] != 0.0 ) ? e_sumET_after[i] / e_sumET_before[i]
                                                       : std::numeric_limits<float>::quiet_NaN();
    std::cout << "-- " << names[i] << " (" << e_ntowers[i] << " towers) --\n"
              << "  HIJING-only sum ET:             " << e_sumET_hijing[i] << " GeV\n"
              << "  HIJING+PYTHIA calibrated ET:     " << e_sumET_before[i] << " GeV\n"
              << "  Scaled/final sum ET:             " << e_sumET_after[i] << " GeV\n"
              << "  Difference after-before:         " << ( e_sumET_after[i] - e_sumET_before[i] ) << " GeV\n"
              << "  Ratio after/before:               " << ratio << "\n"
              << "  Number of towers modified:       " << e_ntowers_modified[i] << "\n"
              << "  Maximum |Delta ET| tower:         " << e_maxAbsDeltaET[i]
              << " GeV  ( ieta=" << e_maxAbsDeltaET_ieta[i] << ", iphi=" << e_maxAbsDeltaET_iphi[i] << " )\n";
  }

  // top-5 towers with the largest correction, across all calorimeters
  std::vector<DebugRow> sorted_rows = m_event_debug_rows;
  std::sort( sorted_rows.begin(), sorted_rows.end(),
             [] ( const DebugRow & a, const DebugRow & b )
             {
               return std::fabs( a.et_after - a.et_before ) > std::fabs( b.et_after - b.et_before );
             } );
  std::cout << "-- Towers with the largest |Delta ET| corrections --\n";
  const int nprint = std::min<int>( 5, static_cast<int>( sorted_rows.size() ) );
  for ( int k = 0; k < nprint; ++k )
  {
    const auto & row = sorted_rows[k];
    std::cout << "  " << row.calo << " ieta=" << row.ieta << " iphi=" << row.iphi
              << " eta=" << row.eta
              << " | E_before=" << row.e_before << " E_hijing=" << row.e_hijing
              << " E_after=" << row.e_after
              << " | ET_before=" << row.et_before << " ET_after=" << row.et_after
              << " Delta_ET=" << ( row.et_after - row.et_before ) << "\n";
  }
  std::cout << "======================================================================\n" << std::endl;

  m_event_debug_rows.clear();
}

int CaloManip::End( PHCompositeNode * /*topNode*/ )
{

  m_tfile -> Close();

  if ( m_debug_mode && m_debug_tfile )
  {
    m_debug_tfile -> cd();
    m_debug_tower_tree -> Write();
    m_debug_event_tree -> Write();
    m_debug_tfile -> Close();
  }

  if ( Verbosity () > 0 )
  {
    std::cout << "CaloManip::End - done" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;

}