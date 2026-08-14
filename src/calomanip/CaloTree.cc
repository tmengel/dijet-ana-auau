#include "CaloTree.h"

#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <phool/getClass.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <ffaobjects/EventHeaderv1.h>


#include <TTree.h>

int CaloTree::Init( PHCompositeNode * /*topNode*/ )
{
  
  PHTFileServer::get().open( m_output_filename, "RECREATE" );

  if ( Verbosity () > 0 ) 
  {
    std::cout << "CaloTree::Init - opening file " << m_output_filename << std::endl;
  } 
  
  m_event_id = -1;
  m_tree = new TTree( "T", "T" );
  m_tree -> Branch( "event_id", &m_event_id, "event_id/I" );
  m_tree -> Branch( "runnumber", &m_runnumber, "runnumber/I" );
  m_tree -> Branch( "evtsequence", &m_evtsequence, "evtsequence/I" );
  m_tree -> Branch( "b", &m_b, "b/F" );
  m_tree -> Branch( "ep_angle", &m_ep_angle, "ep_angle/F" );
  m_tree -> Branch( "ecc", &m_ecc, "ecc/F" );
  m_tree -> Branch( "ncoll", &m_ncoll, "ncoll/F" );
  m_tree -> Branch( "npart", &m_npart, "npart/F" );
  m_tree -> Branch( "psi1", &m_psi1, "psi1/F" );
  m_tree -> Branch( "psi2", &m_psi2, "psi2/F" );
  m_tree -> Branch( "psi3", &m_psi3, "psi3/F" );
  m_tree -> Branch( "ncemc", &k_ncemc, "ncemc/I" );
  m_tree -> Branch( "nhcalin", &k_nhcalin, "nhcalin/I" );
  m_tree -> Branch( "nhcalout", &k_nhcalout, "nhcalout/I" );
  m_tree -> Branch( "cemc_tower_energy", &m_cemc_tower_energy );
  m_tree -> Branch( "cemc_tower_status", &m_cemc_tower_status );
  m_tree -> Branch( "hcalin_tower_energy", &m_hcalin_tower_energy );
  m_tree -> Branch( "hcalin_tower_status", &m_hcalin_tower_status );
  m_tree -> Branch( "hcalout_tower_energy", &m_hcalout_tower_energy );
  m_tree -> Branch( "hcalout_tower_status", &m_hcalout_tower_status );

  if ( Verbosity () > 0 )
  {
    std::cout << "CaloTree::Init - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


int CaloTree::process_event( PHCompositeNode *topNode )
{

  m_event_id++; 

  if( Verbosity() > 1 ) 
  {
    std::cout << "CaloTree::process_event - Process event " << m_event_id << std::endl;
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
    std::cout << PHWHERE << " - runnumber = " << m_runnumber << ", evtsequence = " << m_evtsequence << ", b = " << m_b << ", ep_angle = " << m_ep_angle << ", ecc = " << m_ecc << ", psi2 = " << m_psi2 << ", ncoll = " << m_ncoll << ", npart = " << m_npart << std::endl;
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

  // Store per-channel, in the container's own native channel order -- this
  // is later read back in CaloManip and matched channel-for-channel against
  // the live embedded container of the SAME node type, so no (ieta,iphi)
  // grid-size assumption (retowered 24x64 vs native 96x256) ever needs to
  // be made or re-derived on either side of the two passes.
  k_ncemc = cemc_towers -> size();
  k_nhcalin = hcalin_towers -> size();
  k_nhcalout = hcalout_towers -> size();
  m_cemc_tower_energy.assign( k_ncemc, 0.0 );
  m_cemc_tower_status.assign( k_ncemc, 0 );
  m_hcalin_tower_energy.assign( k_nhcalin, 0.0 );
  m_hcalin_tower_status.assign( k_nhcalin, 0 );
  m_hcalout_tower_energy.assign( k_nhcalout, 0.0 );
  m_hcalout_tower_status.assign( k_nhcalout, 0 );

  for ( auto ich = 0; ich < k_ncemc; ++ich )
  {
    auto * tower = cemc_towers -> get_tower_at_channel( ich );
    if ( !tower )
    {
      continue;
    }
    m_cemc_tower_energy[ich] = tower -> get_energy();
    m_cemc_tower_status[ich] = tower -> get_isGood();
  }
  for ( auto ich = 0; ich < k_nhcalin; ++ich )
  {
    auto * tower = hcalin_towers -> get_tower_at_channel( ich );
    if ( !tower )
    {
      continue;
    }
    m_hcalin_tower_energy[ich] = tower -> get_energy();
    m_hcalin_tower_status[ich] = tower -> get_isGood();
  }
  for ( auto ich = 0; ich < k_nhcalout; ++ich )
  {
    auto * tower = hcalout_towers -> get_tower_at_channel( ich );
    if ( !tower )
    {
      continue;
    }
    m_hcalout_tower_energy[ich] = tower -> get_energy();
    m_hcalout_tower_status[ich] = tower -> get_isGood();
  }

  m_tree -> Fill();

  return Fun4AllReturnCodes::EVENT_OK;

}

int CaloTree::End( PHCompositeNode * /*topNode*/ )
{
  
  PHTFileServer::get().cd(m_output_filename); 
  m_tree -> Write();
  PHTFileServer::get().close();

  if ( Verbosity () > 0 ) 
  {
    std::cout << "CaloTree::EndRun - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;

}







   