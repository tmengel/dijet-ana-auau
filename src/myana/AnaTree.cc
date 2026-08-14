#include "AnaTree.h"

#include <fun4all/PHTFileServer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <ffaobjects/EventHeaderv1.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMapv1.h>

#include <centrality/CentralityInfov2.h>

#include <eventplaneinfo/Eventplaneinfo.h>
#include <eventplaneinfo/EventplaneinfoMap.h>

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/MinimumBiasInfov1.h>


#include <ffarawobjects/Gl1Packetv2.h>

#include <jetbase/Jetv2.h>
#include <jetbase/JetContainerv1.h>

#include <mbd/MbdOutV2.h>

#include <jetbackground/TowerBackgroundv1.h>
#include <jetbackground/TowerRhov1.h>

#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>

AnaTree::AnaTree( const std::string & outputfile )
  : SubsysReco("AnaTree")
  , m_output_filename( outputfile )
{}

int AnaTree::Init( PHCompositeNode * /*topNode*/ )
{
  
    if ( Verbosity () > 0 ) 
    {
        std::cout << "AnaTree::Init - opening file " << m_output_filename << std::endl;
    }

    // init tree
    PHTFileServer::get().open( m_output_filename, "RECREATE" );
    m_tree = new TTree( "T", "T" );
    m_tree -> Branch( "event_id", &m_event_id, "event_id/I" );
    if ( !m_gl1_node.empty() ) 
    {
        m_tree -> Branch( "scaled_triggervec", m_scaled_triggervec, "scaled_triggervec[64]/I" );
        m_tree -> Branch( "live_triggervec", m_live_triggervec, "live_triggervec[64]/I" );
    }
    if ( !m_minbias_node.empty() )
    {
        m_tree -> Branch( "is_minbias", &m_is_minbias, "is_minbias/I" );
    }
    if ( !m_zvrtx_node.empty() )
    {
        m_tree -> Branch( "zvrtx", &m_zvtx, "zvrtx/F" );
    }
    if ( !m_cent_node.empty() )
    {
        m_tree -> Branch( "cent", &m_cent, "cent/I" );
    }
    if ( !m_mbd_node.empty() )
    {
        m_tree -> Branch( "mbd_q_N", &m_mbd_q_N, "mbd_q_N/F" );
        m_tree -> Branch( "mbd_q_S", &m_mbd_q_S, "mbd_q_S/F" );
        m_tree -> Branch( "mbd_t_N", &m_mbd_t_N, "mbd_t_N/F" );
        m_tree -> Branch( "mbd_t_S", &m_mbd_t_S, "mbd_t_S/F" );
    }
    if ( m_save_full_calo )
    {
        m_tree -> Branch( "tower_E", m_tower_E, Form("tower_E[3][%d][%d]/F", k_ieta, k_iphi) );
        m_tree -> Branch( "tower_isgood", m_tower_isgood, Form("tower_isgood[3][%d][%d]/I", k_ieta, k_iphi) );
    }
    if ( m_save_sumeT )
    {
        m_tree -> Branch( "sumeT", m_sumeT, "sumeT[3]/F" );
    }
    if ( m_rho_nodes.size() > 0 )
    {
        m_tree -> Branch( "rho_vals", m_rho_vals, "rho_vals[3]/F" );
        m_tree -> Branch( "rho_sigmas", m_rho_sigmas, "rho_sigmas[3]/F" );
    }
    if ( m_single_rho_nodes.size() > 0 )
    {
       for ( unsigned int i = 0; i < m_single_rho_nodes.size(); ++i )
       {
            m_tree -> Branch( Form("%s_rho_val", m_rho_nicknames[i].c_str()), &m_single_rho_vals[i], Form("%s_rho_val/F", m_rho_nicknames[i].c_str()) );
            m_tree -> Branch( Form("%s_rho_sigma", m_rho_nicknames[i].c_str()), &m_single_rho_sigmas[i], Form("%s_rho_sigma/F", m_rho_nicknames[i].c_str()) );
       }
    }
    if ( !m_towerbkgd_v2_node.empty() )
    {
        m_tree -> Branch( "sub2_v2", &m_sub2_v2, "sub2_v2/F" );
        m_tree -> Branch( "sub2_flowfaliure", &m_sub2_flowfaliure, "sub2_flowfaliure/I" );
        m_tree -> Branch( "sub2_psi2", &m_sub2_psi2, "sub2_psi2/F" );
        m_tree -> Branch( "sub2_towerbkgd_ue", m_sub2_towerbkgd_ue, Form("sub2_towerbkgd_ue[3][%d]/F", k_ieta) );
    }
    if ( !m_sub1jet_node.empty() ) 
    {
        m_tree  -> Branch( "jet_E", &m_sub1_jet_E );
        m_tree  -> Branch( "jet_phi", &m_sub1_jet_phi );
        m_tree  -> Branch( "jet_eta", &m_sub1_jet_eta );
        m_tree  -> Branch( "jet_pT", &m_sub1_jet_pT );
        m_tree  -> Branch( "jet_unsub_pT", &m_sub1_jet_unsub_pT );
        m_tree  -> Branch( "jet_unsub_E", &m_sub1_jet_unsub_E );
        m_tree  -> Branch( "jet_constituent_E", &m_sub1_jet_constituent_E );
        m_tree  -> Branch( "jet_constituent_phi", &m_sub1_jet_constituent_phi );
        m_tree  -> Branch( "jet_constituent_eta", &m_sub1_jet_constituent_eta );
        m_tree  -> Branch( "jet_constituent_pT", &m_sub1_jet_constituent_pT );
        m_tree  -> Branch( "jet_constituent_srcID", &m_sub1_jet_constituent_srcID );
    }
    if ( !m_rhojet_node.empty() ) 
    {
        m_tree  -> Branch( "rho_jet_E", &m_rho_jet_E );
        m_tree  -> Branch( "rho_jet_phi", &m_rho_jet_phi );
        m_tree  -> Branch( "rho_jet_eta", &m_rho_jet_eta );
        m_tree  -> Branch( "rho_jet_pT", &m_rho_jet_pT );
        m_tree  -> Branch( "rho_jet_unsub_pT", &m_rho_jet_unsub_pT );
        m_tree  -> Branch( "rho_jet_unsub_E", &m_rho_jet_unsub_E );
        m_tree  -> Branch( "rho_jet_constituent_E", &m_rho_jet_constituent_E );
        m_tree  -> Branch( "rho_jet_constituent_phi", &m_rho_jet_constituent_phi );
        m_tree  -> Branch( "rho_jet_constituent_eta", &m_rho_jet_constituent_eta );
        m_tree  -> Branch( "rho_jet_constituent_pT", &m_rho_jet_constituent_pT );
        m_tree  -> Branch( "rho_jet_constituent_srcID", &m_rho_jet_constituent_srcID );
    }
    if ( !m_arearhojet_node.empty() ) 
    {
        m_tree  -> Branch( "arearho_jet_E", &m_arearho_jet_E );
        m_tree  -> Branch( "arearho_jet_phi", &m_arearho_jet_phi );
        m_tree  -> Branch( "arearho_jet_eta", &m_arearho_jet_eta );
        m_tree  -> Branch( "arearho_jet_pT", &m_arearho_jet_pT );
        m_tree  -> Branch( "arearho_jet_area", &m_arearho_jet_area );
        m_tree  -> Branch( "arearho_jet_rho", &m_arearho_jet_rho , "arearho_jet_rho/F" );
        m_tree  -> Branch( "arearho_jet_unsub_pT", &m_arearho_jet_unsub_pT );
        m_tree  -> Branch( "arearho_jet_unsub_E", &m_arearho_jet_unsub_E );
        m_tree  -> Branch( "arearho_jet_constituent_E", &m_arearho_jet_constituent_E );
        m_tree  -> Branch( "arearho_jet_constituent_phi", &m_arearho_jet_constituent_phi );
        m_tree  -> Branch( "arearho_jet_constituent_eta", &m_arearho_jet_constituent_eta );
        m_tree  -> Branch( "arearho_jet_constituent_pT", &m_arearho_jet_constituent_pT );
        m_tree  -> Branch( "arearho_jet_constituent_srcID", &m_arearho_jet_constituent_srcID );
    }
    if ( !m_truthjet_node.empty() ) 
    {
        m_tree  -> Branch( "truth_jet_E", &m_truth_jet_E );
        m_tree  -> Branch( "truth_jet_phi", &m_truth_jet_phi );
        m_tree  -> Branch( "truth_jet_eta", &m_truth_jet_eta );
        m_tree  -> Branch( "truth_jet_pT", &m_truth_jet_pT );
    }
    if ( !m_eventhead_node.empty() ) 
    {
        m_tree -> Branch( "b", &m_b, "b/F" );
        m_tree -> Branch( "ep_angle", &m_ep_angle, "ep_angle/F" );
        m_tree -> Branch( "ecc", &m_ecc, "ecc/F" );
        m_tree -> Branch( "psi1", &m_psi1, "psi1/F" );
        m_tree -> Branch( "psi2", &m_psi2, "psi2/F" );
        m_tree -> Branch( "psi3", &m_psi3, "psi3/F" );
        m_tree -> Branch( "ncoll", &m_ncoll, "ncoll/F" );
        m_tree -> Branch( "npart", &m_npart, "npart/F" );
        m_tree -> Branch( "runnumber", &m_runnumber, "runnumber/I" );
        m_tree -> Branch( "evtsequence", &m_evtsequence, "evtsequence/I" );
    }
  
    if ( Verbosity () > 0 )
    {
        std::cout << "AnaTree::Init - done" << std::endl;
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTree::InitRun ( PHCompositeNode */*topNode*/ )
{
    if ( Verbosity () > 0 ) 
    {
        std::cout << "AnaTree::InitRun - initializing run with file " << m_output_filename << std::endl;
    }
    m_event_id = -1;

    return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTree::process_event( PHCompositeNode *topNode )
{
    m_event_id++; 

    if( Verbosity() > 1 ) 
    {
        std::cout << PHWHERE << " Processing event " << m_event_id << std::endl;
    }


    if ( !m_gl1_node.empty() )
    {  
        // get GL1
        memset( m_scaled_triggervec, 0, sizeof(m_scaled_triggervec) );
        memset( m_live_triggervec, 0, sizeof(m_live_triggervec) );
        auto * gl1 = findNode::getClass< Gl1Packetv2 >( topNode, m_gl1_node );
        if( !gl1 ) 
        {
            std::cout << PHWHERE << " No GL1 packet found! Abort." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }
        auto s_triggervec = gl1 -> getScaledVector();
        auto l_triggervec = gl1 -> getLiveVector();
        for ( int i = 0; i < 64; ++i )
        {
            m_scaled_triggervec[i] = (s_triggervec >> i) & 0x1;
            m_live_triggervec[i] = (l_triggervec >> i) & 0x1;
        }
        if ( Verbosity() > 1 ) 
        {
            std::cout << PHWHERE << " - s_triggervec = " << std::hex << s_triggervec << ", l_triggervec = " << l_triggervec << std::dec << std::endl;
        }
    }

    if ( !m_eventhead_node.empty() ) 
    { 
        // get event header info
        m_b = -999;
        m_ep_angle = -999;
        m_ecc = -999;
        m_psi1 = -999;
        m_psi2 = -999;
        m_psi3 = -999;
        m_runnumber = -1;
        m_evtsequence = -1;

        auto * eventhead = findNode::getClass<EventHeader>( topNode, m_eventhead_node );
        if ( !eventhead ) 
        {
            std::cout << PHWHERE << " Input node " << m_eventhead_node << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }
        m_b = eventhead->get_ImpactParameter();
        m_ep_angle = eventhead->get_EventPlaneAngle();
        m_ecc = eventhead->get_eccentricity();
        m_psi1 = eventhead->get_FlowPsiN(1);
        m_psi2 = eventhead->get_FlowPsiN(2);
        m_psi3 = eventhead->get_FlowPsiN(3); 
        // for ( int i = 0; i < 6; ++i )
        // {
        //     m_psi_arr[i] = eventhead->get_FlowPsiN(i+1);
        // }
        m_ncoll = eventhead->get_ncoll();
        m_npart = eventhead->get_npart();
        m_runnumber = eventhead->get_RunNumber();
        m_evtsequence = eventhead->get_EvtSequence();
  
        if ( Verbosity() > 1 ) 
        {
            std::cout << PHWHERE << " - b = " << m_b << ", ep_angle = " << m_ep_angle << ", ecc = " << m_ecc << ", psi2 = " << m_psi2 << ", ncoll = " << m_ncoll << ", npart = " << m_npart << std::endl;
        }

    }
   
    if ( !m_minbias_node.empty() )
    {
        // get minbias
        m_is_minbias = 0;
        auto * minbias = findNode::getClass< MinimumBiasInfov1 >( topNode, m_minbias_node );
        if( !minbias ) 
        {
            std::cout << PHWHERE << " No minbias packet found! Abort." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }

        m_is_minbias = minbias -> isAuAuMinimumBias() ? 1 : 0;
        if ( Verbosity() > 1 ) 
        {
            std::cout << PHWHERE << " - is_minbias = " << m_is_minbias << std::endl;
        }
    }

    if ( !m_cent_node.empty() ) 
    { 
        // get centrality
        m_cent = -1;
        auto * cent_node = findNode::getClass< CentralityInfo >( topNode, m_cent_node );
        if ( !cent_node ) 
        {
            std::cout << PHWHERE << m_cent_node << " node missing, Abort!." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }
        m_cent = static_cast<int>( cent_node -> get_centrality_bin(CentralityInfo::PROP::mbd_NS) );
        if ( Verbosity() > 1 ) 
        {
            std::cout << PHWHERE << "- Centrality = " << m_cent << std::endl;
        }
    }

    if ( !m_zvrtx_node.empty() ) 
    { 
        // get zvtx

        m_zvtx = -999;
        GlobalVertex * vtx { nullptr };
        auto * vertexmap = findNode::getClass<GlobalVertexMap>( topNode, m_zvrtx_node );
        if ( !vertexmap  ) 
        {
            std::cout << PHWHERE << "" << m_zvrtx_node << " node missing, skipping event." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }
        if ( vertexmap->empty() ) 
        {
            std::cout << PHWHERE << "" << m_zvrtx_node << " node has empty vertex map, skipping event." << std::endl;
            return Fun4AllReturnCodes::ABORTEVENT;
        }

        auto vertices = vertexmap -> get_gvtxs_with_type( { GlobalVertex::MBD } );
        if( !vertices.empty() )
        {
            vtx = vertices.at(0);
        }
        else 
        {
            vtx = vertexmap->begin()->second;
        }
        
        if ( vtx )
        {
            m_zvtx = vtx->get_z();
        }
        
        if ( std::isnan(m_zvtx) || std::abs(m_zvtx) > 1e3 )
        {
            static bool z_warning_once = true;
            if ( z_warning_once )
            {
                z_warning_once = false;
                std::cout << PHWHERE << " vertex z is " << m_zvtx << ", skipping event (further warnings will be suppressed)." << std::endl;
            }
            
            return Fun4AllReturnCodes::ABORTEVENT;
            
        }
 
        if ( Verbosity() > 1 ) 
        {
            std::cout << PHWHERE << " - zvtx = " << m_zvtx << std::endl;
        }

    }

    if ( !m_mbd_node.empty() ) 
    { 
        // get mbd info
        m_mbd_q_N = -999;
        m_mbd_q_S = -999;
        m_mbd_t_N = -999;
        m_mbd_t_S = -999;
        auto * mbd_node = findNode::getClass< MbdOutV2 >( topNode, m_mbd_node );
        if ( !mbd_node ) 
        {
            std::cout << PHWHERE << m_mbd_node << " node missing, skipping event." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }
        m_mbd_q_N = mbd_node -> get_q(1);
        m_mbd_q_S = mbd_node -> get_q(0);
        m_mbd_t_N = mbd_node -> get_time(1);
        m_mbd_t_S = mbd_node -> get_time(0);
        if ( Verbosity() > 1 ) 
        {
            std::cout << PHWHERE << " - mbd_q_N = " << m_mbd_q_N << ", mbd_q_S = " << m_mbd_q_S << ", mbd_t_N = " << m_mbd_t_N << ", mbd_t_S = " << m_mbd_t_S << std::endl;
        }

    }

    if ( m_single_rho_nodes.size() > 0 ) 
    {
        for ( unsigned int i = 0; i < m_single_rho_nodes.size(); ++i )
        {
            m_single_rho_vals[i] = -999;
            m_single_rho_sigmas[i] = -999;
            auto * rho_node = findNode::getClass< TowerRhov1 >( topNode, m_single_rho_nodes[i] );
            if ( !rho_node ) 
            {
                std::cout << PHWHERE << " Input node " << m_single_rho_nodes[i] << " Node missing, doing nothing." << std::endl;
                return Fun4AllReturnCodes::ABORTRUN; 
            }
            m_single_rho_vals[i] = rho_node -> get_rho();
            m_single_rho_sigmas[i] = rho_node -> get_sigma();
            if ( Verbosity() > 1 ) 
            {
                std::cout << PHWHERE << " - " << m_single_rho_nodes[i] << "_rho_val = " << m_single_rho_vals[i] << ", " << m_single_rho_nodes[i] << "_rho_sigma = " << m_single_rho_sigmas[i] << std::endl;
            }
        }
    }

    if ( !m_arearhojet_node.empty() ) 
    { 
        m_arearho_jet_E.clear();
        m_arearho_jet_phi.clear();
        m_arearho_jet_eta.clear();
        m_arearho_jet_pT.clear();
        m_arearho_jet_area.clear();
        m_arearho_jet_rho = 0;
        m_arearho_jet_unsub_pT.clear();
        m_arearho_jet_unsub_E.clear();
        m_arearho_jet_constituent_E.clear();
        m_arearho_jet_constituent_phi.clear();
        m_arearho_jet_constituent_eta.clear();
        m_arearho_jet_constituent_pT.clear();
        m_arearho_jet_constituent_srcID.clear();

        auto * jets = findNode::getClass<JetContainer>( topNode, m_arearhojet_node );
        if ( !jets )
        {
            std::cout << PHWHERE << " Input node " << m_arearhojet_node << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN; 
        }

        m_arearho_jet_rho = jets -> get_rho_median();
        auto area_index = jets -> property_index(Jet::PROPERTY::prop_area);
        if ( Verbosity() > 1 ) 
        {
            std::cout << PHWHERE << " - arearho_jet_rho = " << m_arearho_jet_rho << std::endl;
        }

        auto * rho_area_node = findNode::getClass<TowerRhov1>( topNode, "TowerRho_AREA" );
        float rho_from_node = 0;
        if ( rho_area_node )
        {
            rho_from_node = rho_area_node -> get_rho();
            std::cout << PHWHERE << " - arearho_jet_rho = " << m_arearho_jet_rho << ", rho_from_node = " << rho_from_node << std::endl;
        }

        for ( const auto & jet : * jets )
        {


            std::vector<float> constituent_E {};
            std::vector<float> constituent_phi {};
            std::vector<float> constituent_eta {};
            std::vector<float> constituent_pT {};
            std::vector<int> constituent_srcID {};
            for ( const auto & comp : jet -> get_comp_vec() )
            {
                double tower_r = 0.0;
                m_caloid = RawTowerDefs::CalorimeterId::NONE;
                m_towerinfos = nullptr;
                m_towergeom = nullptr;
                if( comp.first == Jet::SRC::HCALIN_TOWERINFO )
                {
                    // const float CALO_RADIUS[3] = {93.5, 127.503, 225.87};
                    m_towerinfos = LoadTowerInfoContainer( topNode, "TOWERINFO_CALIB_HCALIN" );
                    m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALIN" );
                    tower_r = 127.503;
                    m_caloid = RawTowerDefs::CalorimeterId::HCALIN;
                }
                else if ( comp.first == Jet::SRC::HCALOUT_TOWERINFO )
                {
                    m_towerinfos = LoadTowerInfoContainer( topNode, "TOWERINFO_CALIB_HCALOUT" );
                    m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALOUT" );
                    tower_r = 225.87;
                    m_caloid = RawTowerDefs::CalorimeterId::HCALOUT;
                }
                else if ( comp.first == Jet::SRC::CEMC_TOWERINFO_RETOWER )
                {
                    m_towerinfos = LoadTowerInfoContainer( topNode, "TOWERINFO_CALIB_CEMC_RETOWER" );
                    m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALIN" );
                    tower_r = 93.5;
                    m_caloid = RawTowerDefs::CalorimeterId::HCALIN; // use hcalin geometry for cemc towers since we just want eta/phi and the r is only used for calculating unsub pT which will be corrected by UE subtraction
                }
                else 
                {
                    if ( Verbosity() > 3 ) 
                    {
                        std::cout << PHWHERE << " Warning: jet constituent with unknown source " << comp.first << ", skipping." << std::endl;
                    }
                    continue;
                }
                
                auto * tower = m_towerinfos->get_tower_at_channel( comp.second );
                if ( !tower || !tower->get_isGood() ) 
                {
                    if ( Verbosity() > 3 ) 
                    {
                        std::cout << PHWHERE << " Warning: constituent tower with caloid " << comp.first << " and channel " << comp.second << " not found in towerinfo container, skipping." << std::endl;
                    }
                    continue;
                }

                const auto tkey = m_towerinfos->encode_key( comp.second );
                auto ieta = m_towerinfos->getTowerEtaBin(tkey);
                auto iphi = m_towerinfos->getTowerPhiBin(tkey);
                const auto key = RawTowerDefs::encode_towerid( m_caloid, ieta, iphi );
                
                auto * geom = m_towergeom->get_tower_geometry( key );
                if ( !geom )
                {
                    if ( Verbosity() > 3 ) 
                    {
                        std::cout << PHWHERE << " Warning: geometry for tower with caloid " << comp.first << " and channel " << comp.second << " not found, skipping." << std::endl;
                    }
                    continue;
                }

                double tower_z0   = sinh( geom -> get_eta() ) * tower_r;
                double comp_z     = tower_z0 - m_zvtx;
                double comp_eta   =  asinh( comp_z / tower_r );
                double comp_phi   = geom -> get_phi();
                double comp_E     = tower->get_energy();
                double ue         = 0;
                double un_E       = comp_E + ue;
                float un_pT = un_E / cosh( comp_eta );
                constituent_E.push_back(comp_E);
                constituent_phi.push_back(comp_phi);
                constituent_eta.push_back(comp_eta);
                constituent_pT.push_back(un_pT);
                constituent_srcID.push_back(static_cast<int>(comp.first));  

            } // end loop over constituents
            
            auto jet_area = jet -> get_property(area_index);
            m_arearho_jet_E.push_back(jet->get_e());
            m_arearho_jet_phi.push_back(jet->get_phi());
            m_arearho_jet_eta.push_back(jet->get_eta());
            m_arearho_jet_pT.push_back(jet->get_pt() - m_arearho_jet_rho * jet_area);
            m_arearho_jet_area.push_back(jet_area);
            m_arearho_jet_unsub_pT.push_back(jet->get_pt());
            m_arearho_jet_unsub_E.push_back(jet->get_e());
            m_arearho_jet_constituent_E.push_back(constituent_E);
            m_arearho_jet_constituent_phi.push_back(constituent_phi);
            m_arearho_jet_constituent_eta.push_back(constituent_eta);
            m_arearho_jet_constituent_pT.push_back(constituent_pT);
            m_arearho_jet_constituent_srcID.push_back(constituent_srcID);

        } // end loop over jets

        
    }
    
    if ( !m_sub1jet_node.empty() ) 
    { 
        
        m_sub1_jet_E.clear();
        m_sub1_jet_phi.clear();
        m_sub1_jet_eta.clear();
        m_sub1_jet_pT.clear();
        m_sub1_jet_unsub_pT.clear();
        m_sub1_jet_unsub_E.clear();
        m_sub1_jet_constituent_E.clear();
        m_sub1_jet_constituent_phi.clear();
        m_sub1_jet_constituent_eta.clear();
        m_sub1_jet_constituent_pT.clear();
        m_sub1_jet_constituent_srcID.clear();

        memset(m_sub2_towerbkgd_ue, 0, sizeof(m_sub2_towerbkgd_ue));
        m_sub2_v2 = 0.0;
        m_sub2_flowfaliure = 0;
        m_sub2_psi2 = 0.0;

        // get sub1 jet info
        auto * jets = findNode::getClass<JetContainer>( topNode, m_sub1jet_node );
        if ( !jets )
        {
            std::cout << PHWHERE << " Input node " << m_sub1jet_node << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN; 
        }
        
    

        auto * tower_background_sub2 = findNode::getClass<TowerBackgroundv1>(topNode, m_towerbkgd_v2_node);
        if ( !tower_background_sub2  )
        {
            std::cout << PHWHERE << " TowerBackgroundv1 node is missing, skipping." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN; // fatal error
        }
        else 
        {
            // read info
            for ( size_t ilay = 0 ; ilay < 3; ++ilay ) 
            {
                auto this_ue_sub2 = tower_background_sub2->get_UE(ilay);
                for ( size_t ieta = 0; ieta < k_ieta; ++ieta )
                {
                    m_sub2_towerbkgd_ue[ilay][ieta] = this_ue_sub2[ieta];
                }
            }
            m_sub2_v2 = tower_background_sub2->get_v2();
            m_sub2_flowfaliure = tower_background_sub2->get_flow_failure_flag();
            m_sub2_psi2 = tower_background_sub2->get_Psi2();
        }
        

        for ( const auto & jet : * jets )
        {
            
            float unsub_pz = 0;
            float unsub_px = 0; 
            float unsub_py = 0;
            float unsub_E  = 0;
            std::vector<float> constituent_E {};
            std::vector<float> constituent_phi {};
            std::vector<float> constituent_eta {};
            std::vector<float> constituent_pT {};
            std::vector<int> constituent_srcID {};
            for ( const auto & comp : jet -> get_comp_vec() )
            {
                double tower_r = 0.0;
                int layer_idx = -1;
                m_caloid = RawTowerDefs::CalorimeterId::NONE;
                m_towerinfos = nullptr;
                m_towergeom = nullptr;
                if( comp.first == Jet::SRC::HCALIN_TOWERINFO_SUB1 )
                {
                    // const float CALO_RADIUS[3] = {93.5, 127.503, 225.87};
                    m_towerinfos = LoadTowerInfoContainer( topNode, "TOWERINFO_CALIB_HCALIN_SUB1" );
                    m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALIN" );
                    tower_r = 127.503;
                    m_caloid = RawTowerDefs::CalorimeterId::HCALIN;
                    layer_idx = 1;
                }
                else if ( comp.first == Jet::SRC::HCALOUT_TOWERINFO_SUB1 )
                {
                    m_towerinfos = LoadTowerInfoContainer( topNode, "TOWERINFO_CALIB_HCALOUT_SUB1" );
                    m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALOUT" );
                    tower_r = 225.87;
                    m_caloid = RawTowerDefs::CalorimeterId::HCALOUT;
                    layer_idx = 2;
                }
                else if ( comp.first == Jet::SRC::CEMC_TOWERINFO_SUB1 )
                {
                    m_towerinfos = LoadTowerInfoContainer( topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1" );
                    m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALIN" );
                    tower_r = 93.5;
                    layer_idx = 0;
                    m_caloid = RawTowerDefs::CalorimeterId::HCALIN; // use hcalin geometry for cemc towers since we just want eta/phi and the r is only used for calculating unsub pT which will be corrected by UE subtraction
                }
                else 
                {
                    if ( Verbosity() > 3 ) 
                    {
                        std::cout << PHWHERE << " Warning: jet constituent with unknown source " << comp.first << ", skipping." << std::endl;
                    }
                    continue;
                }
                
                auto * tower = m_towerinfos->get_tower_at_channel( comp.second );
                if ( !tower || !tower->get_isGood() ) 
                {
                    if ( Verbosity() > 3 ) 
                    {
                        std::cout << PHWHERE << " Warning: constituent tower with caloid " << comp.first << " and channel " << comp.second << " not found in towerinfo container, skipping." << std::endl;
                    }
                    continue;
                }

                const auto tkey = m_towerinfos->encode_key( comp.second );
                auto ieta = m_towerinfos->getTowerEtaBin(tkey);
                auto iphi = m_towerinfos->getTowerPhiBin(tkey);
                const auto key = RawTowerDefs::encode_towerid( m_caloid, ieta, iphi );
                
                auto * geom = m_towergeom->get_tower_geometry( key );
                if ( !geom )
                {
                    if ( Verbosity() > 3 ) 
                    {
                        std::cout << PHWHERE << " Warning: geometry for tower with caloid " << comp.first << " and channel " << comp.second << " not found, skipping." << std::endl;
                    }
                    continue;
                }

                double tower_z0   = sinh( geom -> get_eta() ) * tower_r;
                double comp_z     = tower_z0 - m_zvtx;
                double comp_eta   =  asinh( comp_z / tower_r );
                double comp_phi   = geom -> get_phi();
                double comp_E     = tower->get_energy();
                double ue         = m_sub2_towerbkgd_ue[layer_idx][ieta];
                double un_E       = comp_E + ue;
                float un_pT = un_E / cosh( comp_eta );
                // float pT = comp_E / cosh( comp_eta );
                float un_px = un_pT * cos( comp_phi );
                float un_py = un_pT * sin( comp_phi );
                float un_pz = un_pT * sinh( comp_eta );
                unsub_px += un_px;
                unsub_py += un_py;
                unsub_pz += un_pz;
                unsub_E  += un_E;

                constituent_E.push_back(comp_E);
                constituent_phi.push_back(comp_phi);
                constituent_eta.push_back(comp_eta);
                constituent_pT.push_back(un_pT);
                constituent_srcID.push_back(static_cast<int>(comp.first));  

            } // end loop over constituents
            
            auto * unsub_jet = new Jetv2();
            unsub_jet->set_px(unsub_px);
            unsub_jet->set_py(unsub_py);
            unsub_jet->set_pz(unsub_pz);
            unsub_jet->set_e(unsub_E);
            
            m_sub1_jet_E.push_back(jet->get_e());
            m_sub1_jet_eta.push_back(jet->get_eta());
            m_sub1_jet_phi.push_back(jet->get_phi());
            m_sub1_jet_pT.push_back(jet->get_pt());
            m_sub1_jet_unsub_pT.push_back(unsub_jet->get_pt());
            m_sub1_jet_unsub_E.push_back(unsub_jet->get_e());
            m_sub1_jet_constituent_E.push_back(constituent_E);
            m_sub1_jet_constituent_phi.push_back(constituent_phi);
            m_sub1_jet_constituent_eta.push_back(constituent_eta);
            m_sub1_jet_constituent_pT.push_back(constituent_pT);
            m_sub1_jet_constituent_srcID.push_back(constituent_srcID);


        } // end loop over jets
    }   

    if ( !m_rhojet_node.empty() )
    { 
        
        m_rho_jet_E.clear();
        m_rho_jet_phi.clear();
        m_rho_jet_eta.clear();
        m_rho_jet_pT.clear();
        m_rho_jet_unsub_pT.clear();
        m_rho_jet_unsub_E.clear();
        m_rho_jet_constituent_E.clear();
        m_rho_jet_constituent_phi.clear();
        m_rho_jet_constituent_eta.clear();
        m_rho_jet_constituent_pT.clear();
        m_rho_jet_constituent_srcID.clear();

        memset(m_rho_vals, 0, sizeof(m_rho_vals));
        memset(m_rho_sigmas, 0, sizeof(m_rho_sigmas));

        // get rho jet info
        auto * jets = findNode::getClass<JetContainer>( topNode, m_rhojet_node );
        if ( !jets )
        {
            std::cout << PHWHERE << " Input node " << m_rhojet_node << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN; 
        }
        
    
        for ( size_t ilay = 0 ; ilay < 3; ++ilay )
        {
            auto * rho_node = findNode::getClass<TowerRhov1>( topNode, m_rho_nodes[ilay] );
            if ( !rho_node  )
            {
                std::cout << PHWHERE << " TowerRhov1 node " << m_rho_nodes[ilay] << " is missing, skipping." << std::endl;
                // return Fun4AllReturnCodes::ABORTRUN; // fatal error
                continue;
            }
            else 
            {
                m_rho_vals[ilay] = rho_node->get_rho();
                m_rho_sigmas[ilay] = rho_node->get_sigma();
            }
        }
        

        for ( const auto & jet : * jets )
        {
            
            float unsub_pz = 0;
            float unsub_px = 0; 
            float unsub_py = 0;
            float unsub_E  = 0;
            std::vector<float> constituent_E {};
            std::vector<float> constituent_phi {};
            std::vector<float> constituent_eta {};
            std::vector<float> constituent_pT {};
            std::vector<int> constituent_srcID {};
            for ( const auto & comp : jet -> get_comp_vec() )
            {
                double tower_r = 0.0;
                int layer_idx = -1;
                m_caloid = RawTowerDefs::CalorimeterId::NONE;
                m_towerinfos = nullptr;
                m_towergeom = nullptr;
                if( comp.first == Jet::SRC::HCALIN_TOWERINFO )
                {
                    // const float CALO_RADIUS[3] = {93.5, 127.503, 225.87};
                    m_towerinfos = LoadTowerInfoContainer( topNode, "TOWERINFO_CALIB_HCALIN" );
                    m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALIN" );
                    tower_r = 127.503;
                    m_caloid = RawTowerDefs::CalorimeterId::HCALIN;
                    layer_idx = 1;
                }
                else if ( comp.first == Jet::SRC::HCALOUT_TOWERINFO )
                {
                    m_towerinfos = LoadTowerInfoContainer( topNode, "TOWERINFO_CALIB_HCALOUT" );
                    m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALOUT" );
                    tower_r = 225.87;
                    m_caloid = RawTowerDefs::CalorimeterId::HCALOUT;
                    layer_idx = 2;
                }
                else if ( comp.first == Jet::SRC::CEMC_TOWERINFO_RETOWER )
                {
                    m_towerinfos = LoadTowerInfoContainer( topNode, "TOWERINFO_CALIB_CEMC_RETOWER" );
                    m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALIN" );
                    tower_r = 93.5;
                    layer_idx = 0;
                    m_caloid = RawTowerDefs::CalorimeterId::HCALIN; // use hcalin geometry for cemc towers since we just want eta/phi and the r is only used for calculating unsub pT which will be corrected by UE subtraction
                }
                else 
                {
                    if ( Verbosity() > 3 ) 
                    {
                        std::cout << PHWHERE << " Warning: jet constituent with unknown source " << comp.first << ", skipping." << std::endl;
                    }
                    continue;
                }
                
                auto * tower = m_towerinfos->get_tower_at_channel( comp.second );
                if ( !tower || !tower->get_isGood() ) 
                {
                    if ( Verbosity() > 3 ) 
                    {
                        std::cout << PHWHERE << " Warning: constituent tower with caloid " << comp.first << " and channel " << comp.second << " not found in towerinfo container, skipping." << std::endl;
                    }
                    continue;
                }

                const auto tkey = m_towerinfos->encode_key( comp.second );
                auto ieta = m_towerinfos->getTowerEtaBin(tkey);
                auto iphi = m_towerinfos->getTowerPhiBin(tkey);
                const auto key = RawTowerDefs::encode_towerid( m_caloid, ieta, iphi );
                
                auto * geom = m_towergeom->get_tower_geometry( key );
                if ( !geom )
                {
                    if ( Verbosity() > 3 ) 
                    {
                        std::cout << PHWHERE << " Warning: geometry for tower with caloid " << comp.first << " and channel " << comp.second << " not found, skipping." << std::endl;
                    }
                    continue;
                }

                double tower_z0   = sinh( geom -> get_eta() ) * tower_r;
                double comp_z     = tower_z0 - m_zvtx;
                double comp_eta   =  asinh( comp_z / tower_r );
                double comp_phi   = geom -> get_phi();
                double comp_E     = tower->get_energy();
                double ue         = m_rho_vals[layer_idx];
                double un_E       = comp_E + ue;
                float un_pT = un_E / cosh( comp_eta );
                float un_px = un_pT * cos( comp_phi );
                float un_py = un_pT * sin( comp_phi );
                float un_pz = un_pT * sinh( comp_eta );
                unsub_px += un_px;
                unsub_py += un_py;
                unsub_pz += un_pz;
                unsub_E  += un_E;

                constituent_E.push_back(comp_E);
                constituent_phi.push_back(comp_phi);
                constituent_eta.push_back(comp_eta);
                constituent_pT.push_back(un_pT);
                constituent_srcID.push_back(static_cast<int>(comp.first));  

            } // end loop over constituents
            
            auto * unsub_jet = new Jetv2();
            unsub_jet->set_px(unsub_px);
            unsub_jet->set_py(unsub_py);
            unsub_jet->set_pz(unsub_pz);
            unsub_jet->set_e(unsub_E);
            
            m_rho_jet_E.push_back(jet->get_e());
            m_rho_jet_eta.push_back(jet->get_eta());
            m_rho_jet_phi.push_back(jet->get_phi());
            m_rho_jet_pT.push_back(jet->get_pt());
            m_rho_jet_unsub_pT.push_back(unsub_jet->get_pt());
            m_rho_jet_unsub_E.push_back(unsub_jet->get_e());
            m_rho_jet_constituent_E.push_back(constituent_E);
            m_rho_jet_constituent_phi.push_back(constituent_phi);
            m_rho_jet_constituent_eta.push_back(constituent_eta);
            m_rho_jet_constituent_pT.push_back(constituent_pT);
            m_rho_jet_constituent_srcID.push_back(constituent_srcID);


        } // end loop over jets
    }

    if ( !m_truthjet_node.empty() ) 
    { 
        // get truth jet info
        m_truth_jet_E.clear();
        m_truth_jet_phi.clear();
        m_truth_jet_eta.clear();
        m_truth_jet_pT.clear();
        auto * truth_jets = findNode::getClass<JetContainer>( topNode, m_truthjet_node );
        if ( !truth_jets )        
        {
            std::cout << PHWHERE << " Input node " << m_truthjet_node << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN; 
        }
        for ( const auto & jet : * truth_jets )
        {
            m_truth_jet_E.push_back(jet->get_e());
            m_truth_jet_eta.push_back(jet->get_eta());
            m_truth_jet_phi.push_back(jet->get_phi());
            m_truth_jet_pT.push_back(jet->get_pt());
        } // end loop over truth jets
    }

    if ( m_save_full_calo )
    {
        memset(m_tower_E, 0, sizeof(m_tower_E));
        memset(m_tower_isgood, 0, sizeof(m_tower_isgood));
        for ( size_t ilay = 0; ilay < 3; ++ilay )   
        {
            m_towerinfos = LoadTowerInfoContainer( topNode, m_fullcalo_nodes[ilay] );
            for ( unsigned int ich = 0; ich < m_towerinfos->size(); ich++ ) 
            {
                auto tower = m_towerinfos->get_tower_at_channel(ich);

                if ( !tower || ! tower->get_isGood() || std::isnan(tower->get_energy() ) )
                {
                    
                    continue; // skip bad towers
                }
                unsigned int key = m_towerinfos -> encode_key(ich);
                int ieta = m_towerinfos -> getTowerEtaBin(key);
                int iphi = m_towerinfos -> getTowerPhiBin(key);
                m_tower_E[ilay][ieta][iphi] = tower->get_energy();
                m_tower_isgood[ilay][ieta][iphi] = 1;

            }
        }
    }
    
    if ( m_save_sumeT )
    {
        memset(m_sumeT, 0, sizeof(m_sumeT));
        auto sumet_vec = SumCaloE( topNode, m_sumeT_nodes );
        for ( size_t ilay = 0; ilay < 3; ++ilay )
        {
            m_sumeT[ilay] = sumet_vec[ilay];
        }
    }
    

    // fill tree
    m_tree->Fill();
    
    return Fun4AllReturnCodes::EVENT_OK;

}

int AnaTree::End( PHCompositeNode * /*topNode*/ )
{
  
  if( Verbosity() > 0 ) 
  {
    std::cout << "AnaTree::EndRun - End run " << std::endl;
    std::cout << "AnaTree::EndRun - Writing to " << m_output_filename << std::endl;
  }

  PHTFileServer::get().cd(m_output_filename); 

  m_tree->Write();
  

  if ( Verbosity() > 0 ) 
  {
    std::cout << "AnaTree::EndRun - Writing run tree" << std::endl;
  }
  PHTFileServer::get().close();
 
  if ( Verbosity () > 0 ) 
  {
    std::cout << "AnaTree::EndRun - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

TowerInfoContainer * AnaTree::LoadTowerInfoContainer(PHCompositeNode *topNode, const std::string &tower_node_name)
{
  auto * tic = findNode::getClass<TowerInfoContainer>(topNode, tower_node_name.c_str());
  if (!tic)
  {
    std::cout << "AnaTree::LoadTowerInfoContainer - cannot find " << tower_node_name << ", exiting" << std::endl;
    exit(1);
  }
  return tic;
}

RawTowerGeomContainer * AnaTree::LoadTowerGeomContainer(PHCompositeNode *topNode, const std::string &tower_geom_node_name)
{
  auto * tg = findNode::getClass<RawTowerGeomContainer>(topNode, tower_geom_node_name.c_str());
  if (!tg)
  {
    std::cout << "AnaTree::LoadTowerGeomContainer - cannot find " << tower_geom_node_name << ", exiting" << std::endl;
    exit(1);
  }
  return tg;
}

std::vector<float> AnaTree::SumCaloE( PHCompositeNode *topNode,  const std::vector< std::string > &towerinfo_nodes )
{
    std::vector<float> sumeT(3, 0.0); 
    for ( size_t ilay = 0; ilay < 3; ++ilay )   
    {
        double tower_r =0.0;
        m_caloid = RawTowerDefs::CalorimeterId::NONE;
        m_towerinfos = LoadTowerInfoContainer( topNode, towerinfo_nodes[ilay] );
        m_towergeom = nullptr;
        if ( ilay == 0 ) 
        {
            m_caloid = RawTowerDefs::CalorimeterId::HCALIN;
            // tower_r = m_calo_r[RawTowerDefs::CalorimeterId::CEMC-1];
            tower_r = 93.5; // use hard coded radius for cemc since we just want eta/phi and the r is only used for calculating unsub pT which will be corrected by UE subtraction
            m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALIN" );
        }
        else if ( ilay == 1 ) 
        {
            m_caloid = RawTowerDefs::CalorimeterId::HCALIN;
            // tower_r = m_calo_r[RawTowerDefs::CalorimeterId::HCALIN-1];
            tower_r = 127.503;
            m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALIN" );
        }
        else if ( ilay == 2 ) 
        {
            m_caloid = RawTowerDefs::CalorimeterId::HCALOUT;
            // tower_r = m_calo_r[RawTowerDefs::CalorimeterId::HCALOUT-1];
            tower_r = 225.87;
            m_towergeom = LoadTowerGeomContainer( topNode, "TOWERGEOM_HCALOUT" );
        }

        for ( unsigned int ich = 0; ich < m_towerinfos->size(); ich++ ) 
        {
            auto tower = m_towerinfos->get_tower_at_channel(ich);
            if ( !tower || ! tower->get_isGood() || std::isnan(tower->get_energy() ) )
            {
                continue; // skip bad towers
            }

            unsigned int key = m_towerinfos -> encode_key(ich);
            int ieta = m_towerinfos -> getTowerEtaBin(key);
            int iphi = m_towerinfos -> getTowerPhiBin(key);
            const RawTowerDefs::keytype geokey = RawTowerDefs::encode_towerid(m_caloid, ieta, iphi);

            auto * tower_geom = m_towergeom -> get_tower_geometry( geokey );
            if ( !tower_geom )
            {
                std::cout << PHWHERE << " Warning: cannot find geometry for tower with calo id " << static_cast<int>(m_caloid) << " and channel " << ich << ", skipping." << std::endl;
                exit(1);
                continue;
            }
            
            double tower_z0   = sinh( tower_geom -> get_eta() ) * tower_r;
            double tower_z    = tower_z0 - m_zvtx;
            double tower_eta  =  asinh( tower_z / tower_r );
            sumeT[ilay] += tower->get_energy() / cosh( tower_eta );
        }
    }

    return sumeT;
}







   