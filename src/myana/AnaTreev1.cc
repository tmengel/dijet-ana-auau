#include "AnaTreev1.h"

#include <fun4all/PHTFileServer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <ffaobjects/EventHeaderv1.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
#include <HepMC/PdfInfo.h>

#pragma GCC diagnostic pop

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

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
#include <TH1D.h>
#include <TMath.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>

AnaTreev1::AnaTreev1( const std::string & outputfile )
  : SubsysReco("AnaTreev1")
  , m_output_filename( outputfile )
{}

int AnaTreev1::Init( PHCompositeNode * /*topNode*/ )
{
  
    if ( Verbosity () > 0 ) 
    {
        std::cout << "AnaTreev1::Init - opening file " << m_output_filename << std::endl;
    }

    // init tree
    PHTFileServer::get().open( m_output_filename, "RECREATE" );

    // h_centrality = new TH1D("h_centrality",";Centrality Bin;Events", 100, -0.5, 99.5);
    // for (int i = 0; i < 100; i++)
    // {
    //   h_jet_spectra[i] = new TH1D(Form("h_jet_spectra_%d", i),";p_{T};Counts", 100, 0, 50);
    //   h_jet_spectra_etacut[i] = new TH1D(Form("h_jet_spectra_etacut_%d",i),";p_{T};Counts", 100, 0, 50);
    // }

    m_tree = new TTree( "T", "T" );
    m_event_id = -1;
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

    if ( !m_mbd_node.empty() )
    {
        m_tree -> Branch( "mbd_q_N", &m_mbd_q_N, "mbd_q_N/F" );
        m_tree -> Branch( "mbd_q_S", &m_mbd_q_S, "mbd_q_S/F" );
        m_tree -> Branch( "mbd_t_N", &m_mbd_t_N, "mbd_t_N/F" );
        m_tree -> Branch( "mbd_t_S", &m_mbd_t_S, "mbd_t_S/F" );
    }

    if ( !m_truth_jet_node.empty() ) 
    {
        m_tree  -> Branch( "truth_jet_E", &m_truth_jet_E );
        m_tree  -> Branch( "truth_jet_phi", &m_truth_jet_phi );
        m_tree  -> Branch( "truth_jet_eta", &m_truth_jet_eta );
        m_tree  -> Branch( "truth_jet_pT", &m_truth_jet_pT );
        if ( !m_phHep_node.empty() ) 
        {
            m_tree -> Branch( "truth_jet_flavor", &m_truth_jet_flavor );
            m_tree -> Branch( "truth_jet_parton_pT", &m_truth_jet_parton_pT );
            m_tree -> Branch( "truth_zvtx", &m_truth_zvtx, "truth_zvtx/F" );
            m_tree -> Branch( "truth_jet_R", &m_truth_jet_R, "truth_jet_R/F" );
            m_tree -> Branch( "truth_jet_maxpT_r04", &m_truth_jet_maxpT_r04, "truth_jet_maxpT_r04/F" );
        }
    }

    if ( !m_cemc_node.empty() )
    {
        m_tree->Branch( "sumeT_cemc", &m_sumeT_cemc, "sumeT_cemc/F" );

        if ( m_save_full_cemc )
        {
            m_tree->Branch(
                "cemc_tower_E",
                m_cemc_tower_E,
                "cemc_tower_E[24][64]/F"
            );

            m_tree->Branch(
                "cemc_tower_isgood",
                m_cemc_tower_isgood,
                "cemc_tower_isgood[24][64]/I"
            );
        }
    }

    if ( !m_ihcal_node.empty() )
    {
        m_tree->Branch( "sumeT_ihcal", &m_sumeT_ihcal, "sumeT_ihcal/F" );

        if ( m_save_full_ihcal )
        {
            m_tree->Branch(
                "ihcal_tower_E",
                m_ihcal_tower_E,
                "ihcal_tower_E[24][64]/F"
            );

            m_tree->Branch(
                "ihcal_tower_isgood",
                m_ihcal_tower_isgood,
                "ihcal_tower_isgood[24][64]/I"
            );
        }
    }

    if ( !m_ohcal_node.empty() )
    {
        m_tree->Branch( "sumeT_ohcal", &m_sumeT_ohcal, "sumeT_ohcal/F" );

        if ( m_save_full_ohcal )
        {
            m_tree->Branch(
                "ohcal_tower_E",
                m_ohcal_tower_E,
                "ohcal_tower_E[24][64]/F"
            );

            m_tree->Branch(
                "ohcal_tower_isgood",
                m_ohcal_tower_isgood,
                "ohcal_tower_isgood[24][64]/I"
            );
        }
    }

    if ( !m_jet_node.empty() )
    {
        m_tree  -> Branch( "jet_R", &m_jet_R, "jet_R/F" );
        m_tree  -> Branch( "jet_E", &m_jet_E );
        m_tree  -> Branch( "jet_phi", &m_jet_phi );
        m_tree  -> Branch( "jet_eta", &m_jet_eta );
        m_tree  -> Branch( "jet_pT", &m_jet_pT );
    }

    if ( !m_sub1_jet_node.empty() )
    {
        m_tree  -> Branch( "sub1_jet_R", &m_sub1_jet_R, "sub1_jet_R/F" );
        m_tree  -> Branch( "sub1_jet_E", &m_sub1_jet_E );
        m_tree  -> Branch( "sub1_jet_phi", &m_sub1_jet_phi );
        m_tree  -> Branch( "sub1_jet_eta", &m_sub1_jet_eta );
        m_tree  -> Branch( "sub1_jet_pT", &m_sub1_jet_pT );
        m_tree  -> Branch( "sub1_jet_unsub_pT", &m_sub1_jet_unsub_pT );
        m_tree  -> Branch( "sub1_jet_unsub_E", &m_sub1_jet_unsub_E );
        m_tree  -> Branch( "sub1_jet_constituent_E", &m_sub1_jet_constituent_E );
        m_tree  -> Branch( "sub1_jet_constituent_phi", &m_sub1_jet_constituent_phi );
        m_tree  -> Branch( "sub1_jet_constituent_eta", &m_sub1_jet_constituent_eta );
        m_tree  -> Branch( "sub1_jet_constituent_pT", &m_sub1_jet_constituent_pT );
        m_tree  -> Branch( "sub1_jet_constituent_srcID", &m_sub1_jet_constituent_srcID );
    }

    if ( !m_towerbkgd_v2_node.empty() )
    {
        m_tree -> Branch( "sub2_v2", &m_sub2_v2, "sub2_v2/F" );
        m_tree -> Branch( "sub2_flowfaliure", &m_sub2_flowfaliure, "sub2_flowfaliure/I" );
        m_tree -> Branch( "sub2_psi2", &m_sub2_psi2, "sub2_psi2/F" );
        m_tree -> Branch( "sub2_towerbkgd_ue", &m_sub2_towerbkgd_ue );
    }
    
    if ( !m_rho_jet_node.empty() ) 
    {
        m_tree  -> Branch( "rho_jet_R", &m_rho_jet_R, "rho_jet_R/F" );
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

    
    if ( m_rho_nodes.size() > 0 )
    {
        m_rho_vals.resize( m_rho_nodes.size(), 0.0 );
        m_rho_sigmas.resize( m_rho_nodes.size(), 0.0 );
        for ( size_t i = 0; i < m_rho_nodes.size(); ++i )
        {
            m_tree -> Branch( Form( "rho_val_%s", m_rho_nodes[i].c_str() ), &m_rho_vals[i], Form( "rho_val_%s/F", m_rho_nodes[i].c_str() ) );
            m_tree -> Branch( Form( "rho_sigma_%s", m_rho_nodes[i].c_str() ), &m_rho_sigmas[i], Form( "rho_sigma_%s/F", m_rho_nodes[i].c_str() ) );
        }
    }

  
    if ( Verbosity () > 0 )
    {
        std::cout << "AnaTreev1::Init - done" << std::endl;
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTreev1::InitRun ( PHCompositeNode */*topNode*/ )
{
    if ( Verbosity () > 0 ) 
    {
        std::cout << "AnaTreev1::InitRun - initializing run with file " << m_output_filename << std::endl;
    }

    m_event_id = -1;

    return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTreev1::process_event( PHCompositeNode *topNode )
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
        // h_centrality ->Fill(m_cent*1.0);
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

    if ( !m_truth_jet_node.empty() ) 
    { 
        // get truth jet info
        m_truth_jet_E.clear();
        m_truth_jet_phi.clear();
        m_truth_jet_eta.clear();
        m_truth_jet_pT.clear();
        m_truth_jet_flavor.clear();
        m_truth_jet_parton_pT.clear();
        m_truth_jet_R = -1;
        m_truth_jet_maxpT_r04 = -1;

        auto * truth_jets = findNode::getClass<JetContainer>( topNode, m_truth_jet_node );
        if ( !truth_jets )        
        {
            std::cout << PHWHERE << " Input node " << m_truth_jet_node << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN; 
        }
        m_truth_jet_R = truth_jets->get_par();

        if ( m_truth_jet_R != 0.4f )
        {
            auto * truth_jets_r04 = findNode::getClass<JetContainer>( topNode, "AntiKt_Truth_r04" );
            if ( truth_jets_r04 )
            {
                float maxpT = -1;
                for ( const auto & jet : * truth_jets_r04 )
                {
                    if ( jet->get_pt() > maxpT )
                    {
                        maxpT = jet->get_pt();
                    }
                }
                m_truth_jet_maxpT_r04 = maxpT;
            }
        }
        else
        {
            float maxpT = -1;
            for ( const auto & jet : * truth_jets )
            {
                if ( jet->get_pt() > maxpT )
                {
                    maxpT = jet->get_pt();
                }
            }
            m_truth_jet_maxpT_r04 = maxpT;
        }

        PHHepMCGenEventMap * geneventmap = nullptr;
        HepMC::GenEvent * hepmc_event = nullptr;
        if ( !m_phHep_node.empty() ) 
        {
            
            m_truth_zvtx = -999;
            geneventmap = findNode::getClass<PHHepMCGenEventMap>( topNode, m_phHep_node );
            if ( !geneventmap )
            {
                std::cout << PHWHERE << " Input node " << m_phHep_node << " Node missing, doing nothing." << std::endl;
                return Fun4AllReturnCodes::ABORTRUN; 
            }
            
            PHHepMCGenEvent * genevt = geneventmap -> get(2);
            if ( genevt ){ hepmc_event = genevt->getEvent(); }
        }

        // walk the shower history (ISR+FSR) forward from a hard-process
        // parton to its last same-flavor copy, i.e. the parton as it
        // enters hadronization. the raw hard-process (status 21/22/23)
        // momentum is defined at the 2->2 vertex, before the backward
        // ISR evolution on the incoming legs finishes reshuffling
        // momentum/frame across the whole event -- matching jets
        // (which live in the final, fully-showered lab frame) against
        // that pre-ISR momentum directly can be systematically
        // mismatched. at each branching, follow the highest-pT
        // daughter that keeps the same |pdg_id|; stop once there's no
        // such daughter left (that particle is the final parton copy).
        auto find_final_parton = []( HepMC::GenParticle * p ) -> HepMC::GenParticle *
        {
            const int pid = abs( p->pdg_id() );
            while ( p && p->end_vertex() )
            {
                HepMC::GenParticle * next = nullptr;
                float next_pt = -1.0;
                for (
                    HepMC::GenVertex::particles_out_const_iterator d = p->end_vertex()->particles_out_const_begin();
                    d != p->end_vertex()->particles_out_const_end();
                    ++d
                )
                {
                    if ( !(*d) || abs( (*d)->pdg_id() ) != pid ) continue;
                    const HepMC::FourVector mom = (*d)->momentum();
                    const float dpt = sqrt( mom.px() * mom.px() + mom.py() * mom.py() );
                    if ( dpt > next_pt )
                    {
                        next_pt = dpt;
                        next = *d;
                    }
                }
                if ( !next ) break; // no same-flavor daughter left -- p is the final copy
                p = next;
            }
            return p;
        };


        for ( const auto & jet : * truth_jets )
        {
            // m_truth_jet_E.push_back(jet->get_e());
            // m_truth_jet_eta.push_back(jet->get_eta());
            // m_truth_jet_phi.push_back(jet->get_phi());
            // m_truth_jet_pT.push_back(jet->get_pt());

            float e = jet->get_e();
            float eta = jet->get_eta();
            float phi = jet->get_phi();
            float pt = jet->get_pt();

            int flavor = -1;
            float max_pt = 0;
            for (
                HepMC::GenEvent::particle_const_iterator p = hepmc_event->particles_begin();
		        p != hepmc_event->particles_end();
                ++p
            )
		    {
                HepMC::GenParticle * particle = *p;
		        if (!particle) continue;

                int pid = abs(particle->pdg_id());
                int status = particle->status();

                // Select outgoing partons from hard scattering (status 23) or
                // partons before hadronization (status 21, 22)
                if (status != 23 && status != 21 && status != 22 ) continue;

                // Only consider quarks (1-6) and gluons (21)
                bool is_quark = (pid >= 1 && pid <= 6);
                bool is_gluon = (pid == 21);
                if (!is_quark && !is_gluon) continue;

                // if (!(pid >= 1 && pid <= 6) && pid != 21) continue;

                HepMC::GenParticle * final_parton = find_final_parton( particle );
                HepMC::FourVector momentum = final_parton->momentum();
                float part_pt = sqrt(momentum.px() * momentum.px() + momentum.py() * momentum.py());
                float part_eta = momentum.eta();
                float part_phi = momentum.phi();
                
                // Require minimum pT for parton matching
                if ( part_pt < 1.0 ) continue;

                // Calculate angular distance between parton and jet
                if ( part_phi > TMath::Pi() )
                {
                    part_phi -= 2*TMath::Pi();
                }

                float deta = fabs(eta - part_eta);
                float dphi = fabs(phi - part_phi);
                if (dphi > TMath::Pi())
                {
                    dphi -= 2*TMath::Pi();
                }

                float dr = sqrt(deta*deta + dphi*dphi);

                // Match closest parton with highest pT
                if (dr < 0.4 && part_pt > max_pt)
                {
                    max_pt = part_pt;
                    flavor = pid;
                }
		    }

            m_truth_jet_E.push_back(e);
            m_truth_jet_eta.push_back(eta);
            m_truth_jet_phi.push_back(phi);
            m_truth_jet_pT.push_back(pt);
            m_truth_jet_flavor.push_back(flavor);
            m_truth_jet_parton_pT.push_back(max_pt);
        }

    }

    if ( !m_sub1_jet_node.empty() ) 
    { 
        
        m_sub1_jet_R = -1;
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

        m_sub1_jet_towerbkgd_ue.clear();
        m_sub1_jet_towerbkgd_v2 = 0.0;
        m_sub1_jet_towerbkgd_flowfaliure = 0;
        m_sub1_jet_towerbkgd_psi2 = 0.0;

        // get sub1 jet info
        auto * jets = findNode::getClass<JetContainer>( topNode, m_sub1_jet_node );
        if ( !jets )
        {
            std::cout << PHWHERE << " Input node " << m_sub1_jet_node << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN; 
        }
        m_sub1_jet_R = jets->get_par();
        
    

        auto * tower_background_sub2 = findNode::getClass<TowerBackgroundv1>(topNode, m_sub1_jet_towerbkgd_node);
        if ( !tower_background_sub2  )
        {
            std::cout << PHWHERE << " TowerBackgroundv1 node is missing, skipping." << std::endl;
            return Fun4AllReturnCodes::ABORTEVENT;
        }
        else 
        {
            // read info
            m_sub1_jet_towerbkgd_v2 = tower_background_sub2->get_v2();
            m_sub1_jet_towerbkgd_flowfaliure = tower_background_sub2->get_flow_failure_flag();
            m_sub1_jet_towerbkgd_psi2 = tower_background_sub2->get_Psi2();
            for ( int i = 0; i < 3; ++i )
            {
                m_sub1_jet_towerbkgd_ue.push_back( tower_background_sub2->get_UE(i) );
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
                else if ( comp.first == Jet::SRC::HCALIN_TOWERINFO )
                {
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
                double ue         = m_sub1_jet_towerbkgd_ue[layer_idx][ieta];
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
    if ( !m_towerbkgd_v2_node.empty() )
    {
        // get tower background info
        m_sub2_v2 = 0.0;
        m_sub2_flowfaliure = 0;
        m_sub2_psi2 = 0.0;
        m_sub2_towerbkgd_ue.clear();

        auto * tower_background_sub2 = findNode::getClass<TowerBackgroundv1>(topNode, m_towerbkgd_v2_node);
        if ( !tower_background_sub2  )
        {
            std::cout << PHWHERE << " TowerBackgroundv1 node is missing, skipping." << std::endl;
            return Fun4AllReturnCodes::ABORTEVENT;
        }
        else 
        {
            // read info
            m_sub2_v2 = tower_background_sub2->get_v2();
            m_sub2_flowfaliure = tower_background_sub2->get_flow_failure_flag();
            m_sub2_psi2 = tower_background_sub2->get_Psi2();
            for ( int i = 0; i < 3; ++i )
            {
                m_sub2_towerbkgd_ue.push_back( tower_background_sub2->get_UE(i) );
            }
        }
    }

   
    if ( !m_rho_jet_node.empty() )
    { 
        
        m_rho_jet_R = -1;
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

        m_rho_jet_cemc_rho = 0.0;
        m_rho_jet_hcalin_rho = 0.0;
        m_rho_jet_hcalout_rho = 0.0;
        
        // get rho jet info
        auto * jets = findNode::getClass<JetContainer>( topNode, m_rho_jet_node );
        if ( !jets )
        {
            std::cout << PHWHERE << " Input node " << m_rho_jet_node << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN; 
        }
        m_rho_jet_R = jets->get_par();
        
    
        auto * cemc_rho_node = findNode::getClass<TowerRhov1>( topNode,m_rho_jet_cemc_rho_node );
        auto * hcalin_rho_node = findNode::getClass<TowerRhov1>( topNode,m_rho_jet_hcalin_rho_node );
        auto * hcalout_rho_node = findNode::getClass<TowerRhov1>( topNode,m_rho_jet_hcalout_rho_node );
        if ( !cemc_rho_node || !hcalin_rho_node || !hcalout_rho_node )
        {
            std::cout << PHWHERE << " One or more TowerRhov1 nodes are missing, skipping." << std::endl;
            return Fun4AllReturnCodes::ABORTEVENT;
        }
        m_rho_jet_cemc_rho = cemc_rho_node->get_rho();
        m_rho_jet_hcalin_rho = hcalin_rho_node->get_rho();
        m_rho_jet_hcalout_rho = hcalout_rho_node->get_rho();

        

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
                    else if ( comp.first == Jet::SRC::HCALIN_TOWERINFO )
                    {
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
                double ue = 0;
                if ( layer_idx == 0 )
                {
                    ue = m_rho_jet_cemc_rho;
                }
                else if ( layer_idx == 1 )
                {
                    ue = m_rho_jet_hcalin_rho;
                }
                else if ( layer_idx == 2 )
                {
                    ue = m_rho_jet_hcalout_rho;
                }
                ue *= cosh(comp_eta);
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

    for ( size_t irho = 0; irho < m_rho_nodes.size(); ++irho )
    {
        auto * rho_node = findNode::getClass<TowerRhov1>( topNode, m_rho_nodes[irho] );
        if ( !rho_node )
        {
            std::cout << PHWHERE << " Input node " << m_rho_nodes[irho] << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN; 
        }
        m_rho_vals[irho] = rho_node->get_rho();
        m_rho_sigmas[irho] = rho_node->get_sigma();
    }

    if ( !m_cemc_node.empty() )
    {
        m_sumeT_cemc = 0.0;
        if ( m_save_full_cemc )
        {
            memset( m_cemc_tower_E, 0, sizeof(m_cemc_tower_E) );
            memset( m_cemc_tower_isgood, 0, sizeof(m_cemc_tower_isgood) );
        }

        auto * cemc_node = findNode::getClass<TowerInfoContainer>( topNode, m_cemc_node );
        if ( !cemc_node )
        {
            std::cout << PHWHERE << " Input node " << m_cemc_node << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }

        for ( unsigned int ich = 0; ich < cemc_node->size(); ich++ )
        {
            auto tower = cemc_node->get_tower_at_channel(ich);
            if ( !tower || !tower->get_isGood() || std::isnan(tower->get_energy()) )
            {
                continue; // skip bad towers
            }
            unsigned int key = cemc_node -> encode_key(ich);
            int ieta = cemc_node -> getTowerEtaBin(key);
            int iphi = cemc_node -> getTowerPhiBin(key);
            m_sumeT_cemc += tower->get_energy();
            if ( m_save_full_cemc )
            {
                m_cemc_tower_E[ieta][iphi] = tower->get_energy();
                m_cemc_tower_isgood[ieta][iphi] = 1;
            }
        }
    }

    if ( !m_ihcal_node.empty() )
    {
        m_sumeT_ihcal = 0.0;
        if ( m_save_full_ihcal )
        {
            memset( m_ihcal_tower_E, 0, sizeof(m_ihcal_tower_E) );
            memset( m_ihcal_tower_isgood, 0, sizeof(m_ihcal_tower_isgood) );
        }

        auto * ihcal_node = findNode::getClass<TowerInfoContainer>( topNode, m_ihcal_node );
        if ( !ihcal_node )
        {
            std::cout << PHWHERE << " Input node " << m_ihcal_node << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }

        for ( unsigned int ich = 0; ich < ihcal_node->size(); ich++ )
        {
            auto tower = ihcal_node->get_tower_at_channel(ich);
            if ( !tower || !tower->get_isGood() || std::isnan(tower->get_energy()) )
            {
                continue; // skip bad towers
            }
            unsigned int key = ihcal_node -> encode_key(ich);
            int ieta = ihcal_node -> getTowerEtaBin(key);
            int iphi = ihcal_node -> getTowerPhiBin(key);
            m_sumeT_ihcal += tower->get_energy();
            if ( m_save_full_ihcal )
            {
                m_ihcal_tower_E[ieta][iphi] = tower->get_energy();
                m_ihcal_tower_isgood[ieta][iphi] = 1;
            }
        }
    }

    if ( !m_ohcal_node.empty() )
    {
        m_sumeT_ohcal = 0.0;
        if ( m_save_full_ohcal )
        {
            memset( m_ohcal_tower_E, 0, sizeof(m_ohcal_tower_E) );
            memset( m_ohcal_tower_isgood, 0, sizeof(m_ohcal_tower_isgood) );
        }

        auto * ohcal_node = findNode::getClass<TowerInfoContainer>( topNode, m_ohcal_node );
        if ( !ohcal_node )
        {
            std::cout << PHWHERE << " Input node " << m_ohcal_node << " Node missing, doing nothing." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }

        for ( unsigned int ich = 0; ich < ohcal_node->size(); ich++ )
        {
            auto tower = ohcal_node->get_tower_at_channel(ich);
            if ( !tower || !tower->get_isGood() || std::isnan(tower->get_energy()) )
            {
                continue; // skip bad towers
            }
            unsigned int key = ohcal_node -> encode_key(ich);
            int ieta = ohcal_node -> getTowerEtaBin(key);
            int iphi = ohcal_node -> getTowerPhiBin(key);
            m_sumeT_ohcal += tower->get_energy();
            if ( m_save_full_ohcal )
            {
                m_ohcal_tower_E[ieta][iphi] = tower->get_energy();
                m_ohcal_tower_isgood[ieta][iphi] = 1;
            }
        }
    }

    
    // fill tree
    m_tree->Fill();
    
    return Fun4AllReturnCodes::EVENT_OK;

}

int AnaTreev1::End( PHCompositeNode * /*topNode*/ )
{
  
  if( Verbosity() > 0 ) 
  {
    std::cout << "AnaTreev1::EndRun - End run " << std::endl;
    std::cout << "AnaTreev1::EndRun - Writing to " << m_output_filename << std::endl;
  }

  PHTFileServer::get().cd(m_output_filename); 

  m_tree->Write();
  

  if ( Verbosity() > 0 ) 
  {
    std::cout << "AnaTreev1::EndRun - Writing run tree" << std::endl;
  }
  PHTFileServer::get().close();
 
  if ( Verbosity () > 0 ) 
  {
    std::cout << "AnaTreev1::EndRun - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

TowerInfoContainer * AnaTreev1::LoadTowerInfoContainer(PHCompositeNode *topNode, const std::string &tower_node_name)
{
  auto * tic = findNode::getClass<TowerInfoContainer>(topNode, tower_node_name.c_str());
  if (!tic)
  {
    std::cout << "AnaTreev1::LoadTowerInfoContainer - cannot find " << tower_node_name << ", exiting" << std::endl;
    exit(1);
  }
  return tic;
}

RawTowerGeomContainer * AnaTreev1::LoadTowerGeomContainer(PHCompositeNode *topNode, const std::string &tower_geom_node_name)
{
  auto * tg = findNode::getClass<RawTowerGeomContainer>(topNode, tower_geom_node_name.c_str());
  if (!tg)
  {
    std::cout << "AnaTreev1::LoadTowerGeomContainer - cannot find " << tower_geom_node_name << ", exiting" << std::endl;
    exit(1);
  }
  return tg;
}








   