
#ifndef _MATCH_STANDALONE_C_
#define _MATCH_STANDALONE_C_

#include <myana/AnaUtils.h>

#include <TBranch.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )


int match_standalone(
    const std::string & infile = "/sphenix/user/tmengel/dijet-ana-auau/macros/truth-matching/rootfiles-merged/08_24_2026_v001/jet10.list",
    const std::string & outfile = "output.root"
)
{
    auto * t = new TChain( "T" );

    auto files = AnaUtils::getFilelist( infile, ".root" );
    for ( const auto & file : files )
    {
        std::cout << "Adding file: " << file << std::endl;
        t -> Add( file.c_str() );
    }

    int nentries = t -> GetEntries();
    if ( nentries < 1 )
    {
        std::cerr << "Error: No entries in TChain." << std::endl;
        return -1;
    }
    std::cout << "Total entries in TChain: " << nentries << std::endl;

    t -> SetBranchStatus( "*", false );

    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * tout = new TTree( "T", "T" );

    auto has_branch = [ t ]( const std::string & name )
    {
        return t -> GetBranch( name.c_str() ) != nullptr;
    };

    // Jet::SRC codes (jetbase/Jet.h): CEMC-layer constituents show up as 25
    // (CEMC_TOWERINFO), 28 (CEMC_TOWERINFO_RETOWER) or 29 (CEMC_TOWERINFO_SUB1);
    // HCALIN as 26 or 30 (HCALIN_TOWERINFO[_SUB1]); HCALOUT as 27 or 31
    // (HCALOUT_TOWERINFO[_SUB1]) -- same mapping as makeMatchedTreesTaggedCaloAuAu.C.
    auto calc_calo_frac = []( 
        const std::vector< float > & compE, 
        const std::vector< int > & compSrc,
        float & emcalfrac, float & ihcalfrac, float & ohcalfrac 
    )
    {
        emcalfrac = 0.0;
        ihcalfrac = 0.0;
        ohcalfrac = 0.0;
        float total = 0.0;
        const size_t n = std::min( compE.size(), compSrc.size() );
        for ( size_t i = 0; i < n; ++i )
        {
            const float e = compE[i];
            total += e;
            const int src = compSrc[i];
            if ( src == 25 || src == 28 || src == 29 ) emcalfrac += e;
            else if ( src == 26 || src == 30 ) ihcalfrac += e;
            else if ( src == 27 || src == 31 ) ohcalfrac += e;
        }
        if ( total > 0.0 )
        {
            emcalfrac /= total;
            ihcalfrac /= total;
            ohcalfrac /= total;
        }
    };

    //----------------------------------------------------------------
    // event_id
    //----------------------------------------------------------------
    int event_id = -1;
    if ( has_branch( "event_id" ) )
    {
        t    -> SetBranchStatus ( "event_id", true );
        t    -> SetBranchAddress( "event_id", &event_id );
        tout -> Branch( "event_id", &event_id, "event_id/I" );
    }

    //----------------------------------------------------------------
    // gl1 trigger vectors
    //----------------------------------------------------------------
    // not read: unused (GL1 bit-10 cut below is disabled, and neither
    // vector is written to output). Re-enable both lines below together
    // with the cut in the event loop if the trigger requirement returns.
    static const int k_gl1_max = 64;
    int scaled_triggervec[k_gl1_max] {};
    int live_triggervec[k_gl1_max] {};
    // if ( has_branch( "scaled_triggervec" ) )
    // {
    //     t -> SetBranchStatus ( "scaled_triggervec", true );
    //     t -> SetBranchAddress( "scaled_triggervec", scaled_triggervec );
    // }

    //----------------------------------------------------------------
    // minbias flag
    //----------------------------------------------------------------
    int is_minbias = 0;
    if ( has_branch( "is_minbias" ) )
    {
        t    -> SetBranchStatus ( "is_minbias", true );
        t    -> SetBranchAddress( "is_minbias", &is_minbias );

        tout -> Branch( "is_minbias", &is_minbias, "is_minbias/I" );
    }

    //----------------------------------------------------------------
    // z-vertex
    //----------------------------------------------------------------
    float zvrtx = 0.0;
    if ( has_branch( "zvrtx" ) )
    {
        t    -> SetBranchStatus ( "zvrtx", true );
        t    -> SetBranchAddress( "zvrtx", &zvrtx );

        tout -> Branch( "zvrtx", &zvrtx, "zvrtx/F" );
    }

    //----------------------------------------------------------------
    // centrality
    //----------------------------------------------------------------
    int cent = -1;
    if ( has_branch( "cent" ) )
    {
        t    -> SetBranchStatus ( "cent", true );
        t    -> SetBranchAddress( "cent", &cent );

        tout -> Branch( "cent", &cent, "cent/I" );
    }

    //----------------------------------------------------------------
    // event header (HepMC)
    //----------------------------------------------------------------
    // only psi2/psi3 are written to output -- b/ep_angle/ecc/psi1/ncoll/
    // npart/runnumber/evtsequence are not read at all.
    float psi2 = 0.0, psi3 = 0.0;
    if ( has_branch( "psi2" ) )
    {
        t -> SetBranchStatus ( "psi2", true );
        t -> SetBranchStatus ( "psi3", true );

        t -> SetBranchAddress( "psi2", &psi2 );
        t -> SetBranchAddress( "psi3", &psi3 );

        tout -> Branch( "psi2", &psi2, "psi2/F" );
        tout -> Branch( "psi3", &psi3, "psi3/F" );
    }

    //----------------------------------------------------------------
    // mbd
    //----------------------------------------------------------------
    // mbd_t_N/mbd_t_S are not read: only mbd_q_N + mbd_q_S feed the
    // derived mbd_q output.
    float mbd_q_N = -999.0, mbd_q_S = -999.0;
    float mbd_q = -999.0;
    if ( has_branch( "mbd_q_N" ) )
    {
        t -> SetBranchStatus ( "mbd_q_N", true );
        t -> SetBranchStatus ( "mbd_q_S", true );

        t -> SetBranchAddress( "mbd_q_N", &mbd_q_N );
        t -> SetBranchAddress( "mbd_q_S", &mbd_q_S );

        tout -> Branch( "mbd_q", &mbd_q, "mbd_q/F" );
    }

    //----------------------------------------------------------------
    // truth jets
    //----------------------------------------------------------------
    std::vector< float > * truth_jet_E   = nullptr;
    std::vector< float > * truth_jet_phi = nullptr;
    std::vector< float > * truth_jet_eta = nullptr;
    std::vector< float > * truth_jet_pT  = nullptr;
    std::vector< int > truth_jet_accept_eta {};
    const bool has_truth_jets = has_branch( "truth_jet_E" );
    if ( has_truth_jets )
    {
        t -> SetBranchStatus ( "truth_jet_E", true );
        t -> SetBranchStatus ( "truth_jet_phi", true );
        t -> SetBranchStatus ( "truth_jet_eta", true );
        t -> SetBranchStatus ( "truth_jet_pT", true );

        t -> SetBranchAddress( "truth_jet_E", &truth_jet_E );
        t -> SetBranchAddress( "truth_jet_phi", &truth_jet_phi );
        t -> SetBranchAddress( "truth_jet_eta", &truth_jet_eta );
        t -> SetBranchAddress( "truth_jet_pT", &truth_jet_pT );

        tout -> Branch( "truth_jet_E", &truth_jet_E );
        tout -> Branch( "truth_jet_phi", &truth_jet_phi );
        tout -> Branch( "truth_jet_eta", &truth_jet_eta );
        tout -> Branch( "truth_jet_pT", &truth_jet_pT );
        tout -> Branch( "truth_jet_accept_eta", &truth_jet_accept_eta );
    }

    // truth jet parton matching (only present if PHHepMCGenEventMap was added)
    std::vector< int >   * truth_jet_flavor    = nullptr;
    std::vector< float > * truth_jet_parton_pT = nullptr;
    std::vector< float > * truth_jet_parton_dr = nullptr;
    float truth_zvtx = 0.0, truth_jet_R = 0.0, truth_jet_maxpT_r04 = 0.0;
    const bool has_truth_jet_partons = has_branch( "truth_jet_flavor" );
    if ( has_truth_jet_partons )
    {
        t -> SetBranchStatus ( "truth_jet_flavor", true );
        t -> SetBranchStatus ( "truth_jet_parton_pT", true );
        t -> SetBranchStatus ( "truth_jet_parton_dr", true );
        t -> SetBranchStatus ( "truth_zvtx", true );
        t -> SetBranchStatus ( "truth_jet_R", true );
        t -> SetBranchStatus ( "truth_jet_maxpT_r04", true );

        t -> SetBranchAddress( "truth_jet_flavor", &truth_jet_flavor );
        t -> SetBranchAddress( "truth_jet_parton_pT", &truth_jet_parton_pT );
        t -> SetBranchAddress( "truth_jet_parton_dr", &truth_jet_parton_dr );
        t -> SetBranchAddress( "truth_zvtx", &truth_zvtx );
        t -> SetBranchAddress( "truth_jet_R", &truth_jet_R );
        t -> SetBranchAddress( "truth_jet_maxpT_r04", &truth_jet_maxpT_r04 );

        tout -> Branch( "truth_jet_flavor", &truth_jet_flavor );
        tout -> Branch( "truth_jet_parton_pT", &truth_jet_parton_pT );
        tout -> Branch( "truth_jet_parton_dr", &truth_jet_parton_dr );
        tout -> Branch( "truth_zvtx", &truth_zvtx, "truth_zvtx/F" );
        tout -> Branch( "truth_jet_R", &truth_jet_R, "truth_jet_R/F" );
        tout -> Branch( "truth_jet_maxpT_r04", &truth_jet_maxpT_r04, "truth_jet_maxpT_r04/F" );
    }

    // re-sorts the truth jet vectors in place, descending by pT, so the
    // output (and any indices into it, e.g. truth_jet_reco_match_idx/
    // jet_truth_match_idx) reflect pT order.
    auto sort_truth_jets_by_pt = [&]()
    {
        const size_t n = truth_jet_pT -> size();
        std::vector< size_t > idx( n );
        std::iota( idx.begin(), idx.end(), 0 );
        std::sort( idx.begin(), idx.end(), [&]( size_t a, size_t b )
        {
            return truth_jet_pT -> at( a ) > truth_jet_pT -> at( b );
        } );

        auto apply_f = [&]( std::vector< float > & v )
        {
            std::vector< float > tmp( n );
            for ( size_t k = 0; k < n; ++k ) tmp[k] = v[ idx[k] ];
            v = std::move( tmp );
        };
        auto apply_i = [&]( std::vector< int > & v )
        {
            std::vector< int > tmp( n );
            for ( size_t k = 0; k < n; ++k ) tmp[k] = v[ idx[k] ];
            v = std::move( tmp );
        };

        apply_f( *truth_jet_E );
        apply_f( *truth_jet_phi );
        apply_f( *truth_jet_eta );
        apply_f( *truth_jet_pT );

        if ( has_truth_jet_partons )
        {
            apply_i( *truth_jet_flavor );
            apply_f( *truth_jet_parton_pT );
            apply_f( *truth_jet_parton_dr );
        }
    };

    //----------------------------------------------------------------
    // calo towers (cemc / ihcal / ohcal)
    //----------------------------------------------------------------
    static const int k_ieta = 24;
    static const int k_iphi = 64;

    // the raw sumeT_cemc/ihcal/ohcal scalar branches AnaTreev1 wrote directly
    // are not read: superseded by the zvrtx-corrected *_calc versions below.
    float cemc_tower_E[k_ieta][k_iphi] {};
    int   cemc_tower_isgood[k_ieta][k_iphi] {};
    if ( has_branch( "cemc_tower_E" ) )
    {
        t -> SetBranchStatus ( "cemc_tower_E", true );
        t -> SetBranchStatus ( "cemc_tower_isgood", true );
        t -> SetBranchAddress( "cemc_tower_E", cemc_tower_E );
        t -> SetBranchAddress( "cemc_tower_isgood", cemc_tower_isgood );

        // tout -> Branch( "cemc_tower_E", cemc_tower_E, "cemc_tower_E[24][64]/F" );
        // tout -> Branch( "cemc_tower_isgood", cemc_tower_isgood, "cemc_tower_isgood[24][64]/I" );
    }

    float ihcal_tower_E[k_ieta][k_iphi] {};
    int   ihcal_tower_isgood[k_ieta][k_iphi] {};
    if ( has_branch( "ihcal_tower_E" ) )
    {
        t -> SetBranchStatus ( "ihcal_tower_E", true );
        t -> SetBranchStatus ( "ihcal_tower_isgood", true );
        t -> SetBranchAddress( "ihcal_tower_E", ihcal_tower_E );
        t -> SetBranchAddress( "ihcal_tower_isgood", ihcal_tower_isgood );

        // tout -> Branch( "ihcal_tower_E", ihcal_tower_E, "ihcal_tower_E[24][64]/F" );
        // tout -> Branch( "ihcal_tower_isgood", ihcal_tower_isgood, "ihcal_tower_isgood[24][64]/I" );
    }

    float ohcal_tower_E[k_ieta][k_iphi] {};
    int   ohcal_tower_isgood[k_ieta][k_iphi] {};
    if ( has_branch( "ohcal_tower_E" ) )
    {
        t -> SetBranchStatus ( "ohcal_tower_E", true );
        t -> SetBranchStatus ( "ohcal_tower_isgood", true );
        t -> SetBranchAddress( "ohcal_tower_E", ohcal_tower_E );
        t -> SetBranchAddress( "ohcal_tower_isgood", ohcal_tower_isgood );

        // tout -> Branch( "ohcal_tower_E", ohcal_tower_E, "ohcal_tower_E[24][64]/F" );
        // tout -> Branch( "ohcal_tower_isgood", ohcal_tower_isgood, "ohcal_tower_isgood[24][64]/I" );
    }

    //----------------------------------------------------------------
    // event-level sum eT, computed from the tower arrays the same way
    // match.C does (AnaUtils::calc_sumeT, zvrtx-corrected).
    //----------------------------------------------------------------
    float sumeT_cemc_calc  = 0.0;
    float sumeT_ihcal_calc = 0.0;
    float sumeT_ohcal_calc = 0.0;
    float sumeT_calc       = 0.0;
    const bool has_sumeT_calc_inputs = has_branch( "zvrtx" )
                                      && has_branch( "cemc_tower_E" )
                                      && has_branch( "ihcal_tower_E" )
                                      && has_branch( "ohcal_tower_E" );
    if ( has_sumeT_calc_inputs )
    {
        tout -> Branch( "sumeT_cemc", &sumeT_cemc_calc, "sumeT_cemc/F" );
        tout -> Branch( "sumeT_ihcal", &sumeT_ihcal_calc, "sumeT_ihcal/F" );
        tout -> Branch( "sumeT_ohcal", &sumeT_ohcal_calc, "sumeT_ohcal/F" );
        tout -> Branch( "sumeT", &sumeT_calc, "sumeT/F" );
    }

    // "reco jets (unsubtracted)" (jet_R/E/eta/phi/pT), the whole sub1
    // (seeded) subtracted jet collection (incl. its constituents), and
    // sub2/towerbkgd_v2 are not read: none of them feed an output branch
    // or a computation below (matching/em-fractions use the rho_jet_*
    // collection instead).
    float jet_R = 0.0;
    std::vector< float > * jet_E   = nullptr;
    std::vector< float > * jet_phi = nullptr;
    std::vector< float > * jet_eta = nullptr;
    std::vector< float > * jet_pT  = nullptr;

    float sub1_jet_R = 0.0;
    std::vector< float > * sub1_jet_E        = nullptr;
    std::vector< float > * sub1_jet_phi      = nullptr;
    std::vector< float > * sub1_jet_eta      = nullptr;
    std::vector< float > * sub1_jet_pT       = nullptr;
    std::vector< float > * sub1_jet_unsub_pT = nullptr;
    std::vector< float > * sub1_jet_unsub_E  = nullptr;

    std::vector< std::vector< float > > * sub1_jet_constituent_E   = nullptr;
    std::vector< std::vector< float > > * sub1_jet_constituent_phi = nullptr;
    std::vector< std::vector< float > > * sub1_jet_constituent_eta = nullptr;
    std::vector< std::vector< float > > * sub1_jet_constituent_pT  = nullptr;
    std::vector< std::vector< int > >   * sub1_jet_constituent_srcID = nullptr;

    float sub2_v2 = 0.0;
    int   sub2_flowfaliure = 0;
    float sub2_psi2 = 0.0;
    std::vector< std::vector< float > > * sub2_towerbkgd_ue = nullptr;

    //----------------------------------------------------------------
    // rho-subtracted jets
    //----------------------------------------------------------------
    float rho_jet_R = 0.0;
    std::vector< float > * rho_jet_E        = nullptr;
    std::vector< float > * rho_jet_phi      = nullptr;
    std::vector< float > * rho_jet_eta      = nullptr;
    std::vector< float > * rho_jet_pT       = nullptr;
    std::vector< float > * rho_jet_unsub_pT = nullptr;
    std::vector< float > * rho_jet_unsub_E  = nullptr;
    std::vector< int > jet_accept_eta {};
    const bool has_rho_jets = has_branch( "rho_jet_R" );
    if ( has_rho_jets )
    {
        t -> SetBranchStatus ( "rho_jet_R", true );
        t -> SetBranchStatus ( "rho_jet_E", true );
        t -> SetBranchStatus ( "rho_jet_phi", true );
        t -> SetBranchStatus ( "rho_jet_eta", true );
        t -> SetBranchStatus ( "rho_jet_pT", true );
        t -> SetBranchStatus ( "rho_jet_unsub_pT", true );
        t -> SetBranchStatus ( "rho_jet_unsub_E", true );

        t -> SetBranchAddress( "rho_jet_R", &rho_jet_R );
        t -> SetBranchAddress( "rho_jet_E", &rho_jet_E );
        t -> SetBranchAddress( "rho_jet_phi", &rho_jet_phi );
        t -> SetBranchAddress( "rho_jet_eta", &rho_jet_eta );
        t -> SetBranchAddress( "rho_jet_pT", &rho_jet_pT );
        t -> SetBranchAddress( "rho_jet_unsub_pT", &rho_jet_unsub_pT );
        t -> SetBranchAddress( "rho_jet_unsub_E", &rho_jet_unsub_E );

        tout -> Branch( "jet_R", &rho_jet_R, "rho_jet_R/F" );
        tout -> Branch( "jet_E", &rho_jet_E );
        tout -> Branch( "jet_phi", &rho_jet_phi );
        tout -> Branch( "jet_eta", &rho_jet_eta );
        tout -> Branch( "jet_pT", &rho_jet_pT );
        tout -> Branch( "jet_unsub_pT", &rho_jet_unsub_pT );
        tout -> Branch( "jet_unsub_E", &rho_jet_unsub_E );
        tout -> Branch( "jet_accept_eta", &jet_accept_eta );
    }

    // only constituent_E/constituent_srcID feed calc_calo_frac --
    // constituent_phi/eta/pT are not read.
    std::vector< std::vector< float > > * rho_jet_constituent_E   = nullptr;
    std::vector< std::vector< int > >   * rho_jet_constituent_srcID = nullptr;
    std::vector< float > rho_jet_cemcfrac {};
    std::vector< float > rho_jet_ihcalfrac{};
    std::vector< float > rho_jet_ohcalfrac{};

    const bool has_rho_jet_constituents = has_branch( "rho_jet_constituent_E" );
    if ( has_rho_jet_constituents )
    {
        t -> SetBranchStatus ( "rho_jet_constituent_E", true );
        t -> SetBranchStatus ( "rho_jet_constituent_srcID", true );

        t -> SetBranchAddress( "rho_jet_constituent_E", &rho_jet_constituent_E );
        t -> SetBranchAddress( "rho_jet_constituent_srcID", &rho_jet_constituent_srcID );

        tout -> Branch( "jet_cemcfrac", &rho_jet_cemcfrac );
        tout -> Branch( "jet_ihcalfrac", &rho_jet_ihcalfrac );
        tout -> Branch( "jet_ohcalfrac", &rho_jet_ohcalfrac );
    }

    // re-sorts the reco (rho-subtracted) jet vectors in place, descending
    // by pT -- including the constituent lists, so the em-fractions (and
    // any indices into this collection, e.g. truth_jet_reco_match_idx/
    // jet_truth_match_idx) stay aligned with the new order. Must run
    // before both the em-fraction calc and the truth/reco matching below.
    auto sort_reco_jets_by_pt = [&]()
    {
        const size_t n = rho_jet_pT -> size();
        std::vector< size_t > idx( n );
        std::iota( idx.begin(), idx.end(), 0 );
        std::sort( idx.begin(), idx.end(), [&]( size_t a, size_t b )
        {
            return rho_jet_pT -> at( a ) > rho_jet_pT -> at( b );
        } );

        auto apply_f = [&]( std::vector< float > & v )
        {
            std::vector< float > tmp( n );
            for ( size_t k = 0; k < n; ++k ) tmp[k] = v[ idx[k] ];
            v = std::move( tmp );
        };
        auto apply_vf = [&]( std::vector< std::vector< float > > & v )
        {
            std::vector< std::vector< float > > tmp( n );
            for ( size_t k = 0; k < n; ++k ) tmp[k] = std::move( v[ idx[k] ] );
            v = std::move( tmp );
        };
        auto apply_vi = [&]( std::vector< std::vector< int > > & v )
        {
            std::vector< std::vector< int > > tmp( n );
            for ( size_t k = 0; k < n; ++k ) tmp[k] = std::move( v[ idx[k] ] );
            v = std::move( tmp );
        };

        apply_f( *rho_jet_E );
        apply_f( *rho_jet_phi );
        apply_f( *rho_jet_eta );
        apply_f( *rho_jet_pT );
        apply_f( *rho_jet_unsub_pT );
        apply_f( *rho_jet_unsub_E );

        if ( has_rho_jet_constituents )
        {
            apply_vf( *rho_jet_constituent_E );
            apply_vi( *rho_jet_constituent_srcID );
        }
    };

    //----------------------------------------------------------------
    // truth jet -> reco jet ("jet_*", i.e. the rho-subtracted collection)
    // index match. No pT/eta selection on either side -- just the
    // highest-pT reco jet within dR < 0.75*0.3 of each truth jet.
    // -1 if no reco jet is within range.
    //----------------------------------------------------------------
    std::vector< int > truth_jet_reco_match_idx {};
    std::vector< int > jet_truth_match_idx {};
    const bool has_truth_reco_match_inputs = has_branch( "truth_jet_pT" ) && has_branch( "rho_jet_pT" );
    if ( has_truth_reco_match_inputs )
    {
        tout -> Branch( "truth_jet_reco_match_idx", &truth_jet_reco_match_idx );
        tout -> Branch( "jet_truth_match_idx", &jet_truth_match_idx );
    }

    //----------------------------------------------------------------
    // generic per-node rho values (rho_val_<node>, rho_sigma_<node>) --
    // branch names are dynamic, so wire these up in a loop instead of
    // one variable per node.
    //----------------------------------------------------------------
    std::vector< std::string > rho_names;
    for ( auto * obj : *t -> GetListOfBranches() )
    {
        std::string name = obj -> GetName();
        if ( name.rfind( "rho_val_", 0 ) == 0 )
        {
            rho_names.push_back( name.substr( std::string( "rho_val_" ).size() ) );
        }
    }

    std::vector< float > rho_vals( rho_names.size(), 0.0 );
    std::vector< float > rho_sigmas( rho_names.size(), 0.0 );
    for ( size_t i = 0; i < rho_names.size(); ++i )
    {
        const std::string val_name   = "rho_val_" + rho_names[i];
        const std::string sigma_name = "rho_sigma_" + rho_names[i];

        t -> SetBranchStatus ( val_name.c_str(), true );
        t -> SetBranchStatus ( sigma_name.c_str(), true );
        t -> SetBranchAddress( val_name.c_str(), &rho_vals[i] );
        t -> SetBranchAddress( sigma_name.c_str(), &rho_sigmas[i] );

        tout -> Branch( val_name.c_str(), &rho_vals[i], ( val_name + "/F" ).c_str() );
        tout -> Branch( sigma_name.c_str(), &rho_sigmas[i], ( sigma_name + "/F" ).c_str() );
    }

    //----------------------------------------------------------------
    // event selection: minbias-triggered, |zvrtx| < 60 cm, GL1 bit 10 fired
    //----------------------------------------------------------------
    const bool has_event_sel_inputs = has_branch( "is_minbias" )
                                     && has_branch( "zvrtx" );
    if ( !has_event_sel_inputs )
    {
        std::cerr << "Warning: is_minbias/zvrtx not all present "
                   << "-- event selection will not be applied." << std::endl;
    }

    //----------------------------------------------------------------
    // event loop
    //----------------------------------------------------------------
    long n_pass = 0;
    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );

        if ( i % ( nentries / 10 ) == 0  && i > 0 )
        {
            std::cout << "Processing entry " << i << " / " << nentries << std::endl;
        }
        
        if ( has_event_sel_inputs )
        {
            if ( !is_minbias ) continue;
            if ( std::abs( zvrtx ) > 60.0 ) continue;
            // if ( !scaled_triggervec[10] ) continue;
        }

        if ( has_truth_jets )
        {
            sort_truth_jets_by_pt();

            const size_t n_truth = truth_jet_pT -> size();
            truth_jet_accept_eta.assign( n_truth, 0 );
            for ( size_t it = 0; it < n_truth; ++it )
            {
                truth_jet_accept_eta[it] = AnaUtils::accept_jet_eta( truth_jet_eta -> at( it ), zvrtx, truth_jet_R ) ? 1 : 0;
            }
        }
        if ( has_rho_jets )
        {
            sort_reco_jets_by_pt();

            const size_t n_reco = rho_jet_pT -> size();
            jet_accept_eta.assign( n_reco, 0 );
            for ( size_t ir = 0; ir < n_reco; ++ir )
            {
                jet_accept_eta[ir] =
                    AnaUtils::accept_jet_eta( rho_jet_eta -> at( ir ), zvrtx, rho_jet_R ) ? 1 : 0;
            }
        }

        //--- derived outputs, computed only for events passing selection ---
        if ( has_sumeT_calc_inputs )
        {
            sumeT_cemc_calc  = AnaUtils::calc_sumeT( AnaUtils::CEMC, zvrtx, cemc_tower_E, cemc_tower_isgood );
            sumeT_ihcal_calc = AnaUtils::calc_sumeT( AnaUtils::HCALIN, zvrtx, ihcal_tower_E, ihcal_tower_isgood );
            sumeT_ohcal_calc = AnaUtils::calc_sumeT( AnaUtils::HCALOUT, zvrtx, ohcal_tower_E, ohcal_tower_isgood );
            sumeT_calc       = sumeT_cemc_calc + sumeT_ihcal_calc + sumeT_ohcal_calc;
        }

        mbd_q = mbd_q_N + mbd_q_S;

        if ( has_rho_jet_constituents )
        {
            const size_t n_jets = rho_jet_constituent_E -> size();
            rho_jet_cemcfrac.assign( n_jets, 0.0 );
            rho_jet_ihcalfrac.assign( n_jets, 0.0 );
            rho_jet_ohcalfrac.assign( n_jets, 0.0 );
            for ( size_t ij = 0; ij < n_jets; ++ij )
            {
                calc_calo_frac(
                    rho_jet_constituent_E -> at( ij ),
                    rho_jet_constituent_srcID -> at( ij ),
                    rho_jet_cemcfrac[ij], rho_jet_ihcalfrac[ij], rho_jet_ohcalfrac[ij]
                );
            }
        }

        if ( has_truth_reco_match_inputs )
        {
            // strict 1-to-1 matching: truth and reco jets are already
            // sorted descending by pT (see sort_truth_jets_by_pt/
            // sort_reco_jets_by_pt above), so walking truth in that order
            // and taking the first not-yet-used reco jet within dR is
            // equivalent to "highest-pT reco jet within range" while
            // guaranteeing no reco jet is claimed by more than one truth
            // jet.
            static const float max_dr = 0.75f * 0.3f;
            const size_t n_truth = truth_jet_pT -> size();
            const size_t n_reco  = rho_jet_pT   -> size();
            truth_jet_reco_match_idx.assign( n_truth, -1 );
            jet_truth_match_idx.assign( n_reco, -1 );

            std::vector< bool > reco_used( n_reco, false );
            for ( size_t it = 0; it < n_truth; ++it )
            {
                for ( size_t ir = 0; ir < n_reco; ++ir )
                {
                    if ( reco_used[ir] ) continue; // already matched

                    const float dr = AnaUtils::calc_dr(
                        truth_jet_eta -> at( it ), truth_jet_phi -> at( it ),
                        rho_jet_eta   -> at( ir ), rho_jet_phi   -> at( ir )
                    );
                    if ( dr > max_dr ) continue; // outside threshold

                    truth_jet_reco_match_idx[it] = static_cast< int >( ir );
                    jet_truth_match_idx[ir] = static_cast< int >( it );
                    reco_used[ir] = true;
                    break; // highest-pT priority: stop after first match
                }
            }
        }

        ++n_pass;
        tout -> Fill();
    }
    std::cout << "Events passing selection: " << n_pass << " / " << nentries << std::endl;

    fout -> cd();
    tout -> Write();
    const long n_written = tout -> GetEntries();
    fout -> Close();

    std::cout << "Wrote " << n_written << " entries to " << outfile << std::endl;

    return 0;
}

#endif
