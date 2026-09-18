
#ifndef _MATCH_DATA_C_
#define _MATCH_DATA_C_

#include <myana/AnaUtils.h>
#include <myana/RhoEtaCalibLookup.h>

#include <TBranch.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

//--------------------------------------------------------------------
// rho eta-shape calibration -- needed to "unsubtract" the stored towers
//
// The cemc/ihcal/ohcal tower arrays in the anatrees come from the MULTSUB
// nodes, i.e. they are already UE-subtracted. SubtractTowersRhov1 built
// them as
//
//     E_sub = E_raw - rho * w( layer, ieta ; zvtx, mbdQ ) * cosh( eta_corr )
//
// with eta_corr the z-vertex-corrected tower eta (AnaUtils::get_corrected_calo_eta,
// same radii/eta table as the tower geometry the module used) and w the
// per-run eta-shape calibration table that Fun4All_Dijets_AuAu.C handed it
// via set_etaCalib_directPath(), rho_calib/calibs/rho_calib_<run>.root in
// CDBTTree format. rho itself is stored per event in the same tree
// (rho_val_TowerRho_MULT_<layer>), so reading the table back recovers
// E_raw exactly -- which is what the sum eT below is computed from.
//--------------------------------------------------------------------
struct RhoEtaCalib
{
    static constexpr int k_n_eta = 24;

    bool loaded { false };
    int n_zvtx { 0 };
    int n_mbdQ { 0 };
    std::vector< float > zvtx_edges {};
    std::vector< float > mbdQ_edges {};
    // one entry per encoded channel, layer 0/1/2 = cemc/ihcal/ohcal
    std::vector< std::array< float, 3 > > w {};

    // same binning convention as SubtractTowersRhov1::find_bin
    static int find_bin( const float val, const std::vector< float > & edges )
    {
        if ( edges.size() < 2 || std::isnan( val ) || val < edges.front() || val >= edges.back() )
        {
            return -1;
        }
        for ( size_t i = 0; i + 1 < edges.size(); ++i )
        {
            if ( val >= edges[i] && val < edges[i + 1] ) return static_cast< int >( i );
        }
        return -1;
    }

    // SubtractTowersRhov1::encode_channel
    int encode_channel( const int ieta, const int izbin, const int imbd ) const
    {
        return izbin * ( n_mbdQ * k_n_eta ) + imbd * k_n_eta + ieta;
    }

    // SubtractTowersRhov1::get_etaWeight -- 1.0 whenever the module itself
    // would have fallen back to a flat rho (no calibration, event outside
    // the calibrated (zvtx,mbdQ) range, or a non-positive stored weight)
    float get_weight( const int layer, const int ieta, const float zvtx, const float mbdQ ) const
    {
        if ( !loaded ) return 1.0;
        if ( layer < 0 || layer > 2 ) return 1.0;
        if ( ieta < 0 || ieta >= k_n_eta ) return 1.0;

        const int izbin = find_bin( zvtx, zvtx_edges );
        const int imbd  = find_bin( mbdQ, mbdQ_edges );
        if ( izbin < 0 || imbd < 0 ) return 1.0;

        const int channel = encode_channel( ieta, izbin, imbd );
        if ( channel < 0 || channel >= static_cast< int >( w.size() ) ) return 1.0;

        const float val = w[ channel ][ layer ];
        if ( !( val > 0 ) || std::isnan( val ) ) return 1.0;
        return val;
    }
};

// reads a CDBTTree-format calibration file by hand (Single tree = scalars
// with a type-letter prefix, Multiple tree = per-channel payload keyed by
// IID) so the macro needs nothing beyond ROOT.
RhoEtaCalib load_rho_eta_calib( const std::string & path )
{
    RhoEtaCalib calib;

    auto * f = TFile::Open( path.c_str(), "READ" );
    if ( !f || f -> IsZombie() )
    {
        std::cerr << "Warning: could not open rho calibration file " << path << std::endl;
        delete f;
        return calib;
    }

    auto * single = dynamic_cast< TTree * >( f -> Get( "Single" ) );
    auto * multi  = dynamic_cast< TTree * >( f -> Get( "Multiple" ) );
    if ( !single || !multi || single -> GetEntries() < 1 )
    {
        std::cerr << "Warning: rho calibration file " << path << " has no Single/Multiple trees" << std::endl;
        f -> Close();
        delete f;
        return calib;
    }

    int n_eta = 0, n_z = 0, n_mbd = 0;
    single -> SetBranchAddress( "In_eta", &n_eta );
    single -> SetBranchAddress( "In_zvtx_bins", &n_z );
    single -> SetBranchAddress( "In_mbdQ_bins", &n_mbd );
    single -> GetEntry( 0 );

    if ( n_eta != RhoEtaCalib::k_n_eta || n_z <= 0 || n_mbd <= 0 )
    {
        std::cerr << "Warning: rho calibration file " << path << " looks invalid (n_eta=" << n_eta
                  << ", n_zvtx_bins=" << n_z << ", n_mbdQ_bins=" << n_mbd << ")" << std::endl;
        f -> Close();
        delete f;
        return calib;
    }

    calib.n_zvtx = n_z;
    calib.n_mbdQ = n_mbd;
    calib.zvtx_edges.assign( n_z + 1, 0.0 );
    calib.mbdQ_edges.assign( n_mbd + 1, 0.0 );
    for ( int i = 0; i <= n_z; ++i )
    {
        single -> SetBranchAddress( Form( "Fzvtx_edge_%d", i ), &calib.zvtx_edges[i] );
    }
    for ( int i = 0; i <= n_mbd; ++i )
    {
        single -> SetBranchAddress( Form( "FmbdQ_edge_%d", i ), &calib.mbdQ_edges[i] );
    }
    single -> GetEntry( 0 );

    int id = -1;
    float w_cemc = 1.0, w_hcalin = 1.0, w_hcalout = 1.0;
    multi -> SetBranchAddress( "IID", &id );
    multi -> SetBranchAddress( "Fw_cemc", &w_cemc );
    multi -> SetBranchAddress( "Fw_hcalin", &w_hcalin );
    multi -> SetBranchAddress( "Fw_hcalout", &w_hcalout );

    calib.w.assign( n_z * n_mbd * n_eta, std::array< float, 3 >{ { 1.0f, 1.0f, 1.0f } } );
    for ( long i = 0; i < multi -> GetEntries(); ++i )
    {
        multi -> GetEntry( i );
        if ( id < 0 || id >= static_cast< int >( calib.w.size() ) ) continue;
        calib.w[ id ] = std::array< float, 3 >{ { w_cemc, w_hcalin, w_hcalout } };
    }

    f -> Close();
    delete f;

    calib.loaded = true;
    return calib;
}

// anatree files are named <...>-<8-digit run>-<5-digit segment>.root; the
// data trees carry no run-number branch, so the per-run calibration has to
// be picked from the file name.
int run_number_from_filename( const std::string & path )
{
    std::string base = path;
    const size_t slash = base.find_last_of( '/' );
    if ( slash != std::string::npos ) base = base.substr( slash + 1 );
    const size_t dot = base.rfind( ".root" );
    if ( dot != std::string::npos ) base = base.substr( 0, dot );

    std::vector< std::string > tokens;
    size_t start = 0;
    while ( start <= base.size() )
    {
        const size_t dash = base.find( '-', start );
        tokens.push_back( base.substr( start, dash - start ) );
        if ( dash == std::string::npos ) break;
        start = dash + 1;
    }
    if ( tokens.size() < 2 ) return -1;

    const std::string & run_tok = tokens[ tokens.size() - 2 ];
    if ( run_tok.empty() ) return -1;
    for ( char c : run_tok )
    {
        if ( !std::isdigit( static_cast< unsigned char >( c ) ) ) return -1;
    }
    return std::atoi( run_tok.c_str() );
}

// Data-only counterpart of truth-matching/match_standalone.C: same event
// selection, sumeT/em-fraction derivations and rho-jet pT sort, but with
// no truth-jet / HepMC-derived branches (data anatree files never carry
// them) and no truth-to-reco matching.
//
// Jets are written out trimmed to pT > jet_pt_min (default 5 GeV). The
// containers hold ~90 jets/event, nearly all of them sub-GeV, and writing
// all of them dominated the output file size.
int match_data(
    const std::string & infile = "/sphenix/user/tmengel/dijet-ana-auau/catalog/data_v004_20260821/data_v004_20260821-000.list",
    const std::string & outfile = "output.root",
    const float jet_pt_min = 3.0 // lower for prob file generation
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
    }

    float ihcal_tower_E[k_ieta][k_iphi] {};
    int   ihcal_tower_isgood[k_ieta][k_iphi] {};
    if ( has_branch( "ihcal_tower_E" ) )
    {
        t -> SetBranchStatus ( "ihcal_tower_E", true );
        t -> SetBranchStatus ( "ihcal_tower_isgood", true );
        t -> SetBranchAddress( "ihcal_tower_E", ihcal_tower_E );
        t -> SetBranchAddress( "ihcal_tower_isgood", ihcal_tower_isgood );
    }

    float ohcal_tower_E[k_ieta][k_iphi] {};
    int   ohcal_tower_isgood[k_ieta][k_iphi] {};
    if ( has_branch( "ohcal_tower_E" ) )
    {
        t -> SetBranchStatus ( "ohcal_tower_E", true );
        t -> SetBranchStatus ( "ohcal_tower_isgood", true );
        t -> SetBranchAddress( "ohcal_tower_E", ohcal_tower_E );
        t -> SetBranchAddress( "ohcal_tower_isgood", ohcal_tower_isgood );
    }

    //----------------------------------------------------------------
    // event-level sum eT, computed from the tower arrays the same way
    // match.C does (AnaUtils::calc_sumeT, zvrtx-corrected).
    //
    // The stored towers are the MULTSUB (rho-subtracted) ones, so they are
    // "unsubtracted" first -- see unsubtract_towers() below -- and the sum
    // eT written out is the raw, unsubtracted one.
    //----------------------------------------------------------------
    float sumeT_cemc_calc  = 0.0;
    float sumeT_ihcal_calc = 0.0;
    float sumeT_ohcal_calc = 0.0;
    float sumeT_calc       = 0.0;
    float cemc_tower_E_unsub[k_ieta][k_iphi] {};
    float ihcal_tower_E_unsub[k_ieta][k_iphi] {};
    float ohcal_tower_E_unsub[k_ieta][k_iphi] {};
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
    // written out per jet, holding only the jets above jet_pt_min -- the
    // input vectors above carry every jet in the container (~90/event, mostly
    // sub-GeV), which is what made the output files unmanageable.
    std::vector< float > jet_E_out {};
    std::vector< float > jet_phi_out {};
    std::vector< float > jet_eta_out {};
    std::vector< float > jet_pT_out {};
    std::vector< float > jet_unsub_pT_out {};
    std::vector< float > jet_unsub_E_out {};
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
        tout -> Branch( "jet_E", &jet_E_out );
        tout -> Branch( "jet_phi", &jet_phi_out );
        tout -> Branch( "jet_eta", &jet_eta_out );
        tout -> Branch( "jet_pT", &jet_pT_out );
        tout -> Branch( "jet_unsub_pT", &jet_unsub_pT_out );
        tout -> Branch( "jet_unsub_E", &jet_unsub_E_out );
        tout -> Branch( "jet_accept_eta", &jet_accept_eta );
    }

    // constituent_E/constituent_srcID feed calc_calo_frac; constituent_E/phi/eta
    // feed the per-jet constituent four-vector sum. constituent_pT is not read --
    // it is the *unsubtracted* tower pT, while the sum below is built from the
    // subtracted constituent energies.
    std::vector< std::vector< float > > * rho_jet_constituent_E   = nullptr;
    std::vector< std::vector< float > > * rho_jet_constituent_phi = nullptr;
    std::vector< std::vector< float > > * rho_jet_constituent_eta = nullptr;
    std::vector< std::vector< int > >   * rho_jet_constituent_srcID = nullptr;
    std::vector< float > rho_jet_cemcfrac {};
    std::vector< float > rho_jet_ihcalfrac{};
    std::vector< float > rho_jet_ohcalfrac{};
    std::vector< float > rho_jet_comp_E {};
    std::vector< float > rho_jet_comp_pT {};

    auto calc_const_fourvec = [&]( const std::vector< float > & E,
                                   const std::vector< float > & phi,
                                   const std::vector< float > & eta ) -> TLorentzVector
    {
        const size_t n = std::min( E.size(), std::min( phi.size(), eta.size() ) );
        double px_tot = 0.0, py_tot = 0.0, pz_tot = 0.0, e_tot = 0.0;
        for ( size_t i = 0; i < n; ++i )
        {
            float e = E[i];
            float pT = e/std::cosh( eta[i] );
            float px = pT * std::cos( phi[i] );
            float py = pT * std::sin( phi[i] );
            float pz = pT * std::sinh( eta[i] );

            px_tot += px;
            py_tot += py;
            pz_tot += pz;
            e_tot += e;
        }
        
        return TLorentzVector( px_tot, py_tot, pz_tot, e_tot );
    };

    const bool has_rho_jet_constituents = has_branch( "rho_jet_constituent_E" )
                                         && has_branch( "rho_jet_constituent_phi" )
                                         && has_branch( "rho_jet_constituent_eta" )
                                         && has_branch( "rho_jet_constituent_srcID" );
    if ( has_rho_jet_constituents )
    {
        t -> SetBranchStatus ( "rho_jet_constituent_E", true );
        t -> SetBranchStatus ( "rho_jet_constituent_phi", true );
        t -> SetBranchStatus ( "rho_jet_constituent_eta", true );
        t -> SetBranchStatus ( "rho_jet_constituent_srcID", true );

        t -> SetBranchAddress( "rho_jet_constituent_E", &rho_jet_constituent_E );
        t -> SetBranchAddress( "rho_jet_constituent_phi", &rho_jet_constituent_phi );
        t -> SetBranchAddress( "rho_jet_constituent_eta", &rho_jet_constituent_eta );
        t -> SetBranchAddress( "rho_jet_constituent_srcID", &rho_jet_constituent_srcID );

        tout -> Branch( "jet_cemcfrac", &rho_jet_cemcfrac );
        tout -> Branch( "jet_ihcalfrac", &rho_jet_ihcalfrac );
        tout -> Branch( "jet_ohcalfrac", &rho_jet_ohcalfrac );
        tout -> Branch( "jet_comp_E", &rho_jet_comp_E );
        tout -> Branch( "jet_comp_pT", &rho_jet_comp_pT );
    }

    // re-sorts the reco (rho-subtracted) jet vectors in place, descending
    // by pT -- including the constituent lists, so the em-fractions stay
    // aligned with the new order. Must run before the em-fraction calc.
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
            apply_vf( *rho_jet_constituent_phi );
            apply_vf( *rho_jet_constituent_eta );
            apply_vi( *rho_jet_constituent_srcID );
        }
    };

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
    // tower "unsubtraction" -- undoes the SubtractTowersRhov1 subtraction
    // that produced the stored MULTSUB towers, so the sum eT above is the
    // raw one. Needs the per-layer rho of the event plus the same per-run
    // eta-shape calibration Fun4All_Dijets_AuAu.C fed the module.
    //----------------------------------------------------------------
    auto rho_index = [&]( const std::string & node ) -> int
    {
        for ( size_t i = 0; i < rho_names.size(); ++i )
        {
            if ( rho_names[i] == node ) return static_cast< int >( i );
        }
        return -1;
    };
    const int irho_cemc  = rho_index( "TowerRho_MULT_CEMC" );
    const int irho_ihcal = rho_index( "TowerRho_MULT_HCALIN" );
    const int irho_ohcal = rho_index( "TowerRho_MULT_HCALOUT" );

    const bool has_unsub_inputs = has_sumeT_calc_inputs
                                 && has_branch( "mbd_q_N" )
                                 && irho_cemc >= 0 && irho_ihcal >= 0 && irho_ohcal >= 0;
    if ( has_sumeT_calc_inputs && !has_unsub_inputs )
    {
        std::cerr << "Warning: rho / mbd_q branches missing -- sum eT will be computed "
                  << "from the subtracted towers as-is." << std::endl;
    }

    // same lookup Fun4All_Dijets_AuAu.C uses: the run's own calibration, else the
    // run2auau default (calibrations/rho_eta/README.md). A file whose run number
    // cannot be parsed gets the default as well.
    std::map< int, RhoEtaCalib > calib_cache;
    auto get_calib = [&]( const int run ) -> const RhoEtaCalib &
    {
        auto it = calib_cache.find( run );
        if ( it != calib_cache.end() ) return it -> second;

        const std::string path = RhoEtaCalibLookup::GetCalibPath(
            run < 0 ? RhoEtaCalibLookup::kRun2AuAuDefaultRun : run );
        return calib_cache.emplace( run, load_rho_eta_calib( path ) ).first -> second;
    };

    // E_raw = E_sub + rho * w * cosh( eta_corr ), for the towers flagged good
    // (bad towers were zeroed by the module and are skipped by calc_sumeT).
    auto unsubtract_towers = [&](
        const AnaUtils::CaloType calo, const int layer, const float rho,
        const RhoEtaCalib & calib, const float zv, const float mbdQ,
        const float E_sub[k_ieta][k_iphi], const int isgood[k_ieta][k_iphi],
        float E_unsub[k_ieta][k_iphi]
    )
    {
        for ( int ieta = 0; ieta < k_ieta; ++ieta )
        {
            const float cosh_eta = std::cosh( AnaUtils::get_corrected_calo_eta( calo, ieta, zv ) );
            const float ue = rho * cosh_eta * calib.get_weight( layer, ieta, zv, mbdQ );
            for ( int iphi = 0; iphi < k_iphi; ++iphi )
            {
                E_unsub[ieta][iphi] = isgood[ieta][iphi] ? E_sub[ieta][iphi] + ue : 0.0f;
            }
        }
    };

    // the chain spans runs, so the calibration follows the current file
    int cur_tree_number = -1;
    const RhoEtaCalib * cur_calib = nullptr;

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

        // indices, into the sorted input vectors, of the jets actually written
        std::vector< size_t > sel_idx {};
        if ( has_rho_jets )
        {
            sort_reco_jets_by_pt();

            // descending pT after the sort, so the jets above threshold are a
            // leading prefix and keep their order -- the same convention
            // calibrate_matched.C's own pT>5 trim relies on.
            const size_t n_reco = rho_jet_pT -> size();
            for ( size_t ir = 0; ir < n_reco; ++ir )
            {
                if ( rho_jet_pT -> at( ir ) <= jet_pt_min ) break;
                sel_idx.push_back( ir );
            }

            const size_t n_sel = sel_idx.size();
            jet_E_out.resize( n_sel );
            jet_phi_out.resize( n_sel );
            jet_eta_out.resize( n_sel );
            jet_pT_out.resize( n_sel );
            jet_unsub_pT_out.resize( n_sel );
            jet_unsub_E_out.resize( n_sel );
            jet_accept_eta.assign( n_sel, 0 );
            for ( size_t k = 0; k < n_sel; ++k )
            {
                const size_t ir = sel_idx[k];
                jet_E_out[k]        = rho_jet_E        -> at( ir );
                jet_phi_out[k]      = rho_jet_phi      -> at( ir );
                jet_eta_out[k]      = rho_jet_eta      -> at( ir );
                jet_pT_out[k]       = rho_jet_pT       -> at( ir );
                jet_unsub_pT_out[k] = rho_jet_unsub_pT -> at( ir );
                jet_unsub_E_out[k]  = rho_jet_unsub_E  -> at( ir );
                jet_accept_eta[k]   =
                    AnaUtils::accept_jet_eta( rho_jet_eta -> at( ir ), zvrtx, rho_jet_R ) ? 1 : 0;
            }
        }
        else if ( has_rho_jet_constituents )
        {
            // no rho_jet_pT to threshold on -- keep every jet
            sel_idx.resize( rho_jet_constituent_E -> size() );
            std::iota( sel_idx.begin(), sel_idx.end(), 0 );
        }

        mbd_q = mbd_q_N + mbd_q_S;

        //--- derived outputs, computed only for events passing selection ---
        if ( has_sumeT_calc_inputs )
        {
            if ( has_unsub_inputs )
            {
                if ( t -> GetTreeNumber() != cur_tree_number )
                {
                    cur_tree_number = t -> GetTreeNumber();
                    auto * cur_file = t -> GetFile();
                    cur_calib = &get_calib( cur_file ? run_number_from_filename( cur_file -> GetName() ) : -1 );
                }

                unsubtract_towers( AnaUtils::CEMC, 0, rho_vals[irho_cemc], *cur_calib, zvrtx, mbd_q,
                                   cemc_tower_E, cemc_tower_isgood, cemc_tower_E_unsub );
                unsubtract_towers( AnaUtils::HCALIN, 1, rho_vals[irho_ihcal], *cur_calib, zvrtx, mbd_q,
                                   ihcal_tower_E, ihcal_tower_isgood, ihcal_tower_E_unsub );
                unsubtract_towers( AnaUtils::HCALOUT, 2, rho_vals[irho_ohcal], *cur_calib, zvrtx, mbd_q,
                                   ohcal_tower_E, ohcal_tower_isgood, ohcal_tower_E_unsub );

                sumeT_cemc_calc  = AnaUtils::calc_sumeT( AnaUtils::CEMC, zvrtx, cemc_tower_E_unsub, cemc_tower_isgood );
                sumeT_ihcal_calc = AnaUtils::calc_sumeT( AnaUtils::HCALIN, zvrtx, ihcal_tower_E_unsub, ihcal_tower_isgood );
                sumeT_ohcal_calc = AnaUtils::calc_sumeT( AnaUtils::HCALOUT, zvrtx, ohcal_tower_E_unsub, ohcal_tower_isgood );
            }
            else
            {
                sumeT_cemc_calc  = AnaUtils::calc_sumeT( AnaUtils::CEMC, zvrtx, cemc_tower_E, cemc_tower_isgood );
                sumeT_ihcal_calc = AnaUtils::calc_sumeT( AnaUtils::HCALIN, zvrtx, ihcal_tower_E, ihcal_tower_isgood );
                sumeT_ohcal_calc = AnaUtils::calc_sumeT( AnaUtils::HCALOUT, zvrtx, ohcal_tower_E, ohcal_tower_isgood );
            }
            sumeT_calc = sumeT_cemc_calc + sumeT_ihcal_calc + sumeT_ohcal_calc;
        }

        if ( has_rho_jet_constituents )
        {
            const size_t n_sel = sel_idx.size();
            rho_jet_cemcfrac.assign( n_sel, 0.0 );
            rho_jet_ihcalfrac.assign( n_sel, 0.0 );
            rho_jet_ohcalfrac.assign( n_sel, 0.0 );
            rho_jet_comp_E.assign( n_sel, 0.0 );
            rho_jet_comp_pT.assign( n_sel, 0.0 );
            for ( size_t k = 0; k < n_sel; ++k )
            {
                const size_t ij = sel_idx[k];
                calc_calo_frac(
                    rho_jet_constituent_E -> at( ij ),
                    rho_jet_constituent_srcID -> at( ij ),
                    rho_jet_cemcfrac[k], rho_jet_ihcalfrac[k], rho_jet_ohcalfrac[k]
                );

                const TLorentzVector p4 = calc_const_fourvec(
                    rho_jet_constituent_E -> at( ij ),
                    rho_jet_constituent_phi -> at( ij ),
                    rho_jet_constituent_eta -> at( ij )
                );

                rho_jet_comp_E[k]  = p4.E();
                rho_jet_comp_pT[k] = p4.Pt();
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
