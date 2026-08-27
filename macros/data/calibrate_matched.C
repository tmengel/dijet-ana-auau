
#ifndef _CALIBRATE_MATCHED_C_
#define _CALIBRATE_MATCHED_C_

#include <myana/AnaUtils.h>

#include <TFile.h>
#include <TProfile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

// Run-by-run (file-by-file) sumeT calibration on top of the merged
// match_data.C output (one of the 119 data_v004_20260821-NNN_matched.root
// chunks). Each chunk is its own calibration unit here -- no trigger info
// survives the merge (match_data.C never wrote scaled_triggervec), so the
// event selection is is_minbias && |zvrtx|<60 only.
//   1st pass -- fill a TProfile of sumeT vs. centrality (2% bins, "s"
//               error option so bin errors are the bin's standard
//               deviation rather than the error on the mean).
//   2nd pass -- reject events whose sumeT falls outside
//               <sumeT>(cent) +/- 3.5*sigma(cent) using this chunk's own
//               profile. Survivors get their jet collection (already
//               pT-sorted by match_data.C) trimmed to pT>5 GeV, which
//               preserves the existing descending-pT order.
int calibrate_matched(
    const std::string & infile = "",
    const std::string & outfile = "output.root"
)
{
    auto * fin = TFile::Open( infile.c_str() );
    if ( !fin || fin -> IsZombie() )
    {
        std::cerr << "Error: could not open input file " << infile << std::endl;
        return -1;
    }
    auto * t = static_cast< TTree * >( fin -> Get( "T" ) );
    if ( !t )
    {
        std::cerr << "Error: no TTree \"T\" in " << infile << std::endl;
        return -1;
    }

    int nentries = t -> GetEntries();
    if ( nentries < 1 )
    {
        std::cerr << "Error: No entries in " << infile << std::endl;
        return -1;
    }
    std::cout << "Total entries: " << nentries << std::endl;

    t -> SetBranchStatus( "*", false );

    auto has_branch = [ t ]( const std::string & name )
    {
        return t -> GetBranch( name.c_str() ) != nullptr;
    };

    //----------------------------------------------------------------
    // event selection inputs: minbias + |zvrtx|<60 (no trigger info
    // survives the merge)
    //----------------------------------------------------------------
    int is_minbias = 0;
    float zvrtx = 0.0;
    if ( !has_branch( "is_minbias" ) || !has_branch( "zvrtx" ) )
    {
        std::cerr << "Error: is_minbias/zvrtx not present in " << infile << std::endl;
        fin -> Close();
        return -1;
    }
    t -> SetBranchStatus ( "is_minbias", true );
    t -> SetBranchStatus ( "zvrtx", true );
    t -> SetBranchAddress( "is_minbias", &is_minbias );
    t -> SetBranchAddress( "zvrtx", &zvrtx );

    auto event_selected = [&]()
    {
        return is_minbias && std::abs( zvrtx ) < 60.0;
    };

    //----------------------------------------------------------------
    // pass-through event-level branches
    //----------------------------------------------------------------
    int event_id = -1;
    if ( has_branch( "event_id" ) )
    {
        t -> SetBranchStatus ( "event_id", true );
        t -> SetBranchAddress( "event_id", &event_id );
    }

    int cent = -1;
    const bool has_cent = has_branch( "cent" );
    if ( has_cent )
    {
        t -> SetBranchStatus ( "cent", true );
        t -> SetBranchAddress( "cent", &cent );
    }

    float mbd_q = -999.0;
    const bool has_mbd = has_branch( "mbd_q" );
    if ( has_mbd )
    {
        t -> SetBranchStatus ( "mbd_q", true );
        t -> SetBranchAddress( "mbd_q", &mbd_q );
    }

    float sumeT_cemc = 0.0, sumeT_ihcal = 0.0, sumeT_ohcal = 0.0, sumeT = 0.0;
    const bool has_sumeT = has_branch( "sumeT" );
    if ( !has_sumeT )
    {
        std::cerr << "Error: sumeT not present in " << infile << std::endl;
        fin -> Close();
        return -1;
    }
    t -> SetBranchStatus ( "sumeT_cemc", true );
    t -> SetBranchStatus ( "sumeT_ihcal", true );
    t -> SetBranchStatus ( "sumeT_ohcal", true );
    t -> SetBranchStatus ( "sumeT", true );
    t -> SetBranchAddress( "sumeT_cemc", &sumeT_cemc );
    t -> SetBranchAddress( "sumeT_ihcal", &sumeT_ihcal );
    t -> SetBranchAddress( "sumeT_ohcal", &sumeT_ohcal );
    t -> SetBranchAddress( "sumeT", &sumeT );

    //----------------------------------------------------------------
    // jet collection (already rho-subtracted and pT-sorted by match_data.C)
    //----------------------------------------------------------------
    float jet_R = 0.0;
    std::vector< float > * jet_E        = nullptr;
    std::vector< float > * jet_phi      = nullptr;
    std::vector< float > * jet_eta      = nullptr;
    std::vector< float > * jet_pT       = nullptr;
    std::vector< float > * jet_unsub_pT = nullptr;
    std::vector< float > * jet_unsub_E  = nullptr;
    std::vector< int >   * jet_accept_eta  = nullptr;
    std::vector< float > * jet_cemcfrac   = nullptr;
    std::vector< float > * jet_ihcalfrac  = nullptr;
    std::vector< float > * jet_ohcalfrac  = nullptr;

    const bool has_jets = has_branch( "jet_pT" );
    if ( has_jets )
    {
        t -> SetBranchStatus ( "jet_R", true );
        t -> SetBranchStatus ( "jet_E", true );
        t -> SetBranchStatus ( "jet_phi", true );
        t -> SetBranchStatus ( "jet_eta", true );
        t -> SetBranchStatus ( "jet_pT", true );
        t -> SetBranchStatus ( "jet_unsub_pT", true );
        t -> SetBranchStatus ( "jet_unsub_E", true );
        t -> SetBranchStatus ( "jet_accept_eta", true );

        t -> SetBranchAddress( "jet_R", &jet_R );
        t -> SetBranchAddress( "jet_E", &jet_E );
        t -> SetBranchAddress( "jet_phi", &jet_phi );
        t -> SetBranchAddress( "jet_eta", &jet_eta );
        t -> SetBranchAddress( "jet_pT", &jet_pT );
        t -> SetBranchAddress( "jet_unsub_pT", &jet_unsub_pT );
        t -> SetBranchAddress( "jet_unsub_E", &jet_unsub_E );
        t -> SetBranchAddress( "jet_accept_eta", &jet_accept_eta );
    }

    const bool has_jet_calo_frac = has_branch( "jet_cemcfrac" );
    if ( has_jet_calo_frac )
    {
        t -> SetBranchStatus ( "jet_cemcfrac", true );
        t -> SetBranchStatus ( "jet_ihcalfrac", true );
        t -> SetBranchStatus ( "jet_ohcalfrac", true );
        t -> SetBranchAddress( "jet_cemcfrac", &jet_cemcfrac );
        t -> SetBranchAddress( "jet_ihcalfrac", &jet_ihcalfrac );
        t -> SetBranchAddress( "jet_ohcalfrac", &jet_ohcalfrac );
    }

    //----------------------------------------------------------------
    // generic per-node rho values (rho_val_<node>, rho_sigma_<node>)
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
    }

    //==================================================================
    // pass 1: calibration -- sumeT vs. centrality profile, 2% bins
    //==================================================================
    static const int    k_cent_nbins   = 50;   // 2% wide bins over [0,100]
    static const double k_cent_lo      = 0.0;
    static const double k_cent_hi      = 100.0;
    static const int    k_min_bin_stat = 20;   // below this, skip the sumeT cut for that bin
    static const double k_nsigma       = 3.5;

    auto * hCal = new TProfile(
        "sumeT_vs_cent", "sumeT vs centrality;centrality [%];sum e_{T} [GeV]",
        k_cent_nbins, k_cent_lo, k_cent_hi, "s"
    );

    long n_cal_events = 0;
    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );
        if ( !event_selected() ) continue;
        if ( !has_cent ) continue;

        hCal -> Fill( cent, sumeT );
        ++n_cal_events;
    }
    std::cout << "Calibration pass: " << n_cal_events << " / " << nentries
               << " events passed event selection" << std::endl;

    std::vector< double > cent_mean( k_cent_nbins, 0.0 );
    std::vector< double > cent_sigma( k_cent_nbins, 0.0 );
    std::vector< bool >   cent_valid( k_cent_nbins, false );
    for ( int b = 0; b < k_cent_nbins; ++b )
    {
        if ( hCal -> GetBinEntries( b + 1 ) < k_min_bin_stat ) continue;
        cent_mean[b]  = hCal -> GetBinContent( b + 1 );
        cent_sigma[b] = hCal -> GetBinError( b + 1 ); // "s" option -> std dev, not error on mean
        cent_valid[b] = true;
    }

    auto cent_bin = [&]( const int c )
    {
        int b = static_cast< int >( c / ( k_cent_hi / k_cent_nbins ) );
        return std::clamp( b, 0, k_cent_nbins - 1 );
    };

    //==================================================================
    // pass 2: apply sumeT calibration cut + jet pT cut, write output
    //==================================================================
    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * tout = new TTree( "T", "T" );

    tout -> Branch( "event_id", &event_id, "event_id/I" );
    tout -> Branch( "is_minbias", &is_minbias, "is_minbias/I" );
    tout -> Branch( "zvrtx", &zvrtx, "zvrtx/F" );
    tout -> Branch( "cent", &cent, "cent/I" );
    if ( has_mbd ) tout -> Branch( "mbd_q", &mbd_q, "mbd_q/F" );
    tout -> Branch( "sumeT_cemc", &sumeT_cemc, "sumeT_cemc/F" );
    tout -> Branch( "sumeT_ihcal", &sumeT_ihcal, "sumeT_ihcal/F" );
    tout -> Branch( "sumeT_ohcal", &sumeT_ohcal, "sumeT_ohcal/F" );
    tout -> Branch( "sumeT", &sumeT, "sumeT/F" );

    std::vector< float > jet_E_out, jet_phi_out, jet_eta_out, jet_pT_out;
    std::vector< float > jet_unsub_pT_out, jet_unsub_E_out;
    std::vector< int >   jet_accept_eta_out;
    std::vector< float > jet_cemcfrac_out, jet_ihcalfrac_out, jet_ohcalfrac_out;
    if ( has_jets )
    {
        tout -> Branch( "jet_R", &jet_R, "jet_R/F" );
        tout -> Branch( "jet_E", &jet_E_out );
        tout -> Branch( "jet_phi", &jet_phi_out );
        tout -> Branch( "jet_eta", &jet_eta_out );
        tout -> Branch( "jet_pT", &jet_pT_out );
        tout -> Branch( "jet_unsub_pT", &jet_unsub_pT_out );
        tout -> Branch( "jet_unsub_E", &jet_unsub_E_out );
        tout -> Branch( "jet_accept_eta", &jet_accept_eta_out );
    }
    if ( has_jet_calo_frac )
    {
        tout -> Branch( "jet_cemcfrac", &jet_cemcfrac_out );
        tout -> Branch( "jet_ihcalfrac", &jet_ihcalfrac_out );
        tout -> Branch( "jet_ohcalfrac", &jet_ohcalfrac_out );
    }
    for ( size_t i = 0; i < rho_names.size(); ++i )
    {
        const std::string val_name   = "rho_val_" + rho_names[i];
        const std::string sigma_name = "rho_sigma_" + rho_names[i];
        tout -> Branch( val_name.c_str(), &rho_vals[i], ( val_name + "/F" ).c_str() );
        tout -> Branch( sigma_name.c_str(), &rho_sigmas[i], ( sigma_name + "/F" ).c_str() );
    }

    static const float k_min_jet_pt = 5.0;

    long n_sumeT_cut = 0;
    long n_written = 0;
    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );
        if ( !event_selected() ) continue;

        if ( has_cent )
        {
            const int b = cent_bin( cent );
            if ( cent_valid[b] && std::abs( sumeT - cent_mean[b] ) > k_nsigma * cent_sigma[b] )
            {
                ++n_sumeT_cut;
                continue;
            }
        }

        if ( has_jets )
        {
            // jet_pT is already sorted descending by match_data.C, so a
            // pT>5 filter that keeps relative order needs no re-sort.
            std::vector< int > idx;
            for ( size_t j = 0; j < jet_pT -> size(); ++j )
            {
                if ( jet_pT -> at( j ) > k_min_jet_pt ) idx.push_back( static_cast< int >( j ) );
            }

            const size_t n_sel = idx.size();
            jet_E_out.resize( n_sel );
            jet_phi_out.resize( n_sel );
            jet_eta_out.resize( n_sel );
            jet_pT_out.resize( n_sel );
            jet_unsub_pT_out.resize( n_sel );
            jet_unsub_E_out.resize( n_sel );
            jet_accept_eta_out.resize( n_sel );
            if ( has_jet_calo_frac )
            {
                jet_cemcfrac_out.resize( n_sel );
                jet_ihcalfrac_out.resize( n_sel );
                jet_ohcalfrac_out.resize( n_sel );
            }

            for ( size_t k = 0; k < n_sel; ++k )
            {
                const int j = idx[k];
                jet_E_out[k]          = jet_E -> at( j );
                jet_phi_out[k]        = jet_phi -> at( j );
                jet_eta_out[k]        = jet_eta -> at( j );
                jet_pT_out[k]         = jet_pT -> at( j );
                jet_unsub_pT_out[k]   = jet_unsub_pT -> at( j );
                jet_unsub_E_out[k]    = jet_unsub_E -> at( j );
                jet_accept_eta_out[k] = jet_accept_eta -> at( j );

                if ( has_jet_calo_frac )
                {
                    jet_cemcfrac_out[k]  = jet_cemcfrac -> at( j );
                    jet_ihcalfrac_out[k] = jet_ihcalfrac -> at( j );
                    jet_ohcalfrac_out[k] = jet_ohcalfrac -> at( j );
                }
            }
        }

        tout -> Fill();
        ++n_written;
    }
    std::cout << "sumeT calibration cut removed " << n_sumeT_cut << " events" << std::endl;
    std::cout << "Wrote " << n_written << " / " << nentries << " entries to " << outfile << std::endl;

    fout -> cd();
    tout -> Write();
    hCal -> Write();
    fout -> Close();
    fin  -> Close();

    return 0;
}

#endif
