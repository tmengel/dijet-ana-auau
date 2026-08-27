
#ifndef _FIND_REFERENCE_CALIB_C_
#define _FIND_REFERENCE_CALIB_C_

#include <cdbobjects/CDBTTree.h>

#include <TFile.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libcdbobjects.so )

// Compares every run-by-run rho eta-shape calibration in calibs/ and picks
// the one whose weight vector is closest to the population mean -- i.e. the
// most "typical" run -- as a candidate default/fallback calibration for runs
// that don't (yet) have their own. All calibration files share the same
// binning (defaults of make_rho_calib.C: n_eta=24, n_zvtx_bins=7,
// n_mbdQ_bins=10 -> 1680 channels x 3 layers = 5040 values per run), so
// per-channel comparison across runs is well-defined without rebinning.
//
// Runs with a large fraction of empty (fallback w=1) channels are excluded
// from being a *candidate* reference (their vector is partly a placeholder,
// not a measurement) but are still included when computing the population
// mean, so they don't bias the target -- they're just not eligible to *be*
// the answer.
int find_reference_calib(
    const std::string & calib_dir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/calibs",
    const double max_fallback_frac = 0.05 // exclude runs with >5% empty channels from candidacy
)
{
    static const int kNCalo = 3;
    static const char * calib_field_names[kNCalo] = { "w_cemc", "w_hcalin", "w_hcalout" };

    // -- discover runs from calibs/rho_calib_<run>.root (skip _qa companions) --
    std::vector<int> runs;
    {
        TSystemDirectory dir( "calibdir", calib_dir.c_str() );
        TList * files = dir.GetListOfFiles();
        if ( files )
        {
            for ( auto * obj : *files )
            {
                TString name = obj->GetName();
                if ( !name.BeginsWith( "rho_calib_" ) || !name.EndsWith( ".root" ) || name.Contains( "_qa" ) ) continue;
                TString numpart = name;
                numpart.ReplaceAll( "rho_calib_", "" );
                numpart.ReplaceAll( ".root", "" );
                if ( numpart.IsDigit() ) runs.push_back( numpart.Atoi() );
            }
        }
    }
    std::sort( runs.begin(), runs.end() );
    std::cout << "Found " << runs.size() << " run calibrations in " << calib_dir << std::endl;
    if ( runs.empty() ) return -1;

    // -- load every run's full weight vector --
    int NETA = -1, n_zbins = -1, n_mbdbins = -1, nchannels = -1;
    std::vector<std::vector<float>> weights; // [run_idx][channel*3 + layer]
    std::vector<int> n_fallback( runs.size(), 0 );
    weights.reserve( runs.size() );

    for ( size_t ir = 0; ir < runs.size(); ++ir )
    {
        const int run = runs[ir];
        const std::string calibfile = Form( "%s/rho_calib_%d.root", calib_dir.c_str(), run );
        auto * cdbttree = new CDBTTree( calibfile );
        cdbttree->LoadCalibrations();

        const int this_neta = cdbttree->GetSingleIntValue( "n_eta" );
        const int this_nz = cdbttree->GetSingleIntValue( "n_zvtx_bins" );
        const int this_nmbd = cdbttree->GetSingleIntValue( "n_mbdQ_bins" );
        if ( NETA < 0 ) { NETA = this_neta; n_zbins = this_nz; n_mbdbins = this_nmbd; nchannels = n_zbins * n_mbdbins * NETA; }
        else if ( this_neta != NETA || this_nz != n_zbins || this_nmbd != n_mbdbins )
        {
            std::cerr << "Warning: run " << run << " has different binning (n_eta=" << this_neta
                      << ", n_zvtx_bins=" << this_nz << ", n_mbdQ_bins=" << this_nmbd
                      << ") -- skipping (not directly comparable)." << std::endl;
            delete cdbttree;
            runs.erase( runs.begin() + ir );
            n_fallback.erase( n_fallback.begin() + ir );
            --ir;
            continue;
        }

        std::vector<float> v( nchannels * kNCalo, 1.0f );
        int nfb = 0;
        for ( int ch = 0; ch < nchannels; ++ch )
        {
            for ( int ic = 0; ic < kNCalo; ++ic )
            {
                const float w = cdbttree->GetFloatValue( ch, calib_field_names[ic] );
                v[ch * kNCalo + ic] = w;
                if ( !( w > 0.0f ) || w == 1.0f ) ++nfb; // fallback/empty-channel heuristic
            }
        }
        weights.push_back( std::move( v ) );
        n_fallback[ir] = nfb;
        delete cdbttree;
    }

    const size_t nruns = runs.size();
    const int ntotal_per_run = nchannels * kNCalo;
    std::cout << "nchannels/run = " << nchannels << " (n_eta=" << NETA << ", n_zvtx_bins=" << n_zbins
              << ", n_mbdQ_bins=" << n_mbdbins << "), " << ntotal_per_run << " values/run (x3 layers)" << std::endl;

    // -- population mean vector across ALL runs (unweighted; fallback runs included so they don't bias candidacy but do count toward "typical") --
    std::vector<double> mean( ntotal_per_run, 0.0 );
    for ( size_t ir = 0; ir < nruns; ++ir )
        for ( int k = 0; k < ntotal_per_run; ++k )
            mean[k] += weights[ir][k];
    for ( int k = 0; k < ntotal_per_run; ++k ) mean[k] /= nruns;

    // -- per-run RMS distance from the population mean --
    std::vector<double> dist( nruns, 0.0 );
    for ( size_t ir = 0; ir < nruns; ++ir )
    {
        double sumsq = 0.0;
        for ( int k = 0; k < ntotal_per_run; ++k )
        {
            const double d = weights[ir][k] - mean[k];
            sumsq += d * d;
        }
        dist[ir] = std::sqrt( sumsq / ntotal_per_run );
    }

    // -- rank candidates (exclude high-fallback runs from being the answer) --
    std::vector<size_t> order( nruns );
    for ( size_t i = 0; i < nruns; ++i ) order[i] = i;
    std::sort( order.begin(), order.end(), [&]( size_t a, size_t b ) { return dist[a] < dist[b]; } );

    std::cout << "\n=== Top 10 most-typical runs (RMS distance from population mean, all runs) ===" << std::endl;
    for ( size_t i = 0; i < std::min<size_t>( 10, nruns ); ++i )
    {
        const size_t ir = order[i];
        const double fb_frac = double( n_fallback[ir] ) / ntotal_per_run;
        std::cout << "  run " << runs[ir] << "  RMS_dist=" << dist[ir]
                  << "  fallback_frac=" << 100.0 * fb_frac << "%" << std::endl;
    }

    std::cout << "\n=== Best candidate (excluding runs with >" << 100.0 * max_fallback_frac << "% fallback channels) ===" << std::endl;
    int best_run = -1; double best_dist = 1e18;
    for ( size_t i = 0; i < nruns; ++i )
    {
        const size_t ir = order[i];
        const double fb_frac = double( n_fallback[ir] ) / ntotal_per_run;
        if ( fb_frac > max_fallback_frac ) continue;
        best_run = runs[ir];
        best_dist = dist[ir];
        break;
    }
    if ( best_run > 0 )
    {
        std::cout << "  -> run " << best_run << "  (RMS_dist=" << best_dist << ")" << std::endl;
        std::cout << "  calibs/rho_calib_" << best_run << ".root" << std::endl;
    }
    else
    {
        std::cout << "  No run passed the fallback-fraction cut." << std::endl;
    }

    // -- also report the worst outliers, for context --
    std::cout << "\n=== 5 least-typical runs (largest RMS distance) ===" << std::endl;
    for ( size_t i = 0; i < 5 && i < nruns; ++i )
    {
        const size_t ir = order[nruns - 1 - i];
        const double fb_frac = double( n_fallback[ir] ) / ntotal_per_run;
        std::cout << "  run " << runs[ir] << "  RMS_dist=" << dist[ir]
                  << "  fallback_frac=" << 100.0 * fb_frac << "%" << std::endl;
    }

    // -- explicit lookup for any specific run of interest (rank among ALL nruns, 1 = most typical) --
    auto report_run = [&]( int run_of_interest ) {
        for ( size_t i = 0; i < nruns; ++i )
        {
            if ( runs[order[i]] != run_of_interest ) continue;
            const size_t ir = order[i];
            const double fb_frac = double( n_fallback[ir] ) / ntotal_per_run;
            std::cout << "\nrun " << run_of_interest << ":  rank " << ( i + 1 ) << " / " << nruns
                      << "  RMS_dist=" << dist[ir] << "  fallback_frac=" << 100.0 * fb_frac << "%" << std::endl;
            return;
        }
        std::cout << "\nrun " << run_of_interest << ": not found among processed calibrations" << std::endl;
    };
    report_run( 54912 );
    report_run( 54590 );

    return 0;
}

#endif
