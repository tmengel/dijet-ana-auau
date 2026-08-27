
#ifndef _MAKE_RHO_CALIB_SIM_C_
#define _MAKE_RHO_CALIB_SIM_C_

#include <myana/AnaUtils.h>

#include <cdbobjects/CDBTTree.h>

#include <TChain.h>
#include <TFile.h>
#include <TProfile.h>
#include <TString.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )
R__LOAD_LIBRARY( libcdbobjects.so )

// Same as make_rho_calib.C, but for the mb_hijing simulation samples (label
// = "mb_hijing_scaled" or "mb_hijing_unsclaed" instead of a numeric run
// number) -- reads the per-event tower statistics tree produced by
// read_tree_sim.C + hadd (one per sample, output/<label>/rho_calib_<label>.root)
// and derives the per-(layer,
// eta ring, zvertex bin, MBD-Q bin) rho calibration weight
//   w_i = <E_tower_i> / (rho * cosh(eta_corr_i))
// (PPG-14 Eq. 1). For every event, its (zvertex, mbdQ) locate a bin from
// zvtx_edges/mbdQ_edges (half-open, [lo,hi), events outside the edges are
// dropped); within that (izbin,imbd) bin, for every eta ring mean_E[ieta]
// (the average energy of the good towers in that ring) is divided by
// rho*cosh(eta_corr) for that event, and the ratio is accumulated into a
// TProfile vs ieta, weighted by ntowers (the number of good towers behind
// that event's mean_E). The calibration in a given (izbin,imbd) bin is
// therefore the ntowers-weighted mean of the ratio over the minbias events
// that fall in it. The centrality proxy is the raw MBD charge sum (mbdQ),
// not a calibrated centrality percentile, to avoid a dependency on a
// separate centrality calibration -- see SubtractTowersRhov1.
//
// Output is written directly in CDBTTree format, matching exactly what
// SubtractTowersRhov1::LoadEtaCalib / get_etaWeight expect:
//   -- single values: n_eta (=24), n_zvtx_bins, n_mbdQ_bins,
//      zvtx_edge_0..n_zvtx_bins, mbdQ_edge_0..n_mbdQ_bins
//   -- per-channel values: w_cemc, w_hcalin, w_hcalout, indexed by
//      channel = izbin*(n_mbdQ_bins*24) + imbd*24 + ieta
//      (SubtractTowersRhov1::encode_channel)
// A companion plain ROOT file (outfile with a "_qa" suffix) holding the
// TProfiles themselves (weight vs ieta, and <eta_corr> vs ieta) is also
// written per (izbin,imbd) bin, for plotting/debugging.
int make_rho_calib_sim(
    const std::string & label = "mb_hijing_scaled",
    const std::string & outdir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/calibs",
    // half-open bin edges, [lo,hi); events outside [edges.front(),edges.back())
    // are dropped from the calibration. Same defaults as the data pipeline
    // (make_rho_calib.C) for direct comparability between data and sim
    // calibrations -- re-tune if the sim zvrtx/mbdQ distributions turn out
    // to differ meaningfully from data's.
    std::vector<float> zvtx_edges = { -60.0, -30.0, -15.0, -5.0, 5.0, 15.0, 30.0, 60.0 },
    std::vector<float> mbdQ_edges = { 0.0, 10.0, 35.0, 85.0, 170.0, 300.0, 475.0, 700.0, 1000.0, 1400.0, 5000.0 }
)
{
    gSystem -> Exec( Form( "mkdir -p %s", outdir.c_str() ) );

    const std::string infile = Form( "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/output/%s/rho_calib_%s.root", label.c_str(), label.c_str() );
    const std::string outfile = Form( "%s/rho_calib_%s.root", outdir.c_str(), label.c_str() );
    const std::string qa_outfile = Form( "%s/rho_calib_%s_qa.root", outdir.c_str(), label.c_str() );
    const bool run_debug = false;

    if ( zvtx_edges.size() < 2 || mbdQ_edges.size() < 2 )
    {
        std::cerr << "Error: zvtx_edges and mbdQ_edges must each have at least 2 entries." << std::endl;
        return -1;
    }
    const int n_zbins = static_cast<int>( zvtx_edges.size() ) - 1;
    const int n_mbdbins = static_cast<int>( mbdQ_edges.size() ) - 1;

    std::cout << "Input file:  " << infile << std::endl;
    std::cout << "Output file: " << outfile << std::endl;
    std::cout << "n_zvtx_bins = " << n_zbins << ", n_mbdQ_bins = " << n_mbdbins << std::endl;

    auto * f = TFile::Open( infile.c_str() );
    if ( !f || f -> IsZombie() )
    {
        std::cerr << "Error: could not open " << infile << std::endl;
        return -1;
    }
    auto * t = (TTree*) f -> Get( "T" );
    if ( !t )
    {
        std::cerr << "Error: no tree \"T\" in " << infile << std::endl;
        return -1;
    }
    const Long64_t nentries = t -> GetEntries();
    if ( nentries < 1 )
    {
        std::cerr << "Error: no entries found in chain." << std::endl;
        return -1;
    }
    std::cout << "Total entries: " << nentries << std::endl;

    t -> SetBranchStatus( "*", false );

    static const int kNCalo = 3; // 0 = cemc, 1 = ihcal, 2 = ohcal
    static const int NETA = 24;

    float zvrtx = 0.0;
    float mbdQ = 0.0;
    float rho[kNCalo];
    float mean_E[kNCalo][NETA];
    int ntowers[kNCalo][NETA];
    float corrected_eta[kNCalo][NETA];

    bool ok = true;
    auto enable = [&]( const char * name, void * addr )
    {
        if ( !t -> GetBranch( name ) )
        {
            std::cerr << "Error: branch \"" << name << "\" not found in input chain." << std::endl;
            ok = false;
            return;
        }
        t -> SetBranchStatus( name, true );
        t -> SetBranchAddress( name, addr );
    };

    enable( "zvrtx", &zvrtx );
    enable( "mbdQ", &mbdQ );
    enable( "rho", rho );
    enable( "mean_E", mean_E );
    enable( "ntowers", ntowers );
    enable( "corrected_eta", corrected_eta );

    if ( !ok )
    {
        std::cerr << "Error: input chain is missing required branches -- this macro needs trees "
                  << "written by read_tree.C." << std::endl;
        return -1;
    }

    static const char * calo_names[kNCalo] = { "cemc", "hcalin", "hcalout" };
    static const char * calib_field_names[kNCalo] = { "w_cemc", "w_hcalin", "w_hcalout" };

    // matches SubtractTowersRhov1::find_bin -- half-open [edges[i],edges[i+1]),
    // -1 if val is outside [edges.front(),edges.back())
    auto find_bin = []( const float val, const std::vector<float> & edges ) -> int
    {
        if ( edges.size() < 2 || std::isnan(val) || val < edges.front() || val >= edges.back() )
        {
            return -1;
        }
        for ( size_t i = 0; i + 1 < edges.size(); ++i )
        {
            if ( val >= edges[i] && val < edges[i + 1] ) return static_cast<int>(i);
        }
        return -1;
    };

    // matches SubtractTowersRhov1::encode_channel
    auto encode_channel = []( const int ieta, const int izbin, const int imbd, const int n_mbd_bins ) -> int
    {
        return izbin * ( n_mbd_bins * NETA ) + imbd * NETA + ieta;
    };

    auto bin_index = [&]( const int izbin, const int imbd, const int ic ) -> int
    {
        return ( izbin * n_mbdbins + imbd ) * kNCalo + ic;
    };

    std::vector<TProfile*> pWeight( n_zbins * n_mbdbins * kNCalo, nullptr );
    std::vector<TProfile*> pEta( n_zbins * n_mbdbins * kNCalo, nullptr );
    for ( int iz = 0; iz < n_zbins; ++iz )
    {
        for ( int imbd = 0; imbd < n_mbdbins; ++imbd )
        {
            for ( int ic = 0; ic < kNCalo; ++ic )
            {
                const int idx = bin_index( iz, imbd, ic );
                pWeight[idx] = new TProfile(
                    Form( "calib_%s_z%d_q%d", calo_names[ic], iz, imbd ),
                    Form( "rho calibration weight, %s, zvtx#in[%.1f,%.1f), mbdQ#in[%.1f,%.1f);i_{#eta};w = <E_{tower}> / (#rho cosh(#eta_{corr}))",
                          calo_names[ic], zvtx_edges[iz], zvtx_edges[iz+1], mbdQ_edges[imbd], mbdQ_edges[imbd+1] ),
                    NETA, -0.5, NETA - 0.5 );
                pEta[idx] = new TProfile(
                    Form( "eta_%s_z%d_q%d", calo_names[ic], iz, imbd ),
                    Form( "mean corrected #eta, %s, zvtx#in[%.1f,%.1f), mbdQ#in[%.1f,%.1f);i_{#eta};<#eta_{corr}>",
                          calo_names[ic], zvtx_edges[iz], zvtx_edges[iz+1], mbdQ_edges[imbd], mbdQ_edges[imbd+1] ),
                    NETA, -0.5, NETA - 0.5 );
            }
        }
    }

    Long64_t n_events_used = 0;
    Long64_t n_events_outside_bins = 0;
    for ( Long64_t i = 0; i < nentries; ++i )
    {
        if ( run_debug && i >= 1000 ) break;

        t -> GetEntry( i );
        if ( nentries >= 10 && i % ( nentries / 10 ) == 0 && i > 0 )
        {
            std::cout << "Processing entry " << i << " / " << nentries << std::endl;
        }

        const int izbin = find_bin( zvrtx, zvtx_edges );
        const int imbd = find_bin( mbdQ, mbdQ_edges );
        if ( izbin < 0 || imbd < 0 )
        {
            ++n_events_outside_bins;
            continue;
        }

        bool used = false;
        for ( int ic = 0; ic < kNCalo; ++ic )
        {
            if ( rho[ic] <= 0.0 ) continue;
            const int idx = bin_index( izbin, imbd, ic );
            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                if ( ntowers[ic][ieta] <= 0 ) continue;
                const float denom = rho[ic] * std::cosh( corrected_eta[ic][ieta] );
                if ( denom <= 0.0 ) continue;
                const float ratio = mean_E[ic][ieta] / denom;
                pWeight[idx] -> Fill( ieta, ratio, ntowers[ic][ieta] );
                pEta[idx] -> Fill( ieta, corrected_eta[ic][ieta], ntowers[ic][ieta] );
                used = true;
            }
        }
        if ( used ) ++n_events_used;
    }
    std::cout << n_events_used << " / " << nentries << " events contributed to the calibration ("
              << n_events_outside_bins << " outside the configured zvtx/mbdQ bins)." << std::endl;

    // -- companion QA file with the raw TProfiles --
    auto * fqa = new TFile( qa_outfile.c_str(), "RECREATE" );
    fqa -> cd();
    for ( auto * p : pWeight ) p -> Write();
    for ( auto * p : pEta ) p -> Write();
    fqa -> Close();
    std::cout << "Wrote " << qa_outfile << std::endl;

    // -- final calibration, in CDBTTree format for SubtractTowersRhov1 --
    auto * cdbttree = new CDBTTree( outfile );
    cdbttree -> SetSingleIntValue( "n_eta", NETA );
    cdbttree -> SetSingleIntValue( "n_zvtx_bins", n_zbins );
    cdbttree -> SetSingleIntValue( "n_mbdQ_bins", n_mbdbins );
    for ( int i = 0; i <= n_zbins; ++i )
    {
        cdbttree -> SetSingleFloatValue( Form( "zvtx_edge_%d", i ), zvtx_edges[i] );
    }
    for ( int i = 0; i <= n_mbdbins; ++i )
    {
        cdbttree -> SetSingleFloatValue( Form( "mbdQ_edge_%d", i ), mbdQ_edges[i] );
    }
    cdbttree -> CommitSingle();

    Long64_t n_empty_channels = 0;
    for ( int iz = 0; iz < n_zbins; ++iz )
    {
        for ( int imbd = 0; imbd < n_mbdbins; ++imbd )
        {
            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                const int channel = encode_channel( ieta, iz, imbd, n_mbdbins );
                for ( int ic = 0; ic < kNCalo; ++ic )
                {
                    const int idx = bin_index( iz, imbd, ic );
                    const int bin = ieta + 1;
                    const double n_eff = pWeight[idx] -> GetBinEntries( bin );
                    // empty (izbin,imbd,ieta) channels fall back to w = 1 (no
                    // correction) -- get_etaWeight() treats any non-positive
                    // weight the same way, but write it explicitly here.
                    float w = 1.0f;
                    if ( n_eff > 0 )
                    {
                        w = static_cast<float>( pWeight[idx] -> GetBinContent( bin ) );
                        if ( !( w > 0.0f ) ) w = 1.0f;
                    }
                    else
                    {
                        ++n_empty_channels;
                    }
                    cdbttree -> SetFloatValue( channel, calib_field_names[ic], w );
                }
            }
        }
    }
    cdbttree -> Commit();
    cdbttree -> WriteCDBTTree();
    delete cdbttree;

    std::cout << "Wrote " << outfile << " (" << n_empty_channels << " / " << ( n_zbins * n_mbdbins * NETA * kNCalo )
              << " channels had no entries and were set to w = 1)" << std::endl;

    return 0;
}

#endif
