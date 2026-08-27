
#ifndef _READ_TREE_C_
#define _READ_TREE_C_

#include <myana/AnaUtils.h>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include <cmath>
#include <iostream>
#include <string>

R__LOAD_LIBRARY( libmyana.so )

// reads per-tower energy/goodness, event centrality, z-vertex, MBD charge,
// sumeT, and rho underlying-event density from Au+Au data trees (written
// by AnaTree with save_full_calo(cemc, hcalin, hcalout, true) and
// save_sum_eT(cemc, hcalin, hcalout, true) enabled -- layer index 0/1/2 =
// CEMC/HCALIN/HCALOUT throughout, matching AnaUtils::CaloType), and writes
// out one flat row per event with the per-(layer, eta ring) tower-energy
// statistics needed to derive the rho calibration weight
//   w = <E_tower> / (rho * cosh(eta_corr))
// (PPG-14 Eq. 1) offline, binned however is convenient after the fact
// (e.g. by zvrtx, mbdQ, cent). events failing the sPHENIX minimum-bias
// requirement or outside |zvrtx| < max_abs_zvrtx are discarded. infile
// may be a file list (read via AnaUtils::getFilelist) or a single ROOT
// file.
int read_tree(
    const std::string & infile = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/lists/54361.list",
    const std::string & outfile = "rho_calib_tree.root",
    const float max_abs_zvrtx = 60.0
)
{
    auto * t = new TChain( "T" );
    auto files = AnaUtils::getFilelist( infile, ".root" );
    if ( !files.empty() )
    {
        for ( const auto & file : files )
        {
            std::cout << "Adding file: " << file << std::endl;
            t -> Add( file.c_str() );
        }
    }
    else
    {
        t -> Add( infile.c_str() );
    }

    const int nentries = t -> GetEntries();
    if ( nentries < 1 )
    {
        std::cerr << "Error: no entries found for " << infile << std::endl;
        return -1;
    }
    std::cout << "Total entries: " << nentries << std::endl;

    t -> SetBranchStatus( "*", false );

    int is_minbias = 0;
    float zvrtx = 0.0;
    int cent = -1;
    float mbd_q_N = 0.0, mbd_q_S = 0.0;

    static const int kNCalo = 3; // 0 = cemc, 1 = ihcal, 2 = ohcal
    static const int NETA = 24;
    static const int NPHI = 64;

    float tower_E[kNCalo][NETA][NPHI];
    int tower_isgood[kNCalo][NETA][NPHI];
    float rho_vals[kNCalo];
    float sumeT[kNCalo];

    bool ok = true;
    auto enable = [&]( const char * name, void * addr )
    {
        if ( !t -> GetBranch( name ) )
        {
            std::cerr << "Error: branch \"" << name << "\" not found in input tree." << std::endl;
            ok = false;
            return;
        }
        t -> SetBranchStatus( name, true );
        t -> SetBranchAddress( name, addr );
    };

    enable( "is_minbias", &is_minbias );
    enable( "zvrtx", &zvrtx );
    enable( "cent", &cent );
    enable( "mbd_q_N", &mbd_q_N );
    enable( "mbd_q_S", &mbd_q_S );
    enable( "tower_E", tower_E );
    enable( "tower_isgood", tower_isgood );
    enable( "rho_vals", rho_vals );
    enable( "sumeT", sumeT );

    if ( !ok )
    {
        std::cerr << "Error: input tree is missing required branches -- this macro needs a tree "
                  << "written by AnaTree with save_full_calo(cemc, hcalin, hcalout, true) and "
                  << "save_sum_eT(cemc, hcalin, hcalout, true) enabled." << std::endl;
        return -1;
    }

    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * tout = new TTree( "T", "per-event tower statistics for rho calibration" );

    int out_cent = -1;
    float out_zvrtx = 0.0;
    float out_mbdQ = 0.0;
    float out_sumeT[kNCalo];
    float out_rho[kNCalo];
    float out_mean_E[kNCalo][NETA];
    int   out_ntowers[kNCalo][NETA];
    float out_corrected_eta[kNCalo][NETA];

    tout -> Branch( "cent", &out_cent, "cent/I" );
    tout -> Branch( "zvrtx", &out_zvrtx, "zvrtx/F" );
    tout -> Branch( "mbdQ", &out_mbdQ, "mbdQ/F" );
    tout -> Branch( "sumeT", out_sumeT, "sumeT[3]/F" );
    tout -> Branch( "rho", out_rho, "rho[3]/F" );
    tout -> Branch( "mean_E", out_mean_E, Form( "mean_E[3][%d]/F", NETA ) );
    tout -> Branch( "ntowers", out_ntowers, Form( "ntowers[3][%d]/I", NETA ) );
    tout -> Branch( "corrected_eta", out_corrected_eta, Form( "corrected_eta[3][%d]/F", NETA ) );

    const AnaUtils::CaloType calo_types[kNCalo] = { AnaUtils::CEMC, AnaUtils::HCALIN, AnaUtils::HCALOUT };

    long n_events_used = 0;

    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );
        if ( nentries >= 10 && i % ( nentries / 10 ) == 0 && i > 0 )
        {
            std::cout << "Processing entry " << i << " / " << nentries << std::endl;
        }

        if ( !is_minbias ) continue;
        if ( std::fabs( zvrtx ) >= max_abs_zvrtx ) continue;

        out_cent = cent;
        out_zvrtx = zvrtx;
        out_mbdQ = mbd_q_N + mbd_q_S;

        for ( int ic = 0; ic < kNCalo; ++ic )
        {
            out_sumeT[ic] = sumeT[ic];
            out_rho[ic] = rho_vals[ic];

            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                out_corrected_eta[ic][ieta] = AnaUtils::get_corrected_calo_eta( calo_types[ic], ieta, zvrtx );

                float sum_E = 0.0;
                int n_good = 0;
                for ( int iphi = 0; iphi < NPHI; ++iphi )
                {
                    if ( !tower_isgood[ic][ieta][iphi] ) continue;
                    sum_E += tower_E[ic][ieta][iphi];
                    ++n_good;
                }
                out_ntowers[ic][ieta] = n_good;
                out_mean_E[ic][ieta] = ( n_good > 0 ) ? sum_E / n_good : 0.0f;
            }
        }

        tout -> Fill();
        ++n_events_used;
    }

    fout -> cd();
    tout -> Write();
    fout -> Close();

    std::cout << n_events_used << " / " << nentries << " events written (minbias, |zvrtx| < "
              << max_abs_zvrtx << " cm)." << std::endl;
    std::cout << "Wrote " << outfile << std::endl;

    return 0;
}

#endif
