
#ifndef _READ_TREE_SIM_C_
#define _READ_TREE_SIM_C_

#include <myana/AnaUtils.h>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include <cmath>
#include <iostream>
#include <string>

R__LOAD_LIBRARY( libmyana.so )

// Same as read_tree.C, but for the mb_hijing simulation CALO_TREE files
// (catalog/mb_hijing_{scaled,unsclaed}/*.list), which use a different branch
// layout than the data anatree: per-layer branches (cemc_tower_E[24][64],
// ihcal_tower_E[24][64], ohcal_tower_E[24][64], and matching *_tower_isgood
// / rho_val_TowerRho_MULT_* / sumeT_* names) instead of the data's combined
// tower_E[3][24][64] / rho_vals[3] / sumeT[3]. This macro maps those onto the
// exact same internal (layer, eta ring) representation and writes an
// IDENTICAL output schema to read_tree.C's, so make_rho_calib_sim.C (the
// only other change needed) can reuse the same calibration-derivation logic
// unchanged. Layer index 0/1/2 = CEMC/HCALIN/HCALOUT, matching
// AnaUtils::CaloType, same as the data pipeline.
int read_tree_sim(
    const std::string & infile = "/sphenix/user/tmengel/dijet-ana-auau/catalog/mb_hijing_scaled/mb_hijing_scaled-000.list",
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

    // sim tree: separate per-layer branches, not combined [3][24][64]/[3] arrays
    float cemc_tower_E[NETA][NPHI], ihcal_tower_E[NETA][NPHI], ohcal_tower_E[NETA][NPHI];
    int cemc_tower_isgood[NETA][NPHI], ihcal_tower_isgood[NETA][NPHI], ohcal_tower_isgood[NETA][NPHI];
    float rho_val_cemc = 0.0, rho_val_ihcal = 0.0, rho_val_ohcal = 0.0;
    float sumeT_cemc = 0.0, sumeT_ihcal = 0.0, sumeT_ohcal = 0.0;

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
    enable( "cemc_tower_E", cemc_tower_E );
    enable( "ihcal_tower_E", ihcal_tower_E );
    enable( "ohcal_tower_E", ohcal_tower_E );
    enable( "cemc_tower_isgood", cemc_tower_isgood );
    enable( "ihcal_tower_isgood", ihcal_tower_isgood );
    enable( "ohcal_tower_isgood", ohcal_tower_isgood );
    enable( "rho_val_TowerRho_MULT_CEMC", &rho_val_cemc );
    enable( "rho_val_TowerRho_MULT_HCALIN", &rho_val_ihcal );
    enable( "rho_val_TowerRho_MULT_HCALOUT", &rho_val_ohcal );
    enable( "sumeT_cemc", &sumeT_cemc );
    enable( "sumeT_ihcal", &sumeT_ihcal );
    enable( "sumeT_ohcal", &sumeT_ohcal );

    if ( !ok )
    {
        std::cerr << "Error: input tree is missing required branches for the mb_hijing sim CALO_TREE layout." << std::endl;
        return -1;
    }

    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * tout = new TTree( "T", "per-event tower statistics for rho calibration (sim)" );

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

        out_sumeT[0] = sumeT_cemc; out_sumeT[1] = sumeT_ihcal; out_sumeT[2] = sumeT_ohcal;
        out_rho[0] = rho_val_cemc; out_rho[1] = rho_val_ihcal; out_rho[2] = rho_val_ohcal;

        for ( int ic = 0; ic < kNCalo; ++ic )
        {
            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                out_corrected_eta[ic][ieta] = AnaUtils::get_corrected_calo_eta( calo_types[ic], ieta, zvrtx );

                float sum_E = 0.0;
                int n_good = 0;
                for ( int iphi = 0; iphi < NPHI; ++iphi )
                {
                    int isgood = 0; float E = 0.0;
                    if ( ic == 0 ) { isgood = cemc_tower_isgood[ieta][iphi]; E = cemc_tower_E[ieta][iphi]; }
                    else if ( ic == 1 ) { isgood = ihcal_tower_isgood[ieta][iphi]; E = ihcal_tower_E[ieta][iphi]; }
                    else { isgood = ohcal_tower_isgood[ieta][iphi]; E = ohcal_tower_E[ieta][iphi]; }
                    if ( !isgood ) continue;
                    sum_E += E;
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
