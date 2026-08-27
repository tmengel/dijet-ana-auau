
#ifndef _MATCH_C_
#define _MATCH_C_

#include <myana/AnaUtils.h>

#include <sPhenixStyle.C>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <vector>
#include <iostream>

R__LOAD_LIBRARY( libmyana.so )

// matches truth jets to rho-subtracted and sub1(seeded)-subtracted reco
// jets and writes out one flat MatchedJetRow per truth/reco relationship:
// matched pairs first, then unmatched truth jets, then unmatched reco
// jets -- separately for each reco collection (reco_type: 0 = rho, 1 = sub1).
int match(
    const std::string & infile = "/sphenix/user/tmengel/dijet-ana-auau/catalog/jet10_hijing_scaled/jet10_hijing_scaled-000.list",
    const std::string & outfile = "output.root",
    const float min_jet_pt = 5.0
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

    t -> SetBranchStatus( "*", false ); // disable all branches

    int event_id = -1;
    int is_minbias = 0;
    float zvrtx = 0.0;
    int cent = -1;
    float psi2 = 0.0;
    float mbd_q_N = 0.0, mbd_q_S = 0.0;

    static const int NETA = 24;
    static const int NPHI = 64;
    float cemc_tower_E[NETA][NPHI];
    int cemc_tower_isgood[NETA][NPHI];
    float ihcal_tower_E[NETA][NPHI];
    int ihcal_tower_isgood[NETA][NPHI];
    float ohcal_tower_E[NETA][NPHI];
    int ohcal_tower_isgood[NETA][NPHI];

    float truth_jet_R = 0.0;
    std::vector<float> * truth_jet_E = nullptr;
    std::vector<float> * truth_jet_eta = nullptr;
    std::vector<float> * truth_jet_phi = nullptr;
    std::vector<float> * truth_jet_pT = nullptr;
    std::vector<int>   * truth_jet_flavor = nullptr;
    std::vector<float> * truth_jet_parton_pT = nullptr;
    std::vector<float> * truth_jet_parton_dr = nullptr;
    float truth_jet_maxpT_r04 = 0.0;

    float rho_jet_R = 0.0;
    std::vector<float> * rho_jet_E = nullptr;
    std::vector<float> * rho_jet_eta = nullptr;
    std::vector<float> * rho_jet_phi = nullptr;
    std::vector<float> * rho_jet_pT = nullptr;
    std::vector<float> * rho_jet_unsub_E = nullptr;
    std::vector<float> * rho_jet_unsub_pT = nullptr;

    float sub1_jet_R = 0.0;
    std::vector<float> * sub1_jet_E = nullptr;
    std::vector<float> * sub1_jet_eta = nullptr;
    std::vector<float> * sub1_jet_phi = nullptr;
    std::vector<float> * sub1_jet_pT = nullptr;
    std::vector<float> * sub1_jet_unsub_E = nullptr;
    std::vector<float> * sub1_jet_unsub_pT = nullptr;

    auto enable = [&]( const char * name, void * addr )
    {
        t -> SetBranchStatus( name, true );
        t -> SetBranchAddress( name, addr );
    };

    enable( "event_id", &event_id );
    enable( "is_minbias", &is_minbias );
    enable( "zvrtx", &zvrtx );
    enable( "cent", &cent );
    enable( "psi2", &psi2 );
    enable( "mbd_q_N", &mbd_q_N );
    enable( "mbd_q_S", &mbd_q_S );

    enable( "cemc_tower_E", cemc_tower_E );
    enable( "cemc_tower_isgood", cemc_tower_isgood );
    enable( "ihcal_tower_E", ihcal_tower_E );
    enable( "ihcal_tower_isgood", ihcal_tower_isgood );
    enable( "ohcal_tower_E", ohcal_tower_E );
    enable( "ohcal_tower_isgood", ohcal_tower_isgood );

    enable( "truth_jet_R", &truth_jet_R );
    enable( "truth_jet_E", &truth_jet_E );
    enable( "truth_jet_eta", &truth_jet_eta );
    enable( "truth_jet_phi", &truth_jet_phi );
    enable( "truth_jet_pT", &truth_jet_pT );
    enable( "truth_jet_flavor", &truth_jet_flavor );
    enable( "truth_jet_parton_pT", &truth_jet_parton_pT );
    enable( "truth_jet_parton_dr", &truth_jet_parton_dr );
    enable( "truth_jet_maxpT_r04", &truth_jet_maxpT_r04 );

    const bool has_rho = t -> GetBranch( "rho_jet_pT" ) != nullptr;
    const bool has_sub1 = t -> GetBranch( "sub1_jet_pT" ) != nullptr && false; // disable sub1 for now

    if ( has_rho )
    {
        enable( "rho_jet_R", &rho_jet_R );
        enable( "rho_jet_E", &rho_jet_E );
        enable( "rho_jet_eta", &rho_jet_eta );
        enable( "rho_jet_phi", &rho_jet_phi );
        enable( "rho_jet_pT", &rho_jet_pT );
        enable( "rho_jet_unsub_E", &rho_jet_unsub_E );
        enable( "rho_jet_unsub_pT", &rho_jet_unsub_pT );
    }
    if ( has_sub1 )
    {
        enable( "sub1_jet_R", &sub1_jet_R );
        enable( "sub1_jet_E", &sub1_jet_E );
        enable( "sub1_jet_eta", &sub1_jet_eta );
        enable( "sub1_jet_phi", &sub1_jet_phi );
        enable( "sub1_jet_pT", &sub1_jet_pT );
        enable( "sub1_jet_unsub_E", &sub1_jet_unsub_E );
        enable( "sub1_jet_unsub_pT", &sub1_jet_unsub_pT );
    }

    if ( !has_rho && !has_sub1 )
    {
        std::cerr << "Error: input tree has neither rho_jet nor sub1_jet branches." << std::endl;
        return -1;
    }

    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * tout = new TTree( "T", "Matched truth/reco jets" );

    AnaUtils::MatchedJetRow row;
    AnaUtils::book_matched_jet_tree( tout, row );

    std::cout << "Minimum jet pT: " << min_jet_pt << std::endl;

    long n_rows = 0;
    long n_matched_rho = 0, n_matched_sub1 = 0;

    for ( int ientry = 0; ientry < nentries; ++ientry )
    {
        t -> GetEntry( ientry );
        if ( nentries >= 10 && ientry % ( nentries / 10 ) == 0 && ientry > 0 )
        {
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        }

        row.event_id = event_id;
        row.cent = cent;
        row.zvrtx = zvrtx;
        row.mbdQ = mbd_q_N + mbd_q_S;
        row.sumeT = AnaUtils::calc_sumeT( AnaUtils::CEMC, zvrtx, cemc_tower_E, cemc_tower_isgood )
                  + AnaUtils::calc_sumeT( AnaUtils::HCALIN, zvrtx, ihcal_tower_E, ihcal_tower_isgood )
                  + AnaUtils::calc_sumeT( AnaUtils::HCALOUT, zvrtx, ohcal_tower_E, ohcal_tower_isgood );
        row.is_minbias = is_minbias;
        row.psi2 = psi2;
        row.truth_jet_maxpt_r04 = truth_jet_maxpT_r04;

        auto truth_sel = AnaUtils::select_jets(
            *truth_jet_pT, *truth_jet_E, *truth_jet_eta,
            min_jet_pt, zvrtx, truth_jet_R, false
        );

        // reco_type 0 = rho, reco_type 1 = sub1
        for ( int reco_type = 0; reco_type < 2; ++reco_type )
        {
            if ( reco_type == 0 && !has_rho ) continue;
            if ( reco_type == 1 && !has_sub1 ) continue;

            std::vector<float> * reco_E   = ( reco_type == 0 ) ? rho_jet_E   : sub1_jet_E;
            std::vector<float> * reco_eta = ( reco_type == 0 ) ? rho_jet_eta : sub1_jet_eta;
            std::vector<float> * reco_phi = ( reco_type == 0 ) ? rho_jet_phi : sub1_jet_phi;
            std::vector<float> * reco_pT  = ( reco_type == 0 ) ? rho_jet_pT  : sub1_jet_pT;
            std::vector<float> * reco_unsub_E  = ( reco_type == 0 ) ? rho_jet_unsub_E  : sub1_jet_unsub_E;
            std::vector<float> * reco_unsub_pT = ( reco_type == 0 ) ? rho_jet_unsub_pT : sub1_jet_unsub_pT;
            const float jet_R = ( reco_type == 0 ) ? rho_jet_R : sub1_jet_R;
            const float max_dr = 0.75 * jet_R;

            auto reco_sel = AnaUtils::select_jets(
                *reco_pT, *reco_E, *reco_eta,
                min_jet_pt, zvrtx, jet_R, true
            );

            auto matches = AnaUtils::match_truth_reco_jets(
                truth_sel, *truth_jet_eta, *truth_jet_phi,
                reco_sel, *reco_eta, *reco_phi,
                max_dr
            );

            for ( const auto & jm : matches )
            {
                row.reco_type = reco_type;
                row.dr = jm.dr;

                const bool has_truth = jm.truth_index >= 0;
                const bool has_reco  = jm.reco_index >= 0;

                if ( has_truth && has_reco ) { row.match_status = 0; if ( reco_type == 0 ) ++n_matched_rho; else ++n_matched_sub1; }
                else if ( has_truth )        { row.match_status = 1; }
                else                         { row.match_status = 2; }

                if ( has_truth )
                {
                    const int ti = jm.truth_index;
                    row.truth_pt = truth_jet_pT->at(ti);
                    row.truth_e = truth_jet_E->at(ti);
                    row.truth_eta = truth_jet_eta->at(ti);
                    row.truth_phi = truth_jet_phi->at(ti);
                    row.truth_flavor = truth_jet_flavor->at(ti);
                    row.truth_parton_pt = truth_jet_parton_pT->at(ti);
                    row.truth_parton_dr = truth_jet_parton_dr->at(ti);
                }
                else
                {
                    row.truth_pt = -999.0;
                    row.truth_e = -999.0;
                    row.truth_eta = -999.0;
                    row.truth_phi = -999.0;
                    row.truth_flavor = -999;
                    row.truth_parton_pt = -999.0;
                    row.truth_parton_dr = -999.0;
                }

                if ( has_reco )
                {
                    const int ri = jm.reco_index;
                    row.reco_pt = reco_pT->at(ri);
                    row.reco_e = reco_E->at(ri);
                    row.reco_eta = reco_eta->at(ri);
                    row.reco_phi = reco_phi->at(ri);
                    row.reco_unsub_e = reco_unsub_E->at(ri);
                    row.reco_unsub_pt = reco_unsub_pT->at(ri);
                }
                else
                {
                    row.reco_pt = -999.0;
                    row.reco_e = -999.0;
                    row.reco_eta = -999.0;
                    row.reco_phi = -999.0;
                    row.reco_unsub_e = -999.0;
                    row.reco_unsub_pt = -999.0;
                }

                tout -> Fill();
                ++n_rows;
            }
        }
    }

    std::cout << "Wrote " << n_rows << " rows (" << n_matched_rho << " matched rho pairs, "
              << n_matched_sub1 << " matched sub1 pairs)" << std::endl;

    fout -> cd();
    tout -> Write();
    fout -> Close();

    std::cout << "Wrote " << outfile << std::endl;

    return 0;
}

#endif
