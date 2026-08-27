
#ifndef _DIJET_PAIR_MATCHING_C_
#define _DIJET_PAIR_MATCHING_C_

#include <myana/AnaUtils.h>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

// Pair-level dijet matching / response bookkeeping, written to sit
// alongside dijet_matching.C rather than replace it, so the two outputs
// can be compared on the same input.
//
// Differences from dijet_matching.C:
//
//  (1) The classification is PAIR-level (one row per event) instead of
//      per-leg. A dijet response (x_J, A_J, pT1-vs-pT2) is a pair-level
//      object -- a per-leg row set has no defined meaning when the two
//      legs land in different categories.
//
//  (2) The reco side is built the way the data analysis builds it: the
//      leading and subleading ACCEPTED reco jets by pT (same convention as
//      makeMatchedTreesTaggedAuAu.C), then the pT and dphi requirements.
//      dijet_matching.C instead only ever looked at each truth leg's own
//      dR match, so a UE jet that outranks a real leg could never appear
//      in its output at all.
//
//  (3) A UE fluctuation that has taken over a dijet leg gets its own
//      category (kUESub) instead of being silently absorbed into Fill or
//      counted as a Miss. Those pairs are excluded from the response
//      matrix AND from fakes AND from misses -- they are the population
//      the subleading-jet efficiency and the inclusive cross-check cover,
//      so counting them here as well would double count them.
//
// The truth dijet is deliberately selected with NO eta cut: the truth
// candidate is always truth jets 0 and 1 (the two hardest truth jets in
// the event, already pT-sorted by match_standalone.C). Keeping the truth
// pair unconditional is what makes it possible to tell a reco jet that is
// unmatched because its truth partner was excluded by the calorimeter
// acceptance from a reco jet that is unmatched because it is a UE
// fluctuation: the per-leg truth_fail bitmask records WHY the truth leg
// failed, and reco_prov records WHAT each reco leg actually is.
//
// Acceptance definitions (unchanged from dijet_matching.C):
//   truth leg 1(2): pT > 14(7) GeV, truth_jet_accept_eta
//   reco  leg 1(2): pT > 20(8) GeV, jet_accept_eta, E > 0
//   both truth and reco pairs additionally require |dphi(1,2)| > 7pi/8
//
// Input: the tree written by match_standalone.C (needs truth_jet_*,
// jet_*, truth_jet_reco_match_idx and jet_truth_match_idx).
//
// The classification itself lives in DijetPair::classify() below, which
// takes plain vectors and no ROOT I/O state, so it can be lifted into
// match_standalone.C as-is to do the matching and the categorization in
// one pass.
namespace DijetPair
{
    // Pair-level category. 0-3 keep the same numbering as the per-leg
    // categories in dijet_matching.C so the two outputs stay comparable;
    // 4 is the new one that keeps UE-substituted pairs out of the other
    // four.
    enum Category
    {
        kFill  = 0, // truth pair in acceptance; reco dijet exists and its
                    // two legs ARE the two truth legs' own matches (in
                    // either order -- see legs_swapped). The only category
                    // that fills the response matrix.
        kMiss  = 1, // truth pair in acceptance; no reco dijet at all.
        kFake  = 2, // truth pair NOT in acceptance; reco dijet exists.
        kSkip  = 3, // truth pair NOT in acceptance; no reco dijet.
        kUESub = 4  // truth pair in acceptance; reco dijet exists, but at
                    // least one leg is not the corresponding truth leg's
                    // match. Neither a Fill, nor a Miss, nor a Fake.
    };

    // What each selected reco dijet leg actually is, from
    // jet_truth_match_idx.
    enum Prov
    {
        kProvNone       = -1, // no reco leg (no reco dijet was formed)
        kProvLead       =  0, // this reco jet is truth leg 1's match
        kProvSub        =  1, // this reco jet is truth leg 2's match
        kProvOtherTruth =  2, // matched to a truth jet, but a softer one
                              // that is not part of the leading truth pair
        kProvUE         =  3  // no truth match at all -> UE fluctuation
    };

    // Bitmask recording WHY a truth leg is out of acceptance. Zero means
    // the leg is in acceptance. kFailDphi is a pair-level condition and so
    // is set on both legs together.
    enum TruthFail
    {
        kFailNone = 0,
        kFailPt   = 1,
        kFailEta  = 2,
        kFailDphi = 4
    };

    struct Config
    {
        float truth_pt_thresh[2] = { 14.0f, 7.0f };
        float reco_pt_thresh[2]  = { 20.0f, 8.0f };
        float min_dphi           = 7.0f * TMath::Pi() / 8.0f;
        // false (default): subleading reco jet is the 2nd-highest-pT
        //   accepted jet, then dphi is required of that pair -- the
        //   convention used in makeMatchedTreesTaggedAuAu.C, and the one
        //   under which a UE jet can actually steal a leg.
        // true: subleading reco jet is instead the highest-pT accepted jet
        //   above threshold that is already back-to-back with the leading
        //   one. Provided as a systematic; it hides part of the UE
        //   substitution by construction.
        bool sublead_back_to_back = false;
    };

    struct Result
    {
        int   category      = kSkip;
        int   truth_in_acc  = 0;
        int   truth_fail[2] = { kFailNone, kFailNone };
        int   truth_idx[2]  = { -1, -1 };
        float truth_dphi    = -999.0f;
        // each truth leg's own dR match, kept even when it is not one of
        // the selected reco dijet legs -- this is what shows "leg 2's real
        // match was at 6 GeV and a 12 GeV UE jet took its place".
        int   truth_match_idx[2] = { -1, -1 };
        float truth_match_pt[2]  = { -999.0f, -999.0f };

        int   reco_pair       = 0;
        int   reco_idx[2]     = { -1, -1 };
        float reco_dphi       = -999.0f;
        int   reco_prov[2]    = { kProvNone, kProvNone };
        int   reco_truth_idx[2] = { -1, -1 }; // truth index each reco leg matched
        int   legs_swapped    = 0;
        int   n_accepted_reco = 0;
    };

    // Pure classification -- no ROOT I/O state, safe to lift into
    // match_standalone.C. Returns false if the event has no truth dijet
    // candidate at all (fewer than 2 truth jets), in which case nothing
    // should be written for it.
    inline bool classify(
        const std::vector< float > & truth_jet_pT,
        const std::vector< float > & truth_jet_phi,
        const std::vector< int >   & truth_jet_accept_eta,
        const std::vector< int >   & truth_jet_reco_match_idx,
        const std::vector< float > & jet_pT,
        const std::vector< float > & jet_phi,
        const std::vector< float > & jet_E,
        const std::vector< int >   & jet_accept_eta,
        const std::vector< int >   & jet_truth_match_idx,
        const Config & cfg,
        Result & r
    )
    {
        if ( truth_jet_pT.size() < 2 ) return false;

        //------------------------------------------------------------
        // truth pair: always jets 0 and 1, NO eta cut at selection time.
        //------------------------------------------------------------
        r.truth_idx[0] = 0;
        r.truth_idx[1] = 1;
        r.truth_dphi = AnaUtils::dphi_wrap( truth_jet_phi.at( 0 ), truth_jet_phi.at( 1 ) );
        const bool truth_dphi_ok = r.truth_dphi > cfg.min_dphi;

        for ( int leg = 0; leg < 2; ++leg )
        {
            const int ti = r.truth_idx[leg];
            int fail = kFailNone;
            if ( !( truth_jet_pT.at( ti ) > cfg.truth_pt_thresh[leg] ) ) fail |= kFailPt;
            if ( !truth_jet_accept_eta.at( ti ) )                        fail |= kFailEta;
            if ( !truth_dphi_ok )                                        fail |= kFailDphi;
            r.truth_fail[leg] = fail;

            const int mi = truth_jet_reco_match_idx.at( ti );
            r.truth_match_idx[leg] = mi;
            r.truth_match_pt[leg]  = ( mi >= 0 ) ? jet_pT.at( mi ) : -999.0f;
        }
        r.truth_in_acc = ( r.truth_fail[0] == kFailNone && r.truth_fail[1] == kFailNone ) ? 1 : 0;

        //------------------------------------------------------------
        // reco pair: built the way the data analysis builds it, from the
        // accepted reco jets in pT order -- NOT from the truth legs'
        // matches. This is what lets a UE jet take a leg.
        //------------------------------------------------------------
        std::vector< int > accepted;
        for ( size_t j = 0; j < jet_pT.size(); ++j )
        {
            if ( jet_accept_eta.at( j ) && jet_E.at( j ) > 0.0f )
            {
                accepted.push_back( static_cast< int >( j ) );
            }
        }
        // jet_pT is pT-descending on input (sort_reco_jets_by_pt in
        // match_standalone.C), so filtering preserves that order.
        for ( int j : accepted )
        {
            if ( jet_pT.at( j ) > cfg.reco_pt_thresh[1] ) ++r.n_accepted_reco;
        }

        int r1 = -1, r2 = -1;
        if ( !accepted.empty() && jet_pT.at( accepted[0] ) > cfg.reco_pt_thresh[0] )
        {
            r1 = accepted[0];
            if ( cfg.sublead_back_to_back )
            {
                for ( size_t k = 1; k < accepted.size(); ++k )
                {
                    const int cand = accepted[k];
                    if ( !( jet_pT.at( cand ) > cfg.reco_pt_thresh[1] ) ) break; // sorted
                    if ( AnaUtils::dphi_wrap( jet_phi.at( r1 ), jet_phi.at( cand ) ) > cfg.min_dphi )
                    {
                        r2 = cand;
                        break;
                    }
                }
            }
            else if ( accepted.size() >= 2
                      && jet_pT.at( accepted[1] ) > cfg.reco_pt_thresh[1]
                      && AnaUtils::dphi_wrap( jet_phi.at( r1 ), jet_phi.at( accepted[1] ) ) > cfg.min_dphi )
            {
                r2 = accepted[1];
            }
        }

        if ( r1 >= 0 && r2 >= 0 )
        {
            r.reco_pair = 1;
            r.reco_idx[0] = r1;
            r.reco_idx[1] = r2;
            r.reco_dphi = AnaUtils::dphi_wrap( jet_phi.at( r1 ), jet_phi.at( r2 ) );

            for ( int leg = 0; leg < 2; ++leg )
            {
                const int ri = r.reco_idx[leg];
                const int ti = jet_truth_match_idx.at( ri );
                r.reco_truth_idx[leg] = ti;
                if      ( ti < 0 )  r.reco_prov[leg] = kProvUE;
                else if ( ti == 0 ) r.reco_prov[leg] = kProvLead;
                else if ( ti == 1 ) r.reco_prov[leg] = kProvSub;
                else                r.reco_prov[leg] = kProvOtherTruth;
            }
        }

        //------------------------------------------------------------
        // pair-level category
        //------------------------------------------------------------
        if ( r.truth_in_acc )
        {
            if ( !r.reco_pair )
            {
                r.category = kMiss;
            }
            else
            {
                // both reco legs must be the two truth legs' own matches;
                // either ordering counts (a resolution-driven pT swap
                // between the two real legs is a genuine response entry,
                // not a UE effect -- legs_swapped records it).
                const bool same_pair =
                       ( r.reco_prov[0] == kProvLead && r.reco_prov[1] == kProvSub )
                    || ( r.reco_prov[0] == kProvSub  && r.reco_prov[1] == kProvLead );
                r.category = same_pair ? kFill : kUESub;
                r.legs_swapped = ( same_pair && r.reco_prov[0] == kProvSub ) ? 1 : 0;
            }
        }
        else
        {
            r.category = r.reco_pair ? kFake : kSkip;
        }

        return true;
    }
}

int dijet_pair_matching(
    const std::string & infile = "output.root",
    const std::string & outfile = "dijet_pair_matching.root",
    const bool sublead_back_to_back = false
)
{
    auto * t = new TChain( "T" );
    // infile can be either a single merged .root file (as written by
    // match_standalone.C) or a .list file listing several of them --
    // AnaUtils::getFilelist always reads its argument as a text file, so a
    // .root file must be added to the chain directly rather than routed
    // through it.
    if ( infile.size() >= 5 && infile.compare( infile.size() - 5, 5, ".root" ) == 0 )
    {
        std::cout << "Adding file: " << infile << std::endl;
        t -> Add( infile.c_str() );
    }
    else
    {
        auto files = AnaUtils::getFilelist( infile, ".root" );
        for ( const auto & file : files )
        {
            std::cout << "Adding file: " << file << std::endl;
            t -> Add( file.c_str() );
        }
    }

    int nentries = t -> GetEntries();
    if ( nentries < 1 )
    {
        std::cerr << "Error: No entries in TChain." << std::endl;
        return -1;
    }
    std::cout << "Total entries in TChain: " << nentries << std::endl;

    t -> SetBranchStatus( "*", false );

    auto has_branch = [ t ]( const std::string & name )
    {
        return t -> GetBranch( name.c_str() ) != nullptr;
    };
    auto enable = [ t ]( const std::string & name, void * addr )
    {
        t -> SetBranchStatus ( name.c_str(), true );
        t -> SetBranchAddress( name.c_str(), addr );
    };

    //----------------------------------------------------------------
    // event-level passthrough
    //----------------------------------------------------------------
    int event_id = -1;
    int cent = -1;
    float zvrtx = 0.0;
    int is_minbias = 0;
    float mbd_q = -999.0;
    float sumeT = -999.0;
    float psi2 = -999.0;

    const bool has_event_id = has_branch( "event_id" );
    const bool has_cent     = has_branch( "cent" );
    const bool has_zvrtx    = has_branch( "zvrtx" );
    const bool has_minbias  = has_branch( "is_minbias" );
    const bool has_mbd_q    = has_branch( "mbd_q" );
    const bool has_sumeT    = has_branch( "sumeT" );
    const bool has_psi2     = has_branch( "psi2" );

    if ( has_event_id ) enable( "event_id", &event_id );
    if ( has_cent )     enable( "cent", &cent );
    if ( has_zvrtx )    enable( "zvrtx", &zvrtx );
    if ( has_minbias )  enable( "is_minbias", &is_minbias );
    if ( has_mbd_q )    enable( "mbd_q", &mbd_q );
    if ( has_sumeT )    enable( "sumeT", &sumeT );
    if ( has_psi2 )     enable( "psi2", &psi2 );

    //----------------------------------------------------------------
    // truth / reco jets and the match_standalone.C matching branches
    //----------------------------------------------------------------
    std::vector< float > * truth_jet_pT  = nullptr;
    std::vector< float > * truth_jet_eta = nullptr;
    std::vector< float > * truth_jet_phi = nullptr;
    std::vector< float > * truth_jet_E   = nullptr;
    std::vector< int >   * truth_jet_accept_eta = nullptr;
    std::vector< int >   * truth_jet_reco_match_idx = nullptr;

    std::vector< float > * jet_pT  = nullptr;
    std::vector< float > * jet_eta = nullptr;
    std::vector< float > * jet_phi = nullptr;
    std::vector< float > * jet_E   = nullptr;
    std::vector< int >   * jet_accept_eta = nullptr;
    std::vector< int >   * jet_truth_match_idx = nullptr;

    const bool has_all = has_branch( "truth_jet_pT" )
                        && has_branch( "truth_jet_accept_eta" )
                        && has_branch( "truth_jet_reco_match_idx" )
                        && has_branch( "jet_pT" )
                        && has_branch( "jet_accept_eta" )
                        && has_branch( "jet_truth_match_idx" );
    if ( !has_all )
    {
        std::cerr << "Error: input tree is missing the truth/reco jet or "
                  << "matching branches written by match_standalone.C "
                  << "(truth_jet_reco_match_idx / jet_truth_match_idx)." << std::endl;
        return -1;
    }

    enable( "truth_jet_pT", &truth_jet_pT );
    enable( "truth_jet_eta", &truth_jet_eta );
    enable( "truth_jet_phi", &truth_jet_phi );
    enable( "truth_jet_E", &truth_jet_E );
    enable( "truth_jet_accept_eta", &truth_jet_accept_eta );
    enable( "truth_jet_reco_match_idx", &truth_jet_reco_match_idx );
    enable( "jet_pT", &jet_pT );
    enable( "jet_eta", &jet_eta );
    enable( "jet_phi", &jet_phi );
    enable( "jet_E", &jet_E );
    enable( "jet_accept_eta", &jet_accept_eta );
    enable( "jet_truth_match_idx", &jet_truth_match_idx );

    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * tout = new TTree( "T", "pair-level dijet matching / response" );

    int   o_event_id = -1, o_cent = -1, o_is_minbias = 0;
    float o_zvrtx = 0.0, o_mbd_q = -999.0, o_sumeT = -999.0;
    float o_psi2 = -999.0, o_dpsi2 = -999.0;

    int   o_category = -1, o_truth_in_acc = 0, o_reco_pair = 0;
    int   o_truth_fail1 = 0, o_truth_fail2 = 0;
    int   o_truth_idx1 = -1, o_truth_idx2 = -1;
    float o_truth_pt1 = -999.0, o_truth_eta1 = -999.0, o_truth_phi1 = -999.0, o_truth_E1 = -999.0;
    float o_truth_pt2 = -999.0, o_truth_eta2 = -999.0, o_truth_phi2 = -999.0, o_truth_E2 = -999.0;
    float o_truth_dphi = -999.0, o_truth_xj = -999.0;
    int   o_truth_match_idx1 = -1, o_truth_match_idx2 = -1;
    float o_truth_match_pt1 = -999.0, o_truth_match_pt2 = -999.0;
    int   o_reco_idx1 = -1, o_reco_idx2 = -1;
    float o_reco_pt1 = -999.0, o_reco_eta1 = -999.0, o_reco_phi1 = -999.0, o_reco_E1 = -999.0;
    float o_reco_pt2 = -999.0, o_reco_eta2 = -999.0, o_reco_phi2 = -999.0, o_reco_E2 = -999.0;
    float o_reco_dphi = -999.0, o_reco_xj = -999.0;
    int   o_reco_prov1 = -1, o_reco_prov2 = -1;
    int   o_reco_truth_idx1 = -1, o_reco_truth_idx2 = -1;
    int   o_legs_swapped = 0, o_n_accepted_reco = 0;

    if ( has_event_id ) tout -> Branch( "event_id", &o_event_id, "event_id/I" );
    if ( has_cent )     tout -> Branch( "cent", &o_cent, "cent/I" );
    if ( has_zvrtx )    tout -> Branch( "zvrtx", &o_zvrtx, "zvrtx/F" );
    if ( has_minbias )  tout -> Branch( "is_minbias", &o_is_minbias, "is_minbias/I" );
    if ( has_mbd_q )    tout -> Branch( "mbd_q", &o_mbd_q, "mbd_q/F" );
    if ( has_sumeT )    tout -> Branch( "sumeT", &o_sumeT, "sumeT/F" );
    if ( has_psi2 )
    {
        tout -> Branch( "psi2", &o_psi2, "psi2/F" );
        // dpsi2 = wrapped angle between the truth leading jet and the
        // 2nd-order event plane; see AnaUtils::get_dpsi2.
        tout -> Branch( "dpsi2", &o_dpsi2, "dpsi2/F" );
    }

    tout -> Branch( "category", &o_category, "category/I" );
    tout -> Branch( "truth_in_acc", &o_truth_in_acc, "truth_in_acc/I" );
    tout -> Branch( "reco_pair", &o_reco_pair, "reco_pair/I" );
    // bitmask: 1 = pT below threshold, 2 = outside truth_jet_accept_eta,
    // 4 = pair fails the dphi requirement (set on both legs together).
    tout -> Branch( "truth_fail1", &o_truth_fail1, "truth_fail1/I" );
    tout -> Branch( "truth_fail2", &o_truth_fail2, "truth_fail2/I" );
    tout -> Branch( "truth_idx1", &o_truth_idx1, "truth_idx1/I" );
    tout -> Branch( "truth_idx2", &o_truth_idx2, "truth_idx2/I" );
    tout -> Branch( "truth_pt1", &o_truth_pt1, "truth_pt1/F" );
    tout -> Branch( "truth_eta1", &o_truth_eta1, "truth_eta1/F" );
    tout -> Branch( "truth_phi1", &o_truth_phi1, "truth_phi1/F" );
    tout -> Branch( "truth_E1", &o_truth_E1, "truth_E1/F" );
    tout -> Branch( "truth_pt2", &o_truth_pt2, "truth_pt2/F" );
    tout -> Branch( "truth_eta2", &o_truth_eta2, "truth_eta2/F" );
    tout -> Branch( "truth_phi2", &o_truth_phi2, "truth_phi2/F" );
    tout -> Branch( "truth_E2", &o_truth_E2, "truth_E2/F" );
    tout -> Branch( "truth_dphi", &o_truth_dphi, "truth_dphi/F" );
    tout -> Branch( "truth_xj", &o_truth_xj, "truth_xj/F" );
    // each truth leg's own dR match, whether or not it ended up being one
    // of the selected reco dijet legs.
    tout -> Branch( "truth_match_idx1", &o_truth_match_idx1, "truth_match_idx1/I" );
    tout -> Branch( "truth_match_idx2", &o_truth_match_idx2, "truth_match_idx2/I" );
    tout -> Branch( "truth_match_pt1", &o_truth_match_pt1, "truth_match_pt1/F" );
    tout -> Branch( "truth_match_pt2", &o_truth_match_pt2, "truth_match_pt2/F" );
    tout -> Branch( "reco_idx1", &o_reco_idx1, "reco_idx1/I" );
    tout -> Branch( "reco_idx2", &o_reco_idx2, "reco_idx2/I" );
    tout -> Branch( "reco_pt1", &o_reco_pt1, "reco_pt1/F" );
    tout -> Branch( "reco_eta1", &o_reco_eta1, "reco_eta1/F" );
    tout -> Branch( "reco_phi1", &o_reco_phi1, "reco_phi1/F" );
    tout -> Branch( "reco_E1", &o_reco_E1, "reco_E1/F" );
    tout -> Branch( "reco_pt2", &o_reco_pt2, "reco_pt2/F" );
    tout -> Branch( "reco_eta2", &o_reco_eta2, "reco_eta2/F" );
    tout -> Branch( "reco_phi2", &o_reco_phi2, "reco_phi2/F" );
    tout -> Branch( "reco_E2", &o_reco_E2, "reco_E2/F" );
    tout -> Branch( "reco_dphi", &o_reco_dphi, "reco_dphi/F" );
    tout -> Branch( "reco_xj", &o_reco_xj, "reco_xj/F" );
    // -1 no leg, 0 truth leg 1's match, 1 truth leg 2's match,
    //  2 a softer truth jet outside the leading pair, 3 no truth match (UE).
    tout -> Branch( "reco_prov1", &o_reco_prov1, "reco_prov1/I" );
    tout -> Branch( "reco_prov2", &o_reco_prov2, "reco_prov2/I" );
    tout -> Branch( "reco_truth_idx1", &o_reco_truth_idx1, "reco_truth_idx1/I" );
    tout -> Branch( "reco_truth_idx2", &o_reco_truth_idx2, "reco_truth_idx2/I" );
    tout -> Branch( "legs_swapped", &o_legs_swapped, "legs_swapped/I" );
    tout -> Branch( "n_accepted_reco", &o_n_accepted_reco, "n_accepted_reco/I" );

    DijetPair::Config cfg;
    cfg.sublead_back_to_back = sublead_back_to_back;

    long n_events = 0;
    long n_cat[5] = { 0, 0, 0, 0, 0 };
    long n_swapped = 0;
    // kUESub breakdown: which leg(s) got taken, and by what
    long n_ue_lead_only = 0, n_ue_sub_only = 0, n_ue_both = 0;
    long n_ue_by_unmatched = 0, n_ue_by_other_truth = 0;
    // kFake breakdown by why the truth pair was out of acceptance
    long n_fake_truth_eta = 0, n_fake_truth_pt = 0, n_fake_truth_dphi = 0;
    long n_fake_legs_are_truth = 0, n_fake_legs_ue = 0;

    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );

        DijetPair::Result r;
        if ( !DijetPair::classify(
                 *truth_jet_pT, *truth_jet_phi, *truth_jet_accept_eta, *truth_jet_reco_match_idx,
                 *jet_pT, *jet_phi, *jet_E, *jet_accept_eta, *jet_truth_match_idx,
                 cfg, r ) )
        {
            continue; // no truth dijet candidate to report
        }

        ++n_events;
        ++n_cat[ r.category ];
        n_swapped += r.legs_swapped;

        if ( r.category == DijetPair::kUESub )
        {
            const bool lead_ok = ( r.reco_prov[0] == DijetPair::kProvLead || r.reco_prov[0] == DijetPair::kProvSub );
            const bool sub_ok  = ( r.reco_prov[1] == DijetPair::kProvLead || r.reco_prov[1] == DijetPair::kProvSub );
            if      ( !lead_ok && !sub_ok ) ++n_ue_both;
            else if ( !lead_ok )            ++n_ue_lead_only;
            else                            ++n_ue_sub_only;

            for ( int leg = 0; leg < 2; ++leg )
            {
                if ( r.reco_prov[leg] == DijetPair::kProvUE )         ++n_ue_by_unmatched;
                if ( r.reco_prov[leg] == DijetPair::kProvOtherTruth ) ++n_ue_by_other_truth;
            }
        }
        else if ( r.category == DijetPair::kFake )
        {
            const int f = r.truth_fail[0] | r.truth_fail[1];
            if ( f & DijetPair::kFailEta )  ++n_fake_truth_eta;
            if ( f & DijetPair::kFailPt )   ++n_fake_truth_pt;
            if ( f & DijetPair::kFailDphi ) ++n_fake_truth_dphi;

            const bool both_are_truth_legs =
                   ( r.reco_prov[0] == DijetPair::kProvLead || r.reco_prov[0] == DijetPair::kProvSub )
                && ( r.reco_prov[1] == DijetPair::kProvLead || r.reco_prov[1] == DijetPair::kProvSub );
            if ( both_are_truth_legs ) ++n_fake_legs_are_truth;
            if ( r.reco_prov[0] == DijetPair::kProvUE || r.reco_prov[1] == DijetPair::kProvUE ) ++n_fake_legs_ue;
        }

        o_event_id = event_id;
        o_cent = cent;
        o_zvrtx = zvrtx;
        o_is_minbias = is_minbias;
        o_mbd_q = mbd_q;
        o_sumeT = sumeT;
        o_psi2 = psi2;
        o_dpsi2 = has_psi2 ? AnaUtils::get_dpsi2( psi2, truth_jet_phi -> at( r.truth_idx[0] ) ) : -999.0f;

        o_category      = r.category;
        o_truth_in_acc  = r.truth_in_acc;
        o_reco_pair     = r.reco_pair;
        o_truth_fail1   = r.truth_fail[0];
        o_truth_fail2   = r.truth_fail[1];
        o_truth_idx1    = r.truth_idx[0];
        o_truth_idx2    = r.truth_idx[1];
        o_truth_dphi    = r.truth_dphi;

        o_truth_pt1  = truth_jet_pT  -> at( r.truth_idx[0] );
        o_truth_eta1 = truth_jet_eta -> at( r.truth_idx[0] );
        o_truth_phi1 = truth_jet_phi -> at( r.truth_idx[0] );
        o_truth_E1   = truth_jet_E   -> at( r.truth_idx[0] );
        o_truth_pt2  = truth_jet_pT  -> at( r.truth_idx[1] );
        o_truth_eta2 = truth_jet_eta -> at( r.truth_idx[1] );
        o_truth_phi2 = truth_jet_phi -> at( r.truth_idx[1] );
        o_truth_E2   = truth_jet_E   -> at( r.truth_idx[1] );
        o_truth_xj   = ( o_truth_pt1 > 0.0f ) ? o_truth_pt2 / o_truth_pt1 : -999.0f;

        o_truth_match_idx1 = r.truth_match_idx[0];
        o_truth_match_idx2 = r.truth_match_idx[1];
        o_truth_match_pt1  = r.truth_match_pt[0];
        o_truth_match_pt2  = r.truth_match_pt[1];

        o_reco_idx1 = r.reco_idx[0];
        o_reco_idx2 = r.reco_idx[1];
        o_reco_prov1 = r.reco_prov[0];
        o_reco_prov2 = r.reco_prov[1];
        o_reco_truth_idx1 = r.reco_truth_idx[0];
        o_reco_truth_idx2 = r.reco_truth_idx[1];
        o_reco_dphi = r.reco_dphi;
        o_legs_swapped = r.legs_swapped;
        o_n_accepted_reco = r.n_accepted_reco;

        if ( r.reco_pair )
        {
            o_reco_pt1  = jet_pT  -> at( r.reco_idx[0] );
            o_reco_eta1 = jet_eta -> at( r.reco_idx[0] );
            o_reco_phi1 = jet_phi -> at( r.reco_idx[0] );
            o_reco_E1   = jet_E   -> at( r.reco_idx[0] );
            o_reco_pt2  = jet_pT  -> at( r.reco_idx[1] );
            o_reco_eta2 = jet_eta -> at( r.reco_idx[1] );
            o_reco_phi2 = jet_phi -> at( r.reco_idx[1] );
            o_reco_E2   = jet_E   -> at( r.reco_idx[1] );
            o_reco_xj   = ( o_reco_pt1 > 0.0f ) ? o_reco_pt2 / o_reco_pt1 : -999.0f;
        }
        else
        {
            o_reco_pt1 = o_reco_eta1 = o_reco_phi1 = o_reco_E1 = -999.0;
            o_reco_pt2 = o_reco_eta2 = o_reco_phi2 = o_reco_E2 = -999.0;
            o_reco_xj = -999.0;
        }

        tout -> Fill();
    }

    std::cout << "\n---- pair-level dijet categories ----" << std::endl;
    std::cout << "Events with a truth dijet candidate: " << n_events << std::endl;
    std::cout << "  Fill  (0): " << n_cat[0] << "   <- the only category that fills the response" << std::endl;
    std::cout << "             of which reco legs pT-swapped: " << n_swapped << std::endl;
    std::cout << "  Miss  (1): " << n_cat[1] << std::endl;
    std::cout << "  Fake  (2): " << n_cat[2] << std::endl;
    std::cout << "             truth pair failed eta: " << n_fake_truth_eta
              << ", pT: " << n_fake_truth_pt
              << ", dphi: " << n_fake_truth_dphi << std::endl;
    std::cout << "             both reco legs are the truth pair: " << n_fake_legs_are_truth
              << ", >=1 leg is an unmatched/UE jet: " << n_fake_legs_ue << std::endl;
    std::cout << "  Skip  (3): " << n_cat[3] << std::endl;
    std::cout << "  UESub (4): " << n_cat[4] << "   <- excluded from Fill/Miss/Fake; covered by"
              << " the subleading-jet efficiency" << std::endl;
    std::cout << "             leading leg taken: " << n_ue_lead_only
              << ", subleading leg taken: " << n_ue_sub_only
              << ", both taken: " << n_ue_both << std::endl;
    std::cout << "             legs taken by an unmatched (UE) jet: " << n_ue_by_unmatched
              << ", by a softer truth jet: " << n_ue_by_other_truth << std::endl;

    fout -> cd();
    tout -> Write();
    const long n_written = tout -> GetEntries();
    fout -> Close();

    std::cout << "\nWrote " << n_written << " rows to " << outfile << std::endl;

    return 0;
}

#endif
