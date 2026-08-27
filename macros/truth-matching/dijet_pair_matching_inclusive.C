
#ifndef _DIJET_PAIR_MATCHING_INCLUSIVE_C_
#define _DIJET_PAIR_MATCHING_INCLUSIVE_C_

#include <myana/AnaUtils.h>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

// Inclusive-pairing variant of dijet_pair_matching.C: instead of only the
// leading truth pair (1,2), the leading truth jet is paired with every
// other truth jet in the event -- (1,2), (1,3), (1,4), ... -- and each
// pairing is categorized with the same pair-level scheme, including the
// kUESub category and the truth_fail / reco_prov bookkeeping.
//
// This is to dijet_matching_inclusive.C what dijet_pair_matching.C is to
// dijet_matching.C. The two additions over dijet_matching_inclusive.C are
// that the LEADING reco leg is now taken from the independently-ranked
// accepted reco population (so a UE fluctuation that has taken the
// leading slot is visible at all), and that such pairings get their own
// category instead of being absorbed into Fill.
//
// RECO SIDE -- per-pairing construction (note the asymmetry):
//   leg 1 = the leading ACCEPTED reco jet in the event (pT > 20,
//           jet_accept_eta, E > 0). Independently ranked, so its
//           provenance is a real measurement: reco_prov1 says whether the
//           leading reco slot is actually the leading truth jet's match, a
//           softer truth jet, or an unmatched UE jet.
//   leg 2 = truth partner j's OWN dR match (truth_jet_reco_match_idx[j]),
//           required to pass pT > 8, jet_accept_eta, E > 0.
//
// Because leg 2 is truth-driven rather than rank-driven, this scheme can
// only ever detect substitution of the LEADING leg. A UE jet taking the
// subleading slot is invisible here by construction -- use
// dijet_pair_matching.C (where both legs are rank-driven) for that.
// The real invariant is that leg 2 is ALWAYS one of the two truth legs'
// own matches -- kProvSub for an ordinary pair, kProvLead for a swapped
// one -- and never kProvUE or kProvOtherTruth. That is exactly why
// subleading substitution cannot show up here.
//
// SWAPPED PAIRS: if the leading accepted reco jet IS the partner's own
// match, the pair is not degenerate -- it is inverted. It is then taken as
// leg 1 = partner's match, leg 2 = the leading truth jet's match, with
// legs_swapped = 1, exactly as dijet_pair_matching.C treats a
// resolution-driven ordering flip between two real legs. If the leading
// truth jet has no match at all in that situation there is no second leg
// and no reco pair is formed.
//
// Acceptance (unchanged; leg 2 uses the subleading threshold set
// regardless of the partner's actual rank in the event):
//   truth leg 1(2): pT > 14(7) GeV, truth_jet_accept_eta
//   reco  leg 1(2): pT > 20(8) GeV, jet_accept_eta, E > 0
//   each truth pair additionally requires |dphi(1,j)| > 7pi/8
//
// Output: one row per (event, pair_partner_rank), for every event with at
// least 2 truth jets. pair_partner_rank identifies the pairing (1, rank),
// so rank==2 is the leading dijet pair (1,2), rank==3 is (1,3), etc. --
// the same convention as dijet_matching_inclusive.C.
//
// The namespace is deliberately distinct from dijet_pair_matching.C's
// DijetPair so that both macros can be loaded in the same ROOT session.
namespace DijetPairInclusive
{
    // Same numbering as DijetPair in dijet_pair_matching.C.
    enum Category
    {
        kFill  = 0, // truth pair in acceptance; reco pair exists and its
                    // two legs ARE the two truth legs' own matches.
        kMiss  = 1, // truth pair in acceptance; no reco pair.
        kFake  = 2, // truth pair NOT in acceptance; reco pair exists.
        kSkip  = 3, // truth pair NOT in acceptance; no reco pair.
        kUESub = 4  // truth pair in acceptance; reco pair exists, but the
                    // leading reco leg is not the leading truth jet's
                    // match -- something else holds the leading slot.
    };

    // What a selected reco leg actually is, from jet_truth_match_idx.
    // "Partner" here means the truth jet j of the current pairing.
    enum Prov
    {
        kProvNone       = -1, // no reco leg
        kProvLead       =  0, // matched to truth jet 0 (the leading truth jet)
        kProvSub        =  1, // matched to truth jet j (this pairing's partner)
        kProvOtherTruth =  2, // matched to some other truth jet
        kProvUE         =  3  // no truth match at all -> UE fluctuation
    };

    // Why a truth leg is out of acceptance; zero means it is in.
    // kFailDphi is a pair-level condition, set on both legs together.
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
    };

    struct Result
    {
        int   category      = kSkip;
        int   truth_in_acc  = 0;
        int   truth_fail[2] = { kFailNone, kFailNone };
        int   truth_idx[2]  = { -1, -1 };
        float truth_dphi    = -999.0f;
        int   truth_match_idx[2] = { -1, -1 };
        float truth_match_pt[2]  = { -999.0f, -999.0f };

        int   reco_pair         = 0;
        int   reco_idx[2]       = { -1, -1 };
        float reco_dphi         = -999.0f;
        int   reco_prov[2]      = { kProvNone, kProvNone };
        int   reco_truth_idx[2] = { -1, -1 };
        int   legs_swapped      = 0;
        int   n_accepted_reco   = 0;
    };

    // Classifies the pairing (truth jet 0, truth jet partner). `r_lead` is
    // the event's leading accepted reco jet index (or -1 if the event has
    // none above the leading threshold) and `n_accepted_reco` its
    // multiplicity above the subleading floor -- both are event-level, so
    // the caller computes them once and reuses them across partners.
    //
    // Pure: no ROOT I/O state, safe to lift into match_standalone.C.
    inline bool classify(
        const int partner,
        const int r_lead,
        const int n_accepted_reco,
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
        if ( partner < 1 || partner >= static_cast< int >( truth_jet_pT.size() ) ) return false;

        r.n_accepted_reco = n_accepted_reco;

        //------------------------------------------------------------
        // truth pair (0, partner) -- NO eta cut at selection time, same
        // as dijet_pair_matching.C: keeping the truth pair unconditional
        // is what lets truth_fail separate "reco jet is unmatched because
        // its truth partner was excluded" from "reco jet is a UE jet".
        //------------------------------------------------------------
        r.truth_idx[0] = 0;
        r.truth_idx[1] = partner;
        r.truth_dphi = AnaUtils::dphi_wrap( truth_jet_phi.at( 0 ), truth_jet_phi.at( partner ) );
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
        // reco pair: leg 1 independently ranked, leg 2 truth-driven.
        //------------------------------------------------------------
        const int m0 = r.truth_match_idx[0];       // leading truth jet's match
        const int mj = r.truth_match_idx[1];       // partner's match

        int leg1 = -1, leg2 = -1;
        if ( r_lead >= 0 )
        {
            if ( r_lead == mj )
            {
                // the partner's own match holds the leading reco slot --
                // an inverted, not a degenerate, pair.
                leg1 = mj;
                leg2 = m0;              // -1 if the leading truth jet has no match
                r.legs_swapped = ( m0 >= 0 ) ? 1 : 0;
            }
            else
            {
                leg1 = r_lead;
                leg2 = mj;
            }
        }

        // leg1 already satisfies the leading requirements by construction
        // (it is the leading accepted reco jet above reco_pt_thresh[0]);
        // leg2 still has to be checked against the subleading set.
        const bool leg2_ok =
               ( leg2 >= 0 )
            && ( jet_pT.at( leg2 ) > cfg.reco_pt_thresh[1] )
            && jet_accept_eta.at( leg2 )
            && ( jet_E.at( leg2 ) > 0.0f );

        if ( leg1 >= 0 && leg2_ok
             && AnaUtils::dphi_wrap( jet_phi.at( leg1 ), jet_phi.at( leg2 ) ) > cfg.min_dphi )
        {
            r.reco_pair = 1;
            r.reco_idx[0] = leg1;
            r.reco_idx[1] = leg2;
            r.reco_dphi = AnaUtils::dphi_wrap( jet_phi.at( leg1 ), jet_phi.at( leg2 ) );

            for ( int leg = 0; leg < 2; ++leg )
            {
                const int ri = r.reco_idx[leg];
                const int ti = jet_truth_match_idx.at( ri );
                r.reco_truth_idx[leg] = ti;
                if      ( ti < 0 )        r.reco_prov[leg] = kProvUE;
                else if ( ti == 0 )       r.reco_prov[leg] = kProvLead;
                else if ( ti == partner ) r.reco_prov[leg] = kProvSub;
                else                      r.reco_prov[leg] = kProvOtherTruth;
            }
        }
        else
        {
            r.legs_swapped = 0; // no pair formed, so nothing was swapped
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
                // both reco legs must be the two truth legs' own matches,
                // in either order (legs_swapped records the order).
                const bool same_pair =
                       ( r.reco_prov[0] == kProvLead && r.reco_prov[1] == kProvSub )
                    || ( r.reco_prov[0] == kProvSub  && r.reco_prov[1] == kProvLead );
                r.category = same_pair ? kFill : kUESub;
            }
        }
        else
        {
            r.category = r.reco_pair ? kFake : kSkip;
        }

        return true;
    }
}

int dijet_pair_matching_inclusive(
    const std::string & infile = "output.root",
    const std::string & outfile = "dijet_pair_matching_inclusive.root",
    const int max_partner_rank = -1   // -1 = every truth jet; else cap the
                                      // pairing at (1, max_partner_rank)
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
    auto * tout = new TTree( "T", "inclusive-pairing dijet matching / response" );

    int   o_event_id = -1, o_cent = -1, o_is_minbias = 0;
    float o_zvrtx = 0.0, o_mbd_q = -999.0, o_sumeT = -999.0;
    float o_psi2 = -999.0, o_dpsi2 = -999.0;

    int   o_pair_partner_rank = -1, o_n_truth_jets = 0;
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
        tout -> Branch( "dpsi2", &o_dpsi2, "dpsi2/F" );
    }

    // pair_partner_rank identifies the pairing (1, pair_partner_rank), so
    // pair_partner_rank==2 is the leading dijet pair (1,2), ==3 is (1,3),
    // etc. -- same convention as dijet_matching_inclusive.C.
    tout -> Branch( "pair_partner_rank", &o_pair_partner_rank, "pair_partner_rank/I" );
    tout -> Branch( "n_truth_jets", &o_n_truth_jets, "n_truth_jets/I" );
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
    // -1 no leg, 0 leading truth jet's match, 1 this pairing's partner's
    //  match, 2 some other truth jet, 3 no truth match (UE).
    // leg 2 is truth-driven, so reco_prov2 is always 1 (ordinary pair) or
    // 0 (swapped pair) when a pair exists -- never 2 or 3.
    tout -> Branch( "reco_prov1", &o_reco_prov1, "reco_prov1/I" );
    tout -> Branch( "reco_prov2", &o_reco_prov2, "reco_prov2/I" );
    tout -> Branch( "reco_truth_idx1", &o_reco_truth_idx1, "reco_truth_idx1/I" );
    tout -> Branch( "reco_truth_idx2", &o_reco_truth_idx2, "reco_truth_idx2/I" );
    tout -> Branch( "legs_swapped", &o_legs_swapped, "legs_swapped/I" );
    tout -> Branch( "n_accepted_reco", &o_n_accepted_reco, "n_accepted_reco/I" );

    DijetPairInclusive::Config cfg;

    long n_events = 0, n_rows = 0;
    long n_cat[5] = { 0, 0, 0, 0, 0 };
    // same, restricted to the leading pairing (rank 2) so it can be lined
    // up against dijet_pair_matching.C on the same input
    long n_cat_r2[5] = { 0, 0, 0, 0, 0 };
    long n_swapped = 0;
    long n_ue_by_unmatched = 0, n_ue_by_other_truth = 0;

    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );
        if ( truth_jet_pT -> size() < 2 ) continue; // no truth pairing possible

        //------------------------------------------------------------
        // event-level reco quantities, computed once and reused across
        // every pairing in this event.
        //------------------------------------------------------------
        int r_lead = -1;
        int n_accepted_reco = 0;
        for ( size_t j = 0; j < jet_pT -> size(); ++j )
        {
            if ( !( jet_accept_eta -> at( j ) && jet_E -> at( j ) > 0.0f ) ) continue;
            // jet_pT is pT-descending on input (sort_reco_jets_by_pt in
            // match_standalone.C), so the first accepted jet is the
            // leading one.
            if ( r_lead < 0 && jet_pT -> at( j ) > cfg.reco_pt_thresh[0] )
            {
                r_lead = static_cast< int >( j );
            }
            if ( jet_pT -> at( j ) > cfg.reco_pt_thresh[1] ) ++n_accepted_reco;
        }

        ++n_events;

        const int n_truth = static_cast< int >( truth_jet_pT -> size() );
        int last_partner = n_truth - 1;
        if ( max_partner_rank > 0 && ( max_partner_rank - 1 ) < last_partner )
        {
            last_partner = max_partner_rank - 1;
        }

        for ( int partner = 1; partner <= last_partner; ++partner )
        {
            DijetPairInclusive::Result r;
            if ( !DijetPairInclusive::classify(
                     partner, r_lead, n_accepted_reco,
                     *truth_jet_pT, *truth_jet_phi, *truth_jet_accept_eta, *truth_jet_reco_match_idx,
                     *jet_pT, *jet_phi, *jet_E, *jet_accept_eta, *jet_truth_match_idx,
                     cfg, r ) )
            {
                continue;
            }

            ++n_rows;
            ++n_cat[ r.category ];
            if ( partner == 1 ) ++n_cat_r2[ r.category ];
            n_swapped += r.legs_swapped;

            if ( r.category == DijetPairInclusive::kUESub )
            {
                if ( r.reco_prov[0] == DijetPairInclusive::kProvUE )         ++n_ue_by_unmatched;
                if ( r.reco_prov[0] == DijetPairInclusive::kProvOtherTruth ) ++n_ue_by_other_truth;
            }

            o_event_id = event_id;
            o_cent = cent;
            o_zvrtx = zvrtx;
            o_is_minbias = is_minbias;
            o_mbd_q = mbd_q;
            o_sumeT = sumeT;
            o_psi2 = psi2;
            o_dpsi2 = has_psi2 ? AnaUtils::get_dpsi2( psi2, truth_jet_phi -> at( 0 ) ) : -999.0f;

            o_pair_partner_rank = partner + 1;
            o_n_truth_jets  = n_truth;
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
    }

    std::cout << "\n---- inclusive-pairing dijet categories ----" << std::endl;
    std::cout << "Events with >= 2 truth jets: " << n_events
              << "   pairings written: " << n_rows << std::endl;
    std::cout << "  all pairings   -- Fill: " << n_cat[0]
              << ", Miss: " << n_cat[1]
              << ", Fake: " << n_cat[2]
              << ", Skip: " << n_cat[3]
              << ", UESub: " << n_cat[4] << std::endl;
    std::cout << "  rank 2 only    -- Fill: " << n_cat_r2[0]
              << ", Miss: " << n_cat_r2[1]
              << ", Fake: " << n_cat_r2[2]
              << ", Skip: " << n_cat_r2[3]
              << ", UESub: " << n_cat_r2[4]
              << "   <- compare against dijet_pair_matching.C" << std::endl;
    std::cout << "  swapped pairs (partner's match held the leading slot): " << n_swapped << std::endl;
    std::cout << "  UESub leading slot held by an unmatched (UE) jet: " << n_ue_by_unmatched
              << ", by a softer truth jet: " << n_ue_by_other_truth << std::endl;
    std::cout << "  NOTE: leg 2 is truth-driven here, so subleading-slot"
              << " substitution is invisible by construction --" << std::endl;
    std::cout << "        use dijet_pair_matching.C for that." << std::endl;

    fout -> cd();
    tout -> Write();
    const long n_written = tout -> GetEntries();
    fout -> Close();

    std::cout << "\nWrote " << n_written << " rows to " << outfile << std::endl;

    return 0;
}

#endif
