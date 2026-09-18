
#ifndef _DIJET_PAIR_MATCHING_INCLUSIVE_V2_C_
#define _DIJET_PAIR_MATCHING_INCLUSIVE_V2_C_

#include <myana/AnaUtils.h>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

// Inclusive-pairing variant of dijet_pair_matching_v2.C: instead of only
// the leading pair (pT1,pT2), the leading jet is paired with every other
// jet -- (pT1,pT2), (pT1,pT3), (pT1,pT4), ... -- and each pairing is
// categorized with exactly the pair-level scheme dijet_pair_matching_v2.C
// applies to the leading pair.
//
// This is to dijet_matching_inclusive.C what dijet_pair_matching_v2.C is
// to dijet_matching.C.
//
// THE SCHEME (identical to dijet_pair_matching_v2.C, applied at each rank)
//
// The pairing rank is the SAME on both sides. For rank `partner`:
//
//   truth pair = truth jets ( 0, partner )
//                pT-sorted on input, no eta cut at selection time
//   reco  pair = ( accepted[0], accepted[partner] )
//                where `accepted` is the event's ACCEPTED, TRUTH-MATCHED
//                reco jets in pT order -- jet_accept_eta, E > 0,
//                pT > 8 GeV, jet_truth_match_idx >= 0
//
// So both legs are rank-driven on the reco side, just as they are in the
// exclusive macro: nothing here is truth-driven, and a softer truth jet
// that outranks the expected leg genuinely takes the slot instead of
// being stepped over. partner==1 is the leading pair, which reproduces
// dijet_pair_matching_v2.C bit for bit; partner==2 is (pT1,pT3),
// partner==3 is (pT1,pT4), and so on.
//
// V2 -- the same two changes as dijet_pair_matching_v2.C:
//
//  (a) BOTH reco legs must be TRUTH-MATCHED (jet_truth_match_idx >= 0).
//      An unmatched UE fluctuation is stepped over when the accepted list
//      is built rather than allowed to take a leg, so kProvUE can never
//      appear and kUESub is never assigned -- the categories in practice
//      are Fill / Miss / Fake / Skip / FakeMiss.
//
//  (b) The reco dphi requirement no longer gates pair FORMATION. The pair
//      is built from the two matched legs, reco_dphi is recorded for every
//      pair, and dphi instead decides reco_in_acc -- whether the formed
//      pair is a reco dijet CANDIDATE. The category keys on that flag, not
//      on reco_pair.
//
//  (c) Five-way categorization, identical to dijet_pair_matching_v2.C:
//
//        truth cand | reco cand | reco IS the truth pair | category
//        -----------+-----------+------------------------+----------
//            no     |    no     |          --            | Skip
//            no     |    yes    |          --            | Fake
//            yes    |    no     |          --            | Miss
//            yes    |    yes    |          yes           | Fill
//            yes    |    yes    |          no            | FakeMiss
//
//      "reco IS the truth pair" means the two reco legs are the two truth
//      legs' own dR matches, in either order -- legs_swapped records a
//      resolution-driven ordering flip between two real legs, which is a
//      Fill, not a failure.
//
//      kFakeMiss is the population where both candidates are in
//      acceptance but the reco candidate is built from the wrong jets: it
//      is a fake AND a miss at once, and filling the response from it
//      would bias the subleading response. Either slot can be the wrong
//      one here, since both are rank-driven.
//
// Acceptance (leg 2 uses the subleading threshold set regardless of the
// pairing's actual rank in the event):
//   truth leg 1(2): pT > 14(7) GeV, truth_jet_accept_eta
//   reco  leg 1(2): pT > 20(8) GeV, jet_accept_eta, E > 0, truth-matched
//   the truth pair requires |dphi(1,j)| > 7pi/8; the reco pair's dphi
//   decides reco_in_acc rather than whether the pair is formed
//
// The truth pair is deliberately selected with NO eta cut: keeping it
// unconditional is what makes it possible to tell a reco jet that is
// unmatched because its truth partner was excluded by the calorimeter
// acceptance from a reco jet that is unmatched because it is a UE
// fluctuation. The per-leg truth_fail bitmask records WHY a truth leg
// failed, and reco_prov records WHAT each reco leg actually is.
//
// Output: one row per (event, pair_partner_rank), for every event with at
// least 2 truth jets. pair_partner_rank identifies the pairing (1, rank),
// so rank==2 is the leading dijet pair (pT1,pT2), rank==3 is (pT1,pT3),
// etc. -- the same convention as dijet_matching_inclusive.C.
//
// The namespace is deliberately distinct from dijet_pair_matching.C's
// DijetPair, dijet_pair_matching_v2.C's DijetPairV2 and
// dijet_pair_matching_inclusive.C's DijetPairInclusive, so that all of
// them can be loaded in the same ROOT session.
namespace DijetPairInclusiveV2
{
    // Same numbering as DijetPairV2 in dijet_pair_matching_v2.C.
    enum Category
    {
        kFill  = 0, // truth dijet candidate AND reco dijet candidate both
                    // in acceptance, and the reco candidate IS the truth
                    // candidate's own two matches (either order -- see
                    // legs_swapped). The only category that fills the
                    // response matrix.
        kMiss  = 1, // truth dijet candidate in acceptance; no reco dijet
                    // candidate in acceptance.
        kFake  = 2, // truth dijet candidate NOT in acceptance; reco dijet
                    // candidate in acceptance.
        kSkip  = 3, // neither candidate in acceptance.
        kUESub = 4, // never assigned in v2 (both legs are truth-matched by
                    // construction); kept so the category numbering stays
                    // comparable with the v1 macros.
        kFakeMiss = 5 // both candidates in acceptance, but the reco
                      // candidate is NOT the truth candidate's matches.
                      // Counts as BOTH a fake (that reco pair) and a miss
                      // (this truth pair) -- filling the response from it
                      // would bias the subleading response.
    };

    // What a selected reco leg actually is, from jet_truth_match_idx.
    // "Partner" here means the truth jet `partner` of the current pairing.
    enum Prov
    {
        kProvNone       = -1, // no reco leg
        kProvLead       =  0, // matched to truth jet 0 (the leading truth jet)
        kProvSub        =  1, // matched to truth jet `partner`
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
        // 1 once the formed pair also satisfies the reco dphi requirement
        // -- i.e. a reco dijet CANDIDATE exists, not merely two legs.
        int   reco_in_acc       = 0;
        int   reco_prov[2]      = { kProvNone, kProvNone };
        int   reco_truth_idx[2] = { -1, -1 };
        int   legs_swapped      = 0;
        int   n_accepted_reco   = 0;
    };

    // Builds the event-level list of ACCEPTED, TRUTH-MATCHED reco jets in
    // pT order. This is exactly the selection dijet_pair_matching_v2.C
    // applies inside its own classify(); it is lifted out here because
    // every pairing in the event reuses the same list, and because the
    // reco pairing rank indexes straight into it.
    //
    // jet_pT is pT-descending on input (sort_reco_jets_by_pt in
    // match_standalone.C), so `accepted` is pT-descending too: accepted[0]
    // is the leading accepted matched reco jet, accepted[k] the k-th.
    inline void build_accepted_reco(
        const std::vector< float > & jet_pT,
        const std::vector< float > & jet_E,
        const std::vector< int >   & jet_accept_eta,
        const std::vector< int >   & jet_truth_match_idx,
        const Config & cfg,
        std::vector< int > & accepted
    )
    {
        accepted.clear();
        for ( size_t j = 0; j < jet_pT.size(); ++j )
        {
            if (
                jet_accept_eta.at( j )
                && jet_E.at( j ) > 0.0f
                && jet_pT.at( j ) > cfg.reco_pt_thresh[1]
                && jet_truth_match_idx.at( j ) >= 0
            )
            {
                accepted.push_back( static_cast< int >( j ) );
            }
        }
    }

    // Classifies the pairing of rank `partner` -- truth jets (0, partner)
    // against reco jets (accepted[0], accepted[partner]). partner==1 is
    // the leading pair (pT1,pT2) and reproduces dijet_pair_matching_v2.C
    // exactly; partner==2 is (pT1,pT3), partner==3 is (pT1,pT4), and so
    // on. BOTH sides use the same rank, so the two legs are rank-driven on
    // the reco side just as they are in the exclusive macro -- nothing
    // here is truth-driven.
    //
    // `accepted` is the event-level list from build_accepted_reco(),
    // computed once by the caller and reused across every pairing.
    //
    // Returns false when this event cannot form the pairing at all (fewer
    // than partner+1 truth jets), in which case nothing should be written.
    // Pure otherwise: no ROOT I/O state, safe to lift into
    // match_standalone.C.
    inline bool classify(
        const int partner,
        const std::vector< int >   & accepted,
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
        if ( partner < 1 || partner >= static_cast< int >( truth_jet_pT.size() ) ) return false;

        r.n_accepted_reco = static_cast< int >( accepted.size() );

        //------------------------------------------------------------
        // truth pair (0, partner) -- NO eta cut at selection time, same
        // as dijet_pair_matching_v2.C: keeping the truth pair
        // unconditional is what lets truth_fail separate "reco jet is
        // unmatched because its truth partner was excluded" from "reco jet
        // is a UE jet".
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

        const bool has_truth_dijet_candidate = ( r.truth_in_acc == 1 );

        //------------------------------------------------------------
        // reco pair: built the way the data analysis builds it, from the
        // accepted reco jets in pT order -- NOT from the truth legs'
        // matches. The pairing rank is the SAME on both sides, so the
        // reco counterpart of the truth pairing (0, partner) is
        // (accepted[0], accepted[partner]).
        //------------------------------------------------------------
        int r1 = -1, r2 = -1;
        if ( !accepted.empty() && jet_pT.at( accepted[0] ) > cfg.reco_pt_thresh[0] )
        {
            r1 = accepted[0];
        }
        if ( static_cast< int >( accepted.size() ) > partner
             && jet_pT.at( accepted[partner] ) > cfg.reco_pt_thresh[1] )
        {
            r2 = accepted[partner];
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
                if      ( ti < 0 )        r.reco_prov[leg] = kProvUE;
                else if ( ti == 0 )       r.reco_prov[leg] = kProvLead;
                else if ( ti == partner ) r.reco_prov[leg] = kProvSub;
                else                      r.reco_prov[leg] = kProvOtherTruth;
            }
            r.reco_in_acc = ( r.reco_dphi >= cfg.min_dphi );
        }

        const bool has_reco_dijet_candidate  = ( r.reco_in_acc == 1 );
        const bool reco_is_matched_to_truth  = (r.reco_prov[0] == kProvLead && r.reco_prov[1] == kProvSub) || (r.reco_prov[0] == kProvSub && r.reco_prov[1] == kProvLead);

        //------------------------------------------------------------
        // pair-level category
        //------------------------------------------------------------

        if ( !has_truth_dijet_candidate && !has_reco_dijet_candidate ) r.category = kSkip;
        if ( !has_truth_dijet_candidate && has_reco_dijet_candidate ) r.category = kFake;
        if ( has_truth_dijet_candidate )
        {
            if ( !has_reco_dijet_candidate )
            {
                r.category = kMiss;
            }
            else
            {
                if ( reco_is_matched_to_truth )
                {
                    r.category = kFill;
                }
                else
                {
                    r.category = kFakeMiss;
                }
                r.legs_swapped = ( r.reco_prov[0] == kProvSub && r.reco_prov[1] == kProvLead ) ? 1 : 0;
            }

        }
        return true;
    }
}

int dijet_pair_matching_inclusive_v2(
    const std::string & infile = "output.root",
    const std::string & outfile = "dijet_pair_matching_inclusive_v2.root",
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
    int   o_category = -1, o_truth_in_acc = 0, o_reco_pair = 0, o_reco_in_acc = 0;
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
    // reco_pair alone is two legs above threshold; reco_in_acc is the
    // reco dijet CANDIDATE, i.e. those two legs plus the dphi
    // requirement. The category keys on reco_in_acc, so without it the
    // tree cannot reproduce its own categorization.
    tout -> Branch( "reco_in_acc", &o_reco_in_acc, "reco_in_acc/I" );
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
    // Both legs are rank-driven, so either can come back as 2; 3 cannot
    // appear, since both legs are required to be truth-matched.
    tout -> Branch( "reco_prov1", &o_reco_prov1, "reco_prov1/I" );
    tout -> Branch( "reco_prov2", &o_reco_prov2, "reco_prov2/I" );
    tout -> Branch( "reco_truth_idx1", &o_reco_truth_idx1, "reco_truth_idx1/I" );
    tout -> Branch( "reco_truth_idx2", &o_reco_truth_idx2, "reco_truth_idx2/I" );
    tout -> Branch( "legs_swapped", &o_legs_swapped, "legs_swapped/I" );
    tout -> Branch( "n_accepted_reco", &o_n_accepted_reco, "n_accepted_reco/I" );

    DijetPairInclusiveV2::Config cfg;

    long n_events = 0, n_rows = 0;
    // one slot per Category, kFakeMiss (5) included -- sizing this to the
    // number of categories is load-bearing, ++n_cat[r.category] indexes it
    // with the raw enum value.
    long n_cat[6] = { 0, 0, 0, 0, 0, 0 };
    // same, restricted to the leading pairing (rank 2) so it can be lined
    // up against dijet_pair_matching_v2.C on the same input
    long n_cat_r2[6] = { 0, 0, 0, 0, 0, 0 };
    long n_swapped = 0;
    long n_ue_by_unmatched = 0, n_ue_by_other_truth = 0;
    // kFakeMiss breakdown, the same one dijet_pair_matching_v2.C prints:
    // both legs are rank-driven here, so either can be the wrong one.
    long n_fm_lead_wrong = 0, n_fm_sub_wrong = 0, n_fm_both_wrong = 0;
    long n_fm_legs_other_truth = 0;
    long n_fm_truth_lead_unmatched = 0, n_fm_truth_sub_unmatched = 0;

    // reused across events so the per-event rebuild does not reallocate
    std::vector< int > accepted;

    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );
        if ( truth_jet_pT -> size() < 2 ) continue; // no truth pairing possible

        //------------------------------------------------------------
        // event-level reco jet list, computed once and reused across every
        // pairing in this event. Identical to the selection inside
        // dijet_pair_matching_v2.C's classify(), so the two macros'
        // n_accepted_reco mean the same thing; the pairing rank indexes
        // straight into it, accepted[0] being the leading accepted
        // matched reco jet.
        //------------------------------------------------------------
        DijetPairInclusiveV2::build_accepted_reco(
            *jet_pT, *jet_E, *jet_accept_eta, *jet_truth_match_idx, cfg, accepted );

        ++n_events;

        const int n_truth = static_cast< int >( truth_jet_pT -> size() );
        int last_partner = n_truth - 1;
        if ( max_partner_rank > 0 && ( max_partner_rank - 1 ) < last_partner )
        {
            last_partner = max_partner_rank - 1;
        }

        for ( int partner = 1; partner <= last_partner; ++partner )
        {
            DijetPairInclusiveV2::Result r;
            if ( !DijetPairInclusiveV2::classify(
                     partner, accepted,
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

            if ( r.category == DijetPairInclusiveV2::kUESub )
            {
                if ( r.reco_prov[0] == DijetPairInclusiveV2::kProvUE )         ++n_ue_by_unmatched;
                if ( r.reco_prov[0] == DijetPairInclusiveV2::kProvOtherTruth ) ++n_ue_by_other_truth;
            }
            else if ( r.category == DijetPairInclusiveV2::kFakeMiss )
            {
                // both candidates are in acceptance but the reco pair is
                // built from the wrong jets. Both legs are rank-driven, so
                // record which slot is wrong and what is sitting in it --
                // the same breakdown dijet_pair_matching_v2.C prints.
                const bool lead_ok = ( r.reco_prov[0] == DijetPairInclusiveV2::kProvLead
                                    || r.reco_prov[0] == DijetPairInclusiveV2::kProvSub );
                const bool sub_ok  = ( r.reco_prov[1] == DijetPairInclusiveV2::kProvLead
                                    || r.reco_prov[1] == DijetPairInclusiveV2::kProvSub );
                if      ( !lead_ok && !sub_ok ) ++n_fm_both_wrong;
                else if ( !lead_ok )            ++n_fm_lead_wrong;
                else                            ++n_fm_sub_wrong;

                for ( int leg = 0; leg < 2; ++leg )
                {
                    if ( r.reco_prov[leg] == DijetPairInclusiveV2::kProvOtherTruth ) ++n_fm_legs_other_truth;
                }

                // "no reco match at all" means the truth leg was never
                // reconstructed; otherwise it was simply outranked.
                if ( r.truth_match_idx[0] < 0 ) ++n_fm_truth_lead_unmatched;
                if ( r.truth_match_idx[1] < 0 ) ++n_fm_truth_sub_unmatched;
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
            o_reco_in_acc   = r.reco_in_acc;
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
              << ", UESub: " << n_cat[4]
              << ", FakeMiss: " << n_cat[5] << std::endl;
    std::cout << "  rank 2 only    -- Fill: " << n_cat_r2[0]
              << ", Miss: " << n_cat_r2[1]
              << ", Fake: " << n_cat_r2[2]
              << ", Skip: " << n_cat_r2[3]
              << ", UESub: " << n_cat_r2[4]
              << ", FakeMiss: " << n_cat_r2[5]
              << "   <- compare against dijet_pair_matching_v2.C" << std::endl;
    std::cout << "  Fill pairs with the reco legs pT-swapped: " << n_swapped << std::endl;
    std::cout << "  UESub legs taken by an unmatched (UE) jet: " << n_ue_by_unmatched
              << ", by a softer truth jet: " << n_ue_by_other_truth << std::endl;
    std::cout << "  FakeMiss wrong leg -- leading: " << n_fm_lead_wrong
              << ", subleading: " << n_fm_sub_wrong
              << ", both: " << n_fm_both_wrong << std::endl;
    std::cout << "           legs held by a softer truth jet: " << n_fm_legs_other_truth
              << " | truth leg with no reco match at all -- leading: " << n_fm_truth_lead_unmatched
              << ", partner: " << n_fm_truth_sub_unmatched << std::endl;
    std::cout << "  categorized: "
              << ( n_cat[0] + n_cat[1] + n_cat[2] + n_cat[3] + n_cat[4] + n_cat[5] )
              << " / " << n_rows << std::endl;
    std::cout << "  misses (Miss + FakeMiss): " << ( n_cat[1] + n_cat[5] )
              << " | fakes (Fake + FakeMiss): " << ( n_cat[2] + n_cat[5] ) << std::endl;
    std::cout << "  NOTE: both legs are rank-driven, on the truth side and on the"
              << " reco side, at the SAME pairing rank --" << std::endl;
    std::cout << "        rank 2 is therefore bit-identical to dijet_pair_matching_v2.C."
              << std::endl;

    fout -> cd();
    tout -> Write();
    const long n_written = tout -> GetEntries();
    fout -> Close();

    std::cout << "\nWrote " << n_written << " rows to " << outfile << std::endl;

    return 0;
}

#endif
