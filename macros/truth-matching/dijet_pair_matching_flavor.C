
#ifndef _DIJET_PAIR_MATCHING_FLAVOR_C_
#define _DIJET_PAIR_MATCHING_FLAVOR_C_

// Pulls in the DijetPair namespace (Category / Prov / TruthFail / Config /
// Result / classify) so the classification lives in exactly one place --
// edit dijet_pair_matching.C and this macro follows automatically. The
// include guard there keeps the extra dijet_pair_matching() entry point
// harmless if both macros are loaded in the same ROOT session.
#include "dijet_pair_matching.C"

#include <myana/AnaUtils.h>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

// Flavor-split version of dijet_pair_matching.C: the exclusive leading
// truth pair (truth jets 0 and 1) is categorized with the pair-level
// Fill / Miss / Fake / Skip / UESub scheme, then routed to one of two
// output files according to the hard-parton flavor of the two legs.
//
// This is to dijet_pair_matching.C what dijet_matching_flavor.C is to
// dijet_matching.C. The flavor tagging is unchanged from
// dijet_matching_flavor.C (truth_jet_flavor = PDG id of the matched
// parton: 1-6 quark, 21 gluon, -1 unmatched -- the leading_dijet_flavor.C
// convention); events whose leading dijet is unmatched or "other" go to
// neither output.
//
// What IS new relative to dijet_matching_flavor.C:
//   - one row per event (pair-level) instead of one row per leg, since a
//     dijet response is a pair-level object;
//   - the reco dijet is the leading/subleading ACCEPTED reco jets, so a UE
//     fluctuation that has taken a leg is visible;
//   - kUESub keeps those pairs out of Fill, Miss AND Fake;
//   - truth_fail / reco_prov bookkeeping is written out.
// See dijet_pair_matching.C for the full description of all of the above.
//
// Flavor is taken from the two leading TRUTH jets, which is the same pair
// the classification uses -- so the flavor tag is well defined for every
// row including Fake and Skip, where the reco legs may be different jets
// entirely. reco_prov1/2 tell you what the reco legs actually were.
int dijet_pair_matching_flavor(
    const std::string & infile = "output.root",
    const std::string & outfile_qq    = "dijet_pair_matching_qq.root",
    const std::string & outfile_qg_gg = "dijet_pair_matching_qg_gg.root",
    const bool sublead_back_to_back = false
)
{
    // same convention as leading_dijet_flavor.C / dijet_matching_flavor.C
    auto is_qq_qg_gg = []( const int flavor1, const int flavor2, bool & is_qq ) -> bool
    {
        const bool isq1 = ( flavor1 >= 1 && flavor1 <= 6 );
        const bool isq2 = ( flavor2 >= 1 && flavor2 <= 6 );
        const bool isg1 = ( flavor1 == 21 );
        const bool isg2 = ( flavor2 == 21 );

        if ( isq1 && isq2 ) { is_qq = true;  return true; }
        if ( ( isq1 && isg2 ) || ( isg1 && isq2 ) ) { is_qq = false; return true; }
        if ( isg1 && isg2 ) { is_qq = false; return true; }
        return false; // unmatched (flavor<=0) or "other" -- neither qq nor qg/gg
    };

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
    // truth / reco jets, flavor, and the match_standalone.C matching
    //----------------------------------------------------------------
    std::vector< float > * truth_jet_pT  = nullptr;
    std::vector< float > * truth_jet_eta = nullptr;
    std::vector< float > * truth_jet_phi = nullptr;
    std::vector< float > * truth_jet_E   = nullptr;
    std::vector< int >   * truth_jet_accept_eta = nullptr;
    std::vector< int >   * truth_jet_flavor = nullptr;
    std::vector< int >   * truth_jet_reco_match_idx = nullptr;

    std::vector< float > * jet_pT  = nullptr;
    std::vector< float > * jet_eta = nullptr;
    std::vector< float > * jet_phi = nullptr;
    std::vector< float > * jet_E   = nullptr;
    std::vector< int >   * jet_accept_eta = nullptr;
    std::vector< int >   * jet_truth_match_idx = nullptr;

    const bool has_all = has_branch( "truth_jet_pT" )
                        && has_branch( "truth_jet_accept_eta" )
                        && has_branch( "truth_jet_flavor" )
                        && has_branch( "truth_jet_reco_match_idx" )
                        && has_branch( "jet_pT" )
                        && has_branch( "jet_accept_eta" )
                        && has_branch( "jet_truth_match_idx" );
    if ( !has_all )
    {
        std::cerr << "Error: input tree is missing the truth/reco jet, flavor, or "
                  << "matching branches written by match_standalone.C "
                  << "(truth_jet_flavor / truth_jet_reco_match_idx / jet_truth_match_idx)." << std::endl;
        return -1;
    }

    enable( "truth_jet_pT", &truth_jet_pT );
    enable( "truth_jet_eta", &truth_jet_eta );
    enable( "truth_jet_phi", &truth_jet_phi );
    enable( "truth_jet_E", &truth_jet_E );
    enable( "truth_jet_accept_eta", &truth_jet_accept_eta );
    enable( "truth_jet_flavor", &truth_jet_flavor );
    enable( "truth_jet_reco_match_idx", &truth_jet_reco_match_idx );
    enable( "jet_pT", &jet_pT );
    enable( "jet_eta", &jet_eta );
    enable( "jet_phi", &jet_phi );
    enable( "jet_E", &jet_E );
    enable( "jet_accept_eta", &jet_accept_eta );
    enable( "jet_truth_match_idx", &jet_truth_match_idx );

    //----------------------------------------------------------------
    // one identically-schemed pair-level output tree per flavor category
    //----------------------------------------------------------------
    struct Output
    {
        TFile * file = nullptr;
        TTree * tree = nullptr;

        int   event_id = -1, cent = -1, is_minbias = 0;
        float zvrtx = 0.0, mbd_q = -999.0, sumeT = -999.0;
        float psi2 = -999.0, dpsi2 = -999.0;

        int   category = -1, truth_in_acc = 0, reco_pair = 0;
        int   truth_fail1 = 0, truth_fail2 = 0;
        int   truth_idx1 = -1, truth_idx2 = -1;
        float truth_pt1 = -999.0, truth_eta1 = -999.0, truth_phi1 = -999.0, truth_E1 = -999.0;
        float truth_pt2 = -999.0, truth_eta2 = -999.0, truth_phi2 = -999.0, truth_E2 = -999.0;
        float truth_dphi = -999.0, truth_xj = -999.0;
        int   truth_flavor1 = -1, truth_flavor2 = -1, is_qq = 0;
        int   truth_match_idx1 = -1, truth_match_idx2 = -1;
        float truth_match_pt1 = -999.0, truth_match_pt2 = -999.0;
        int   reco_idx1 = -1, reco_idx2 = -1;
        float reco_pt1 = -999.0, reco_eta1 = -999.0, reco_phi1 = -999.0, reco_E1 = -999.0;
        float reco_pt2 = -999.0, reco_eta2 = -999.0, reco_phi2 = -999.0, reco_E2 = -999.0;
        float reco_dphi = -999.0, reco_xj = -999.0;
        int   reco_prov1 = -1, reco_prov2 = -1;
        int   reco_truth_idx1 = -1, reco_truth_idx2 = -1;
        int   legs_swapped = 0, n_accepted_reco = 0;

        long  n_events = 0;
        long  n_cat[5] = { 0, 0, 0, 0, 0 };
        long  n_swapped = 0;
        long  n_ue_by_unmatched = 0, n_ue_by_other_truth = 0;
    };

    // Books into an already-constructed Output, rather than returning one
    // by value -- the branch addresses have to point at the object that
    // actually outlives this call.
    auto book_output = [&]( Output & out, const std::string & filename )
    {
        out.file = new TFile( filename.c_str(), "RECREATE" );
        out.tree = new TTree( "T", "pair-level dijet matching / response (flavor-tagged)" );

        if ( has_event_id ) out.tree -> Branch( "event_id", &out.event_id, "event_id/I" );
        if ( has_cent )     out.tree -> Branch( "cent", &out.cent, "cent/I" );
        if ( has_zvrtx )    out.tree -> Branch( "zvrtx", &out.zvrtx, "zvrtx/F" );
        if ( has_minbias )  out.tree -> Branch( "is_minbias", &out.is_minbias, "is_minbias/I" );
        if ( has_mbd_q )    out.tree -> Branch( "mbd_q", &out.mbd_q, "mbd_q/F" );
        if ( has_sumeT )    out.tree -> Branch( "sumeT", &out.sumeT, "sumeT/F" );
        if ( has_psi2 )
        {
            out.tree -> Branch( "psi2", &out.psi2, "psi2/F" );
            out.tree -> Branch( "dpsi2", &out.dpsi2, "dpsi2/F" );
        }

        out.tree -> Branch( "category", &out.category, "category/I" );
        out.tree -> Branch( "truth_in_acc", &out.truth_in_acc, "truth_in_acc/I" );
        out.tree -> Branch( "reco_pair", &out.reco_pair, "reco_pair/I" );
        // bitmask: 1 = pT below threshold, 2 = outside truth_jet_accept_eta,
        // 4 = pair fails the dphi requirement (set on both legs together).
        out.tree -> Branch( "truth_fail1", &out.truth_fail1, "truth_fail1/I" );
        out.tree -> Branch( "truth_fail2", &out.truth_fail2, "truth_fail2/I" );
        out.tree -> Branch( "truth_idx1", &out.truth_idx1, "truth_idx1/I" );
        out.tree -> Branch( "truth_idx2", &out.truth_idx2, "truth_idx2/I" );
        out.tree -> Branch( "truth_pt1", &out.truth_pt1, "truth_pt1/F" );
        out.tree -> Branch( "truth_eta1", &out.truth_eta1, "truth_eta1/F" );
        out.tree -> Branch( "truth_phi1", &out.truth_phi1, "truth_phi1/F" );
        out.tree -> Branch( "truth_E1", &out.truth_E1, "truth_E1/F" );
        out.tree -> Branch( "truth_pt2", &out.truth_pt2, "truth_pt2/F" );
        out.tree -> Branch( "truth_eta2", &out.truth_eta2, "truth_eta2/F" );
        out.tree -> Branch( "truth_phi2", &out.truth_phi2, "truth_phi2/F" );
        out.tree -> Branch( "truth_E2", &out.truth_E2, "truth_E2/F" );
        out.tree -> Branch( "truth_dphi", &out.truth_dphi, "truth_dphi/F" );
        out.tree -> Branch( "truth_xj", &out.truth_xj, "truth_xj/F" );
        // PDG id of each truth leg's matched hard parton: 1-6 quark,
        // 21 gluon. is_qq is 1 for the qq output, 0 for the qg/gg one.
        out.tree -> Branch( "truth_flavor1", &out.truth_flavor1, "truth_flavor1/I" );
        out.tree -> Branch( "truth_flavor2", &out.truth_flavor2, "truth_flavor2/I" );
        out.tree -> Branch( "is_qq", &out.is_qq, "is_qq/I" );
        out.tree -> Branch( "truth_match_idx1", &out.truth_match_idx1, "truth_match_idx1/I" );
        out.tree -> Branch( "truth_match_idx2", &out.truth_match_idx2, "truth_match_idx2/I" );
        out.tree -> Branch( "truth_match_pt1", &out.truth_match_pt1, "truth_match_pt1/F" );
        out.tree -> Branch( "truth_match_pt2", &out.truth_match_pt2, "truth_match_pt2/F" );
        out.tree -> Branch( "reco_idx1", &out.reco_idx1, "reco_idx1/I" );
        out.tree -> Branch( "reco_idx2", &out.reco_idx2, "reco_idx2/I" );
        out.tree -> Branch( "reco_pt1", &out.reco_pt1, "reco_pt1/F" );
        out.tree -> Branch( "reco_eta1", &out.reco_eta1, "reco_eta1/F" );
        out.tree -> Branch( "reco_phi1", &out.reco_phi1, "reco_phi1/F" );
        out.tree -> Branch( "reco_E1", &out.reco_E1, "reco_E1/F" );
        out.tree -> Branch( "reco_pt2", &out.reco_pt2, "reco_pt2/F" );
        out.tree -> Branch( "reco_eta2", &out.reco_eta2, "reco_eta2/F" );
        out.tree -> Branch( "reco_phi2", &out.reco_phi2, "reco_phi2/F" );
        out.tree -> Branch( "reco_E2", &out.reco_E2, "reco_E2/F" );
        out.tree -> Branch( "reco_dphi", &out.reco_dphi, "reco_dphi/F" );
        out.tree -> Branch( "reco_xj", &out.reco_xj, "reco_xj/F" );
        // -1 no leg, 0 truth leg 1's match, 1 truth leg 2's match,
        //  2 a softer truth jet outside the leading pair, 3 no truth match (UE).
        out.tree -> Branch( "reco_prov1", &out.reco_prov1, "reco_prov1/I" );
        out.tree -> Branch( "reco_prov2", &out.reco_prov2, "reco_prov2/I" );
        out.tree -> Branch( "reco_truth_idx1", &out.reco_truth_idx1, "reco_truth_idx1/I" );
        out.tree -> Branch( "reco_truth_idx2", &out.reco_truth_idx2, "reco_truth_idx2/I" );
        out.tree -> Branch( "legs_swapped", &out.legs_swapped, "legs_swapped/I" );
        out.tree -> Branch( "n_accepted_reco", &out.n_accepted_reco, "n_accepted_reco/I" );
    };

    Output out_qq, out_qg_gg;
    book_output( out_qq, outfile_qq );
    book_output( out_qg_gg, outfile_qg_gg );

    DijetPair::Config cfg;
    cfg.sublead_back_to_back = sublead_back_to_back;

    long n_events_unmatched_flavor = 0;

    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );

        // truth dijet candidate = leading two truth jets overall, no
        // acceptance cut at selection time. truth jets are pT-sorted on
        // input, so this is jets 0 and 1.
        if ( truth_jet_pT -> size() < 2 ) continue;

        bool is_qq = false;
        const int flavor1 = truth_jet_flavor -> at( 0 );
        const int flavor2 = truth_jet_flavor -> at( 1 );
        if ( !is_qq_qg_gg( flavor1, flavor2, is_qq ) )
        {
            ++n_events_unmatched_flavor; // unmatched/"other" -- neither output
            continue;
        }
        Output & out = is_qq ? out_qq : out_qg_gg;

        DijetPair::Result r;
        if ( !DijetPair::classify(
                 *truth_jet_pT, *truth_jet_phi, *truth_jet_accept_eta, *truth_jet_reco_match_idx,
                 *jet_pT, *jet_phi, *jet_E, *jet_accept_eta, *jet_truth_match_idx,
                 cfg, r ) )
        {
            continue;
        }

        ++out.n_events;
        ++out.n_cat[ r.category ];
        out.n_swapped += r.legs_swapped;

        if ( r.category == DijetPair::kUESub )
        {
            for ( int leg = 0; leg < 2; ++leg )
            {
                if ( r.reco_prov[leg] == DijetPair::kProvUE )         ++out.n_ue_by_unmatched;
                if ( r.reco_prov[leg] == DijetPair::kProvOtherTruth ) ++out.n_ue_by_other_truth;
            }
        }

        out.event_id = event_id;
        out.cent = cent;
        out.zvrtx = zvrtx;
        out.is_minbias = is_minbias;
        out.mbd_q = mbd_q;
        out.sumeT = sumeT;
        out.psi2 = psi2;
        out.dpsi2 = has_psi2 ? AnaUtils::get_dpsi2( psi2, truth_jet_phi -> at( r.truth_idx[0] ) ) : -999.0f;

        out.category      = r.category;
        out.truth_in_acc  = r.truth_in_acc;
        out.reco_pair     = r.reco_pair;
        out.truth_fail1   = r.truth_fail[0];
        out.truth_fail2   = r.truth_fail[1];
        out.truth_idx1    = r.truth_idx[0];
        out.truth_idx2    = r.truth_idx[1];
        out.truth_dphi    = r.truth_dphi;

        out.truth_pt1  = truth_jet_pT  -> at( r.truth_idx[0] );
        out.truth_eta1 = truth_jet_eta -> at( r.truth_idx[0] );
        out.truth_phi1 = truth_jet_phi -> at( r.truth_idx[0] );
        out.truth_E1   = truth_jet_E   -> at( r.truth_idx[0] );
        out.truth_pt2  = truth_jet_pT  -> at( r.truth_idx[1] );
        out.truth_eta2 = truth_jet_eta -> at( r.truth_idx[1] );
        out.truth_phi2 = truth_jet_phi -> at( r.truth_idx[1] );
        out.truth_E2   = truth_jet_E   -> at( r.truth_idx[1] );
        out.truth_xj   = ( out.truth_pt1 > 0.0f ) ? out.truth_pt2 / out.truth_pt1 : -999.0f;

        out.truth_flavor1 = flavor1;
        out.truth_flavor2 = flavor2;
        out.is_qq = is_qq ? 1 : 0;

        out.truth_match_idx1 = r.truth_match_idx[0];
        out.truth_match_idx2 = r.truth_match_idx[1];
        out.truth_match_pt1  = r.truth_match_pt[0];
        out.truth_match_pt2  = r.truth_match_pt[1];

        out.reco_idx1 = r.reco_idx[0];
        out.reco_idx2 = r.reco_idx[1];
        out.reco_prov1 = r.reco_prov[0];
        out.reco_prov2 = r.reco_prov[1];
        out.reco_truth_idx1 = r.reco_truth_idx[0];
        out.reco_truth_idx2 = r.reco_truth_idx[1];
        out.reco_dphi = r.reco_dphi;
        out.legs_swapped = r.legs_swapped;
        out.n_accepted_reco = r.n_accepted_reco;

        if ( r.reco_pair )
        {
            out.reco_pt1  = jet_pT  -> at( r.reco_idx[0] );
            out.reco_eta1 = jet_eta -> at( r.reco_idx[0] );
            out.reco_phi1 = jet_phi -> at( r.reco_idx[0] );
            out.reco_E1   = jet_E   -> at( r.reco_idx[0] );
            out.reco_pt2  = jet_pT  -> at( r.reco_idx[1] );
            out.reco_eta2 = jet_eta -> at( r.reco_idx[1] );
            out.reco_phi2 = jet_phi -> at( r.reco_idx[1] );
            out.reco_E2   = jet_E   -> at( r.reco_idx[1] );
            out.reco_xj   = ( out.reco_pt1 > 0.0f ) ? out.reco_pt2 / out.reco_pt1 : -999.0f;
        }
        else
        {
            out.reco_pt1 = out.reco_eta1 = out.reco_phi1 = out.reco_E1 = -999.0;
            out.reco_pt2 = out.reco_eta2 = out.reco_phi2 = out.reco_E2 = -999.0;
            out.reco_xj = -999.0;
        }

        out.tree -> Fill();
    }

    std::cout << "\n---- pair-level dijet categories, flavor split ----" << std::endl;
    std::cout << "Events with unmatched/other leading-dijet flavor (written to neither output): "
              << n_events_unmatched_flavor << std::endl;

    for ( Output * out : { &out_qq, &out_qg_gg } )
    {
        const bool qq = ( out == &out_qq );
        std::cout << ( qq ? "qq   : " : "qg/gg: " ) << out -> n_events << " events"
                  << " -- Fill: " << out -> n_cat[0]
                  << ", Miss: " << out -> n_cat[1]
                  << ", Fake: " << out -> n_cat[2]
                  << ", Skip: " << out -> n_cat[3]
                  << ", UESub: " << out -> n_cat[4] << std::endl;
        std::cout << "         swapped Fills: " << out -> n_swapped
                  << " | UESub legs taken by an unmatched (UE) jet: " << out -> n_ue_by_unmatched
                  << ", by a softer truth jet: " << out -> n_ue_by_other_truth << std::endl;

        out -> file -> cd();
        out -> tree -> Write();
        const long n_written = out -> tree -> GetEntries();
        out -> file -> Close();

        std::cout << "         wrote " << n_written << " rows to "
                  << ( qq ? outfile_qq : outfile_qg_gg ) << std::endl;
    }

    return 0;
}

#endif
