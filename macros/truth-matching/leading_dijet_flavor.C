
#ifndef _LEADING_DIJET_FLAVOR_C_
#define _LEADING_DIJET_FLAVOR_C_

#include <myana/AnaUtils.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

// for each event in the tree written by match.C, take the truth jets
// (matched or not -- reco matching status is irrelevant here), rank by
// pT, and classify the leading+subleading truth dijet by the sum of
// their flavor codes:
//   flavor_sum < 0             -> unmatched
//   0 < flavor_sum < 13        -> qq
//   21 < flavor_sum < 28       -> qg
//   flavor_sum == 42           -> gg
//   otherwise                  -> other
// only events with lead pt in [min_pt1, max_pt1), sublead pt in
// [min_pt2, max_pt2), and dphi(lead,sublead) >= min_dphi are classified.
// if require_lead_dijet_matched is true, both the leading and subleading
// truth jet must additionally have been matched to a reco jet of
// reco_type (match_status == 0) -- the truth-level pT ranking that picks
// out "leading"/"subleading" is unchanged, this only adds a pass/fail
// requirement on top of it.
// reco_type just needs to name one collection present in the file ("rho"
// or "sub1"): match.C runs the same truth jet selection for both, so
// either one gives each truth jet exactly once per event. the resulting
// category fractions are saved as a TH1F "h_flavor_fraction" in outfile.
int leading_dijet_flavor(
    const std::string & infile =  "/sphenix/user/tmengel/dijet-ana-auau/macros/truth-matching/rootfiles/08_23_2026_v001/jet30_hijing_scaled_all.root",
    const std::string & outfile = "dijet_flavor_fractions.root",
    const std::string & reco_type = "rho",
    const float min_pt1 = 30.0,
    const float max_pt1 = 43.2,
    const float min_pt2 = 10.1,
    const float max_pt2 = 999.0,
    const float min_dphi = (7.0/8.0)*TMath::Pi(),
    const bool require_lead_dijet_matched = false
)
{
    int want_reco_type = -1;
    if ( reco_type == "rho" )       want_reco_type = 0;
    else if ( reco_type == "sub1" ) want_reco_type = 1;
    else
    {
        std::cerr << "Error: reco_type must be \"rho\" or \"sub1\", got \"" << reco_type << "\"." << std::endl;
        return -1;
    }

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

    AnaUtils::MatchedJetRow row;
    AnaUtils::read_matched_jet_tree( t, row );

    const int nentries = t -> GetEntries();

    // match.C always writes every row of an event contiguously, so a new
    // event is detected the moment event_id changes between rows.
    int current_event_id = -1;
    bool have_event = false;

    struct TruthJetInfo { float pt; float phi; int flavor; int match_status; };
    std::vector<TruthJetInfo> truth_jets;

    long n_events_seen = 0;
    long n_events_with_dijet = 0;

    enum FlavorCategory { kUnmatched = 0, kQQ = 1, kQG = 2, kGG = 3, kOther = 4, kNCategories = 5 };
    const char * category_labels[kNCategories] = { "unmatched", "qq", "qg", "gg", "other" };
    long n_category[kNCategories] = { 0, 0, 0, 0, 0 };

    auto categorize = []( const int flavor1, const int flavor2 ) -> FlavorCategory
    {
        // const int flavor_sum = flavor1 + flavor2;
        // if ( flavor_sum < 0 )                     return kUnmatched;
        // if ( flavor_sum > 0  && flavor_sum < 13 ) return kQQ;
        // if ( flavor_sum > 21 && flavor_sum < 28 ) return kQG;
        // if ( flavor_sum == 42 )                   return kGG;
        // return kOther;
        bool isq1 = (flavor1 >= 1 && flavor1 <= 6);
        bool isq2 = (flavor2 >= 1 && flavor2 <= 6);
        bool isg1 = (flavor1 == 21);
        bool isg2 = (flavor2 == 21);
        if ( isq1 && isq2 ) return kQQ;
        if ( (isq1 && isg2) || (isg1 && isq2) ) return kQG;
        if ( isg1 && isg2 ) return kGG;
        if ( flavor1 <= 0 || flavor2 <= 0 ) return kUnmatched;
        return kOther;
    };
   

    auto process_event = [&]()
    {
        ++n_events_seen;
        if ( truth_jets.size() < 2 ) return;

        std::sort( truth_jets.begin(), truth_jets.end(),
                   []( const TruthJetInfo & a, const TruthJetInfo & b ) { return a.pt > b.pt; } );

        const TruthJetInfo & lead = truth_jets[0];
        const TruthJetInfo & sublead = truth_jets[1];
        ++n_events_with_dijet;

        if ( lead.pt < min_pt1 || lead.pt >= max_pt1 ) return;
        if ( sublead.pt < min_pt2 || sublead.pt >= max_pt2 ) return;
        const float dphi = AnaUtils::dphi_wrap( lead.phi, sublead.phi );
        if ( dphi < min_dphi ) return;
        if ( require_lead_dijet_matched && ( lead.match_status != 0 || sublead.match_status != 0 ) ) return;

        int lf = lead.flavor;
        if ( lf == 0 ) lf = - 999;
        int sf = sublead.flavor;
        if ( sf == 0 ) sf = -999;
        const int flavor_sum = lf + sf;
        // const int flavor_sum = lead.flavor + sublead.flavor;
        ++n_category[ categorize( lead.flavor, sublead.flavor ) ];
        if ( categorize( lead.flavor, sublead.flavor ) == kOther )
        {
            std::cout << "Event " << current_event_id << " lead/sublead flavor sum = " << flavor_sum
                      << " (lead flavor = " << lead.flavor << ", sublead flavor = " << sublead.flavor << ")"
                      << std::endl;
        }
    };

    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );

        if ( !have_event || row.event_id != current_event_id )
        {
            if ( have_event ) process_event();
            truth_jets.clear();
            current_event_id = row.event_id;
            have_event = true;
        }

        if ( row.reco_type != want_reco_type ) continue;
        if ( row.match_status == 2 ) continue; // unmatched reco row -- no truth jet here

        truth_jets.push_back( { row.truth_pt, row.truth_phi, row.truth_flavor, row.match_status } );
    }
    if ( have_event ) process_event();

    long n_passing_cuts = 0;
    for ( int i = 0; i < kNCategories; ++i ) n_passing_cuts += n_category[i];
    // don't count unmatched 
    n_passing_cuts -= n_category[kUnmatched];
    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * h_flavor_fraction = new TH1F(
        "h_flavor_fraction", "leading+subleading truth dijet flavor;;fraction of flavor-matched dijet pairs",
        kNCategories, 0, kNCategories
    );
    for ( int i = 0; i < kNCategories; ++i )
    {
        h_flavor_fraction -> GetXaxis() -> SetBinLabel( i + 1, category_labels[i] );
        const double fraction = ( n_passing_cuts > 0 ) ? static_cast<double>( n_category[i] ) / n_passing_cuts : 0.0;
        h_flavor_fraction -> SetBinContent( i + 1, fraction );
    }
    fout -> cd();
    h_flavor_fraction -> Write();
    fout -> Close();

    std::cout << n_events_with_dijet << " / " << n_events_seen << " events had a leading+subleading truth dijet, "
              << n_passing_cuts << " passed the pT1/pT2/dphi cuts." << std::endl;
    for ( int i = 0; i < kNCategories; ++i )
    {
        const double pct = ( n_passing_cuts > 0 ) ? 100.0 * n_category[i] / n_passing_cuts : 0.0;
        std::cout << "  " << category_labels[i] << ": " << n_category[i] << " (" << pct << "%)" << std::endl;
    }
    std::cout << "Wrote " << outfile << std::endl;

    return 0;
}

#endif
