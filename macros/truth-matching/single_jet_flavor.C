
#ifndef _SINGLE_JET_FLAVOR_C_
#define _SINGLE_JET_FLAVOR_C_

#include <myana/AnaUtils.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

#include <iostream>
#include <string>

R__LOAD_LIBRARY( libmyana.so )

// for each truth jet (matched or not -- reco matching status is
// irrelevant unless require_jet_matched is set) in the tree written by
// match.C, classify by its own flavor code:
//   flavor <= 0       -> unmatched (no hard parton found in the jet)
//   1 <= flavor <= 6  -> quark
//   flavor == 21      -> gluon
//   otherwise         -> other
// only jets with pt in [min_pt, max_pt) are classified. if
// require_jet_matched is true, the jet must additionally have been
// matched to a reco jet of reco_type (match_status == 0). reco_type just
// needs to name one collection present in the file ("rho" or "sub1"):
// match.C runs the same truth jet selection for both, so either one
// gives each truth jet exactly once (no double counting). the resulting
// category fractions are saved as a TH1F "h_jet_flavor_fraction" in
// outfile.
int single_jet_flavor(
    const std::string & infile = "/sphenix/user/tmengel/dijet-ana-auau/macros/truth-matching/rootfiles/08_23_2026_v001/jet30_hijing_scaled_all.root",
    const std::string & outfile = "single_jet_flavor_fractions.root",
    const std::string & reco_type = "rho",
    const float min_pt = 30.0,
    const float max_pt = 43.2,
    const bool require_jet_matched = true
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

    enum FlavorCategory { kUnmatched = 0, kQuark = 1, kGluon = 2, kOther = 3, kNCategories = 4 };
    const char * category_labels[kNCategories] = { "unmatched", "quark", "gluon", "other" };
    long n_category[kNCategories] = { 0, 0, 0, 0 };

    auto categorize = []( const int flavor ) -> FlavorCategory
    {
        if ( flavor <= 0 )                return kUnmatched;
        if ( flavor >= 1 && flavor <= 6 ) return kQuark;
        if ( flavor == 21 )               return kGluon;
        return kOther;
    };

    long n_jets_considered = 0;
    long n_parton_matched = 0;
    double sum_parton_dr = 0.0;

    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );

        if ( row.reco_type != want_reco_type ) continue;
        if ( row.match_status == 2 ) continue; // unmatched reco row -- no truth jet here
        if ( row.truth_pt < min_pt || row.truth_pt >= max_pt ) continue;
        if ( require_jet_matched && row.match_status != 0 ) continue;

        ++n_jets_considered;
        const FlavorCategory category = categorize( row.truth_flavor );
        ++n_category[ category ];

        if ( category != kUnmatched ) // truth_parton_dr is only meaningful when a parton was found
        {
            ++n_parton_matched;
            sum_parton_dr += row.truth_parton_dr;
        }
    }

    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * h_jet_flavor_fraction = new TH1F(
        "h_jet_flavor_fraction", "single truth jet flavor;;fraction of truth jets",
        kNCategories, 0, kNCategories
    );
    for ( int i = 0; i < kNCategories; ++i )
    {
        h_jet_flavor_fraction -> GetXaxis() -> SetBinLabel( i + 1, category_labels[i] );
        const double fraction = ( n_jets_considered > 0 ) ? static_cast<double>( n_category[i] ) / n_jets_considered : 0.0;
        h_jet_flavor_fraction -> SetBinContent( i + 1, fraction );
    }
    fout -> cd();
    h_jet_flavor_fraction -> Write();
    fout -> Close();

    std::cout << n_jets_considered << " truth jets passed the pT/match cuts." << std::endl;
    for ( int i = 0; i < kNCategories; ++i )
    {
        const double pct = ( n_jets_considered > 0 ) ? 100.0 * n_category[i] / n_jets_considered : 0.0;
        std::cout << "  " << category_labels[i] << ": " << n_category[i] << " (" << pct << "%)" << std::endl;
    }
    const double mean_parton_dr = ( n_parton_matched > 0 ) ? sum_parton_dr / n_parton_matched : 0.0;
    std::cout << "  <dR(jet,parton)> = " << mean_parton_dr << " (over " << n_parton_matched << " parton-matched jets)" << std::endl;
    std::cout << "Wrote " << outfile << std::endl;

    return 0;
}

#endif
