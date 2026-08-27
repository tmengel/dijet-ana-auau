#ifndef _DRAW_PLOTS_C_
#define _DRAW_PLOTS_C_

#include <myana/AnaUtils.h>

#include <sPhenixStyle.C>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TF1.h>
#include <TStyle.h>
#include <TString.h>

#include <iostream>
#include <string>

R__LOAD_LIBRARY( libmyana.so )

// draws the histograms/graphs written by parse_tree.C (random cone
// background fluctuation, rho- vs sub1-subtracted jet kinematics, and the
// rho<->sub1 jet matching) into a single multi-page PDF, sPHENIX-style.
int draw_plots(
    const std::string & infile = "parse_tree.root",
    const std::string & outfile = "parse_tree_plots.pdf",
    const std::string & label = "Au+Au #sqrt{s_{NN}} = 200 GeV",
    const std::string & sphenix_tag = "#it{#bf{sPHENIX}} Internal"
)
{
    SetsPhenixStyle();
    gStyle -> SetOptStat( 0 );

    auto * f = TFile::Open( infile.c_str() );
    if ( !f || f -> IsZombie() )
    {
        std::cerr << "Error: could not open " << infile << std::endl;
        return -1;
    }

    bool ok = true;
    auto get_h1 = [&]( const char * name ) -> TH1D*
    {
        auto * h = (TH1D*) f -> Get( name );
        if ( !h ) { std::cerr << "Error: missing histogram \"" << name << "\"" << std::endl; ok = false; }
        return h;
    };
    auto get_h2 = [&]( const char * name ) -> TH2D*
    {
        auto * h = (TH2D*) f -> Get( name );
        if ( !h ) { std::cerr << "Error: missing histogram \"" << name << "\"" << std::endl; ok = false; }
        return h;
    };
    auto get_g = [&]( const char * name ) -> TGraphErrors*
    {
        auto * g = (TGraphErrors*) f -> Get( name );
        if ( !g ) { std::cerr << "Error: missing graph \"" << name << "\"" << std::endl; ok = false; }
        return g;
    };

    auto * h_rc_pt_all   = get_h1( "h_rc_pt_all" );
    auto * g_rc_mean_cent  = get_g( "g_rc_mean_cent" );
    auto * g_rc_sigma_cent = get_g( "g_rc_sigma_cent" );

    auto * h_rc_pt_sub1_all     = get_h1( "h_rc_pt_sub1_all" );
    auto * g_rc_sub1_mean_cent  = get_g( "g_rc_sub1_mean_cent" );
    auto * g_rc_sub1_sigma_cent = get_g( "g_rc_sub1_sigma_cent" );

    auto * h_rc_pt_corr_tower_vs_sub1 = get_h2( "h_rc_pt_corr_tower_vs_sub1" );
    auto * h_rc_dpt_sub1_minus_tower  = get_h1( "h_rc_dpt_sub1_minus_tower" );

    // RAW/unsubtracted sanity check (rho_jet's constituent_E is not
    // per-tower rho-subtracted -- see parse_tree.C for why), not a
    // subtraction cross-check like the sub1 one above.
    auto * h_rc_pt_rhoconst_all          = get_h1( "h_rc_pt_rhoconst_all" );
    auto * g_rc_rhoconst_mean_cent       = get_g( "g_rc_rhoconst_mean_cent" );
    auto * g_rc_rhoconst_sigma_cent      = get_g( "g_rc_rhoconst_sigma_cent" );
    auto * h_rc_pt_corr_tower_vs_rhoconst = get_h2( "h_rc_pt_corr_tower_vs_rhoconst" );
    auto * h_rc_dpt_rhoconst_minus_tower  = get_h1( "h_rc_dpt_rhoconst_minus_tower" );

    // manual (from-scratch, starting from the RAW towers above) re-derivations
    // of the sub1 and rho subtractions, cross-checking the two "baked-in"
    // ones (h_rc_pt_sub1_all, h_rc_pt_all) independently.
    auto * h_rc_pt_manualsub1_all    = get_h1( "h_rc_pt_manualsub1_all" );
    auto * g_rc_manualsub1_mean_cent  = get_g( "g_rc_manualsub1_mean_cent" );
    auto * g_rc_manualsub1_sigma_cent = get_g( "g_rc_manualsub1_sigma_cent" );

    auto * h_rc_pt_rhouncalib_all     = get_h1( "h_rc_pt_rhouncalib_all" );
    auto * g_rc_rhouncalib_mean_cent  = get_g( "g_rc_rhouncalib_mean_cent" );
    auto * g_rc_rhouncalib_sigma_cent = get_g( "g_rc_rhouncalib_sigma_cent" );

    auto * h_rc_pt_rhocalib_all     = get_h1( "h_rc_pt_rhocalib_all" );
    auto * g_rc_rhocalib_mean_cent  = get_g( "g_rc_rhocalib_mean_cent" );
    auto * g_rc_rhocalib_sigma_cent = get_g( "g_rc_rhocalib_sigma_cent" );

    auto * h_rc_pt_corr_manualsub1_vs_sub1const = get_h2( "h_rc_pt_corr_manualsub1_vs_sub1const" );
    auto * h_rc_dpt_manualsub1_minus_sub1const  = get_h1( "h_rc_dpt_manualsub1_minus_sub1const" );
    auto * h_rc_pt_corr_rhocalib_vs_tower       = get_h2( "h_rc_pt_corr_rhocalib_vs_tower" );
    auto * h_rc_dpt_rhocalib_minus_tower        = get_h1( "h_rc_dpt_rhocalib_minus_tower" );
    auto * h_rc_pt_corr_rhouncalib_vs_rhocalib  = get_h2( "h_rc_pt_corr_rhouncalib_vs_rhocalib" );

    // seed-protected variants: a tower with eT > sqrt(2)*rho' is left
    // unsubtracted (treated as a jet fragment, not background).
    auto * h_rc_pt_rhouncalib_seed_all     = get_h1( "h_rc_pt_rhouncalib_seed_all" );
    auto * g_rc_rhouncalib_seed_mean_cent  = get_g( "g_rc_rhouncalib_seed_mean_cent" );
    auto * g_rc_rhouncalib_seed_sigma_cent = get_g( "g_rc_rhouncalib_seed_sigma_cent" );

    auto * h_rc_pt_rhocalib_seed_all     = get_h1( "h_rc_pt_rhocalib_seed_all" );
    auto * g_rc_rhocalib_seed_mean_cent  = get_g( "g_rc_rhocalib_seed_mean_cent" );
    auto * g_rc_rhocalib_seed_sigma_cent = get_g( "g_rc_rhocalib_seed_sigma_cent" );

    auto * h_rc_dpt_rhouncalib_seed_minus_rhouncalib = get_h1( "h_rc_dpt_rhouncalib_seed_minus_rhouncalib" );
    auto * h_rc_dpt_rhocalib_seed_minus_rhocalib     = get_h1( "h_rc_dpt_rhocalib_seed_minus_rhocalib" );
    auto * h_rc_pt_corr_rhocalib_seed_vs_tower       = get_h2( "h_rc_pt_corr_rhocalib_seed_vs_tower" );
    auto * h_rc_dpt_rhocalib_seed_minus_tower        = get_h1( "h_rc_dpt_rhocalib_seed_minus_tower" );

    auto * h_rho_jet_pt   = get_h1( "h_rho_jet_pt" );
    auto * h_sub1_jet_pt  = get_h1( "h_sub1_jet_pt" );
    auto * h_rho_jet_eta  = get_h1( "h_rho_jet_eta" );
    auto * h_sub1_jet_eta = get_h1( "h_sub1_jet_eta" );
    auto * h_rho_jet_phi  = get_h1( "h_rho_jet_phi" );
    auto * h_sub1_jet_phi = get_h1( "h_sub1_jet_phi" );
    auto * h_njets_rho    = get_h1( "h_njets_rho" );
    auto * h_njets_sub1   = get_h1( "h_njets_sub1" );

    auto * h_match_pt_rho_vs_sub1 = get_h2( "h_match_pt_rho_vs_sub1" );
    auto * h_match_dpt       = get_h1( "h_match_dpt" );
    auto * h_match_dpt_relpt = get_h1( "h_match_dpt_relpt" );
    auto * h_match_dr        = get_h1( "h_match_dr" );
    auto * h_match_deta      = get_h1( "h_match_deta" );
    auto * h_match_dphi      = get_h1( "h_match_dphi" );

    if ( !ok )
    {
        std::cerr << "Error: one or more required objects were missing from " << infile << ", aborting." << std::endl;
        return -1;
    }

    const int kRho         = kAzure + 2;
    const int kSub1        = kRed + 1;
    const int kRhoConst    = kGreen + 2;
    const int kManualSub1  = kMagenta + 1;
    const int kRhoUncalib  = kOrange + 7;
    const int kRhoCalib    = kViolet + 1;
    const int kRhoUncalibSeed = kSpring + 4;
    const int kRhoCalibSeed   = kPink + 6;

    auto draw_header = [&]()
    {
        AnaUtils::myText( 0.16, 0.955, kBlack, sphenix_tag.c_str(), 0.045 );
        if ( !label.empty() )
        {
            TLatex tl;
            tl.SetNDC();
            tl.SetTextAlign( 31 );
            tl.SetTextSize( 0.028 );
            tl.DrawLatex( 0.95, 0.955, label.c_str() );
        }
    };

    auto style_pad = [&]( TVirtualPad * p )
    {
        p -> SetTopMargin( 0.10 );
        p -> SetLeftMargin( 0.14 );
        p -> SetRightMargin( 0.05 );
        p -> SetBottomMargin( 0.13 );
    };

    auto * c = new TCanvas( "c_random_cones_plots", "c_random_cones_plots", 800, 700 );

    const std::string book  = outfile + "(";
    const std::string mid   = outfile;
    const std::string close = outfile + ")";

    //--------------------------------------------------------------------
    // page 1: random cone pT distribution
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        gPad -> SetLogy();
        c -> SetTopMargin( 0.14 );

        h_rc_pt_all -> SetLineColor( kRho );
        h_rc_pt_all -> SetLineWidth( 2 );
        h_rc_pt_sub1_all -> SetLineColor( kSub1 );
        h_rc_pt_sub1_all -> SetLineWidth( 2 );

        h_rc_pt_all -> Scale( 1.0 / h_rc_pt_all -> Integral() );
        h_rc_pt_sub1_all -> Scale( 1.0 / h_rc_pt_sub1_all -> Integral() );

        h_rc_pt_all -> GetXaxis() -> SetTitle( "random cone p_{T} [GeV]" );
        h_rc_pt_all -> GetYaxis() -> SetTitle( "cones" );
        h_rc_pt_all -> GetYaxis() -> SetTitleOffset( 1.5 );
        // h_rc_pt_all -> SetMaximum( 1.5 * std::max( h_rc_pt_all -> GetMaximum(), h_rc_pt_sub1_all -> GetMaximum() ) );
        h_rc_pt_all -> GetYaxis() -> SetRangeUser( 1e-3, 1.5 * std::max( h_rc_pt_all -> GetMaximum(), h_rc_pt_sub1_all -> GetMaximum() ) );
        h_rc_pt_all -> GetXaxis() -> SetRangeUser( -15, 15 );
        h_rc_pt_all -> Draw( "hist" );
        h_rc_pt_sub1_all -> Draw( "hist SAME" );

        auto * fgaus = new TF1( "fgaus", "gaus",
            h_rc_pt_all -> GetMean() - 2.0 * h_rc_pt_all -> GetRMS(),
            h_rc_pt_all -> GetMean() + 2.0 * h_rc_pt_all -> GetRMS()
        );
        fgaus -> SetLineColor( kRho );
        fgaus -> SetLineStyle( 2 );
        h_rc_pt_all -> Fit( fgaus, "RQ0" );
        fgaus -> Draw( "SAME" );

        auto * fgaus_sub1 = new TF1( "fgaus_sub1", "gaus",
            h_rc_pt_sub1_all -> GetMean() - 2.0 * h_rc_pt_sub1_all -> GetRMS(),
            h_rc_pt_sub1_all -> GetMean() + 2.0 * h_rc_pt_sub1_all -> GetRMS()
        );
        fgaus_sub1 -> SetLineColor( kSub1 );
        fgaus_sub1 -> SetLineStyle( 2 );
        h_rc_pt_sub1_all -> Fit( fgaus_sub1, "RQ0" );
        fgaus_sub1 -> Draw( "SAME" );

        auto * leg = new TLegend( 0.56, 0.60, 0.94, 0.80 );
        leg -> SetFillStyle( 0 );
        leg -> SetBorderSize( 0 );
        leg -> SetTextSize( 0.024 );
        leg -> AddEntry( h_rc_pt_all, "tower (rho)", "l" );
        leg -> AddEntry( "", Form( "Mean=%.3f  RMS=%.3f", h_rc_pt_all -> GetMean(), h_rc_pt_all -> GetRMS() ), "" );
        leg -> AddEntry( h_rc_pt_sub1_all, "sub1 constituents", "l" );
        leg -> AddEntry( "", Form( "Mean=%.3f  RMS=%.3f", h_rc_pt_sub1_all -> GetMean(), h_rc_pt_sub1_all -> GetRMS() ), "" );
        leg -> Draw();

        draw_header();
        c -> RedrawAxis();
        c -> Print( book.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2: random cone mean/sigma vs centrality
    //--------------------------------------------------------------------
    {
        c -> Clear();
        gPad -> SetLogy( 0 );
        c -> Divide( 2, 1 );

        c -> cd( 1 );
        style_pad( gPad );
        g_rc_mean_cent -> SetMarkerStyle( 20 );
        g_rc_mean_cent -> SetMarkerColor( kRho );
        g_rc_mean_cent -> SetLineColor( kRho );
        g_rc_sub1_mean_cent -> SetMarkerStyle( 21 );
        g_rc_sub1_mean_cent -> SetMarkerColor( kSub1 );
        g_rc_sub1_mean_cent -> SetLineColor( kSub1 );

        auto * mg_mean = new TMultiGraph( "mg_rc_mean_cent", ";centrality [%];<p_{T}^{cone}> [GeV]" );
        mg_mean -> Add( g_rc_mean_cent, "P" );
        mg_mean -> Add( g_rc_sub1_mean_cent, "P" );
        mg_mean -> Draw( "A" );
        auto * l0 = new TLine( 0, 0, 90, 0 );
        l0 -> SetLineStyle( 2 );
        l0 -> SetLineColor( kGray + 2 );
        l0 -> Draw();
        mg_mean -> Draw( "P SAME" );

        auto * leg_mean = new TLegend( 0.30, 0.16, 0.94, 0.30 );
        leg_mean -> SetFillStyle( 0 );
        leg_mean -> SetBorderSize( 0 );
        leg_mean -> SetTextSize( 0.03 );
        leg_mean -> AddEntry( g_rc_mean_cent, "tower (rho)", "p" );
        leg_mean -> AddEntry( g_rc_sub1_mean_cent, "sub1 constituents", "p" );
        leg_mean -> Draw();

        c -> cd( 2 );
        style_pad( gPad );
        g_rc_sigma_cent -> SetMarkerStyle( 20 );
        g_rc_sigma_cent -> SetMarkerColor( kRho );
        g_rc_sigma_cent -> SetLineColor( kRho );
        g_rc_sub1_sigma_cent -> SetMarkerStyle( 21 );
        g_rc_sub1_sigma_cent -> SetMarkerColor( kSub1 );
        g_rc_sub1_sigma_cent -> SetLineColor( kSub1 );

        auto * mg_sigma = new TMultiGraph( "mg_rc_sigma_cent", ";centrality [%];#sigma(p_{T}^{cone}) [GeV]" );
        mg_sigma -> Add( g_rc_sigma_cent, "P" );
        mg_sigma -> Add( g_rc_sub1_sigma_cent, "P" );
        mg_sigma -> Draw( "A" );

        auto * leg_sigma = new TLegend( 0.30, 0.76, 0.94, 0.90 );
        leg_sigma -> SetFillStyle( 0 );
        leg_sigma -> SetBorderSize( 0 );
        leg_sigma -> SetTextSize( 0.03 );
        leg_sigma -> AddEntry( g_rc_sigma_cent, "tower (rho)", "p" );
        leg_sigma -> AddEntry( g_rc_sub1_sigma_cent, "sub1 constituents", "p" );
        leg_sigma -> Draw();

        c -> cd( 0 );
        draw_header();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2b: random cone pT, tower(rho) vs sub1-constituent, cone-by-cone
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetRightMargin( 0.14 );
        c -> SetLogy( 0 );
        c -> SetLogz();

        h_rc_pt_corr_tower_vs_sub1 -> GetXaxis() -> SetTitle( "tower (rho-subtracted) cone p_{T} [GeV]" );
        h_rc_pt_corr_tower_vs_sub1 -> GetYaxis() -> SetTitle( "sub1-subtracted constituent cone p_{T} [GeV]" );
        h_rc_pt_corr_tower_vs_sub1 -> GetYaxis() -> SetTitleOffset( 1.5 );
        h_rc_pt_corr_tower_vs_sub1 -> Draw( "COLZ" );

        const double rc_xmin = h_rc_pt_corr_tower_vs_sub1 -> GetXaxis() -> GetXmin();
        const double rc_xmax = h_rc_pt_corr_tower_vs_sub1 -> GetXaxis() -> GetXmax();
        auto * rc_ldiag = new TLine( rc_xmin, rc_xmin, rc_xmax, rc_xmax );
        rc_ldiag -> SetLineStyle( 2 );
        rc_ldiag -> SetLineColor( kBlack );
        rc_ldiag -> SetLineWidth( 2 );
        rc_ldiag -> Draw();

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetLogy();

        h_rc_dpt_sub1_minus_tower -> SetLineColor( kRho );
        h_rc_dpt_sub1_minus_tower -> SetLineWidth( 2 );
        h_rc_dpt_sub1_minus_tower -> GetXaxis() -> SetTitle( "p_{T}^{cone,sub1} - p_{T}^{cone,tower} [GeV]" );
        h_rc_dpt_sub1_minus_tower -> GetYaxis() -> SetTitle( "cones" );
        h_rc_dpt_sub1_minus_tower -> GetYaxis() -> SetTitleOffset( 1.5 );
        h_rc_dpt_sub1_minus_tower -> Draw( "hist" );
        auto * rc_l0 = new TLine( 0, 0, 0, h_rc_dpt_sub1_minus_tower -> GetMaximum() );
        rc_l0 -> SetLineStyle( 2 );
        rc_l0 -> SetLineColor( kGray + 2 );
        rc_l0 -> Draw();
        h_rc_dpt_sub1_minus_tower -> Draw( "hist SAME" );

        TLatex tl_rc;
        tl_rc.SetNDC();
        tl_rc.SetTextSize( 0.03 );
        tl_rc.DrawLatex( 0.60, 0.80, Form( "Mean = %.3f #pm %.3f", h_rc_dpt_sub1_minus_tower -> GetMean(), h_rc_dpt_sub1_minus_tower -> GetMeanError() ) );
        tl_rc.DrawLatex( 0.60, 0.74, Form( "RMS  = %.3f #pm %.3f", h_rc_dpt_sub1_minus_tower -> GetRMS(), h_rc_dpt_sub1_minus_tower -> GetRMSError() ) );

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2c: random cone pT, RAW/unsubtracted via rho_jet constituents.
    // NOT a subtraction cross-check (rho_jet's constituent_E has no rho
    // term subtracted -- the rho correction is jet-level only, see
    // parse_tree.C) -- this is an independent sanity check on the typical
    // raw underlying-event energy density in a cone.
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        gPad -> SetLogy();
        c -> SetTopMargin( 0.14 );

        h_rc_pt_all -> SetLineColor( kRho );
        h_rc_pt_rhoconst_all -> SetLineColor( kRhoConst );
        h_rc_pt_rhoconst_all -> SetLineWidth( 2 );

        h_rc_pt_rhoconst_all -> GetXaxis() -> SetTitle( "random cone p_{T} [GeV]" );
        h_rc_pt_rhoconst_all -> GetYaxis() -> SetTitle( "cones" );
        h_rc_pt_rhoconst_all -> GetYaxis() -> SetTitleOffset( 1.5 );
        h_rc_pt_rhoconst_all -> Draw( "hist" );
        h_rc_pt_all -> Draw( "hist SAME" );

        auto * leg = new TLegend( 0.50, 0.62, 0.94, 0.80 );
        leg -> SetFillStyle( 0 );
        leg -> SetBorderSize( 0 );
        leg -> SetTextSize( 0.026 );
        leg -> AddEntry( h_rc_pt_rhoconst_all, "RAW (rho_jet constituents)", "l" );
        leg -> AddEntry( "", Form( "Mean=%.3f  RMS=%.3f", h_rc_pt_rhoconst_all -> GetMean(), h_rc_pt_rhoconst_all -> GetRMS() ), "" );
        leg -> AddEntry( h_rc_pt_all, "tower (rho-subtracted)", "l" );
        leg -> AddEntry( "", Form( "Mean=%.3f  RMS=%.3f", h_rc_pt_all -> GetMean(), h_rc_pt_all -> GetRMS() ), "" );
        leg -> Draw();

        TLatex tl_note;
        tl_note.SetNDC();
        tl_note.SetTextSize( 0.024 );
        tl_note.SetTextColor( kGray + 2 );
        tl_note.DrawLatex( 0.16, 0.20, "rho_jet constituents are NOT per-tower subtracted -- sanity check only, not a subtraction cross-check" );

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    {
        c -> Clear();
        gPad -> SetLogy( 0 );
        c -> Divide( 2, 1 );

        c -> cd( 1 );
        style_pad( gPad );
        g_rc_rhoconst_mean_cent -> SetMarkerStyle( 22 );
        g_rc_rhoconst_mean_cent -> SetMarkerColor( kRhoConst );
        g_rc_rhoconst_mean_cent -> SetLineColor( kRhoConst );
        g_rc_rhoconst_mean_cent -> SetTitle( ";centrality [%];<p_{T}^{cone}> [GeV] (RAW)" );
        g_rc_rhoconst_mean_cent -> Draw( "AP" );

        c -> cd( 2 );
        style_pad( gPad );
        g_rc_rhoconst_sigma_cent -> SetMarkerStyle( 22 );
        g_rc_rhoconst_sigma_cent -> SetMarkerColor( kRhoConst );
        g_rc_rhoconst_sigma_cent -> SetLineColor( kRhoConst );
        g_rc_rhoconst_sigma_cent -> SetTitle( ";centrality [%];#sigma(p_{T}^{cone}) [GeV] (RAW)" );
        g_rc_rhoconst_sigma_cent -> Draw( "AP" );

        c -> cd( 0 );
        draw_header();
        c -> Print( mid.c_str() );
    }

    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetRightMargin( 0.14 );
        c -> SetLogz();

        h_rc_pt_corr_tower_vs_rhoconst -> GetXaxis() -> SetTitle( "tower (rho-subtracted) cone p_{T} [GeV]" );
        h_rc_pt_corr_tower_vs_rhoconst -> GetYaxis() -> SetTitle( "RAW (rho_jet constituent) cone p_{T} [GeV]" );
        h_rc_pt_corr_tower_vs_rhoconst -> GetYaxis() -> SetTitleOffset( 1.5 );
        h_rc_pt_corr_tower_vs_rhoconst -> Draw( "COLZ" );

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2d: all 5 subtracted random-cone reconstructions overlaid --
    // the 3 manual (from-scratch, from RAW towers) ones plus the 2
    // pre-existing ("baked-in") ones. RAW/unsubtracted is intentionally
    // excluded here (already shown on its own page above) since its scale
    // (mean~9, extending to >100 GeV) would swamp these ~O(2 GeV) ones.
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        gPad -> SetLogy();
        c -> SetTopMargin( 0.14 );

        // "manual" (from-scratch) series are dashed so they stay visible
        // even where they sit almost exactly on top of the baked-in ones
        // (which is the point -- see the validation pages below).
        struct Series { TH1D * h; int color; int style; const char * label; };
        Series all5[5] = {
            { h_rc_pt_all,            kRho,        1, "tower (rho, baked-in)" },
            { h_rc_pt_sub1_all,       kSub1,       1, "sub1 constituents (baked-in)" },
            { h_rc_pt_manualsub1_all, kManualSub1, 2, "manual sub1" },
            { h_rc_pt_rhouncalib_all, kRhoUncalib, 2, "manual rho, uncalibrated" },
            { h_rc_pt_rhocalib_all,   kRhoCalib,   2, "manual rho, calibrated" }
        };

        double ymax = 0;
        for ( auto & s : all5 ) ymax = std::max( ymax, s.h -> GetMaximum() );

        for ( int i = 0; i < 5; ++i )
        {
            all5[i].h -> SetLineColor( all5[i].color );
            all5[i].h -> SetLineStyle( all5[i].style );
            all5[i].h -> SetLineWidth( 2 );
            all5[i].h -> Scale( 1.0 / all5[i].h -> Integral() );
            all5[i].h -> GetXaxis() -> SetTitle( "random cone p_{T} [GeV]" );
            all5[i].h -> GetXaxis() -> SetRangeUser( -30, 30 );
            all5[i].h -> GetYaxis() -> SetTitle( "cones" );
            all5[i].h -> GetYaxis() -> SetTitleOffset( 1.5 );
            all5[i].h -> SetMaximum( 40.0 * ymax ); // big log-y headroom so the legend strip below stays clear of all 5 curves
            all5[i].h -> Draw( i == 0 ? "hist" : "hist SAME" );
        }

        // one compact row across the top of the frame -- with 40x headroom
        // above the peak, this stays clear of every curve at every x.
        auto * leg = new TLegend( 0.16, 0.79, 0.94, 0.855 );
        leg -> SetFillStyle( 0 );
        leg -> SetBorderSize( 0 );
        leg -> SetTextSize( 0.020 );
        leg -> SetNColumns( 3 );
        for ( auto & s : all5 ) leg -> AddEntry( s.h, s.label, "l" );
        leg -> Draw();

        TLatex tl_stats;
        tl_stats.SetNDC();
        tl_stats.SetTextSize( 0.020 );
        tl_stats.SetTextColor( kGray + 2 );
        double y = 0.60;
        for ( auto & s : all5 )
        {
            tl_stats.SetTextColor( s.color );
            tl_stats.DrawLatex( 0.16, y, Form( "#mu=%.2f  #sigma=%.2f", s.h -> GetMean(), s.h -> GetRMS() ) );
            y -= 0.045;
        }

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2e: mean/sigma vs centrality for all 5 subtracted reconstructions
    //--------------------------------------------------------------------
    {
        c -> Clear();
        gPad -> SetLogy( 0 );
        c -> Divide( 2, 1 );

        struct GSeries { TGraphErrors * g; int color; int marker; const char * label; };
        GSeries mean5[5] = {
            { g_rc_mean_cent,           kRho,        20, "tower (rho)" },
            { g_rc_sub1_mean_cent,      kSub1,       21, "sub1 const" },
            { g_rc_manualsub1_mean_cent, kManualSub1, 22, "manual sub1" },
            { g_rc_rhouncalib_mean_cent, kRhoUncalib, 23, "rho uncalib" },
            { g_rc_rhocalib_mean_cent,   kRhoCalib,   33, "rho calib" }
        };
        GSeries sigma5[5] = {
            { g_rc_sigma_cent,           kRho,        20, "tower (rho)" },
            { g_rc_sub1_sigma_cent,      kSub1,       21, "sub1 const" },
            { g_rc_manualsub1_sigma_cent, kManualSub1, 22, "manual sub1" },
            { g_rc_rhouncalib_sigma_cent, kRhoUncalib, 23, "rho uncalib" },
            { g_rc_rhocalib_sigma_cent,   kRhoCalib,   33, "rho calib" }
        };

        c -> cd( 1 );
        style_pad( gPad );
        auto * mg_mean5 = new TMultiGraph( "mg_rc_mean5_cent", ";centrality [%];<p_{T}^{cone}> [GeV]" );
        for ( auto & s : mean5 )
        {
            s.g -> SetMarkerStyle( s.marker );
            s.g -> SetMarkerColor( s.color );
            s.g -> SetLineColor( s.color );
            mg_mean5 -> Add( s.g, "P" );
        }
        mg_mean5 -> Draw( "A" );
        auto * l0 = new TLine( 0, 0, 90, 0 );
        l0 -> SetLineStyle( 2 );
        l0 -> SetLineColor( kGray + 2 );
        l0 -> Draw();
        mg_mean5 -> Draw( "P SAME" );

        auto * leg_mean5 = new TLegend( 0.30, 0.14, 0.94, 0.30 );
        leg_mean5 -> SetFillStyle( 0 );
        leg_mean5 -> SetBorderSize( 0 );
        leg_mean5 -> SetTextSize( 0.026 );
        leg_mean5 -> SetNColumns( 2 );
        for ( auto & s : mean5 ) leg_mean5 -> AddEntry( s.g, s.label, "p" );
        leg_mean5 -> Draw();

        c -> cd( 2 );
        style_pad( gPad );
        auto * mg_sigma5 = new TMultiGraph( "mg_rc_sigma5_cent", ";centrality [%];#sigma(p_{T}^{cone}) [GeV]" );
        for ( auto & s : sigma5 )
        {
            s.g -> SetMarkerStyle( s.marker );
            s.g -> SetMarkerColor( s.color );
            s.g -> SetLineColor( s.color );
            mg_sigma5 -> Add( s.g, "P" );
        }
        mg_sigma5 -> Draw( "A" );

        auto * leg_sigma5 = new TLegend( 0.30, 0.70, 0.94, 0.86 );
        leg_sigma5 -> SetFillStyle( 0 );
        leg_sigma5 -> SetBorderSize( 0 );
        leg_sigma5 -> SetTextSize( 0.026 );
        leg_sigma5 -> SetNColumns( 2 );
        for ( auto & s : sigma5 ) leg_sigma5 -> AddEntry( s.g, s.label, "p" );
        leg_sigma5 -> Draw();

        c -> cd( 0 );
        draw_header();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2f: validation -- manual sub1 (RAW - <UE>) vs the pre-existing
    // sub1_jet constituent cone. These should (and do) agree essentially
    // exactly, since both ultimately implement the same per-tower <UE>
    // subtraction -- this validates the raw-tower recovery, layer
    // decoding, and eta-ring lookup used for all 3 manual reconstructions.
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetRightMargin( 0.14 );
        c -> SetLogz();

        h_rc_pt_corr_manualsub1_vs_sub1const -> GetXaxis() -> SetTitle( "manual sub1 (RAW - <UE>) cone p_{T} [GeV]" );
        h_rc_pt_corr_manualsub1_vs_sub1const -> GetYaxis() -> SetTitle( "sub1_jet constituent cone p_{T} [GeV]" );
        h_rc_pt_corr_manualsub1_vs_sub1const -> GetYaxis() -> SetTitleOffset( 1.5 );
        h_rc_pt_corr_manualsub1_vs_sub1const -> Draw( "COLZ" );

        const double xmin1 = h_rc_pt_corr_manualsub1_vs_sub1const -> GetXaxis() -> GetXmin();
        const double xmax1 = h_rc_pt_corr_manualsub1_vs_sub1const -> GetXaxis() -> GetXmax();
        auto * ldiag1 = new TLine( xmin1, xmin1, xmax1, xmax1 );
        ldiag1 -> SetLineStyle( 2 );
        ldiag1 -> SetLineColor( kBlack );
        ldiag1 -> SetLineWidth( 2 );
        ldiag1 -> Draw();

        TLatex tl_corr1;
        tl_corr1.SetNDC();
        tl_corr1.SetTextSize( 0.03 );
        tl_corr1.DrawLatex( 0.18, 0.82, Form( "corr = %.4f", h_rc_pt_corr_manualsub1_vs_sub1const -> GetCorrelationFactor() ) );

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2g: validation -- manual calibrated-rho vs the pre-existing
    // tower-grid rho cone. Close (but not exact) agreement validates the
    // from-scratch rho_calib calibration path independently of the
    // pre-baked cemc/ihcal/ohcal_tower_E arrays.
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetRightMargin( 0.14 );
        c -> SetLogz();

        h_rc_pt_corr_rhocalib_vs_tower -> GetXaxis() -> SetTitle( "manual calibrated-rho cone p_{T} [GeV]" );
        h_rc_pt_corr_rhocalib_vs_tower -> GetYaxis() -> SetTitle( "tower (rho-subtracted) cone p_{T} [GeV]" );
        h_rc_pt_corr_rhocalib_vs_tower -> GetYaxis() -> SetTitleOffset( 1.5 );
        h_rc_pt_corr_rhocalib_vs_tower -> Draw( "COLZ" );

        const double xmin2 = h_rc_pt_corr_rhocalib_vs_tower -> GetXaxis() -> GetXmin();
        const double xmax2 = h_rc_pt_corr_rhocalib_vs_tower -> GetXaxis() -> GetXmax();
        auto * ldiag2 = new TLine( xmin2, xmin2, xmax2, xmax2 );
        ldiag2 -> SetLineStyle( 2 );
        ldiag2 -> SetLineColor( kBlack );
        ldiag2 -> SetLineWidth( 2 );
        ldiag2 -> Draw();

        TLatex tl_corr2;
        tl_corr2.SetNDC();
        tl_corr2.SetTextSize( 0.03 );
        tl_corr2.DrawLatex( 0.18, 0.82, Form( "corr = %.4f", h_rc_pt_corr_rhocalib_vs_tower -> GetCorrelationFactor() ) );

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2h: the two validation delta distributions side by side
    //--------------------------------------------------------------------
    {
        c -> Clear();
        c -> Divide( 2, 1 );

        auto draw_delta = [&]( TVirtualPad * p, TH1D * h, const char * xtitle )
        {
            p -> cd();
            p -> SetTopMargin( 0.10 );
            p -> SetLeftMargin( 0.18 );
            p -> SetRightMargin( 0.05 );
            p -> SetBottomMargin( 0.16 );
            h -> SetLineColor( kManualSub1 );
            h -> SetLineWidth( 2 );
            h -> GetXaxis() -> SetTitle( xtitle );
            h -> GetYaxis() -> SetTitle( "cones" );
            h -> GetXaxis() -> SetLabelSize( 0.05 );
            h -> GetXaxis() -> SetTitleSize( 0.05 );
            h -> GetYaxis() -> SetLabelSize( 0.05 );
            h -> GetYaxis() -> SetTitleSize( 0.05 );
            h -> GetYaxis() -> SetTitleOffset( 1.7 );
            h -> Draw( "hist" );
            auto * l0 = new TLine( 0, 0, 0, h -> GetMaximum() );
            l0 -> SetLineStyle( 2 );
            l0 -> SetLineColor( kGray + 2 );
            l0 -> Draw();
            h -> Draw( "hist SAME" );

            TLatex tl;
            tl.SetNDC();
            tl.SetTextSize( 0.04 );
            tl.DrawLatex( 0.24, 0.82, Form( "Mean = %.3g", h -> GetMean() ) );
            tl.DrawLatex( 0.24, 0.75, Form( "RMS  = %.3g", h -> GetRMS() ) );
        };

        draw_delta( c -> cd(1), h_rc_dpt_manualsub1_minus_sub1const, "p_{T}^{manual sub1} - p_{T}^{sub1 const} [GeV]" );
        draw_delta( c -> cd(2), h_rc_dpt_rhocalib_minus_tower, "p_{T}^{manual rho-calib} - p_{T}^{tower} [GeV]" );

        c -> cd( 0 );
        draw_header();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2i: effect of the eta/zvtx/mbdQ calibration -- uncalibrated vs
    // calibrated manual rho subtraction, same cone
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetRightMargin( 0.14 );
        c -> SetLogz();

        h_rc_pt_corr_rhouncalib_vs_rhocalib -> GetXaxis() -> SetTitle( "uncalibrated-rho cone p_{T} [GeV]" );
        h_rc_pt_corr_rhouncalib_vs_rhocalib -> GetYaxis() -> SetTitle( "calibrated-rho cone p_{T} [GeV]" );
        h_rc_pt_corr_rhouncalib_vs_rhocalib -> GetYaxis() -> SetTitleOffset( 1.5 );
        h_rc_pt_corr_rhouncalib_vs_rhocalib -> Draw( "COLZ" );

        const double xmin3 = h_rc_pt_corr_rhouncalib_vs_rhocalib -> GetXaxis() -> GetXmin();
        const double xmax3 = h_rc_pt_corr_rhouncalib_vs_rhocalib -> GetXaxis() -> GetXmax();
        auto * ldiag3 = new TLine( xmin3, xmin3, xmax3, xmax3 );
        ldiag3 -> SetLineStyle( 2 );
        ldiag3 -> SetLineColor( kBlack );
        ldiag3 -> SetLineWidth( 2 );
        ldiag3 -> Draw();

        TLatex tl_corr3;
        tl_corr3.SetNDC();
        tl_corr3.SetTextSize( 0.03 );
        tl_corr3.DrawLatex( 0.18, 0.82, Form( "corr = %.4f", h_rc_pt_corr_rhouncalib_vs_rhocalib -> GetCorrelationFactor() ) );

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2j: effect of seed/outlier protection -- a tower with
    // eT > sqrt(2)*rho' is left unsubtracted (treated as a jet fragment
    // rather than background). Overlays each rho variant against its
    // seed-protected counterpart.
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        gPad -> SetLogy();
        c -> SetTopMargin( 0.14 );

        struct Series { TH1D * h; int color; int style; const char * label; };
        Series all4[4] = {
            { h_rc_pt_rhouncalib_all,      kRhoUncalib,     1, "uncalibrated rho" },
            { h_rc_pt_rhouncalib_seed_all, kRhoUncalibSeed, 2, "uncalibrated rho, seed-protected" },
            { h_rc_pt_rhocalib_all,        kRhoCalib,       1, "calibrated rho" },
            { h_rc_pt_rhocalib_seed_all,   kRhoCalibSeed,   2, "calibrated rho, seed-protected" }
        };

        double ymax = 0;
        for ( auto & s : all4 ) ymax = std::max( ymax, s.h -> GetMaximum() );

        for ( int i = 0; i < 4; ++i )
        {
            all4[i].h -> SetLineColor( all4[i].color );
            all4[i].h -> SetLineStyle( all4[i].style );
            all4[i].h -> SetLineWidth( 2 );
            all4[i].h -> GetXaxis() -> SetTitle( "random cone p_{T} [GeV]" );
            all4[i].h -> GetXaxis() -> SetRangeUser( -30, 30 );
            all4[i].h -> GetYaxis() -> SetTitle( "cones" );
            all4[i].h -> GetYaxis() -> SetTitleOffset( 1.5 );
            all4[i].h -> SetMaximum( 40.0 * ymax );
            all4[i].h -> Draw( i == 0 ? "hist" : "hist SAME" );
        }

        auto * leg = new TLegend( 0.16, 0.79, 0.94, 0.855 );
        leg -> SetFillStyle( 0 );
        leg -> SetBorderSize( 0 );
        leg -> SetTextSize( 0.020 );
        leg -> SetNColumns( 2 );
        for ( auto & s : all4 ) leg -> AddEntry( s.h, s.label, "l" );
        leg -> Draw();

        TLatex tl_stats;
        tl_stats.SetNDC();
        tl_stats.SetTextSize( 0.020 );
        double y = 0.60;
        for ( auto & s : all4 )
        {
            tl_stats.SetTextColor( s.color );
            tl_stats.DrawLatex( 0.16, y, Form( "#mu=%.2f  #sigma=%.2f", s.h -> GetMean(), s.h -> GetRMS() ) );
            y -= 0.045;
        }

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2k: mean/sigma vs centrality, rho variants with/without seed
    // protection
    //--------------------------------------------------------------------
    {
        c -> Clear();
        gPad -> SetLogy( 0 );
        c -> Divide( 2, 1 );

        struct GSeries { TGraphErrors * g; int color; int marker; const char * label; };
        GSeries mean4[4] = {
            { g_rc_rhouncalib_mean_cent,      kRhoUncalib,     23, "rho uncalib" },
            { g_rc_rhouncalib_seed_mean_cent, kRhoUncalibSeed, 26, "rho uncalib, seed" },
            { g_rc_rhocalib_mean_cent,        kRhoCalib,       33, "rho calib" },
            { g_rc_rhocalib_seed_mean_cent,   kRhoCalibSeed,   27, "rho calib, seed" }
        };
        GSeries sigma4[4] = {
            { g_rc_rhouncalib_sigma_cent,      kRhoUncalib,     23, "rho uncalib" },
            { g_rc_rhouncalib_seed_sigma_cent, kRhoUncalibSeed, 26, "rho uncalib, seed" },
            { g_rc_rhocalib_sigma_cent,        kRhoCalib,       33, "rho calib" },
            { g_rc_rhocalib_seed_sigma_cent,   kRhoCalibSeed,   27, "rho calib, seed" }
        };

        c -> cd( 1 );
        style_pad( gPad );
        auto * mg_mean4 = new TMultiGraph( "mg_rc_mean4_cent", ";centrality [%];<p_{T}^{cone}> [GeV]" );
        for ( auto & s : mean4 )
        {
            s.g -> SetMarkerStyle( s.marker );
            s.g -> SetMarkerColor( s.color );
            s.g -> SetLineColor( s.color );
            mg_mean4 -> Add( s.g, "P" );
        }
        mg_mean4 -> Draw( "A" );
        auto * l0 = new TLine( 0, 0, 90, 0 );
        l0 -> SetLineStyle( 2 );
        l0 -> SetLineColor( kGray + 2 );
        l0 -> Draw();
        mg_mean4 -> Draw( "P SAME" );

        auto * leg_mean4 = new TLegend( 0.16, 0.70, 0.94, 0.86 );
        leg_mean4 -> SetFillStyle( 0 );
        leg_mean4 -> SetBorderSize( 0 );
        leg_mean4 -> SetTextSize( 0.026 );
        for ( auto & s : mean4 ) leg_mean4 -> AddEntry( s.g, s.label, "p" );
        leg_mean4 -> Draw();

        c -> cd( 2 );
        style_pad( gPad );
        auto * mg_sigma4 = new TMultiGraph( "mg_rc_sigma4_cent", ";centrality [%];#sigma(p_{T}^{cone}) [GeV]" );
        for ( auto & s : sigma4 )
        {
            s.g -> SetMarkerStyle( s.marker );
            s.g -> SetMarkerColor( s.color );
            s.g -> SetLineColor( s.color );
            mg_sigma4 -> Add( s.g, "P" );
        }
        mg_sigma4 -> Draw( "A" );

        auto * leg_sigma4 = new TLegend( 0.16, 0.70, 0.94, 0.86 );
        leg_sigma4 -> SetFillStyle( 0 );
        leg_sigma4 -> SetBorderSize( 0 );
        leg_sigma4 -> SetTextSize( 0.026 );
        for ( auto & s : sigma4 ) leg_sigma4 -> AddEntry( s.g, s.label, "p" );
        leg_sigma4 -> Draw();

        c -> cd( 0 );
        draw_header();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2l: does seed protection bring the manual calibrated-rho cone
    // closer to, or further from, the pre-existing tower-grid rho cone?
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetRightMargin( 0.14 );
        c -> SetLogz();

        h_rc_pt_corr_rhocalib_seed_vs_tower -> GetXaxis() -> SetTitle( "manual rho-calib (seed) cone p_{T} [GeV]" );
        h_rc_pt_corr_rhocalib_seed_vs_tower -> GetYaxis() -> SetTitle( "tower (rho-subtracted) cone p_{T} [GeV]" );
        h_rc_pt_corr_rhocalib_seed_vs_tower -> GetYaxis() -> SetTitleOffset( 1.5 );
        h_rc_pt_corr_rhocalib_seed_vs_tower -> Draw( "COLZ" );

        const double xmin4 = h_rc_pt_corr_rhocalib_seed_vs_tower -> GetXaxis() -> GetXmin();
        const double xmax4 = h_rc_pt_corr_rhocalib_seed_vs_tower -> GetXaxis() -> GetXmax();
        auto * ldiag4 = new TLine( xmin4, xmin4, xmax4, xmax4 );
        ldiag4 -> SetLineStyle( 2 );
        ldiag4 -> SetLineColor( kBlack );
        ldiag4 -> SetLineWidth( 2 );
        ldiag4 -> Draw();

        TLatex tl_corr4;
        tl_corr4.SetNDC();
        tl_corr4.SetTextSize( 0.03 );
        tl_corr4.DrawLatex( 0.18, 0.82, Form( "corr = %.4f", h_rc_pt_corr_rhocalib_seed_vs_tower -> GetCorrelationFactor() ) );

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 2m: seed-protection deltas -- how much each rho variant shifts
    // when outlier towers are excluded from subtraction
    //--------------------------------------------------------------------
    {
        c -> Clear();
        c -> Divide( 2, 1 );

        auto draw_delta = [&]( TVirtualPad * p, TH1D * h, const char * xtitle, const int color )
        {
            p -> cd();
            p -> SetTopMargin( 0.10 );
            p -> SetLeftMargin( 0.18 );
            p -> SetRightMargin( 0.05 );
            p -> SetBottomMargin( 0.16 );
            h -> SetLineColor( color );
            h -> SetLineWidth( 2 );
            h -> GetXaxis() -> SetTitle( xtitle );
            h -> GetYaxis() -> SetTitle( "cones" );
            h -> GetXaxis() -> SetLabelSize( 0.05 );
            h -> GetXaxis() -> SetTitleSize( 0.05 );
            h -> GetYaxis() -> SetLabelSize( 0.05 );
            h -> GetYaxis() -> SetTitleSize( 0.05 );
            h -> GetYaxis() -> SetTitleOffset( 1.7 );
            h -> Draw( "hist" );
            auto * l0 = new TLine( 0, 0, 0, h -> GetMaximum() );
            l0 -> SetLineStyle( 2 );
            l0 -> SetLineColor( kGray + 2 );
            l0 -> Draw();
            h -> Draw( "hist SAME" );

            TLatex tl;
            tl.SetNDC();
            tl.SetTextSize( 0.04 );
            tl.DrawLatex( 0.24, 0.82, Form( "Mean = %.3g", h -> GetMean() ) );
            tl.DrawLatex( 0.24, 0.75, Form( "RMS  = %.3g", h -> GetRMS() ) );
        };

        draw_delta( c -> cd(1), h_rc_dpt_rhouncalib_seed_minus_rhouncalib, "p_{T}^{uncalib,seed} - p_{T}^{uncalib} [GeV]", kRhoUncalibSeed );
        draw_delta( c -> cd(2), h_rc_dpt_rhocalib_seed_minus_rhocalib, "p_{T}^{calib,seed} - p_{T}^{calib} [GeV]", kRhoCalibSeed );

        c -> cd( 0 );
        draw_header();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 3: rho vs sub1 jet pT spectra
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetLogy();

        h_rho_jet_pt -> SetLineColor( kRho );
        h_rho_jet_pt -> SetLineWidth( 2 );
        h_sub1_jet_pt -> SetLineColor( kSub1 );
        h_sub1_jet_pt -> SetLineWidth( 2 );

        h_rho_jet_pt -> GetXaxis() -> SetTitle( "jet p_{T} [GeV]" );
        h_rho_jet_pt -> GetYaxis() -> SetTitle( "jets" );
        h_rho_jet_pt -> GetYaxis() -> SetTitleOffset( 1.5 );
        h_rho_jet_pt -> GetXaxis() -> SetRangeUser( 0, 40 );
        h_rho_jet_pt -> SetMaximum( 1.5 * std::max( h_rho_jet_pt -> GetMaximum(), h_sub1_jet_pt -> GetMaximum() ) );
        h_rho_jet_pt -> Draw( "hist" );
        h_sub1_jet_pt -> Draw( "hist SAME" );

        auto * leg = new TLegend( 0.42, 0.72, 0.94, 0.86 );
        leg -> SetFillStyle( 0 );
        leg -> SetBorderSize( 0 );
        leg -> SetTextSize( 0.03 );
        leg -> AddEntry( h_rho_jet_pt, "rho-subtracted", "l" );
        leg -> AddEntry( h_sub1_jet_pt, "sub1 (seeded) subtracted", "l" );
        leg -> Draw();

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 4: rho vs sub1 jet eta
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetLogy( 0 );

        h_rho_jet_eta -> SetLineColor( kRho );
        h_rho_jet_eta -> SetLineWidth( 2 );
        h_sub1_jet_eta -> SetLineColor( kSub1 );
        h_sub1_jet_eta -> SetLineWidth( 2 );

        h_rho_jet_eta -> GetXaxis() -> SetTitle( "#eta_{jet}" );
        h_rho_jet_eta -> GetYaxis() -> SetTitle( "jets" );
        h_rho_jet_eta -> GetYaxis() -> SetTitleOffset( 1.5 );
        // h_rho_jet_eta -> SetMaximum( 1.7 * std::max( h_rho_jet_eta -> GetMaximum(), h_sub1_jet_eta -> GetMaximum() ) );
        // h_rho_jet_eta -> SetMinimum( 0.0 );

        h_rho_jet_eta -> Scale( 1.0 / h_rho_jet_eta -> Integral() );
        h_sub1_jet_eta -> Scale( 1.0 / h_sub1_jet_eta -> Integral() );

        h_rho_jet_eta -> SetMaximum( 1.7 * std::max( h_rho_jet_eta -> GetMaximum(), h_sub1_jet_eta -> GetMaximum() ) );
        h_rho_jet_eta -> SetMinimum( 0.0 );

        h_rho_jet_eta -> Draw( "hist" );
        h_sub1_jet_eta -> Draw( "hist SAME" );

        auto * leg = new TLegend( 0.28, 0.66, 0.72, 0.82 );
        leg -> SetFillStyle( 0 );
        leg -> SetBorderSize( 0 );
        leg -> SetTextSize( 0.03 );
        leg -> AddEntry( h_rho_jet_eta, "rho-subtracted", "l" );
        leg -> AddEntry( h_sub1_jet_eta, "sub1 (seeded) subtracted", "l" );
        leg -> Draw();

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 5: rho vs sub1 jet phi
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetLogy( 0 );

        h_rho_jet_phi -> SetLineColor( kRho );
        h_rho_jet_phi -> SetLineWidth( 2 );
        h_sub1_jet_phi -> SetLineColor( kSub1 );
        h_sub1_jet_phi -> SetLineWidth( 2 );

        h_rho_jet_phi -> GetXaxis() -> SetTitle( "#phi_{jet}" );
        h_rho_jet_phi -> GetYaxis() -> SetTitle( "jets" );
        h_rho_jet_phi -> GetYaxis() -> SetTitleOffset( 1.5 );

        h_rho_jet_phi -> Scale( 1.0 / h_rho_jet_phi -> Integral() );
        h_sub1_jet_phi -> Scale( 1.0 / h_sub1_jet_phi -> Integral() );

        h_rho_jet_phi -> SetMaximum( 1.7 * std::max( h_rho_jet_phi -> GetMaximum(), h_sub1_jet_phi -> GetMaximum() ) );
        h_rho_jet_phi -> SetMinimum( 0.0 );
        h_rho_jet_phi -> Draw( "hist" );
        h_sub1_jet_phi -> Draw( "hist SAME" );

        auto * leg = new TLegend( 0.28, 0.66, 0.72, 0.82 );
        leg -> SetFillStyle( 0 );
        leg -> SetBorderSize( 0 );
        leg -> SetTextSize( 0.03 );
        leg -> AddEntry( h_rho_jet_phi, "rho-subtracted", "l" );
        leg -> AddEntry( h_sub1_jet_phi, "sub1 (seeded) subtracted", "l" );
        leg -> Draw();

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 6: selected jet multiplicity
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetLogy();

        h_njets_rho -> SetLineColor( kRho );
        h_njets_rho -> SetLineWidth( 2 );
        h_njets_sub1 -> SetLineColor( kSub1 );
        h_njets_sub1 -> SetLineWidth( 2 );

        h_njets_rho -> GetXaxis() -> SetTitle( "N_{jets} (selected)" );
        h_njets_rho -> GetYaxis() -> SetTitle( "events" );
        h_njets_rho -> GetYaxis() -> SetTitleOffset( 1.5 );
        h_njets_rho -> SetMaximum( 3.0 * std::max( h_njets_rho -> GetMaximum(), h_njets_sub1 -> GetMaximum() ) );
        h_njets_rho -> Draw( "hist" );
        h_njets_sub1 -> Draw( "hist SAME" );

        auto * leg = new TLegend( 0.42, 0.72, 0.94, 0.86 );
        leg -> SetFillStyle( 0 );
        leg -> SetBorderSize( 0 );
        leg -> SetTextSize( 0.03 );
        leg -> AddEntry( h_njets_rho, "rho-subtracted", "l" );
        leg -> AddEntry( h_njets_sub1, "sub1 (seeded) subtracted", "l" );
        leg -> Draw();

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 7: matched jet pT correlation
    //--------------------------------------------------------------------
    {
        c -> Clear();
        style_pad( c );
        c -> SetTopMargin( 0.14 );
        c -> SetRightMargin( 0.14 );
        c -> SetLogy( 0 );
        c -> SetLogz();

        h_match_pt_rho_vs_sub1 -> GetXaxis() -> SetTitle( "rho-subtracted jet p_{T} [GeV]" );
        h_match_pt_rho_vs_sub1 -> GetYaxis() -> SetTitle( "sub1-subtracted jet p_{T} [GeV]" );
        h_match_pt_rho_vs_sub1 -> GetYaxis() -> SetTitleOffset( 1.5 );
        h_match_pt_rho_vs_sub1 -> GetXaxis() -> SetRangeUser( 0, 40 );
        h_match_pt_rho_vs_sub1 -> GetYaxis() -> SetRangeUser( 0, 40 );
        h_match_pt_rho_vs_sub1 -> Draw( "COLZ" );

        const double xmax = h_match_pt_rho_vs_sub1 -> GetXaxis() -> GetXmax();
        auto * ldiag = new TLine( 0, 0, xmax, xmax );
        ldiag -> SetLineStyle( 2 );
        ldiag -> SetLineColor( kBlack );
        ldiag -> SetLineWidth( 2 );
        ldiag -> Draw();

        draw_header();
        c -> RedrawAxis();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 8: matched jet delta-pT
    //--------------------------------------------------------------------
    // style + draw a single 1D histogram in a sub-pad of a Divide(ncols,1)
    // canvas: margins, label/title sizes, and the mean/RMS box are all
    // scaled up to compensate for the narrower pad.
    auto style_divided_pad = [&]( TVirtualPad * p, const int ncols )
    {
        p -> cd();
        const double left = ( ncols >= 3 ) ? 0.30 : 0.20;
        p -> SetTopMargin( 0.10 );
        p -> SetLeftMargin( left );
        p -> SetRightMargin( 0.05 );
        p -> SetBottomMargin( 0.16 );
    };

    auto draw_1d_with_zero_line = [&]( TVirtualPad * p, TH1D * h, const char * xtitle, const int ncols )
    {
        style_divided_pad( p, ncols );
        gPad -> SetLogy( 1 );
        const double fsize = ( ncols >= 3 ) ? 0.065 : 0.05;
        const double left = ( ncols >= 3 ) ? 0.30 : 0.20;
        h -> SetLineColor( kRho );
        h -> SetLineWidth( 2 );
        h -> GetXaxis() -> SetTitle( xtitle );
        h -> GetYaxis() -> SetTitle( "pairs" );
        h -> GetXaxis() -> SetLabelSize( fsize );
        h -> GetXaxis() -> SetTitleSize( fsize );
        h -> GetXaxis() -> SetTitleOffset( 1.0 );
        h -> GetYaxis() -> SetLabelSize( fsize );
        h -> GetYaxis() -> SetTitleSize( fsize );
        h -> GetYaxis() -> SetTitleOffset( ( ncols >= 3 ) ? 2.3 : 1.7 );
        h -> Draw( "hist" );
        auto * l0 = new TLine( 0, 0, 0, h -> GetMaximum() );
        l0 -> SetLineStyle( 2 );
        l0 -> SetLineColor( kGray + 2 );
        l0 -> Draw();
        h -> Draw( "hist SAME" );

        TLatex tl;
        tl.SetNDC();
        tl.SetTextSize( fsize * 0.7 );
        tl.DrawLatex( left + 0.02, 0.82, Form( "Mean = %.3f", h -> GetMean() ) );
        tl.DrawLatex( left + 0.02, 0.75, Form( "RMS  = %.3f", h -> GetRMS() ) );
    };

    {
        c -> Clear();
        c -> Divide( 2, 1 );
        draw_1d_with_zero_line( c -> cd(1), h_match_dpt, "p_{T}^{sub1} - p_{T}^{rho} [GeV]", 2 );
        draw_1d_with_zero_line( c -> cd(2), h_match_dpt_relpt, "(p_{T}^{sub1} - p_{T}^{rho}) / p_{T}^{rho}", 2 );
        c -> cd( 0 );
        draw_header();
        c -> Print( mid.c_str() );
    }

    //--------------------------------------------------------------------
    // page 9: matched jet dR, deta, dphi
    //--------------------------------------------------------------------
    {
        c -> Clear();
        c -> Divide( 3, 1 );

        style_divided_pad( c -> cd(1), 3 );
        gPad -> SetLogy( 1 );
        h_match_dr -> SetLineColor( kRho );
        h_match_dr -> SetLineWidth( 2 );
        h_match_dr -> GetXaxis() -> SetTitle( "#Delta R(rho, sub1)" );
        h_match_dr -> GetYaxis() -> SetTitle( "pairs" );
        h_match_dr -> GetXaxis() -> SetLabelSize( 0.065 );
        h_match_dr -> GetXaxis() -> SetTitleSize( 0.065 );
        h_match_dr -> GetXaxis() -> SetTitleOffset( 1.0 );
        h_match_dr -> GetYaxis() -> SetLabelSize( 0.065 );
        h_match_dr -> GetYaxis() -> SetTitleSize( 0.065 );
        h_match_dr -> GetYaxis() -> SetTitleOffset( 2.3 );
        h_match_dr -> Draw( "hist" );

        draw_1d_with_zero_line( c -> cd(2), h_match_deta, "#eta^{sub1} - #eta^{rho}", 3 );
        draw_1d_with_zero_line( c -> cd(3), h_match_dphi, "#phi^{sub1} - #phi^{rho}", 3 );

        c -> cd( 0 );
        draw_header();
        c -> Print( close.c_str() );
    }

    f -> Close();

    std::cout << "Wrote " << outfile << std::endl;

    return 0;
}

#endif // _DRAW_PLOTS_C_
