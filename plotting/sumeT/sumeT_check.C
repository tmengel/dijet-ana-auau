#ifndef _SUMET_CHECK_C_
#define _SUMET_CHECK_C_

#include <sPhenixStyle.C>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

#include <iostream>
#include <string>

// quick check of the event-level "sumeT" branch (MBD sum eT) between the
// unscaled and scaled hijing samples.
int sumeT_check(
    const std::string & infile_unscaled = "rootfiles/08_23_2026_v001/jet30_hijing_unscaled_all.root",
    const std::string & infile_scaled   = "rootfiles/08_23_2026_v001/jet30_hijing_scaled_all.root",
    const std::string & infile_data     = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/output/54912/rho_calib_54912.root",
    const std::string & outfile         = "sumeT_check.pdf",
    const std::string & label           = "Hijing, jet30"
)
{
    SetsPhenixStyle();
    gStyle -> SetOptStat( 0 );

    auto * f_unscaled = TFile::Open( infile_unscaled.c_str() );
    if ( !f_unscaled || f_unscaled -> IsZombie() )
    {
        std::cerr << "Error: could not open " << infile_unscaled << std::endl;
        return -1;
    }
    auto * f_scaled = TFile::Open( infile_scaled.c_str() );
    if ( !f_scaled || f_scaled -> IsZombie() )
    {
        std::cerr << "Error: could not open " << infile_scaled << std::endl;
        return -1;
    }
    auto * f_data = TFile::Open( infile_data.c_str() );
    if ( !f_data || f_data -> IsZombie() )
    {
        std::cerr << "Error: could not open " << infile_data << std::endl;
        return -1;
    }

    auto * t_unscaled = (TTree*) f_unscaled -> Get( "T" );
    auto * t_scaled   = (TTree*) f_scaled   -> Get( "T" );
    auto * t_data     = (TTree*) f_data     -> Get( "T" );
    if ( !t_unscaled || !t_scaled || !t_data )
    {
        std::cerr << "Error: no tree \"T\" found in one of the input files" << std::endl;
        return -1;
    }

    const int nbins   = 100;
    const double xmin = 0.0;
    const double xmax = 2000.0;

    auto * h_unscaled = new TH1F( "h_sumeT_unscaled", ";sum e_{T} (MBD) [GeV];events", nbins, xmin, xmax );
    auto * h_scaled   = new TH1F( "h_sumeT_scaled",   ";sum e_{T} (MBD) [GeV];events", nbins, xmin, xmax );
    auto * h_data     = new TH1F( "h_sumeT_data",     ";sum e_{T} (MBD) [GeV];events", nbins, xmin, xmax );
    h_unscaled -> Sumw2();
    h_scaled   -> Sumw2();
    h_data     -> Sumw2();

    // sumeT is an event-level quantity duplicated across every jet row
    // sharing the same event_id, so only fill once per distinct event_id.
    auto fill_per_event = [] ( TTree * t, TH1F * h )
    {
        int event_id = 0;
        float sumeT  = 0.0;
        t -> SetBranchStatus( "*", 0 );
        t -> SetBranchStatus( "event_id", 1 );
        t -> SetBranchStatus( "sumeT", 1 );
        t -> SetBranchAddress( "event_id", &event_id );
        t -> SetBranchAddress( "sumeT", &sumeT );

        int last_event_id = -1;
        const Long64_t n = t -> GetEntries();
        for ( Long64_t i = 0; i < n; ++i )
        {
            t -> GetEntry( i );
            if ( event_id == last_event_id ) continue;
            last_event_id = event_id;
            h -> Fill( sumeT );
        }
        t -> ResetBranchAddresses();
        t -> SetBranchStatus( "*", 1 );
    };

    auto fill_per_event_scaled = [] ( TTree * t, TH1F * h , const float scale_factor = 1.11 )
    {
        int event_id = 0;
        float sumeT  = 0.0;
        t -> SetBranchStatus( "*", 0 );
        t -> SetBranchStatus( "event_id", 1 );
        t -> SetBranchStatus( "sumeT", 1 );
        t -> SetBranchAddress( "event_id", &event_id );
        t -> SetBranchAddress( "sumeT", &sumeT );

        int last_event_id = -1;
        const Long64_t n = t -> GetEntries();
        for ( Long64_t i = 0; i < n; ++i )
        {
            t -> GetEntry( i );
            if ( event_id == last_event_id ) continue;
            last_event_id = event_id;
            h -> Fill( sumeT * scale_factor );
        }
        t -> ResetBranchAddresses();
        t -> SetBranchStatus( "*", 1 );
    };

    fill_per_event_scaled( t_unscaled, h_unscaled, 1.00 );
    fill_per_event( t_scaled, h_scaled );

    // real data (rho_calib output, run 54912): already one row per event, and
    // sumeT is split by calo layer ([0]=cemc, [1]=hcalin, [2]=hcalout), so
    // sum the three layers to get the same total sum eT as the sim branch.
    {
        float sumeT_layers[3] = { 0.0, 0.0, 0.0 };
        t_data -> SetBranchStatus( "*", 0 );
        t_data -> SetBranchStatus( "sumeT", 1 );
        t_data -> SetBranchAddress( "sumeT", sumeT_layers );

        const Long64_t n = t_data -> GetEntries();
        for ( Long64_t i = 0; i < n; ++i )
        {
            t_data -> GetEntry( i );
            h_data -> Fill( sumeT_layers[0] + sumeT_layers[1] + sumeT_layers[2] );
        }
        t_data -> ResetBranchAddresses();
        t_data -> SetBranchStatus( "*", 1 );
    }

    std::cout << "unscaled: " << h_unscaled -> GetEntries() << " unique events" << std::endl;
    std::cout << "scaled:   " << h_scaled   -> GetEntries() << " unique events" << std::endl;
    std::cout << "data:     " << h_data     -> GetEntries() << " events" << std::endl;

    if ( h_unscaled -> Integral() > 0 ) h_unscaled -> Scale( 1.0 / h_unscaled -> Integral() );
    if ( h_scaled   -> Integral() > 0 ) h_scaled   -> Scale( 1.0 / h_scaled   -> Integral() );
    if ( h_data     -> Integral() > 0 ) h_data     -> Scale( 1.0 / h_data     -> Integral() );

    h_unscaled -> SetLineColor( kBlack );
    h_unscaled -> SetLineWidth( 2 );
    h_unscaled -> SetMarkerStyle( 20 );
    h_unscaled -> SetMarkerColor( kBlack );

    h_scaled -> SetLineColor( kAzure + 2 );
    h_scaled -> SetLineWidth( 2 );
    h_scaled -> SetMarkerStyle( 21 );
    h_scaled -> SetMarkerColor( kAzure + 2 );

    h_data -> SetLineColor( kRed + 1 );
    h_data -> SetLineWidth( 2 );
    h_data -> SetMarkerStyle( 22 );
    h_data -> SetMarkerColor( kRed + 1 );

    const double ymax = 1.35 * std::max( { h_unscaled -> GetMaximum(), h_scaled -> GetMaximum(), h_data -> GetMaximum() } );

    // ratio of each simulation histogram to data
    auto * h_ratio_unscaled = (TH1F*) h_unscaled -> Clone( "h_sumeT_ratio_unscaled" );
    auto * h_ratio_scaled   = (TH1F*) h_scaled   -> Clone( "h_sumeT_ratio_scaled" );
    h_ratio_unscaled -> Divide( h_data );
    h_ratio_scaled   -> Divide( h_data );

    auto * c = new TCanvas( "c_sumeT_check", "c_sumeT_check", 600, 650 );

    auto * p_main = new TPad( "p_main", "p_main", 0.0, 0.30, 1.0, 1.0 );
    p_main -> SetTopMargin( 0.14 );
    p_main -> SetLeftMargin( 0.14 );
    p_main -> SetRightMargin( 0.05 );
    p_main -> SetBottomMargin( 0.02 );
    p_main -> SetLogy( 1 );
    p_main -> Draw();

    auto * p_ratio = new TPad( "p_ratio", "p_ratio", 0.0, 0.0, 1.0, 0.30 );
    p_ratio -> SetTopMargin( 0.02 );
    p_ratio -> SetLeftMargin( 0.14 );
    p_ratio -> SetRightMargin( 0.05 );
    p_ratio -> SetBottomMargin( 0.35 );
    p_ratio -> Draw();

    p_main -> cd();

    auto * hframe = p_main -> DrawFrame( xmin, 1e-8, xmax, ymax, ";;fraction of events" );
    hframe -> GetYaxis() -> SetTitleOffset( 1.5 );
    hframe -> GetXaxis() -> SetLabelSize( 0 );

    h_unscaled -> Draw( "HIST E SAME" );
    h_scaled   -> Draw( "HIST E SAME" );
    h_data     -> Draw( "HIST E SAME" );

    auto * leg = new TLegend( 0.55, 0.66, 0.93, 0.83 );
    leg -> SetBorderSize( 0 );
    leg -> SetFillStyle( 0 );
    leg -> SetTextSize( 0.03 );
    leg -> AddEntry( h_unscaled, "unscaled", "l" );
    leg -> AddEntry( h_scaled, "scaled", "l" );
    leg -> AddEntry( h_data, "data (run 54912)", "l" );
    leg -> Draw();

    TLatex tl_sphenix;
    tl_sphenix.SetNDC();
    tl_sphenix.SetTextColor( kBlack );
    tl_sphenix.SetTextSize( 0.045 );
    tl_sphenix.DrawLatex( 0.14, 0.955, "#it{#bf{sPHENIX}} Simulation" );

    if ( !label.empty() )
    {
        TLatex tl_label;
        tl_label.SetNDC();
        tl_label.SetTextAlign( 31 );
        tl_label.SetTextSize( 0.028 );
        tl_label.DrawLatex( 0.95, 0.955, label.c_str() );
    }

    p_main -> RedrawAxis();

    p_ratio -> cd();

    auto * hratio_frame = p_ratio -> DrawFrame( xmin, 0.0, xmax, 2.0, ";sum e_{T} (MBD) [GeV];sim / data" );
    hratio_frame -> GetXaxis() -> SetTitleSize( 0.11 );
    hratio_frame -> GetXaxis() -> SetLabelSize( 0.09 );
    hratio_frame -> GetYaxis() -> SetTitleSize( 0.11 );
    hratio_frame -> GetYaxis() -> SetLabelSize( 0.09 );
    hratio_frame -> GetYaxis() -> SetTitleOffset( 0.55 );
    hratio_frame -> GetYaxis() -> SetNdivisions( 505 );

    TLine unity_line( xmin, 1.0, xmax, 1.0 );
    unity_line.SetLineStyle( 2 );
    unity_line.SetLineColor( kGray + 2 );
    unity_line.Draw();

    h_ratio_unscaled -> Draw( "E SAME" );
    h_ratio_scaled   -> Draw( "E SAME" );

    p_ratio -> RedrawAxis();

    c -> SaveAs( outfile.c_str() );

    std::cout << "Wrote " << outfile << std::endl;

    return 0;
}

#endif
