
#ifndef _PLOT_REFERENCE_CALIB_C_
#define _PLOT_REFERENCE_CALIB_C_

#include <cdbobjects/CDBTTree.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TString.h>
#include <TSystemDirectory.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libcdbobjects.so )

// Companion plots for find_reference_calib.C's pick of run 54590 as the
// default/reference rho eta-shape calibration: (1) per-run RMS distance
// from the population mean, sorted, showing 54590 sitting at the minimum
// and the two low-statistics outliers (54853/54859) sitting far above
// everyone else; (2) per-layer <weight> vs ieta (averaged over zvtx/mbdQ
// bins), population mean +/- RMS band vs run 54590's own profile, plus the
// two outlier runs for contrast.
int plot_reference_calib(
    const std::string & calib_dir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/calibs",
    const int highlight_run = 54590,
    const int compare_run = 54912,
    const std::string & outdir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib"
)
{
    static const int kNCalo = 3;
    static const char * calo_names[kNCalo] = { "cemc", "hcalin", "hcalout" };
    static const char * calib_field_names[kNCalo] = { "w_cemc", "w_hcalin", "w_hcalout" };

    std::vector<int> runs;
    {
        TSystemDirectory dir( "calibdir", calib_dir.c_str() );
        TList * files = dir.GetListOfFiles();
        for ( auto * obj : *files )
        {
            TString name = obj->GetName();
            if ( !name.BeginsWith( "rho_calib_" ) || !name.EndsWith( ".root" ) || name.Contains( "_qa" ) ) continue;
            TString numpart = name;
            numpart.ReplaceAll( "rho_calib_", "" );
            numpart.ReplaceAll( ".root", "" );
            if ( numpart.IsDigit() ) runs.push_back( numpart.Atoi() );
        }
    }
    std::sort( runs.begin(), runs.end() );
    const size_t nruns = runs.size();
    std::cout << "Found " << nruns << " runs" << std::endl;

    int NETA = -1, n_zbins = -1, n_mbdbins = -1, nchannels = -1;
    std::vector<std::vector<float>> weights; // [run][channel*3+layer]
    std::vector<int> n_fallback;
    weights.reserve( nruns );
    n_fallback.reserve( nruns );

    for ( size_t ir = 0; ir < runs.size(); ++ir )
    {
        auto * cdbttree = new CDBTTree( Form( "%s/rho_calib_%d.root", calib_dir.c_str(), runs[ir] ) );
        cdbttree->LoadCalibrations();
        if ( NETA < 0 )
        {
            NETA = cdbttree->GetSingleIntValue( "n_eta" );
            n_zbins = cdbttree->GetSingleIntValue( "n_zvtx_bins" );
            n_mbdbins = cdbttree->GetSingleIntValue( "n_mbdQ_bins" );
            nchannels = n_zbins * n_mbdbins * NETA;
        }
        std::vector<float> v( nchannels * kNCalo, 1.0f );
        int nfb = 0;
        for ( int ch = 0; ch < nchannels; ++ch )
            for ( int ic = 0; ic < kNCalo; ++ic )
            {
                const float w = cdbttree->GetFloatValue( ch, calib_field_names[ic] );
                v[ch * kNCalo + ic] = w;
                if ( !( w > 0.0f ) || w == 1.0f ) ++nfb;
            }
        weights.push_back( std::move( v ) );
        n_fallback.push_back( nfb );
        delete cdbttree;
    }

    const int ntotal = nchannels * kNCalo;
    std::vector<double> mean( ntotal, 0.0 );
    for ( size_t ir = 0; ir < nruns; ++ir ) for ( int k = 0; k < ntotal; ++k ) mean[k] += weights[ir][k];
    for ( int k = 0; k < ntotal; ++k ) mean[k] /= nruns;

    std::vector<double> dist( nruns, 0.0 );
    for ( size_t ir = 0; ir < nruns; ++ir )
    {
        double s = 0.0;
        for ( int k = 0; k < ntotal; ++k ) { const double d = weights[ir][k] - mean[k]; s += d * d; }
        dist[ir] = std::sqrt( s / ntotal );
    }

    // ============ Plot 1: sorted RMS distance, highlight best + worst ============
    std::vector<size_t> order( nruns );
    for ( size_t i = 0; i < nruns; ++i ) order[i] = i;
    std::sort( order.begin(), order.end(), [&]( size_t a, size_t b ) { return dist[a] < dist[b]; } );

    auto * gAll = new TGraph( nruns );
    auto * gHighlight = new TGraph( 0 );
    int highlight_rank = -1;
    for ( size_t i = 0; i < nruns; ++i )
    {
        gAll->SetPoint( i, i, dist[order[i]] );
        if ( runs[order[i]] == highlight_run )
        {
            highlight_rank = i;
            gHighlight->SetPoint( 0, i, dist[order[i]] );
        }
    }

    auto * c1 = new TCanvas( "c1", "rms distance", 1000, 700 );
    c1->SetLogy();
    gAll->SetTitle( Form( "run-by-run rho calibration: RMS distance from population mean (%zu runs);rank (sorted, most-typical first);RMS(w_{run} - w_{mean})", nruns ) );
    gAll->SetMarkerStyle( 20 );
    gAll->SetMarkerSize( 0.9 );
    gAll->SetMarkerColor( kAzure + 2 );
    gAll->Draw( "AP" );
    gHighlight->SetMarkerStyle( 29 );
    gHighlight->SetMarkerSize( 2.5 );
    gHighlight->SetMarkerColor( kRed + 1 );
    gHighlight->Draw( "P SAME" );

    auto * leg1 = new TLegend( 0.15, 0.72, 0.55, 0.85 );
    leg1->AddEntry( gAll, "all runs", "p" );
    leg1->AddEntry( gHighlight, Form( "run %d (rank %d = best)", highlight_run, highlight_rank + 1 ), "p" );
    leg1->Draw();
    c1->SaveAs( Form( "%s/reference_calib_rms_distance.png", outdir.c_str() ) );
    c1->SaveAs( Form( "%s/reference_calib_rms_distance.pdf", outdir.c_str() ) );

    // ============ Plot 2: <weight> vs ieta per layer, mean+-RMS band vs 54590 vs outliers ============
    // find the 2 worst (outlier) runs among those actually included
    int compare_rank = -1;
    for ( size_t i = 0; i < nruns; ++i ) if ( runs[order[i]] == compare_run ) { compare_rank = i; break; }
    const int outlier1 = runs[order[nruns - 1]];
    const int outlier2 = runs[order[nruns - 2]];
    int idx_highlight = -1, idx_out1 = -1, idx_out2 = -1, idx_compare = -1;
    for ( size_t ir = 0; ir < nruns; ++ir )
    {
        if ( runs[ir] == highlight_run ) idx_highlight = ir;
        if ( runs[ir] == outlier1 ) idx_out1 = ir;
        if ( runs[ir] == outlier2 ) idx_out2 = ir;
        if ( runs[ir] == compare_run ) idx_compare = ir;
    }

    auto * c2 = new TCanvas( "c2", "eta shape comparison", 1800, 600 );
    c2->Divide( 3, 1 );
    auto * c3 = new TCanvas( "c3", "reference vs compare", 1800, 600 );
    c3->Divide( 3, 1 );

    auto channel_of = [&]( int ieta, int izbin, int imbd ) { return izbin * ( n_mbdbins * NETA ) + imbd * NETA + ieta; };

    for ( int ic = 0; ic < kNCalo; ++ic )
    {
        // per-run per-ieta average over (izbin,imbd), skipping exact-1.0 fallback entries
        auto avg_over_eta = [&]( const std::vector<float> & v ) {
            std::vector<double> out( NETA, 0.0 );
            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                double s = 0.0; int n = 0;
                for ( int iz = 0; iz < n_zbins; ++iz )
                    for ( int imbd = 0; imbd < n_mbdbins; ++imbd )
                    {
                        const int ch = channel_of( ieta, iz, imbd );
                        const float w = v[ch * kNCalo + ic];
                        if ( w == 1.0f ) continue; // skip fallback
                        s += w; ++n;
                    }
                out[ieta] = ( n > 0 ) ? s / n : 1.0;
            }
            return out;
        };

        // population mean/RMS of the per-run eta-profile (across runs)
        std::vector<std::vector<double>> per_run_profile( nruns );
        for ( size_t ir = 0; ir < nruns; ++ir ) per_run_profile[ir] = avg_over_eta( weights[ir] );

        auto * gBand = new TGraphErrors( NETA );
        auto * gRef = new TGraph( NETA );
        auto * gOut1 = new TGraph( NETA );
        auto * gOut2 = new TGraph( NETA );
        for ( int ieta = 0; ieta < NETA; ++ieta )
        {
            double s = 0.0;
            for ( size_t ir = 0; ir < nruns; ++ir ) s += per_run_profile[ir][ieta];
            const double m = s / nruns;
            double s2 = 0.0;
            for ( size_t ir = 0; ir < nruns; ++ir ) { const double d = per_run_profile[ir][ieta] - m; s2 += d * d; }
            const double rms = std::sqrt( s2 / nruns );
            gBand->SetPoint( ieta, ieta, m );
            gBand->SetPointError( ieta, 0.0, rms );
            gRef->SetPoint( ieta, ieta, per_run_profile[idx_highlight][ieta] );
            gOut1->SetPoint( ieta, ieta, per_run_profile[idx_out1][ieta] );
            gOut2->SetPoint( ieta, ieta, per_run_profile[idx_out2][ieta] );
        }

        auto * gCmp = new TGraph( NETA );
        for ( int ieta = 0; ieta < NETA; ++ieta ) gCmp->SetPoint( ieta, ieta, per_run_profile[idx_compare][ieta] );

        // per-panel auto range from the population band + reference run (outliers may clip -- that's fine, they're shown only for contrast)
        double ymin = 1e18, ymax = -1e18;
        for ( int ieta = 0; ieta < NETA; ++ieta )
        {
            double x, y, ey;
            gBand->GetPoint( ieta, x, y ); ey = gBand->GetErrorY( ieta );
            ymin = std::min( { ymin, y - ey, per_run_profile[idx_highlight][ieta] } );
            ymax = std::max( { ymax, y + ey, per_run_profile[idx_highlight][ieta] } );
        }
        const double pad = 0.15 * ( ymax - ymin );
        ymin -= pad; ymax += pad;

        c2->cd( ic + 1 );
        gBand->SetTitle( Form( "%s;i_{#eta};<weight> (avg over zvtx/mbdQ bins)", calo_names[ic] ) );
        gBand->SetFillColorAlpha( kGray + 1, 0.5 );
        gBand->SetFillStyle( 1001 );
        gBand->SetLineColor( kGray + 2 );
        gBand->GetYaxis()->SetRangeUser( ymin, ymax );
        gBand->Draw( "A3" );
        auto * gBandLine = (TGraph*) gBand->Clone();
        gBandLine->SetLineColor( kBlack );
        gBandLine->SetLineWidth( 2 );
        gBandLine->Draw( "L SAME" );

        gRef->SetLineColor( kRed + 1 );
        gRef->SetLineWidth( 3 );
        gRef->SetMarkerColor( kRed + 1 );
        gRef->SetMarkerStyle( 20 );
        gRef->Draw( "LP SAME" );

        gOut1->SetLineColor( kAzure + 2 );
        gOut1->SetLineStyle( 2 );
        gOut1->SetLineWidth( 2 );
        gOut1->Draw( "L SAME" );

        gOut2->SetLineColor( kSpring + 4 );
        gOut2->SetLineStyle( 2 );
        gOut2->SetLineWidth( 2 );
        gOut2->Draw( "L SAME" );

        if ( ic == 0 )
        {
            auto * leg2 = new TLegend( 0.13, 0.68, 0.62, 0.89 );
            leg2->AddEntry( gBandLine, "population mean (122 runs)", "l" );
            leg2->AddEntry( (TObject*) gBand, "#pm1 RMS band", "f" );
            leg2->AddEntry( gRef, Form( "run %d (reference)", highlight_run ), "lp" );
            leg2->AddEntry( gOut1, Form( "run %d (outlier, low-stat)", outlier1 ), "l" );
            leg2->AddEntry( gOut2, Form( "run %d (outlier, low-stat)", outlier2 ), "l" );
            leg2->Draw();
        }

        // -- c3: reference (54590) vs compare (54912) vs population mean, no outliers --
        double ymin3 = 1e18, ymax3 = -1e18;
        for ( int ieta = 0; ieta < NETA; ++ieta )
        {
            double x, y, ey;
            gBand->GetPoint( ieta, x, y ); ey = gBand->GetErrorY( ieta );
            ymin3 = std::min( { ymin3, y - ey, per_run_profile[idx_highlight][ieta], per_run_profile[idx_compare][ieta] } );
            ymax3 = std::max( { ymax3, y + ey, per_run_profile[idx_highlight][ieta], per_run_profile[idx_compare][ieta] } );
        }
        const double pad3 = 0.15 * ( ymax3 - ymin3 );
        ymin3 -= pad3; ymax3 += pad3;

        c3->cd( ic + 1 );
        auto * gBand3 = (TGraphErrors*) gBand->Clone();
        gBand3->SetTitle( Form( "%s;i_{#eta};<weight> (avg over zvtx/mbdQ bins)", calo_names[ic] ) );
        gBand3->GetYaxis()->SetRangeUser( ymin3, ymax3 );
        gBand3->Draw( "A3" );
        auto * gBandLine3 = (TGraph*) gBand3->Clone();
        gBandLine3->SetLineColor( kBlack );
        gBandLine3->SetLineWidth( 2 );
        gBandLine3->Draw( "L SAME" );

        gRef->Draw( "LP SAME" );

        gCmp->SetLineColor( kAzure + 2 );
        gCmp->SetLineWidth( 3 );
        gCmp->SetMarkerColor( kAzure + 2 );
        gCmp->SetMarkerStyle( 21 );
        gCmp->Draw( "LP SAME" );

        if ( ic == 0 )
        {
            auto * leg3 = new TLegend( 0.13, 0.68, 0.62, 0.89 );
            leg3->AddEntry( gBandLine3, "population mean (122 runs)", "l" );
            leg3->AddEntry( (TObject*) gBand3, "#pm1 RMS band", "f" );
            leg3->AddEntry( gRef, Form( "run %d (rank 1, reference)", highlight_run ), "lp" );
            leg3->AddEntry( gCmp, Form( "run %d (rank %d)", compare_run, compare_rank + 1 ), "lp" );
            leg3->Draw();
        }
    }
    c2->SaveAs( Form( "%s/reference_calib_eta_shape.png", outdir.c_str() ) );
    c2->SaveAs( Form( "%s/reference_calib_eta_shape.pdf", outdir.c_str() ) );
    c3->SaveAs( Form( "%s/reference_calib_vs_compare.png", outdir.c_str() ) );
    c3->SaveAs( Form( "%s/reference_calib_vs_compare.pdf", outdir.c_str() ) );

    std::cout << "Wrote reference_calib_rms_distance.{png,pdf}, reference_calib_eta_shape.{png,pdf}, "
              << "and reference_calib_vs_compare.{png,pdf} in " << outdir << std::endl;
    return 0;
}

#endif
