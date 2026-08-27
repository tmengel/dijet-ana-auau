
#ifndef _FINAL_CROSSCHECK_LOCALRHO_C_
#define _FINAL_CROSSCHECK_LOCALRHO_C_

#include <myana/AnaUtils.h>
#include <cdbobjects/CDBTTree.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )
R__LOAD_LIBRARY( libcdbobjects.so )

// Attempt at a width-minimizing background estimator that is still "a rho
// method" -- i.e. still a robust local-density estimator subtracted
// directly from the tower, just localized much further than the standard
// single flat rho per calo layer per event.
//
// Standard rho: one MEDIAN-based density for the whole calo layer, per
// event -- captures event-to-event multiplicity/centrality variation but
// nothing about eta shape (hence needing the separately-derived w(ieta)
// calibration to fix bias, per the earlier cross-check).
//
// "local rho" (this file): the MEDIAN of that SAME ring's (up to 64) good
// tower energies, computed independently PER (event, layer, eta ring). This
// keeps the same jet-robust median philosophy as rho (a jet landing in a
// few of the 64 phi slots barely moves the median), but is local enough to
// self-adjust to whatever the true instantaneous background level is in
// that specific ring/event -- including the eta-dependent geometric effects
// that otherwise require the external w(ieta) calibration -- with NO cosh(eta)
// or w(ieta) needed at all: the median IS the predicted tower energy
// directly, in the same units, for that ring.
//
// This differs from the earlier failed "event-by-event" attempt in two
// structural ways: (1) it's a per-RING statistic (24 independent estimates
// per event, not one whole-event scalar diluting/mixing all rings
// together), and (2) it's a MEDIAN (robust to a single jet-hit tower) used
// ADDITIVELY (predicted = median, residual = raw - median), not a MEAN used
// as a fragile multiplicative denominator.
int final_crosscheck_localrho(
    const int run = 54590,
    const int nfiles_max = 2,
    const std::string & anatree_dir = "/sphenix/tg/tg01/jets/tmengel/jstgtf03/data/v001_20260720/anatree",
    const std::string & calib_dir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/calibs",
    const std::string & outdir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib",
    const float max_abs_zvrtx = 60.0
)
{
    static const int kNCalo = 3;
    static const int NETA = 24;
    static const int NPHI = 64;
    static const char * calo_names[kNCalo] = { "cemc", "hcalin", "hcalout" };
    static const char * calib_field_names[kNCalo] = { "w_cemc", "w_hcalin", "w_hcalout" };
    const AnaUtils::CaloType calo_types[kNCalo] = { AnaUtils::CEMC, AnaUtils::HCALIN, AnaUtils::HCALOUT };

    auto * cdbttree = new CDBTTree( Form( "%s/rho_calib_%d.root", calib_dir.c_str(), run ) );
    cdbttree->LoadCalibrations();
    const int n_zbins = cdbttree->GetSingleIntValue( "n_zvtx_bins" );
    const int n_mbdbins = cdbttree->GetSingleIntValue( "n_mbdQ_bins" );
    std::vector<float> zvtx_edges( n_zbins + 1 );
    for ( int i = 0; i <= n_zbins; ++i ) zvtx_edges[i] = cdbttree->GetSingleFloatValue( Form( "zvtx_edge_%d", i ) );
    std::vector<float> mbdQ_edges( n_mbdbins + 1 );
    for ( int i = 0; i <= n_mbdbins; ++i ) mbdQ_edges[i] = cdbttree->GetSingleFloatValue( Form( "mbdQ_edge_%d", i ) );

    auto find_bin = []( const float val, const std::vector<float> & edges ) -> int
    {
        if ( edges.size() < 2 || std::isnan(val) || val < edges.front() || val >= edges.back() ) return -1;
        for ( size_t i = 0; i + 1 < edges.size(); ++i ) if ( val >= edges[i] && val < edges[i + 1] ) return static_cast<int>(i);
        return -1;
    };
    auto encode_channel = [&]( const int ieta, const int izbin, const int imbd ) -> int
    { return izbin * ( n_mbdbins * NETA ) + imbd * NETA + ieta; };
    auto w_of = [&]( int ic, int ieta, int izbin, int imbd ) -> float
    { return cdbttree->GetFloatValue( encode_channel( ieta, izbin, imbd ), calib_field_names[ic] ); };

    auto * ch = new TChain( "T" );
    int nfiles_added = 0;
    for ( int seg = 0; seg < 999 && nfiles_added < nfiles_max; ++seg )
    {
        const std::string fn = Form( "%s/anatree_run2auau_pro001_pcdb001_v001-%08d-%05d.root", anatree_dir.c_str(), run, seg );
        if ( gSystem->AccessPathName( fn.c_str() ) ) continue;
        ch->Add( fn.c_str() );
        ++nfiles_added;
    }
    std::cout << "Chained " << nfiles_added << " anatree file(s) for run " << run << std::endl;

    ch->SetBranchStatus( "*", false );
    int is_minbias = 0;
    float zvrtx = 0, mbd_q_N = 0, mbd_q_S = 0;
    float rho_vals[kNCalo];
    float tower_E[kNCalo][NETA][NPHI];
    int tower_isgood[kNCalo][NETA][NPHI];
    float sub2_ue[kNCalo][NETA];
    auto enable = [&]( const char * name, void * addr ) { ch->SetBranchStatus( name, true ); ch->SetBranchAddress( name, addr ); };
    enable( "is_minbias", &is_minbias ); enable( "zvrtx", &zvrtx );
    enable( "mbd_q_N", &mbd_q_N ); enable( "mbd_q_S", &mbd_q_S );
    enable( "rho_vals", rho_vals ); enable( "tower_E", tower_E );
    enable( "tower_isgood", tower_isgood ); enable( "sub2_towerbkgd_ue", sub2_ue );

    const Long64_t nentries = ch->GetEntries();
    std::cout << "Total entries in chain: " << nentries << std::endl;

    static const int kNMethod = 4; // 0=uncal 1=cal 2=sub1 3=local_rho(median)
    const char * method_names[kNMethod] = { "uncal (rho)", "cal (rho, static w)", "sub1 (iterative UE)", "local rho (per-ring median)" };

    struct Acc { double sum = 0, sumsq = 0; long n = 0;
        void add( double v ) { sum += v; sumsq += v * v; ++n; }
        double mean() const { return n > 0 ? sum / n : 0.0; }
        double rms() const { return n > 0 ? std::sqrt( std::max( 0.0, sumsq / n - mean() * mean() ) ) : 0.0; } };

    std::vector<std::vector<Acc>> accRaw( kNCalo, std::vector<Acc>( NETA ) );
    std::vector<std::vector<std::vector<Acc>>> accResid( kNMethod, std::vector<std::vector<Acc>>( kNCalo, std::vector<Acc>( NETA ) ) );

    std::vector<float> ring_vals;
    ring_vals.reserve( NPHI );

    long n_events_used = 0;
    for ( Long64_t i = 0; i < nentries; ++i )
    {
        ch->GetEntry( i );
        if ( nentries >= 10 && i % ( nentries / 10 ) == 0 && i > 0 )
            std::cout << "processing " << i << " / " << nentries << std::endl;

        if ( !is_minbias ) continue;
        if ( std::fabs( zvrtx ) >= max_abs_zvrtx ) continue;
        const float mbdQ = mbd_q_N + mbd_q_S;
        const int izbin = find_bin( zvrtx, zvtx_edges );
        const int imbd = find_bin( mbdQ, mbdQ_edges );
        if ( izbin < 0 || imbd < 0 ) continue;
        ++n_events_used;

        for ( int ic = 0; ic < kNCalo; ++ic )
        {
            if ( rho_vals[ic] <= 0.0 ) continue;
            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                const double eta = AnaUtils::get_corrected_calo_eta( calo_types[ic], ieta, zvrtx );
                const double coshE = std::cosh( eta );
                const double w = w_of( ic, ieta, izbin, imbd );
                const double pred_uncal = rho_vals[ic] * coshE;
                const double pred_cal = pred_uncal * w;
                const double pred_sub1 = sub2_ue[ic][ieta];

                // -- local rho: median of this ring's good towers, this event --
                ring_vals.clear();
                for ( int iphi = 0; iphi < NPHI; ++iphi )
                    if ( tower_isgood[ic][ieta][iphi] ) ring_vals.push_back( tower_E[ic][ieta][iphi] );
                double pred_local = 0.0;
                if ( !ring_vals.empty() )
                {
                    std::sort( ring_vals.begin(), ring_vals.end() );
                    const size_t n = ring_vals.size();
                    pred_local = ( n % 2 == 1 ) ? ring_vals[n / 2] : 0.5 * ( ring_vals[n / 2 - 1] + ring_vals[n / 2] );
                }

                for ( int iphi = 0; iphi < NPHI; ++iphi )
                {
                    if ( !tower_isgood[ic][ieta][iphi] ) continue;
                    const double E_raw = tower_E[ic][ieta][iphi];
                    accRaw[ic][ieta].add( E_raw );
                    accResid[0][ic][ieta].add( E_raw - pred_uncal );
                    accResid[1][ic][ieta].add( E_raw - pred_cal );
                    accResid[2][ic][ieta].add( E_raw - pred_sub1 );
                    accResid[3][ic][ieta].add( E_raw - pred_local );
                }
            }
        }
    }
    std::cout << n_events_used << " / " << nentries << " events used." << std::endl;

    auto fixpad = []( double left = 0.17 ) { gPad->SetLeftMargin( left ); gPad->SetRightMargin( 0.04 ); };
    auto fixaxisG = []( TGraph * g ) { g->GetYaxis()->SetTitleOffset( 1.6 ); g->GetYaxis()->SetLabelSize( 0.035 ); };
    const int plot_start = 1; // drop uncal from the plots (still in console summary)

    auto * cBias = new TCanvas( "c_bias_lr", "fractional bias (local rho)", 1800, 700 );
    cBias->Divide( 3, 1 );
    auto * cWidth = new TCanvas( "c_width_lr", "fractional width (local rho)", 1800, 700 );
    cWidth->Divide( 3, 1 );
    int colors[kNMethod] = { kAzure + 2, kRed + 1, (int) kSpring + 4, kMagenta + 1 };

    for ( int ic = 0; ic < kNCalo; ++ic )
    {
        std::vector<TGraph*> gBias( kNMethod ), gWidth( kNMethod );
        double biasMin = 1e18, biasMax = -1e18, widthMax = 0.0;
        for ( int m = 0; m < kNMethod; ++m )
        {
            gBias[m] = new TGraph( NETA );
            gWidth[m] = new TGraph( NETA );
            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                const double eref = accRaw[ic][ieta].mean();
                const double b = ( eref != 0 ) ? accResid[m][ic][ieta].mean() / eref : 0.0;
                const double wdt = ( eref != 0 ) ? accResid[m][ic][ieta].rms() / eref : 0.0;
                gBias[m]->SetPoint( ieta, ieta, b );
                gWidth[m]->SetPoint( ieta, ieta, wdt );
                if ( m >= plot_start ) { biasMin = std::min( biasMin, b ); biasMax = std::max( biasMax, b ); widthMax = std::max( widthMax, wdt ); }
            }
            gBias[m]->SetLineColor( colors[m] ); gBias[m]->SetLineWidth( 2 );
            gWidth[m]->SetLineColor( colors[m] ); gWidth[m]->SetLineWidth( 2 );
        }

        cBias->cd( ic + 1 );
        fixpad();
        const double bpad = 0.15 * ( biasMax - biasMin );
        gBias[plot_start]->SetTitle( Form( "%s (full-tower stats);i_{#eta};mean(E_{raw}-pred) / <E_{raw}>  (fractional bias)", calo_names[ic] ) );
        gBias[plot_start]->GetYaxis()->SetRangeUser( biasMin - bpad, biasMax + bpad );
        gBias[plot_start]->Draw( "AL" );
        fixaxisG( gBias[plot_start] );
        for ( int m = plot_start + 1; m < kNMethod; ++m ) gBias[m]->Draw( "L SAME" );
        auto * lz = new TLine( -0.5, 0.0, NETA - 0.5, 0.0 ); lz->SetLineStyle( 2 ); lz->Draw( "SAME" );
        if ( ic == 0 )
        {
            auto * leg = new TLegend( 0.13, 0.68, 0.6, 0.89 );
            for ( int m = plot_start; m < kNMethod; ++m ) leg->AddEntry( gBias[m], method_names[m], "l" );
            leg->Draw();
        }

        cWidth->cd( ic + 1 );
        fixpad();
        gWidth[plot_start]->SetTitle( Form( "%s (full-tower stats);i_{#eta};RMS(E_{raw}-pred) / <E_{raw}>  (fractional width)", calo_names[ic] ) );
        gWidth[plot_start]->GetYaxis()->SetRangeUser( 0.0, 1.15 * widthMax );
        gWidth[plot_start]->Draw( "AL" );
        fixaxisG( gWidth[plot_start] );
        for ( int m = plot_start + 1; m < kNMethod; ++m ) gWidth[m]->Draw( "L SAME" );
        if ( ic == 0 )
        {
            auto * leg = new TLegend( 0.13, 0.68, 0.6, 0.89 );
            for ( int m = plot_start; m < kNMethod; ++m ) leg->AddEntry( gWidth[m], method_names[m], "l" );
            leg->Draw();
        }

        std::cout << "\n[" << calo_names[ic] << "] eta-averaged |fractional bias| and fractional width by method:" << std::endl;
        for ( int m = 0; m < kNMethod; ++m )
        {
            double sb = 0, sw = 0;
            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                double x, b, w2; gBias[m]->GetPoint( ieta, x, b ); gWidth[m]->GetPoint( ieta, x, w2 );
                sb += std::fabs( b ); sw += w2;
            }
            std::cout << "  " << method_names[m] << ":  <|bias|>=" << sb / NETA << "  <width>=" << sw / NETA << std::endl;
        }
    }
    cBias->SaveAs( Form( "%s/crosscheck_localrho_bias.png", outdir.c_str() ) );
    cBias->SaveAs( Form( "%s/crosscheck_localrho_bias.pdf", outdir.c_str() ) );
    cWidth->SaveAs( Form( "%s/crosscheck_localrho_width.png", outdir.c_str() ) );
    cWidth->SaveAs( Form( "%s/crosscheck_localrho_width.pdf", outdir.c_str() ) );

    std::cout << "\nWrote crosscheck_localrho_{bias,width}.{png,pdf} in " << outdir << std::endl;
    return 0;
}

#endif
