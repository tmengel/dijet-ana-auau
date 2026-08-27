
#ifndef _FINAL_CROSSCHECK_C_
#define _FINAL_CROSSCHECK_C_

#include <cdbobjects/CDBTTree.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libcdbobjects.so )

// Final cross-check of the rho eta-shape calibration on run 54590 (the
// reference run picked by find_reference_calib.C). Compares, per event and
// per (layer,ieta): the raw tower energy mean_E[ieta] against four
// background predictions --
//   uncal   = rho * cosh(eta)                                  (no eta calibration)
//   cal     = rho * cosh(eta) * w(ieta,izbin,imbd)              (static run calibration)
//   evt     = rho * cosh(eta) * mean_E[ieta] / <mean_E>_evtavg  (literal event-by-event ask --
//             NOTE: algebraically this reduces to a rescaling of mean_E[ieta] itself, i.e. it
//             uses the ring's OWN observed value as part of its own prediction. Included as asked,
//             flagged as circular in the writeup.)
//   hybrid  = cal * (<mean_E>_evtavg / <<mean_E>>_populationRef) (proposed alternative: keep the
//             robust static SHAPE w(ieta) but let the OVERALL normalization float event-by-event,
//             using only the event's own other rings -- not circular the way "evt" is.)
// Reports both the BIAS (mean residual, i.e. closure) and the WIDTH (RMS residual, i.e.
// event-by-event fluctuation) vs ieta for all four, plus raw E +/- sigma and sliced (z-bin,
// mbdQ-bin) views of the uncalibrated & calibrated ratios.
int final_crosscheck(
    const int run = 54590,
    const std::string & tree_dir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/output",
    const std::string & calib_dir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/calibs",
    const std::string & outdir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib"
)
{
    static const int kNCalo = 3;
    static const int NETA = 24;
    static const char * calo_names[kNCalo] = { "cemc", "hcalin", "hcalout" };
    static const char * calib_field_names[kNCalo] = { "w_cemc", "w_hcalin", "w_hcalout" };

    // -- load calibration --
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

    // -- open per-event tree --
    auto * f = TFile::Open( Form( "%s/%d/rho_calib_%d.root", tree_dir.c_str(), run, run ) );
    auto * t = (TTree*) f->Get( "T" );
    const Long64_t nentries = t->GetEntries();
    std::cout << "run " << run << ": " << nentries << " entries" << std::endl;

    t->SetBranchStatus( "*", false );
    float zvrtx = 0, mbdQ = 0;
    float rho[kNCalo];
    float mean_E[kNCalo][NETA];
    int ntowers[kNCalo][NETA];
    float corrected_eta[kNCalo][NETA];
    auto enable = [&]( const char * name, void * addr ) { t->SetBranchStatus( name, true ); t->SetBranchAddress( name, addr ); };
    enable( "zvrtx", &zvrtx ); enable( "mbdQ", &mbdQ ); enable( "rho", rho );
    enable( "mean_E", mean_E ); enable( "ntowers", ntowers ); enable( "corrected_eta", corrected_eta );

    // ===================== PASS 1: population reference <<mean_E>> per (layer,izbin,imbd) =====================
    std::vector<std::vector<std::vector<double>>> refSum( kNCalo, std::vector<std::vector<double>>( n_zbins, std::vector<double>( n_mbdbins, 0.0 ) ) );
    std::vector<std::vector<std::vector<long>>> refN( kNCalo, std::vector<std::vector<long>>( n_zbins, std::vector<long>( n_mbdbins, 0 ) ) );

    for ( Long64_t i = 0; i < nentries; ++i )
    {
        t->GetEntry( i );
        const int izbin = find_bin( zvrtx, zvtx_edges );
        const int imbd = find_bin( mbdQ, mbdQ_edges );
        if ( izbin < 0 || imbd < 0 ) continue;
        for ( int ic = 0; ic < kNCalo; ++ic )
        {
            if ( rho[ic] <= 0.0 ) continue;
            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                if ( ntowers[ic][ieta] <= 0 ) continue;
                refSum[ic][izbin][imbd] += mean_E[ic][ieta];
                ++refN[ic][izbin][imbd];
            }
        }
    }
    std::vector<std::vector<std::vector<double>>> refAvg( kNCalo, std::vector<std::vector<double>>( n_zbins, std::vector<double>( n_mbdbins, 1.0 ) ) );
    for ( int ic = 0; ic < kNCalo; ++ic )
        for ( int iz = 0; iz < n_zbins; ++iz )
            for ( int imbd = 0; imbd < n_mbdbins; ++imbd )
                if ( refN[ic][iz][imbd] > 0 ) refAvg[ic][iz][imbd] = refSum[ic][iz][imbd] / refN[ic][iz][imbd];
    std::cout << "Pass 1 done (population reference <<mean_E>> per layer/izbin/imbd)." << std::endl;

    // ===================== PASS 2: accumulate raw E, 4 method predictions & residuals =====================
    static const int kNMethod = 4; // 0=uncal 1=cal 2=evt(literal) 3=hybrid
    const char * method_names[kNMethod] = { "uncal", "cal", "evt", "hybrid" };

    // marginalized-over-everything accumulators, per (layer, ieta): raw E, and residual per method
    struct Acc { double sum = 0, sumsq = 0; long n = 0;
        void add( double v ) { sum += v; sumsq += v * v; ++n; }
        double mean() const { return n > 0 ? sum / n : 0.0; }
        double rms() const { return n > 0 ? std::sqrt( std::max( 0.0, sumsq / n - mean() * mean() ) ) : 0.0; } };

    std::vector<std::vector<Acc>> accRaw( kNCalo, std::vector<Acc>( NETA ) );
    std::vector<std::vector<std::vector<Acc>>> accResid( kNMethod, std::vector<std::vector<Acc>>( kNCalo, std::vector<Acc>( NETA ) ) );

    // sliced (z-bin only, mbdQ marginalized) and (mbdQ-bin only, z marginalized) ratio accumulators, uncal & cal
    std::vector<std::vector<std::vector<Acc>>> accRatioUncal_z( kNCalo, std::vector<std::vector<Acc>>( n_zbins, std::vector<Acc>( NETA ) ) );
    std::vector<std::vector<std::vector<Acc>>> accRatioCal_z( kNCalo, std::vector<std::vector<Acc>>( n_zbins, std::vector<Acc>( NETA ) ) );
    std::vector<std::vector<std::vector<Acc>>> accRatioUncal_mbd( kNCalo, std::vector<std::vector<Acc>>( n_mbdbins, std::vector<Acc>( NETA ) ) );
    std::vector<std::vector<std::vector<Acc>>> accRatioCal_mbd( kNCalo, std::vector<std::vector<Acc>>( n_mbdbins, std::vector<Acc>( NETA ) ) );

    for ( Long64_t i = 0; i < nentries; ++i )
    {
        t->GetEntry( i );
        if ( nentries >= 10 && i % ( nentries / 10 ) == 0 && i > 0 )
            std::cout << "pass 2: " << i << " / " << nentries << std::endl;

        const int izbin = find_bin( zvrtx, zvtx_edges );
        const int imbd = find_bin( mbdQ, mbdQ_edges );
        if ( izbin < 0 || imbd < 0 ) continue;

        for ( int ic = 0; ic < kNCalo; ++ic )
        {
            if ( rho[ic] <= 0.0 ) continue;

            // this event's <mean_E>_evtavg over all valid ieta (for methods evt & hybrid)
            double evtSum = 0.0; int evtN = 0;
            for ( int ieta = 0; ieta < NETA; ++ieta )
                if ( ntowers[ic][ieta] > 0 ) { evtSum += mean_E[ic][ieta]; ++evtN; }
            if ( evtN == 0 ) continue;
            const double evtAvg = evtSum / evtN;

            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                if ( ntowers[ic][ieta] <= 0 ) continue;
                const double E_raw = mean_E[ic][ieta];
                const double eta = corrected_eta[ic][ieta];
                const double coshE = std::cosh( eta );
                const double w = w_of( ic, ieta, izbin, imbd );

                const double pred_uncal = rho[ic] * coshE;
                const double pred_cal = rho[ic] * coshE * w;
                const double pred_evt = ( evtAvg > 0.0 ) ? rho[ic] * coshE * ( E_raw / evtAvg ) : pred_cal;
                const double refA = refAvg[ic][izbin][imbd];
                const double pred_hybrid = ( refA > 0.0 ) ? pred_cal * ( evtAvg / refA ) : pred_cal;

                accRaw[ic][ieta].add( E_raw );
                accResid[0][ic][ieta].add( E_raw - pred_uncal );
                accResid[1][ic][ieta].add( E_raw - pred_cal );
                accResid[2][ic][ieta].add( E_raw - pred_evt );
                accResid[3][ic][ieta].add( E_raw - pred_hybrid );

                if ( pred_uncal > 0 ) accRatioUncal_z[ic][izbin][ieta].add( E_raw / pred_uncal );
                if ( pred_cal > 0 ) accRatioCal_z[ic][izbin][ieta].add( E_raw / pred_cal );
                if ( pred_uncal > 0 ) accRatioUncal_mbd[ic][imbd][ieta].add( E_raw / pred_uncal );
                if ( pred_cal > 0 ) accRatioCal_mbd[ic][imbd][ieta].add( E_raw / pred_cal );
            }
        }
    }
    std::cout << "Pass 2 done." << std::endl;

    // consistent fix for the y-axis title getting clipped/overlapping the tick labels
    // in narrow multi-pad canvases -- wider left margin + pushed-out title offset.
    auto fixpad = []( double left = 0.17 ) { gPad->SetLeftMargin( left ); gPad->SetRightMargin( 0.04 ); };
    auto fixaxis = []( TH1 * h ) { h->GetYaxis()->SetTitleOffset( 1.6 ); h->GetYaxis()->SetLabelSize( 0.035 ); };
    auto fixaxisG = []( TGraph * g ) { g->GetYaxis()->SetTitleOffset( 1.6 ); g->GetYaxis()->SetLabelSize( 0.035 ); };

    // ===================== PLOT 1: raw E +/- sigma vs eta =====================
    {
        auto * c = new TCanvas( "c_rawE", "raw E +/- sigma", 1800, 700 );
        c->Divide( 3, 1 );
        for ( int ic = 0; ic < kNCalo; ++ic )
        {
            c->cd( ic + 1 );
            fixpad();
            auto * g = new TGraphErrors( NETA );
            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                g->SetPoint( ieta, ieta, accRaw[ic][ieta].mean() );
                g->SetPointError( ieta, 0.0, accRaw[ic][ieta].rms() );
            }
            g->SetTitle( Form( "%s;i_{#eta};raw tower <E> #pm #sigma (GeV)", calo_names[ic] ) );
            g->SetLineColor( kAzure + 2 );
            g->SetMarkerColor( kAzure + 2 );
            g->SetMarkerStyle( 20 );
            g->SetFillColorAlpha( kAzure + 2, 0.25 );
            g->Draw( "A3" );
            fixaxisG( g );
            auto * gl = (TGraph*) g->Clone(); gl->SetLineWidth( 2 ); gl->Draw( "LP SAME" );
        }
        c->SaveAs( Form( "%s/crosscheck_raw_E_sigma.png", outdir.c_str() ) );
        c->SaveAs( Form( "%s/crosscheck_raw_E_sigma.pdf", outdir.c_str() ) );
    }

    // ===================== PLOT 2 & 3: sliced ratio views (z-bin and mbdQ-bin), uncal vs cal =====================
    auto make_sliced_plot = [&]( const std::string & tag, const std::string & label,
                                  std::vector<std::vector<std::vector<Acc>>> & accUncalSliced,
                                  std::vector<std::vector<std::vector<Acc>>> & accCalSliced,
                                  int nslices, const std::vector<float> & edges ) {
        auto * c = new TCanvas( Form( "c_%s", tag.c_str() ), label.c_str(), 1800, 1200 );
        c->Divide( 3, 2 );
        for ( int ic = 0; ic < kNCalo; ++ic )
        {
            for ( int row = 0; row < 2; ++row ) // row 0 = uncal, row 1 = cal
            {
                c->cd( row * 3 + ic + 1 );
                fixpad();
                double ymax = 0.0, ymin = 1e18;
                std::vector<TGraph*> graphs;
                for ( int isl = 0; isl < nslices; ++isl )
                {
                    auto * g = new TGraph( NETA );
                    for ( int ieta = 0; ieta < NETA; ++ieta )
                    {
                        const double v = ( row == 0 ) ? accUncalSliced[ic][isl][ieta].mean() : accCalSliced[ic][isl][ieta].mean();
                        g->SetPoint( ieta, ieta, v );
                        if ( v > 0 ) { ymax = std::max( ymax, v ); ymin = std::min( ymin, v ); }
                    }
                    const double frac = double( isl ) / std::max( 1, nslices - 1 );
                    const int color = TColor::GetColor( (float) frac, (float) ( 0.2 ), (float) ( 1.0 - frac ) );
                    g->SetLineColor( color );
                    g->SetLineWidth( 2 );
                    graphs.push_back( g );
                }
                auto * frame = graphs[0];
                frame->SetTitle( Form( "%s -- %s;i_{#eta};E_{raw}/pred (%s)", calo_names[ic], label.c_str(), row == 0 ? "uncalibrated" : "calibrated" ) );
                const double pad = 0.15 * ( ymax - ymin );
                frame->GetYaxis()->SetRangeUser( std::max( 0.0, ymin - pad ), ymax + pad );
                frame->Draw( "AL" );
                fixaxisG( frame );
                for ( size_t k = 1; k < graphs.size(); ++k ) graphs[k]->Draw( "L SAME" );
                auto * line = new TLine( -0.5, 1.0, NETA - 0.5, 1.0 );
                line->SetLineStyle( 2 ); line->Draw( "SAME" );
            }
        }
        c->SaveAs( Form( "%s/crosscheck_%s.png", outdir.c_str(), tag.c_str() ) );
        c->SaveAs( Form( "%s/crosscheck_%s.pdf", outdir.c_str(), tag.c_str() ) );
    };
    make_sliced_plot( "zslices", "z-bin slices (colored low#rightarrowhigh |z|... actually edge order)", accRatioUncal_z, accRatioCal_z, n_zbins, zvtx_edges );
    make_sliced_plot( "mbdslices", "mbdQ-bin slices (colored peripheral#rightarrowcentral)", accRatioUncal_mbd, accRatioCal_mbd, n_mbdbins, mbdQ_edges );

    // ===================== PLOT 4: bias (mean residual / <E_raw>) vs ieta, 4 methods =====================
    // PLOT 5: width (RMS residual / <E_raw>) vs ieta, 4 methods
    {
        auto * cBias = new TCanvas( "c_bias", "fractional bias", 1800, 700 );
        cBias->Divide( 3, 1 );
        auto * cWidth = new TCanvas( "c_width", "fractional width", 1800, 700 );
        cWidth->Divide( 3, 1 );

        int colors[kNMethod] = { kAzure + 2, kRed + 1, (int) kSpring + 4, (int) kOrange + 7 };
        std::string labels[kNMethod] = { "uncalibrated", "calibrated (static w)", "event-by-event (literal)", "hybrid (static shape, event scale)" };

        for ( int ic = 0; ic < kNCalo; ++ic )
        {
            cBias->cd( ic + 1 );
            cWidth->cd( ic + 1 );
        }
        // method 0 (uncalibrated/raw rho) is excluded from the plots -- its much larger
        // bias/width dominates the y-range and hides the differences between the correction
        // methods themselves. Still computed & reported in the console summary below.
        const int plot_start = 1;
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
            gBias[plot_start]->SetTitle( Form( "%s;i_{#eta};mean(E_{raw}-pred) / <E_{raw}>  (fractional bias)", calo_names[ic] ) );
            gBias[plot_start]->GetYaxis()->SetRangeUser( biasMin - bpad, biasMax + bpad );
            gBias[plot_start]->Draw( "AL" );
            fixaxisG( gBias[plot_start] );
            for ( int m = plot_start + 1; m < kNMethod; ++m ) gBias[m]->Draw( "L SAME" );
            auto * lz = new TLine( -0.5, 0.0, NETA - 0.5, 0.0 ); lz->SetLineStyle( 2 ); lz->Draw( "SAME" );
            if ( ic == 0 )
            {
                auto * leg = new TLegend( 0.13, 0.68, 0.55, 0.89 );
                for ( int m = plot_start; m < kNMethod; ++m ) leg->AddEntry( gBias[m], labels[m].c_str(), "l" );
                leg->Draw();
            }

            cWidth->cd( ic + 1 );
            fixpad();
            gWidth[plot_start]->SetTitle( Form( "%s;i_{#eta};RMS(E_{raw}-pred) / <E_{raw}>  (fractional width)", calo_names[ic] ) );
            gWidth[plot_start]->GetYaxis()->SetRangeUser( 0.0, 1.15 * widthMax );
            gWidth[plot_start]->Draw( "AL" );
            fixaxisG( gWidth[plot_start] );
            for ( int m = plot_start + 1; m < kNMethod; ++m ) gWidth[m]->Draw( "L SAME" );
            if ( ic == 0 )
            {
                auto * leg = new TLegend( 0.13, 0.68, 0.55, 0.89 );
                for ( int m = plot_start; m < kNMethod; ++m ) leg->AddEntry( gWidth[m], labels[m].c_str(), "l" );
                leg->Draw();
            }

            // console summary
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
        cBias->SaveAs( Form( "%s/crosscheck_bias.png", outdir.c_str() ) );
        cBias->SaveAs( Form( "%s/crosscheck_bias.pdf", outdir.c_str() ) );
        cWidth->SaveAs( Form( "%s/crosscheck_width.png", outdir.c_str() ) );
        cWidth->SaveAs( Form( "%s/crosscheck_width.pdf", outdir.c_str() ) );
    }

    std::cout << "\nWrote crosscheck_{raw_E_sigma,zslices,mbdslices,bias,width}.{png,pdf} in " << outdir << std::endl;
    return 0;
}

#endif
