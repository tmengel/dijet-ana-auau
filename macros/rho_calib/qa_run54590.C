
#ifndef _QA_RUN54590_C_
#define _QA_RUN54590_C_

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TString.h>
#include <TTree.h>

#include <iostream>
#include <string>

// Sanity-check plots for the reference-calibration run pick (54590): cent,
// zvrtx, mbdQ, and per-layer sumeT distributions, overlaid against a second
// run (54912, the run used throughout the earlier closure test/Fun4All
// checks) as a normal-run baseline, to confirm 54590 isn't hiding anything
// pathological (truncated zvtx range, double-peaked mbdQ, bad centrality
// coding, absurd sumeT, etc.) despite being the "most typical" calibration
// by eta-weight shape.
int qa_run54590(
    const int run_ref = 54590,
    const int run_compare = 54912,
    const std::string & tree_dir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/output",
    const std::string & outplot = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/qa_run54590.pdf",
    const std::string & outpng = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/qa_run54590.png"
)
{
    auto * fRef = TFile::Open( Form( "%s/%d/rho_calib_%d.root", tree_dir.c_str(), run_ref, run_ref ) );
    auto * fCmp = TFile::Open( Form( "%s/%d/rho_calib_%d.root", tree_dir.c_str(), run_compare, run_compare ) );
    auto * tRef = (TTree*) fRef->Get( "T" );
    auto * tCmp = (TTree*) fCmp->Get( "T" );
    std::cout << "run " << run_ref << ": " << tRef->GetEntries() << " entries" << std::endl;
    std::cout << "run " << run_compare << ": " << tCmp->GetEntries() << " entries" << std::endl;

    auto * c = new TCanvas( "c", "run QA", 1600, 1200 );
    c->Divide( 3, 3 );

    auto draw_pair = [&]( int pad, const char * var, const char * cut, const char * title,
                           int nbins, double lo, double hi, bool logy ) {
        c->cd( pad );
        if ( logy ) gPad->SetLogy();
        auto * hRef = new TH1F( Form( "h_%s_ref", var ), title, nbins, lo, hi );
        auto * hCmp = new TH1F( Form( "h_%s_cmp", var ), title, nbins, lo, hi );
        tRef->Draw( Form( "%s>>h_%s_ref", var, var ), cut, "goff" );
        tCmp->Draw( Form( "%s>>h_%s_cmp", var, var ), cut, "goff" );
        if ( hRef->Integral() > 0 ) hRef->Scale( 1.0 / hRef->Integral() );
        if ( hCmp->Integral() > 0 ) hCmp->Scale( 1.0 / hCmp->Integral() );
        hRef->SetLineColor( kRed + 1 );
        hRef->SetLineWidth( 2 );
        hRef->SetStats( 0 );
        hCmp->SetLineColor( kAzure + 2 );
        hCmp->SetLineWidth( 2 );
        hCmp->SetLineStyle( 2 );
        hCmp->SetStats( 0 );
        const double ymax = 1.3 * std::max( hRef->GetMaximum(), hCmp->GetMaximum() );
        hRef->GetYaxis()->SetRangeUser( logy ? 1e-5 : 0.0, ymax );
        hRef->Draw( "HIST" );
        hCmp->Draw( "HIST SAME" );
        if ( pad == 1 )
        {
            auto * leg = new TLegend( 0.5, 0.7, 0.88, 0.88 );
            leg->AddEntry( hRef, Form( "run %d (reference)", run_ref ), "l" );
            leg->AddEntry( hCmp, Form( "run %d (compare)", run_compare ), "l" );
            leg->Draw();
        }
        std::cout << var << ": ref mean=" << hRef->GetMean() << " rms=" << hRef->GetRMS()
                  << " | cmp mean=" << hCmp->GetMean() << " rms=" << hCmp->GetRMS() << std::endl;
    };

    draw_pair( 1, "cent", "cent>=0", "centrality;cent (%);events (norm.)", 50, 0, 100, false );
    draw_pair( 2, "zvrtx", "", "z-vertex;z_{vtx} (cm);events (norm.)", 60, -60, 60, false );
    draw_pair( 3, "mbdQ", "mbdQ>0", "MBD charge sum;mbdQ;events (norm.)", 100, 0, 2000, true );
    draw_pair( 4, "sumeT[0]", "sumeT[0]>0", "sumeT CEMC;sumeT (GeV);events (norm.)", 100, 0, 200, true );
    draw_pair( 5, "sumeT[1]", "sumeT[1]>0", "sumeT HCALIN;sumeT (GeV);events (norm.)", 100, 0, 100, true );
    draw_pair( 6, "sumeT[2]", "sumeT[2]>0", "sumeT HCALOUT;sumeT (GeV);events (norm.)", 100, 0, 150, true );
    draw_pair( 7, "rho[0]", "rho[0]>0", "#rho CEMC;#rho;events (norm.)", 100, 0, 5, true );
    draw_pair( 8, "rho[1]", "rho[1]>0", "#rho HCALIN;#rho;events (norm.)", 100, 0, 3, true );
    draw_pair( 9, "rho[2]", "rho[2]>0", "#rho HCALOUT;#rho;events (norm.)", 100, 0, 3, true );

    c->SaveAs( outplot.c_str() );
    c->SaveAs( outpng.c_str() );
    std::cout << "Wrote " << outplot << " and " << outpng << std::endl;
    return 0;
}

#endif
