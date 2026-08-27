
#ifndef _CLOSURE_TEST_C_
#define _CLOSURE_TEST_C_

#include <cdbobjects/CDBTTree.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TProfile.h>
#include <TString.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libcdbobjects.so )

// Closure test for the rho eta-shape calibration produced by make_rho_calib.C.
// Re-reads the same per-event tower-statistics tree that the calibration was
// derived from (output/<run>/rho_calib_<run>.root) and, for every event and
// eta ring, computes the raw ratio w_raw = <E_tower> / (rho*cosh(eta_corr))
// alongside the calibrated ratio w_raw / w_calib(ieta,izbin,imbd) taken from
// the CDBTTree written by make_rho_calib.C. If the calibration closes, the
// calibrated <E_tower>/(rho*cosh(eta)) profile vs ieta should be flat at 1,
// in contrast to the raw (uncalibrated) profile which reproduces the eta
// shape seen in the calib QA file.
int closure_test(
    const int run_number = 54912,
    const std::string & tree_dir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/output",
    const std::string & calib_dir = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/calibs",
    const std::string & outplot = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/closure_test_54912.pdf",
    const std::string & outfile = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/closure_test_54912.root"
)
{
    const std::string treefile = Form( "%s/%d/rho_calib_%d.root", tree_dir.c_str(), run_number, run_number );
    const std::string calibfile = Form( "%s/rho_calib_%d.root", calib_dir.c_str(), run_number );

    std::cout << "Tree file:  " << treefile << std::endl;
    std::cout << "Calib file: " << calibfile << std::endl;

    // -- load calibration (same conventions as make_rho_calib.C / SubtractTowersRhov1) --
    auto * cdbttree = new CDBTTree( calibfile );
    cdbttree->LoadCalibrations();

    const int NETA = cdbttree->GetSingleIntValue( "n_eta" );
    const int n_zbins = cdbttree->GetSingleIntValue( "n_zvtx_bins" );
    const int n_mbdbins = cdbttree->GetSingleIntValue( "n_mbdQ_bins" );
    std::cout << "n_eta=" << NETA << " n_zvtx_bins=" << n_zbins << " n_mbdQ_bins=" << n_mbdbins << std::endl;

    std::vector<float> zvtx_edges( n_zbins + 1 );
    for ( int i = 0; i <= n_zbins; ++i ) zvtx_edges[i] = cdbttree->GetSingleFloatValue( Form( "zvtx_edge_%d", i ) );
    std::vector<float> mbdQ_edges( n_mbdbins + 1 );
    for ( int i = 0; i <= n_mbdbins; ++i ) mbdQ_edges[i] = cdbttree->GetSingleFloatValue( Form( "mbdQ_edge_%d", i ) );

    static const int kNCalo = 3;
    static const char * calo_names[kNCalo] = { "cemc", "hcalin", "hcalout" };
    static const char * calib_field_names[kNCalo] = { "w_cemc", "w_hcalin", "w_hcalout" };

    auto find_bin = []( const float val, const std::vector<float> & edges ) -> int
    {
        if ( edges.size() < 2 || std::isnan(val) || val < edges.front() || val >= edges.back() ) return -1;
        for ( size_t i = 0; i + 1 < edges.size(); ++i ) if ( val >= edges[i] && val < edges[i + 1] ) return static_cast<int>(i);
        return -1;
    };
    auto encode_channel = [&]( const int ieta, const int izbin, const int imbd ) -> int
    {
        return izbin * ( n_mbdbins * NETA ) + imbd * NETA + ieta;
    };

    // -- open the per-event tree --
    auto * f = TFile::Open( treefile.c_str() );
    if ( !f || f->IsZombie() ) { std::cerr << "Error: could not open " << treefile << std::endl; return -1; }
    auto * t = (TTree*) f->Get( "T" );
    const Long64_t nentries = t->GetEntries();
    std::cout << "Total entries: " << nentries << std::endl;

    t->SetBranchStatus( "*", false );
    float zvrtx = 0.0, mbdQ = 0.0;
    float rho[kNCalo];
    float mean_E[kNCalo][24];
    int ntowers[kNCalo][24];
    float corrected_eta[kNCalo][24];
    auto enable = [&]( const char * name, void * addr ) { t->SetBranchStatus( name, true ); t->SetBranchAddress( name, addr ); };
    enable( "zvrtx", &zvrtx );
    enable( "mbdQ", &mbdQ );
    enable( "rho", rho );
    enable( "mean_E", mean_E );
    enable( "ntowers", ntowers );
    enable( "corrected_eta", corrected_eta );

    std::vector<TProfile*> pRaw( kNCalo, nullptr );
    std::vector<TProfile*> pCalib( kNCalo, nullptr );
    for ( int ic = 0; ic < kNCalo; ++ic )
    {
        pRaw[ic] = new TProfile( Form( "raw_%s", calo_names[ic] ),
            Form( "%s;i_{#eta};<E_{tower}> / (#rho cosh(#eta_{corr}))", calo_names[ic] ), NETA, -0.5, NETA - 0.5 );
        pCalib[ic] = new TProfile( Form( "calibrated_%s", calo_names[ic] ),
            Form( "%s (calibrated);i_{#eta};[<E_{tower}> / (#rho cosh(#eta_{corr}))] / w_{calib}", calo_names[ic] ), NETA, -0.5, NETA - 0.5 );
    }

    Long64_t n_used = 0;
    for ( Long64_t i = 0; i < nentries; ++i )
    {
        t->GetEntry( i );
        const int izbin = find_bin( zvrtx, zvtx_edges );
        const int imbd = find_bin( mbdQ, mbdQ_edges );
        if ( izbin < 0 || imbd < 0 ) continue;

        bool used = false;
        for ( int ic = 0; ic < kNCalo; ++ic )
        {
            if ( rho[ic] <= 0.0 ) continue;
            for ( int ieta = 0; ieta < NETA; ++ieta )
            {
                if ( ntowers[ic][ieta] <= 0 ) continue;
                const float denom = rho[ic] * std::cosh( corrected_eta[ic][ieta] );
                if ( denom <= 0.0 ) continue;
                const float raw_ratio = mean_E[ic][ieta] / denom;

                const int channel = encode_channel( ieta, izbin, imbd );
                const float w = cdbttree->GetFloatValue( channel, calib_field_names[ic] );

                pRaw[ic]->Fill( ieta, raw_ratio, ntowers[ic][ieta] );
                if ( w > 0.0f ) pCalib[ic]->Fill( ieta, raw_ratio / w, ntowers[ic][ieta] );
                used = true;
            }
        }
        if ( used ) ++n_used;
    }
    std::cout << n_used << " / " << nentries << " events used." << std::endl;

    // -- plots + summary --
    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * c = new TCanvas( "c", "closure test", 1800, 600 );
    c->Divide( 3, 1 );
    for ( int ic = 0; ic < kNCalo; ++ic )
    {
        c->cd( ic + 1 );
        pRaw[ic]->SetLineColor( kRed + 1 );
        pRaw[ic]->SetMarkerColor( kRed + 1 );
        pRaw[ic]->SetMarkerStyle( 20 );
        pCalib[ic]->SetLineColor( kBlue + 1 );
        pCalib[ic]->SetMarkerColor( kBlue + 1 );
        pCalib[ic]->SetMarkerStyle( 21 );

        pRaw[ic]->SetStats( 0 );
        pCalib[ic]->SetStats( 0 );

        const double ymax = 1.5 * pRaw[ic]->GetMaximum();
        pRaw[ic]->GetYaxis()->SetRangeUser( 0, ymax );
        pRaw[ic]->Draw( "PE" );
        pCalib[ic]->Draw( "PE SAME" );

        auto * line = new TLine( -0.5, 1.0, NETA - 0.5, 1.0 );
        line->SetLineStyle( 2 );
        line->SetLineColor( kBlack );
        line->Draw( "SAME" );

        auto * leg = new TLegend( 0.15, 0.75, 0.5, 0.88 );
        leg->AddEntry( pRaw[ic], "raw (uncalibrated)", "pl" );
        leg->AddEntry( pCalib[ic], "calibrated / w", "pl" );
        leg->Draw();

        pRaw[ic]->Write();
        pCalib[ic]->Write();

        // RMS spread of the calibrated profile across eta rings -- closure metric
        double sum = 0, sum2 = 0; int n = 0;
        for ( int b = 1; b <= NETA; ++b )
        {
            if ( pCalib[ic]->GetBinEntries( b ) <= 0 ) continue;
            const double v = pCalib[ic]->GetBinContent( b );
            sum += v; sum2 += v * v; ++n;
        }
        const double mean = ( n > 0 ) ? sum / n : 0.0;
        const double rms = ( n > 0 ) ? std::sqrt( std::max( 0.0, sum2 / n - mean * mean ) ) : 0.0;
        std::cout << calo_names[ic] << ": calibrated mean=" << mean << "  eta-to-eta RMS=" << rms
                  << " (" << 100.0 * rms / mean << "%)" << std::endl;
    }
    c->Write();
    c->SaveAs( outplot.c_str() );
    fout->Close();

    std::cout << "Wrote " << outfile << " and " << outplot << std::endl;
    return 0;
}

#endif
