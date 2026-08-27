
#ifndef _COMPARE_RHO_SUB1_ETA_C_
#define _COMPARE_RHO_SUB1_ETA_C_

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TRatioPlot.h>
#include <TString.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

// Compares the jet eta distribution reconstructed with the rho-subtraction
// method (rho_jet_*, calibrated with the new run-54912 eta-shape weights
// via Fun4All_Dijets_AuAu.C) against the Sub1/iterative-subtraction method
// (sub1_jet_*) from the same AnaTreev1 output, to check the rho calibration
// hasn't distorted the jet eta shape relative to the established Sub1
// baseline.
int compare_rho_sub1_eta(
    const std::string & infile = "/sphenix/user/tmengel/dijet-ana-auau/macros/data/rho_calibrated_54912_00000.root",
    const double pt_min = 7.0,
    const std::string & outplot = "/sphenix/user/tmengel/dijet-ana-auau/macros/data/rho_vs_sub1_eta_54912_00000.pdf",
    const std::string & outfile = "/sphenix/user/tmengel/dijet-ana-auau/macros/data/rho_vs_sub1_eta_54912_00000.root"
)
{
    auto * f = TFile::Open( infile.c_str() );
    if ( !f || f->IsZombie() ) { std::cerr << "Error: could not open " << infile << std::endl; return -1; }
    auto * t = (TTree*) f->Get( "T" );
    const Long64_t nentries = t->GetEntries();
    std::cout << "Total entries: " << nentries << std::endl;

    std::vector<float> * rho_eta = nullptr; std::vector<float> * rho_pT = nullptr;
    std::vector<float> * sub1_eta = nullptr; std::vector<float> * sub1_pT = nullptr;
    t->SetBranchAddress( "rho_jet_eta", &rho_eta );
    t->SetBranchAddress( "rho_jet_pT", &rho_pT );
    t->SetBranchAddress( "sub1_jet_eta", &sub1_eta );
    t->SetBranchAddress( "sub1_jet_pT", &sub1_pT );

    const int NBINS = 11;
    auto * hRho = new TH1F( "h_rho_jet_eta", Form( "jet #eta, p_{T} > %.0f GeV;jet #eta;jets (normalized)", pt_min ), NBINS, -1.1, 1.1 );
    auto * hSub1 = new TH1F( "h_sub1_jet_eta", "sub1", NBINS, -1.1, 1.1 );

    long n_rho = 0, n_sub1 = 0;
    for ( Long64_t i = 0; i < nentries; ++i )
    {
        t->GetEntry( i );
        if ( rho_eta )
        {
            for ( size_t j = 0; j < rho_eta->size(); ++j )
            {
                if ( (*rho_pT)[j] < pt_min ) continue;
                hRho->Fill( (*rho_eta)[j] );
                ++n_rho;
            }
        }
        if ( sub1_eta )
        {
            for ( size_t j = 0; j < sub1_eta->size(); ++j )
            {
                if ( (*sub1_pT)[j] < pt_min ) continue;
                hSub1->Fill( (*sub1_eta)[j] );
                ++n_sub1;
            }
        }
    }
    std::cout << "jets with pT > " << pt_min << " GeV: rho=" << n_rho << "  sub1=" << n_sub1 << std::endl;

    if ( hRho->Integral() > 0 ) hRho->Scale( 1.0 / hRho->Integral() );
    if ( hSub1->Integral() > 0 ) hSub1->Scale( 1.0 / hSub1->Integral() );

    hRho->SetLineColor( kBlue + 1 );
    hRho->SetMarkerColor( kBlue + 1 );
    hRho->SetMarkerStyle( 20 );
    hRho->SetStats( 0 );
    hSub1->SetLineColor( kRed + 1 );
    hSub1->SetMarkerColor( kRed + 1 );
    hSub1->SetMarkerStyle( 21 );
    hSub1->SetStats( 0 );

    auto * c = new TCanvas( "c", "rho vs sub1 jet eta", 900, 800 );
    auto * rp = new TRatioPlot( hRho, hSub1 );
    rp->Draw();
    rp->GetLowerRefYaxis()->SetTitle( "rho / sub1" );
    rp->GetUpperRefYaxis()->SetTitle( "jets (normalized)" );

    c->cd();
    auto * leg = new TLegend( 0.15, 0.75, 0.5, 0.88 );
    leg->AddEntry( hRho, "rho method (new eta calib)", "pl" );
    leg->AddEntry( hSub1, "sub1 method", "pl" );
    rp->GetUpperPad()->cd();
    leg->Draw();

    c->SaveAs( outplot.c_str() );

    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    hRho->Write();
    hSub1->Write();

    // simple closure metric: max fractional deviation of the ratio from 1
    // in bins with enough entries in both histograms
    double max_dev = 0.0; int max_dev_bin = -1;
    for ( int b = 1; b <= NBINS; ++b )
    {
        const double r = hRho->GetBinContent( b );
        const double s = hSub1->GetBinContent( b );
        if ( r <= 0 || s <= 0 ) continue;
        const double dev = std::fabs( r / s - 1.0 );
        if ( dev > max_dev ) { max_dev = dev; max_dev_bin = b; }
    }
    std::cout << "Max |rho/sub1 - 1| = " << max_dev
              << " at eta bin center " << hRho->GetBinCenter( max_dev_bin )
              << " (" << 100.0 * max_dev << "%)" << std::endl;

    fout->Close();
    std::cout << "Wrote " << outfile << " and " << outplot << std::endl;
    return 0;
}

#endif
