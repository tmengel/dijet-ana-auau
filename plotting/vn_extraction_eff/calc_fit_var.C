#include <sPhenixStyle.C>


int calc_fit_var( 
    const std::string & infile = "/sphenix/user/tmengel/JSTG-TF-03/macros/jetv2_reweight/dphi_plots/rootfiles_MC_dpsi2_8x8_cent/run31_MB_scaled111_jets_reweighted.root",
    // const std::string & infile = "/sphenix/user/tmengel/JSTG-TF-03/macros/jetv2_reweight/dphi_plots/rootfiles_MC_dpsi2_8x8_cent/run2auau_ana509_2024p022_v001_grl.root",
    const std::string & outfile = "dphi_mc.pdf"
) 
{

    auto * fin = TFile::Open( infile.c_str(), "READ" );
    if ( !fin || !fin->IsOpen() ) { std::cerr << "Error: Could not open input file: " << infile << std::endl; return -1; }
    auto * h2_dphi_cent_incl_w = (TH2F*) fin->Get("h2_dphi_cent_incl_w") -> Clone("h2_dphi_cent_incl_w");
    h2_dphi_cent_incl_w -> SetDirectory(0); // detach from file
    fin -> Close();

    std::cout << "Calculating fit variable..." << std::endl;

    const double cent_min = 0;
    const double cent_max = 10;
    const int cent_bin_i = h2_dphi_cent_incl_w -> GetYaxis() -> FindBin( cent_min );
    const int cent_bin_f = h2_dphi_cent_incl_w -> GetYaxis() -> FindBin( cent_max );
    h2_dphi_cent_incl_w -> GetYaxis() -> SetRange( cent_bin_i, cent_bin_f );

    auto * h1_dphi_cent_incl_w = (TH1F*) h2_dphi_cent_incl_w -> ProjectionX("h1_dphi_cent_incl_w", cent_bin_i, cent_bin_f );
    h1_dphi_cent_incl_w -> SetDirectory(0);
    h1_dphi_cent_incl_w -> GetXaxis() -> SetTitle( "#Delta#phi" );
    h1_dphi_cent_incl_w -> GetYaxis() -> SetTitle( "Weighted Counts" );

    SetsPhenixStyle();

    
    const int ndphi_bins = h1_dphi_cent_incl_w -> GetNbinsX();
    
    std::cout << "Number of dphi bins: " << ndphi_bins << std::endl;
    // std::cout << "Combs" << 
    // calculate the number of fit ranges
    int n_fit_ranges = 0;
    const int nvarbins = 12;
    for ( int i = 0; i < nvarbins; ++i )
    {
        for ( int j = ndphi_bins - nvarbins ; j < ndphi_bins; ++j )
        {
            n_fit_ranges++;
        }
    }
    std::cout << "Number of fit ranges: " << n_fit_ranges << std::endl;

    // TF1 * f1s[ n_fit_ranges ];
    // int fit_range_i[ n_fit_ranges ];
    // int fit_range_f[ n_fit_ranges ];
    // double v22s[ n_fit_ranges ];
    // double v23s[ n_fit_ranges ];
    std::vector<TF1*> f1s;
    std::vector<int> fit_range_i;
    std::vector<int> fit_range_f;
    std::vector<double> v22s;
    std::vector<double> v23s;
    std::vector<double> fit_range_ix;
    std::vector<double> fit_range_fx;
    for ( int i = 0; i < nvarbins; ++i )
    {
        for ( int j = ndphi_bins - nvarbins ; j < ndphi_bins; ++j )
        {
            // int range_index =  (i * (ndphi_bins - 1)) - (i * (i + 1) / 2) + (j - i - 1);
            // fit_range_i[ range_index ] = i;
            // fit_range_f[ range_index ] = j;
            // f1s[ range_index ] = new TF1( Form("f1_%d", range_index), "[0] * ( 1 + 2 * [1] * cos(2.0 * x) + 2 * [2] * cos(3.0 * x))", h1_dphi_cent_incl_w->GetBinLowEdge(i+1), h1_dphi_cent_incl_w->GetBinLowEdge(j+1) );
            // h1_dphi_cent_incl_w -> Fit( f1s[ range_index ], "0RlQ", "", h1_dphi_cent_incl_w->GetBinLowEdge(i+1), h1_dphi_cent_incl_w->GetBinLowEdge(j+1) );
            // v22s[ range_index ] = f1s[ range_index ] -> GetParameter( 1 );
            // v23s[ range_index ] = f1s[ range_index ] -> GetParameter( 2 );

            auto * f1 = new TF1( Form("f1_%d_%d", i, j), "[0] * ( 1 + 2 * [1] * cos(2.0 * x) + 2 * [2] * cos(3.0 * x))", h1_dphi_cent_incl_w->GetBinLowEdge(i+1), h1_dphi_cent_incl_w->GetBinLowEdge(j+1) );
            h1_dphi_cent_incl_w -> Fit( f1, "0RlQ", "", h1_dphi_cent_incl_w->GetBinLowEdge(i+1), h1_dphi_cent_incl_w->GetBinLowEdge(j+1) );
            double v22 = f1 -> GetParameter( 1 );
            double v23 = f1 -> GetParameter( 2 );
            f1s.push_back( f1 );
            fit_range_i.push_back( i );
            fit_range_f.push_back( j );
            v22s.push_back( v22 );
            v23s.push_back( v23 );
            fit_range_ix.push_back( h1_dphi_cent_incl_w->GetBinLowEdge(i+1) );
            fit_range_fx.push_back( h1_dphi_cent_incl_w->GetBinLowEdge(j+1) );
            // std::cout << "Fit range: " << h1_d

        }
    }


    auto * f1 = new TF1( "f1", "[0] * ( 1 + 2 * [1] * cos(2.0 * x) + 2 * [2] * cos(3.0 * x))", 0, TMath::Pi() );

    // fit full range
    h1_dphi_cent_incl_w -> Fit( f1, "0RlQ",  "", 0, TMath::Pi() );
    double v22_full = f1 -> GetParameter( 1 );
    double v23_full = f1 -> GetParameter( 2 );
    double v22_full_err = f1 -> GetParError( 1 );
    double v23_full_err = f1 -> GetParError( 2 );
    std::cout << "Full range fit: v22 = " << v22_full << " +/- " << v22_full_err << ", v23 = " << v23_full << " +/- " << v23_full_err << std::endl;

    
    TH2D * h2_v22 = new TH2D( "h2_v22", "v22 vs fit range", ndphi_bins, 0, ndphi_bins, ndphi_bins, 0, ndphi_bins );
    TH2D * h2_v23 = new TH2D( "h2_v23", "v23 vs fit range", ndphi_bins, 0, ndphi_bins, ndphi_bins, 0, ndphi_bins );
    TH1D * h1_v22_var = new TH1D( "h1_v22_var", "v22 variance vs fit range", 100, -1, 1 );
    TH1D * h1_v23_var = new TH1D( "h1_v23_var", "v23 variance vs fit range", 100, -1, 1 );
    for ( size_t idx = 0; idx < f1s.size(); ++idx )
    {
        int i = fit_range_i[ idx ];
        int j = fit_range_f[ idx ];
        double v22_var = v22s[ idx ] - v22_full;
        double v23_var = v23s[ idx ] - v23_full;
        if (v22_full !=0 )
        {
            v22_var /= v22_full;
            h2_v22 -> SetBinContent( i + 1, j + 1, abs(v22_var) );
            h1_v22_var -> Fill( v22_var );
        }
        if (v23_full !=0 )
        {
            v23_var /= v23_full;
            h2_v23 -> SetBinContent( i + 1, j + 1, abs(v23_var) );
            h1_v23_var -> Fill( v23_var );
        }
    }

    TH2D * h2_v22_var_dphix = new TH2D( "h2_v22_var_dphix", "v22 variance vs fit range (dphi)", nvarbins, fit_range_ix.front(), fit_range_ix.back(), nvarbins, fit_range_fx.front(), fit_range_fx.back() );
    for ( size_t idx = 0; idx < f1s.size(); ++idx )
    {
        double v22_var = v22s[ idx ] - v22_full;
        if (v22_full !=0 )
        {
            v22_var /= v22_full;
            int i = h2_v22_var_dphix -> GetXaxis() -> FindBin( fit_range_ix[ idx ] );
            int j = h2_v22_var_dphix -> GetYaxis() -> FindBin( fit_range_fx[ idx ] );
            h2_v22_var_dphix -> SetBinContent( i, j, abs(v22_var) );
            // std::cout << "
            // h2_v22_var_dphix -> SetBinContent( idx + 1, 1, abs(v22_var) );
        }
    }

    std::cout << "Writing output to " << outfile << std::endl;


    auto * c = new TCanvas();
    c -> cd();
    h1_dphi_cent_incl_w -> Draw("hist");
    // for ( int i = 0; i < n_fit_ranges; ++i )
    // {
    //     f1s[ i ] -> SetLineColor( kGray );
    //     // f1s[ i ] -> Draw("same");
    // }
    f1 -> SetLineColor( kRed );
    f1 -> Draw("same");

    c -> SaveAs( outfile.c_str() );

    gPad -> SetRightMargin( 0.15 );
    // gPad -> SetLogz( 1 );
    h2_v22 -> GetZaxis() -> SetTitle( "#delta v_{2,2} / v_{2,2}" );
    h2_v22 -> GetZaxis() -> SetTitleOffset( 1.5 );
    h2_v22 -> GetYaxis() -> SetTitle( "Fit Range End Bin" );
    h2_v22 -> GetXaxis() -> SetTitle( "Fit Range Start Bin" );
    // h2_v22 -> GetXaxis() -> SetRangeUser( 1, ndphi_bins-1 );
    // h2_v22 -> GetYaxis() -> SetRangeUser( 1, ndphi_bins-1 );
    // h2_v22 -> GetZaxis() -> SetRangeUser( -1
    // h2_v22 -> GetZaxis() -> SetRangeUser( -1, 1 );

    h2_v22 -> Draw("colz");
    c -> SaveAs( "v22_var.png" );

    c -> Clear();
    h2_v22_var_dphix -> GetZaxis() -> SetTitle( "#delta v_{2,2} / v_{2,2}" );
    h2_v22_var_dphix -> GetZaxis() -> SetTitleOffset( 1.5 );
    // h2_v22_var_dphix -> GetXaxis() -> SetRangeUser( 0, TMath::Pi() );
    // h2_v22_var_dphix -> GetYaxis() -> SetRangeUser( 0, TMath::Pi() );
    h2_v22_var_dphix -> GetYaxis() -> SetTitle( "Fit Range End #Delta#phi" );
    h2_v22_var_dphix -> GetXaxis() -> SetTitle( "Fit Range Start #Delta#phi" );
    // gPad -> SetRightMargin( 0.
    h2_v22_var_dphix -> Draw("colz");
    c -> SaveAs( "v22_var_dphix.png" );

    std::cout << "Writing output to calc_fit_var.root" << std::endl;

    auto fout = TFile::Open( "calc_fit_var.root", "RECREATE" );
    fout -> cd();
    h1_dphi_cent_incl_w -> Write();
    f1 -> Write();
    h1_v22_var -> Write();
    h1_v23_var -> Write();
    h2_v22 -> Write();
    h2_v23 -> Write();
    // for ( int i = 0; i < n_fit_ranges; ++i )
    // {
    //     // f1s[ i ] -> Write();
    //     if ( f1s[i] ) { f1s[ i ] -> Write(); }
    // }
    fout -> Close();
    std::cout << "Output written to calc_fit_var.root" << std::endl;


   

    return 0;

}
