// Compares per-tower calibrated energy distributions (good towers only) for
// all three calorimeters, before vs after the no-noise waveform fit. The two
// input files are produced by test_waveform.sh, which runs
// Fun4All_UEScaling_Pass1.C twice on the same events:
//   do_waveform_fit=false -> infile1 : read the existing, already-built
//       calibrated towers directly from DST_CALO_CLUSTER (no waveform sim)
//   do_waveform_fit=true  -> infile2 : rebuild towers from G4Hits through
//       CaloWaveformSim (NOISE_NONE) + CaloTowerBuilder (TEMPLATE fit)

#include <sPhenixStyle.C>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>

#include <vector>
#include <string>
#include <iostream>

int plot_tower_energy_comp(
    const std::string & infile1 = "CALO_TREE_nominal_sHijing_0_20fm-00000031-00000.root",
    const std::string & infile2 = "CALO_TREE_noNoise_sHijing_0_20fm-00000031-00000.root",
    const std::string & outdir  = "plots"
)
{
    SetsPhenixStyle();
    gStyle -> SetOptStat( 0 );
    gStyle -> SetOptFit( 0 );

    gSystem -> Exec( Form( "mkdir -p %s", outdir.c_str() ) );

    TFile * f1 = TFile::Open( infile1.c_str(), "READ" );
    TFile * f2 = TFile::Open( infile2.c_str(), "READ" );
    if ( !f1 || f1 -> IsZombie() || !f2 || f2 -> IsZombie() )
    {
        std::cout << "plot_tower_energy_comp: cannot open input file(s)" << std::endl;
        return 1;
    }
    TTree * t1 = (TTree*) f1 -> Get( "T" );
    TTree * t2 = (TTree*) f2 -> Get( "T" );

    const std::vector< std::string > tower_types = { "cemc", "hcalin", "hcalout" };

    std::vector<float> * energy1 = nullptr;
    std::vector<int>   * status1 = nullptr;
    std::vector<float> * energy2 = nullptr;
    std::vector<int>   * status2 = nullptr;

    for ( const auto & det : tower_types )
    {
        t1 -> SetBranchAddress( ( det + "_tower_energy" ).c_str(), &energy1 );
        t1 -> SetBranchAddress( ( det + "_tower_status" ).c_str(), &status1 );
        t2 -> SetBranchAddress( ( det + "_tower_energy" ).c_str(), &energy2 );
        t2 -> SetBranchAddress( ( det + "_tower_status" ).c_str(), &status2 );

        // TH1F * h1 = new TH1F( ( det + "_before" ).c_str(), ";Tower energy [GeV];Towers", 100, -1.0, 5.0 );
        // TH1F * h2 = new TH1F( ( det + "_after" ).c_str(), "", 100, -1.0, 5.0 );
        TH1F * h1 = nullptr;
        TH1F * h2 = nullptr;
        if ( det == "cemc" )
        {
            h1 = new TH1F( ( det + "_before" ).c_str(), ";Tower energy [GeV];Towers", 100, -1.0, 5.0 );
            h2 = new TH1F( ( det + "_after" ).c_str(), "", 100, -1.0, 5.0 );
        }
        else if ( det == "hcalin" )
        {
            h1 = new TH1F( ( det + "_before" ).c_str(), ";Tower energy [GeV];Towers", 100, -0.1, 1.5 );
            h2 = new TH1F( ( det + "_after" ).c_str(), "", 100, -0.1, 1.5 );
        }
        else if ( det == "hcalout" )
        {
            h1 = new TH1F( ( det + "_before" ).c_str(), ";Tower energy [GeV];Towers", 100, -0.1, 1.5 );
            h2 = new TH1F( ( det + "_after" ).c_str(), "", 100, -0.1, 1.5 );
        }

        double sum1 = 0, sum2 = 0;
        int n1 = 0, n2 = 0;

        std::vector<double> sumeT1, sumeT2; // per-event total good-tower energy

        Long64_t nentries1 = t1 -> GetEntries();
        sumeT1.reserve( nentries1 );
        for ( Long64_t ievt = 0; ievt < nentries1; ++ievt )
        {
            t1 -> GetEntry( ievt );
            double evt_sum1 = 0;
            for ( size_t ich = 0; ich < energy1 -> size(); ++ich )
            {
                if ( status1 -> at( ich ) != 1 )
                {
                    continue; // only good towers
                }
                h1 -> Fill( energy1 -> at( ich ) );
                sum1 += energy1 -> at( ich );
                evt_sum1 += energy1 -> at( ich );
                n1++;
            }
            sumeT1.push_back( evt_sum1 );
        }

        Long64_t nentries2 = t2 -> GetEntries();
        sumeT2.reserve( nentries2 );
        for ( Long64_t ievt = 0; ievt < nentries2; ++ievt )
        {
            t2 -> GetEntry( ievt );
            double evt_sum2 = 0;
            for ( size_t ich = 0; ich < energy2 -> size(); ++ich )
            {
                if ( status2 -> at( ich ) != 1 )
                {
                    continue; // only good towers
                }
                h2 -> Fill( energy2 -> at( ich ) );
                sum2 += energy2 -> at( ich );
                evt_sum2 += energy2 -> at( ich );
                n2++;
            }
            sumeT2.push_back( evt_sum2 );
        }

        std::cout << det << ": before (nominal) sum=" << sum1 << " GeV over " << n1
                   << " good towers (" << nentries1 << " events), mean="
                   << ( n1 > 0 ? sum1 / n1 : 0.0 ) << " GeV" << std::endl;
        std::cout << det << ": after (no-noise fit) sum=" << sum2 << " GeV over " << n2
                   << " good towers (" << nentries2 << " events), mean="
                   << ( n2 > 0 ? sum2 / n2 : 0.0 ) << " GeV" << std::endl;

        h1 -> SetLineColor( kGray + 2 );
        h1 -> SetFillColorAlpha( kGray, 0.5 );
        h2 -> SetLineColor( kAzure + 2 );
        h2 -> SetFillColorAlpha( kAzure + 2, 0.35 );

        TCanvas * cv = new TCanvas( ( "cv_" + det ).c_str(), "", 700, 600 );
        cv -> SetLogy();
        h1 -> SetTitle( Form( "%s good-tower energy: before vs after no-noise waveform fit", det.c_str() ) );
        h1 -> Draw( "HIST" );
        h2 -> Draw( "HIST SAME" );

        TLegend * leg = new TLegend( 0.45, 0.72, 0.88, 0.86 );
        leg -> AddEntry( h1, "before (nominal)", "f" );
        leg -> AddEntry( h2, "after (no-noise fit)", "f" );
        leg -> Draw();

        cv -> SaveAs( Form( "%s/tower_energy_comp_%s.png", outdir.c_str(), det.c_str() ) );
        cv -> SaveAs( Form( "%s/tower_energy_comp_%s.pdf", outdir.c_str(), det.c_str() ) );

        // per-event total good-tower energy (sumeT) distribution
        double lo = 0.0;
        double hi = 1.0;
        for ( double v : sumeT1 ) { if ( v > hi ) hi = v; if ( v < lo ) lo = v; }
        for ( double v : sumeT2 ) { if ( v > hi ) hi = v; if ( v < lo ) lo = v; }
        hi *= 1.1;

        TH1F * hs1 = new TH1F( ( det + "_sumeT_before" ).c_str(), ";#Sigma E_{tower} per event [GeV];Events", 100, lo, hi );
        TH1F * hs2 = new TH1F( ( det + "_sumeT_after" ).c_str(), "", 100, lo, hi );
        for ( double v : sumeT1 ) { hs1 -> Fill( v ); }
        for ( double v : sumeT2 ) { hs2 -> Fill( v ); }

        std::cout << det << ": before (nominal) mean sumeT/event="
                   << hs1 -> GetMean() << " GeV" << std::endl;
        std::cout << det << ": after (no-noise fit) mean sumeT/event="
                   << hs2 -> GetMean() << " GeV" << std::endl;

        hs1 -> SetLineColor( kGray + 2 );
        hs1 -> SetFillColorAlpha( kGray, 0.5 );
        hs2 -> SetLineColor( kAzure + 2 );
        hs2 -> SetFillColorAlpha( kAzure + 2, 0.35 );

        TCanvas * cvs = new TCanvas( ( "cv_sumeT_" + det ).c_str(), "", 700, 600 );
        cvs -> SetLogy();
        hs1 -> SetTitle( Form( "%s per-event #Sigma E_{tower}: before vs after no-noise waveform fit", det.c_str() ) );
        hs1 -> Draw( "HIST" );
        hs2 -> Draw( "HIST SAME" );

        TLegend * legs = new TLegend( 0.45, 0.72, 0.88, 0.86 );
        legs -> AddEntry( hs1, "before (nominal)", "f" );
        legs -> AddEntry( hs2, "after (no-noise fit)", "f" );
        legs -> Draw();

        cvs -> SaveAs( Form( "%s/tower_sumeT_comp_%s.png", outdir.c_str(), det.c_str() ) );
        cvs -> SaveAs( Form( "%s/tower_sumeT_comp_%s.pdf", outdir.c_str(), det.c_str() ) );
    }

    return 0;
}
