#include <iostream>

#include <sPhenixStyle.C>

#include <myana/AnaUtils.h>

R__LOAD_LIBRARY( libmyana.so )

int makeplot(
    const char* inputFile = "probs.root",
    const char* outputName = "probabilities_AA_r03")
{
    constexpr int nCents = 4;
    SetsPhenixStyle();

    const char* centLabels[nCents] = {
        "0-10%",
        "10-30%",
        "30-50%",
        "50-90%"
    };

    const int colors[nCents] = {
        kRed + 1,
        kBlue + 1,
        kGreen + 2,
        kOrange + 7
    };

    auto* fin = TFile::Open(inputFile, "READ");

    if (!fin || fin->IsZombie())
    {
        std::cerr << "Could not open " << inputFile << std::endl;
        return 1;
    }

    TH1D* hPt[nCents] = {};
    TH1D* hEff[nCents] = {};

    for (int icent = 0; icent < nCents; ++icent)
    {
        hPt[icent] = dynamic_cast<TH1D*>(
            fin->Get(Form("h1_pt_reco_cent%d", icent))
        );

        hEff[icent] = dynamic_cast<TH1D*>(
            fin->Get(Form("h1_prob_cent%d", icent))
        );

        if (!hPt[icent] || !hEff[icent])
        {
            std::cerr
                << "Missing histogram(s) for centrality index "
                << icent
                << "\nExpected:"
                << "\n  h1_pt_reco_cent" << icent
                << "\n  h1_prob_cent" << icent
                << std::endl;

            fin->Close();
            return 1;
        }

        // Detach histograms so they remain valid after closing the file.
        hPt[icent]->SetDirectory(nullptr);
        hEff[icent]->SetDirectory(nullptr);

        hPt[icent]->SetLineColor(colors[icent]);
        hPt[icent]->SetMarkerColor(colors[icent]);
        hPt[icent]->SetLineWidth(2);

        hEff[icent]->SetLineColor(colors[icent]);
        hEff[icent]->SetMarkerColor(colors[icent]);
        hEff[icent]->SetLineWidth(2);
    }

    fin->Close();


    auto * c = new TCanvas("c", "c", 1600, 1200);
    c -> cd();
    gPad -> SetLogy();
    gPad -> SetLeftMargin(0.18);
    gPad -> SetRightMargin(0.03);
    gPad -> SetBottomMargin(0.15);

    hPt[0] -> Draw("HIST");
    hPt[0] -> GetXaxis() -> SetTitle("p_{T} [GeV]");
    hPt[0] -> GetYaxis() -> SetTitle("#frac{1}{N_{evt}} #frac{dN_{jet}}{dp_{T}} [GeV^{-1}]");
    hPt[0] -> GetXaxis() -> SetRangeUser(5.0, 30.0);
    hPt[0] -> GetYaxis() -> SetRangeUser(1.e-7, 5.0);
    hPt[0] -> GetXaxis() -> SetTitleSize(0.050);
    hPt[0] -> GetYaxis() -> SetTitleSize(0.050);
    hPt[0] -> GetXaxis() -> SetLabelSize(0.040);
    hPt[0] -> GetYaxis() -> SetLabelSize(0.040);
    hPt[0] -> GetYaxis() -> SetTitleOffset(1.45);
    hPt[0] -> GetXaxis() -> SetTitleOffset(1.15);
    hPt[0] -> SetLineWidth(2);
    // hPt[0] -> SetMinimum(1.e-7);
    // hPt[0] -> SetMaximum(5.0);

    for (int icent = 1; icent < nCents; ++icent)
    {
        hPt[icent] -> Draw("HIST SAME");
    }

    auto * legend = new TLegend(0.75, 0.72, 0.92, 0.87);
    legend -> SetTextSize(0.037);
    legend -> SetBorderSize(0);
    legend -> SetFillStyle(0);
    legend -> SetFillColor(0);
    for (int icent = 0; icent < nCents; ++icent)
    {
        legend -> AddEntry(
            hPt[icent],
            centLabels[icent],
            "l"
        );
    }
    legend -> Draw("SAME");

    std::vector<std::string> tags = {
        "#it{#bf{sPHENIX}} Internal",
        "Au+Au#kern[0.1]{#sqrt{s_{_{NN}}}} = 200 GeV",
        "R = 0.3 , anti-k_{t}"
    };
    double xx = 0.40;
    double yy = 0.87;
    double dy = 0.05;
    for (const auto& tag : tags)
    {
        AnaUtils::myText(xx, yy, kBlack, tag.c_str(), 0.04);
        yy -= dy;
    }

    c -> SaveAs(Form("%s.pdf", outputName));

    c -> cd();
    gPad -> SetLogy(0);

    hEff[0]->SetTitle("");
    hEff[0]->GetXaxis()->SetTitle("p_{T} [GeV]");
    hEff[0]->GetYaxis()->SetTitle("Subleading efficiency");

    hEff[0] -> GetXaxis() -> SetRangeUser(5.0, 30.0);
    hEff[0] -> GetYaxis() -> SetRangeUser(0.0, 1.15);
    hEff[0] -> GetXaxis() -> SetTitleSize(0.050);
    hEff[0] -> GetYaxis() -> SetTitleSize(0.050);
    hEff[0] -> GetXaxis() -> SetLabelSize(0.040);
    hEff[0] -> GetYaxis() -> SetLabelSize(0.040);
    hEff[0] -> GetYaxis() -> SetTitleOffset(1.45);
    hEff[0] -> GetXaxis() -> SetTitleOffset(1.15);
    hEff[0] -> SetLineWidth(2);
   

    hEff[0]->Draw("HIST");

    for (int icent = 1; icent < nCents; ++icent)
    {
        hEff[icent]->Draw("HIST SAME");
    }

    auto* unityLine = new TLine(5.0, 1.0, 30.0, 1.0);
    unityLine->SetLineColor(kBlack);
    unityLine->SetLineStyle(2);
    unityLine->SetLineWidth(2);
    unityLine->Draw("SAME");

    auto legend2 = new TLegend(0.75, 0.42, 0.92, 0.57);
    legend2 -> SetTextSize(0.037);
    legend2 -> SetBorderSize(0);
    legend2 -> SetFillStyle(0);
    legend2 -> SetFillColor(0);
    for (int icent = 0; icent < nCents; ++icent)
    {
        legend2 -> AddEntry(
            hEff[icent],
            centLabels[icent],
            "l"
        );
    }
    legend2 -> Draw("SAME");

    yy = 0.67;
    xx = 0.60;
    for (const auto& tag : tags)
    {
        AnaUtils::myText(xx, yy, kBlack, tag.c_str(), 0.04);
        yy -= dy;
    }
    c -> SaveAs(Form("%s_eff.pdf", outputName));

    // std::cout
    //     << "Created:"
    //     << "\n  " << outputName << ".pdf"
    //     << "\n  " << outputName << ".png"
    //     << "\n  " << outputName << ".root"
    //     << std::endl;

    return 0;
}
