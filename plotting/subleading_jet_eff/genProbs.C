#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"

namespace
{
  constexpr int kNCents = 4;

  constexpr double kCentBins[kNCents + 1] =
  {
    0.0, 10.0, 30.0, 50.0, 90.0
  };

  // Detector and analysis geometry.
  //
  // sPHENIX calorimeter coverage is approximately |eta| < 1.1.
  // Therefore Delta eta_calo = 2.2.
  constexpr double kEtaCaloMax = 1.1;
  constexpr double kDeltaEtaCalo = 2.0 * kEtaCaloMax;

  constexpr double kJetR = 0.3;

  // For an away-side requirement such as
  //
  //   Delta phi_12 > 3*pi/4,
  //
  // the allowed subleading region is centered at pi and has half-width pi/4.
  //
  // Its total azimuthal width is therefore 2*DeltaPhi.
  constexpr double kDeltaPhi = TMath::Pi() / 8.0;

  constexpr int kNPtBins = 45;
  constexpr double kPtMin = 5.0;
  constexpr double kPtMax = 50.0;

  int GetCentralityBin(const double centrality)
  {
    for (int icent = 0; icent < kNCents; ++icent)
    {
      if (centrality >= kCentBins[icent] &&
          centrality < kCentBins[icent + 1])
      {
        return icent;
      }
    }

    return -1;
  }
}

int genProbs()
{
  const double jetArea = TMath::Pi() * kJetR * kJetR;
  const double caloArea = 2.0 * TMath::Pi() * kDeltaEtaCalo;

  const double subleadingSearchArea =
    2.0 * kDeltaPhi * kDeltaEtaCalo;

  const double A1 =
    (caloArea - jetArea - subleadingSearchArea) / caloArea;

  const double A2 =
    (caloArea - jetArea) / caloArea;

  if (A1 < 0.0 || A1 > 1.0 || A2 < 0.0 || A2 > 1.0)
  {
    std::cerr
      << "Invalid geometric acceptance factors:"
      << "\n  A1 = " << A1
      << "\n  A2 = " << A2
      << std::endl;

    return 1;
  }

  std::cout
    << "Geometric acceptance:"
    << "\n  A_calo = " << caloArea
    << "\n  pi R^2 = " << jetArea
    << "\n  away-side search area = " << subleadingSearchArea
    << "\n  A1 = " << A1
    << "\n  A2 = " << A2
    << std::endl;

  auto* fin = TFile::Open(
    "/sphenix/user/tmengel/JetUESub-JSTG-TF03/macros/data/"
    "rootfiles/run2auau_rho_jet.root",
    "READ"
  );

  if (!fin || fin->IsZombie())
  {
    std::cerr << "Could not open input file." << std::endl;
    return 1;
  }

  auto* t = dynamic_cast<TTree*>(fin->Get("T"));

  if (!t)
  {
    std::cerr << "Could not find TTree named T." << std::endl;
    fin->Close();
    return 1;
  }

  const Long64_t nentries = t->GetEntries();

  std::cout
    << "Total entries in TTree: "
    << nentries
    << std::endl;

  t->SetBranchStatus("*", false);

  std::vector<float>* jet_pT = nullptr;
  std::vector<float>* jet_eta = nullptr;
  std::vector<float>* jet_phi = nullptr;

  float sumeT = 0.0;
  int cent = -1;

  float rho_cemc = 0.0;
  float rho_hcalin = 0.0;
  float rho_hcalout = 0.0;

  t->SetBranchStatus("jet_pT", true);
  t->SetBranchStatus("jet_eta", true);
  t->SetBranchStatus("jet_phi", true);
  t->SetBranchStatus("sumeT_all", true);
  t->SetBranchStatus("cent", true);
  t->SetBranchStatus("rho_cemc", true);
  t->SetBranchStatus("rho_hcalin", true);
  t->SetBranchStatus("rho_hcalout", true);

  t->SetBranchAddress("jet_pT", &jet_pT);
  t->SetBranchAddress("jet_eta", &jet_eta);
  t->SetBranchAddress("jet_phi", &jet_phi);
  t->SetBranchAddress("sumeT_all", &sumeT);
  t->SetBranchAddress("cent", &cent);
  t->SetBranchAddress("rho_cemc", &rho_cemc);
  t->SetBranchAddress("rho_hcalin", &rho_hcalin);
  t->SetBranchAddress("rho_hcalout", &rho_hcalout);

  TH1D* h1_pt_counts[kNCents] = {};
  TH1D* h1_pt_reco[kNCents] = {};
  TH1D* h1_prob[kNCents] = {};

  double nEvents[kNCents] = {};

  for (int icent = 0; icent < kNCents; ++icent)
  {
    h1_pt_counts[icent] = new TH1D(
      Form("h1_pt_counts_cent%d", icent),
      Form(
        "%.0f-%.0f%%;p_{T}^{jet} [GeV];Jet counts",
        kCentBins[icent],
        kCentBins[icent + 1]
      ),
      kNPtBins,
      kPtMin,
      kPtMax
    );

    h1_pt_counts[icent]->Sumw2();

    h1_prob[icent] = new TH1D(
      Form("h1_prob_cent%d", icent),
      Form(
        "%.0f-%.0f%%;p_{T,2} [GeV];Subleading efficiency",
        kCentBins[icent],
        kCentBins[icent + 1]
      ),
      kNPtBins,
      kPtMin,
      kPtMax
    );

    h1_prob[icent]->Sumw2();
  }

  const Long64_t progressStep =
    std::max<Long64_t>(1, nentries / 10);

  // Require the full R=0.3 jet to be contained inside the calorimeter.
  const double jetEtaMax = kEtaCaloMax - kJetR;

  for (Long64_t ientry = 0; ientry < nentries; ++ientry)
  {
    t->GetEntry(ientry);

    if (ientry > 0 && ientry % progressStep == 0)
    {
      std::cout
        << "Processing entry "
        << ientry
        << " / "
        << nentries
        << std::endl;
    }

    const int icent = GetCentralityBin(cent);

    if (icent < 0)
    {
      continue;
    }

    // This assumes each TTree entry is one accepted minimum-bias event.
    //
    // Add any event-quality, trigger, or minimum-bias selection before
    // incrementing nEvents if such selections are required.
    nEvents[icent] += 1.0;

    if (!jet_pT || !jet_eta)
    {
      continue;
    }

    if (jet_pT->size() != jet_eta->size())
    {
      std::cerr
        << "jet_pT and jet_eta size mismatch in entry "
        << ientry
        << std::endl;
      continue;
    }

    for (std::size_t ijet = 0; ijet < jet_pT->size(); ++ijet)
    {
      const double pt = jet_pT->at(ijet);
      const double eta = jet_eta->at(ijet);

      if (!std::isfinite(pt) || !std::isfinite(eta))
      {
        continue;
      }

      if (std::abs(eta) >= jetEtaMax)
      {
        continue;
      }

      // Do not impose an upper-pT cut. ROOT stores pT >= 50 GeV
      // in the overflow bin, which is included in lambda.
      if (pt < kPtMin)
      {
        continue;
      }

      h1_pt_counts[icent]->Fill(pt);
    }
  }

  for (int icent = 0; icent < kNCents; ++icent)
  {
    std::cout
      << "Centrality "
      << kCentBins[icent]
      << "-"
      << kCentBins[icent + 1]
      << "%: N_evt = "
      << nEvents[icent]
      << std::endl;

    if (nEvents[icent] <= 0.0)
    {
      std::cerr
        << "No events found for centrality index "
        << icent
        << std::endl;
      continue;
    }

    const int nBins = h1_pt_counts[icent]->GetNbinsX();

    for (int ibin = 1; ibin <= nBins; ++ibin)
    {
      /*
       * For bin [x,x+1]:
       *
       * lambda =
       *   A1 * integral_x^{x+1} dN/(Nevt dpT) dpT
       * + A2 * integral_{x+1}^{infinity} dN/(Nevt dpT) dpT.
       *
       * Since the histogram contains raw counts, the integrals are:
       *
       *   N_current / Nevt
       *   N_above   / Nevt.
       *
       * Bin nBins+1 is the ROOT overflow bin and represents jets
       * above the histogram maximum.
       */

      const double nCurrent =
        h1_pt_counts[icent]->GetBinContent(ibin);

      const double nAbove =
        h1_pt_counts[icent]->Integral(
          ibin + 1,
          nBins + 1
        );

      const double lambda =
        (A1 * nCurrent + A2 * nAbove) / nEvents[icent];

      const double efficiency = std::exp(-lambda);

      /*
       * Assuming independent Poisson uncertainties on N_current
       * and N_above:
       *
       * Var(lambda) =
       *   [A1^2 N_current + A2^2 N_above] / Nevt^2.
       *
       * Since epsilon = exp(-lambda):
       *
       * sigma_epsilon = epsilon * sigma_lambda.
       */
      const double lambdaVariance =
        (
          A1 * A1 * nCurrent +
          A2 * A2 * nAbove
        ) /
        (
          nEvents[icent] *
          nEvents[icent]
        );

      const double efficiencyError =
        efficiency * std::sqrt(lambdaVariance);

      h1_prob[icent]->SetBinContent(ibin, efficiency);
      h1_prob[icent]->SetBinError(ibin, efficiencyError);
    }

    // Make the published/displayed spectrum:
    //
    //   (1/Nevt) dNjet/dpT.
    h1_pt_reco[icent] = dynamic_cast<TH1D*>(
      h1_pt_counts[icent]->Clone(
        Form("h1_pt_reco_cent%d", icent)
      )
    );

    h1_pt_reco[icent]->SetTitle(
      Form(
        "%.0f-%.0f%%;"
        "p_{T}^{jet} [GeV];"
        "#frac{1}{N_{evt}} #frac{dN_{jet}}{dp_{T}} [GeV^{-1}]",
        kCentBins[icent],
        kCentBins[icent + 1]
      )
    );

    h1_pt_reco[icent]->Scale(
      1.0 / nEvents[icent],
      "width"
    );
  }

  auto* fout = TFile::Open("probs.root", "RECREATE");

  if (!fout || fout->IsZombie())
  {
    std::cerr << "Could not create probs.root." << std::endl;
    fin->Close();
    return 1;
  }

  fout->cd();

  auto* h1_event_count = new TH1D(
    "h1_event_count",
    ";Centrality bin;N_{evt}",
    kNCents,
    0,
    kNCents
  );

  for (int icent = 0; icent < kNCents; ++icent)
  {
    h1_event_count->GetXaxis()->SetBinLabel(
      icent + 1,
      Form(
        "%.0f-%.0f%%",
        kCentBins[icent],
        kCentBins[icent + 1]
      )
    );

    h1_event_count->SetBinContent(
      icent + 1,
      nEvents[icent]
    );

    h1_pt_counts[icent]->Write();

    if (h1_pt_reco[icent])
    {
      h1_pt_reco[icent]->Write();
    }

    h1_prob[icent]->Write();
  }

  h1_event_count->Write();

  // Save geometry parameters in a small histogram for reproducibility.
  auto* h1_geometry = new TH1D(
    "h1_geometry",
    ";Parameter;Value",
    6,
    0,
    6
  );

  h1_geometry->GetXaxis()->SetBinLabel(1, "R");
  h1_geometry->GetXaxis()->SetBinLabel(2, "DeltaEtaCalo");
  h1_geometry->GetXaxis()->SetBinLabel(3, "DeltaPhi");
  h1_geometry->GetXaxis()->SetBinLabel(4, "Acalo");
  h1_geometry->GetXaxis()->SetBinLabel(5, "A1");
  h1_geometry->GetXaxis()->SetBinLabel(6, "A2");

  h1_geometry->SetBinContent(1, kJetR);
  h1_geometry->SetBinContent(2, kDeltaEtaCalo);
  h1_geometry->SetBinContent(3, kDeltaPhi);
  h1_geometry->SetBinContent(4, caloArea);
  h1_geometry->SetBinContent(5, A1);
  h1_geometry->SetBinContent(6, A2);

  h1_geometry->Write();

  fout->Close();
  fin->Close();

  std::cout << "Wrote spectra and efficiencies to probs.root."
            << std::endl;

  return 0;
}