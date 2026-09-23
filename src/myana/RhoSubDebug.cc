#include "RhoSubDebug.h"

#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetInput.h>
#include <jetbase/TowerJetInput.h>

#include <jetbackground/TowerBackground.h>
#include <jetbackground/TowerRho.h>
#include <jetbackground/TowerRhov1.h>

#include <mbd/MbdOutV2.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  TowerJetInput *mk_input(Jet::SRC src, const std::string &prefix)
  {
    auto *in = new TowerJetInput(src, prefix);
    in->set_GlobalVertexType(GlobalVertex::MBD);
    return in;
  }

  double dphi_wrap(double a, double b)
  {
    double d = a - b;
    while (d > M_PI) { d -= 2 * M_PI; }
    while (d < -M_PI) { d += 2 * M_PI; }
    return d;
  }

  double dR(double eta1, double phi1, double eta2, double phi2)
  {
    const double de = eta1 - eta2;
    const double dp = dphi_wrap(phi1, phi2);
    return std::sqrt((de * de) + (dp * dp));
  }
}  // namespace

RhoSubDebug::RhoSubDebug(const std::string &name, const std::string &outfile)
  : SubsysReco(name)
  , m_outfilename(outfile)
{
  std::ostringstream nullstream;
  fastjet::ClusterSequence::set_fastjet_banner_stream(&nullstream);
  fastjet::ClusterSequence::print_banner();
  fastjet::ClusterSequence::set_fastjet_banner_stream(&std::cout);
}

RhoSubDebug::~RhoSubDebug()
{
  for (auto &v : m_coll_inputs)
  {
    for (auto *in : v) { delete in; }
    v.clear();
  }
}

int RhoSubDebug::Init(PHCompositeNode * /*topNode*/)
{
  m_file = new TFile(m_outfilename.c_str(), "RECREATE");
  book_histograms();
  return Fun4AllReturnCodes::EVENT_OK;
}

void RhoSubDebug::book_histograms()
{
  m_file->cd();

  const char *cnames[NCOLL + 1] = {"all", "emcal", "ihcal", "ohcal", "ext"};
  const char *ctitles[NCOLL + 1] = {"internal all 3 layers", "internal EMCal only",
                                    "internal IHCal only", "internal OHCal only",
                                    "external JetReco seeds"};
  for (int c = 0; c <= NCOLL; ++c)
  {
    h_kt_pt[c] = new TH1F(Form("h_kt_pt_%s", cnames[c]),
                          Form("%s;k_{T} jet p_{T} [GeV];jets", ctitles[c]), 120, 0, 60);
    h_kt_eta[c] = new TH1F(Form("h_kt_eta_%s", cnames[c]),
                           Form("%s;k_{T} jet #eta;jets", ctitles[c]), 80, -1.6, 1.6);
    h_kt_nconst[c] = new TH1F(Form("h_kt_nconst_%s", cnames[c]),
                              Form("%s;k_{T} jet N_{const};jets", ctitles[c]), 100, 0, 400);
    h_kt_njet[c] = new TH1F(Form("h_kt_njet_%s", cnames[c]),
                            Form("%s;k_{T} jets / event;events", ctitles[c]), 80, 0, 80);
  }
  for (int c = 0; c < NCOLL; ++c)
  {
    h_kt_ptovern[c] = new TH1F(Form("h_kt_ptovern_%s", cnames[c]),
                               Form("%s;k_{T} jet p_{T} / N_{const} [GeV];jets", ctitles[c]),
                               120, 0, 0.6);
  }

  h_kt_dR_vs_pt = new TH2F("h_kt_dR_vs_pt",
      "k_{T} jet axis: flipped vs unflipped;k_{T} jet p_{T} [GeV];#DeltaR(fastjet axis, physical axis)",
      60, 0, 30, 100, 0, 2.0);
  p_kt_dR_vs_pt = new TProfile("p_kt_dR_vs_pt",
      "k_{T} jet axis: flipped vs unflipped;k_{T} jet p_{T} [GeV];#LT#DeltaR#GT", 60, 0, 30);
  h_kt_deta_vs_pt = new TH2F("h_kt_deta_vs_pt",
      "k_{T} jet axis: flipped vs unflipped;k_{T} jet p_{T} [GeV];#eta_{fastjet} - #eta_{physical}",
      60, 0, 30, 100, -1.5, 1.5);
  h_kt_dR = new TH1F("h_kt_dR",
      "k_{T} jet axis: flipped vs unflipped;#DeltaR(fastjet axis, physical axis);jets", 200, 0, 2.0);
  h_kt_match_dR = new TH1F("h_kt_match_dR",
      "internal (physical axis) matched to external seed;#DeltaR;jets", 100, 0, 0.02);
  h_kt_match_dpt = new TH1F("h_kt_match_dpt",
      "internal vs matched external seed;p_{T}^{int} - p_{T}^{ext} [GeV];jets", 101, -0.05, 0.05);
  h_kt_acc_pt_all = new TH1F("h_kt_acc_pt_all",
      "all internal k_{T} jets;k_{T} jet p_{T} [GeV];jets", 60, 0, 30);
  h_kt_acc_pt_flip = new TH1F("h_kt_acc_pt_flip",
      "acceptance decision differs between axes;k_{T} jet p_{T} [GeV];jets", 60, 0, 30);

  for (int l = 0; l < NLAYER; ++l)
  {
    h_rho[l] = new TH1F(Form("h_rho_%s", layer_tag(l)),
                        Form("%s;#rho [GeV];events", layer_name(l)), 120, 0, 0.6);
    h_tow_dE[l] = new TH1F(Form("h_tow_dE_%s", layer_tag(l)),
        Form("%s;E_{MULTSUB} - E_{SUB1} [GeV];towers", layer_name(l)), 200, -1e-6, 1e-6);
    p_tow_dE_vs_ieta[l] = new TProfile(Form("p_tow_dE_vs_ieta_%s", layer_tag(l)),
        Form("%s;i#eta;#LTE_{MULTSUB} - E_{SUB1}#GT [GeV]", layer_name(l)), 24, -0.5, 23.5);
    p_ue_vs_ieta[l] = new TProfile(Form("p_ue_vs_ieta_%s", layer_tag(l)),
        Form("%s;i#eta;#LTUE#GT [GeV]", layer_name(l)), 24, -0.5, 23.5);
    p_flat_vs_ieta[l] = new TProfile(Form("p_flat_vs_ieta_%s", layer_tag(l)),
        Form("%s;i#eta;#LT#rho cosh#eta#GT [GeV]", layer_name(l)), 24, -0.5, 23.5);
    p_rmMULT_vs_ieta[l] = new TProfile(Form("p_rmMULT_vs_ieta_%s", layer_tag(l)),
        Form("%s;i#eta;#LTremoved#GT [GeV]", layer_name(l)), 24, -0.5, 23.5);
    p_rmSUB1_vs_ieta[l] = new TProfile(Form("p_rmSUB1_vs_ieta_%s", layer_tag(l)),
        Form("%s;i#eta;#LTremoved#GT [GeV]", layer_name(l)), 24, -0.5, 23.5);
    h_ue_minus_rmSUB1[l] = new TH1F(Form("h_ue_minus_rmSUB1_%s", layer_tag(l)),
        Form("%s;UE - removed_{SUB1} [GeV];i#eta strips", layer_name(l)), 200, -1e-6, 1e-6);
    h_rmMULT_minus_rmSUB1[l] = new TH1F(Form("h_rmMULT_minus_rmSUB1_%s", layer_tag(l)),
        Form("%s;removed_{MULTSUB} - removed_{SUB1} [GeV];i#eta strips", layer_name(l)),
        200, -1e-6, 1e-6);
    p_w_vs_ieta[l] = new TProfile(Form("p_w_vs_ieta_%s", layer_tag(l)),
        Form("%s;i#eta;#LTw#GT = UE / (#rho cosh#eta)", layer_name(l)), 24, -0.5, 23.5);
  }

  h_tow_dE_all = new TH1F("h_tow_dE_all",
      "all layers;E_{MULTSUB} - E_{SUB1} [GeV];towers", 200, -1e-6, 1e-6);

  h_jet_pt_rho = new TH1F("h_jet_pt_rho",
      "anti-k_{T} R=0.3 from MULTSUB towers;jet p_{T} [GeV];jets", 100, 0, 50);
  h_jet_pt_sub1 = new TH1F("h_jet_pt_sub1",
      "anti-k_{T} R=0.3 from SUB1 towers;jet p_{T} [GeV];jets", 100, 0, 50);
  h_jet_dpt = new TH1F("h_jet_dpt",
      "matched jets;p_{T}^{MULTSUB} - p_{T}^{SUB1} [GeV];jets", 201, -1e-4, 1e-4);
  // The narrow histogram above resolves the float-rounding core but sends the
  // rare anti-kT reclustering tail into its overflow, which silently biases any
  // RMS read off it. This one is wide enough to contain the tail.
  h_jet_dpt_wide = new TH1F("h_jet_dpt_wide",
      "matched jets;p_{T}^{MULTSUB} - p_{T}^{SUB1} [GeV];jets", 401, -0.5, 0.5);
  h_jet_unmatched_pt = new TH1F("h_jet_unmatched_pt",
      "jets present in one path but not the other;jet p_{T} [GeV];jets", 100, 0, 10);
  h_jet_njet_rho = new TH1F("h_jet_njet_rho", "MULTSUB;jets / event;events", 60, 0, 60);
  h_jet_njet_sub1 = new TH1F("h_jet_njet_sub1", "SUB1;jets / event;events", 60, 0, 60);
  h_jet_pt_corr = new TH2F("h_jet_pt_corr",
      "jet-by-jet;p_{T}^{MULTSUB} [GeV];p_{T}^{SUB1} [GeV]", 100, 0, 50, 100, 0, 50);

  const char *jtag[2] = {"rho", "sub1"};
  const char *jname[2] = {"Rho1 (MULTSUB towers)", "Rho2 (SUB1 towers)"};
  for (int i = 0; i < 2; ++i)
  {
    h_jet_eta[i] = new TH1F(Form("h_jet_eta_%s", jtag[i]),
        Form("%s;anti-k_{T} jet #eta;jets", jname[i]), 60, -1.2, 1.2);
    h_jet_phi[i] = new TH1F(Form("h_jet_phi_%s", jtag[i]),
        Form("%s;anti-k_{T} jet #phi;jets", jname[i]), 64, -M_PI, M_PI);
    h_jet_etaphi[i] = new TH2F(Form("h_jet_etaphi_%s", jtag[i]),
        Form("%s;anti-k_{T} jet #eta;anti-k_{T} jet #phi", jname[i]),
        48, -1.2, 1.2, 64, -M_PI, M_PI);
  }

  // truth-matched response, in variable pT bins to keep each slice populated
  const double ptbins[] = {5, 8, 11, 15, 20, 25, 30, 40, 60};
  const int nptbins = (sizeof(ptbins) / sizeof(double)) - 1;
  h_truth_pt = new TH1F("h_truth_pt", "truth jets in acceptance;p_{T}^{truth} [GeV];jets",
                        nptbins, ptbins);
  h_truth_eta = new TH1F("h_truth_eta", "truth jets;#eta^{truth};jets", 60, -1.2, 1.2);
  for (int i = 0; i < 2; ++i)
  {
    h_truth_pt_matched[i] = new TH1F(Form("h_truth_pt_matched_%s", jtag[i]),
        Form("%s;p_{T}^{truth} [GeV];matched jets", jname[i]), nptbins, ptbins);
    h_jes[i] = new TH2F(Form("h_jes_%s", jtag[i]),
        Form("%s;p_{T}^{truth} [GeV];p_{T}^{reco} / p_{T}^{truth}", jname[i]),
        nptbins, ptbins, 120, 0.0, 2.4);
  }

  for (int l = 0; l < NLAYER; ++l)
  {
    h_stored_pt[l] = new TH1F(Form("h_stored_pt_%s", layer_tag(l)),
        Form("stored %s;k_{T} jet p_{T} [GeV];jets", layer_name(l)), 120, 0, 60);
    h_stored_eta[l] = new TH1F(Form("h_stored_eta_%s", layer_tag(l)),
        Form("stored %s;k_{T} jet #eta;jets", layer_name(l)), 80, -1.6, 1.6);
    h_stored_nconst[l] = new TH1F(Form("h_stored_nconst_%s", layer_tag(l)),
        Form("stored %s;k_{T} jet N_{const};jets", layer_name(l)), 100, 0, 400);
    h_stored_ptovern[l] = new TH1F(Form("h_stored_ptovern_%s", layer_tag(l)),
        Form("stored %s;k_{T} jet p_{T}/N_{const} [GeV];jets", layer_name(l)), 120, 0, 0.6);
    h_stored_njet[l] = new TH1F(Form("h_stored_njet_%s", layer_tag(l)),
        Form("stored %s;k_{T} jets / event;events", layer_name(l)), 60, 0, 60);
    h_stored_rho[l] = new TH1F(Form("h_stored_rho_%s", layer_tag(l)),
        Form("stored %s;JetContainer::get_rho_median() [GeV];events", layer_name(l)), 120, 0, 0.6);
    h_stored_rho_minus_node[l] = new TH1F(Form("h_stored_rho_minus_node_%s", layer_tag(l)),
        Form("stored %s;container #rho - TowerRho node #rho [GeV];events", layer_name(l)),
        201, -1e-6, 1e-6);
  }

  const char *seedname[NCOLL] = {"all 3 layers", "EMCal", "IHCal", "OHCal"};
  for (int c = 0; c < NCOLL; ++c)
  {
    h_seed_etaphi_rho[c] = new TH2F(Form("h_seed_etaphi_rho_%s", coll_tag(c)),
        Form("DetermineTowerRho seeds, %s;#eta;#phi", seedname[c]), 48, -1.2, 1.2, 64, -M_PI, M_PI);
    h_seed_pt_rho[c] = new TH1F(Form("h_seed_pt_rho_%s", coll_tag(c)),
        Form("DetermineTowerRho seeds, %s;seed p_{T} [GeV];seeds", seedname[c]), 120, 0, 60);
    h_seed_eta_rho[c] = new TH1F(Form("h_seed_eta_rho_%s", coll_tag(c)),
        Form("DetermineTowerRho seeds, %s;seed #eta;seeds", seedname[c]), 60, -1.2, 1.2);
  }
  h_seed_etaphi_dtb = new TH2F("h_seed_etaphi_dtb",
      "DetermineTowerBackgroundv1 seeds;#eta;#phi", 48, -1.2, 1.2, 64, -M_PI, M_PI);
  h_seed_pt_dtb = new TH1F("h_seed_pt_dtb",
      "DetermineTowerBackgroundv1 seeds;seed p_{T} [GeV];seeds", 120, 0, 60);
  h_seed_eta_dtb = new TH1F("h_seed_eta_dtb",
      "DetermineTowerBackgroundv1 seeds;seed #eta;seeds", 60, -1.2, 1.2);
  h_seed_nseed_dtb = new TH1F("h_seed_nseed_dtb",
      "DetermineTowerBackgroundv1;seeds tagged per event;events", 10, -0.5, 9.5);
  h_seed_match_dR = new TH1F("h_seed_match_dR",
      "DTBv1 seed #rightarrow nearest all-layer DetermineTowerRho seed;#DeltaR;seeds", 100, 0, 3.2);
  h_seed_match_dR_emcal = new TH1F("h_seed_match_dR_emcal",
      "DTBv1 seed #rightarrow nearest EMCal DetermineTowerRho seed;#DeltaR;seeds", 100, 0, 3.2);
  h_seed_nshared = new TH1F("h_seed_nshared",
      "per event;background-method seeds also DetermineTowerRho seeds;events", 10, -0.5, 9.5);
  h_kt_etaphi_ext = new TH2F("h_kt_etaphi_ext",
      "external k_{T} jets (all);#eta;#phi", 64, -1.6, 1.6, 64, -M_PI, M_PI);
  for (int l = 0; l < NLAYER; ++l)
  {
    h_stored_etaphi[l] = new TH2F(Form("h_stored_etaphi_%s", layer_tag(l)),
        Form("stored %s k_{T} jets;#eta;#phi", layer_name(l)), 64, -1.6, 1.6, 64, -M_PI, M_PI);
  }

  // seeds vs the two leading truth jets
  const char *leadtag[2] = {"lead", "sublead"};
  const char *leadname[2] = {"leading", "subleading"};
  for (int i = 0; i < 2; ++i)
  {
    h_st_truth_pt[i] = new TH1F(Form("h_st_truth_pt_%s", leadtag[i]),
        Form("%s truth jet (R=0.4, |#eta|<1.1);p_{T}^{truth} [GeV];truth jets", leadname[i]), 40, 0, 80);
    h_st_truth_eta[i] = new TH1F(Form("h_st_truth_eta_%s", leadtag[i]),
        Form("%s truth jet;#eta^{truth};truth jets", leadname[i]), 22, -1.1, 1.1);
    for (int m = 0; m < NSEEDMETH; ++m)
    {
      h_st_tag_pt[m][i] = new TH1F(Form("h_st_tag_pt_%s_%s", seedmeth_tag(m), leadtag[i]),
          Form("%s truth jet found by %s;p_{T}^{truth} [GeV];truth jets", leadname[i], seedmeth_name(m).c_str()), 40, 0, 80);
      h_st_tag_eta[m][i] = new TH1F(Form("h_st_tag_eta_%s_%s", seedmeth_tag(m), leadtag[i]),
          Form("%s truth jet found by %s;#eta^{truth};truth jets", leadname[i], seedmeth_name(m).c_str()), 22, -1.1, 1.1);
    }
  }
  for (int m = 0; m < NSEEDMETH; ++m)
  {
    h_st_nfound[m] = new TH1F(Form("h_st_nfound_%s", seedmeth_tag(m)),
        Form("%s;leading truth jets found by the seeds;events", seedmeth_name(m).c_str()), 3, -0.5, 2.5);
    h_st_seed_dR[m] = new TH1F(Form("h_st_seed_dR_%s", seedmeth_tag(m)),
        Form("%s;#DeltaR(seed, nearest leading truth jet);seeds", seedmeth_name(m).c_str()), 64, 0, 3.2);
    h_st_seed_pt_all[m] = new TH1F(Form("h_st_seed_pt_all_%s", seedmeth_tag(m)),
        Form("%s;seed p_{T} [GeV];seeds", seedmeth_name(m).c_str()), 30, 0, 60);
    h_st_seed_pt_matched[m] = new TH1F(Form("h_st_seed_pt_matched_%s", seedmeth_tag(m)),
        Form("%s, seeds matched to a leading truth jet;seed p_{T} [GeV];seeds", seedmeth_name(m).c_str()), 30, 0, 60);
  }

  h_tow_nexact_diff = new TH1F("h_tow_nexact_diff",
      "per event;towers not bit-identical (MULTSUB vs SUB1);events", 50, -0.5, 49.5);
  h_jet_ncount_diff = new TH1F("h_jet_ncount_diff",
      "per event, all jets (no p_{T} cut);N_{jets}^{Rho1} - N_{jets}^{Rho2};events", 11, -5.5, 5.5);
  h_saved_kt_ncount_diff = new TH1F("h_saved_kt_ncount_diff",
      "per event, summed over layers;#Sigma |N_{stored} - N_{replicated}|;events", 11, -0.5, 10.5);

  h_saved_kt_dpt = new TH1F("h_saved_kt_dpt",
      "DetermineTowerRho saved jets vs replication;p_{T}^{saved} - p_{T}^{replica} [GeV];jets",
      201, -1e-3, 1e-3);
  h_saved_kt_dR = new TH1F("h_saved_kt_dR",
      "DetermineTowerRho saved jets vs replication;#DeltaR;jets", 100, 0, 0.02);
}

int RhoSubDebug::InitRun(PHCompositeNode * /*topNode*/)
{
  if (m_abs_jet_eta_range < 0)
  {
    m_abs_jet_eta_range = m_abs_input_eta_range - m_par;
  }

  // collection 0 = all three layers together (what the external JetReco kT pass sees)
  m_coll_inputs[0].push_back(mk_input(Jet::CEMC_TOWERINFO_RETOWER, "TOWERINFO_CALIB"));
  m_coll_inputs[0].push_back(mk_input(Jet::HCALIN_TOWERINFO, "TOWERINFO_CALIB"));
  m_coll_inputs[0].push_back(mk_input(Jet::HCALOUT_TOWERINFO, "TOWERINFO_CALIB"));
  // collections 1..3 = one per DetermineTowerRho instance in the macro
  m_coll_inputs[1].push_back(mk_input(Jet::CEMC_TOWERINFO_RETOWER, "TOWERINFO_CALIB"));
  m_coll_inputs[2].push_back(mk_input(Jet::HCALIN_TOWERINFO, "TOWERINFO_CALIB"));
  m_coll_inputs[3].push_back(mk_input(Jet::HCALOUT_TOWERINFO, "TOWERINFO_CALIB"));

  std::cout << "RhoSubDebug::InitRun - jet eta range = " << m_abs_jet_eta_range
            << ", tower eta range = " << m_abs_input_eta_range
            << ", kT R = " << m_par << ", omit n hardest = " << m_omit_nhardest
            << ", writing " << m_outfilename << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

float RhoSubDebug::grab_zvtx(PHCompositeNode *topNode)
{
  float z = 0.0;
  auto *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap || vertexmap->empty()) { return 0.0; }
  auto vertices = vertexmap->get_gvtxs_with_type({m_vertex_type});
  if (!vertices.empty() && vertices.at(0)) { z = vertices.at(0)->get_z(); }
  if (std::isnan(z) || std::abs(z) > 1e3) { z = 0.0; }
  return z;
}

float RhoSubDebug::grab_mbdQ(PHCompositeNode *topNode)
{
  auto *mbd = findNode::getClass<MbdOutV2>(topNode, "MbdOut");
  if (!mbd) { return -1.0; }
  return mbd->get_q(0) + mbd->get_q(1);
}

// Mirrors DetermineTowerRho::process_event's clustering stage. Applies only the
// "omit N hardest" half of the selector so that the eta acceptance can be studied.
std::vector<RhoSubDebug::KtJet>
RhoSubDebug::cluster_like_rho(PHCompositeNode *topNode, const std::vector<JetInput *> &inputs)
{
  std::vector<Jet *> particles{};
  for (auto *input : inputs)
  {
    std::vector<Jet *> const parts = input->get_input(topNode);
    for (auto *part : parts)
    {
      particles.push_back(part);
      particles.back()->set_id(particles.size() - 1);
    }
  }

  std::vector<fastjet::PseudoJet> pjs{};
  for (unsigned int ip = 0; ip < particles.size(); ++ip)
  {
    float e = particles[ip]->get_e();
    if (e == 0.) { continue; }
    float px = particles[ip]->get_px();
    float py = particles[ip]->get_py();
    float pz = particles[ip]->get_pz();
    if (e < 0)
    {
      // make energy = +1 MeV for purposes of clustering
      const float r = 0.001 / e;
      e *= r; px *= r; py *= r; pz *= r;
    }
    fastjet::PseudoJet pj(px, py, pz, e);
    if (std::abs(pj.eta()) > m_abs_input_eta_range) { continue; }
    pj.set_user_index(ip);
    pjs.push_back(pj);
  }

  fastjet::JetDefinition jet_def(fastjet::kt_algorithm, m_par, fastjet::E_scheme, fastjet::Best);
  fastjet::ClusterSequence cs(pjs, jet_def);
  auto fjs = cs.inclusive_jets();

  // Reproduce DetermineTowerRho::get_jet_selector() exactly. The real selector is
  //   (!SelectorNHardest(N)) * SelectorAbsEtaMax(max)
  // and fastjet's operator* applies its RIGHT operand first, so the N hardest are
  // dropped from the jets that already passed the eta cut -- not from all jets.
  const auto in_acc = fastjet::SelectorAbsEtaMax(m_abs_jet_eta_range)(fjs);
  const auto used_fj = (!fastjet::SelectorNHardest(m_omit_nhardest))(in_acc);
  auto same = [](const fastjet::PseudoJet &a, const fastjet::PseudoJet &b)
  { return a.cluster_hist_index() == b.cluster_hist_index(); };
  auto contains = [&same](const std::vector<fastjet::PseudoJet> &v, const fastjet::PseudoJet &j)
  { return std::any_of(v.begin(), v.end(), [&](const fastjet::PseudoJet &x) { return same(x, j); }); };

  std::vector<KtJet> out{};
  out.reserve(fjs.size());
  for (auto &fj : fjs)
  {
    // rebuild the jet from the ORIGINAL, unflipped constituents -- this is what
    // FastJetAlgoSub stores into the external seed container
    double px = 0, py = 0, pz = 0, e = 0;
    int n = 0;
    for (auto &c : fj.constituents())
    {
      auto *p = particles[c.user_index()];
      px += p->get_px();
      py += p->get_py();
      pz += p->get_pz();
      e += p->get_e();
      n++;
    }
    KtJet k{};
    k.pt_uf = std::sqrt((px * px) + (py * py));
    k.eta_uf = (k.pt_uf > 0) ? std::asinh(pz / k.pt_uf) : 0.0;
    k.phi_uf = std::atan2(py, px);
    k.pt_fj = fj.perp();
    k.eta_fj = fj.eta();
    k.phi_fj = fj.phi_std();
    k.e_fj = e;
    k.nconst = n;
    // SelectorAbsEtaMax keeps |eta| <= max
    k.pass_eta_fj = (std::abs(k.eta_fj) <= m_abs_jet_eta_range);
    k.pass_eta_uf = (std::abs(k.eta_uf) <= m_abs_jet_eta_range);
    k.used = contains(used_fj, fj);
    k.seed = contains(in_acc, fj) && !k.used;
    out.push_back(k);
  }

  std::sort(out.begin(), out.end(),
            [](const KtJet &a, const KtJet &b) { return a.pt_uf > b.pt_uf; });

  for (auto *p : particles) { delete p; }
  return out;
}

void RhoSubDebug::compare_kt_seeds(PHCompositeNode *topNode)
{
  const bool verbose = (m_evt <= m_nprint);

  // ---- external seeds ----
  auto *ext = findNode::getClass<JetContainer>(topNode, m_ext_kt_node);
  std::vector<KtJet> ext_jets{};
  if (ext)
  {
    for (auto *j : *ext)
    {
      if (!j) { continue; }
      KtJet k{};
      k.pt_uf = j->get_pt();
      k.eta_uf = j->get_eta();
      k.phi_uf = j->get_phi();
      k.pt_fj = k.pt_uf; k.eta_fj = k.eta_uf; k.phi_fj = k.phi_uf;
      k.e_fj = j->get_e();
      k.nconst = static_cast<int>(j->size_comp());
      k.pass_eta_fj = k.pass_eta_uf = true;
      k.used = true;
      k.seed = false;
      ext_jets.push_back(k);
    }
    std::sort(ext_jets.begin(), ext_jets.end(),
              [](const KtJet &a, const KtJet &b) { return a.pt_uf > b.pt_uf; });
    h_kt_njet[NCOLL]->Fill(ext_jets.size());
    for (const auto &k : ext_jets)
    {
      h_kt_pt[NCOLL]->Fill(k.pt_uf);
      h_kt_eta[NCOLL]->Fill(k.eta_uf);
      h_kt_nconst[NCOLL]->Fill(k.nconst);
    }
  }

  // ---- internal collections ----
  std::vector<KtJet> coll[NCOLL];
  for (int c = 0; c < NCOLL; ++c)
  {
    coll[c] = cluster_like_rho(topNode, m_coll_inputs[c]);

    int npass = 0;
    for (const auto &k : coll[c])
    {
      // exactly the jets DetermineTowerRho feeds into the median
      if (!k.used) { continue; }
      npass++;
      h_kt_pt[c]->Fill(k.pt_uf);
      h_kt_eta[c]->Fill(k.eta_uf);
      h_kt_nconst[c]->Fill(k.nconst);
      if (k.nconst > 0) { h_kt_ptovern[c]->Fill(k.pt_uf / k.nconst); }
    }
    h_kt_njet[c]->Fill(npass);
  }

  // ---- axis convention: flipped (what the selector sees) vs physical ----
  for (const auto &k : coll[0])
  {
    const double d = dR(k.eta_fj, k.phi_fj, k.eta_uf, k.phi_uf);
    h_kt_dR->Fill(d);
    h_kt_dR_vs_pt->Fill(k.pt_uf, d);
    p_kt_dR_vs_pt->Fill(k.pt_uf, d);
    h_kt_deta_vs_pt->Fill(k.pt_uf, k.eta_fj - k.eta_uf);
    h_kt_acc_pt_all->Fill(k.pt_uf);
    m_n_kt_jets++;
    if (k.pass_eta_fj != k.pass_eta_uf)
    {
      h_kt_acc_pt_flip->Fill(k.pt_uf);
      m_n_axis_flip_accept++;
    }
  }

  validate_saved_kt(topNode, coll);
  compare_seeds(topNode, coll);

  // ---- match internal (physical axis) to external seeds ----
  int nmatch = 0;
  for (const auto &a : coll[0])
  {
    if (!a.used) { continue; }
    double best = 1e9;
    const KtJet *bm = nullptr;
    for (const auto &b : ext_jets)
    {
      const double d = dR(a.eta_uf, a.phi_uf, b.eta_uf, b.phi_uf);
      if (d < best) { best = d; bm = &b; }
    }
    if (bm)
    {
      h_kt_match_dR->Fill(std::min(best, 0.0199));
      if (best < 0.01)
      {
        nmatch++;
        h_kt_match_dpt->Fill(a.pt_uf - bm->pt_uf);
      }
    }
  }

  if (verbose)
  {
    std::cout << "\n================ [1] kT SEED COMPARISON (event " << m_evt << ") ================\n";
    std::cout << "  external " << m_ext_kt_node << " : " << ext_jets.size() << " jets\n";
    for (int c = 0; c < NCOLL; ++c)
    {
      int npass = 0;
      for (const auto &k : coll[c]) { if (k.used) { npass++; } }
      std::cout << "  " << std::setw(24) << coll_name(c) << " : " << npass
                << " jets pass |eta_fastjet| < " << m_abs_jet_eta_range
                << "  (of " << coll[c].size() << " after omitting the " << m_omit_nhardest
                << " hardest)\n";
    }
    int nflip = 0;
    for (const auto &k : coll[0]) { if (k.pass_eta_fj != k.pass_eta_uf) { nflip++; } }
    std::cout << "  ==> " << nmatch << " internal all-layer jets matched an external seed"
              << " within dR<0.01 (physical axis)\n";
    std::cout << "  ==> " << nflip << " jets where the flipped and physical axes disagree"
              << " on the |eta|<" << m_abs_jet_eta_range << " acceptance\n";
  }
}

// Cross-check the JetContainers DetermineTowerRho writes via add_method(..., jet_node)
// against this module's independent replication of the same clustering.
void RhoSubDebug::validate_saved_kt(PHCompositeNode *topNode, const std::vector<KtJet> coll[NCOLL])
{
  const bool verbose = (m_evt <= m_nprint);
  long count_diff_sum = 0;
  bool any_layer = false;
  for (int c = 1; c < NCOLL; ++c)
  {
    if (m_saved_kt_node[c].empty()) { continue; }
    auto *jc = findNode::getClass<JetContainer>(topNode, m_saved_kt_node[c]);
    if (!jc)
    {
      if (verbose)
      {
        std::cout << "  !! saved kT node " << m_saved_kt_node[c] << " NOT FOUND\n";
      }
      continue;
    }

    // the replication, restricted to the jets that actually entered rho
    std::vector<const KtJet *> ref{};
    for (const auto &k : coll[c])
    {
      if (k.used && k.nconst > 0) { ref.push_back(&k); }
    }

    // rho as recorded by the container itself, cross-checked against the TowerRho node
    const int lay = c - 1;  // collections 1..3 map onto layers 0..2
    const std::string rho_nodes[NLAYER] = {m_cemc_rho_node, m_ihcal_rho_node, m_ohcal_rho_node};
    auto *rnode = findNode::getClass<TowerRhov1>(topNode, rho_nodes[lay]);
    const double node_rho = rnode ? rnode->get_rho() : 0.0;
    h_stored_rho[lay]->Fill(jc->get_rho_median());
    const double rho_diff = jc->get_rho_median() - node_rho;
    h_stored_rho_minus_node[lay]->Fill(rho_diff);
    m_max_stored_rho_diff = std::max(m_max_stored_rho_diff, std::abs(rho_diff));

    // the container also holds the omitted seeds, flagged with prop_SeedItr = 1;
    // only the jets that entered rho are compared with the replication
    const bool has_seed_flag = jc->has_property(Jet::PROPERTY::prop_SeedItr);
    const Jet::PROPERTY seed_idx = has_seed_flag ? jc->property_index(Jet::PROPERTY::prop_SeedItr)
                                                 : Jet::PROPERTY::no_property;

    long nsaved = 0;
    int nmatch = 0;
    double worst = 0.0;
    std::vector<char> ref_used(ref.size(), 0);
    for (auto *j : *jc)
    {
      if (!j) { continue; }
      if (has_seed_flag && j->get_property(seed_idx) > 0.5F) { continue; }
      nsaved++;
      m_n_saved_kt++;

      // distributions of the stored object itself
      const double jpt = j->get_pt();
      const double jeta = j->get_eta();
      const int jn = static_cast<int>(j->size_comp());
      h_stored_pt[lay]->Fill(jpt);
      h_stored_eta[lay]->Fill(jeta);
      h_stored_nconst[lay]->Fill(jn);
      if (jn > 0) { h_stored_ptovern[lay]->Fill(jpt / jn); }
      // the stored axis is the physical one, so jets beyond the nominal jet
      // acceptance are exactly the axis-convention leakage
      if (std::abs(jeta) > m_abs_jet_eta_range) { m_n_stored_outside_eta++; }
      h_stored_etaphi[lay]->Fill(jeta, j->get_phi());

      // exclusive match, so the leftovers on either side are real disagreements
      double best = 1e9;
      long bk = -1;
      for (size_t r = 0; r < ref.size(); ++r)
      {
        if (ref_used[r]) { continue; }
        const double d = dR(j->get_eta(), j->get_phi(), ref[r]->eta_uf, ref[r]->phi_uf);
        if (d < best) { best = d; bk = static_cast<long>(r); }
      }
      if (bk >= 0)
      {
        h_saved_kt_dR->Fill(std::min(best, 0.0199));
        if (best < 0.01)
        {
          ref_used[bk] = 1;
          nmatch++;
          m_n_saved_kt_matched++;
          const double dpt = j->get_pt() - ref[bk]->pt_uf;
          h_saved_kt_dpt->Fill(dpt);
          worst = std::max(worst, std::abs(dpt));
          m_max_saved_kt_dpt = std::max(m_max_saved_kt_dpt, std::abs(dpt));
        }
      }
    }
    for (const char u : ref_used) { if (!u) { m_n_replica_unmatched++; } }
    count_diff_sum += std::abs(nsaved - static_cast<long>(ref.size()));  // abs: layers must not cancel

    h_stored_njet[lay]->Fill(nsaved);

    if (verbose)
    {
      std::cout << "  saved " << std::setw(24) << m_saved_kt_node[c] << " : " << nsaved
                << " jets (replication expects " << ref.size() << "), " << nmatch
                << " matched within dR<0.01, worst |dpt| = " << worst << " GeV"
                << ", rho_median stored = " << jc->get_rho_median() << "\n";
    }
    any_layer = true;
  }
  if (any_layer)
  {
    h_saved_kt_ncount_diff->Fill(count_diff_sum);
    if (count_diff_sum != 0) { m_n_evt_saved_countdiff++; }
  }
}

// Compare the jets each path EXCLUDES from its background estimate.
//  - DetermineTowerRho drops the N hardest kT jets inside |eta_jet| <= max, separately
//    for each rho instance (one per layer here), ranked by the flipped-axis pT.
//  - DetermineTowerBackgroundv1 tags the N external kT jets (all layers, no eta cut)
//    with the highest recomputed pT as prop_SeedItr == 1 and masks their towers.
// These are different definitions; this shows where they land and how often they agree.
void RhoSubDebug::compare_seeds(PHCompositeNode *topNode, const std::vector<KtJet> coll[NCOLL])
{
  // seed positions per method, for the truth comparison. DetermineTowerRho seeds use
  // the physical axis so they sit where the jet actually is.
  std::vector<SeedPos> seeds[NSEEDMETH];
  for (int c = 0; c < NCOLL; ++c)
  {
    for (const auto &k : coll[c])
    {
      if (!k.seed) { continue; }
      h_seed_etaphi_rho[c]->Fill(k.eta_uf, k.phi_uf);
      h_seed_pt_rho[c]->Fill(k.pt_uf);
      h_seed_eta_rho[c]->Fill(k.eta_uf);
      seeds[c + 1].push_back({k.pt_uf, k.eta_uf, k.phi_uf});
    }
  }

  if (auto *ext = findNode::getClass<JetContainer>(topNode, m_ext_kt_node))
  {
    for (auto *j : *ext)
    {
      if (j) { h_kt_etaphi_ext->Fill(j->get_eta(), j->get_phi()); }
    }
  }

  // optional first iteration of an iterative background method (seed method 5)
  if (!m_bkgd_seed_node_it1.empty())
  {
    if (auto *it1 = findNode::getClass<JetContainer>(topNode, m_bkgd_seed_node_it1))
    {
      const auto idx1 = it1->property_index(Jet::PROPERTY::prop_SeedItr);
      for (auto *j : *it1)
      {
        if (j && std::abs(j->get_property(idx1) - m_bkgd_seed_it1) < 0.5)
        {
          seeds[5].push_back({j->get_pt(), j->get_eta(), j->get_phi()});
        }
      }
    }
  }

  const std::string &seednode = m_bkgd_seed_node.empty() ? m_ext_kt_node : m_bkgd_seed_node;
  auto *bkgd = findNode::getClass<JetContainer>(topNode, seednode);
  if (!bkgd)
  {
    match_seeds_to_truth(topNode, seeds);
    return;
  }
  const auto idx_itr = bkgd->property_index(Jet::PROPERTY::prop_SeedItr);

  int nseed = 0, nshared = 0;
  for (auto *j : *bkgd)
  {
    if (!j) { continue; }
    if (std::abs(j->get_property(idx_itr) - m_bkgd_seed_itr) >= 0.5) { continue; }

    nseed++;
    m_n_seed_dtb++;
    seeds[0].push_back({j->get_pt(), j->get_eta(), j->get_phi()});
    h_seed_etaphi_dtb->Fill(j->get_eta(), j->get_phi());
    h_seed_pt_dtb->Fill(j->get_pt());
    h_seed_eta_dtb->Fill(j->get_eta());

    // nearest DetermineTowerRho seed: the all-layer replica is built from the same
    // towers as the external container, so it is the like-for-like comparison
    double best = 1e9, best_em = 1e9;
    for (const auto &k : coll[0])
    {
      if (k.seed) { best = std::min(best, dR(j->get_eta(), j->get_phi(), k.eta_uf, k.phi_uf)); }
    }
    for (const auto &k : coll[1])
    {
      if (k.seed) { best_em = std::min(best_em, dR(j->get_eta(), j->get_phi(), k.eta_uf, k.phi_uf)); }
    }
    if (best < 1e8) { h_seed_match_dR->Fill(std::min(best, 3.19)); }
    if (best_em < 1e8) { h_seed_match_dR_emcal->Fill(std::min(best_em, 3.19)); }
    if (best < m_seed_share_dR) { nshared++; m_n_seed_dtb_shared++; }
  }
  h_seed_nseed_dtb->Fill(nseed);
  h_seed_nshared->Fill(std::min(nshared, 9));
  if (nseed > 0 && nshared == nseed) { m_n_evt_seeds_identical++; }

  match_seeds_to_truth(topNode, seeds);
}

// Do the jets each method excludes from its background estimate correspond to the
// hard scatter? Compares every method's seeds with the two leading truth jets
// (highest pT among truth jets inside the calorimeter acceptance).
//   efficiency: fraction of leading truth jets with a seed within m_seed_truth_dR
//   purity:     fraction of seeds within m_seed_truth_dR of a leading truth jet
void RhoSubDebug::match_seeds_to_truth(PHCompositeNode *topNode, const std::vector<SeedPos> seeds[NSEEDMETH])
{
  if (m_truth_seed_node.empty()) { return; }
  auto *truth = findNode::getClass<JetContainer>(topNode, m_truth_seed_node);
  if (!truth) { return; }

  SeedPos lead[2] = {{-1, 0, 0}, {-1, 0, 0}};
  for (auto *tj : *truth)
  {
    if (!tj || std::abs(tj->get_eta()) >= m_truth_seed_absetamax) { continue; }
    const SeedPos t{tj->get_pt(), tj->get_eta(), tj->get_phi()};
    if (t.pt > lead[0].pt)
    {
      lead[1] = lead[0];
      lead[0] = t;
    }
    else if (t.pt > lead[1].pt)
    {
      lead[1] = t;
    }
  }
  const int ntruth = (lead[0].pt > 0) + (lead[1].pt > 0);
  if (ntruth == 0) { return; }
  m_st_events++;

  for (int i = 0; i < ntruth; ++i)
  {
    h_st_truth_pt[i]->Fill(lead[i].pt);
    h_st_truth_eta[i]->Fill(lead[i].eta);
  }

  for (int m = 0; m < NSEEDMETH; ++m)
  {
    int nfound = 0;
    for (int i = 0; i < ntruth; ++i)
    {
      bool found = false;
      for (const auto &sd : seeds[m])
      {
        if (dR(sd.eta, sd.phi, lead[i].eta, lead[i].phi) < m_seed_truth_dR) { found = true; break; }
      }
      if (found)
      {
        ++nfound;
        h_st_tag_pt[m][i]->Fill(lead[i].pt);
        h_st_tag_eta[m][i]->Fill(lead[i].eta);
      }
    }
    h_st_nfound[m]->Fill(nfound);
    m_st_nfound[m][nfound]++;

    for (const auto &sd : seeds[m])
    {
      double best = 1e9;
      for (int i = 0; i < ntruth; ++i) { best = std::min(best, dR(sd.eta, sd.phi, lead[i].eta, lead[i].phi)); }
      h_st_seed_dR[m]->Fill(std::min(best, 3.19));
      h_st_seed_pt_all[m]->Fill(sd.pt);
      m_st_nseed[m]++;
      if (best < m_seed_truth_dR)
      {
        h_st_seed_pt_matched[m]->Fill(sd.pt);
        m_st_nseed_matched[m]++;
      }
    }
  }
}

// Truth-matched jet energy scale / resolution for both subtraction paths.
void RhoSubDebug::truth_match(PHCompositeNode *topNode)
{
  if (m_truth_jet_node.empty()) { return; }
  auto *truth = findNode::getClass<JetContainer>(topNode, m_truth_jet_node);
  if (!truth) { return; }

  JetContainer *reco[2] = {findNode::getClass<JetContainer>(topNode, m_jet_node_rho),
                           findNode::getClass<JetContainer>(topNode, m_jet_node_sub1)};

  for (auto *tj : *truth)
  {
    if (!tj) { continue; }
    const double tpt = tj->get_pt();
    const double teta = tj->get_eta();
    const double tphi = tj->get_phi();
    if (tpt < m_truth_ptmin) { continue; }
    h_truth_eta->Fill(teta);
    if (std::abs(teta) > m_truth_absetamax) { continue; }
    h_truth_pt->Fill(tpt);

    for (int i = 0; i < 2; ++i)
    {
      if (!reco[i]) { continue; }
      // Highest-pT reco jet inside R (the anti-kT R here is 0.3), above a soft
      // threshold. Matching to the *closest* jet instead lets a stray soft jet win
      // and smears the response down towards zero.
      double bestpt = -1;
      for (auto *rj : *reco[i])
      {
        if (!rj) { continue; }
        if (rj->get_pt() < m_reco_match_ptmin) { continue; }
        if (dR(teta, tphi, rj->get_eta(), rj->get_phi()) > m_match_dR) { continue; }
        if (rj->get_pt() > bestpt) { bestpt = rj->get_pt(); }
      }
      if (bestpt > 0)
      {
        h_truth_pt_matched[i]->Fill(tpt);
        h_jes[i]->Fill(tpt, bestpt / tpt);
      }
    }
  }
}

void RhoSubDebug::compare_towers(PHCompositeNode *topNode)
{
  const bool verbose = (m_evt <= m_nprint);
  if (verbose)
  {
    std::cout << "\n================ [2] TOWER COMPARISON MULTSUB vs SUB1 (event " << m_evt
              << ") ================\n";
  }

  const std::vector<std::string> bases = {"TOWERINFO_CALIB_CEMC_RETOWER",
                                          "TOWERINFO_CALIB_HCALIN",
                                          "TOWERINFO_CALIB_HCALOUT"};
  long evt_exact_diff = 0;
  for (int lay = 0; lay < NLAYER; ++lay)
  {
    const std::string &base = bases.at(lay);
    auto *raw = findNode::getClass<TowerInfoContainer>(topNode, base);
    auto *mul = findNode::getClass<TowerInfoContainer>(topNode, "MULTSUB_" + base);
    auto *sb1 = findNode::getClass<TowerInfoContainer>(topNode, base + "_SUB1");
    if (!raw || !mul || !sb1)
    {
      if (verbose) { std::cout << "  !! missing node for layer " << layer_name(lay) << "\n"; }
      continue;
    }

    double maxdiff = 0.0, sumdiff = 0.0;
    int ndiff = 0;
    const int n = static_cast<int>(raw->size());
    for (int ich = 0; ich < n; ++ich)
    {
      auto *tm = mul->get_tower_at_channel(ich);
      auto *ts = sb1->get_tower_at_channel(ich);

      // Bitwise check on exactly what TowerJetInput reads: the float energy and
      // the good-tower flag. Any mismatch here can change the jet collection.
      const float fm = tm->get_energy();
      const float fs = ts->get_energy();
      const bool same_e = (fm == fs) || (std::isnan(fm) && std::isnan(fs));
      const bool same_good = (tm->get_isGood() == ts->get_isGood());
      if (!same_e || !same_good)
      {
        evt_exact_diff++;
        m_n_tower_exact_differ++;
        if (m_n_printed_mismatch < 20)
        {
          m_n_printed_mismatch++;
          const unsigned int kk = raw->encode_key(ich);
          std::cout << std::setprecision(9) << "RhoSubDebug TOWER-MISMATCH evt=" << m_evt
                    << " layer=" << layer_name(lay) << " ieta=" << raw->getTowerEtaBin(kk)
                    << " iphi=" << raw->getTowerPhiBin(kk)
                    << " raw=" << raw->get_tower_at_channel(ich)->get_energy()
                    << " multsub=" << fm << " sub1=" << fs
                    << " good(multsub/sub1)=" << tm->get_isGood() << "/" << ts->get_isGood()
                    << " zvtx=" << m_vtxz << " mbdQ=" << grab_mbdQ(topNode) << std::endl;
        }
      }

      const double em = fm;
      const double es = fs;
      const double d = em - es;
      const unsigned int k = raw->encode_key(ich);
      h_tow_dE[lay]->Fill(d);
      h_tow_dE_all->Fill(d);
      p_tow_dE_vs_ieta[lay]->Fill(raw->getTowerEtaBin(k), d);
      sumdiff += std::abs(d);
      if (std::abs(d) > 1e-4) { ndiff++; m_n_tower_disagree++; }
      maxdiff = std::max(maxdiff, std::abs(d));
      m_n_tower_pairs++;
      m_sum_tower_absdiff += std::abs(d);
      m_max_tower_absdiff = std::max(m_max_tower_absdiff, std::abs(d));
    }

    if (verbose)
    {
      std::cout << "  " << std::setw(15) << layer_name(lay) << " : " << n << " towers, "
                << ndiff << " differ by >1e-4 GeV, <|dE|>=" << std::scientific
                << std::setprecision(4) << (n > 0 ? sumdiff / n : 0.0)
                << ", max|dE|=" << maxdiff << std::fixed << "\n";
    }
  }
  h_tow_nexact_diff->Fill(static_cast<double>(std::min(evt_exact_diff, 49L)));
  if (evt_exact_diff > 0) { m_n_evt_tower_exact_differ++; }
}

void RhoSubDebug::compare_ue(PHCompositeNode *topNode)
{
  const bool verbose = (m_evt <= m_nprint);

  auto *bkgd = findNode::getClass<TowerBackground>(topNode, m_bkgd_node);
  if (!bkgd) { return; }

  double rho[NLAYER] = {0, 0, 0};
  const std::string rho_nodes[NLAYER] = {m_cemc_rho_node, m_ihcal_rho_node, m_ohcal_rho_node};
  for (int l = 0; l < NLAYER; ++l)
  {
    auto *r = findNode::getClass<TowerRhov1>(topNode, rho_nodes[l]);
    if (r) { rho[l] = r->get_rho(); }
    h_rho[l]->Fill(rho[l]);
  }

  if (verbose)
  {
    std::cout << "\n================ [3] UE / TowerBackground vs ACTUAL SUBTRACTION (event "
              << m_evt << ") ================\n";
    std::cout << "  zvtx = " << m_vtxz << " cm;  rho_emcal=" << rho[0]
              << "  rho_ihcal=" << rho[1] << "  rho_ohcal=" << rho[2] << "\n";
  }

  const std::vector<std::string> bases = {"TOWERINFO_CALIB_CEMC_RETOWER",
                                          "TOWERINFO_CALIB_HCALIN",
                                          "TOWERINFO_CALIB_HCALOUT"};
  auto *geom_ih = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  auto *geom_oh = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  auto *geom_em = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if (!geom_ih || !geom_oh || !geom_em) { return; }

  const double r_em = geom_em->get_tower_geometry(
      RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, 0, 0))->get_center_radius();
  const double r_ih = geom_ih->get_tower_geometry(
      RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, 0, 0))->get_center_radius();
  const double r_oh = geom_oh->get_tower_geometry(
      RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, 0, 0))->get_center_radius();

  for (int lay = 0; lay < NLAYER; ++lay)
  {
    const std::string &base = bases.at(lay);
    auto *raw = findNode::getClass<TowerInfoContainer>(topNode, base);
    auto *mul = findNode::getClass<TowerInfoContainer>(topNode, "MULTSUB_" + base);
    auto *sb1 = findNode::getClass<TowerInfoContainer>(topNode, base + "_SUB1");
    if (!raw || !mul || !sb1) { continue; }

    RawTowerGeomContainer *geom = (lay == 2) ? geom_oh : geom_ih;
    RawTowerDefs::CalorimeterId cid =
        (lay == 2) ? RawTowerDefs::CalorimeterId::HCALOUT : RawTowerDefs::CalorimeterId::HCALIN;
    const double rad = (lay == 0) ? r_em : ((lay == 1) ? r_ih : r_oh);

    const auto &ue = bkgd->get_UE(lay);
    const int neta = geom->get_etabins();

    std::vector<double> rm_mult(neta, 0.0), rm_sub1(neta, 0.0);
    std::vector<int> ngood(neta, 0);
    const int n = static_cast<int>(raw->size());
    for (int ich = 0; ich < n; ++ich)
    {
      auto *t = raw->get_tower_at_channel(ich);
      if (!t->get_isGood()) { continue; }
      const unsigned int k = raw->encode_key(ich);
      const int ieta = raw->getTowerEtaBin(k);
      if (ieta < 0 || ieta >= neta) { continue; }
      rm_mult.at(ieta) += t->get_energy() - mul->get_tower_at_channel(ich)->get_energy();
      rm_sub1.at(ieta) += t->get_energy() - sb1->get_tower_at_channel(ich)->get_energy();
      ngood.at(ieta)++;
    }

    for (int ieta = 0; ieta < neta; ++ieta)
    {
      if (ngood.at(ieta) == 0) { continue; }
      const double eta0 = geom->get_tower_geometry(
          RawTowerDefs::encode_towerid(cid, ieta, 0))->get_eta();
      const double eta_corr = std::asinh(((std::sinh(eta0) * rad) - m_vtxz) / rad);
      const double flat = rho[lay] * std::cosh(eta_corr);
      const double rmM = rm_mult.at(ieta) / ngood.at(ieta);
      const double rmS = rm_sub1.at(ieta) / ngood.at(ieta);
      const double ueval = (ieta < static_cast<int>(ue.size())) ? ue.at(ieta) : 0.0;

      p_ue_vs_ieta[lay]->Fill(ieta, ueval);
      p_flat_vs_ieta[lay]->Fill(ieta, flat);
      p_rmMULT_vs_ieta[lay]->Fill(ieta, rmM);
      p_rmSUB1_vs_ieta[lay]->Fill(ieta, rmS);
      h_ue_minus_rmSUB1[lay]->Fill(ueval - rmS);
      h_rmMULT_minus_rmSUB1[lay]->Fill(rmM - rmS);
      if (flat != 0.0) { p_w_vs_ieta[lay]->Fill(ieta, ueval / flat); }
      m_max_ue_absdiff = std::max(m_max_ue_absdiff, std::abs(rmM - rmS));
    }
  }
}

void RhoSubDebug::compare_jets(PHCompositeNode *topNode)
{
  auto *ja = findNode::getClass<JetContainer>(topNode, m_jet_node_rho);
  auto *jb = findNode::getClass<JetContainer>(topNode, m_jet_node_sub1);
  if (!ja || !jb) { return; }

  // ---- full collections, no pT threshold, bitwise ----
  // Identical towers into identically configured JetReco instances must give the
  // same jets in the same order, so compare index by index with no matching at all.
  const long na = static_cast<long>(ja->size());
  const long nb = static_cast<long>(jb->size());
  m_n_jets_all_rho += na;
  m_n_jets_all_sub1 += nb;
  h_jet_ncount_diff->Fill(static_cast<double>(std::max(-5L, std::min(5L, na - nb))));
  bool exact = (na == nb);
  if (na != nb) { m_n_evt_jet_allcount_differ++; }
  for (long i = 0; exact && i < na; ++i)
  {
    auto *x = ja->get_jet(i);
    auto *y = jb->get_jet(i);
    if (!x || !y || x->get_px() != y->get_px() || x->get_py() != y->get_py() ||
        x->get_pz() != y->get_pz() || x->get_e() != y->get_e())
    {
      exact = false;
    }
  }
  if (!exact)
  {
    m_n_evt_jet_exact_differ++;
    if (m_n_printed_mismatch < 40)
    {
      m_n_printed_mismatch++;
      std::cout << "RhoSubDebug JET-MISMATCH evt=" << m_evt << " N(Rho1)=" << na
                << " N(Rho2)=" << nb << " zvtx=" << m_vtxz << std::endl;
    }
  }

  struct J { double pt, eta, phi; };
  std::vector<J> A{}, B{};
  for (auto *j : *ja) { if (j && j->get_pt() > 1.0) { A.push_back({j->get_pt(), j->get_eta(), j->get_phi()}); } }
  for (auto *j : *jb) { if (j && j->get_pt() > 1.0) { B.push_back({j->get_pt(), j->get_eta(), j->get_phi()}); } }
  std::sort(A.begin(), A.end(), [](const J &x, const J &y) { return x.pt > y.pt; });
  std::sort(B.begin(), B.end(), [](const J &x, const J &y) { return x.pt > y.pt; });

  h_jet_njet_rho->Fill(A.size());
  h_jet_njet_sub1->Fill(B.size());
  for (const auto &a : A)
  {
    h_jet_pt_rho->Fill(a.pt);
    h_jet_eta[0]->Fill(a.eta);
    h_jet_phi[0]->Fill(a.phi);
    h_jet_etaphi[0]->Fill(a.eta, a.phi);
  }
  for (const auto &b : B)
  {
    h_jet_pt_sub1->Fill(b.pt);
    h_jet_eta[1]->Fill(b.eta);
    h_jet_phi[1]->Fill(b.phi);
    h_jet_etaphi[1]->Fill(b.eta, b.phi);
  }

  // exclusive geometric matching, so a jet cannot be claimed twice and the
  // leftovers on each side are genuine collection differences
  std::vector<char> used(B.size(), 0);
  m_n_jets_rho += static_cast<long>(A.size());
  m_n_jets_sub1 += static_cast<long>(B.size());
  if (A.size() != B.size()) { m_n_evt_jetcount_differ++; }

  for (const auto &a : A)
  {
    double best = 1e9;
    long bk = -1;
    for (size_t k = 0; k < B.size(); ++k)
    {
      if (used[k]) { continue; }
      const double d = dR(a.eta, a.phi, B[k].eta, B[k].phi);
      if (d < best) { best = d; bk = static_cast<long>(k); }
    }
    if (bk >= 0 && best < 0.05)
    {
      used[bk] = 1;
      const double dpt = a.pt - B[bk].pt;
      h_jet_dpt->Fill(dpt);
      h_jet_dpt_wide->Fill(dpt);
      h_jet_pt_corr->Fill(a.pt, B[bk].pt);
      m_n_jets_matched++;
      m_max_jet_dpt = std::max(m_max_jet_dpt, std::abs(dpt));
    }
    else
    {
      h_jet_unmatched_pt->Fill(a.pt);
      m_n_jets_unmatched++;
    }
  }
  for (size_t k = 0; k < B.size(); ++k)
  {
    if (!used[k])
    {
      h_jet_unmatched_pt->Fill(B[k].pt);
      m_n_jets_unmatched++;
    }
  }

  if (m_evt <= m_nprint)
  {
    std::cout << "\n================ [4] JET COMPARISON (event " << m_evt << ") ================\n";
    std::cout << "  " << m_jet_node_rho << " (MULTSUB towers): " << A.size() << " jets pt>1\n";
    std::cout << "  " << m_jet_node_sub1 << " (SUB1 towers)   : " << B.size() << " jets pt>1\n";
    if (!A.empty() && !B.empty())
    {
      std::cout << "  leading: pt_MULTSUB=" << A[0].pt << "  pt_SUB1=" << B[0].pt
                << "  dpt=" << (A[0].pt - B[0].pt) << " GeV\n";
    }
  }
}

int RhoSubDebug::process_event(PHCompositeNode *topNode)
{
  m_evt++;
  m_vtxz = grab_zvtx(topNode);
  compare_kt_seeds(topNode);
  compare_towers(topNode);
  compare_ue(topNode);
  compare_jets(topNode);
  truth_match(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int RhoSubDebug::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "\n################ RhoSubDebug SUMMARY over " << m_evt << " events ################\n";
  std::cout << "  tower pairs compared      : " << m_n_tower_pairs << "\n";
  std::cout << "  pairs differing >1e-4 GeV : " << m_n_tower_disagree << "  ("
            << (m_n_tower_pairs > 0 ? 100.0 * m_n_tower_disagree / m_n_tower_pairs : 0.0) << " %)\n";
  std::cout << "  mean |E_MULTSUB - E_SUB1| : "
            << (m_n_tower_pairs > 0 ? m_sum_tower_absdiff / m_n_tower_pairs : 0.0) << " GeV\n";
  std::cout << "  max  |E_MULTSUB - E_SUB1| : " << m_max_tower_absdiff << " GeV\n";
  std::cout << "  max per-ieta |UE_MULTSUB - UE_SUB1| : " << m_max_ue_absdiff << " GeV\n";
  std::cout << "  kT jets examined          : " << m_n_kt_jets << "\n";
  std::cout << "  kT jets whose acceptance depends on the axis convention : "
            << m_n_axis_flip_accept << "  ("
            << (m_n_kt_jets > 0 ? 100.0 * m_n_axis_flip_accept / m_n_kt_jets : 0.0) << " %)\n";
  std::cout << "  kT jets saved by DetermineTowerRho : " << m_n_saved_kt << ", of which "
            << m_n_saved_kt_matched << " matched the replication (max |dpt| = "
            << m_max_saved_kt_dpt << " GeV)\n";
  std::cout << "  stored jets beyond |eta| < " << m_abs_jet_eta_range << " : "
            << m_n_stored_outside_eta << "  ("
            << (m_n_saved_kt > 0 ? 100.0 * m_n_stored_outside_eta / m_n_saved_kt : 0.0)
            << " %)  <- axis-convention leakage, visible in the stored object\n";
  std::cout << "  max |container rho - TowerRho node rho| : " << m_max_stored_rho_diff
            << " GeV\n";
  std::cout << "  --- bitwise agreement (no thresholds) ---\n";
  std::cout << "  towers not bit-identical        : " << m_n_tower_exact_differ << " in "
            << m_n_evt_tower_exact_differ << " events\n";
  std::cout << "  jets (all) Rho1 = " << m_n_jets_all_rho << ",  Rho2 = " << m_n_jets_all_sub1
            << ",  difference = " << (m_n_jets_all_rho - m_n_jets_all_sub1) << "\n";
  std::cout << "  events with different jet COUNT : " << m_n_evt_jet_allcount_differ << "\n";
  std::cout << "  events with any jet not bit-identical : " << m_n_evt_jet_exact_differ << "\n";
  std::cout << "  --- stored DetermineTowerRho jets vs replication ---\n";
  std::cout << "  events with a per-layer count mismatch : " << m_n_evt_saved_countdiff
            << ",  replica jets with no stored partner : " << m_n_replica_unmatched << "\n";
  std::cout << "  --- seeds ---\n";
  std::cout << "  " << m_label_bkgd << " seeds = " << m_n_seed_dtb << ", of which also all-layer DetermineTowerRho seeds (dR<" << m_seed_share_dR << ") = "
            << m_n_seed_dtb_shared << "  ("
            << (m_n_seed_dtb > 0 ? 100.0 * m_n_seed_dtb_shared / m_n_seed_dtb : 0.0) << " %)\n";
  std::cout << "  events where every " << m_label_bkgd << " seed is a DetermineTowerRho seed : "
            << m_n_evt_seeds_identical << "\n";
  if (m_st_events > 0)
  {
    std::cout << "  --- seeds vs the two leading truth jets (" << m_truth_seed_node << ", |eta|<"
              << m_truth_seed_absetamax << ", dR<" << m_seed_truth_dR << "), " << m_st_events << " events ---\n";
    for (int m = 0; m < NSEEDMETH; ++m)
    {
      if (m_st_nseed[m] == 0) { continue; }  // method not configured in this job
      std::cout << "  " << std::setw(30) << seedmeth_name(m) << " : found 0/1/2 leading truth jets in "
                << std::fixed << std::setprecision(1)
                << 100.0 * m_st_nfound[m][0] / m_st_events << "% / "
                << 100.0 * m_st_nfound[m][1] / m_st_events << "% / "
                << 100.0 * m_st_nfound[m][2] / m_st_events << "% of events;  seed purity "
                << (m_st_nseed[m] > 0 ? 100.0 * m_st_nseed_matched[m] / m_st_nseed[m] : 0.0) << "%\n";
    }
    std::cout << std::defaultfloat;
  }
  std::cout << "  --- jet collections (pT > 1 GeV, geometric match) ---\n";
  std::cout << "  jets  Rho1 = " << m_n_jets_rho << ",  Rho2 = " << m_n_jets_sub1
            << ",  difference = " << (m_n_jets_rho - m_n_jets_sub1) << "\n";
  std::cout << "  matched = " << m_n_jets_matched << ",  unmatched (either side) = "
            << m_n_jets_unmatched << "  ("
            << (m_n_jets_rho > 0 ? 100.0 * m_n_jets_unmatched / m_n_jets_rho : 0.0) << " %)\n";
  std::cout << "  events where the jet COUNT differs : " << m_n_evt_jetcount_differ
            << "  (" << (m_evt > 0 ? 100.0 * m_n_evt_jetcount_differ / m_evt : 0.0) << " %)\n";
  std::cout << "  max |dpt| among matched jets : " << m_max_jet_dpt << " GeV\n";
  std::cout << "############################################################\n";

  if (m_file)
  {
    m_file->cd();
    TNamed("label_path1", m_label_path1.c_str()).Write();
    TNamed("label_path2", m_label_path2.c_str()).Write();
    TNamed("label_bkgd_seeds", m_label_bkgd.c_str()).Write();
    TNamed("label_bkgd_seeds_it1", m_label_bkgd_it1.c_str()).Write();
    TNamed("seed_share_dR", Form("%g", m_seed_share_dR)).Write();
    // the actual cut this run used, so plotting can draw it correctly instead of
    // assuming the historical 0.7 default
    TNamed("jet_abs_eta", Form("%g", m_abs_jet_eta_range)).Write();
    TNamed("input_abs_eta", Form("%g", m_abs_input_eta_range)).Write();
    m_file->Write();
    m_file->Close();
    m_file = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
