#ifndef MYANA_RHOSUBDEBUG_H
#define MYANA_RHOSUBDEBUG_H

#include <fun4all/SubsysReco.h>

#include <jetbase/Jet.h>
#include <globalvertex/GlobalVertex.h>

#include <string>
#include <vector>

class PHCompositeNode;
class JetInput;
class TFile;
class TH1F;
class TH2F;
class TProfile;

/// \class RhoSubDebug
///
/// \brief Purely diagnostic module comparing the two rho-subtraction paths:
///   (1) DetermineTowerRho -> SubtractTowersRhov1       -> MULTSUB_TOWERINFO_CALIB_*
///   (2) DetermineTowerRho -> DetermineTowerBackgroundv1 -> SubtractTowers
///                                                       -> TOWERINFO_CALIB_*_SUB1
///
/// Writes nothing to the node tree; fills a standalone histogram file. Four checks:
///   [1] external kT seeds vs the kT clustering DetermineTowerRho does internally,
///       including the effect of the negative-energy flip on the jet axis
///   [2] MULTSUB_TOWERINFO_CALIB_* vs TOWERINFO_CALIB_*_SUB1, tower by tower
///   [3] TowerBackground UE[layer][ieta] vs the amount actually removed by each path
///   [4] the resulting jet collections
class RhoSubDebug : public SubsysReco
{
 public:
  RhoSubDebug(const std::string &name = "RhoSubDebug",
              const std::string &outfile = "rhosub_debug_hists.root");
  ~RhoSubDebug() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_external_kt_node(const std::string &n) { m_ext_kt_node = n; }

  /// Jet-level |eta| acceptance this diagnostic's internal replication of
  /// DetermineTowerRho's clustering uses (drives pass_eta_fj/uf, the acceptance-flip
  /// study, and the seed pool). Must be set to whatever value the real
  /// DetermineTowerRho instances in the macro were given via set_jet_abs_eta(), or
  /// this diagnostic silently reimposes DetermineTowerRho's own default
  /// (tower_abs_eta - R) instead of tracking the real configuration. Left unset
  /// (default), that default resolution still runs in InitRun, unchanged.
  void set_jet_abs_eta(float abseta) { m_abs_jet_eta_range = abseta; }

  /// Seeds of the background method being compared with DetermineTowerRho: jets in
  /// `node` whose prop_SeedItr equals `itr`. Defaults to the external kT container
  /// with itr = 1 (DetermineTowerBackgroundv1). For the iterative HIJetReco method use
  /// AntiKt_TowerInfo_HIRecoSeedsSub_r02 with itr = 2, and optionally its first
  /// iteration via set_bkgd_seed_node_iter1(AntiKt_TowerInfo_HIRecoSeedsRaw_r02, 1).
  void set_bkgd_seed_node(const std::string &node, int itr) { m_bkgd_seed_node = node; m_bkgd_seed_itr = itr; }
  void set_bkgd_seed_node_iter1(const std::string &node, int itr) { m_bkgd_seed_node_it1 = node; m_bkgd_seed_it1 = itr; }
  /// dR within which a background-method seed counts as the same object as a
  /// DetermineTowerRho seed (0.01 when both are R=0.4 kT jets)
  void set_seed_share_dR(float dRmax) { m_seed_share_dR = dRmax; }
  /// Labels written into the histogram file for the plotting macro
  void set_labels(const std::string &path1, const std::string &path2, const std::string &bkgd_seeds,
                  const std::string &bkgd_seeds_it1 = "")
  {
    m_label_path1 = path1; m_label_path2 = path2; m_label_bkgd = bkgd_seeds; m_label_bkgd_it1 = bkgd_seeds_it1;
  }
  void set_background_node(const std::string &n) { m_bkgd_node = n; }
  void set_rho_nodes(const std::string &c, const std::string &i, const std::string &o)
  {
    m_cemc_rho_node = c; m_ihcal_rho_node = i; m_ohcal_rho_node = o;
  }
  void set_jet_nodes(const std::string &a, const std::string &b)
  {
    m_jet_node_rho = a; m_jet_node_sub1 = b;
  }
  void set_vertex_type(GlobalVertex::VTXTYPE t) { m_vertex_type = t; }
  void set_nprint(int n) { m_nprint = n; }

  /// JetContainers written by DetermineTowerRho itself (add_method(..., jet_node)).
  /// When set, the saved jets are validated against this module's own replication
  /// of the internal clustering.
  void set_saved_kt_nodes(const std::string &c, const std::string &i, const std::string &o)
  {
    m_saved_kt_node[1] = c; m_saved_kt_node[2] = i; m_saved_kt_node[3] = o;
  }
  void set_truth_jet_node(const std::string &n) { m_truth_jet_node = n; }
  /// Truth jets used to judge the seeds each background method excludes. Should
  /// have the same R as the kT seeds (0.4). Empty (default) disables the block.
  void set_truth_seed_node(const std::string &n) { m_truth_seed_node = n; }
  void set_seed_truth_cuts(float dRmax, float truth_absetamax)
  {
    m_seed_truth_dR = dRmax; m_truth_seed_absetamax = truth_absetamax;
  }

  /// seed-defining methods compared against truth:
  /// 0 = DetermineTowerBackgroundv1, 1 = DetermineTowerRho all-layer replica,
  /// 2/3/4 = DetermineTowerRho EMCal / IHCal / OHCal instance,
  /// 5 = background method's first iteration (only if set_bkgd_seed_node_iter1)
  static const int NSEEDMETH = 6;
  /// truth jets are required above this pT and inside this |eta| for the JES/JER
  void set_truth_cuts(float ptmin, float absetamax)
  {
    m_truth_ptmin = ptmin; m_truth_absetamax = absetamax;
  }
  /// truth<->reco matching: highest-pT reco jet above reco_ptmin within dRmax
  void set_match_cuts(float reco_ptmin, float dRmax)
  {
    m_reco_match_ptmin = reco_ptmin; m_match_dR = dRmax;
  }

  /// number of layers indexed as 0=EMCal(retower), 1=IHCal, 2=OHCal
  static const int NLAYER = 3;
  /// internal kT collections: 0=all three layers, 1=EMCal, 2=IHCal, 3=OHCal
  static const int NCOLL = 4;

 private:
  /// One kT jet, carrying BOTH axis conventions.
  ///  *_fj  : from the fastjet PseudoJet, i.e. built from the negative-energy-flipped
  ///          constituents. This is what DetermineTowerRho's eta selector actually sees.
  ///  *_uf  : rebuilt by summing the ORIGINAL (unflipped) constituent momenta, which is
  ///          what FastJetAlgoSub stores into the external seed container.
  struct KtJet
  {
    double pt_uf, eta_uf, phi_uf;
    double pt_fj, eta_fj, phi_fj;
    double e_fj;
    int nconst;
    bool pass_eta_fj;  // |eta_fj| <= m_abs_jet_eta_range  <- the cut as implemented
    bool pass_eta_uf;  // |eta_uf| <= m_abs_jet_eta_range  <- the cut on the physical axis
    // DetermineTowerRho's selector is (!SelectorNHardest(N)) * SelectorAbsEtaMax(max).
    // fastjet applies the right-hand selector FIRST: eta cut, then drop the N hardest
    // of what survived. So among jets passing the eta cut:
    bool seed;         // one of the N hardest in acceptance -- omitted from rho
    bool used;         // in acceptance and not a seed -- enters the rho median
  };

  void compare_kt_seeds(PHCompositeNode *topNode);
  void compare_towers(PHCompositeNode *topNode);
  void compare_ue(PHCompositeNode *topNode);
  void compare_jets(PHCompositeNode *topNode);
  void validate_saved_kt(PHCompositeNode *topNode, const std::vector<KtJet> coll[NCOLL]);
  void compare_seeds(PHCompositeNode *topNode, const std::vector<KtJet> coll[NCOLL]);
  struct SeedPos { double pt, eta, phi; };
  void match_seeds_to_truth(PHCompositeNode *topNode, const std::vector<SeedPos> seeds[NSEEDMETH]);
  float grab_mbdQ(PHCompositeNode *topNode);
  void truth_match(PHCompositeNode *topNode);

  /// Replicate DetermineTowerRho's internal clustering. Applies only the
  /// "omit N hardest" part of the selector, so the eta acceptance can be studied.
  std::vector<KtJet> cluster_like_rho(PHCompositeNode *topNode,
                                      const std::vector<JetInput *> &inputs);

  float grab_zvtx(PHCompositeNode *topNode);
  void book_histograms();

  std::string m_outfilename    {"rhosub_debug_hists.root"};
  std::string m_ext_kt_node    {"Kt_TowerInfo_HIRecoSeedsRaw_r04"};
  std::string m_bkgd_node      {"TowerInfoBackground_Rho"};
  std::string m_cemc_rho_node  {"TowerRho_MULT_CEMC"};
  std::string m_ihcal_rho_node {"TowerRho_MULT_HCALIN"};
  std::string m_ohcal_rho_node {"TowerRho_MULT_HCALOUT"};
  std::string m_jet_node_rho   {"AntiKt_TowerInfo_r03_Rho1"};
  std::string m_jet_node_sub1  {"AntiKt_TowerInfo_r03_Rho2"};
  std::string m_saved_kt_node[NCOLL] {};   // [0] unused; [1..3] = per-layer nodes
  std::string m_truth_jet_node {""};
  std::string m_truth_seed_node {""};
  std::string m_bkgd_seed_node {""};      // empty = m_ext_kt_node
  int m_bkgd_seed_itr {1};
  std::string m_bkgd_seed_node_it1 {""};
  int m_bkgd_seed_it1 {1};
  float m_seed_share_dR {0.01};
  std::string m_label_path1 {"Rho1 (MULTSUB)"};
  std::string m_label_path2 {"Rho2 (SUB1)"};
  std::string m_label_bkgd {"DetermineTowerBackgroundv1"};
  std::string m_label_bkgd_it1 {""};
  float m_seed_truth_dR {0.4};
  float m_truth_seed_absetamax {1.1};
  float m_truth_ptmin {5.0};
  float m_truth_absetamax {0.8};
  float m_reco_match_ptmin {3.0};
  float m_match_dR {0.3};

  GlobalVertex::VTXTYPE m_vertex_type {GlobalVertex::MBD};

  // DetermineTowerRho settings that we replicate
  float m_par                 {0.4};
  float m_abs_input_eta_range {1.1};
  float m_abs_jet_eta_range   {-999.0};  // resolved in InitRun to eta_range - par
  unsigned int m_omit_nhardest{2};

  std::vector<JetInput *> m_coll_inputs[NCOLL] {};

  int m_evt {0};
  int m_nprint {3};

  // running summaries
  double m_max_tower_absdiff {0.0};
  double m_sum_tower_absdiff {0.0};
  long   m_n_tower_pairs {0};
  long   m_n_tower_disagree {0};
  double m_max_ue_absdiff {0.0};
  long   m_n_axis_flip_accept {0};   // jets where the two eta conventions disagree on acceptance
  long   m_n_kt_jets {0};
  long   m_n_saved_kt {0};
  long   m_n_saved_kt_matched {0};
  double m_max_saved_kt_dpt {0.0};
  long   m_n_stored_outside_eta {0};
  double m_max_stored_rho_diff {0.0};
  // jet-collection agreement between the two paths
  long   m_n_jets_rho {0};
  long   m_n_jets_sub1 {0};
  long   m_n_jets_matched {0};
  long   m_n_jets_unmatched {0};
  long   m_n_evt_jetcount_differ {0};
  double m_max_jet_dpt {0.0};
  // bitwise agreement, all jets / all towers, no thresholds
  long   m_n_tower_exact_differ {0};
  long   m_n_evt_tower_exact_differ {0};
  long   m_n_jets_all_rho {0};
  long   m_n_jets_all_sub1 {0};
  long   m_n_evt_jet_allcount_differ {0};
  long   m_n_evt_jet_exact_differ {0};
  long   m_n_printed_mismatch {0};
  // stored DetermineTowerRho jets vs replication, both directions
  long   m_n_evt_saved_countdiff {0};
  long   m_n_replica_unmatched {0};
  // seeds
  long   m_n_seed_dtb {0};
  long   m_n_seed_dtb_shared {0};
  long   m_n_evt_seeds_identical {0};

  float m_vtxz {0.0};

  TFile *m_file {nullptr};

  // ---- [1] kT seeds ----
  TH1F *h_kt_pt[NCOLL + 1] {};      // last slot = external
  TH1F *h_kt_eta[NCOLL + 1] {};
  TH1F *h_kt_nconst[NCOLL + 1] {};
  TH1F *h_kt_njet[NCOLL + 1] {};
  TH1F *h_kt_ptovern[NCOLL] {};     // pT / nconst -- the quantity rho is the median of
  TH2F *h_kt_dR_vs_pt {nullptr};    // flipped vs unflipped axis
  TProfile *p_kt_dR_vs_pt {nullptr};
  TH2F *h_kt_deta_vs_pt {nullptr};
  TH1F *h_kt_dR {nullptr};
  TH1F *h_kt_match_dR {nullptr};    // internal(unflipped) vs external seed
  TH1F *h_kt_match_dpt {nullptr};
  TH1F *h_kt_acc_pt_all {nullptr};  // denominator for the acceptance-flip fraction
  TH1F *h_kt_acc_pt_flip {nullptr};
  TH1F *h_rho[NLAYER] {};

  // ---- [2] towers ----
  TH1F *h_tow_dE[NLAYER] {};
  TProfile *p_tow_dE_vs_ieta[NLAYER] {};
  TH1F *h_tow_dE_all {nullptr};

  // ---- [3] UE ----
  TProfile *p_ue_vs_ieta[NLAYER] {};
  TProfile *p_flat_vs_ieta[NLAYER] {};
  TProfile *p_rmMULT_vs_ieta[NLAYER] {};
  TProfile *p_rmSUB1_vs_ieta[NLAYER] {};
  TH1F *h_ue_minus_rmSUB1[NLAYER] {};
  TH1F *h_rmMULT_minus_rmSUB1[NLAYER] {};
  TProfile *p_w_vs_ieta[NLAYER] {};   // implied eta weight UE / (rho cosh eta)

  // ---- [4] jets ----
  TH1F *h_jet_pt_rho {nullptr};
  TH1F *h_jet_pt_sub1 {nullptr};
  TH1F *h_jet_dpt {nullptr};
  TH1F *h_jet_dpt_wide {nullptr};   // same quantity, range wide enough not to clip the tail
  TH1F *h_jet_unmatched_pt {nullptr};
  TH1F *h_jet_njet_rho {nullptr};
  TH1F *h_jet_njet_sub1 {nullptr};
  TH2F *h_jet_pt_corr {nullptr};
  TH1F *h_jet_eta[2] {};
  TH1F *h_jet_phi[2] {};
  TH2F *h_jet_etaphi[2] {};

  // ---- [5] truth-matched JES / JER ----
  TH1F *h_truth_pt {nullptr};
  TH1F *h_truth_pt_matched[2] {};
  TH2F *h_jes[2] {};
  TH1F *h_truth_eta {nullptr};

  // ---- saved-jet validation ----
  TH1F *h_saved_kt_dpt {nullptr};
  TH1F *h_saved_kt_dR {nullptr};

  // ---- seeds: jets excluded from the background by each path ----
  // DetermineTowerRho: the N hardest kT jets inside its eta acceptance, per rho instance
  // DetermineTowerBackgroundv1: external kT jets it tagged prop_SeedItr == 1
  TH2F *h_seed_etaphi_rho[NCOLL] {};   // [0] all-layer replica, [1..3] per layer
  TH1F *h_seed_pt_rho[NCOLL] {};
  TH1F *h_seed_eta_rho[NCOLL] {};
  TH2F *h_seed_etaphi_dtb {nullptr};
  TH1F *h_seed_pt_dtb {nullptr};
  TH1F *h_seed_eta_dtb {nullptr};
  TH1F *h_seed_nseed_dtb {nullptr};
  TH1F *h_seed_match_dR {nullptr};     // each DTBv1 seed -> nearest all-layer DetermineTowerRho seed
  TH1F *h_seed_match_dR_emcal {nullptr};
  TH1F *h_seed_nshared {nullptr};      // per event: DTBv1 seeds that are also DetermineTowerRho seeds
  TH2F *h_kt_etaphi_ext {nullptr};
  TH2F *h_stored_etaphi[NLAYER] {};

  // ---- seeds vs the two leading truth jets (|eta| < m_truth_seed_absetamax) ----
  // [i] = leading (0) / subleading (1) truth jet
  TH1F *h_st_truth_pt[2] {};
  TH1F *h_st_truth_eta[2] {};
  TH1F *h_st_tag_pt[NSEEDMETH][2] {};    // truth jet has a seed of this method within dR
  TH1F *h_st_tag_eta[NSEEDMETH][2] {};
  TH1F *h_st_nfound[NSEEDMETH] {};       // per event: leading truth jets found (0,1,2)
  TH1F *h_st_seed_dR[NSEEDMETH] {};      // each seed -> nearest of the two leading truth jets
  TH1F *h_st_seed_pt_all[NSEEDMETH] {};
  TH1F *h_st_seed_pt_matched[NSEEDMETH] {};
  long m_st_events {0};
  long m_st_nfound[NSEEDMETH][3] {};
  long m_st_nseed[NSEEDMETH] {};
  long m_st_nseed_matched[NSEEDMETH] {};

  // ---- exact (bitwise) agreement ----
  TH1F *h_tow_nexact_diff {nullptr};   // per event: towers whose energies are not bit-identical
  TH1F *h_jet_ncount_diff {nullptr};   // per event: N_jets(Rho1) - N_jets(Rho2), all jets, no pT cut
  TH1F *h_saved_kt_ncount_diff {nullptr};  // per event: stored - replicated, summed over layers

  // ---- distributions filled directly FROM the stored DetermineTowerRho containers ----
  TH1F *h_stored_pt[NLAYER] {};
  TH1F *h_stored_eta[NLAYER] {};
  TH1F *h_stored_nconst[NLAYER] {};
  TH1F *h_stored_ptovern[NLAYER] {};
  TH1F *h_stored_njet[NLAYER] {};
  TH1F *h_stored_rho[NLAYER] {};
  TH1F *h_stored_rho_minus_node[NLAYER] {};

  static const char *layer_name(int l)
  {
    return (l == 0) ? "EMCAL(retower)" : ((l == 1) ? "IHCAL" : "OHCAL");
  }
  static const char *layer_tag(int l)
  {
    return (l == 0) ? "emcal" : ((l == 1) ? "ihcal" : "ohcal");
  }
  static const char *seedmeth_tag(int m)
  {
    static const char *t[NSEEDMETH] = {"dtb", "rho_all", "rho_emcal", "rho_ihcal", "rho_ohcal", "dtb_it1"};
    return t[m];
  }
  std::string seedmeth_name(int m) const
  {
    static const char *t[NSEEDMETH] = {"", "DetermineTowerRho all layers", "DetermineTowerRho EMCal",
                                       "DetermineTowerRho IHCal", "DetermineTowerRho OHCal", ""};
    if (m == 0) { return m_label_bkgd; }
    if (m == 5) { return m_label_bkgd_it1.empty() ? std::string("first iteration") : m_label_bkgd_it1; }
    return t[m];
  }
  static const char *coll_tag(int c)
  {
    return (c == 0) ? "all" : ((c == 1) ? "emcal" : ((c == 2) ? "ihcal" : "ohcal"));
  }
  static const char *coll_name(int c)
  {
    return (c == 0) ? "internal all 3 layers"
                    : ((c == 1) ? "internal EMCal only"
                                : ((c == 2) ? "internal IHCal only" : "internal OHCal only"));
  }
};

#endif
