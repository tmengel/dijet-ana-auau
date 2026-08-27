#ifndef _PARSE_TREE_C_
#define _PARSE_TREE_C_

#include <myana/AnaUtils.h>

#include <cdbobjects/CDBTTree.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TString.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )
R__LOAD_LIBRARY( libcdbobjects.so )

// 1) throws random cones on top of the (already rho-subtracted) calorimeter
//    towers to measure the residual background fluctuation left after
//    subtraction, binned in centrality. The same cone positions are also
//    summed over the constituents of every sub1(seeded)-subtracted jet in
//    the event (i.e. the same underlying towers, reached via the jet
//    constituent lists instead of the tower grid), so the two subtraction
//    methods' residual fluctuations can be compared cone-by-cone. A third
//    cone is summed over rho_jet's constituents as well -- but since
//    rho_jet's rho subtraction is applied only at the jet level (not per
//    tower), its constituent_E is RAW/unsubtracted, so that third cone is
//    just an independent sanity check on the raw underlying-event energy
//    density, not a subtraction cross-check like the sub1 one.
// 2) compares general kinematics of the rho-subtracted (rho_jet) and
//    seeded/iterative-subtracted (sub1_jet) jet collections, and greedily
//    dR-matches the two collections to see how each method treats the same
//    underlying jet (pT correlation, dR/deta/dphi of the match).
// 3) starting from the RAW towers recovered via the rho_jet constituent
//    lists (see (1) above), independently re-derives three subtracted cone
//    energies by hand: a "manual sub1" using the per-(layer,eta ring) <UE>
//    from the TowerBackgroundv1 object (sub2_towerbkgd_ue), an
//    uncalibrated-rho subtraction using the raw per-event rho density, and
//    a calibrated-rho subtraction using the eta/zvtx/mbdQ-binned weights
//    from the rho_calib pipeline (macros/rho_calib/make_rho_calib.C). These
//    are compared against each other and against the two "baked-in"
//    subtractions from (1) (sub1_jet constituents and the tower grid) as a
//    from-scratch cross-check of both the sub1 and rho subtraction chains.
//
// events failing the sPHENIX minimum-bias requirement or outside
// |zvrtx| < max_abs_zvrtx are discarded, matching the convention used in
// macros/rho_calib/read_tree.C.
int parse_tree(
    const std::string & infile = "/sphenix/tg/tg01/jets/tmengel/jstgtf03/data/v003_20260821/anatree/anatree_run2auau_pro001_pcdb001_v001-00054912-00000.root",
    const std::string & outfile = "parse_tree.root",
    const std::string & rho_calib_file = "/sphenix/user/tmengel/dijet-ana-auau/macros/rho_calib/calibs/rho_calib_54912.root",
    const float max_abs_zvrtx = 60.0,
    const float min_jet_pt = 5.0,
    const int n_cones_per_event = 2,
    const int rng_seed = 12345
)
{
    auto * f = TFile::Open( infile.c_str() );
    if ( !f || f -> IsZombie() )
    {
        std::cerr << "Error: could not open " << infile << std::endl;
        return -1;
    }

    auto * t = (TTree*) f -> Get( "T" );
    if ( !t )
    {
        std::cerr << "Error: no tree \"T\" in " << infile << std::endl;
        return -1;
    }

    const int nentries = t -> GetEntries();
    if ( nentries < 1 )
    {
        std::cerr << "Error: no entries found in tree." << std::endl;
        return -1;
    }
    std::cout << "Total entries: " << nentries << std::endl;

    bool ok = true;
    auto enable = [&]( const char * name, void * addr )
    {
        if ( !t -> GetBranch( name ) )
        {
            std::cerr << "Error: branch \"" << name << "\" not found in input chain." << std::endl;
            ok = false;
            return;
        }
        t -> SetBranchStatus( name, true );
        t -> SetBranchAddress( name, addr );
    };


    int scaled_triggervec[64];
    int is_minbias;
    float zvrtx;
    int cent;
    float mbd_q_N, mbd_q_S;
    enable( "scaled_triggervec", scaled_triggervec );
    enable( "is_minbias", &is_minbias );
    enable( "zvrtx", &zvrtx );
    enable( "cent", &cent );
    enable( "mbd_q_N", &mbd_q_N );
    enable( "mbd_q_S", &mbd_q_S );

    // these are actually the rho subtracted towers
    // the have been subtracted by E' = E - rho' cosh(eta_corr)
    // where rho' = rho * w(eta) and eta_corr = eta - delta_eta
    float sumeT_cemc, sumeT_ihcal, sumeT_ohcal;
    float cemc_tower_E[24][64], ihcal_tower_E[24][64], ohcal_tower_E[24][64];
    int cemc_tower_isgood[24][64], ihcal_tower_isgood[24][64], ohcal_tower_isgood[24][64];
    enable( "sumeT_cemc", &sumeT_cemc );
    enable( "sumeT_ihcal", &sumeT_ihcal );
    enable( "sumeT_ohcal", &sumeT_ohcal );
    enable( "cemc_tower_E", cemc_tower_E );
    enable( "ihcal_tower_E", ihcal_tower_E );
    enable( "ohcal_tower_E", ohcal_tower_E );
    enable( "cemc_tower_isgood", cemc_tower_isgood );
    enable( "ihcal_tower_isgood", ihcal_tower_isgood );
    enable( "ohcal_tower_isgood", ohcal_tower_isgood );

    float sub1_jet_R;
    std::vector<float> * sub1_jet_E = nullptr;
    std::vector<float> * sub1_jet_pT = nullptr;
    std::vector<float> * sub1_jet_phi = nullptr;
    std::vector<float> * sub1_jet_eta = nullptr;
    std::vector<float> * sub1_jet_unsub_pT = nullptr;
    std::vector<float> * sub1_jet_unsub_E = nullptr;
    std::vector< std::vector<float> > * sub1_jet_constituent_E = nullptr;
    std::vector< std::vector<float> > * sub1_jet_constituent_pT = nullptr;
    std::vector< std::vector<float> > * sub1_jet_constituent_phi = nullptr;
    std::vector< std::vector<float> > * sub1_jet_constituent_eta = nullptr;
    std::vector< std::vector<int> > * sub1_jet_constituent_srcID = nullptr;

    std::vector< std::vector<float > > * sub2_towerbkgd_ue = nullptr;

    enable( "sub1_jet_R", &sub1_jet_R );
    enable( "sub1_jet_E", &sub1_jet_E );
    enable( "sub1_jet_pT", &sub1_jet_pT );
    enable( "sub1_jet_phi", &sub1_jet_phi );
    enable( "sub1_jet_eta", &sub1_jet_eta );
    enable( "sub1_jet_unsub_pT", &sub1_jet_unsub_pT );
    enable( "sub1_jet_unsub_E", &sub1_jet_unsub_E );
    enable( "sub1_jet_constituent_E", &sub1_jet_constituent_E );
    enable( "sub1_jet_constituent_eta", &sub1_jet_constituent_eta );
    enable( "sub1_jet_constituent_phi", &sub1_jet_constituent_phi );
    enable( "sub2_towerbkgd_ue", &sub2_towerbkgd_ue );

    float rho_jet_R;
    std::vector<float> * rho_jet_E = nullptr;
    std::vector<float> * rho_jet_pT = nullptr;
    std::vector<float> * rho_jet_phi = nullptr;
    std::vector<float> * rho_jet_eta = nullptr;
    std::vector<float> * rho_jet_unsub_pT = nullptr;
    std::vector<float> * rho_jet_unsub_E = nullptr;
    std::vector< std::vector<float> > * rho_jet_constituent_E = nullptr;
    std::vector< std::vector<float> > * rho_jet_constituent_pT = nullptr;
    std::vector< std::vector<float> > * rho_jet_constituent_phi = nullptr;
    std::vector< std::vector<float> > * rho_jet_constituent_eta = nullptr;
    std::vector< std::vector<int> > * rho_jet_constituent_srcID = nullptr;

    enable( "rho_jet_R", &rho_jet_R );
    enable( "rho_jet_E", &rho_jet_E );
    enable( "rho_jet_pT", &rho_jet_pT );
    enable( "rho_jet_phi", &rho_jet_phi );
    enable( "rho_jet_eta", &rho_jet_eta );
    enable( "rho_jet_unsub_pT", &rho_jet_unsub_pT );
    enable( "rho_jet_unsub_E", &rho_jet_unsub_E );
    enable( "rho_jet_constituent_E", &rho_jet_constituent_E );
    enable( "rho_jet_constituent_eta", &rho_jet_constituent_eta );
    enable( "rho_jet_constituent_phi", &rho_jet_constituent_phi );
    enable( "rho_jet_constituent_srcID", &rho_jet_constituent_srcID );

    // uncalibrated per-event, per-layer rho density (ET/area-like units --
    // see the ue*=cosh(eta) convention in AnaTreev1::rho_jet's constituent
    // loop), used for the manual rho subtraction below.
    float rho_val_cemc, rho_val_hcalin, rho_val_hcalout;
    enable( "rho_val_TowerRho_MULT_CEMC", &rho_val_cemc );
    enable( "rho_val_TowerRho_MULT_HCALIN", &rho_val_hcalin );
    enable( "rho_val_TowerRho_MULT_HCALOUT", &rho_val_hcalout );

    if ( !ok )
    {
        std::cerr << "Error: one or more required branches were missing, aborting." << std::endl;
        return -1;
    }

    //--------------------------------------------------------------------
    // rho calibration (eta/zvtx/mbdQ-binned weights, see
    // macros/rho_calib/make_rho_calib.C / SubtractTowersRhov1)
    //--------------------------------------------------------------------
    auto * calib_tree = new CDBTTree( rho_calib_file );
    calib_tree -> LoadCalibrations();
    const int calib_n_eta     = calib_tree -> GetSingleIntValue( "n_eta" );
    const int calib_n_zbins   = calib_tree -> GetSingleIntValue( "n_zvtx_bins" );
    const int calib_n_mbdbins = calib_tree -> GetSingleIntValue( "n_mbdQ_bins" );
    if ( calib_n_eta != 24 || calib_n_zbins < 1 || calib_n_mbdbins < 1 )
    {
        std::cerr << "Error: could not load a valid rho calibration from " << rho_calib_file << std::endl;
        return -1;
    }
    std::vector<float> calib_zvtx_edges( calib_n_zbins + 1 );
    for ( int i = 0; i <= calib_n_zbins; ++i )
    {
        calib_zvtx_edges[i] = calib_tree -> GetSingleFloatValue( Form( "zvtx_edge_%d", i ) );
    }
    std::vector<float> calib_mbdQ_edges( calib_n_mbdbins + 1 );
    for ( int i = 0; i <= calib_n_mbdbins; ++i )
    {
        calib_mbdQ_edges[i] = calib_tree -> GetSingleFloatValue( Form( "mbdQ_edge_%d", i ) );
    }
    const char * calib_field_names[3] = { "w_cemc", "w_hcalin", "w_hcalout" };

    // matches SubtractTowersRhov1::find_bin -- half-open [edges[i],edges[i+1]),
    // -1 if val is outside [edges.front(),edges.back())
    auto find_bin = []( const float val, const std::vector<float> & edges ) -> int
    {
        if ( edges.size() < 2 || std::isnan(val) || val < edges.front() || val >= edges.back() ) return -1;
        for ( size_t i = 0; i + 1 < edges.size(); ++i )
        {
            if ( val >= edges[i] && val < edges[i + 1] ) return static_cast<int>(i);
        }
        return -1;
    };
    // matches SubtractTowersRhov1::encode_channel
    auto encode_channel = []( const int ieta, const int izbin, const int imbd, const int n_mbd_bins ) -> int
    {
        return izbin * ( n_mbd_bins * 24 ) + imbd * 24 + ieta;
    };
    // rho_jet clusters from plain raw calibrated towers, not the "_SUB1"
    // containers (confirmed empirically: rho_jet_constituent_srcID only
    // ever takes the values 26/27/28, not 29/30/31). Jet::SRC values:
    // CEMC_TOWERINFO_RETOWER=28, HCALIN_TOWERINFO=26, HCALOUT_TOWERINFO=27
    // -> layer {0=CEMC,1=HCALIN,2=HCALOUT}.
    auto decode_layer = []( const int srcID ) -> int
    {
        if ( srcID == 28 ) return 0; // CEMC (retowered)
        if ( srcID == 26 ) return 1; // HCALIN
        if ( srcID == 27 ) return 2; // HCALOUT
        return -1;
    };

    std::cout << "Loaded rho calibration " << rho_calib_file << " (n_eta=" << calib_n_eta
              << ", n_zvtx_bins=" << calib_n_zbins << ", n_mbdQ_bins=" << calib_n_mbdbins << ")" << std::endl;

    //--------------------------------------------------------------------
    // histograms
    //--------------------------------------------------------------------

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // centrality bins: 0-10, 10-20, ..., 80-90 (percent)
    const int NCENT = 9;
    const float cent_edges[NCENT + 1] = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 };
    auto find_cent_bin = [&]( const float c ) -> int
    {
        if ( c < cent_edges[0] || c >= cent_edges[NCENT] ) return -1;
        for ( int i = 0; i < NCENT; ++i )
        {
            if ( c >= cent_edges[i] && c < cent_edges[i + 1] ) return i;
        }
        return -1;
    };

    // -------------------- random cones --------------------
    // "tower": cones summed over the rho-subtracted calorimeter tower grid
    // (cemc/ihcal/ohcal_tower_E). "sub1": the same cone positions, but summed
    // over the constituents of every sub1(seeded)-subtracted jet in the event
    // instead -- i.e. the same underlying (sub1-subtracted) towers, just
    // reached via the jet constituent lists rather than the tower grid.
    TH1D * h_rc_pt_all = new TH1D( "h_rc_pt_all", "Random cone p_{T} (rho-subtracted towers);cone p_{T} [GeV];cones", 200, -50, 50 );
    TH1D * h_rc_pt_cent[NCENT];
    TH1D * h_rc_pt_sub1_all = new TH1D( "h_rc_pt_sub1_all", "Random cone p_{T} (sub1-subtracted jet constituents);cone p_{T} [GeV];cones", 200, -50, 50 );
    TH1D * h_rc_pt_sub1_cent[NCENT];
    // "rho_const": same recipe as "sub1", but summed over the constituents of
    // every rho_jet instead. Unlike sub1_jet, rho_jet's rho subtraction is
    // applied only at the jet level (E_jet = E_unsub - rho*Area); its
    // constituent_E is the RAW, unsubtracted tower energy (see
    // AnaTreev1::rho_jet constituent loop, comp_E = tower->get_energy()
    // with no rho term subtracted). So this is NOT a subtraction cross-check
    // like the sub1 one above -- it's an independent sanity check on the
    // typical raw (unsubtracted) underlying-event energy in a cone.
    TH1D * h_rc_pt_rhoconst_all = new TH1D( "h_rc_pt_rhoconst_all", "Random cone p_{T}, RAW/unsubtracted (rho_jet constituents);cone p_{T} [GeV];cones", 200, -50, 150 );
    TH1D * h_rc_pt_rhoconst_cent[NCENT];
    // manual re-derivations from the same RAW towers, cross-checking the
    // sub1 and rho subtraction chains independently of what's baked into
    // sub1_jet_constituent_E / cemc_tower_E:
    //   manualsub1  : E_raw - <UE>(layer,eta ring)                 [TowerBackgroundv1]
    //   rhouncalib  : E_raw - rho(layer) * cosh(eta)                [raw event rho, no eta weight]
    //   rhocalib    : E_raw - rho(layer) * w(mbdQ,zvtx,eta) * cosh(eta)  [rho_calib CDB weights]
    TH1D * h_rc_pt_manualsub1_all = new TH1D( "h_rc_pt_manualsub1_all", "Random cone p_{T}, manual sub1 (RAW - <UE>);cone p_{T} [GeV];cones", 200, -50, 50 );
    TH1D * h_rc_pt_manualsub1_cent[NCENT];
    TH1D * h_rc_pt_rhouncalib_all = new TH1D( "h_rc_pt_rhouncalib_all", "Random cone p_{T}, uncalibrated rho (RAW - #rho#upointcosh#eta);cone p_{T} [GeV];cones", 200, -50, 50 );
    TH1D * h_rc_pt_rhouncalib_cent[NCENT];
    TH1D * h_rc_pt_rhocalib_all = new TH1D( "h_rc_pt_rhocalib_all", "Random cone p_{T}, calibrated rho (RAW - #rho#upoint w#upointcosh#eta);cone p_{T} [GeV];cones", 200, -50, 50 );
    TH1D * h_rc_pt_rhocalib_cent[NCENT];
    // same as rhouncalib/rhocalib, but a tower with eT > sqrt(2)*rho' is left
    // unsubtracted (treated as a jet fragment/seed, not background).
    TH1D * h_rc_pt_rhouncalib_seed_all = new TH1D( "h_rc_pt_rhouncalib_seed_all", "Random cone p_{T}, uncalibrated rho, seed-protected;cone p_{T} [GeV];cones", 200, -50, 50 );
    TH1D * h_rc_pt_rhouncalib_seed_cent[NCENT];
    TH1D * h_rc_pt_rhocalib_seed_all = new TH1D( "h_rc_pt_rhocalib_seed_all", "Random cone p_{T}, calibrated rho, seed-protected;cone p_{T} [GeV];cones", 200, -50, 50 );
    TH1D * h_rc_pt_rhocalib_seed_cent[NCENT];
    for ( int i = 0; i < NCENT; ++i )
    {
        h_rc_pt_cent[i] = new TH1D(
            Form( "h_rc_pt_cent%d_%d", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            Form( "Random cone p_{T}, cent %d-%d%%;cone p_{T} [GeV];cones", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            200, -50, 50
        );
        h_rc_pt_sub1_cent[i] = new TH1D(
            Form( "h_rc_pt_sub1_cent%d_%d", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            Form( "Random cone p_{T} (sub1 constituents), cent %d-%d%%;cone p_{T} [GeV];cones", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            200, -50, 50
        );
        h_rc_pt_rhoconst_cent[i] = new TH1D(
            Form( "h_rc_pt_rhoconst_cent%d_%d", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            Form( "Random cone p_{T}, RAW/unsubtracted (rho_jet constituents), cent %d-%d%%;cone p_{T} [GeV];cones", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            200, -50, 150
        );
        h_rc_pt_manualsub1_cent[i] = new TH1D(
            Form( "h_rc_pt_manualsub1_cent%d_%d", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            Form( "Random cone p_{T}, manual sub1, cent %d-%d%%;cone p_{T} [GeV];cones", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            200, -50, 50
        );
        h_rc_pt_rhouncalib_cent[i] = new TH1D(
            Form( "h_rc_pt_rhouncalib_cent%d_%d", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            Form( "Random cone p_{T}, uncalibrated rho, cent %d-%d%%;cone p_{T} [GeV];cones", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            200, -50, 50
        );
        h_rc_pt_rhocalib_cent[i] = new TH1D(
            Form( "h_rc_pt_rhocalib_cent%d_%d", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            Form( "Random cone p_{T}, calibrated rho, cent %d-%d%%;cone p_{T} [GeV];cones", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            200, -50, 50
        );
        h_rc_pt_rhouncalib_seed_cent[i] = new TH1D(
            Form( "h_rc_pt_rhouncalib_seed_cent%d_%d", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            Form( "Random cone p_{T}, uncalibrated rho, seed-protected, cent %d-%d%%;cone p_{T} [GeV];cones", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            200, -50, 50
        );
        h_rc_pt_rhocalib_seed_cent[i] = new TH1D(
            Form( "h_rc_pt_rhocalib_seed_cent%d_%d", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            Form( "Random cone p_{T}, calibrated rho, seed-protected, cent %d-%d%%;cone p_{T} [GeV];cones", (int) cent_edges[i], (int) cent_edges[i + 1] ),
            200, -50, 50
        );
    }

    // -------------------- random cone: tower vs sub1-constituent comparison --------------------
    TH2D * h_rc_pt_corr_tower_vs_sub1 = new TH2D(
        "h_rc_pt_corr_tower_vs_sub1", "random cones;tower (rho-subtracted) cone p_{T} [GeV];sub1-subtracted constituent cone p_{T} [GeV]",
        100, -50, 50, 100, -50, 50
    );
    TH1D * h_rc_dpt_sub1_minus_tower = new TH1D(
        "h_rc_dpt_sub1_minus_tower", "random cones;p_{T}^{cone,sub1} - p_{T}^{cone,tower} [GeV];cones", 100, -30, 30
    );

    // -------------------- random cone: tower(rho-subtracted) vs rho_jet-constituent(RAW) sanity check --------------------
    TH2D * h_rc_pt_corr_tower_vs_rhoconst = new TH2D(
        "h_rc_pt_corr_tower_vs_rhoconst", "random cones;tower (rho-subtracted) cone p_{T} [GeV];RAW/unsubtracted (rho_jet constituent) cone p_{T} [GeV]",
        100, -50, 50, 100, -50, 150
    );
    TH1D * h_rc_dpt_rhoconst_minus_tower = new TH1D(
        "h_rc_dpt_rhoconst_minus_tower", "random cones;p_{T}^{cone,RAW} - p_{T}^{cone,tower(rho)} [GeV];cones", 100, -30, 130
    );

    // -------------------- validation: manual (from-scratch) vs baked-in subtractions --------------------
    TH2D * h_rc_pt_corr_manualsub1_vs_sub1const = new TH2D(
        "h_rc_pt_corr_manualsub1_vs_sub1const", "random cones;manual sub1 (RAW - <UE>) cone p_{T} [GeV];sub1_jet constituent cone p_{T} [GeV]",
        100, -50, 50, 100, -50, 50
    );
    TH1D * h_rc_dpt_manualsub1_minus_sub1const = new TH1D(
        "h_rc_dpt_manualsub1_minus_sub1const", "random cones;p_{T}^{cone,manual sub1} - p_{T}^{cone,sub1 const} [GeV];cones", 100, -20, 20
    );

    TH2D * h_rc_pt_corr_rhocalib_vs_tower = new TH2D(
        "h_rc_pt_corr_rhocalib_vs_tower", "random cones;manual calibrated-rho cone p_{T} [GeV];tower (rho-subtracted) cone p_{T} [GeV]",
        100, -50, 50, 100, -50, 50
    );
    TH1D * h_rc_dpt_rhocalib_minus_tower = new TH1D(
        "h_rc_dpt_rhocalib_minus_tower", "random cones;p_{T}^{cone,manual rho-calib} - p_{T}^{cone,tower} [GeV];cones", 100, -20, 20
    );

    TH2D * h_rc_pt_corr_rhouncalib_vs_rhocalib = new TH2D(
        "h_rc_pt_corr_rhouncalib_vs_rhocalib", "random cones;uncalibrated-rho cone p_{T} [GeV];calibrated-rho cone p_{T} [GeV]",
        100, -50, 50, 100, -50, 50
    );

    // -------------------- effect of seed/outlier protection (eT > sqrt(2)*rho' left unsubtracted) --------------------
    TH1D * h_rc_dpt_rhouncalib_seed_minus_rhouncalib = new TH1D(
        "h_rc_dpt_rhouncalib_seed_minus_rhouncalib", "random cones;p_{T}^{cone,uncalib rho,seed} - p_{T}^{cone,uncalib rho} [GeV];cones", 100, -10, 30
    );
    TH1D * h_rc_dpt_rhocalib_seed_minus_rhocalib = new TH1D(
        "h_rc_dpt_rhocalib_seed_minus_rhocalib", "random cones;p_{T}^{cone,calib rho,seed} - p_{T}^{cone,calib rho} [GeV];cones", 100, -10, 30
    );
    TH2D * h_rc_pt_corr_rhocalib_seed_vs_tower = new TH2D(
        "h_rc_pt_corr_rhocalib_seed_vs_tower", "random cones;manual calibrated-rho, seed-protected cone p_{T} [GeV];tower (rho-subtracted) cone p_{T} [GeV]",
        100, -50, 50, 100, -50, 50
    );
    TH1D * h_rc_dpt_rhocalib_seed_minus_tower = new TH1D(
        "h_rc_dpt_rhocalib_seed_minus_tower", "random cones;p_{T}^{cone,manual rho-calib,seed} - p_{T}^{cone,tower} [GeV];cones", 100, -20, 20
    );

    // -------------------- rho vs sub1 jet kinematics --------------------
    TH1D * h_rho_jet_pt   = new TH1D( "h_rho_jet_pt", "rho-subtracted jets;p_{T} [GeV];jets", 100, 0, 100 );
    TH1D * h_sub1_jet_pt  = new TH1D( "h_sub1_jet_pt", "sub1-subtracted jets;p_{T} [GeV];jets", 100, 0, 100 );
    TH1D * h_rho_jet_eta  = new TH1D( "h_rho_jet_eta", "rho-subtracted jets;#eta;jets", 44, -1.1, 1.1 );
    TH1D * h_sub1_jet_eta = new TH1D( "h_sub1_jet_eta", "sub1-subtracted jets;#eta;jets", 44, -1.1, 1.1 );
    TH1D * h_rho_jet_phi  = new TH1D( "h_rho_jet_phi", "rho-subtracted jets;#phi;jets", 64, -TMath::Pi(), TMath::Pi() );
    TH1D * h_sub1_jet_phi = new TH1D( "h_sub1_jet_phi", "sub1-subtracted jets;#phi;jets", 64, -TMath::Pi(), TMath::Pi() );
    TH2D * h_rho_jet_pt_eta  = new TH2D( "h_rho_jet_pt_eta", "rho-subtracted jets;#eta;p_{T} [GeV]", 44, -1.1, 1.1, 100, 0, 100 );
    TH2D * h_sub1_jet_pt_eta = new TH2D( "h_sub1_jet_pt_eta", "sub1-subtracted jets;#eta;p_{T} [GeV]", 44, -1.1, 1.1, 100, 0, 100 );

    TH1D * h_njets_rho  = new TH1D( "h_njets_rho", "selected jet multiplicity;N_{jets};events", 20, 0, 20 );
    TH1D * h_njets_sub1 = new TH1D( "h_njets_sub1", "selected jet multiplicity;N_{jets};events", 20, 0, 20 );

    // -------------------- rho <-> sub1 matching --------------------
    TH2D * h_match_pt_rho_vs_sub1 = new TH2D(
        "h_match_pt_rho_vs_sub1", "matched jets;rho-subtracted p_{T} [GeV];sub1-subtracted p_{T} [GeV]",
        100, 0, 100, 100, 0, 100
    );
    TH1D * h_match_dpt        = new TH1D( "h_match_dpt", "matched jets;p_{T}^{sub1} - p_{T}^{rho} [GeV];pairs", 100, -30, 30 );
    TH1D * h_match_dpt_relpt  = new TH1D( "h_match_dpt_relpt", "matched jets;(p_{T}^{sub1} - p_{T}^{rho}) / p_{T}^{rho};pairs", 100, -2, 2 );
    TH1D * h_match_dr         = new TH1D( "h_match_dr", "matched jets;#Delta R(rho, sub1);pairs", 60, 0, 0.6 );
    TH1D * h_match_deta       = new TH1D( "h_match_deta", "matched jets;#eta^{sub1} - #eta^{rho};pairs", 60, -0.3, 0.3 );
    TH1D * h_match_dphi       = new TH1D( "h_match_dphi", "matched jets;#phi^{sub1} - #phi^{rho};pairs", 60, -0.3, 0.3 );

    long n_matched = 0, n_rho_only = 0, n_sub1_only = 0;
    long n_cones_thrown = 0;

    TRandom3 rnd( rng_seed );

    // sum of transverse energy (already rho-subtracted) of good towers within
    // dR < R of (cone_eta, cone_phi), using the z-vertex corrected tower eta
    // for each calorimeter layer.
    auto sum_cone_pt = [&]( const float cone_eta, const float cone_phi, const float R ) -> float
    {
        float sum = 0.0;
        struct Layer { AnaUtils::CaloType calo; const float (*E)[64]; const int (*isgood)[64]; };
        Layer layers[3] = {
            { AnaUtils::CEMC,    cemc_tower_E,  cemc_tower_isgood  },
            { AnaUtils::HCALIN,  ihcal_tower_E, ihcal_tower_isgood },
            { AnaUtils::HCALOUT, ohcal_tower_E, ohcal_tower_isgood }
        };
        for ( const auto & layer : layers )
        {
            for ( int ieta = 0; ieta < 24; ++ieta )
            {
                const float tower_eta = AnaUtils::get_corrected_calo_eta( layer.calo, ieta, zvrtx );
                if ( AnaUtils::deta_abs( tower_eta, cone_eta ) > R ) continue; // cheap pre-cut
                const float inv_cosh_eta = 1.0f / std::cosh( tower_eta );
                for ( int iphi = 0; iphi < 64; ++iphi )
                {
                    if ( !layer.isgood[ieta][iphi] ) continue;
                    const float tower_phi = AnaUtils::get_calo_phi( layer.calo, iphi );
                    if ( AnaUtils::calc_dr( tower_eta, tower_phi, cone_eta, cone_phi ) >= R ) continue;
                    sum += layer.E[ieta][iphi] * inv_cosh_eta;
                }
            }
        }
        return sum;
    };

    // sum of transverse energy of sub1(seeded)-subtracted jet constituents
    // within dR < R of (cone_eta, cone_phi). Looping over every jet's
    // constituent list (regardless of jet pT) recovers essentially every
    // sub1-subtracted tower in the event -- the anti-kt input -- just via
    // the jet constituent bookkeeping instead of the tower grid.
    auto sum_cone_pt_sub1 = [&]( const float cone_eta, const float cone_phi, const float R ) -> float
    {
        float sum = 0.0;
        const size_t njets = sub1_jet_constituent_E -> size();
        for ( size_t ij = 0; ij < njets; ++ij )
        {
            const auto & cE   = sub1_jet_constituent_E   -> at(ij);
            const auto & cEta = sub1_jet_constituent_eta -> at(ij);
            const auto & cPhi = sub1_jet_constituent_phi -> at(ij);
            const size_t ncon = cE.size();
            for ( size_t ic = 0; ic < ncon; ++ic )
            {
                if ( AnaUtils::deta_abs( cEta[ic], cone_eta ) > R ) continue; // cheap pre-cut
                if ( AnaUtils::calc_dr( cEta[ic], cPhi[ic], cone_eta, cone_phi ) >= R ) continue;
                sum += cE[ic] / std::cosh( cEta[ic] );
            }
        }
        return sum;
    };

    // one combined pass over the RAW rho_jet constituents in a cone,
    // computing the raw sum (same quantity as h_rc_pt_rhoconst_all) plus
    // all three manual subtractions at once. The compute_cone_sums lambda
    // is (re)built fresh inside the event loop below, since it captures
    // per-event lookup tables (eta ring positions, calibration weights,
    // raw rho).
    struct ConeSums
    {
        float raw              = 0.0f;
        float sub1             = 0.0f;
        float rho_uncalib      = 0.0f;
        float rho_calib        = 0.0f;
        // same as rho_uncalib/rho_calib, but a tower is left unsubtracted
        // (kept at its raw eT) if its eT exceeds sqrt(2)*rho' -- treating it
        // as a jet fragment ("seed") rather than background, the way seeded
        // iterative subtraction methods do.
        float rho_uncalib_seed = 0.0f;
        float rho_calib_seed   = 0.0f;
    };

    std::cout << "Minimum jet pT: " << min_jet_pt << std::endl;

    for ( int ientry = 0; ientry < nentries; ++ientry )
    {
        t -> GetEntry( ientry );
        if ( nentries >= 10 && ientry % ( nentries / 10 ) == 0 && ientry > 0 )
        {
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        }

        if ( !is_minbias ) continue;
        if ( std::fabs( zvrtx ) >= max_abs_zvrtx ) continue;

        const int icent = find_cent_bin( (float) cent );

        //----------------------------------------------------------------
        // per-event setup for the manual (from-scratch) subtractions:
        // z-vertex-corrected eta ring positions (to map a constituent's
        // continuous eta back to its ieta ring), the calibration bin
        // (izbin,imbd) and per-(layer,ieta) weights, and the raw per-layer
        // rho density.
        //----------------------------------------------------------------
        float eta_ring[3][24];
        for ( int layer = 0; layer < 3; ++layer )
        {
            for ( int ieta = 0; ieta < 24; ++ieta )
            {
                eta_ring[layer][ieta] = AnaUtils::get_corrected_calo_eta( (AnaUtils::CaloType) layer, ieta, zvrtx );
            }
        }

        const float mbdQ = mbd_q_N + mbd_q_S;
        const int izbin = find_bin( zvrtx, calib_zvtx_edges );
        const int imbd  = find_bin( mbdQ, calib_mbdQ_edges );
        const bool calib_valid = ( izbin >= 0 && imbd >= 0 );

        float w_calib[3][24] = {};
        if ( calib_valid )
        {
            for ( int layer = 0; layer < 3; ++layer )
            {
                for ( int ieta = 0; ieta < 24; ++ieta )
                {
                    const int channel = encode_channel( ieta, izbin, imbd, calib_n_mbdbins );
                    w_calib[layer][ieta] = calib_tree -> GetFloatValue( channel, calib_field_names[layer] );
                }
            }
        }

        const float rho_layer_val[3] = { rho_val_cemc, rho_val_hcalin, rho_val_hcalout };

        // sum RAW rho_jet constituents within dR < R of (cone_eta, cone_phi),
        // decoding each one's calo layer from srcID and its eta ring from
        // eta_ring[][], and applying all three manual subtractions in the
        // same pass.
        auto compute_cone_sums = [&]( const float cone_eta, const float cone_phi, const float R ) -> ConeSums
        {
            static const float kSqrt2 = std::sqrt( 2.0f );
            ConeSums s;
            const size_t njets = rho_jet_constituent_E -> size();
            for ( size_t ij = 0; ij < njets; ++ij )
            {
                const auto & cE   = rho_jet_constituent_E     -> at(ij);
                const auto & cEta = rho_jet_constituent_eta   -> at(ij);
                const auto & cPhi = rho_jet_constituent_phi   -> at(ij);
                const auto & cSrc = rho_jet_constituent_srcID -> at(ij);
                const size_t ncon = cE.size();
                for ( size_t ic = 0; ic < ncon; ++ic )
                {
                    const float eta = cEta[ic];
                    if ( AnaUtils::deta_abs( eta, cone_eta ) > R ) continue; // cheap pre-cut
                    const float phi = cPhi[ic];
                    if ( AnaUtils::calc_dr( eta, phi, cone_eta, cone_phi ) >= R ) continue;

                    const int layer = decode_layer( cSrc[ic] );
                    if ( layer < 0 ) continue;

                    int ieta = -1;
                    float best_d = 1e9f;
                    for ( int j = 0; j < 24; ++j )
                    {
                        const float d = std::fabs( eta_ring[layer][j] - eta );
                        if ( d < best_d ) { best_d = d; ieta = j; }
                    }

                    const float E_raw = cE[ic];
                    const float inv_cosh = 1.0f / std::cosh( eta );
                    const float eT = E_raw * inv_cosh;
                    s.raw         += eT;
                    s.sub1        += ( E_raw - sub2_towerbkgd_ue -> at(layer).at(ieta) ) * inv_cosh;

                    const float rho_uncalib_val = rho_layer_val[layer];
                    s.rho_uncalib += ( E_raw - rho_uncalib_val * std::cosh( eta ) ) * inv_cosh;
                    // leave the tower unsubtracted if it looks like a jet
                    // fragment rather than background
                    s.rho_uncalib_seed += ( eT > kSqrt2 * rho_uncalib_val ) ? eT : ( E_raw - rho_uncalib_val * std::cosh( eta ) ) * inv_cosh;

                    if ( calib_valid )
                    {
                        const float rho_calib_val = rho_layer_val[layer] * w_calib[layer][ieta];
                        s.rho_calib += ( E_raw - rho_calib_val * std::cosh( eta ) ) * inv_cosh;
                        s.rho_calib_seed += ( eT > kSqrt2 * rho_calib_val ) ? eT : ( E_raw - rho_calib_val * std::cosh( eta ) ) * inv_cosh;
                    }
                }
            }
            return s;
        };

        //----------------------------------------------------------------
        // 1) random cones on the rho-subtracted towers
        //----------------------------------------------------------------
        const float cone_R = ( rho_jet_R > 0 ) ? rho_jet_R : ( ( sub1_jet_R > 0 ) ? sub1_jet_R : 0.4f );
        const int max_attempts = 50;
        for ( int icone = 0; icone < n_cones_per_event; ++icone )
        {
            float cone_eta = 0, cone_phi = 0;
            bool placed = false;
            for ( int iattempt = 0; iattempt < max_attempts; ++iattempt )
            {
                cone_eta = rnd.Uniform( -1.1, 1.1 );
                cone_phi = rnd.Uniform( 0.0, 2.0 * TMath::Pi() );
                if ( !AnaUtils::accept_jet_eta( cone_eta, zvrtx, cone_R ) ) continue;

                // no veto on proximity to a jet -- cones are allowed to land
                // on top of jets.
                placed = true;
                break;
            }
            if ( !placed ) continue;

            const float cone_pt = sum_cone_pt( cone_eta, cone_phi, cone_R );
            const float cone_pt_sub1 = sum_cone_pt_sub1( cone_eta, cone_phi, cone_R );
            const ConeSums cs = compute_cone_sums( cone_eta, cone_phi, cone_R );
            const float cone_pt_rhoconst = cs.raw;
            ++n_cones_thrown;
            h_rc_pt_all -> Fill( cone_pt );
            if ( icent >= 0 ) h_rc_pt_cent[icent] -> Fill( cone_pt );
            h_rc_pt_sub1_all -> Fill( cone_pt_sub1 );
            if ( icent >= 0 ) h_rc_pt_sub1_cent[icent] -> Fill( cone_pt_sub1 );
            h_rc_pt_rhoconst_all -> Fill( cone_pt_rhoconst );
            if ( icent >= 0 ) h_rc_pt_rhoconst_cent[icent] -> Fill( cone_pt_rhoconst );
            h_rc_pt_corr_tower_vs_sub1 -> Fill( cone_pt, cone_pt_sub1 );
            h_rc_dpt_sub1_minus_tower  -> Fill( cone_pt_sub1 - cone_pt );
            h_rc_pt_corr_tower_vs_rhoconst -> Fill( cone_pt, cone_pt_rhoconst );
            h_rc_dpt_rhoconst_minus_tower  -> Fill( cone_pt_rhoconst - cone_pt );

            h_rc_pt_manualsub1_all -> Fill( cs.sub1 );
            if ( icent >= 0 ) h_rc_pt_manualsub1_cent[icent] -> Fill( cs.sub1 );
            h_rc_pt_rhouncalib_all -> Fill( cs.rho_uncalib );
            if ( icent >= 0 ) h_rc_pt_rhouncalib_cent[icent] -> Fill( cs.rho_uncalib );
            h_rc_pt_rhouncalib_seed_all -> Fill( cs.rho_uncalib_seed );
            if ( icent >= 0 ) h_rc_pt_rhouncalib_seed_cent[icent] -> Fill( cs.rho_uncalib_seed );

            h_rc_pt_corr_manualsub1_vs_sub1const -> Fill( cs.sub1, cone_pt_sub1 );
            h_rc_dpt_manualsub1_minus_sub1const  -> Fill( cs.sub1 - cone_pt_sub1 );
            h_rc_dpt_rhouncalib_seed_minus_rhouncalib -> Fill( cs.rho_uncalib_seed - cs.rho_uncalib );

            if ( calib_valid )
            {
                h_rc_pt_rhocalib_all -> Fill( cs.rho_calib );
                if ( icent >= 0 ) h_rc_pt_rhocalib_cent[icent] -> Fill( cs.rho_calib );
                h_rc_pt_rhocalib_seed_all -> Fill( cs.rho_calib_seed );
                if ( icent >= 0 ) h_rc_pt_rhocalib_seed_cent[icent] -> Fill( cs.rho_calib_seed );

                h_rc_pt_corr_rhocalib_vs_tower -> Fill( cs.rho_calib, cone_pt );
                h_rc_dpt_rhocalib_minus_tower  -> Fill( cs.rho_calib - cone_pt );
                h_rc_pt_corr_rhouncalib_vs_rhocalib -> Fill( cs.rho_uncalib, cs.rho_calib );

                h_rc_dpt_rhocalib_seed_minus_rhocalib -> Fill( cs.rho_calib_seed - cs.rho_calib );
                h_rc_pt_corr_rhocalib_seed_vs_tower   -> Fill( cs.rho_calib_seed, cone_pt );
                h_rc_dpt_rhocalib_seed_minus_tower    -> Fill( cs.rho_calib_seed - cone_pt );
            }
        }

        //----------------------------------------------------------------
        // 2) rho vs sub1 jet kinematics
        //----------------------------------------------------------------
        auto rho_sel = AnaUtils::select_jets(
            *rho_jet_pT, *rho_jet_E, *rho_jet_eta,
            min_jet_pt, zvrtx, rho_jet_R, true
        );
        auto sub1_sel = AnaUtils::select_jets(
            *sub1_jet_pT, *sub1_jet_E, *sub1_jet_eta,
            min_jet_pt, zvrtx, sub1_jet_R, true
        );

        h_njets_rho  -> Fill( (double) rho_sel.size() );
        h_njets_sub1 -> Fill( (double) sub1_sel.size() );

        for ( const int i : rho_sel )
        {
            h_rho_jet_pt     -> Fill( rho_jet_pT->at(i) );
            h_rho_jet_eta    -> Fill( rho_jet_eta->at(i) );
            h_rho_jet_phi    -> Fill( rho_jet_phi->at(i) );
            h_rho_jet_pt_eta -> Fill( rho_jet_eta->at(i), rho_jet_pT->at(i) );
        }
        for ( const int i : sub1_sel )
        {
            h_sub1_jet_pt     -> Fill( sub1_jet_pT->at(i) );
            h_sub1_jet_eta    -> Fill( sub1_jet_eta->at(i) );
            h_sub1_jet_phi    -> Fill( sub1_jet_phi->at(i) );
            h_sub1_jet_pt_eta -> Fill( sub1_jet_eta->at(i), sub1_jet_pT->at(i) );
        }

        //----------------------------------------------------------------
        // 3) match rho <-> sub1 jets (same underlying jet, two subtractions)
        //----------------------------------------------------------------
        const float max_dr = 0.75f * std::max( rho_jet_R, sub1_jet_R );
        auto matches = AnaUtils::match_truth_reco_jets(
            rho_sel, *rho_jet_eta, *rho_jet_phi,
            sub1_sel, *sub1_jet_eta, *sub1_jet_phi,
            max_dr
        );

        for ( const auto & jm : matches )
        {
            const bool has_rho  = jm.truth_index >= 0;
            const bool has_sub1 = jm.reco_index >= 0;
            if ( has_rho && has_sub1 )
            {
                ++n_matched;
                const float rho_pt  = rho_jet_pT->at( jm.truth_index );
                const float sub1_pt = sub1_jet_pT->at( jm.reco_index );
                h_match_pt_rho_vs_sub1 -> Fill( rho_pt, sub1_pt );
                h_match_dpt            -> Fill( sub1_pt - rho_pt );
                if ( rho_pt > 0 ) h_match_dpt_relpt -> Fill( ( sub1_pt - rho_pt ) / rho_pt );
                h_match_dr    -> Fill( jm.dr );
                h_match_deta  -> Fill( sub1_jet_eta->at( jm.reco_index ) - rho_jet_eta->at( jm.truth_index ) );
                h_match_dphi  -> Fill( AnaUtils::dphi_wrap( rho_jet_phi->at( jm.truth_index ), sub1_jet_phi->at( jm.reco_index ) ) );
            }
            else if ( has_rho )  { ++n_rho_only; }
            else                 { ++n_sub1_only; }
        }
    }

    std::cout << "Random cones: thrown " << n_cones_thrown << std::endl;
    std::cout << "rho<->sub1 matching: " << n_matched << " matched, "
              << n_rho_only << " rho-only, " << n_sub1_only << " sub1-only" << std::endl;

    //--------------------------------------------------------------------
    // mean/sigma of the random cone pT distribution vs centrality
    //--------------------------------------------------------------------
    auto * g_rc_mean_cent  = new TGraphErrors( NCENT );
    g_rc_mean_cent  -> SetName( "g_rc_mean_cent" );
    g_rc_mean_cent  -> SetTitle( "Random cone p_{T} mean vs centrality;centrality [%];<p_{T}^{cone}> [GeV]" );
    auto * g_rc_sigma_cent = new TGraphErrors( NCENT );
    g_rc_sigma_cent -> SetName( "g_rc_sigma_cent" );
    g_rc_sigma_cent -> SetTitle( "Random cone p_{T} width vs centrality;centrality [%];#sigma(p_{T}^{cone}) [GeV]" );

    auto * g_rc_sub1_mean_cent  = new TGraphErrors( NCENT );
    g_rc_sub1_mean_cent  -> SetName( "g_rc_sub1_mean_cent" );
    g_rc_sub1_mean_cent  -> SetTitle( "Random cone p_{T} mean vs centrality (sub1 constituents);centrality [%];<p_{T}^{cone}> [GeV]" );
    auto * g_rc_sub1_sigma_cent = new TGraphErrors( NCENT );
    g_rc_sub1_sigma_cent -> SetName( "g_rc_sub1_sigma_cent" );
    g_rc_sub1_sigma_cent -> SetTitle( "Random cone p_{T} width vs centrality (sub1 constituents);centrality [%];#sigma(p_{T}^{cone}) [GeV]" );

    auto * g_rc_rhoconst_mean_cent  = new TGraphErrors( NCENT );
    g_rc_rhoconst_mean_cent  -> SetName( "g_rc_rhoconst_mean_cent" );
    g_rc_rhoconst_mean_cent  -> SetTitle( "Random cone p_{T} mean vs centrality, RAW/unsubtracted (rho_jet constituents);centrality [%];<p_{T}^{cone}> [GeV]" );
    auto * g_rc_rhoconst_sigma_cent = new TGraphErrors( NCENT );
    g_rc_rhoconst_sigma_cent -> SetName( "g_rc_rhoconst_sigma_cent" );
    g_rc_rhoconst_sigma_cent -> SetTitle( "Random cone p_{T} width vs centrality, RAW/unsubtracted (rho_jet constituents);centrality [%];#sigma(p_{T}^{cone}) [GeV]" );

    auto * g_rc_manualsub1_mean_cent  = new TGraphErrors( NCENT );
    g_rc_manualsub1_mean_cent  -> SetName( "g_rc_manualsub1_mean_cent" );
    g_rc_manualsub1_mean_cent  -> SetTitle( "Random cone p_{T} mean vs centrality, manual sub1;centrality [%];<p_{T}^{cone}> [GeV]" );
    auto * g_rc_manualsub1_sigma_cent = new TGraphErrors( NCENT );
    g_rc_manualsub1_sigma_cent -> SetName( "g_rc_manualsub1_sigma_cent" );
    g_rc_manualsub1_sigma_cent -> SetTitle( "Random cone p_{T} width vs centrality, manual sub1;centrality [%];#sigma(p_{T}^{cone}) [GeV]" );

    auto * g_rc_rhouncalib_mean_cent  = new TGraphErrors( NCENT );
    g_rc_rhouncalib_mean_cent  -> SetName( "g_rc_rhouncalib_mean_cent" );
    g_rc_rhouncalib_mean_cent  -> SetTitle( "Random cone p_{T} mean vs centrality, uncalibrated rho;centrality [%];<p_{T}^{cone}> [GeV]" );
    auto * g_rc_rhouncalib_sigma_cent = new TGraphErrors( NCENT );
    g_rc_rhouncalib_sigma_cent -> SetName( "g_rc_rhouncalib_sigma_cent" );
    g_rc_rhouncalib_sigma_cent -> SetTitle( "Random cone p_{T} width vs centrality, uncalibrated rho;centrality [%];#sigma(p_{T}^{cone}) [GeV]" );

    auto * g_rc_rhocalib_mean_cent  = new TGraphErrors( NCENT );
    g_rc_rhocalib_mean_cent  -> SetName( "g_rc_rhocalib_mean_cent" );
    g_rc_rhocalib_mean_cent  -> SetTitle( "Random cone p_{T} mean vs centrality, calibrated rho;centrality [%];<p_{T}^{cone}> [GeV]" );
    auto * g_rc_rhocalib_sigma_cent = new TGraphErrors( NCENT );
    g_rc_rhocalib_sigma_cent -> SetName( "g_rc_rhocalib_sigma_cent" );
    g_rc_rhocalib_sigma_cent -> SetTitle( "Random cone p_{T} width vs centrality, calibrated rho;centrality [%];#sigma(p_{T}^{cone}) [GeV]" );

    auto * g_rc_rhouncalib_seed_mean_cent  = new TGraphErrors( NCENT );
    g_rc_rhouncalib_seed_mean_cent  -> SetName( "g_rc_rhouncalib_seed_mean_cent" );
    g_rc_rhouncalib_seed_mean_cent  -> SetTitle( "Random cone p_{T} mean vs centrality, uncalibrated rho, seed-protected;centrality [%];<p_{T}^{cone}> [GeV]" );
    auto * g_rc_rhouncalib_seed_sigma_cent = new TGraphErrors( NCENT );
    g_rc_rhouncalib_seed_sigma_cent -> SetName( "g_rc_rhouncalib_seed_sigma_cent" );
    g_rc_rhouncalib_seed_sigma_cent -> SetTitle( "Random cone p_{T} width vs centrality, uncalibrated rho, seed-protected;centrality [%];#sigma(p_{T}^{cone}) [GeV]" );

    auto * g_rc_rhocalib_seed_mean_cent  = new TGraphErrors( NCENT );
    g_rc_rhocalib_seed_mean_cent  -> SetName( "g_rc_rhocalib_seed_mean_cent" );
    g_rc_rhocalib_seed_mean_cent  -> SetTitle( "Random cone p_{T} mean vs centrality, calibrated rho, seed-protected;centrality [%];<p_{T}^{cone}> [GeV]" );
    auto * g_rc_rhocalib_seed_sigma_cent = new TGraphErrors( NCENT );
    g_rc_rhocalib_seed_sigma_cent -> SetName( "g_rc_rhocalib_seed_sigma_cent" );
    g_rc_rhocalib_seed_sigma_cent -> SetTitle( "Random cone p_{T} width vs centrality, calibrated rho, seed-protected;centrality [%];#sigma(p_{T}^{cone}) [GeV]" );

    for ( int i = 0; i < NCENT; ++i )
    {
        const float bin_center = 0.5f * ( cent_edges[i] + cent_edges[i + 1] );
        const float bin_half_width = 0.5f * ( cent_edges[i + 1] - cent_edges[i] );
        g_rc_mean_cent  -> SetPoint( i, bin_center, h_rc_pt_cent[i] -> GetMean() );
        g_rc_mean_cent  -> SetPointError( i, bin_half_width, h_rc_pt_cent[i] -> GetMeanError() );
        g_rc_sigma_cent -> SetPoint( i, bin_center, h_rc_pt_cent[i] -> GetRMS() );
        g_rc_sigma_cent -> SetPointError( i, bin_half_width, h_rc_pt_cent[i] -> GetRMSError() );

        g_rc_sub1_mean_cent  -> SetPoint( i, bin_center, h_rc_pt_sub1_cent[i] -> GetMean() );
        g_rc_sub1_mean_cent  -> SetPointError( i, bin_half_width, h_rc_pt_sub1_cent[i] -> GetMeanError() );
        g_rc_sub1_sigma_cent -> SetPoint( i, bin_center, h_rc_pt_sub1_cent[i] -> GetRMS() );
        g_rc_sub1_sigma_cent -> SetPointError( i, bin_half_width, h_rc_pt_sub1_cent[i] -> GetRMSError() );

        g_rc_rhoconst_mean_cent  -> SetPoint( i, bin_center, h_rc_pt_rhoconst_cent[i] -> GetMean() );
        g_rc_rhoconst_mean_cent  -> SetPointError( i, bin_half_width, h_rc_pt_rhoconst_cent[i] -> GetMeanError() );
        g_rc_rhoconst_sigma_cent -> SetPoint( i, bin_center, h_rc_pt_rhoconst_cent[i] -> GetRMS() );
        g_rc_rhoconst_sigma_cent -> SetPointError( i, bin_half_width, h_rc_pt_rhoconst_cent[i] -> GetRMSError() );

        g_rc_manualsub1_mean_cent  -> SetPoint( i, bin_center, h_rc_pt_manualsub1_cent[i] -> GetMean() );
        g_rc_manualsub1_mean_cent  -> SetPointError( i, bin_half_width, h_rc_pt_manualsub1_cent[i] -> GetMeanError() );
        g_rc_manualsub1_sigma_cent -> SetPoint( i, bin_center, h_rc_pt_manualsub1_cent[i] -> GetRMS() );
        g_rc_manualsub1_sigma_cent -> SetPointError( i, bin_half_width, h_rc_pt_manualsub1_cent[i] -> GetRMSError() );

        g_rc_rhouncalib_mean_cent  -> SetPoint( i, bin_center, h_rc_pt_rhouncalib_cent[i] -> GetMean() );
        g_rc_rhouncalib_mean_cent  -> SetPointError( i, bin_half_width, h_rc_pt_rhouncalib_cent[i] -> GetMeanError() );
        g_rc_rhouncalib_sigma_cent -> SetPoint( i, bin_center, h_rc_pt_rhouncalib_cent[i] -> GetRMS() );
        g_rc_rhouncalib_sigma_cent -> SetPointError( i, bin_half_width, h_rc_pt_rhouncalib_cent[i] -> GetRMSError() );

        g_rc_rhocalib_mean_cent  -> SetPoint( i, bin_center, h_rc_pt_rhocalib_cent[i] -> GetMean() );
        g_rc_rhocalib_mean_cent  -> SetPointError( i, bin_half_width, h_rc_pt_rhocalib_cent[i] -> GetMeanError() );
        g_rc_rhocalib_sigma_cent -> SetPoint( i, bin_center, h_rc_pt_rhocalib_cent[i] -> GetRMS() );
        g_rc_rhocalib_sigma_cent -> SetPointError( i, bin_half_width, h_rc_pt_rhocalib_cent[i] -> GetRMSError() );

        g_rc_rhouncalib_seed_mean_cent  -> SetPoint( i, bin_center, h_rc_pt_rhouncalib_seed_cent[i] -> GetMean() );
        g_rc_rhouncalib_seed_mean_cent  -> SetPointError( i, bin_half_width, h_rc_pt_rhouncalib_seed_cent[i] -> GetMeanError() );
        g_rc_rhouncalib_seed_sigma_cent -> SetPoint( i, bin_center, h_rc_pt_rhouncalib_seed_cent[i] -> GetRMS() );
        g_rc_rhouncalib_seed_sigma_cent -> SetPointError( i, bin_half_width, h_rc_pt_rhouncalib_seed_cent[i] -> GetRMSError() );

        g_rc_rhocalib_seed_mean_cent  -> SetPoint( i, bin_center, h_rc_pt_rhocalib_seed_cent[i] -> GetMean() );
        g_rc_rhocalib_seed_mean_cent  -> SetPointError( i, bin_half_width, h_rc_pt_rhocalib_seed_cent[i] -> GetMeanError() );
        g_rc_rhocalib_seed_sigma_cent -> SetPoint( i, bin_center, h_rc_pt_rhocalib_seed_cent[i] -> GetRMS() );
        g_rc_rhocalib_seed_sigma_cent -> SetPointError( i, bin_half_width, h_rc_pt_rhocalib_seed_cent[i] -> GetRMSError() );
    }

    //--------------------------------------------------------------------
    // write output
    //--------------------------------------------------------------------
    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    fout -> cd();

    h_rc_pt_all -> Write();
    for ( int i = 0; i < NCENT; ++i ) h_rc_pt_cent[i] -> Write();
    g_rc_mean_cent  -> Write();
    g_rc_sigma_cent -> Write();

    h_rc_pt_sub1_all -> Write();
    for ( int i = 0; i < NCENT; ++i ) h_rc_pt_sub1_cent[i] -> Write();
    g_rc_sub1_mean_cent  -> Write();
    g_rc_sub1_sigma_cent -> Write();

    h_rc_pt_corr_tower_vs_sub1 -> Write();
    h_rc_dpt_sub1_minus_tower  -> Write();

    h_rc_pt_rhoconst_all -> Write();
    for ( int i = 0; i < NCENT; ++i ) h_rc_pt_rhoconst_cent[i] -> Write();
    g_rc_rhoconst_mean_cent  -> Write();
    g_rc_rhoconst_sigma_cent -> Write();

    h_rc_pt_corr_tower_vs_rhoconst -> Write();
    h_rc_dpt_rhoconst_minus_tower  -> Write();

    h_rc_pt_manualsub1_all -> Write();
    for ( int i = 0; i < NCENT; ++i ) h_rc_pt_manualsub1_cent[i] -> Write();
    g_rc_manualsub1_mean_cent  -> Write();
    g_rc_manualsub1_sigma_cent -> Write();

    h_rc_pt_rhouncalib_all -> Write();
    for ( int i = 0; i < NCENT; ++i ) h_rc_pt_rhouncalib_cent[i] -> Write();
    g_rc_rhouncalib_mean_cent  -> Write();
    g_rc_rhouncalib_sigma_cent -> Write();

    h_rc_pt_rhocalib_all -> Write();
    for ( int i = 0; i < NCENT; ++i ) h_rc_pt_rhocalib_cent[i] -> Write();
    g_rc_rhocalib_mean_cent  -> Write();
    g_rc_rhocalib_sigma_cent -> Write();

    h_rc_pt_corr_manualsub1_vs_sub1const -> Write();
    h_rc_dpt_manualsub1_minus_sub1const  -> Write();

    h_rc_pt_corr_rhocalib_vs_tower -> Write();
    h_rc_dpt_rhocalib_minus_tower  -> Write();

    h_rc_pt_corr_rhouncalib_vs_rhocalib -> Write();

    h_rc_pt_rhouncalib_seed_all -> Write();
    for ( int i = 0; i < NCENT; ++i ) h_rc_pt_rhouncalib_seed_cent[i] -> Write();
    g_rc_rhouncalib_seed_mean_cent  -> Write();
    g_rc_rhouncalib_seed_sigma_cent -> Write();

    h_rc_pt_rhocalib_seed_all -> Write();
    for ( int i = 0; i < NCENT; ++i ) h_rc_pt_rhocalib_seed_cent[i] -> Write();
    g_rc_rhocalib_seed_mean_cent  -> Write();
    g_rc_rhocalib_seed_sigma_cent -> Write();

    h_rc_dpt_rhouncalib_seed_minus_rhouncalib -> Write();
    h_rc_dpt_rhocalib_seed_minus_rhocalib     -> Write();
    h_rc_pt_corr_rhocalib_seed_vs_tower       -> Write();
    h_rc_dpt_rhocalib_seed_minus_tower        -> Write();

    h_rho_jet_pt     -> Write();
    h_sub1_jet_pt    -> Write();
    h_rho_jet_eta    -> Write();
    h_sub1_jet_eta   -> Write();
    h_rho_jet_phi    -> Write();
    h_sub1_jet_phi   -> Write();
    h_rho_jet_pt_eta  -> Write();
    h_sub1_jet_pt_eta -> Write();
    h_njets_rho  -> Write();
    h_njets_sub1 -> Write();

    h_match_pt_rho_vs_sub1 -> Write();
    h_match_dpt       -> Write();
    h_match_dpt_relpt -> Write();
    h_match_dr        -> Write();
    h_match_deta      -> Write();
    h_match_dphi      -> Write();

    fout -> Close();

    std::cout << "Wrote " << outfile << std::endl;

    return 0;
}


#endif // _PARSE_TREE_C_
