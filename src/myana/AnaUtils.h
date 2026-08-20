#ifndef _ANAUTILS_H_
#define _ANAUTILS_H_

#include <vector>
#include <string>
#include <array>

class TTree;

namespace AnaUtils
{
    float get_dpsi2( const float psi2, const float phi );
    
    float deta_abs( const float eta1, const float eta2 ) ;
    float dphi_wrap( const float phi1, const float phi2 ) ;

    float calc_dr( const float eta1, const float phi1, const float eta2, const float phi2 ) ;

    bool accept_jet_eta( const float eta, const float zvrtx, const float jet_R  = 0.3 ) ;
    float correct_calo_eta( const float eta0, const float zvrtx, const float R ) ;

    enum CaloType { CEMC = 0, HCALIN = 1, HCALOUT = 2 };
    inline static constexpr std::array<float, 24> eta = {
        -1.05416667f, -0.96249998f, -0.87083334f, -0.77916664f, 
        -0.68750000f, -0.59583336f, -0.50416666f, -0.41249999f, 
        -0.32083333f, -0.22916667f, -0.13750000f, -0.04583333f, 
        0.04583333f, 0.13750000f, 0.22916667f, 0.32083333f, 
        0.41249999f, 0.50416666f, 0.59583336f, 0.68750000f, 
        0.77916664f, 0.87083334f, 0.96249998f, 1.05416667f
    };
    inline static constexpr std::array<float,64> hcalin_phi  = {
        6.27865314f, 0.09364238f, 0.19181715f, 0.28999192f, 0.38816670f, 
        0.48634145f, 0.58451623f, 0.68269098f, 0.78086579f, 0.87904054f, 
        0.97721529f, 1.07539010f, 1.17356479f, 1.27173960f, 1.36991441f, 
        1.46808910f, 1.56626391f, 1.66443872f, 1.76261342f, 1.86078823f, 
        1.95896304f, 2.05713773f, 2.15531254f, 2.25348735f, 2.35166216f, 
        2.44983697f, 2.54801154f, 2.64618635f, 2.74436116f, 2.84253597f, 
        2.94071078f, 3.03888559f, 3.13706017f, 3.23523498f, 3.33340979f, 
        3.43158460f, 3.52975941f, 3.62793422f, 3.72610879f, 3.82428360f, 
        3.92245841f, 4.02063322f, 4.11880779f, 4.21698284f, 4.31515741f, 
        4.41333246f, 4.51150703f, 4.60968161f, 4.70785666f, 4.80603123f, 
        4.90420628f, 5.00238085f, 5.10055590f, 5.19873047f, 5.29690504f, 
        5.39508009f, 5.49325466f, 5.59142971f, 5.68960428f, 5.78777885f, 
        5.88595390f, 5.98412848f, 6.08230352f, 6.18047810f
    };
    inline static constexpr std::array<float,64> hcalout_phi = {
        6.25815964f, 0.07314893f, 0.17132370f, 0.26949847f, 0.36767325f, 
        0.46584800f, 0.56402278f, 0.66219753f, 0.76037234f, 0.85854709f, 
        0.95672184f, 1.05489659f, 1.15307140f, 1.25124621f, 1.34942091f, 
        1.44759572f, 1.54577053f, 1.64394522f, 1.74212003f, 1.84029484f, 
        1.93846953f, 2.03664422f, 2.13481903f, 2.23299384f, 2.33116865f, 
        2.42934346f, 2.52751827f, 2.62569284f, 2.72386765f, 2.82204247f, 
        2.92021728f, 3.01839209f, 3.11656690f, 3.21474147f, 3.31291628f, 
        3.41109109f, 3.50926590f, 3.60744071f, 3.70561552f, 3.80379009f, 
        3.90196490f, 4.00013971f, 4.09831429f, 4.19648933f, 4.29466391f, 
        4.39283895f, 4.49101353f, 4.58918858f, 4.68736315f, 4.78553772f, 
        4.88371277f, 4.98188734f, 5.08006239f, 5.17823696f, 5.27641153f, 
        5.37458658f, 5.47276115f, 5.57093620f, 5.66911077f, 5.76728582f, 
        5.86546040f, 5.96363497f, 6.06181002f, 6.15998459f
    };

    float get_calo_eta( const int ieta );
    float get_calo_phi( const CaloType calo, const int iphi );
    float get_calo_r( const CaloType calo );
    float get_corrected_calo_eta( const CaloType calo, const int ieta, const float zvrtx );

    // sum of transverse energy (E/cosh(eta), using the z-vertex corrected
    // tower eta) over towers flagged good in tower_isgood[ieta][iphi]
    float calc_sumeT(
        const CaloType calo,
        const float zvrtx,
        const float tower_E[24][64],
        const int tower_isgood[24][64]
    );

    void myText( double x, double y, int color, const char * text, const float size = 0.03 );

    double flow_func( double * x, double * par );

    bool phi_top ( const float phi );
    bool phi_bottom ( const float phi );
    bool phi_west ( const float phi );
    bool phi_east ( const float phi );
    bool phi_vert ( const float phi );
    bool phi_horz ( const float phi );
    bool eta_neg ( const float eta );
    bool eta_pos ( const float eta );

    bool in_plane( const float psi2, const float phi_jet ) ;
    bool mid_plane( const float psi2, const float phi_jet ) ;
    bool out_of_plane( const float psi2, const float phi_jet ) ;   


    

   std::vector< std::string > getFilelist( const std::string & inlist , const std::string & ext = ".root" );

   //--------------------------------------------------------------------
   // truth-reco jet matching
   //--------------------------------------------------------------------

   // indices (into the original, unfiltered jet vectors) of jets passing
   // basic kinematic cuts, sorted by descending pT. the eta cut uses
   // accept_jet_eta(), i.e. the z-vertex- and jet-radius-dependent
   // calorimeter acceptance window, not a flat |eta| cut.
   std::vector<int> select_jets(
       const std::vector<float> & pt,
       const std::vector<float> & e,
       const std::vector<float> & eta,
       const float min_pt,
       const float zvrtx,
       const float jet_R,
       const bool require_e_positive = false
   );

   // one row per truth/reco jet relationship: a matched pair has both
   // indices set, an unmatched truth jet has reco_index == -1, and an
   // unmatched reco jet has truth_index == -1. Indices refer back into
   // the original (unfiltered) truth/reco jet vectors passed in.
   struct JetMatch
   {
       int truth_index { -1 };
       int reco_index  { -1 };
       float dr        { -1.0f };
   };

   // greedy nearest-neighbor matching (globally ascending dR) between the
   // truth jets listed in truth_indices and the reco jets listed in
   // reco_indices. Returns matched pairs first, then unmatched truth jets,
   // then unmatched reco jets -- ready to be written out in that order.
   std::vector<JetMatch> match_truth_reco_jets(
       const std::vector<int> & truth_indices,
       const std::vector<float> & truth_eta,
       const std::vector<float> & truth_phi,
       const std::vector<int> & reco_indices,
       const std::vector<float> & reco_eta,
       const std::vector<float> & reco_phi,
       const float max_dr
   );

   // flat, one-row-per-jet output format for a matched truth/reco jet tree.
   // match_status: 0 = matched pair, 1 = unmatched truth, 2 = unmatched reco.
   // reco_type: 0 = rho-subtracted jet, 1 = sub1(seeded)-subtracted jet.
   // Unfilled fields (e.g. truth_* for an unmatched reco row) are set to -999.
   struct MatchedJetRow
   {
       // event level
       int   event_id    { -999 };
       int   cent        { -999 };
       float zvrtx       { -999.0f };
       float mbdQ        { -999.0f };
       float sumeT       { -999.0f };
       int   is_minbias  { -999 };
       float psi2        { -999.0f };
       float truth_jet_maxpt_r04 { -999.0f };

       // matching info
       int   reco_type    { -999 };
       int   match_status { -999 };
       float dr           { -999.0f };

       // truth jet
       float truth_pt     { -999.0f };
       float truth_e      { -999.0f };
       float truth_eta    { -999.0f };
       float truth_phi    { -999.0f };
       int   truth_flavor { -999 };

       // matched reco jet
       float reco_pt        { -999.0f };
       float reco_e         { -999.0f };
       float reco_eta       { -999.0f };
       float reco_phi       { -999.0f };
       float reco_unsub_e   { -999.0f };
       float reco_unsub_pt  { -999.0f };
   };

   // creates one branch per MatchedJetRow member on tree, bound to row.
   // call once, right after constructing the (empty) output TTree.
   void book_matched_jet_tree( TTree * tree, MatchedJetRow & row );

   // binds tree's branches to row for reading, e.g. when analyzing the
   // output of a macro built on book_matched_jet_tree(). call once, right
   // after opening the tree, then GetEntry() to fill row.
   void read_matched_jet_tree( TTree * tree, MatchedJetRow & row );

} // namespace MyAna

#endif // _ANAUTILS_H_