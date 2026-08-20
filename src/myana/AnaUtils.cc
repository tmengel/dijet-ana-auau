#include "AnaUtils.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <TMath.h>
#include <TLatex.h>
#include <TTree.h>

std::vector< std::string > AnaUtils::getFilelist( const std::string & inlist , const std::string & ext )
{
    std::vector<std::string> files;
    std::ifstream infile( inlist );
    if ( !infile.is_open() )
    {
        std::cerr << "Error: could not open input file list " << inlist << std::endl;
        return files;
    }

    std::string line;
    while ( std::getline( infile, line ) )
    {
        if ( line.empty() ) 
        {
            continue;
        }
        if ( ext != "" && line.find( ext ) == std::string::npos )
        {
             continue;
        }
        files.push_back( line );
    }
    infile.close();
    return files;
}

float AnaUtils::deta_abs(const float eta1, const float eta2)
{
    return fabs(eta2 - eta1);
}

float AnaUtils::dphi_wrap(const float phi1, const float phi2)
{
    float dphi = phi2-phi1;
    if (dphi > TMath::Pi())
    {
        dphi -= 2 * TMath::Pi();
    }
    if (dphi < -TMath::Pi())
    {
        dphi += 2 * TMath::Pi();
    }
    return fabs(dphi);
}

float AnaUtils::calc_dr(const float eta1, const float phi1, const float eta2, const float phi2)
{
    float dphi = dphi_wrap(phi1, phi2);
    float deta = deta_abs(eta1, eta2);
    return TMath::Sqrt(dphi*dphi + deta*deta);
}

float AnaUtils::get_dpsi2( const float psi2, const float phi_jet )
{
    float psi2_mod = psi2 - TMath::Pi();
    if ( psi2 < 0 )
    {
        psi2_mod = psi2 + TMath::Pi();
    }
    float dA = dphi_wrap( phi_jet, psi2 );
    float dB = dphi_wrap( phi_jet, psi2_mod );
    float dmin = std::min( dA, dB );
    return fabs(dmin);
}

float AnaUtils::correct_calo_eta( const float eta0, const float zvrtx, const float R )
{
    double z0 = sinh(eta0) * R;
    double z1 = z0 - zvrtx;
    double eta1 = asinh( z1 / R );
    return eta1;
}

bool AnaUtils::accept_jet_eta( const float eta, const float zvrtx, const float jet_R )
{
    const float CALO_ABS_Z[3] = {130.23, 170.299, 301.683};
    const float CALO_RADIUS[3] = {93.5, 127.503, 225.87};
    float min_eta = -999, MAX_eta = 999;
    for (int i=0; i<3; i++) 
    {
       float min_z = -1.0*CALO_ABS_Z[i];
       float max_z = CALO_ABS_Z[i];
       float z_l = min_z - zvrtx;
       float z_h = max_z - zvrtx;
       float r = CALO_RADIUS[i];
       float eta_min = asinh( z_l / r );
       float eta_max = asinh( z_h / r );
       if ( eta_min > min_eta ) { min_eta = eta_min; }
       if ( eta_max < MAX_eta ) { MAX_eta = eta_max; }
    }
    min_eta += jet_R;
    MAX_eta -= jet_R;
    return (eta >= min_eta && eta <= MAX_eta);
}

float AnaUtils::get_calo_eta( const int ieta )
{
    if ( ieta < 0 || ieta >= 24 )
    {
        std::cerr << "Error: invalid ieta " << ieta << ", must be in [0, 23]" << std::endl;
        return 0.0;
    }
    return eta[ieta];
}
float AnaUtils::get_calo_phi( const CaloType calo, const int iphi )
{
    if ( iphi < 0 || iphi >= 64 )
    {
        std::cerr << "Error: invalid iphi " << iphi << ", must be in [0, 63]" << std::endl;
        return 0.0;
    }
    if ( calo == CaloType::HCALOUT )
    {
        return hcalout_phi[iphi];
    }

    return hcalin_phi[iphi];
}
float AnaUtils::get_calo_r( const CaloType calo )
{
    if ( calo == CaloType::CEMC )
    {
        return 93.5;
    }
    else if ( calo == CaloType::HCALIN )
    {
        return 127.503;
    }
    else if ( calo == CaloType::HCALOUT )
    {
        return 225.87;
    }
    else
    {
        std::cerr << "Error: invalid calo type " << static_cast<int>(calo) << std::endl;
        return 0.0;
    }
}
float AnaUtils::get_corrected_calo_eta( const CaloType calo, const int ieta, const float zvrtx )
{
    float R = get_calo_r( calo );
    float eta0 = get_calo_eta( ieta );
    return correct_calo_eta( eta0, zvrtx, R );
}

float AnaUtils::calc_sumeT(
    const CaloType calo,
    const float zvrtx,
    const float tower_E[24][64],
    const int tower_isgood[24][64]
)
{
    float sum = 0.0;
    for ( int ieta = 0; ieta < 24; ++ieta )
    {
        const float eta_corr = get_corrected_calo_eta( calo, ieta, zvrtx );
        const float inv_cosh_eta = 1.0f / std::cosh( eta_corr );
        for ( int iphi = 0; iphi < 64; ++iphi )
        {
            if ( !tower_isgood[ieta][iphi] ) continue;
            sum += tower_E[ieta][iphi] * inv_cosh_eta;
        }
    }
    return sum;
}

double AnaUtils::flow_func( double * x, double * par )
{
   double a = par[0];
   double bsum = 1.0;
    for ( int i = 1; i < 5; ++i )
    {
        bsum += 2.0* par[i] * cos( static_cast<float>(i) * x[0] );
    }
    return a * bsum;
}

bool AnaUtils::phi_top ( const float phi )
{
    return ( phi >= TMath::Pi()/4 && phi < 3*TMath::Pi()/4 );
}

bool AnaUtils::phi_bottom ( const float phi )
{
    return ( phi >= -3*TMath::Pi()/4 && phi < -TMath::Pi()/4 );
}

bool AnaUtils::phi_west ( const float phi )
{
    return ( std::fabs(phi) < TMath::Pi()/4 );
}

bool AnaUtils::phi_east ( const float phi )
{
    return ( std::fabs(phi) >= 3*TMath::Pi()/4 );
}

bool AnaUtils::phi_vert ( const float phi )
{
    return ( (phi >= TMath::Pi()/4 && phi < 3*TMath::Pi()/4) || (phi >= -3*TMath::Pi()/4 && phi < -TMath::Pi()/4) );
}

bool AnaUtils::phi_horz ( const float phi )
{
    return ( std::fabs(phi) < 3*TMath::Pi()/4 && std::fabs(phi) >= TMath::Pi()/4 );
}

bool AnaUtils::eta_neg ( const float eta )
{
    return ( eta < 0 );
}

bool AnaUtils::eta_pos ( const float eta )
{
    return ( eta >= 0 );
}

bool AnaUtils::in_plane( const float psi2, const float phi_jet ) 
{
    float dpsi2 = get_dpsi2( psi2, phi_jet );
    return ( dpsi2 < TMath::Pi()/6.0f ); // in-plane if within 30 degrees of event plane
}

bool AnaUtils::mid_plane( const float psi2, const float phi_jet ) 
{
    float dpsi2 = get_dpsi2( psi2, phi_jet );
    return ( dpsi2 >= TMath::Pi()/6.0f && dpsi2 < TMath::Pi()/3.0f ); // mid-plane if between 30 and 60 degrees from event plane
}

bool AnaUtils::out_of_plane( const float psi2, const float phi_jet ) 
{
    float dpsi2 = get_dpsi2( psi2, phi_jet );
    return ( dpsi2 >= TMath::Pi()/3.0f ); // out-of-plane if greater than 60 degrees from event plane
}  

void AnaUtils::myText( double x, double y, int color, const char * text, const float size )
{
    TLatex * t = new TLatex();
    t -> SetNDC();
    t -> SetTextSize(size);
    t -> SetTextColor(color);
    t -> DrawLatex(x, y, text);
}

std::vector<int> AnaUtils::select_jets(
    const std::vector<float> & pt,
    const std::vector<float> & e,
    const std::vector<float> & eta,
    const float min_pt,
    const float zvrtx,
    const float jet_R,
    const bool require_e_positive
)
{
    std::vector<int> selected;
    for ( size_t i = 0; i < pt.size(); ++i )
    {
        if ( pt[i] < min_pt ) continue;
        if ( !accept_jet_eta( eta[i], zvrtx, jet_R ) ) continue;
        if ( require_e_positive && e[i] < 0.0 ) continue;
        selected.push_back( static_cast<int>(i) );
    }
    std::sort( selected.begin(), selected.end(), [&]( int a, int b ) { return pt[a] > pt[b]; } );
    return selected;
}

std::vector<AnaUtils::JetMatch> AnaUtils::match_truth_reco_jets(
    const std::vector<int> & truth_indices,
    const std::vector<float> & truth_eta,
    const std::vector<float> & truth_phi,
    const std::vector<int> & reco_indices,
    const std::vector<float> & reco_eta,
    const std::vector<float> & reco_phi,
    const float max_dr
)
{
    struct Candidate { int truth_index; int reco_index; float dr; };
    std::vector<Candidate> candidates;
    candidates.reserve( truth_indices.size() * reco_indices.size() );
    for ( const int ti : truth_indices )
    {
        for ( const int ri : reco_indices )
        {
            const float dr = calc_dr( truth_eta[ti], truth_phi[ti], reco_eta[ri], reco_phi[ri] );
            if ( dr < max_dr )
            {
                candidates.push_back( { ti, ri, dr } );
            }
        }
    }
    std::sort( candidates.begin(), candidates.end(), [] ( const Candidate & a, const Candidate & b ) { return a.dr < b.dr; } );

    std::vector<bool> truth_used( truth_indices.empty() ? 0 : *std::max_element( truth_indices.begin(), truth_indices.end() ) + 1, false );
    std::vector<bool> reco_used( reco_indices.empty() ? 0 : *std::max_element( reco_indices.begin(), reco_indices.end() ) + 1, false );

    std::vector<JetMatch> matched;
    for ( const auto & c : candidates )
    {
        if ( truth_used[c.truth_index] || reco_used[c.reco_index] ) continue;
        truth_used[c.truth_index] = true;
        reco_used[c.reco_index] = true;
        matched.push_back( { c.truth_index, c.reco_index, c.dr } );
    }

    std::vector<JetMatch> result = matched;
    for ( const int ti : truth_indices )
    {
        if ( !truth_used[ti] ) result.push_back( { ti, -1, -1.0f } );
    }
    for ( const int ri : reco_indices )
    {
        if ( !reco_used[ri] ) result.push_back( { -1, ri, -1.0f } );
    }
    return result;
}

void AnaUtils::book_matched_jet_tree( TTree * tree, MatchedJetRow & row )
{
    tree -> Branch( "event_id", &row.event_id, "event_id/I" );
    tree -> Branch( "cent", &row.cent, "cent/I" );
    tree -> Branch( "zvrtx", &row.zvrtx, "zvrtx/F" );
    tree -> Branch( "mbdQ", &row.mbdQ, "mbdQ/F" );
    tree -> Branch( "sumeT", &row.sumeT, "sumeT/F" );
    tree -> Branch( "is_minbias", &row.is_minbias, "is_minbias/I" );
    tree -> Branch( "psi2", &row.psi2, "psi2/F" );
    tree -> Branch( "truth_jet_maxpt_r04", &row.truth_jet_maxpt_r04, "truth_jet_maxpt_r04/F" );

    tree -> Branch( "reco_type", &row.reco_type, "reco_type/I" );
    tree -> Branch( "match_status", &row.match_status, "match_status/I" );
    tree -> Branch( "dr", &row.dr, "dr/F" );

    tree -> Branch( "truth_pt", &row.truth_pt, "truth_pt/F" );
    tree -> Branch( "truth_e", &row.truth_e, "truth_e/F" );
    tree -> Branch( "truth_eta", &row.truth_eta, "truth_eta/F" );
    tree -> Branch( "truth_phi", &row.truth_phi, "truth_phi/F" );
    tree -> Branch( "truth_flavor", &row.truth_flavor, "truth_flavor/I" );

    tree -> Branch( "reco_pt", &row.reco_pt, "reco_pt/F" );
    tree -> Branch( "reco_e", &row.reco_e, "reco_e/F" );
    tree -> Branch( "reco_eta", &row.reco_eta, "reco_eta/F" );
    tree -> Branch( "reco_phi", &row.reco_phi, "reco_phi/F" );
    tree -> Branch( "reco_unsub_e", &row.reco_unsub_e, "reco_unsub_e/F" );
    tree -> Branch( "reco_unsub_pt", &row.reco_unsub_pt, "reco_unsub_pt/F" );
}

void AnaUtils::read_matched_jet_tree( TTree * tree, MatchedJetRow & row )
{
    tree -> SetBranchAddress( "event_id", &row.event_id );
    tree -> SetBranchAddress( "cent", &row.cent );
    tree -> SetBranchAddress( "zvrtx", &row.zvrtx );
    tree -> SetBranchAddress( "mbdQ", &row.mbdQ );
    tree -> SetBranchAddress( "sumeT", &row.sumeT );
    tree -> SetBranchAddress( "is_minbias", &row.is_minbias );
    tree -> SetBranchAddress( "psi2", &row.psi2 );
    tree -> SetBranchAddress( "truth_jet_maxpt_r04", &row.truth_jet_maxpt_r04 );

    tree -> SetBranchAddress( "reco_type", &row.reco_type );
    tree -> SetBranchAddress( "match_status", &row.match_status );
    tree -> SetBranchAddress( "dr", &row.dr );

    tree -> SetBranchAddress( "truth_pt", &row.truth_pt );
    tree -> SetBranchAddress( "truth_e", &row.truth_e );
    tree -> SetBranchAddress( "truth_eta", &row.truth_eta );
    tree -> SetBranchAddress( "truth_phi", &row.truth_phi );
    tree -> SetBranchAddress( "truth_flavor", &row.truth_flavor );

    tree -> SetBranchAddress( "reco_pt", &row.reco_pt );
    tree -> SetBranchAddress( "reco_e", &row.reco_e );
    tree -> SetBranchAddress( "reco_eta", &row.reco_eta );
    tree -> SetBranchAddress( "reco_phi", &row.reco_phi );
    tree -> SetBranchAddress( "reco_unsub_e", &row.reco_unsub_e );
    tree -> SetBranchAddress( "reco_unsub_pt", &row.reco_unsub_pt );
}

