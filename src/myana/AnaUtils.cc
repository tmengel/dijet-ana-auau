#include "AnaUtils.h"

#include <fstream>
#include <iostream>
#include <cmath>

#include <TMath.h>
#include <TLatex.h>

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

