#ifndef _FUN4ALL_UESCALING_PASS2_C_
#define _FUN4ALL_UESCALING_PASS2_C_

#include <GlobalVariables.C>

#include <G4_CEmc_Spacal.C>
#include <G4_HcalIn_ref.C>
#include <G4_HcalOut_ref.C>
#include <HIJetReco.C>
#include <G4_Input.C>

#include <Calo_Calib.C>

#include <ffamodules/CDBInterface.h>

#include <phool/recoConsts.h>

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllUtils.h>

#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/CaloTowerStatus.h>

#include <calowaveformsim/CaloWaveformSim.h>

#include <calomanip/CaloTree.h>
#include <calomanip/CaloManip.h>

R__LOAD_LIBRARY( libfun4all.so )
R__LOAD_LIBRARY( libffamodules.so )
R__LOAD_LIBRARY( libCaloWaveformSim.so )
R__LOAD_LIBRARY( libcalo_reco.so )
R__LOAD_LIBRARY( libcalomanip.so )


void Fun4All_UEScaling_Pass2 (
    const int nEvents               = 10,
    const int run_number            = 31,
    const int segment               = 0,
    const int jet_flag              = 10,
    const std::string & embfile     = "CALO_TREE_noNoise_waveformFit_sHijing_0_20fm-00000031-00000.root",
    const std::string & outfile     = "DST_CALO_CLUSTER_pythia8_Jet10_scaled11perc_sHijing_0_20fm-00000031-00000.root"
)
{
    
    std::cout << "Fun4All_UEScaling_Pass2" << std::endl;

    Enable::VERBOSITY        = 0;

    const std::string & cdbtag = "MDC2";
    Enable::CDB = true;
        
    auto * se = Fun4AllServer::instance();
    se -> Verbosity( 1 );

    auto * rc = recoConsts::instance();
    rc -> set_StringFlag( "CDB_GLOBALTAG", cdbtag );
    rc -> set_uint64Flag( "TIMESTAMP", run_number );

    auto * cdb = CDBInterface::instance();
    cdb -> Verbosity( Enable::VERBOSITY );
    
    auto * flag = new FlagHandler();
    se -> registerSubsystem( flag );

    Input::VERBOSITY = 1;
    Input::READHITS = false;
    // INPUTREADHITS::filename[0] = Form( "G4Hits_sHijing_0_20fm-%010d-%06d.root", run_number, segment );
    InputInit();
    InputRegister();

    for ( const auto & DSTTPYE : { "DST_CALO_CLUSTER" , "DST_GLOBAL",  "DST_MBD_EPD", "DST_TRUTH_JET"} )
    {
        std::string infile = Form( "%s_pythia8_Jet%d_sHijing_0_20fm-%010d-%06d.root", DSTTPYE, jet_flag, run_number, segment );
        std::cout << "\tAdding input file: " << infile << std::endl;
        auto input = new Fun4AllDstInputManager( Form( "DSTINPUT_%s", DSTTPYE ) );
        input -> AddFile( infile );
        input -> Verbosity( Enable::VERBOSITY );
        se -> registerInputManager( input );
    }

    InputManagers();

    auto * ingeom = new Fun4AllRunNodeInputManager( "DST_GEO" );
    ingeom -> AddFile( CDBInterface::instance() -> getUrl( "calo_geo" ) );
    se -> registerInputManager( ingeom );

    Process_Calo_Calib( ); 
    
    auto * cm  = new CaloManip( embfile.c_str() );
    cm -> add_event_header ( "EventHeader" );
    cm -> add_cemc_node ( "TOWERINFO_CALIB_CEMC" );
    cm -> add_hcalin_node ( "TOWERINFO_CALIB_HCALIN" );
    cm -> add_hcalout_node ( "TOWERINFO_CALIB_HCALOUT" );
    cm -> set_scale_factor ( 1.11 );
    cm -> Verbosity ( Enable::VERBOSITY  );
    se -> registerSubsystem( cm );

    auto * rcemc = new RetowerCEMC( ); 
    rcemc -> set_towerinfo( true );
    rcemc -> set_frac_cut( 1.0 );
    rcemc -> set_do_rescale( false );
    rcemc -> set_towerNodePrefix( HIJETS::tower_prefix );
    rcemc -> Verbosity( Enable::VERBOSITY );
    se -> registerSubsystem( rcemc );

    auto * out = new Fun4AllDstOutputManager( "DSTOUTPUT", outfile );
    // if these are left then process_calo_calib will overwrite the overlayed enregies
    out -> StripNode( "TOWERS_CEMC" );
    out -> StripNode( "TOWERS_HCALIN" );
    out -> StripNode( "TOWERS_HCALOUT" );
    out -> Verbosity( Enable::VERBOSITY );
    se -> registerOutputManager( out );

    se -> run( nEvents );
    se -> End( );   
    se -> PrintTimer( );

    CDBInterface::instance() -> Print();

    delete se;

    std::cout << "Done4All" << std::endl;

    gSystem -> Exit( 0 );

}

#endif // _FUN4ALL_UESCALING_PASS2_C_