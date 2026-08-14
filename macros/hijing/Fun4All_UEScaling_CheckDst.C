#ifndef _FUN4ALL_UESCALING_PASS1_C_
#define _FUN4ALL_UESCALING_PASS1_C_

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

void Fun4All_UEScaling_CheckDst( 
    const int nEvents                   = -1,
    const std::string & infile          = "DST_SCALED_jet10_hijing31_pass2-00000.root",
    const std::string & outfile         = "CALO_TREE_DST_CHECK_sHijing_0_20fm-00000031-00000.root",
    const bool nominal                    = false
)
{
    std::cout << "Fun4All_UEScaling_Pass1" << std::endl;

    Enable::VERBOSITY        = 0;

    const std::string & cdbtag = "MDC2";
    Enable::CDB = true;
        
    const CaloWaveformSim::NoiseType noise_type = CaloWaveformSim::NOISE_NONE;

    const int run_number                    = 31;
    const int segment                       = 0;
    const int jet_flag                      = 10;

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

    if ( nominal )
    {
        for ( const auto & DSTTPYE : { "DST_CALO_CLUSTER" , "DST_GLOBAL",  "DST_MBD_EPD", "DST_TRUTH_JET"} )
        {
            std::string infile = Form( "%s_pythia8_Jet%d_sHijing_0_20fm-%010d-%06d.root", DSTTPYE, jet_flag, run_number, segment );
            std::cout << "\tAdding input file: " << infile << std::endl;
            auto input = new Fun4AllDstInputManager( Form( "DSTINPUT_%s", DSTTPYE ) );
            input -> AddFile( infile );
            input -> Verbosity( Enable::VERBOSITY );
            se -> registerInputManager( input );
        }

    }
    else
    {
        auto * input = new Fun4AllDstInputManager( "DSTINPUT" );
        input -> AddFile( infile );
        input -> Verbosity( Enable::VERBOSITY );
        se -> registerInputManager( input );
    }
    // auto * input = new Fun4AllDstInputManager( "DSTINPUT" );
    // input -> AddFile( infile );
    // input -> Verbosity( Enable::VERBOSITY );
    // se -> registerInputManager( input );


    auto * ingeom = new Fun4AllRunNodeInputManager( "DST_GEO" );
    ingeom -> AddFile( CDBInterface::instance() -> getUrl( "calo_geo" ) );
    se -> registerInputManager( ingeom );

    if ( nominal )
    {
        Process_Calo_Calib( );
    }

    auto * rcemc = new RetowerCEMC( );
    rcemc -> set_towerinfo( true );
    rcemc -> set_frac_cut( 1.0 );
    rcemc -> set_do_rescale( false );
    rcemc -> set_towerNodePrefix( HIJETS::tower_prefix );
    rcemc -> Verbosity( Enable::VERBOSITY );
    se -> registerSubsystem( rcemc );

    auto * out = new CaloTree( outfile );
    out -> add_event_header ( "EventHeader" );
    out -> add_cemc_node ( "TOWERINFO_CALIB_CEMC_RETOWER" );
    out -> add_hcalin_node ( "TOWERINFO_CALIB_HCALIN" );
    out -> add_hcalout_node ( "TOWERINFO_CALIB_HCALOUT" );
    out -> Verbosity( Enable::VERBOSITY  );
    se -> registerSubsystem( out );

    se -> run( nEvents );
    se -> End( );   
    se -> PrintTimer( );

    CDBInterface::instance() -> Print();

    delete se;

    std::cout << "Done4All" << std::endl;

    gSystem -> Exit( 0 );

}

#endif // _FUN4ALL_UESCALING_PASS1_C_