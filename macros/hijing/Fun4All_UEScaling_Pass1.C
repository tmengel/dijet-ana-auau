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

void Fun4All_UEScaling_Pass1 ( 
    const int nEvents                   = 10,
    const int run_number                = 31,
    const int segment                   = 0,
    const std::string & outfile         = "CALO_TREE_noNoise_sHijing_0_20fm-00000031-00000.root",
    const bool do_waveform_fit          = true
)
{
    std::cout << "Fun4All_UEScaling_Pass1" << std::endl;

    Enable::VERBOSITY        = 0;

    const std::string & cdbtag = "MDC2";
    Enable::CDB = true;
        
    const CaloWaveformSim::NoiseType noise_type = CaloWaveformSim::NOISE_NONE;

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
    
    Input::READHITS = true && do_waveform_fit;
    INPUTREADHITS::filename[0] = Form( "G4Hits_sHijing_0_20fm-%010d-%06d.root", run_number, segment );
    
    InputInit();
    
    InputRegister();
    
    for ( const auto & DSTTPYE : { "DST_GLOBAL",  "DST_MBD_EPD" , "DST_CALO_CLUSTER" } ) 
    {
        if ( std::string( DSTTPYE ).find( "CALO_CLUSTER" ) != std::string::npos && do_waveform_fit ) 
        {
            continue; // skip the calo cluster input if we are doing waveform fit, since we will build it from hits
        }
        std::string infile = Form( "%s_sHijing_0_20fm-%010d-%06d.root", DSTTPYE, run_number, segment );
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

    if ( do_waveform_fit )
    {
        auto * cws_hcalout = new CaloWaveformSim( "HCALOUTWAVEFORMSIM" );
        cws_hcalout -> set_detector_type( CaloTowerDefs::HCALOUT );
        cws_hcalout -> set_detector( "HCALOUT" );
        cws_hcalout -> set_nsamples( 12 );
        cws_hcalout -> set_timewidth( 0.2 );
        cws_hcalout -> set_peakpos( 6 );
        cws_hcalout -> set_noise_type( noise_type );
        se -> registerSubsystem( cws_hcalout );

        auto * cws_hcalin = new CaloWaveformSim( "HCALINWAVEFORMSIM" );
        cws_hcalin -> set_detector_type( CaloTowerDefs::HCALIN );
        cws_hcalin -> set_detector( "HCALIN" );
        cws_hcalin -> set_nsamples( 12 );
        cws_hcalin -> set_timewidth( 0.2 );
        cws_hcalin -> set_peakpos( 6 );
        cws_hcalin -> set_noise_type( noise_type );
        se -> registerSubsystem( cws_hcalin );

        auto * cws_cemc = new CaloWaveformSim( "CEMCWAVEFORMSIM" );
        cws_cemc -> set_detector_type( CaloTowerDefs::CEMC );
        cws_cemc -> set_detector( "CEMC" );
        cws_cemc -> set_nsamples( 12 );
        cws_cemc -> set_timewidth( 0.2 );
        cws_cemc -> set_peakpos( 6 );
        cws_cemc -> set_noise_type( noise_type );
        cws_cemc -> get_light_collection_model().load_data_file(
            std::string( getenv("CALIBRATIONROOT") ) +
            std::string( "/CEMC/LightCollection/Prototype3Module.xml" ),
            "data_grid_light_guide_efficiency", "data_grid_fiber_trans" 
        );
        se -> registerSubsystem( cws_cemc );

        auto * ctb_hcalout = new CaloTowerBuilder( "HCALOUTTOWERBUILDER" );
        ctb_hcalout -> set_detector_type( CaloTowerDefs::HCALOUT );
        ctb_hcalout -> set_nsamples( 12 );
        ctb_hcalout -> set_dataflag( false );
        ctb_hcalout -> set_processing_type( CaloWaveformProcessing::TEMPLATE );
        // match our current ZS threshold ~7ADC for hcal
        ctb_hcalout -> set_softwarezerosuppression( true, 7 );
        se -> registerSubsystem( ctb_hcalout );

        auto * ctb_hcalin = new CaloTowerBuilder( "HCALINTOWERBUILDER" );
        ctb_hcalin -> set_detector_type( CaloTowerDefs::HCALIN );
        ctb_hcalin -> set_nsamples( 12 );
        ctb_hcalin -> set_dataflag( false );
        ctb_hcalin -> set_processing_type( CaloWaveformProcessing::TEMPLATE );
        ctb_hcalin -> set_softwarezerosuppression( true, 7 );
        se -> registerSubsystem( ctb_hcalin );

        auto * ctb_cemc = new CaloTowerBuilder( "CEMCTOWERBUILDER" );
        ctb_cemc -> set_detector_type( CaloTowerDefs::CEMC );
        ctb_cemc -> set_nsamples( 12 );
        ctb_cemc -> set_dataflag( false );
        ctb_cemc -> set_processing_type( CaloWaveformProcessing::TEMPLATE );
        // match our current ZS threshold ~14ADC for emcal
        ctb_cemc -> set_softwarezerosuppression( true, 14 );
        se -> registerSubsystem( ctb_cemc );

        auto * cts_emcal = new CaloTowerStatus( "CEMCSTATUS" );
        cts_emcal -> set_detector_type( CaloTowerDefs::CEMC );
        se -> registerSubsystem( cts_emcal );

        auto * cts_hcalin = new CaloTowerStatus( "HCALINSTATUS" );
        cts_hcalin -> set_detector_type( CaloTowerDefs::HCALIN );
        se -> registerSubsystem( cts_hcalin );

        auto * cts_hcalout = new CaloTowerStatus( "HCALOUTSTATUS" );
        cts_hcalout -> set_detector_type( CaloTowerDefs::HCALOUT );
        se -> registerSubsystem( cts_hcalout );

        auto * ctc_emcal = new CaloTowerCalib( "CEMCCALIB" );
        ctc_emcal -> set_detector_type( CaloTowerDefs::CEMC );
        ctc_emcal -> set_outputNodePrefix( "TOWERINFO_CALIB_" );
        se -> registerSubsystem( ctc_emcal );

        auto * ctc_hcalin = new CaloTowerCalib( "HCALINCALIB" );
        ctc_hcalin -> set_detector_type( CaloTowerDefs::HCALIN );
        ctc_hcalin -> set_outputNodePrefix( "TOWERINFO_CALIB_" );
        se -> registerSubsystem( ctc_hcalin );

        auto * ctc_hcalout = new CaloTowerCalib( "HCALOUTCALIB" );
        ctc_hcalout -> set_detector_type( CaloTowerDefs::HCALOUT );
        ctc_hcalout -> set_outputNodePrefix( "TOWERINFO_CALIB_" );
        se -> registerSubsystem( ctc_hcalout );
    }

    Process_Calo_Calib( );

    // auto * rcemc = new RetowerCEMC( );
    // rcemc -> set_towerinfo( true );
    // rcemc -> set_frac_cut( 1.0 );
    // rcemc -> set_do_rescale( false );
    // rcemc -> set_towerNodePrefix( HIJETS::tower_prefix );
    // rcemc -> Verbosity( Enable::VERBOSITY );
    // se -> registerSubsystem( rcemc );

    auto * out = new CaloTree( outfile );
    out -> add_event_header ( "EventHeader" );
    out -> add_cemc_node ( "TOWERINFO_CALIB_CEMC" );
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