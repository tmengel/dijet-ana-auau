#ifndef _FUN4ALL_DIJETS_AUAU_C_
#define _FUN4ALL_DIJETS_AUAU_C_

#include <GlobalVariables.C>

#include <G4_CEmc_Spacal.C>
#include <G4_HcalIn_ref.C>
#include <G4_HcalOut_ref.C>
#include <G4_Input.C>

#include <HIJetReco.C>

#include <Calo_Calib.C>

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllUtils.h>

#include <ffamodules/CDBInterface.h>

#include <phool/recoConsts.h>

#include <mbd/MbdReco.h>

#include <epd/EpdReco.h>

#include <zdcinfo/ZdcReco.h>

#include <globalvertex/GlobalVertexReco.h>

#include <centrality/CentralityReco.h>

#include <calotrigger/MinimumBiasClassifier.h>

#include <jetbase/FastJetOptions.h>
#include <jetbase/JetCalib.h>

#include <calotrigger/TriggerRunInfoReco.h>

#include <jetbackground/DetermineTowerRho.h>
#include <jetbackground/TowerRho.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowersRho.h>
#include <jetbackground/SubtractTowersRhov1.h>

#include <myana/EventSelector.h>
#include <myana/MinBiasCut.h>
#include <myana/ZVertexCut.h>
#include <myana/TriggerSelect.h>
#include <myana/AnaTreev1.h>
#include <myana/RhoEtaCalibLookup.h>

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <mbd/MbdEvent.h>
#include <mbd/MbdCalib.h>

#include <string>

R__LOAD_LIBRARY( libfun4all.so )
R__LOAD_LIBRARY( libffamodules.so )
R__LOAD_LIBRARY( libmbd.so )
R__LOAD_LIBRARY( libepd.so )
R__LOAD_LIBRARY( libzdcinfo.so )
R__LOAD_LIBRARY( libglobalvertex.so )
R__LOAD_LIBRARY( libcentrality.so )
R__LOAD_LIBRARY( libcalotrigger.so )
R__LOAD_LIBRARY( libjetbase.so )
R__LOAD_LIBRARY( libjetbackground.so )

R__LOAD_LIBRARY( libmyana.so )


TowerJetInput * GetTowerInput(
    const Jet::SRC src,
    const std::string & prefix = "TOWERINFO_CALIB"
);

void Fun4All_Dijets_AuAu (
    const int nEvents                    = 3,
    const std::string & infile_calo      = "DST_CALOFITTING_run2auau_pro001_pcdb001_v001-00054912-00000.root",
    const std::string & infile_zdc       = "DST_ZDC_RAW_run2auau_pro001_pcdb001_v001-00054912-00000.root",
    const std::string & infile_sepd      = "DST_SEPD_RAW_run2auau_pro001_pcdb001_v001-00054912-00000.root",
    const std::string & outfile_ana      = "output_ana.root",
    const std::string & cdbtag           = "newcdbtag"
)
{
    
    Enable::VERBOSITY = 0;

    auto runseg    = Fun4AllUtils::GetRunSegment( infile_calo );
    int run_number = runseg.first;
    int segment    = runseg.second;

    auto * se = Fun4AllServer::instance( );
    se -> Verbosity( 1 );

    auto * rc = recoConsts::instance( );
    rc -> set_StringFlag( "CDB_GLOBALTAG", cdbtag );
    rc -> set_uint64Flag( "TIMESTAMP", run_number );

    auto * cdb = CDBInterface::instance();
    cdb -> Verbosity( Enable::VERBOSITY );

    auto * flag = new FlagHandler();
    se -> registerSubsystem( flag );

    int ifile = 0;
    for ( const auto & infile : { infile_calo, infile_zdc, infile_sepd } )
    {
        std::cout << "Input file " << ifile << ": " << infile << std::endl;
        auto * input = new Fun4AllDstInputManager( Form( "DSTINPUT_%d", ifile ) );
        input -> AddFile( infile );
        input -> Verbosity( Enable::VERBOSITY );
        se -> registerInputManager( input );
        ++ifile;
    }

    CaloTowerDefs::BuilderType buildertype = CaloTowerDefs::kPRDFTowerv4;
    
    auto * ingeom = new Fun4AllRunNodeInputManager( "DST_GEO" );
    auto geoLocation = CDBInterface::instance() -> getUrl( "calo_geo" );
    ingeom -> AddFile( geoLocation );
    se -> registerInputManager( ingeom );

    auto * trig = new TriggerRunInfoReco( );
    se -> registerSubsystem(trig);

    auto * mbd = new MbdReco( );
    se -> registerSubsystem( mbd );

    auto * sepdbuilder = new CaloTowerBuilder( "SEPDBUILDER" );
    sepdbuilder -> set_detector_type( CaloTowerDefs::SEPD );
    sepdbuilder -> set_builder_type( buildertype );
    sepdbuilder -> set_processing_type( CaloWaveformProcessing::TEMPLATE );
    sepdbuilder -> set_nsamples( 12 );
    sepdbuilder -> set_offlineflag( );
    se -> registerSubsystem( sepdbuilder );

    auto * sepd = new EpdReco( );
    se -> registerSubsystem(sepd);

    auto * caZDC = new CaloTowerBuilder( "ZDCBUILDER" );
    caZDC -> set_detector_type( CaloTowerDefs::ZDC );
    caZDC -> set_builder_type( buildertype );
    caZDC -> set_processing_type( CaloWaveformProcessing::FAST );
    caZDC -> set_nsamples( 16 );
    caZDC -> set_offlineflag( );
    se -> registerSubsystem( caZDC );

    auto * zdc = new ZdcReco( );
    zdc -> set_zdc1_cut( 0.0 );
    zdc -> set_zdc2_cut( 0.0 );
    se -> registerSubsystem( zdc );
    
    auto * zvtrx = new GlobalVertexReco( );
    se -> registerSubsystem( zvtrx );

    Process_Calo_Calib( ); 

    auto * epreco = new EventPlaneReco( );
    epreco -> set_inputNode( "TOWERINFO_CALIB_SEPD" );
    epreco -> set_EventPlaneInfoNodeName( "EventplaneinfoMap" );
    se -> registerSubsystem( epreco );

    auto * mb = new MinimumBiasClassifier( );
    se -> registerSubsystem( mb );

    auto * cr = new CentralityReco( );
    cr -> Verbosity( Enable::VERBOSITY );
    se -> registerSubsystem( cr );

    auto * rcemc = new RetowerCEMC( ); 
    rcemc -> set_towerinfo( true );
    rcemc -> set_frac_cut( 1.0 );
    rcemc -> set_do_rescale( false );
    rcemc -> set_towerNodePrefix( "TOWERINFO_CALIB" );
    se -> registerSubsystem( rcemc );

    auto * es = new EventSelector( );
    es -> Verbosity( Enable::VERBOSITY );
    // auto * tcut = new TriggerSelect();
    // tcut -> SetPacket( 14001 );
    // tcut -> SelectTrigger( 10 );
    // es -> AddCut( tcut );
    auto * mbcut = new MinBiasCut( );
    mbcut -> SetNodeName( "MinimumBiasInfo" );
    es -> AddCut( mbcut );
    es -> PrintCuts( );
    se -> registerSubsystem( es );

    auto * tjr = new JetReco( "TowerJetReco" );
    for ( const auto & src : { Jet::CEMC_TOWERINFO_RETOWER, Jet::HCALIN_TOWERINFO, Jet::HCALOUT_TOWERINFO } )
    {
        tjr -> add_input( GetTowerInput( src , "TOWERINFO_CALIB" ) );
    }
    tjr -> add_algo( HIJETS::GetFJAlgo( 0.2 ), "AntiKt_TowerInfo_HIRecoSeedsRaw_r02" );
    for ( const auto & R : { 0.3 } )
    {
        tjr -> add_algo( HIJETS::GetFJAlgo( R ), Form( "AntiKt_TowerInfo_r0%d", static_cast<int>( R * 10 ) ) );
    }
    tjr -> set_algo_node( "ANTIKT" );
    tjr -> set_input_node( "TOWER" );
    tjr -> Verbosity( 0 );
    se -> registerSubsystem( tjr );

    auto * dtb1 = new DetermineTowerBackground( "DetermineTowerBackground1" );
    dtb1 -> SetBackgroundOutputName( "TowerInfoBackground_Sub1" );
    dtb1 -> set_towerNodePrefix( "TOWERINFO_CALIB" );
    dtb1 -> SetFlow( 0 );
    dtb1 -> SetSeedJetPt( 5.0 );
    dtb1 -> SetSeedJetD( 3.0 );
    dtb1 -> SetSeedType( 0 );
    dtb1 -> SetSeedMaxConst( 3.0 );
    dtb1 -> Verbosity( 0 );
    dtb1 -> UseReweighting( true );
    se -> registerSubsystem( dtb1 );

    auto * casj1 = new CopyAndSubtractJets( "CopyAndSubtractJets" );
    casj1 -> set_towerNodePrefix( "TOWERINFO_CALIB" );
    casj1 -> SetFlowModulation( 0 );
    casj1 -> Verbosity( 0 );
    casj1 -> set_towerinfo( true );
    se -> registerSubsystem( casj1 );

    auto * dtb2 = new DetermineTowerBackground( "DetermineTowerBackground_Sub1" );
    dtb2 -> SetBackgroundOutputName( "TowerInfoBackground_Sub2" );
    dtb2 -> set_towerNodePrefix( "TOWERINFO_CALIB" );
    dtb2 -> SetFlow( 0 );
    dtb2 -> SetSeedJetPt( 7.0 );
    dtb2 -> SetSeedType( 1 );
    dtb2 -> Verbosity( 0 );
    dtb2 -> UseReweighting( true );
    se -> registerSubsystem( dtb2 );

    auto * st1 = new SubtractTowers( "SubtractTowers" );
    st1 -> set_towerNodePrefix( "TOWERINFO_CALIB" );
    st1 -> SetFlowModulation( 0 );
    st1 -> Verbosity( 0 );
    st1 -> set_towerinfo( true );
    se -> registerSubsystem( st1 );

    tjr = new JetReco( "TowerJetReco_Sub1" );
    for ( const auto & src : { Jet::CEMC_TOWERINFO_SUB1 , Jet::HCALIN_TOWERINFO_SUB1, Jet::HCALOUT_TOWERINFO_SUB1 } )
    {
        tjr -> add_input( GetTowerInput( src , "TOWERINFO_CALIB" ) );
    }
    for ( const auto & R : { 0.3 } )
    {
        tjr -> add_algo( HIJETS::GetFJAlgo( R ), Form( "AntiKt_TowerInfo_r0%d_Sub1", static_cast<int>( R * 10 ) ) );
    }
    tjr -> set_algo_node( "ANTIKT" );
    tjr -> set_input_node( "TOWER" );
    tjr -> Verbosity( 0 );
    se -> registerSubsystem( tjr );
    
    auto * trc = new DetermineTowerRho( "DetermineTowerRho_CEMC_Mult" );
    trc -> add_method( TowerRho::Method::MULT, "TowerRho_MULT_CEMC" );
    trc -> add_tower_input( GetTowerInput( Jet::CEMC_TOWERINFO_RETOWER ) );
    se -> registerSubsystem( trc );

    trc = new DetermineTowerRho( "DetermineTowerRho_HCALIN_Mult" );
    trc -> add_method( TowerRho::Method::MULT, "TowerRho_MULT_HCALIN" );
    trc -> add_tower_input( GetTowerInput( Jet::HCALIN_TOWERINFO ) );
    se -> registerSubsystem( trc );  
    
    trc = new DetermineTowerRho( "DetermineTowerRho_HCALOUT_Mult" );
    trc -> add_method( TowerRho::Method::MULT, "TowerRho_MULT_HCALOUT" );
    trc -> add_tower_input( GetTowerInput( Jet::HCALOUT_TOWERINFO ) );
    se -> registerSubsystem( trc );

    const std::string rho_eta_calib_path = RhoEtaCalibLookup::GetCalibPath( run_number );

    auto * subrho = new SubtractTowersRhov1(  "SubtractTowersRho_CEMC_Mult" );
    subrho -> set_rhoNode("TowerRho_MULT_CEMC");
    subrho -> add_targetTowerNode( "TOWERINFO_CALIB_CEMC_RETOWER" );
    subrho -> set_etaCalib_directPath( rho_eta_calib_path );
    subrho -> set_subSuffix( "MULTSUB" );
    subrho -> set_globalVertexType( GlobalVertex::MBD );
    se -> registerSubsystem( subrho );

    subrho = new SubtractTowersRhov1( "SubtractTowersRho_HCALIN_Mult" );
    subrho -> set_rhoNode("TowerRho_MULT_HCALIN");
    subrho -> add_targetTowerNode( "TOWERINFO_CALIB_HCALIN" );
    subrho -> set_etaCalib_directPath( rho_eta_calib_path );
    subrho -> set_subSuffix( "MULTSUB" );
    subrho -> set_globalVertexType( GlobalVertex::MBD );
    se -> registerSubsystem( subrho );

    subrho = new SubtractTowersRhov1( "SubtractTowersRho_HCALOUT_Mult" );
    subrho -> set_rhoNode("TowerRho_MULT_HCALOUT");
    subrho -> add_targetTowerNode( "TOWERINFO_CALIB_HCALOUT" );
    subrho -> set_etaCalib_directPath( rho_eta_calib_path );
    subrho -> set_subSuffix( "MULTSUB" );
    subrho -> set_globalVertexType( GlobalVertex::MBD );
    se -> registerSubsystem( subrho );

    tjr = new JetReco( "TowerJetReco_Rho" );
    // these are the unsub tower types but have been updated w/ multsub
    for ( const auto & src : { Jet::CEMC_TOWERINFO_RETOWER , Jet::HCALIN_TOWERINFO, Jet::HCALOUT_TOWERINFO } )
    {
        tjr -> add_input( GetTowerInput( src , "MULTSUB_TOWERINFO_CALIB" ) );
    }
    for ( const auto & R : { 0.3 } )
    {
        tjr -> add_algo( HIJETS::GetFJAlgo( R ), Form( "AntiKt_TowerInfo_r0%d_Rho", static_cast<int>( R * 10 ) ) );
    }
    tjr -> set_algo_node( "ANTIKT" );
    tjr -> set_input_node( "TOWER" );
    tjr -> Verbosity( 0 );
    se -> registerSubsystem( tjr );

    // JES calibration (legacy z-vertex + eta dependent TF1s, CDB JES_Calib_Default)
    // for both subtracted jet collections; output keeps the input jet order
    for ( const auto & R : { 0.3 } )
    {
        const int cone = static_cast<int>( R * 10 );
        for ( const std::string & jettype : { "Sub1", "Rho" } )
        {
            auto * jetCalib = new JetCalib( Form( "JetCalib_r%02d_%s", cone, jettype.c_str() ) );
            jetCalib -> set_InputNode( Form( "AntiKt_TowerInfo_r%02d_%s", cone, jettype.c_str() ) );
            jetCalib -> set_OutputNode( Form( "AntiKt_TowerInfo_r%02d_%s_calib", cone, jettype.c_str() ) );
            jetCalib -> set_JetRadius( R );
            jetCalib -> set_UseEMfracCalib( false );
            jetCalib -> set_ApplyZvrtxDependentCalib( true );
            jetCalib -> set_ApplyEtaDependentCalib( true );
            jetCalib -> Verbosity( Enable::VERBOSITY );
            se -> registerSubsystem( jetCalib );
        }
    }

    auto * anaout = new AnaTreev1( outfile_ana );
    anaout -> Verbosity( Enable::VERBOSITY );
    anaout -> add_zvrtx_node( "GlobalVertexMap" );
    anaout -> add_cent_node( "CentralityInfo" );
    anaout -> add_gl1_node( "14001" );
    anaout -> add_mbd_node( "MbdOut" );
    anaout -> add_minbias_node( "MinimumBiasInfo" );
    anaout -> add_sub1_jet_node( "AntiKt_TowerInfo_r03_Sub1", "TowerInfoBackground_Sub2" );
    anaout -> add_sub1_jet_calib_node( "AntiKt_TowerInfo_r03_Sub1_calib" );
    anaout -> add_towerbkgd_v2_node( "TowerInfoBackground_Sub2" );
    anaout -> add_rho_jet_node( "AntiKt_TowerInfo_r03_Rho" , "TowerRho_MULT_CEMC", "TowerRho_MULT_HCALIN", "TowerRho_MULT_HCALOUT", "MULTSUB_TOWERINFO_CALIB" );
    anaout -> add_rho_jet_calib_node( "AntiKt_TowerInfo_r03_Rho_calib" );
    for ( const auto & rho_node : { "TowerRho_MULT_CEMC", "TowerRho_MULT_HCALIN", "TowerRho_MULT_HCALOUT" } )
    {
        anaout -> add_rho_nodes( rho_node );
    }
    anaout -> add_jet_node( "AntiKt_TowerInfo_r03" );
    anaout -> add_cemc_node(  "TOWERINFO_CALIB_CEMC_RETOWER" , true );
    anaout -> add_ihcal_node( "TOWERINFO_CALIB_HCALIN" , true );
    anaout -> add_ohcal_node( "TOWERINFO_CALIB_HCALOUT" , true );

    se -> registerSubsystem( anaout );

    se -> run( nEvents );

    se -> End( );   

    se -> PrintTimer( );

    CDBInterface::instance() -> Print();

    delete se;

    std::cout << "Done4All" << std::endl;

    gSystem -> Exit( 0 );

}

TowerJetInput * GetTowerInput(
    const Jet::SRC src,
    const std::string & prefix
)
{
    auto * input = new TowerJetInput( src, prefix );
    input -> set_GlobalVertexType( GlobalVertex::MBD );
    return input;
}

#endif