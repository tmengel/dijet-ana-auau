#ifndef _FUN4ALL_UESCALING_PASS3_C_
#define _FUN4ALL_UESCALING_PASS3_C_

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

#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/DetermineTowerBackgroundv1.h>
#include <jetbackground/DetermineTowerRho.h>
#include <jetbackground/TowerRho.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersRhov1.h>

#include <myana/EventSelector.h>
#include <myana/MinBiasCut.h>
#include <myana/ZVertexCut.h>
#include <myana/TriggerSelect.h>
#include <myana/AnaTreev1.h>
#include <myana/RhoEtaCalibLookup.h>

R__LOAD_LIBRARY( libfun4all.so )
R__LOAD_LIBRARY( libffamodules.so )

R__LOAD_LIBRARY( libcalotrigger.so )
R__LOAD_LIBRARY( libcentrality.so )

R__LOAD_LIBRARY( libmbd.so )
R__LOAD_LIBRARY( libepd.so )
R__LOAD_LIBRARY( libzdcinfo.so )
R__LOAD_LIBRARY( libglobalvertex.so )

R__LOAD_LIBRARY( libjetbackground.so )
R__LOAD_LIBRARY( libjetbase.so )

R__LOAD_LIBRARY( libmyana.so )

TowerJetInput * GetTowerInput( 
    const Jet::SRC src, 
    const std::string & prefix = "TOWERINFO_CALIB" 
);

void Fun4All_UEScaling_NoScale (
    const int nEvents               = 10,
    const std::string & outfile     = "TREE_SCALED_jet10_hijing31_pass3.root"
)
{
    
    std::cout << "Fun4All_UEScaling_Pass3" << std::endl;

    Enable::VERBOSITY = 0;

    const std::string & cdbtag = "MDC2";
    const int run_number = 31;
    const int segment = 0;
    const int jet_flag = 20;
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

    for ( const auto & DSTTPYE : { "DST_CALO_CLUSTER" , "DST_GLOBAL",  "DST_MBD_EPD",  "DST_TRUTH_JET"})
    {
        std::string infile = Form( "%s_pythia8_Jet%d_sHijing_0_20fm-%010d-%06d.root", DSTTPYE, jet_flag, run_number, segment );
        std::cout << "\tAdding input file: " << infile << std::endl;
        auto input = new Fun4AllDstInputManager( Form( "DSTINPUT_%s", DSTTPYE ) );
        input -> AddFile( infile );
        input -> Verbosity( Enable::VERBOSITY );
        se -> registerInputManager( input );
    }

    auto * rcemc = new RetowerCEMC( "RetowerCEMC" );
    rcemc -> set_towerinfo( true );
    rcemc -> set_frac_cut( 1.0 );
    rcemc -> set_do_rescale( false );
    rcemc -> set_towerNodePrefix( HIJETS::tower_prefix );
    rcemc -> Verbosity( Enable::VERBOSITY );
    se -> registerSubsystem( rcemc );

    auto * mb = new MinimumBiasClassifier(  "MinimumBiasClassifier" );
    mb -> setIsSim( true );
    mb -> setOverwriteScale( "/sphenix/user/dlis/Projects/centrality/cdb/calibrations/scales/cdb_centrality_scale_1.root" );
    mb -> setOverwriteVtx( "/sphenix/user/dlis/Projects/centrality/cdb/calibrations/vertexscales/cdb_centrality_vertex_scale_1.root" );
    mb -> Verbosity( Enable::VERBOSITY );
    se -> registerSubsystem( mb );

    auto * cr = new CentralityReco( "CentralityReco" );
    cr -> setOverwriteScale( "/sphenix/user/dlis/Projects/centrality/cdb/calibrations/scales/cdb_centrality_scale_1.root" );
    cr -> setOverwriteVtx( "/sphenix/user/dlis/Projects/centrality/cdb/calibrations/vertexscales/cdb_centrality_vertex_scale_1.root" );
    cr -> setOverwriteDivs( "/sphenix/user/dlis/Projects/centrality/cdb/calibrations/divs/cdb_centrality_1.root" );
    cr -> Verbosity( Enable::VERBOSITY );
    se -> registerSubsystem( cr );

    auto * es = new EventSelector( "EventSelector" );
    es -> Verbosity( Enable::VERBOSITY );
    auto * mbc = new MinBiasCut();
    mbc -> SetNodeName( "MinimumBiasInfo" );
    es -> AddCut( mbc );
    auto * zvc = new ZVertexCut( 60.0, -60.0 );
    zvc -> SetNodeName( "GlobalVertexMap" );
    es -> AddCut( zvc );
    es -> PrintCuts();
    se -> registerSubsystem( es );

    const std::string sub1_seed_raw = "AntiKt_TowerInfo_HIRecoSeedsRaw_r02";
    const std::string sub1_seed_sub = "AntiKt_TowerInfo_HIRecoSeedsSub_r02";
    const std::string rho_kt_seed = "Kt_TowerInfo_HIRecoSeedsRaw_r04";
    std::vector< float > jetRs = {  0.3 };
    // akT seeds and raw aKT jets
    auto * tjr = new JetReco( "TowerJetReco_Raw_akT" );
    for ( const auto & src : { Jet::CEMC_TOWERINFO_RETOWER , Jet::HCALIN_TOWERINFO, Jet::HCALOUT_TOWERINFO } )
    {  
        tjr -> add_input( GetTowerInput( src , HIJETS::tower_prefix ) );
    }
    tjr -> add_algo( HIJETS::GetFJAlgo( 0.2 ), sub1_seed_raw );
    for ( const auto & R : jetRs )
    {
        tjr -> add_algo( HIJETS::GetFJAlgo( R ), Form( "AntiKt_TowerInfo_r0%d", static_cast<int>( R * 10 ) ) );
    }
    tjr -> set_algo_node( "ANTIKT" );
    tjr -> set_input_node( "TOWER" );
    tjr -> Verbosity( 0 );
    se -> registerSubsystem( tjr );

    tjr = new JetReco( "TowerJetReco_Raw_kT" );
    for ( const auto & src : { Jet::CEMC_TOWERINFO_RETOWER , Jet::HCALIN_TOWERINFO, Jet::HCALOUT_TOWERINFO } )
    {   
        tjr -> add_input( GetTowerInput( src , HIJETS::tower_prefix ) );
    }
    FastJetOptions kt_fj_opts({Jet::KT, JET_R, 0.4, VERBOSITY, static_cast<float>(Enable::HIJETS_VERBOSITY)});  
    tjr -> add_algo( new FastJetAlgoSub(kt_fj_opts), rho_kt_seed );
    tjr -> set_algo_node( "KT" );
    tjr -> set_input_node( "TOWER" );
    tjr -> Verbosity( 0 );
    se -> registerSubsystem( tjr );

    // now rho
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

    
    // rho eta-shape calibration for this run, or the dataset default
    // (calibrations/rho_eta/README.md)
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
    for ( const auto & src : { Jet::CEMC_TOWERINFO_RETOWER , Jet::HCALIN_TOWERINFO, Jet::HCALOUT_TOWERINFO } )
    {
        tjr -> add_input( GetTowerInput( src , "MULTSUB_TOWERINFO_CALIB" ) );
    }
    for ( const auto & R : { 0.3 } )
    {
        tjr -> add_algo( HIJETS::GetFJAlgo( R ), Form( "AntiKt_TowerInfo_r0%d_Rho1", static_cast<int>( R * 10 ) ) );
    }
    tjr -> set_algo_node( "ANTIKT" );
    tjr -> set_input_node( "TOWER" );
    // tjr -> Verbosity( 10 );
    se -> registerSubsystem( tjr );

    auto * dtb = new DetermineTowerBackgroundv1( "DetermineTowerBackground_Sub1" );
    dtb -> SetBackgroundOutputName( "TowerInfoBackground_Rho" );
    dtb -> SetVertexType( GlobalVertex::MBD );
    dtb -> SetFlowMode( 0 );
    dtb -> SetPsi2Mode( 0 );
    dtb -> SetNOmitSeeds( 2 );
    dtb -> SetEtaCalib_DirectPath(rho_eta_calib_path);
    dtb -> SetSeedJetName( rho_kt_seed );
    dtb -> SetCEMC_RhoNode("TowerRho_MULT_CEMC");
    dtb -> SetIHCAL_RhoNode("TowerRho_MULT_HCALIN");
    dtb -> SetOHCAL_RhoNode("TowerRho_MULT_HCALOUT");
    dtb -> Verbosity( 0 );
    se -> registerSubsystem( dtb );
    
    auto * subtower = new SubtractTowers( "SubtractTowers" );
    subtower -> set_inputTowerBackgroundNode( "TowerInfoBackground_Rho" );
    subtower -> set_towerNodePrefix( "TOWERINFO_CALIB" );
    subtower -> SetFlowModulation( 0 );
    subtower -> Verbosity( 0 );
    subtower -> set_towerinfo( true );
    se -> registerSubsystem( subtower );

    tjr = new JetReco( "TowerJetReco_Sub1" );
    for (const auto & src : { Jet::CEMC_TOWERINFO_SUB1, Jet::HCALIN_TOWERINFO_SUB1, Jet::HCALOUT_TOWERINFO_SUB1 })
    {
        tjr -> add_input( GetTowerInput( src , "TOWERINFO_CALIB" ) );
    }
    for ( const auto & R : { 0.3 } )
    {
        tjr -> add_algo( HIJETS::GetFJAlgo( R ), Form( "AntiKt_TowerInfo_r0%d_Rho2", static_cast<int>( R * 10 ) ) );
    }
    tjr -> set_algo_node( "ANTIKT" );
    tjr -> set_input_node( "TOWER" );
    // tjr -> Verbosity( 10 );
    se -> registerSubsystem( tjr );

    auto * anaout = new AnaTreev1( outfile );
    anaout -> Verbosity( Enable::VERBOSITY );
    anaout -> add_zvrtx_node( "GlobalVertexMap" );
    anaout -> add_cent_node( "CentralityInfo" );
    anaout -> add_mbd_node( "MbdOut" );
    anaout -> add_minbias_node( "MinimumBiasInfo" );
    anaout -> add_event_header( "EventHeader" );

    if ( jetId > 0 ) 
    {
        // anaout -> add_phHep_node( "PHHepMCGenEventMap" );
        anaout -> add_truth_jet_node( "AntiKt_Truth_r03" );
    }
    
    anaout -> add_sub1_jet_node( "AntiKt_TowerInfo_r03_Rho2", "TowerInfoBackground_Rho" );
    anaout -> add_towerbkgd_v2_node( "TowerInfoBackground_Rho" );

    anaout -> add_rho_jet_node( "AntiKt_TowerInfo_r03_Rho1" , "TowerRho_MULT_CEMC", "TowerRho_MULT_HCALIN", "TowerRho_MULT_HCALOUT" );
    for ( const auto & rho_node : { "TowerRho_MULT_CEMC", "TowerRho_MULT_HCALIN", "TowerRho_MULT_HCALOUT" } )
    {
        anaout -> add_rho_nodes( rho_node );
    }
    anaout -> add_jet_node( rho_kt_seed );
    anaout -> add_cemc_node( "MULTSUB_TOWERINFO_CALIB_CEMC_RETOWER" , true );
    anaout -> add_ihcal_node( "MULTSUB_TOWERINFO_CALIB_HCALIN" , true );
    anaout -> add_ohcal_node( "MULTSUB_TOWERINFO_CALIB_HCALOUT" , true );
   
    
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