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

#include <globalvertex/GlobalVertexReco.h>

#include <centrality/CentralityReco.h>

#include <calotrigger/MinimumBiasClassifier.h>
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

// For runs with a genuine gap in their own MBD calibration (confirmed via
// macros/data/check_mbd_cdb.C -- MBD_QFIT/MBD_T0CORR resolve to nothing at
// that run's timestamp, which starves MinimumBiasClassifier and makes
// MinBiasCut reject every event): re-downloads just those two calibration
// components, on MbdReco's own MbdCalib instance, from a neighboring
// run's already-resolved file. This never touches CDBInterface -- the two
// file paths are looked up once via a transient TIMESTAMP flip that is
// restored immediately, and Download_Gains()/Download_T0Corr() take a
// literal file path, not a CDB domain lookup. Every other subsystem's
// calibration (which for these runs already resolves fine per-run) is
// completely unaffected.
//
// Runs after MbdReco's own InitRun in the same se->run() call (Fun4All
// calls InitRun on registered subsystems in registration order), so this
// must be registered immediately after MbdReco, before se->run() is
// called -- see registration site below.
class MbdCalibOverride : public SubsysReco
{
 public:
    MbdCalibOverride(
        MbdReco * mbdreco,
        const std::string & qfit_url,
        const std::string & t0corr_url
    )
        : SubsysReco( "MbdCalibOverride" )
        , m_mbdreco( mbdreco )
        , m_qfit_url( qfit_url )
        , m_t0corr_url( t0corr_url )
    {}

    int InitRun( PHCompositeNode * ) override
    {
        auto * mbdevent = m_mbdreco -> GetMbdEvent();
        auto * calib = mbdevent ? mbdevent -> GetCalib() : nullptr;
        if ( !calib )
        {
            std::cout << "MbdCalibOverride: MbdReco has no MbdCalib yet -- "
                       << "did registration order change? Nothing overridden." << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }

        std::cout << "MbdCalibOverride: overriding MBD_QFIT -> " << m_qfit_url << std::endl;
        std::cout << "MbdCalibOverride: overriding MBD_T0CORR -> " << m_t0corr_url << std::endl;
        calib -> Download_Gains( m_qfit_url );
        calib -> Download_T0Corr( m_t0corr_url );

        return Fun4AllReturnCodes::EVENT_OK;
    }

 private:
    MbdReco * m_mbdreco{ nullptr };
    std::string m_qfit_url;
    std::string m_t0corr_url;
};

void Fun4All_Dijets_AuAu (
    const int nEvents                    = 100,
    const std::string & infile_calo      = "DST_CALOFITTING_run2auau_pro001_pcdb001_v001-00054912-00000.root",
    const std::string & infile_zdc       = "DST_ZDC_RAW_run2auau_pro001_pcdb001_v001-00054912-00000.root",
    const std::string & infile_sepd      = "DST_SEPD_RAW_run2auau_pro001_pcdb001_v001-00054912-00000.root",
    const std::string & outfile_ana      = "output_ana.root",
    const std::string & cdbtag           = "newcdbtag",
    // known-good run to borrow MBD_QFIT/MBD_T0CORR from (see
    // MbdCalibOverride above) when run_number's own MBD calibration has
    // the gap described above. -1 = disabled (default; MbdReco uses its
    // normal per-run CDB calibration, same as before this option existed).
    const int neighbor_mbd_run           = -1
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

    if ( neighbor_mbd_run > 0 )
    {
        // transient TIMESTAMP flip, restored immediately -- CDBInterface
        // never leaves its normal live-lookup mode, so every other
        // subsystem's own per-run calibration is unaffected.
        rc -> set_uint64Flag( "TIMESTAMP", neighbor_mbd_run );
        const std::string qfit_url   = cdb -> getUrl( "MBD_QFIT" );
        const std::string t0corr_url = cdb -> getUrl( "MBD_T0CORR" );
        rc -> set_uint64Flag( "TIMESTAMP", run_number );

        se -> registerSubsystem( new MbdCalibOverride( mbd, qfit_url, t0corr_url ) );
    }

    auto * sepdbuilder = new CaloTowerBuilder( "SEPDBUILDER" );
    sepdbuilder -> set_detector_type( CaloTowerDefs::SEPD );
    sepdbuilder -> set_builder_type( buildertype );
    sepdbuilder -> set_processing_type( CaloWaveformProcessing::TEMPLATE );
    sepdbuilder -> set_nsamples( 12 );
    sepdbuilder -> set_offlineflag( );
    se -> registerSubsystem( sepdbuilder );

    auto * sepd = new EpdReco( );
    se -> registerSubsystem(sepd);

    auto * caZDC = new CaloTowerBuilder("ZDCBUILDER");
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
    auto * tcut = new TriggerSelect();
    tcut -> SetPacket( 14001 );
    tcut -> SelectTrigger( 10 );
    es -> AddCut( tcut );
    auto * mbcut = new MinBiasCut( );
    mbcut -> SetNodeName( "MinimumBiasInfo" );
    es -> AddCut( mbcut );
    es -> PrintCuts( );
    se -> registerSubsystem( es );

    auto * tjr = new JetReco( "TowerJetReco_Emb" );
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

    auto * dtb1 = new DetermineTowerBackground( "DetermineTowerBackground1_Emb" );
    dtb1 -> SetBackgroundOutputName( "TowerInfoBackground_Sub1" );
    dtb1 -> set_towerNodePrefix( "TOWERINFO_CALIB" );
    dtb1 -> SetFlow( 0 );
    dtb1 -> SetSeedJetPt( 5.0 );
    dtb1 -> SetSeedJetD( 3.0 );
    dtb1 -> SetSeedMaxConst( 3.0 );
    dtb1 -> Verbosity( 0 );
    dtb1 -> UseReweighting( true );
    se -> registerSubsystem( dtb1 );

    auto * casj1 = new CopyAndSubtractJets( "CopyAndSubtractJets_Emb" );
    casj1 -> set_towerNodePrefix( "TOWERINFO_CALIB" );
    casj1 -> SetFlowModulation( 0 );
    casj1 -> Verbosity( 0 );
    casj1 -> set_towerinfo( true );
    se -> registerSubsystem( casj1 );

    auto * dtb2 = new DetermineTowerBackground( "DetermineTowerBackground_Emb_Sub1" );
    dtb2 -> SetBackgroundOutputName( "TowerInfoBackground_Sub2" );
    dtb2 -> set_towerNodePrefix( "TOWERINFO_CALIB" );
    dtb2 -> SetFlow( 0 );
    dtb2 -> SetSeedJetPt( 7.0 );
    dtb2 -> Verbosity( 0 );
    dtb2 -> UseReweighting( true );
    se -> registerSubsystem( dtb2 );

    auto * st1 = new SubtractTowers( "SubtractTowers_Emb" );
    st1 -> set_towerNodePrefix( "TOWERINFO_CALIB" );
    st1 -> SetFlowModulation( 0 );
    st1 -> Verbosity( 0 );
    st1 -> set_towerinfo( true );
    se -> registerSubsystem( st1 );

    tjr = new JetReco( "TowerJetReco_Emb_Sub1" );
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
    
    
    auto * trc = new DetermineTowerRho( "DetermineTowerRho_CEMC_MultEmb" );
    trc -> add_method( TowerRho::Method::MULT, "TowerRho_MULT_CEMC" );
    trc -> add_tower_input( GetTowerInput( Jet::CEMC_TOWERINFO_RETOWER ) );
    se -> registerSubsystem( trc );

    trc = new DetermineTowerRho( "DetermineTowerRho_HCALIN_MultEmb" );
    trc -> add_method( TowerRho::Method::MULT, "TowerRho_MULT_HCALIN" );
    trc -> add_tower_input( GetTowerInput( Jet::HCALIN_TOWERINFO ) );
    se -> registerSubsystem( trc );  
    
    trc = new DetermineTowerRho( "DetermineTowerRho_HCALOUT_MultEmb" );
    trc -> add_method( TowerRho::Method::MULT, "TowerRho_MULT_HCALOUT" );
    trc -> add_tower_input( GetTowerInput( Jet::HCALOUT_TOWERINFO ) );
    se -> registerSubsystem( trc );

    // rho eta-shape calibration for this run, or the dataset default
    // (calibrations/rho_eta/README.md)
    const std::string rho_eta_calib_path = RhoEtaCalibLookup::GetCalibPath( run_number );

    auto * subrho = new SubtractTowersRhov1(  "SubtractTowersRho_CEMC_MultEmb" );
    subrho -> set_rhoNode("TowerRho_MULT_CEMC");
    subrho -> add_targetTowerNode( "TOWERINFO_CALIB_CEMC_RETOWER" );
    subrho -> set_etaCalib_directPath( rho_eta_calib_path );
    subrho -> set_subSuffix( "MULTSUB" );
    subrho -> set_globalVertexType( GlobalVertex::MBD );
    se -> registerSubsystem( subrho );

    subrho = new SubtractTowersRhov1( "SubtractTowersRho_HCALIN_MultEmb" );
    subrho -> set_rhoNode("TowerRho_MULT_HCALIN");
    subrho -> add_targetTowerNode( "TOWERINFO_CALIB_HCALIN" );
    subrho -> set_etaCalib_directPath( rho_eta_calib_path );
    subrho -> set_subSuffix( "MULTSUB" );
    subrho -> set_globalVertexType( GlobalVertex::MBD );
    se -> registerSubsystem( subrho );

    subrho = new SubtractTowersRhov1( "SubtractTowersRho_HCALOUT_MultEmb" );
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
        tjr -> add_algo( HIJETS::GetFJAlgo( R ), Form( "AntiKt_TowerInfo_r0%d_Rho2", static_cast<int>( R * 10 ) ) );
    }
    tjr -> set_algo_node( "ANTIKT" );
    tjr -> set_input_node( "TOWER" );
    tjr -> Verbosity( 0 );
    se -> registerSubsystem( tjr );

    auto * anaout = new AnaTreev1( outfile_ana );
    anaout -> Verbosity( Enable::VERBOSITY );
    anaout -> add_zvrtx_node( "GlobalVertexMap" );
    anaout -> add_cent_node( "CentralityInfo" );
    anaout -> add_gl1_node( "14001" );
    anaout -> add_mbd_node( "MbdOut" );
    anaout -> add_minbias_node( "MinimumBiasInfo" );
    anaout -> add_sub1_jet_node( "AntiKt_TowerInfo_r03_Sub1", "TowerInfoBackground_Sub2" );
    anaout -> add_towerbkgd_v2_node( "TowerInfoBackground_Sub2" );
    anaout -> add_rho_jet_node( "AntiKt_TowerInfo_r03_Rho2" , "TowerRho_MULT_CEMC", "TowerRho_MULT_HCALIN", "TowerRho_MULT_HCALOUT" );
    for ( const auto & rho_node : { "TowerRho_MULT_CEMC", "TowerRho_MULT_HCALIN", "TowerRho_MULT_HCALOUT" } )
    {
        anaout -> add_rho_nodes( rho_node );
    }
    anaout -> add_jet_node( "AntiKt_TowerInfo_r03" );
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

void LoadNeighborMbdCalib(
    const int run_number,
    const int neighbor_run,
    const std::string & cdbtag
)
{
    // domains confirmed missing at the affected runs' own timestamp by
    // macros/data/check_mbd_cdb.C -- everything else (MBD_SAMPMAX,
    // MBD_SLEWCORR, MBD_TIMECORR, MBD_TQ_T0, MBD_TT_T0, and every non-MBD
    // domain) already resolves fine for these runs and must keep coming
    // from run_number itself, not the neighbor.
    static const std::vector<std::string> mbd_gap_domains = { "MBD_QFIT", "MBD_T0CORR" };

    auto * rc  = recoConsts::instance();
    auto * cdb = CDBInterface::instance();

    rc -> set_StringFlag( "CDB_GLOBALTAG", cdbtag );

    const std::string owndump  = Form( "/tmp/cdb_own_%d.txt", run_number );
    const std::string neighdump = Form( "/tmp/cdb_neighbor_%d_for_%d.txt", neighbor_run, run_number );
    const std::string hybrid    = Form( "/tmp/cdb_hybrid_%d.txt", run_number );

    // run_number's own full resolution (missing the gap domains, present
    // for everything else)
    rc -> set_uint64Flag( "TIMESTAMP", run_number );
    cdb -> getUrl( "calo_geo" );
    cdb -> DumpCalibrations( owndump );

    // neighbor_run's full resolution, used only to source the gap domains
    rc -> set_uint64Flag( "TIMESTAMP", neighbor_run );
    cdb -> getUrl( "calo_geo" );
    cdb -> DumpCalibrations( neighdump );

    std::ofstream out( hybrid );
    std::ifstream own_in( owndump );
    std::string line;
    while ( std::getline( own_in, line ) ) out << line << "\n";
    own_in.close();

    std::ifstream neigh_in( neighdump );
    while ( std::getline( neigh_in, line ) )
    {
        for ( const auto & domain : mbd_gap_domains )
        {
            if ( line.rfind( domain + " ", 0 ) == 0 ) out << line << "\n";
        }
    }
    neigh_in.close();
    out.close();

    rc -> set_uint64Flag( "TIMESTAMP", run_number );
    cdb -> ReadCalibrationsFromFile( hybrid );

    std::cout << "LoadNeighborMbdCalib: run " << run_number << " borrowing ";
    for ( const auto & domain : mbd_gap_domains ) std::cout << domain << " ";
    std::cout << "from neighbor run " << neighbor_run
               << "; all other calibrations remain run " << run_number << "'s own." << std::endl;

    std::remove( owndump.c_str() );
    std::remove( neighdump.c_str() );
    std::remove( hybrid.c_str() );
}

#endif