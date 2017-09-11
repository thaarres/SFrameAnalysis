
// Local include(s):
#include "../include/VVanalysis.h"
#include "../include/VVAnalysisTools.h"

// External include(s):
// #include "../GoodRunsLists/include/TGoodRunsListReader.h"


#include <TMath.h>
#include "TGraph.h"
#include <iostream>
#include <fstream>

ClassImp( VVanalysis );

VVanalysis::VVanalysis()
  : SCycleBase()
   , m_jetAK4               ( this )
   , m_jetAK8               ( this )
   , m_jetAK8Puppi          ( this )
   , m_eventInfo            ( this )
   , m_electron             ( this )
   , m_muon                 ( this )
   , m_missingEt            ( this )
   , m_genParticle          ( this )
   , m_pileupReweightingTool( this, "pileupReweighting" )
   , m_bTaggingScaleTool    ( this )
   , m_xSec                 ()
   , m_allEvents( "allEvents", this )
   , m_foundLepton( "passedLoose", this )
   , m_passedPuppi( "passedPuppi", this )
   , m_passedAK4( "passedDeta", this )
   , m_passedMET( "passedMjj", this )
   , m_passedEvents( "passedEvents", this )
   , m_test( "test", this ) 
   , m_genjetAK8( this )
 {
   m_logger << INFO << "Initializing" << SLogger::endmsg;
   SetLogName( GetName() );
   
   DeclareProperty( "JetAK8Name"      ,         m_jetAK8Name              = "jetAK8"        );
   DeclareProperty( "JetAK4Name"      ,         m_jetAK4Name              = "jetAK4"        );
   DeclareProperty( "GenJetAK8Name"   ,         m_genjetAK8Name           = "genJetAK8"     );
   DeclareProperty( "JetAK8PuppiName" ,         m_jetAK8PuppiName         = "jetAK8_puppi"  );
   DeclareProperty( "GenParticleName" ,         m_genParticleName         = "genParticle"   );
   DeclareProperty( "ElectronName"    ,         m_electronName            = "el"            );
   DeclareProperty( "MuonName"        ,         m_muonName                = "mu"            );
   DeclareProperty( "MissingEtName",            m_missingEtName            = "MET" );
   DeclareProperty( "RecoTreeName"    ,         m_recoTreeName            = "tree" );
   DeclareProperty( "PUPPIJEC"        ,         m_PUPPIJEC                = "weights/puppiCorr.root" );
   // DeclareProperty( "JSONName"        ,         m_jsonName                = std::string (std::getenv("SFRAME_DIR")) + "/../GoodRunsLists/JSON/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.txt" );
   DeclareProperty( "IsData"          ,         m_isData                  = false );
   DeclareProperty( "Channel"         ,         m_Channel                 = "el" );
   
   
   nSumGenWeights           = 0;
   
}

VVanalysis::~VVanalysis() {
  
  m_logger << INFO << "Done!" << SLogger::endmsg;

}

void VVanalysis::BeginCycle() throw( SError ) { //Block used to perform an initial configuration of the cycle (ead some local file likegood data ranges)
  
  // // Load the GRL:
  // if (m_isData) {
  //   m_logger << INFO << "Loading GoodRunsList from file: " << m_jsonName << SLogger::endmsg;
  //   Root::TGoodRunsListReader reader( TString( m_jsonName.c_str() ) );
  //   if( ! reader.Interpret() ) {
  //     m_logger << FATAL << "Couldn't read in the GRL!" << SLogger::endmsg;
  //     throw SError( ( "Couldn't read in file: " + m_jsonName ).c_str(), SError::SkipCycle );
  //   }
  //   m_grl = reader.GetMergedGoodRunsList();
  //   m_grl.Summary();
  //   m_grl.SetName( "MyGoodRunsList" );
  //
  //   // Add it as a configuration object, so that the worker nodes will receive it:
  //   AddConfigObject( &m_grl );
  // }

   return;

}

void VVanalysis::EndCycle() throw( SError ) { //Any finalisation steps should be done here. (Closure of some helper files opened by the user code for instance.)

   return;

}

void VVanalysis::BeginInputFile( const SInputData& ) throw( SError ) { //For each new input file the user has to connect his input variables. This has to be performed in this function.
  
  m_logger << INFO << "Connecting input variables" << SLogger::endmsg;
  if (m_isData) {
    
    m_jetAK8Puppi .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetSubstructure, (m_jetAK8PuppiName + "_").c_str() );
    m_jetAK8      .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetSubstructure|Ntuple::JetSoftdropSubjets, (m_jetAK8Name + "_").c_str() );
    m_jetAK4      .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetTruth, (m_jetAK4Name + "_").c_str() );
    m_eventInfo   .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::EventInfoBasic|Ntuple::EventInfoTrigger|Ntuple::EventInfoMETFilters, "" );
    m_jetAK4      .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetTruth, (m_jetAK4Name + "_").c_str() );
    m_electron    .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::ElectronBasic|Ntuple::ElectronID, (m_electronName + "_").c_str() );
    m_muon        .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::MuonBasic|Ntuple::MuonID|Ntuple::MuonIsolation, (m_muonName + "_").c_str() );
    m_missingEt    .ConnectVariables( m_recoTreeName.c_str(), Ntuple::MissingEtBasic, (m_missingEtName + "_").c_str() );
  
  }
  else {
    m_jetAK8Puppi .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetSubstructure, (m_jetAK8PuppiName + "_").c_str() );
    m_jetAK8      .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetSubstructure|Ntuple::JetTruth|Ntuple::JetSoftdropSubjets|Ntuple::JetSoftdropSubjetsTruth, (m_jetAK8Name + "_").c_str() );
    m_jetAK4      .ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetTruth, (m_jetAK4Name + "_").c_str() );
    m_eventInfo   .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::EventInfoBasic|Ntuple::EventInfoTrigger|Ntuple::EventInfoMETFilters|Ntuple::EventInfoTruth, "" );
    m_genParticle .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::GenParticleBasic, (m_genParticleName + "_").c_str() );
    m_genjetAK8   .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::GenJet, (m_genjetAK8Name + "_").c_str() );
    m_electron    .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::ElectronBasic|Ntuple::ElectronID|Ntuple::ElectronSuperCluster, (m_electronName + "_").c_str() );
    m_muon        .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::MuonBasic|Ntuple::MuonID|Ntuple::MuonIsolation, (m_muonName + "_").c_str() );
    m_missingEt   .ConnectVariables(    m_recoTreeName.c_str(), Ntuple::MissingEtBasic, (m_missingEtName + "_").c_str() );
  }
  // Unused collections
  // m_jetAK4.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetTruth, (m_jetAK4Name + "_").c_str() );
  // m_missingEt.ConnectVariables(    m_recoTreeName.c_str(), Ntuple::MissingEtBasic, (m_missingEtName + "_").c_str() );
  
   
  m_logger << INFO << "Connecting input variables completed" << SLogger::endmsg;
  
  m_logger << INFO << "Initializing inputs for puppi softdop mass corrections" << SLogger::endmsg; 
  TFile* file = TFile::Open( m_PUPPIJEC.c_str(),"READ");
  if(file->IsZombie()) throw SError(SError::StopExecution);
  m_puppisd_corr.push_back( dynamic_cast<TF1*>(file->Get("puppiJECcorr_gen")));
  m_puppisd_corr.push_back( dynamic_cast<TF1*>(file->Get("puppiJECcorr_reco_0eta1v3")));
  m_puppisd_corr.push_back( dynamic_cast<TF1*>(file->Get("puppiJECcorr_reco_1v3eta2v5")));
  
  // Get cross section
  b_xSec                        = 1.;
  TString infile = TString(this->GetHistInputFile()->GetName());
  b_xSec  = m_xSec.getLumiWeight( infile ); 
  m_logger << INFO << "Cross section set to " << b_xSec << " for file " << infile <<  SLogger::endmsg;
  if( b_xSec < 0.) throw SError( SError::SkipFile );

  
  return;

}


void VVanalysis::BeginInputData( const SInputData& id ) throw( SError ) { //called once before processing each of the input data types. If you need to initialise output objects (histograms, etc.) before the event-by-event execution, you should do that here. Also the declaration of the output variables has to be done here.
  
  if (!m_isData) m_pileupReweightingTool.BeginInputData( id );
  

  // Declare output variables:
  DeclareVariable( b_lumi                   , "lumi"  );
  DeclareVariable( b_event                  , "evt"  );
  DeclareVariable( b_run                    , "run"  );
  DeclareVariable( b_weight                 , "weight"  );
  DeclareVariable( b_weightGen              , "genWeight");
  DeclareVariable( b_weightPU               , "puWeight");
  DeclareVariable( b_xSec                   , "xsec");
  DeclareVariable( m_o_nJ                   , "njj");
  DeclareVariable( m_o_mpuppisoftdrop       , "Whadr_softDrop_mass");
  DeclareVariable( m_o_mgensoftdrop         , "Whadr_gen_softDrop_mass");
  DeclareVariable( m_o_tau1                 , "Whadr_tau1");  
  DeclareVariable( m_o_tau2                 , "Whadr_tau2");  
  DeclareVariable( m_o_tau2                 , "Whadr_tau3");  
  DeclareVariable( m_o_tau21                , "Whadr_tau21"); 
  DeclareVariable( m_o_pt                   , "Whadr_pt");
  DeclareVariable( m_o_eta                  , "Whadr_eta");
  DeclareVariable( m_o_genpt                , "Whadr_gen_pt");
  DeclareVariable( m_o_csv                  , "Whadr_csv");
  DeclareVariable( m_o_lep_pt               , "lep_pt");
  DeclareVariable( m_o_Wlep_pt              , "Wlep_pt");
  DeclareVariable( m_o_lep_eta              , "lep_eta");
  DeclareVariable( m_o_lep_phi              , "lep_phi");
  
  DeclareVariable( Flag_goodVertices              , "Flag_goodVertices");
  DeclareVariable( Flag_globalTightHalo2016Filter , "Flag_CSCTightHaloFilter");
  DeclareVariable( Flag_HBHENoiseFilter           , "Flag_HBHENoiseFilter");
  DeclareVariable( Flag_HBHENoiseIsoFilter        , "Flag_HBHENoiseIsoFilter");
  DeclareVariable( Flag_eeBadScFilter             , "Flag_eeBadScFilter");
  DeclareVariable( Flag_badChargedHadronFilter    , "Flag_badChargedHadronFilter");
  DeclareVariable( Flag_badMuonFilter             , "Flag_badMuonFilter");
  DeclareVariable( Flag_ECALDeadCell              , "Flag_ECALDeadCell");
  DeclareVariable( HLTMu50           , "HLTMu50");
  DeclareVariable( HLTMu45_eta2p1            , "HLTMu45_eta2p1");
  DeclareVariable( HLTEle105_CaloIdVT_GsfTrkIdT       , "HLTEle105_CaloIdVT_GsfTrkIdT");
  DeclareVariable( HLTEle115_CaloIdVT_GsfTrkIdT       , "HLTEle115_CaloIdVT_GsfTrkIdT");
  DeclareVariable( HLT_all                         , "HLT_all");
  DeclareVariable( m_o_nLeptons           , "nOtherLeptons");
  DeclareVariable( mergedVTruth                , "mergedVTruth");

  // // Declare the output histograms:
  // Book( TH1F( "Mhist"     , "Mjj", 100, 0.0, 3000. ) );
  
 
  m_test->resize( 6, 0 ); // Reserve six entries in vector for counting entries:
  TH1* genEvents = Book(TH1F("genEvents", "number of generated Events",2,0,2));
              
   return;
}

void VVanalysis::EndInputData( const SInputData& ) throw( SError ) {
  
  // static const Int_t n = 5;
 //  Float_t x_array[ n ] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
 //  Float_t y_array[ n ] = { 0.0, 2.0, 4.0, 6.0, 8.0 };
 //  TGraph mygraph( n, x_array, y_array );
 //  mygraph.SetName( "MyGraph" );
 //  WriteObj( mygraph, "graph_dir" );
    
//   std::ofstream myfile;
//   myfile.open ("debug.txt");
//   std::sort (debugEvent.begin(), debugEvent.end());
//   for(int i=0;i< debugEvent.size();i++){
//     myfile << debugEvent.at(i) << "\n";
//   }
//   myfile.close();
   return;

}



void VVanalysis::BeginMasterInputData( const SInputData& ) throw( SError ){

   return;

}

void VVanalysis::EndMasterInputData( const SInputData& ) throw( SError ){ //this is a good place to print some summaries, do some final calculations on the created histograms (for instance fitting them), etc
    
  m_logger << INFO << "Number of all processed events      :  "<< *m_allEvents   << "   " << ( m_test->size() > 0 ? ( *m_test )[ 0 ] : 0 )<< SLogger::endmsg;
  m_logger << INFO << "Number of events passing lepton req.:  "<< *m_foundLepton << "   " << ( m_test->size() > 1 ? ( *m_test )[ 1 ] : 0 )<< SLogger::endmsg;
  m_logger << INFO << "Found AK8 PUPPI jet                 :  "<< *m_passedPuppi << "   " << ( m_test->size() > 2 ? ( *m_test )[ 2 ] : 0 )<< SLogger::endmsg;
  m_logger << INFO << "Found at least one b-tagged AK4     :  "<< *m_passedAK4  << "   " << ( m_test->size() > 4 ? ( *m_test )[ 4 ] : 0 )<< SLogger::endmsg;
  m_logger << INFO << "Number of events passing MET cut    :  "<< *m_passedMET   << "   " << ( m_test->size() > 3 ? ( *m_test )[ 3 ] : 0 )<< SLogger::endmsg;
  m_logger << INFO << "Number of events passing selection  :  "<< *m_passedEvents<< "   " << ( m_test->size() > 5 ? ( *m_test )[ 5 ] : 0 )<< SLogger::endmsg;
  
  m_logger << INFO << "Number of generated Events (weighted) : "<< nSumGenWeights << SLogger::endmsg;
  
   return;

}
  
  


void VVanalysis::ExecuteEvent( const SInputData&, Double_t weight) throw( SError ) { //This is the main analysis function that is called for each event. It receives the weight of the event, as it is calculated by the framework from the luminosities and generator cuts defined in the XML configuration.
  
  // float gW = (m_eventInfo.genEventWeight < 0) ? -1 : 1;
  float gW = m_eventInfo.genEventWeight;  //needed for Herwig
  nSumGenWeights += gW;
  Hist( "genEvents" )->Fill(gW);
  
  clearBranches();
  
  if (!m_isData) {
    b_weight= getEventWeight();
  }
  // Only take generated hadronic events
  // if( m_isSignal && !SignalIsHad( m_genParticle, m_Channel)) throw SError( SError::SkipEvent);

  ++m_allEvents;
  ( *m_test )[ 0 ]++;
  if( m_jetAK8.N < 1 ) throw SError( SError::SkipEvent );
  //-------------Find lepton-------------//
  // electrons: HEEP - pt>120 GeV - |η| <1.442 or 1.56<|η| <2.5
  // muons: HighPt ID - pt>53 GeV - |η| <2.1 - relIso03<0.1
  std::vector<UZH::Electron> goodElectrons = FindGoodLeptons(m_electron);
  std::vector<UZH::Muon>     goodMuons     = FindGoodLeptons(m_muon);
  
  if( (m_Channel.find("el")!= std::string::npos && goodElectrons.size()<1) || (m_Channel.find("mu")!= std::string::npos && goodMuons.size()<1) ) throw SError( SError::SkipEvent);
  ++m_foundLepton;
  ( *m_test )[ 1 ]++;
  //-------------Find good AK8-------------//
  // ak08: Loose ID - pt>100 GeV - |η| <2.4 - ΔR(lepton)>1
  std::vector<UZH::Jet> goodFatJets;
  std::vector<UZH::Jet> goodGenJets;
  std::vector<int> puppiMatch;
  for ( int i = 0; i < (m_jetAK8.N); ++i ) {
    UZH::Jet myjet( &m_jetAK8, i );
    if ( i == 0 && !myjet.IDLoose()) break;
    if (! (myjet.pt() > 200       )) continue;
    if (! (fabs(myjet.eta()) < 2.5)) continue;
    if (! (myjet.IDTight()        )) continue;
    //Match to puppi jet
    float dRmin = 99.;
    if(m_jetAK8Puppi.pt->size() < 1) throw SError( SError::SkipEvent);
    for ( int ii = 0; ii < abs((m_jetAK8Puppi.pt->size())); ++ii ) {
      UZH::Jet mypuppijet( &m_jetAK8Puppi, ii ); 
      float dR = myjet.DeltaR(mypuppijet);
      if ( dR > dRmin ) continue;
      bool samePuppiJet =0;
      for(int m=0;m<abs(puppiMatch.size());m++)
      {
        if(ii == puppiMatch.at(m)) samePuppiJet=1;
      }
      if (samePuppiJet) continue;
      puppiMatch.push_back(ii);
      dRmin = dR;   
      myjet.puppi_softdropmass= ApplyPuppiSoftdropMassCorrections(mypuppijet,m_puppisd_corr,m_isData);//mypuppijet.softdrop_mass();
      myjet.puppi_tau1        = mypuppijet.tau1();
      myjet.puppi_tau2        = mypuppijet.tau2();
      myjet.puppi_tau3        = mypuppijet.tau3();
    }
    if(! FoundNoLeptonOverlap(goodElectrons,goodMuons,myjet.tlv(), 1.0 ) ) continue;
    goodFatJets.push_back(myjet);
  }
  if(goodFatJets.size()<1) throw SError( SError::SkipEvent );
  //-------------Make sure we have a PUPPI match-------------// 
  if( goodFatJets[0].puppi_softdropmass == -99  || goodFatJets[0].puppi_tau1 == -99 || goodFatJets[0].puppi_tau2 == -99) throw SError( SError::SkipEvent );
  ++m_passedPuppi;
  ( *m_test )[ 2 ]++;
 
  //Match to gen jet
  if(!m_isData){
    for( unsigned int i=0; i < goodFatJets.size(); ++i){
      float dRmin = 99.;
      int jetIdx = 99;
      for ( int j = 0; j < (m_genjetAK8.N); ++j ) {
        UZH::Jet genJet( &m_genjetAK8, j );
        if ( genJet.pt() < 50. ) continue;
        float dR = (genJet.tlv()).DeltaR((goodFatJets.at(i)).tlv());
        if ( dR > dRmin ) continue;
        dRmin = dR;
        jetIdx = j;
      }
      UZH::Jet selectedJet( &m_genjetAK8, jetIdx );
      goodGenJets.push_back(selectedJet);
    }
  }
  //-------------Find AK4-------------//
  // ak04: Loose ID - pt>30 GeV - |η| <2.4 - ΔR(lepton)>0.3 - ΔR(ak08)>0.8
  std::vector<UZH::Jet> goodJetsAK4 = FindGoodJetsAK4(m_jetAK4,goodElectrons,goodMuons,goodFatJets[0].tlv());
   if(goodJetsAK4.size()<1) throw SError( SError::SkipEvent );
  ++m_passedAK4;
  ( *m_test )[ 3 ]++;

  
  //-------------Find MET-------------//
  // electron events: met>80 GeV
  // muon events: met>40 GeV
  float m_metCut = 40.0;
  if (m_Channel.find("el")!= std::string::npos) m_metCut==80.0;
  bool foundMet=false;
  UZH::MissingEt goodMet( &m_missingEt, 0 );
  if (goodMet.et() > m_metCut) foundMet=true;
  if (!foundMet) throw SError( SError::SkipEvent);
  ++m_passedMET;
  ( *m_test )[ 4 ]++;
//
//
//
//   //-------------Find leptonic W-------------//
//   //sum of lepton passing selections and met,pt>200 GeV
  TLorentzVector leptonCand_;

  if      (m_Channel.find("el")!= std::string::npos)  leptonCand_ = goodElectrons[0].tlv();
  else if (m_Channel.find("mu")!= std::string::npos)  leptonCand_ = goodMuons[0].tlv();
  else{
    std::cout<< "Channel not defined!! Please type in Channel==mu or Channel=el "<< std::endl;
    throw SError( SError::SkipEvent);
  }

 TLorentzVector p4nu = NuMomentum( leptonCand_.Px(), leptonCand_.Py(), leptonCand_.Pz(), leptonCand_.Pt(), leptonCand_.E(), goodMet.corrPx(), goodMet.corrPy() );

 if( (p4nu+leptonCand_).Pt() < 200. ) throw SError( SError::SkipEvent);
  
// Double_t WmassLep    = TMath::Sqrt( 2*METCand_.at(0).p4.Pt()*leptonCand_.at(0).p4.Pt()*(1-cos(leptonCand_.at(0).p4.DeltaPhi(METCand_.at(0).p4))));


  ++m_passedEvents;
  ( *m_test )[ 5 ]++;
  

  // Fill tree
  m_o_Wlep_pt         = (p4nu+leptonCand_).Pt()   ;
  m_o_nLeptons        = goodElectrons.size()+goodMuons.size();
  m_o_nJ              = goodFatJets.size();
  m_o_mpuppisoftdrop  = goodFatJets[0].puppi_softdropmass;
  m_o_tau1            = goodFatJets[0].puppi_tau1;
  m_o_tau2            = goodFatJets[0].puppi_tau2;
  m_o_tau3            = goodFatJets[0].puppi_tau3;
  m_o_csv             = goodFatJets[0].csv();
  m_o_tau21           = goodFatJets[0].puppi_tau2/goodFatJets[0].puppi_tau1;
  m_o_pt              = goodFatJets[0].tlv().Pt();
  m_o_eta             = goodFatJets[0].tlv().Eta();
  m_o_lep_pt          = leptonCand_.Pt()   ;
  m_o_lep_eta         = leptonCand_.Eta()   ;
  m_o_lep_phi         = leptonCand_.Phi()   ;
 
  
  
  bool isfired = false;
  for( std::map<std::string,bool>::iterator it = (m_eventInfo.trigDecision)->begin(); it != (m_eventInfo.trigDecision)->end(); ++it){
    if( (it->first).find("HLT_Mu50_v")                          != std::string::npos) HLTMu50 = it->second;
    if( (it->first).find("HLT_Mu45_eta2p1_v")                   != std::string::npos) HLTMu45_eta2p1 = it->second;
    if( (it->first).find("HLT_Ele105_CaloIdVT_GsfTrkIdT_v")     != std::string::npos) HLTEle105_CaloIdVT_GsfTrkIdT  = it->second;
    if( (it->first).find("HLT_Ele115_CaloIdVT_GsfTrkIdT_v")     != std::string::npos) HLTEle115_CaloIdVT_GsfTrkIdT  = it->second;
    if      (m_Channel.find("el")!= std::string::npos && (HLTEle105_CaloIdVT_GsfTrkIdT or HLTEle115_CaloIdVT_GsfTrkIdT)) isfired = true;
    else if (m_Channel.find("mu")!= std::string::npos && (HLTMu50 or HLTMu45_eta2p1 )   )                                isfired = true;
  }
  
  b_lumi                          = m_eventInfo.lumiBlock;
  b_event                         = m_eventInfo.eventNumber;
  b_run                           = m_eventInfo.runNumber;
  Flag_goodVertices               = m_eventInfo.PV_filter;
  Flag_globalTightHalo2016Filter  = m_eventInfo.passFilter_CSCHalo;
  Flag_HBHENoiseFilter            = m_eventInfo.passFilter_HBHELoose;
  Flag_HBHENoiseIsoFilter         = m_eventInfo.passFilter_HBHEIso;
  Flag_eeBadScFilter              = m_eventInfo.passFilter_EEBadSc;
  Flag_badChargedHadronFilter     = m_eventInfo.passFilter_chargedHadronTrackResolution;
  Flag_badMuonFilter              = m_eventInfo.passFilter_muonBadTrack;
  Flag_ECALDeadCell               = m_eventInfo.passFilter_ECALDeadCell;
  HLT_all                          = isfired;
  

  if(!m_isData){
    m_o_mgensoftdrop          = goodGenJets[0].softdropmass();
    m_o_genpt                 = goodGenJets[0].tlv().Pt();
  }
  
  std::vector<UZH::GenParticle> GenQuarks  = FindGeneratedQuarks(m_genParticle, m_isData);
  mergedVTruth                =  isMergedVJet(goodFatJets[0].tlv(),GenQuarks) ;
  //

   
  
  // }
  // TString StatusString = "Before_Cuts_";
  // FillValidationHists( GENERAL, StatusString );
  // FillValidationHists( ELECTRON, StatusString );
  // FillValidationHists( MUON, StatusString );
  // FillValidationHists( JETS, StatusString );
 
  return;

}

 

void VVanalysis::FillValidationHists( ValHistsType ht, const TString& status ) {  

   if( ht == GENERAL ){
      TString suffix = status + "General_";
   } 
   else if( ht == ELECTRON ){
      TString suffix = status + "Electron_";
   } 
   else if( ht == MUON ){
      TString suffix = status + "Muon_";
   } 
   else if( ht == JETS ){
      TString suffix = status + "Jets_";
      // Book( TH1F( suffix + "mjj", "mjj", 300, 0.0, 3000.0 ) )->Fill( m_o_mjj );
   } 
   else {
      SError error( SError::StopExecution );
      error << GetName() << ": Validation histogram type"  << ht << " not supported";
      throw error;
   }

   return;
}

double VVanalysis::getEventWeight() {
  
  double weight = 1.;
  for( unsigned int v = 0; v < (m_eventInfo.actualIntPerXing)->size(); ++v ){
    if ( (*m_eventInfo.bunchCrossing)[v] == 0 ) {
      b_weightPU = m_pileupReweightingTool.getPileUpweight( (*m_eventInfo.actualIntPerXing)[v] );
      m_logger << VERBOSE << "Weight: " << b_weightPU << " for true: " << (*m_eventInfo.actualIntPerXing)[v] << SLogger::endmsg;
      break;
    }
  }
  // b_weightGen = (m_eventInfo.genEventWeight < 0) ? -1 : 1;
  b_weightGen = m_eventInfo.genEventWeight; //needed for Herwig
  weight *= b_weightPU*b_weightGen;
  
  return weight;
  
}


void VVanalysis::clearBranches() {
  b_weight                      = 1.;
  b_weightGen                   = 1.;
  b_weightPU                    = 1.;
  b_weightBtag                  = 1.;
  m_o_mpuppisoftdrop            = -99;
  m_o_mgensoftdrop              = -99;
  m_o_tau1                      = -99;
  m_o_tau2                      = -99;
  m_o_tau21                     = -99;
  m_o_lep_pt                    = -99;
  m_o_Wlep_pt                   = -99;
  m_o_lep_eta                   = -99;
  m_o_lep_phi                   = -99;
  
  Flag_goodVertices             = -99;
  Flag_globalTightHalo2016Filter= -99;
  Flag_HBHENoiseFilter          = -99;
  Flag_HBHENoiseIsoFilter       = -99;
  Flag_eeBadScFilter            = -99;
  Flag_badChargedHadronFilter   = -99;
  Flag_badMuonFilter            = -99;
  Flag_ECALDeadCell             = -99;
  HLTMu50          = -99;
  HLTMu45_eta2p1           = -99;
  HLTEle105_CaloIdVT_GsfTrkIdT      = -99;
  HLTEle115_CaloIdVT_GsfTrkIdT      = -99;
  nLeptonOverlap                = -99;
  mergedVTruth               = -99;
}
