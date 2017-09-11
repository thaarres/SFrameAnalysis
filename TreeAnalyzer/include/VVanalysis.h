// Dear emacs, this is -*- c++ -*-

#ifndef VVanalysis_H
#define VVanalysis_H

// STL include(s):
#include <vector>
#include <string>

// SFrame include(s):
#include "core/include/SCycleBase.h"
#include "plug-ins/include/SSummedVar.h"


// ROOT include(s):
#include <TBits.h>
#include <TFile.h>

// External include(s):
#include "../NtupleVariables/include/JetNtupleObject.h"
#include "../NtupleVariables/include/Jet.h"
#include "../NtupleVariables/include/EventInfoNtupleObject.h"
#include "../NtupleVariables/include/ElectronNtupleObject.h"
#include "../NtupleVariables/include/Electron.h"
#include "../NtupleVariables/include/MuonNtupleObject.h"
#include "../NtupleVariables/include/Muon.h"
#include "../NtupleVariables/include/MissingEtNtupleObject.h"
#include "../NtupleVariables/include/MissingEt.h"
#include "../NtupleVariables/include/GenParticleNtupleObject.h"
#include "../NtupleVariables/include/GenParticle.h"
// #include "../GoodRunsLists/include/TGoodRunsList.h"
#include "../PileupReweightingTool/include/PileupReweightingTool.h"
#include "../BTaggingTools/include/BTaggingScaleTool.h"
#include "../TreeAnalyzer/include/LumiWeight.h"

class TH1D;
class TH2D;
class TRandom3;
class TBits;
namespace UZH {
  class Jet;
  class Electron;
  class Muon;
  class MissingEt;
  class GenParticle;
}

/**
 *   @short Analyse flat UZH ntuple
 *
 *
 *  @author Thea Arrestad
 * @version $Revision: 1 $
 */
class VVanalysis : public SCycleBase {

public:
   /// Default constructor
   VVanalysis();
   /// Default destructor
   ~VVanalysis();

   /// Function called at the beginning/end of the cycle
   virtual void BeginCycle() throw( SError );
   virtual void EndCycle() throw( SError );
   
   /// Function called at the beginning of/finished processing a new input data
   virtual void BeginInputData( const SInputData& ) throw( SError );
   virtual void EndInputData  ( const SInputData& ) throw( SError );
   
   ///  Function called before/after processing each InputData block
   virtual void BeginMasterInputData( const SInputData& ) throw( SError );
   virtual void EndMasterInputData  ( const SInputData& ) throw( SError );
   


   /// Function called after opening each new input file
   virtual void BeginInputFile( const SInputData& ) throw( SError );

   /// Function called for every event
   virtual void ExecuteEvent( const SInputData&, Double_t ) throw( SError );
   
/// Function to obtain event weights for MC
   virtual double getEventWeight();
   
private:
   
  // Input variable objects:
  Ntuple::JetNtupleObject         m_jetAK4;       ///< jet container
  Ntuple::JetNtupleObject         m_jetAK8;       ///< jet container
  Ntuple::JetNtupleObject         m_jetAK8Puppi;  ///< jet container
  Ntuple::EventInfoNtupleObject   m_eventInfo;    ///< event info container
  Ntuple::ElectronNtupleObject    m_electron;     ///< electron container
  Ntuple::MuonNtupleObject        m_muon;         ///< muon container
  Ntuple::MissingEtNtupleObject   m_missingEt;    ///< missing E_T container
  Ntuple::GenParticleNtupleObject m_genParticle;  ///< gen particle container
  
  
  // Further objects
  // Root::TGoodRunsList m_grl;
  PileupReweightingTool m_pileupReweightingTool;
  BTaggingScaleTool     m_bTaggingScaleTool;
  LumiWeight            m_xSec;
  
  // Some counters:
  //
  SSummedVar< Int_t > m_allEvents; //!
  SSummedVar< Int_t > m_foundLepton; //!
  SSummedVar< Int_t > m_passedPuppi; //!
  SSummedVar< Int_t > m_passedAK4; //!
  SSummedVar< Int_t > m_passedMET; //!
  SSummedVar< Int_t > m_passedEvents; //!
  SSummedVar< std::vector< Int_t > > m_test; //!
  
  Ntuple::JetNtupleObject         m_genjetAK8;    ///< jet container
  
  float     nSumGenWeights;
  
  // The output variables
  float     m_o_mpuppisoftdrop  ; 
  float     m_o_mgensoftdrop    ; 
  float     m_o_tau1            ; 
  float     m_o_tau2            ; 
  float     m_o_tau3            ; 
  float     m_o_tau21           ;
  float     m_o_csv             ;  
  float     m_o_genpt           ;
  float     m_o_pt              ;
  float     m_o_eta             ;
  int       m_o_nJ             ;
  int       m_o_nLeptons        ;
  float     m_o_Wlep_pt              ;
  float     m_o_lep_pt              ;
  float     m_o_lep_eta             ;
  float     m_o_lep_phi             ;
  
  
  int       Flag_goodVertices             ; 
  int       Flag_globalTightHalo2016Filter; 
  int       Flag_HBHENoiseFilter          ; 
  int       Flag_HBHENoiseIsoFilter       ; 
  int       Flag_eeBadScFilter            ; 
  int       Flag_badChargedHadronFilter   ; 
  int       Flag_badMuonFilter            ;    
  int       Flag_ECALDeadCell             ;                                     
  float     HLTMu50          ; 
  float     HLTMu45_eta2p1           ; 
  float     HLTEle105_CaloIdVT_GsfTrkIdT      ; 
  float     HLTEle115_CaloIdVT_GsfTrkIdT      ; 
  bool      HLT_all                        ;
  int       nLeptonOverlap                ; 
  int       mergedVTruth          ; 
  double    b_xSec                        ;
  int       b_event                       ;
  int       b_lumi                        ;
  int       b_run                         ;
  float     b_weight                      ;
  float     b_weightGen                   ;
  float     b_weightPU                    ;
  float     b_weightBtag                  ;                 
  
  
 
  
  //naming
  std::string m_jetAK8Name;
  std::string m_jetAK4Name;
  std::string m_genjetAK8Name;
  std::string m_muonName;
  std::string m_electronName;
  std::string m_missingEtName;
  std::string m_jetAK8PuppiName;
  std::string m_genParticleName;
  std::string m_PUPPIJEC;

  std::vector<TF1*> m_puppisd_corr      ;

  std::vector<int> debugEvent;
  
  // Names of the input/output trees:
  std::string m_recoTreeName;       ///< name of tree with reconstructed objects in ntuple
  std::string m_outputTreeName;    ///< name of output tree
  
  enum ValHistsType { GENERAL, ELECTRON, MUON, JETS };
  void FillValidationHists( ValHistsType, const TString& status );
  void clearBranches( );
  
  std::string m_jsonName;
  std::string m_Channel;
  bool        m_isData;
  bool        m_isSignal;
  
  
  

   // Macro adding the functions for dictionary generation
   ClassDef( VVanalysis, 0 );

}; // class VVanalysis

#endif // VVanalysis_H

