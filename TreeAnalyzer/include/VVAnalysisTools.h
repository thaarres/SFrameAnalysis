#ifndef VVANALYSISTOOLS_H
#define VVANALYSISTOOLS_H

#include "../include/VVanalysis.h"
#include <TF1.h>

bool isMergedVJet(TLorentzVector goodFatJet, std::vector<UZH::GenParticle> VdecayProducts);
bool isMergedTopJet(TLorentzVector goodFatJet, std::vector<UZH::GenParticle> TopdecayProducts);

std::vector<UZH::GenParticle> FindGeneratedQuarks(Ntuple::GenParticleNtupleObject m_genParticle, bool m_isData, bool findW);

float ApplyPuppiSoftdropMassCorrections(UZH::Jet puppiJet,std::vector<TF1*> m_puppisd_corr, bool m_isData);

bool FoundNoLeptonOverlap(std::vector<UZH::Electron> goodEle, std::vector<UZH::Muon> goodMu, TLorentzVector Jet, float dRmax);

std::vector<UZH::Electron> FindGoodLeptons(Ntuple::ElectronNtupleObject m_electrons);
std::vector<UZH::Muon> FindGoodLeptons(Ntuple::MuonNtupleObject m_muons);
std::vector<UZH::Jet>  FindGoodJetsAK4(Ntuple::JetNtupleObject m_jetAK4,std::vector<UZH::Electron> goodElectrons,std::vector<UZH::Muon> goodMuons,TLorentzVector Jet);

bool SignalIsHad( Ntuple::GenParticleNtupleObject data_ , std::string m_Channel);

std::vector<UZH::Jet> SortAfterPuppiSDMass(std::vector <UZH::Jet> jets);
std::vector<UZH::Jet> SortAfterTau21(std::vector <UZH::Jet> jets);

void PrintEvent(std::vector<UZH::Jet> );

TLorentzVector NuMomentum( float leptonPx, float leptonPy, float leptonPz, float leptonPt, float leptonE, float metPx, float metPy );
  
#endif
