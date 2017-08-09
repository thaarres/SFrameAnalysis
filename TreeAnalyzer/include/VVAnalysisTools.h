#ifndef VVANALYSISTOOLS_H
#define VVANALYSISTOOLS_H

#include "../include/VVanalysis.h"
#include <TF1.h>

bool isMergedVJet(TLorentzVector goodFatJet, std::vector<UZH::GenParticle> VdecayProducts);

void FindGeneratedQuarks(Ntuple::GenParticleNtupleObject m_genParticle, bool m_isData , std::vector<UZH::GenParticle> GenQuarks);

float ApplyPuppiSoftdropMassCorrections(UZH::Jet puppiJet,std::vector<TF1*> m_puppisd_corr, bool m_isData);

#endif
