#ifndef VVANALYSISTOOLS_H
#define VVANALYSISTOOLS_H

#include "../include/VVanalysis.h"

bool isMergedVJet(TLorentzVector goodFatJet, std::vector<UZH::GenParticle> VdecayProducts);

std::vector<UZH::GenParticle> FindGeneratedQuarks(Ntuple::GenParticleNtupleObject m_genParticle, bool m_isData);



#endif
