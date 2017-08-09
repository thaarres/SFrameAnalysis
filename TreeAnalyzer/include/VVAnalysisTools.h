#ifndef VVANALYSISTOOLS_H
#define VVANALYSISTOOLS_H

#include "../include/VVanalysis.h"

bool isMergedVJet(TLorentzVector goodFatJet, std::vector<UZH::GenParticle> VdecayProducts);

void FindGeneratedQuarks(Ntuple::GenParticleNtupleObject m_genParticle, bool m_isData , std::vector<UZH::GenParticle> GenQuarks);



#endif
