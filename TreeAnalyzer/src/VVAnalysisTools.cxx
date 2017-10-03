#include "../include/VVAnalysisTools.h"
#include "../include/VVanalysis.h"
#include "include/EquationSolver.h"

bool isMergedVJet(TLorentzVector goodFatJet, std::vector<UZH::GenParticle> VdecayProducts) {
    int associatedQuarks=0;
    for(unsigned int i=0;i< VdecayProducts.size();i++)
    {
      TLorentzVector q = VdecayProducts[i].tlv();
      if (goodFatJet.DeltaR(q) < 0.8) associatedQuarks +=1;

    }
    if (associatedQuarks ==2) return 1;
    else return 0;   
}


std::vector<UZH::GenParticle> FindGeneratedQuarks(Ntuple::GenParticleNtupleObject m_genParticle,bool m_isData)
{
  std::vector<UZH::GenParticle> GenQuarks;
  if(m_isData) return GenQuarks;
  for(int i=0;i< m_genParticle.N;i++)
  {
    bool containsV=0;
    if(TMath::Abs(m_genParticle.pdgId->at(i)) > 6 or TMath::Abs(m_genParticle.pdgId->at(i)) <1) continue;
       
    for(int ii=0;ii< m_genParticle.nMoth->at(i);ii++)
    {
           
      if(TMath::Abs(m_genParticle.mother->at(i).at(ii)) == 24 or TMath::Abs(m_genParticle.mother->at(i).at(ii)) == 23) containsV = 1;  

    }
    if (!containsV) continue;
    UZH::GenParticle MyQuark( &m_genParticle, i );
    GenQuarks.push_back(MyQuark);
  }
  return GenQuarks;    
}


float ApplyPuppiSoftdropMassCorrections(UZH::Jet puppiJet,std::vector<TF1*> m_puppisd_corr, bool m_isData){
 float genCorr =1;
 if(!m_isData) genCorr = m_puppisd_corr[0]->Eval(puppiJet.pt());
 float recoCorr = 1;
 if( fabs(puppiJet.eta()) <= 1.3) recoCorr = m_puppisd_corr[1]->Eval(puppiJet.pt());
 else if (fabs(puppiJet.eta()) > 1.3) recoCorr = m_puppisd_corr[2]->Eval(puppiJet.pt());
    
 return puppiJet.softdrop_mass()*genCorr*recoCorr;
}


bool FoundNoLeptonOverlap(std::vector<UZH::Electron> goodEle, std::vector<UZH::Muon> goodMu, TLorentzVector Jet, float dRmax){
  for( unsigned int e =0;e< goodEle.size(); e++)
  { 
    if(Jet.DeltaR(goodEle.at(e).tlv()) < dRmax) return 0;  
  }
  for( unsigned int m =0;m< goodMu.size(); m++)
  { 
    if(Jet.DeltaR(goodMu.at(m).tlv()) < dRmax) return 0;   
  }
  return 1;   
}

std::vector<UZH::Electron> FindGoodLeptons(Ntuple::ElectronNtupleObject m_electrons){
    std::vector<UZH::Electron> goodEle;
    for(int i=0;i< m_electrons.N;i++)
    {
        UZH::Electron ele(&m_electrons,i);
        if( ! ele.isHeepElectron()  ) continue;
        if( ele.pt() <= 120.) continue;
        goodEle.push_back(ele);
    }
    return goodEle;
}

std::vector<UZH::Muon> FindGoodLeptons(Ntuple::MuonNtupleObject m_muons){
    std::vector<UZH::Muon> goodMu;
    for(int i=0;i< m_muons.N;i++)
    {
        UZH::Muon mu(&m_muons,i);
        if( mu.pt() <= 53.) continue;
        if( fabs( mu.eta() ) >= 2.1 ) continue;
        if( ! mu.isHighPtMuon()  ) continue; // //isHighPtMuon
        if( mu.trackIso()/mu.pt() >= 0.1 ) continue;
        goodMu.push_back(mu);
    }
    return goodMu;
}

std::vector<UZH::Jet> FindGoodJetsAK4(Ntuple::JetNtupleObject m_jetAK4, std::vector<UZH::Electron> goodElectrons, std::vector<UZH::Muon> goodMuons, TLorentzVector Jet){
    std::vector<UZH::Jet> goodJet;
    for(int i=0;i< m_jetAK4.N;i++)
    {
        UZH::Jet jet(&m_jetAK4,i);
        if( jet.csv() <= 0.89) continue;
        if( jet.pt() <= 30.) continue;
        if( fabs( jet.eta() ) >= 2.4 ) continue;
        if( ! jet.IDLoose()  ) continue; 
        if( !(FoundNoLeptonOverlap(goodElectrons,goodMuons,Jet, 0.3 ) ) ) continue;
        if(Jet.DeltaR(jet.tlv()) < 0.8) continue;  
        goodJet.push_back(jet);
    }
    return goodJet;
}




bool SignalIsHad( Ntuple::GenParticleNtupleObject data_ , std::string m_Channel) {
  
  bool isSignal = false;
  int nVs = 0;
  
  for( int p = 0; p < (data_.N); ++p ){
    bool isHad = false;
    if( (*data_.pt).at(p) < 0.01) continue;
    if( fabs((*data_.pdgId).at(p)) != 24 and fabs((*data_.pdgId).at(p)) != 23  and fabs((*data_.pdgId).at(p)) != 25) continue;
    if ((*data_.dau)[p].size() < 2) continue;
    for( unsigned int d = 0; d < (*data_.dau)[p].size(); ++d ) {
      if( fabs((*data_.dau)[p][d]) > 9) continue;
      isHad = true;
    }
    if (!isHad) continue;
    nVs ++;
  }
  
  if (m_Channel.find("qV")!=std::string::npos && nVs>0) isSignal = true;
  if (m_Channel.find("VV")!=std::string::npos && nVs>1) isSignal = true;
  return isSignal;
}


std::vector<UZH::Jet> SortAfterPuppiSDMass(std::vector <UZH::Jet> jets){
    std::vector<float> mass;
    std::map<float,int> map;
    for(unsigned int i=0;i< jets.size(); i++)
    {
      mass.push_back(jets[i].puppi_softdropmass);
      map[jets[i].puppi_softdropmass] = i;
    }
    std::sort(mass.begin(),mass.end());
    std::vector<UZH::Jet> sorted;
    for(unsigned int i=0;i<mass.size();i++)
    {
      sorted.push_back(jets.at(map[mass.at(mass.size()-i-1)]));
    }
    return sorted;
}

std::vector<UZH::Jet> SortAfterTau21(std::vector <UZH::Jet> jets){
    std::vector<float> tau21;
    std::map<float,int> map;
    for(unsigned int i=0;i< jets.size(); i++)
    {
      tau21.push_back(jets[i].puppi_tau2/jets[i].puppi_tau1);
      map[jets[i].puppi_tau2/jets[i].puppi_tau1] = i;
    }
    std::sort(tau21.begin(),tau21.end(),std::greater<float>() );
    std::vector<UZH::Jet> sorted;
    for(unsigned int i=0;i<tau21.size();i++)
    {
      sorted.push_back(jets.at(map[tau21.at(tau21.size()-i-1)]));
    }
    return sorted;
}



void PrintEvent(std::vector<UZH::Jet> jets  )
{
    std::cout << "D_eta " << fabs(jets.at(0).tlv().Eta()-jets.at(1).tlv().Eta()) << std::endl;
    std::cout << "Mjj   " << (jets.at(0).tlv()+jets.at(1).tlv()).M() <<std::endl;
    for(unsigned int i=0;i<jets.size(); i++)
    {
        std::cout << "jet pt " << jets.at(i).tlv().Pt() << " eta " << jets.at(i).tlv().Eta() << " phi " << jets.at(i).tlv().Phi() << " e " << jets.at(i).tlv().E() << " puppi sd mass "<< jets.at(i).puppi_softdropmass << std::endl;
    
    }
}

TLorentzVector NuMomentum( float leptonPx, float leptonPy, float leptonPz, float leptonPt, float leptonE, float metPx, float metPy ){

  double  mW = 80.399;

  TLorentzVector result;

  //  double Wmt = sqrt(pow(Lepton.et()+MET.pt(),2) - pow(Lepton.px()+metPx,2) - pow(leptonPy+metPy,2) );

  double MisET2 = (metPx * metPx + metPy * metPy);
  double mu = (mW * mW) / 2 + metPx * leptonPx + metPy * leptonPy;
  double a  = (mu * leptonPz) / (leptonE * leptonE - leptonPz * leptonPz);
  double a2 = TMath::Power(a, 2);
  double b  = (TMath::Power(leptonE, 2.) * (MisET2) - TMath::Power(mu, 2.)) / (TMath::Power(leptonE, 2) - TMath::Power(leptonPz, 2));
  double pz1(0), pz2(0), pznu(0);
  int nNuSol(0);

  TLorentzVector p4nu_rec;
  TLorentzVector p4W_rec;
  TLorentzVector p4b_rec;
  TLorentzVector p4Top_rec;
  TLorentzVector p4lep_rec;

  p4lep_rec.SetPxPyPzE(leptonPx, leptonPy, leptonPz, leptonE);

  TLorentzVector p40_rec(0, 0, 0, 0);

  if (a2 - b > 0 )
    {
      //if(!usePositiveDeltaSolutions_)
      //  {
      //  result.push_back(p40_rec);
      //  return result;
      //  }
      double root = sqrt(a2 - b);
      pz1 = a + root;
      pz2 = a - root;
      nNuSol = 2;

      //    if(usePzPlusSolutions_)pznu = pz1;
      //    if(usePzMinusSolutions_)pznu = pz2;
      //if(usePzAbsValMinimumSolutions_){
      pznu = pz1;
      if (fabs(pz1) > fabs(pz2)) pznu = pz2;
      //}


      double Enu = sqrt(MisET2 + pznu * pznu);

      p4nu_rec.SetPxPyPzE(metPx, metPy, pznu, Enu);

      //    result =.push_back(p4nu_rec);
      result = p4nu_rec;

    }
  else
    {

      // if(!useNegativeDeltaSolutions_){
      //result.push_back(p40_rec);
      //  return result;
      //    }
      //    double xprime = sqrt(mW;


      double ptlep = leptonPt, pxlep = leptonPx, pylep = leptonPy, metpx = metPx, metpy = metPy;

      double EquationA = 1;
      double EquationB = -3 * pylep * mW / (ptlep);
      double EquationC = mW * mW * (2 * pylep * pylep) / (ptlep * ptlep) + mW * mW - 4 * pxlep * pxlep * pxlep * metpx / (ptlep * ptlep) - 4 * pxlep * pxlep * pylep * metpy / (ptlep * ptlep);
      double EquationD = 4 * pxlep * pxlep * mW * metpy / (ptlep) - pylep * mW * mW * mW / ptlep;

      std::vector<long double> solutions = EquationSolve<long double>((long double)EquationA, (long double)EquationB, (long double)EquationC, (long double)EquationD);

      std::vector<long double> solutions2 = EquationSolve<long double>((long double)EquationA, -(long double)EquationB, (long double)EquationC, -(long double)EquationD);


      double deltaMin = 14000 * 14000;
      double zeroValue = -mW * mW / (4 * pxlep);
      double minPx = 0;
      double minPy = 0;

      //    std::cout<<"a "<<EquationA << " b " << EquationB  <<" c "<< EquationC <<" d "<< EquationD << std::endl;

      //  if(usePxMinusSolutions_){
      for ( int i = 0; i < (int)solutions.size(); ++i)
        {
	  if (solutions[i] < 0 ) continue;
	  double p_x = (solutions[i] * solutions[i] - mW * mW) / (4 * pxlep);
	  double p_y = ( mW * mW * pylep + 2 * pxlep * pylep * p_x - mW * ptlep * solutions[i]) / (2 * pxlep * pxlep);
	  double Delta2 = (p_x - metpx) * (p_x - metpx) + (p_y - metpy) * (p_y - metpy);

	  //      std:://cout<<"intermediate solution1 met x "<<metpx << " min px " << p_x  <<" met y "<<metpy <<" min py "<< p_y << std::endl;

	  if (Delta2 < deltaMin && Delta2 > 0)
            {
	      deltaMin = Delta2;
	      minPx = p_x;
	      minPy = p_y;
            }
	  //     std:://cout<<"solution1 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl;
        }

      //    }

      //if(usePxPlusSolutions_){
      for ( int i = 0; i < (int)solutions2.size(); ++i)
        {
	  if (solutions2[i] < 0 ) continue;
	  double p_x = (solutions2[i] * solutions2[i] - mW * mW) / (4 * pxlep);
	  double p_y = ( mW * mW * pylep + 2 * pxlep * pylep * p_x + mW * ptlep * solutions2[i]) / (2 * pxlep * pxlep);
	  double Delta2 = (p_x - metpx) * (p_x - metpx) + (p_y - metpy) * (p_y - metpy);
	  //  std:://cout<<"intermediate solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl;
	  if (Delta2 < deltaMin && Delta2 > 0)
            {
	      deltaMin = Delta2;
	      minPx = p_x;
	      minPy = p_y;
            }
	  //  std:://cout<<"solution2 met x "<<metpx << " min px " << minPx  <<" mult;
	}
      //}

      double pyZeroValue = ( mW * mW * pxlep + 2 * pxlep * pylep * zeroValue);
      double delta2ZeroValue = (zeroValue - metpx) * (zeroValue - metpx) + (pyZeroValue - metpy) * (pyZeroValue - metpy);

      if (deltaMin == 14000 * 14000)return result;
      //    else std:://cout << " test " << std::endl;

      if (delta2ZeroValue < deltaMin)
        {
	  deltaMin = delta2ZeroValue;
	  minPx = zeroValue;
	  minPy = pyZeroValue;
        }

      //    std:://cout<<" MtW2 from min py and min px "<< sqrt((minPy*minPy+minPx*minPx))*ptlep*2 -2*(pxlep*minPx + pylep*minPy)  <<std::endl;
      ///    ////Y part

      double mu_Minimum = (mW * mW) / 2 + minPx * pxlep + minPy * pylep;
      double a_Minimum  = (mu_Minimum * leptonPz) / (leptonE * leptonE - leptonPz * leptonPz);
      pznu = a_Minimum;

      //if(!useMetForNegativeSolutions_){
      double Enu = sqrt(minPx * minPx + minPy * minPy + pznu * pznu);
      p4nu_rec.SetPxPyPzE(minPx, minPy, pznu , Enu);

      //    }
      //    else{
      //      pznu = a;
      //      double Enu = sqrt(metpx*metpx+metpy*metpy + pznu*pznu);
      //      p4nu_rec.SetPxPyPzE(metpx, metpy, pznu , Enu);
      //    }

      //      result.push_back(p4nu_rec);
      result = p4nu_rec;
    }
  return result;
}












