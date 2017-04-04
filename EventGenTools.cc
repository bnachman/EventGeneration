#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  

#include "EventGenTools.h"
#include "myFastJetBase.h"

#include "TRandom3.h"
#include "TError.h"
#include "TVector3.h"

using namespace std;

// Constructor 
EventGenTools::EventGenTools(){
    m_test = 0;
}

fastjet::PseudoJet EventGenTools::Add(fastjet::PseudoJet jet){
  fastjet::PseudoJet out(0.,0.,0.,0.);
  TRandom3 *rand = new TRandom3(0);
  double myang = rand->Uniform(0,2*3.141);
  double myrad = rand->Uniform(0,1);
  myrad = sqrt(myrad);
  double eta = jet.eta()+0.4*myrad*cos(myang);
  double phi = jet.phi()+0.4*myrad*sin(myang);
  fastjet::PseudoJet out2(0.,0.,0.,0.);
  out2.reset_momentum_PtYPhiM(0.5,eta,phi,0.);
  out = jet;
  out+=out2;
  return out;
}

fastjet::PseudoJet EventGenTools::Angles(fastjet::PseudoJet jet){
  fastjet::PseudoJet out(0.,0.,0.,0.);
  TRandom3 *rand = new TRandom3(0);
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    double phi_smear = rand->Gaus(jet.constituents()[i].phi(),0.005);
    double eta_smear = rand->Gaus(jet.constituents()[i].eta(),0.005);
    fastjet::PseudoJet out2;
    //out2.reset_momentum_PtYPhiM(jet.constituents()[i].e()/cosh(eta_smear),eta_smear,phi_smear,0.);
    out2.reset_momentum_PtYPhiM(jet.constituents()[i].pt(),eta_smear,phi_smear,0.);
    out+=out2;
  }
  return out;
}

fastjet::PseudoJet EventGenTools::Drop(fastjet::PseudoJet jet){
  fastjet::PseudoJet out(0.,0.,0.,0.);
  TRandom3 *rand = new TRandom3(0);
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    double r = 0.25*exp(-2*jet.constituents()[i].e());
    double flip = rand->Uniform(0,1);
    if (( r < flip) && (jet.constituents()[i].e() <2.5))
      continue;
    out+=jet.constituents()[i];
  }
  return out;
}

fastjet::PseudoJet EventGenTools::Scale(fastjet::PseudoJet jet){
  fastjet::PseudoJet out(0.,0.,0.,0.);
  TRandom3 *rand = new TRandom3(0);
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    double fact = 0.05*(1+1.5/jet.pt());
    double fact_smear = rand->Gaus(1,fact);
    out+=jet.constituents()[i]*fact_smear;
  }
  return out;
}

double EventGenTools::width(fastjet::PseudoJet jet){
  double out=0.;
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    out+= (jet.constituents()[i].pt()/jet.pt())*pow(jet.constituents()[i].delta_R(jet),2);
  }
  return sqrt(out);
}

int EventGenTools::ntrack(fastjet::PseudoJet jet,double pT){
  int out=0.;
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    if (fabs(jet.constituents()[i].user_info<MyUserInfo>().charge()) > 0.1 && jet.constituents()[i].pt() > pT){
      //std::cout << jet.constituents()[i].user_info<MyUserInfo>().charge() << std::endl;
      out++;
    }
  }
  return out;
}

double EventGenTools::JetCharge(fastjet::PseudoJet jet,double kappa){
  //Returns the jet charge with weighting factor kappa
  double charge=0.;
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    charge+=jet.constituents()[i].user_info<MyUserInfo>().charge()*pow(jet.constituents()[i].pt(),kappa);
  }
  return charge/pow(jet.pt(),kappa);
}

bool EventGenTools::IsBHadron(int pdgId){
  int abs_pdgId = abs(pdgId);
  int abs_pdgId_mod10k = (abs(pdgId)%10000);
  if( (abs_pdgId_mod10k>=500 && abs_pdgId_mod10k<600) /*mesons*/  ||
      (abs_pdgId>=5000      && abs_pdgId<6000)      /*baryons*/   )
    return true;

  return false;
}

bool EventGenTools::IsCHadron(int pdgId){
  int abs_pdgId = abs(pdgId);
  int abs_pdgId_mod10k = (abs(pdgId)%10000);
  if( (abs_pdgId_mod10k>=400 && abs_pdgId_mod10k<500) /*mesons*/  ||
      (abs_pdgId>=4000      && abs_pdgId<5000)      /*baryons*/   )
    return true;

  return false;
}

bool EventGenTools::Btag(fastjet::PseudoJet jet,vector<fastjet::PseudoJet> bhadrons,vector<fastjet::PseudoJet> chadrons,double jetrad,double b,double c,double uds){

  TRandom3 *rand = new TRandom3(0);

  int foundb=0;
  int foundc=0;
  
  for (unsigned int i=0; i<bhadrons.size(); i++){
    if (bhadrons[i].delta_R(jet)<jetrad){
      foundb=1;
    }
  }

  for (unsigned int i=0; i<chadrons.size(); i++){
    if (chadrons[i].delta_R(jet)<jetrad){
      if (foundb!=1) foundc=1;
    }
  }

  if (foundb==1){
    double flip = rand->Uniform(0.,1.);
    if (flip < b){
      delete rand;
      return true;
    }
  }
  if (foundc==1){
    double flip= rand->Uniform(0.,1.);
    if (flip < 1./c){
      delete rand;
      return true;
    }
  }
  double flip= rand->Uniform(0.,1.);
  if (flip < 1./uds){
    delete rand;
    return true;
  }
  
  delete rand;
  return false;
}

bool EventGenTools::BosonMatch(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> Bosons, double jetrad, int BosonID ){

  for (unsigned int i=0; i<Bosons.size(); i++){
      if (Bosons[i].user_info<MyUserInfo>().pdg_id() != BosonID) continue;
      if (Bosons[i].delta_R(jet)<jetrad){
        return true;
      }
  }
  return false;
}

bool EventGenTools::IsIsolated(Pythia8::Particle* particle, Pythia8::Pythia* pythia8, float rel_iso, float ConeSize){
    float sumpT=0;
    fastjet::PseudoJet part(particle->px(), particle->py(), particle->pz(),particle->e() );
    for (unsigned int ip=0; ip<pythia8->event.size(); ++ip){
        if (!pythia8->event[ip].isFinal() )      continue;
        if (fabs(pythia8->event[ip].id())  ==12) continue;
        if (fabs(pythia8->event[ip].id())  ==14) continue;
        if (fabs(pythia8->event[ip].id())  ==16) continue;
        if (pythia8->event[ip].pT()       < 0.5) continue;
        if (&pythia8->event[ip] == particle    ) continue; //same particle
        fastjet::PseudoJet p(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e() );
        if(p.delta_R(part)>ConeSize)             continue;
        sumpT+=p.pt();
    }
    if(sumpT/part.pt()>rel_iso) return false;
    else return true;
}
