#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "Math/ProbFunc.h"

#include "TVector3.h"

#include "EventGenAnalysis.h"
#include "EventGenTools.h"

#include "myFastJetBase.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"

#include "Pythia8/Pythia.h"

//Alternative methods
#include "fastjet/contrib/SoftKiller.hh"
#include "fastjet/contrib/VariableRPlugin.hh"

using namespace std;
using namespace fastjet;

// Constructor 
EventGenAnalysis::EventGenAnalysis(){
    if(fDebug) cout << "EventGenAnalysis::EventGenAnalysis Start " << endl;
    ftest = 0;
    fDebug = false;
    fOutName = "test.txt";
    tool = new EventGenTools();
    
    //myfile.open ("example.txt");

    // jet def 
    m_jet_def               = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);//1.0);
    m_jet_def_forrho=  new fastjet::JetDefinition(fastjet::kt_algorithm,0.4);
    if(fDebug) cout << "EventGenAnalysis::EventGenAnalysis End " << endl;
}

// Destructor 
EventGenAnalysis::~EventGenAnalysis(){
    delete tool;
    delete m_jet_def;
}

// Begin method
void EventGenAnalysis::Begin(){
   myfile.open (fOutName.c_str());
   return;
}

// End
void EventGenAnalysis::End(){
    myfile.close();
    return;
}

// Analyze
void EventGenAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia_MB, int NPV){
    if(fDebug) cout << "EventGenAnalysis::AnalyzeEvent Begin " << endl;

    // -------------------------
    if (!pythia8->next()) return;
    if(fDebug) cout << "EventGenAnalysis::AnalyzeEvent Event Number " << ievt << endl;

    // new event-----------------------
    fTEventNumber = ievt;
    std::vector <fastjet::PseudoJet>           particlesForJets;
    std::vector <fastjet::PseudoJet>           HS;
    std::vector <double> alpha_vals;

    for (int iPU = 0; iPU <= NPV; ++iPU) {
      for (int i = 0; i < pythia_MB->event.size(); ++i) {
        if (!pythia_MB->event[i].isFinal()    ) continue;
        if (fabs(pythia_MB->event[i].id())==12) continue;
        if (fabs(pythia_MB->event[i].id())==13) continue;
        if (fabs(pythia_MB->event[i].id())==14) continue;
        if (fabs(pythia_MB->event[i].id())==16) continue;

	fastjet::PseudoJet p(pythia_MB->event[i].px(), pythia_MB->event[i].py(), pythia_MB->event[i].pz(),pythia_MB->event[i].e() );
	p.reset_PtYPhiM(p.pt(), p.rapidity(), p.phi(), 0.);
	p.set_user_info(new MyUserInfo(pythia_MB->event[i].id(),i,iPU,pythia_MB->event[i].charge(),0.,0.));
	particlesForJets.push_back(p);
      }
      if (!pythia_MB->next()) continue;
    }

    // Particle loop -----------------------------------------------------------
    double quark_pT = -1;
    //double quark_pT =  1000000000;
    int ndaughts = 0;

    while (quark_pT < 95){
      ndaughts = 0;
      quark_pT = 1000000000.;
      for (unsigned int ip=0; ip<pythia8->event.size(); ++ip){
	if (abs(pythia8->event[ip].id()) == 25){
	  vector<int> daugthers = pythia8->event.daughterList(ip);
	  for (unsigned int j=0; j<daugthers.size(); j++){
	    if (pythia8->event[daugthers[j]].id()==25) continue;
	    ndaughts++;
	    if (pythia8->event[daugthers[j]].pT() < quark_pT) quark_pT = pythia8->event[daugthers[j]].pT();
	  }
	}
	if (ndaughts > 1) break;
      }
      if (quark_pT < 95) pythia8->next();
    }
    //std::cout << ndaughts << " " << quark_pT << std::endl;

    for (unsigned int ip=0; ip<pythia8->event.size(); ++ip){
      
      fastjet::PseudoJet p(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e());
      p.reset_PtYPhiM(p.pt(), p.rapidity(), p.phi(), 0.);
      p.set_user_info(new MyUserInfo(pythia8->event[ip].id(),ip,-1,pythia8->event[ip].charge(),0.,0.));

      // particles for jets --------------
      if (!pythia8->event[ip].isFinal() )      continue;
      if (fabs(pythia8->event[ip].id())  ==12) continue;
      if (fabs(pythia8->event[ip].id())  ==13) continue;
      if (fabs(pythia8->event[ip].id())  ==14) continue;
      if (fabs(pythia8->event[ip].id())  ==16) continue;
      particlesForJets.push_back(p);
      HS.push_back(p);
      
    } // end particle loop -----------------------------------------------

    if (quark_pT < 95) return; //as in the cleansing paper.

    fastjet::AreaDefinition area_def(fastjet::active_area);
    fastjet::ClusterSequenceArea csa(particlesForJets,*m_jet_def_forrho,area_def);
    fastjet::RangeDefinition range(4.0);
    double median_pt = csa.median_pt_per_unit_area(range);

    //std::cout << HS.size() << std::endl;
    fastjet::ClusterSequence hs(HS, *m_jet_def);
    vector<fastjet::PseudoJet> myHSJets = fastjet::sorted_by_pt(hs.inclusive_jets(25.0));
    //std::cout << myHSJets.size() << std::endl;
    if (myHSJets.size() < 2) return;

    //PUPPI; https://arxiv.org/pdf/1407.6013.pdf 
    double R0 = 0.3;
    double Rmin = 0.02;
    double wcut = 0.1;
    double pTcut = 0.1 + NPV * 0.007;
    
    for (unsigned int i=0; i<particlesForJets.size(); i++){
      //if (i%100==0) std::cout << i << " " << particlesForJets.size() << std::endl;
      bool ischargedi = (fabs(particlesForJets[i].user_info<MyUserInfo>().charge()) > 0.1);
      int whichvertexi =particlesForJets[i].user_info<MyUserInfo>().ipileup();
      double alphai = 0.;
      for (unsigned int j=0; j<particlesForJets.size(); j++){
	bool ischargedj = (fabs(particlesForJets[j].user_info<MyUserInfo>().charge()) > 0.1);
	int whichvertexj = particlesForJets[j].user_info<MyUserInfo>().ipileup();
	if (!ischargedj) continue;
	if (whichvertexj > -1) continue;
	if (particlesForJets[j].delta_R(particlesForJets[i]) > R0) continue;
	if (particlesForJets[j].delta_R(particlesForJets[i]) < Rmin) continue;
	double xi = particlesForJets[j].pt()/particlesForJets[j].delta_R(particlesForJets[i]);
	alphai+=xi;
      }
      if (alphai > 0){
	particlesForJets[i].set_user_info(new MyUserInfo(particlesForJets[i].user_info<MyUserInfo>().pdg_id(),particlesForJets[i].user_info<MyUserInfo>().pythia_id(),particlesForJets[i].user_info<MyUserInfo>().ipileup() ,particlesForJets[i].user_info<MyUserInfo>().charge(),log(alphai),0.));
	if (ischargedi && whichvertexi > -1) alpha_vals.push_back(log(alphai));
      }
    }
    std::sort(alpha_vals.begin(), alpha_vals.end());
    double median_alpha = 0.;
    if (alpha_vals.size() > 1) median_alpha = alpha_vals[int(alpha_vals.size()*0.5)];
    double left_RMS = 0.;
    double left_RMS_counter = 0.;
    for (unsigned int i=0; i<alpha_vals.size(); i++){
      if (alpha_vals[i] < median_alpha){
	left_RMS+=pow(alpha_vals[i]-median_alpha,2);
	left_RMS_counter+=1;
      }
    }
    if (left_RMS_counter > 0) left_RMS = left_RMS/left_RMS_counter;
    else left_RMS = 0.1;
    //std::cout << "median and left RMS " << median_alpha << " " << left_RMS << std::endl;
    for (unsigned int i=0; i<particlesForJets.size(); i++){
      double alpha = particlesForJets[i].user_info<MyUserInfo>().puppi();
      if (alpha != 0){
	if (alpha < median_alpha){
	  particlesForJets[i].set_user_info(new MyUserInfo(particlesForJets[i].user_info<MyUserInfo>().pdg_id(),particlesForJets[i].user_info<MyUserInfo>().pythia_id(),particlesForJets[i].user_info<MyUserInfo>().ipileup() ,particlesForJets[i].user_info<MyUserInfo>().charge(),0.,0.));
	}
	else{
	  double chi2 = pow(alpha-median_alpha,2)/left_RMS;
	  double w = ROOT::Math::chisquared_cdf(chi2,1);
	  if (w < wcut) w = 0.;
	  if (w*particlesForJets[i].pt() < pTcut) w = 0.; 
	  particlesForJets[i].set_user_info(new MyUserInfo(particlesForJets[i].user_info<MyUserInfo>().pdg_id(),particlesForJets[i].user_info<MyUserInfo>().pythia_id(),particlesForJets[i].user_info<MyUserInfo>().ipileup() ,particlesForJets[i].user_info<MyUserInfo>().charge(),w,0.));
	}
      } 
    }

    //std::cout << "done with puppi " << std::endl;

    //SOFTKILLER
    double grid_size = 0.4;
    double rapmax = 5.0;
    contrib::SoftKiller soft_killer(rapmax, grid_size);
    double pt_threshold;
    vector<PseudoJet> soft_killed_event;
    soft_killer.apply(particlesForJets, soft_killed_event, pt_threshold);

    for (unsigned int i=0; i<particlesForJets.size(); i++){
      if(std::find(soft_killed_event.begin(), soft_killed_event.end(), particlesForJets[i]) != soft_killed_event.end()) {
	//std::cout << "Found it !!! " << std::endl;
	particlesForJets[i].set_user_info(new MyUserInfo(particlesForJets[i].user_info<MyUserInfo>().pdg_id(),particlesForJets[i].user_info<MyUserInfo>().pythia_id(),particlesForJets[i].user_info<MyUserInfo>().ipileup() ,particlesForJets[i].user_info<MyUserInfo>().charge(),particlesForJets[i].user_info<MyUserInfo>().puppi(),1));
      } else {
	//std::cout << "NO CAN DO! " << std::endl;
      }
    }
    myfile << NPV << " " ;
    for (unsigned int i=0; i<particlesForJets.size(); i++){
      bool is_charged = (fabs(particlesForJets[i].user_info<MyUserInfo>().charge()) > 0.1);
      int vertex = particlesForJets[i].user_info<MyUserInfo>().ipileup();
      double pt = particlesForJets[i].pt();
      double eta = particlesForJets[i].eta();
      double phi = particlesForJets[i].phi();
      double m = particlesForJets[i].m();
      double puppi_weight = (is_charged && vertex < 0) ? 1. : particlesForJets[i].user_info<MyUserInfo>().puppi();
      double sk_weight = particlesForJets[i].user_info<MyUserInfo>().sk();
      if(std::find(soft_killed_event.begin(), soft_killed_event.end(), particlesForJets[i]) != soft_killed_event.end()) {
	sk_weight = 1;
      }
      myfile << pt << " " << eta << " " << phi << " " << is_charged << " " << vertex << " " << puppi_weight << " " << sk_weight << " ";
    }
    myfile << endl;
    //std::cout << "done with event " << std::endl;
    
    if(fDebug) cout << "EventGenAnalysis::AnalyzeEvent End " << endl;
    return;
}

// declate branches
void EventGenAnalysis::DeclareBranches(){
   return;
}


// resets vars
void EventGenAnalysis::ResetBranches(){

}
