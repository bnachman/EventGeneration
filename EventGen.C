#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "TString.h"
#include "TSystem.h"
#include "TError.h"
//#include "TPythia8.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"

#include "fastjet/PseudoJet.hh"  
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

#include "EventGenTools.h"
#include "EventGenAnalysis.h"

#include "boost/program_options.hpp"

using std::cout;
using std::endl;
using std::string;
using std::map;
using namespace std;
namespace po = boost::program_options;

int getSeed(int seed){
  if (seed > -1) return seed;
  int timeSeed = time(NULL);
  return abs(((timeSeed*181)*((getpid()-83)*359))%104729);
}

int main(int argc, char* argv[]){
    // argument parsing  ------------------------
    cout << "Called as: ";
    for(int ii=0; ii<argc; ++ii){
        cout << argv[ii] << " ";
    }
    cout << endl;

    // agruments 
    int nEvents = 0;
    int fDebug  = 0;
    string outName = "test.txt";
    int    seed      =-1;
    int npu = 0;

    po::options_description desc("Allowed options");
    desc.add_options() //100000
      ("help", "produce help message")
      ("NEvents", po::value<int>(&nEvents)->default_value(10000) ,    "Number of Events ")
      ("Debug",   po::value<int>(&fDebug) ->default_value(0) ,     "Debug flag")
      ("OutFile", po::value<string>(&outName)->default_value("test.txt"), "output file name")
      ("Seed",      po::value<int>(&seed)->default_value(-1), "seed. -1 means random seed")
      ("npu",      po::value<int>(&npu)->default_value(0), "npu")
      ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
        return 1;
    }
    //------

    seed = getSeed(seed);
    Pythia8::Pythia* pythia8b = new Pythia8::Pythia();
    pythia8b->readString("Random:setSeed = on");
    std::stringstream ss; ss << "Random:seed = " << seed;
    cout << ss.str() << endl;
    pythia8b->readString(ss.str());

    EventGenAnalysis * analysis1 = new EventGenAnalysis();

    //Signal
    pythia8b->readString("HiggsSM:ffbar2H = on");
    pythia8b->readString("25:m0=500");
    pythia8b->readString("25:doForceWidth=on");
    pythia8b->readString("25:mWidth=0.1");
    pythia8b->readString("25:onMode = off");
    pythia8b->readString("25:onIfAny = 1 2");

    //pythia8b->readString("PartonLevel:MPI = off");
    //pythia8b->readString("PartonLevel:ISR = off");
    //pythia8b->readString("PartonLevel:FSR = off");

    //Background
    //pythia8b->readString("HardQCD:all = on"); 
    //pythia8b->readString("PhaseSpace:pTHatMin  = 200");
    //pythia8b->readString("PhaseSpace:pTHatMax  = 220");

   pythia8b->readString("Beams:idA = 2212");
   pythia8b->readString("Beams:idB = 2212");
   pythia8b->readString("Beams:eCM = 13000");
   pythia8b->init();

   analysis1->SetOutName(outName.c_str());
   analysis1->Begin();

   Pythia8::Pythia* pythia_MB = new Pythia8::Pythia();
   pythia_MB->readString("Random:setSeed = on");
   ss.clear(); ss.str(""); ss << "Random:seed = " << seed+1;
   cout << ss.str() << endl;
   pythia_MB->readString(ss.str());
   pythia_MB->readString("SoftQCD:nonDiffractive = on");
   pythia_MB->readString("HardQCD:all = off");
   pythia_MB->readString("PhaseSpace:pTHatMin  = .1");
   pythia_MB->readString("PhaseSpace:pTHatMax  = 20000");
   pythia_MB->readString("Beams:idA = 2212");
   pythia_MB->readString("Beams:idB = 2212");
   pythia_MB->readString("Beams:eCM = 13000");
   pythia_MB->init();

   cout << "running on " << nEvents << " events " << endl;
   for (Int_t iev = 0; iev < nEvents; iev++) {
     if (iev%100==0) cout << iev << " " << nEvents << endl;
     analysis1->AnalyzeEvent(iev, pythia8b, pythia_MB, npu);//200);
   }
   analysis1->End();
   delete analysis1;

   pythia8b->stat();
   
   // that was it
   delete pythia8b;
   delete pythia_MB;
   //delete analysis;
   //delete analysis2;
   return 0;
}
