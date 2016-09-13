#include<string>
#include<set>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <random>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "Math/Vector4D.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSRecoJet.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSSamplingSection.hh"

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

using namespace std;

int main(int argc, char** argv) {
    
  std::cout << "argc: " << argc << std::endl;
  std::cout << "TChain files: " << argv[1] << std::endl;
  std::cout << "Output file tag: " << argv[2] << std::endl;
  TChain* tree = new TChain("HGCSSTree");
  tree->Add(argv[1]);

  std::vector<HGCSSSimHit> * hitVec = 0;
  std::vector<HGCSSGenParticle> * hadVec = 0;
  HGCSSEvent *ev = 0;
  tree->SetBranchAddress("HGCSSEvent", &ev);
  tree->SetBranchAddress("HGCSSSimHitVec", &hitVec);
  tree->SetBranchAddress("HGCSSHadAction", &hadVec);
	
  unsigned nEvts = tree->GetEntries();

  char outFileName[256];
  sprintf(outFileName,"simpleHCalDigis_%s.root",argv[2]);
  TFile hfile(outFileName,"RECREATE");
  TTree t1("hcal_digi", "");

  // const int nL = 25;
  // const int hcalL1 = 42;

  // Int_t nLayers,layerPEs[nL],layerNum[nL],layerSensDep[nL];
  // Float_t layerHFlux[nL],layerNFlux[nL],layerEFlux[nL],layerGFlux[nL],layerMFlux[nL];
	
  const int nLh = 15;
  const int nLe = 0;
  int event,seedx,seedy,seedz; 
  int nLayersH = nLh; 
  int nLayersE = nLe; 
  int nhitsE, nhitsH;
  int nHadrons;
  float henergy[nLh],hcalLayerPEs[nLh],hcalLayerNum[nLh];
  float eenergy[nLe];
  // std::vector<float> allHits_e;
  float allHits_e[500];

  float hpar_mass[100];
  //float hpar_px[100], hpar_py[100], hpar_pz[100];
  //float hpar_x[100], hpar_y[100], hpar_z[100];
  float hpar_id[100];
	
  t1.Branch("event",&event,"event/I");
  t1.Branch("seedx",&seedx,"seedx/I");
  t1.Branch("seedy",&seedy,"seedy/I");
  t1.Branch("seedz",&seedz,"seedz/I");

  t1.Branch("nhitsE",&nhitsE,"nhitsE/I");
  t1.Branch("nhitsH",&nhitsH,"nhitsH/I");
  t1.Branch("nLayersE", &nLayersE, "nLayersE/I");
  t1.Branch("eenergy", &eenergy, "eenergy[nLayersE]/F");
  t1.Branch("nLayersH", &nLayersH, "nLayersH/I");
  t1.Branch("henergy", &henergy, "henergy[nLayersH]/F");
  t1.Branch("hcalLayerPEs", &hcalLayerPEs, "hcalLayerPEs[nLayersH]/F");
  t1.Branch("hcalLayerNum", &hcalLayerNum, "hcalLayerNum[nLayersH]/F");

  t1.Branch("allHits_e", &allHits_e, "allHits_e[500]/F");
  t1.Branch("nHadrons", &nHadrons, "nHadrons/I");
  t1.Branch("hpar_mass", &hpar_mass, "hpar_mass[100]/F");
  //t1.Branch("hpar_px", &hpar_px, "hpar_px[100]/F");
  //t1.Branch("hpar_py", &hpar_py, "hpar_py[100]/F");
  //t1.Branch("hpar_pz", &hpar_pz, "hpar_pz[100]/F");
  //t1.Branch("hpar_x", &hpar_x, "hpar_x[100]/F");
  //t1.Branch("hpar_y", &hpar_y, "hpar_y[100]/F");
  //t1.Branch("hpar_z", &hpar_z, "hpar_z[100]/F");
  t1.Branch("hpar_id", &hpar_id, "hpar_id[100]/F");

  //t1.Branch("layerNum", &layerNum, "layerNum[nLayers]/I");
  //t1.Branch("layerSensDep", &layerSensDep, "layerSensDep[nLayers]/I");
  //t1.Branch("layerHFlux", &layerHFlux, "layerHFlux[nLayers]/F");
  //t1.Branch("layerNFlux", &layerNFlux, "layerNFlux[nLayers]/F");
  //t1.Branch("layerEFlux", &layerEFlux, "layerEFlux[nLayers]/F");
  //t1.Branch("layerGFlux", &layerGFlux, "layerGFlux[nLayers]/F");
  //t1.Branch("layerMFlux", &layerMFlux, "layerMFlux[nLayers]/F");

  int firstHCALLayer = 2; 
  int firstECALLayer = 1; 

  std::cout << "number of Evts = " << nEvts << std::endl;
	
  int nhitsTrue = 0;
  for (unsigned ievt(0); ievt < nEvts; ++ievt) { //loop on entries
    tree->GetEntry(ievt);
	
    // event info
    event = ev->eventNumber();
    seedx = ev->seeds().x();
    seedy = ev->seeds().y();
    seedz = ev->seeds().z();

    nhitsTrue = 0;
    nhitsE = 0;
    nhitsH = 0;
    // allHits_e.clear();
    for (int i = 0; i < 500;++i) allHits_e[i] = 0.;
    for (int i = 0; i < 100;++i){
      hpar_mass[i] = 0.;
      //hpar_px[i] = 0.;
      //hpar_py[i] = 0.;
      //hpar_pz[i] = 0.;
      //hpar_x[i] = 0.;
      //hpar_y[i] = 0.;
      //hpar_z[i] = 0.;
      hpar_id[i] = 0.;
    }


    // initialize
    for (int iL(0); iL < nLayersE; ++iL){ 
      eenergy[iL] = 0;
    }
    for (int iL(0); iL < nLayersH; ++iL){ 
      henergy[iL] = 0;
    }

    for (unsigned ihit(0); ihit < hitVec->size(); ++ihit){

      // ECAL
      if (hitVec->at(ihit).energy() > 0 && hitVec->at(ihit).layer() < nLayersE && hitVec->at(ihit).layer() > 0){ 
				
	int necalLayer = hitVec->at(ihit).layer() - firstECALLayer;
	eenergy[necalLayer] += hitVec->at(ihit).energy();

	allHits_e[nhitsTrue] = hitVec->at(ihit).energy();

	nhitsE++;
	nhitsTrue++;
      }

      // HCAL
      if (hitVec->at(ihit).energy() > 0 && hitVec->at(ihit).layer() >firstECALLayer+nLayersE ){ 

	int nhcalLayer = hitVec->at(ihit).layer() - firstHCALLayer;
	henergy[nhcalLayer] += hitVec->at(ihit).energy();

	allHits_e[nhitsTrue] = hitVec->at(ihit).energy();
				
	nhitsH++;
	nhitsTrue++;
				
      }
      for (int iL(0); iL < nLayersH; ++iL){ 

	//// DIGITIZED INFORMATION
	float MeVperMIP = 1.40;
	float PEperMIP = 13.5*6./4.;
	float depEnergy = henergy[iL];
	float meanPE = depEnergy/MeVperMIP*PEperMIP;
	std::default_random_engine generator;
	std::poisson_distribution<int> distribution(meanPE);
	hcalLayerPEs[iL] = distribution(generator);
	hcalLayerNum[iL] = iL;
			
      }

    }

    nHadrons = hadVec->size();
    for (unsigned ihad(0); ihad < hadVec->size(); ++ihad ){
      hpar_mass[ihad] = hadVec->at(ihad).mass();
      //hpar_px[ihad] = hadVec->at(ihad).px();
      //hpar_py[ihad] = hadVec->at(ihad).py();
      //hpar_pz[ihad] = hadVec->at(ihad).pz();
      //hpar_x[ihad] = hadVec->at(ihad).x();
      //hpar_y[ihad] = hadVec->at(ihad).y();
      //hpar_z[ihad] = hadVec->at(ihad).z();
      hpar_id[ihad] = hadVec->at(ihad).pdgid();		
    }

    t1.Fill();

    // 	//// DIGITIZED INFORMATION
    // 	float MeVperMIP = 1.40;
    // 	float PEperMIP = 13.5*6./4.;
    // 	float depEnergy = sec.sensDep();
    // 	float meanPE = depEnergy/MeVperMIP*PEperMIP;
    // 	std::default_random_engine generator;
    // 	std::poisson_distribution<int> distribution(meanPE);
    // 	layerPEs[j-hcalL1] = distribution(generator);
    // 	layerNum[j-hcalL1] = j-hcalL1;
    // 	layerSensDep[j-hcalL1] = sec.sensDep();
    // 	/*
    // 	std::cout << "- - - - - " << j-42 << " - - - - - - - " << std::endl;
    // 	std::cout << "depEnergy: " << depEnergy << std::endl;
    // 	std::cout << "meanPE: " << meanPE << std::endl;
    // 	std::cout << "PE: " << layerPEs[j-43] << std::endl;
    // 	*/

    // }

    // t1.Fill();
  }

  t1.Write();
  hfile.Close();
  return 1;
}
