#ifndef _samplingsection_hh_
#define _samplingsection_hh_

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Colour.hh"
#include "G4Track.hh"
//#include "EventAction.hh"
#include "G4RunManager.hh"

#include <iomanip>
#include <vector>
#include "G4SiHit.hh"
//class EventAction;
class SamplingSection {
public:
	//CTOR
	SamplingSection(std::vector<std::pair <G4double,std::string>> iEle);


	//DTOR
	~SamplingSection() {};

	inline void setNumberOfSectors(const unsigned nSec) {
		n_sectors = nSec;
		sublayer_vol.clear();
		for (unsigned ie(0); ie<n_elements*n_sectors; ++ie) {
			sublayer_vol.push_back(0);
		}
	};

	//
	std::pair<G4bool,G4bool> add(G4double parentKE, G4double depositRawE,G4VPhysicalVolume* vol,G4Track *lTrack,const G4ThreeVector & position, G4int pdgId, G4bool isForward);

	inline bool isSensitiveElement(const unsigned & aEle) {
		if (aEle < n_elements &&
				(ele_name[aEle] == "Si" || ele_name[aEle] == "Scintillator")
		) return true;
		return false;
	};

	inline unsigned getSensitiveLayerIndex(std::string astr) {
		if (astr.find("_")== astr.npos) return 0;
		size_t pos = astr.find("phys");
		if (pos != astr.npos && pos>1) {
			unsigned idx = 0;
			std::istringstream(astr.substr(pos-1,1))>>idx;
			return idx;
		}
		return 0;
	};
	template <typename T>
	inline void clearVec( T & t ) {
	    T tmp;
	    using std::swap;
	    t.swap( tmp );
	}

	inline G4Colour g4Colour(const unsigned & aEle) {
		if (isSensitiveElement(aEle)) return G4Colour::Red();
		if (ele_name[aEle] == "Cu") return G4Colour::Black();
		if (isAbsorberElement(aEle)) return G4Colour::Gray();
		if (ele_name[aEle] == "PCB") return G4Colour::Blue();
		if (ele_name[aEle] == "Air") return G4Colour::Cyan();
		return G4Colour::Yellow();
	};

	inline bool isAbsorberElement(const unsigned & aEle) {
		if (aEle < n_elements &&
				(
						ele_name[aEle] == "Pb" || ele_name[aEle] == "Cu" ||
						ele_name[aEle] == "W" || ele_name[aEle] == "Brass" ||
						ele_name[aEle] == "Fe" || ele_name[aEle] == "Steel" ||
						ele_name[aEle] == "SSteel" || ele_name[aEle] == "Al" ||
						ele_name[aEle] == "WCu" || ele_name[aEle] == "NeutMod"
				)
		) return true;
		return false;
	};

	//reset
	inline void resetCounters()
	{
		sublayer_RawDep.clear();
		sublayer_RawDep.resize(n_elements,0);

		sens_gamKinFlux.clear();
		sens_neutronKinFlux.clear();
		sens_eleKinFlux.clear();
		sens_muKinFlux.clear();
		sens_hadKinFlux.clear();

		sens_gamKinFlux.resize(n_elements,0);
		sens_neutronKinFlux.resize(n_elements,0);
		sens_eleKinFlux.resize(n_elements,0);
		sens_muKinFlux.resize(n_elements,0);
		sens_hadKinFlux.resize(n_elements,0);

		Etracks.clear();
		Gtracks.clear();
		Htracks.clear();
		Mtracks.clear();
		Ntracks.clear();

		Etracks.resize(n_sens_elements);
		Gtracks.resize(n_sens_elements);
		Htracks.resize(n_sens_elements);
		Mtracks.resize(n_sens_elements);
		Ntracks.resize(n_sens_elements);

		//clearVec(sublayer_RawDep);
		for (unsigned idx(0); idx<n_sens_elements; ++idx) {
			//sens_HitVec[idx].clear();
			clearVec(sens_HitVec[idx]);
		}

	}

	inline unsigned getLayer(std::string astr) {
		
		size_t pos = astr.find("phys");
		unsigned num = 0;
		if (astr.find("_")== astr.npos) {
			std::string truncated = astr.substr(0,pos);
			size_t last_index = truncated.find_last_not_of("0123456789");
			num = std::atoi(truncated.substr(last_index+1).c_str());

		}
		else{
			pos = astr.find("_");
			std::string truncated = astr.substr(0,pos);
			size_t last_index = truncated.find_last_not_of("0123456789");
			num = std::atoi(truncated.substr(last_index+1).c_str());
		}

		return num;

	};

	inline G4double getKinGam() {
		double val = 0;
		for (unsigned ie(0); ie < n_sens_elements; ++ie) { val += sens_gamKinFlux[ie]; }
		return val;
	};
	inline G4double getKinNeutron() {
		double val = 0;
		for (unsigned ie(0); ie < n_sens_elements; ++ie) { val += sens_neutronKinFlux[ie]; }
		return val;
	};
	inline G4double getKinEle() {
		double val = 0;
		for (unsigned ie(0); ie < n_sens_elements; ++ie) { val += sens_eleKinFlux[ie]; }
		return val;
	};
	inline G4double getKinHadron() {
		double val = 0;
		for (unsigned ie(0); ie < n_sens_elements; ++ie) { val += sens_hadKinFlux[ie]; }
		return val;
	};
	inline G4double getKinMuon() {
		double val = 0;
		for (unsigned ie(0); ie < n_sens_elements; ++ie) { val += sens_muKinFlux[ie]; }
		return val;
	};	
	const G4SiHitVec & getSiHitVec(const unsigned & idx) const;

	G4double getTotalSensE();

	//
	void report(bool header=false);

	//members
	unsigned n_elements;
	unsigned n_sectors;
	unsigned n_sens_elements;

	//EventAction *eventAction_;

	std::vector<G4double> sublayer_thick;
	std::vector<std::string> ele_name;
	std::vector<G4double> sublayer_X0;
	std::vector<G4double> sublayer_dEdx;
	std::vector<G4double> sublayer_L0;
	std::vector<G4double> sublayer_RawDep;
	std::vector<G4double> sublayer_PrimaryDep;
	std::vector<G4double> sublayer_dl;
	std::vector<G4VPhysicalVolume*> sublayer_vol;

	std::vector<G4double> sens_gamKinFlux, sens_neutronKinFlux, sens_eleKinFlux, sens_muKinFlux, sens_hadKinFlux;
	std::vector<std::vector<unsigned int>> Etracks,Gtracks,Mtracks,Ntracks,Htracks;

	std::vector<G4SiHitVec> sens_HitVec;
	G4double Total_thick;
	bool hasScintillator;
	unsigned sens_HitVec_size_max;

};

#endif
