#include "G4VPhysicalVolume.hh"
#include "SamplingSection.hh"
//

SamplingSection::SamplingSection(std::vector<std::pair <G4double,std::string>> iEle) {

	std::vector<G4double> aThicknessVec;std::vector<std::string> aMaterialVec;
	for (unsigned i = 0; i < iEle.size(); i++)
	{
		aThicknessVec.push_back(iEle.at(i).first);
		aMaterialVec.push_back(iEle.at(i).second);

	}
	Total_thick = 0;
	n_sens_elements=0;
	n_elements=0;
	n_sectors=0;
	sublayer_thick.clear();
	ele_name.clear();
	sublayer_X0.clear();
	sublayer_L0.clear();
	sublayer_vol.clear();
	hasScintillator = false;
	for (unsigned ie(0); ie<aThicknessVec.size(); ++ie) {
		//consider only material with some non-0 width...
		if (aThicknessVec[ie]>0) {
			sublayer_thick.push_back(aThicknessVec[ie]);
			ele_name.push_back(aMaterialVec[ie]);
			if (aMaterialVec[ie]== "Scintillator") hasScintillator = true;
			sublayer_X0.push_back(0);
			sublayer_dEdx.push_back(0);
			sublayer_L0.push_back(0);
			sublayer_vol.push_back(0);
			Total_thick+=aThicknessVec[ie];
			++n_elements;
			//the following method check the total size...
			//so incrementing first.
			if (isSensitiveElement(n_elements-1)) {
				G4SiHitVec lVec;
				sens_HitVec.push_back(lVec);
				++n_sens_elements;
			}
		}
	}
	sens_HitVec_size_max = 0;

	resetCounters();
	//eventAction_ =(EventAction*) G4RunManager::GetRunManager()->GetUserEventAction();
	//std::cout << " -- End of sampling section initialisation. Input " << aThicknessVec.size() << " elements, constructing " << n_elements << " elements with " << n_sens_elements << " sensitive elements." << std::endl;

}

std::pair<G4bool,G4bool> SamplingSection::add(G4double parentKE, G4double depositRawE, G4VPhysicalVolume* vol,G4Track *lTrack, const G4ThreeVector & position, G4int pdgId, G4bool isForward) 
{
	std::string lstr = vol->GetName();
	bool breakSwitch = false;
	bool isSens = false;
	for (unsigned ie(0); ie < n_elements * n_sectors; ++ie) {
		if (breakSwitch) return std::make_pair(breakSwitch,isSens);
		if (sublayer_vol[ie] && lstr == sublayer_vol[ie]->GetName()) {
			
			breakSwitch = true;
			unsigned idx = getSensitiveLayerIndex(lstr);
			unsigned eleidx = ie % n_elements;
			
			if (isSensitiveElement(ie)){
				isSens = isSensitiveElement(eleidx);
				sublayer_RawDep[eleidx] += depositRawE;
				G4SiHit lHit;
				lHit.parentKE = lTrack->GetKineticEnergy()* MeV;
				lHit.energyDep = depositRawE;

				lHit.hit_x = position.x();
				lHit.hit_y = position.y();
				lHit.hit_z = position.z();

				lHit.layer = getLayer(lstr);
				lHit.trackId = lTrack->GetTrackID();
				lHit.pdgId = lTrack->GetDefinition()->GetPDGEncoding();

				sens_HitVec[idx].push_back(lHit);

				// sampling section information
				if (abs(pdgId) == 22){
					if (isForward){
						//Prevent twice counting any photon in any given sub-layer.
						unsigned int trackLoc  = std::find(Gtracks[idx].begin(),Gtracks[idx].end(), lHit.trackId) - Gtracks[idx].begin();
						if (trackLoc == Gtracks[idx].size()){
							sens_gamKinFlux[idx] += parentKE;
							Gtracks[idx].push_back(lHit.trackId);
						}
					}
				}
				else if (abs(pdgId) == 11){
					if (isForward){
						unsigned int trackLoc  = std::find(Etracks[idx].begin(),Etracks[idx].end(), lHit.trackId) - Etracks[idx].begin();
						if (trackLoc == Etracks[idx].size()){
							sens_eleKinFlux[idx] += parentKE;
							Etracks[idx].push_back(lHit.trackId);
						}
					}
				}				
				else if (abs(pdgId) == 13){
					if (isForward){
						unsigned int trackLoc  = std::find(Mtracks[idx].begin(),Mtracks[idx].end(), lHit.trackId) - Mtracks[idx].begin();
						if (trackLoc == Mtracks[idx].size()){
							sens_muKinFlux[idx] += parentKE;
							Mtracks[idx].push_back(lHit.trackId);
						}
					}
				}				
				else if (abs(pdgId) == 2112){
					if (isForward){
						unsigned int trackLoc  = std::find(Ntracks[idx].begin(),Ntracks[idx].end(), lHit.trackId) - Ntracks[idx].begin();
						if (trackLoc == Ntracks[idx].size()){
							sens_neutronKinFlux[idx] += parentKE;
							Ntracks[idx].push_back(lHit.trackId);
						}
					}
				}				
				else{
					if (isForward){
						unsigned int trackLoc  = std::find(Htracks[idx].begin(),Htracks[idx].end(), lHit.trackId) - Htracks[idx].begin();
						if (trackLoc == Htracks[idx].size()){
							sens_hadKinFlux[idx] += parentKE;
							Htracks[idx].push_back(lHit.trackId);
						}
					}
				}				

			}
		} //if in right material
	} //loop on available materials

	return std::make_pair(breakSwitch,isSens);
}

//
void SamplingSection::report(bool header) {
	if (header)
		G4cout
				<< "E/[MeV]\t  Si\tAbsorber\tTotal\tSi g frac\tSi e frac\tSi mu frac\tSi had frac\tSi <t> \t nG4SiHits"
				<< G4endl;
		G4cout << std::setprecision(3) << "\t  " << getTotalSensE()
		<< G4endl;
	}

const G4SiHitVec & SamplingSection::getSiHitVec(const unsigned & idx) const {
	return sens_HitVec[idx];
}
G4double SamplingSection::getTotalSensE() {
	double etot = 0;
	for (unsigned ie(0); ie < n_elements; ++ie) {
		if (isSensitiveElement(ie))
			etot += sublayer_RawDep[ie];
	}
	return etot;
}
