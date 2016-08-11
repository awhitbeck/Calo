#include "HGCSSDetector.hh"
#include <iostream>

HGCSSDetector &theDetector() {

	static HGCSSDetector lDet;
	static bool firstDet = true;
	if (firstDet)
		std::cout << " -- Created detector static object." << std::endl;
	firstDet = false;
	return lDet;
}

void HGCSSDetector::buildDetector(const unsigned versionNumber) {

	reset();
	initialiseIndices(versionNumber);
	HGCSSSubDetector ECAL;
	ECAL.type       = DetectorEnum::ECAL;
	ECAL.name       = "ECAL";
	ECAL.layerIdMin = indices_[0];
	ECAL.layerIdMax = indices_[1];
	ECAL.mipWeight  = 1. / 0.0548; //mip for 200um si
	ECAL.absWeight  = 1.; //ratio of abs dedx
	ECAL.gevWeight  = 1.0;
	ECAL.gevOffset  = 0.0;
	ECAL.isSi       = true;
	ECAL.radiusLim  = 0;

	if (ECAL.nLayers() > 0)
		theDetector().addSubdetector(ECAL);

	HGCSSSubDetector HCAL;
	HCAL.type       = DetectorEnum::HCAL;
	HCAL.name       = "HCAL";
	HCAL.layerIdMin = indices_[1];
	HCAL.layerIdMax = indices_[2];
	HCAL.mipWeight  = 1. / 0.0548; //mip for 200um si
	HCAL.absWeight  = 8.001 / 5.848; //ratio of abs dedx
        HCAL.gevWeight  = 1.0;
	HCAL.gevOffset  = 0.0;
	HCAL.isScint    = true;
	HCAL.radiusLim  = 0;

	if (HCAL.nLayers() > 0) 
		theDetector().addSubdetector(HCAL);
        
	finishInitialisation();
}

const HGCSSSubDetector &HGCSSDetector::subDetectorByLayer(const unsigned aLayer) {

	unsigned section = getSection(aLayer);
	return subdets_[section];
}

unsigned HGCSSDetector::getSection(const unsigned aLayer) const {

	if (aLayer >= nLayers_) {
		std::cerr << " -- Error ! Trying to access layer " << aLayer
				<< " outside of range. nLayers = " << nLayers_ << std::endl;
		exit(1);
	}
	return section_[aLayer];
}

void HGCSSDetector::addSubdetector(const HGCSSSubDetector &adet) {

	subdets_.push_back(adet);
	enumMap_[adet.type] = subdets_.size() - 1;
}

void HGCSSDetector::finishInitialisation() {

	nSections_ = subdets_.size();
	const unsigned lastEle = indices_.size() - 1;
	nLayers_ = indices_[lastEle];

	unsigned secIndex[lastEle];
	unsigned iS(0);
	for (unsigned i(0); i < lastEle; i++) {
		secIndex[i] = iS;
		if (indices_[i] < indices_[i + 1])
			iS += 1;
	}

	//initialise layer-section conversion
	section_.resize(nLayers_, 0);
	for (unsigned iL(0); iL < nLayers_; ++iL) {
		for (unsigned i(0); i < lastEle; ++i) {
			if (iL >= indices_[i] && iL < indices_[i + 1])
				section_[iL] = secIndex[i];
		}
	}
	printDetector(std::cout);
}

const HGCSSSubDetector & HGCSSDetector::subDetectorByEnum(DetectorEnum adet) {

	if (enumMap_.find(adet) == enumMap_.end()) {
		std::cerr
				<< " -- Error ! Trying to access subdetector enum not present in this detector: "
				<< adet << std::endl;
		exit(1);
	}
	return subdets_[enumMap_[adet]];
}

void HGCSSDetector::reset() {

	subdets_.clear();
	enumMap_.clear();
	indices_.clear();
	section_.clear();
}

void HGCSSDetector::printDetector(std::ostream &aOs) const {

	std::cout << " -------------------------- " << std::endl
			<< " -- Detector information -- " << std::endl
			<< " -------------------------- " << std::endl << " - nSections = "
			<< nSections_ << std::endl << " - nLayers = " << nLayers_
			<< std::endl << " - detNames = ";
	for (unsigned i(0); i < nSections_; ++i) {
		std::cout << " " << detName(i);
	}
	std::cout << std::endl;
	std::cout << " -------------------------- " << std::endl;
}
