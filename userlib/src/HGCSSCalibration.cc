#include "HGCSSCalibration.hh"
#include <sstream>
#include <iostream>
#include <cmath>

HGCSSCalibration::HGCSSCalibration(const unsigned nSi) {

	vtx_x_ = 0;
	vtx_y_ = 0;
	vtx_z_ = 0;
	nSiLayers_ = nSi;
}

HGCSSCalibration::~HGCSSCalibration() {}

double HGCSSCalibration::MeVToMip(const unsigned layer,
		const bool absWeight) const {
	if (layer < theDetector().nLayers())
		return theDetector().subDetectorByLayer(layer).mipWeight
				* (absWeight ?
						theDetector().subDetectorByLayer(layer).absWeight : 1.0);
	return 1;
}

