#include "HGCSSCalibration.hh"
#include <sstream>
#include <iostream>
#include <cmath>

HGCSSCalibration::~HGCSSCalibration() {}

double HGCSSCalibration::MeVToMip(const unsigned layer,
		const bool absWeight) const {
	if (layer < theDetector().nLayers())
		return theDetector().subDetectorByLayer(layer).mipWeight
				* (absWeight ?
						theDetector().subDetectorByLayer(layer).absWeight : 1.0);
	return 1;
}

