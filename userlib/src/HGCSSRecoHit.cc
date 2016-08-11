#include "HGCSSRecoHit.hh"

#include <iomanip>
#include <cmath>
#include <stdlib.h>

HGCSSRecoHit::HGCSSRecoHit(const HGCSSSimHit &aSimHit,
		const bool isScintillator, const HGCSSGeometryConversion &aGeom) {

	energy_    = aSimHit.energy();
	adcCounts_ = 0;
	zpos_      = aSimHit.get_z();

	layer_     = aSimHit.layer();
	noiseFrac_ = 0;

	std::pair<double, double> xy = aSimHit.get_xy(isScintillator, aGeom);
	xpos_ = xy.first;
	ypos_ = xy.second;

}

double HGCSSRecoHit::theta() const {
	return 2 * atan(exp(-1. * eta()));
}

double HGCSSRecoHit::eta() const {
	return position().eta();
}

double HGCSSRecoHit::phi() const {
	return position().phi();
}

void HGCSSRecoHit::Add(const HGCSSSimHit &aSimHit) {
	energy_ += aSimHit.energy();
}

void HGCSSRecoHit::Print(std::ostream &aOs) const {
	aOs << "====================================" << std::endl << " = Layer "
			<< layer_
			<< std::endl << " = Energy " << energy_ << " noiseFrac "
			<< noiseFrac_ << std::endl << " = Digi E " << adcCounts_
			<< " adcCounts." << std::endl
			<< "====================================" << std::endl;

}
