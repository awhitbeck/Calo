#ifndef HGCSSCalibration_h
#define HGCSSCalibration_h

#include <string>
#include <vector>
#include "TH2D.h"

#include "HGCSSDetector.hh"

class HGCSSCalibration {

public:
	HGCSSCalibration() {};

	~HGCSSCalibration();

	inline void setVertex(const double &x, const double &y,	const double &z) {

		vtx_x_ = x;
		vtx_y_ = y;
		vtx_z_ = z;

	}
	;

	double MeVToMip(const unsigned layer, const bool absWeight = false) const;

private:

	double   vtx_x_;
	double   vtx_y_;
	double   vtx_z_;
};

#endif
