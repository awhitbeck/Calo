#ifndef _hgcssrecohit_hh_
#define _hgcssrecohit_hh_

#include <iomanip>
#include <vector>
#include <cmath>
#include "Rtypes.h"
#include <sstream>
#include "TMath.h"

#include "HGCSSSimHit.hh"
#include "HGCSSGeometryConversion.hh"
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"

class HGCSSRecoHit {

public:
	HGCSSRecoHit() :
			energy_(0), adcCounts_(0), xpos_(0), ypos_(0), zpos_(0), layer_(0),
			noiseFrac_(0) {
	}
	;

	HGCSSRecoHit(const HGCSSSimHit &aSimHit, const bool isScintillator,
			const HGCSSGeometryConversion &aGeom);

	~HGCSSRecoHit() {};

	inline double energy() const {
		return energy_;
	}
	;

	double eta() const;
	double theta() const;
	double phi() const;

	inline ROOT::Math::XYZPoint position() const {
		return ROOT::Math::XYZPoint(xpos_ / 10., ypos_ / 10., zpos_ / 10.);
	}
	;

	inline double E() const {
		return energy_;
	}
	;

	inline double pt() const {
		return energy_ / cosh(eta());
	}
	;

	inline double px() const {
		return pt() * cos(phi());
	}
	;

	inline double py() const {
		return pt() * sin(phi());
	}
	;

	inline double pz() const {
		return pt() * sinh(eta());
	}
	;

	inline void energy(const double &energy) {
		energy_ = energy;
	}
	;

	inline void x(const double &pos) {
		xpos_ = pos;
	}
	;

	inline void y(const double &pos) {
		ypos_ = pos;
	}
	;

	inline void z(const double &pos) {
		zpos_ = pos;
	}
	;

	inline unsigned adcCounts() const {
		return adcCounts_;
	}
	;

	inline void adcCounts(const unsigned &adcCounts) {
		adcCounts_ = adcCounts;
	}
	;

	inline unsigned layer() const {
		return layer_;
	}
	;

	inline void layer(const unsigned &layer) {
		layer_ = layer;
	}
	;

	inline double noiseFraction() const {
		return noiseFrac_;
	}
	;

	inline void noiseFraction(const double &aFrac) {
		noiseFrac_ = aFrac;
	}
	;

	void Add(const HGCSSSimHit &aSimHit);

	inline double get_x() const {
		return xpos_;
	}
	;

	inline double get_y() const {
		return ypos_;
	}
	;

	inline double get_z() const {
		return zpos_;
	}
	;

	void Print(std::ostream &aOs) const;

private:

	double energy_;
	unsigned adcCounts_;
	double xpos_;
	double ypos_;
	double zpos_;
	unsigned layer_;
	double noiseFrac_;

ClassDef(HGCSSRecoHit,1);

};

typedef std::vector<HGCSSRecoHit> HGCSSRecoHitVec;

#endif
