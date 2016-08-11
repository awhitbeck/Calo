#ifndef HGCSSDigitisation_h
#define HGCSSDigitisation_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "TRandom3.h"
#include "TH2D.h"
#include "HGCSSDetector.hh"

class HGCSSDigitisation {

public:

	HGCSSDigitisation() :
		seed_(0), npe_(9), crossTalk_(0.25), ipXtalk_(0.025), nTotal_(1156), sigmaPix_(3) {

		rndm_.SetSeed(seed_);
		//noise_[DetectorEnum::ECAL]      = 0.12;
		//noise_[DetectorEnum::HCAL]      = 0.12;
		maxADC_[DetectorEnum::ECAL]       = 65535; // 16-bit
		maxADC_[DetectorEnum::HCAL]       = 65535;
		mipToADC_[DetectorEnum::ECAL]     = 8.2; //ADC per mips.
		mipToADC_[DetectorEnum::HCAL]     = 10;
		gainSmearing_[DetectorEnum::ECAL] = 0.02; //2% intercalibration
		gainSmearing_[DetectorEnum::HCAL] = 0.02;

		Print(std::cout);

	}
	;

	~HGCSSDigitisation() {};

	inline void setIntercalibrationFactor(const unsigned icFactor) {
		gainSmearing_[DetectorEnum::ECAL] = icFactor / 100.;
		gainSmearing_[DetectorEnum::HCAL] = icFactor / 100.;
	}
	;

	inline void setRandomSeed(const unsigned aSeed) {
		seed_ = aSeed;
		rndm_.SetSeed(seed_);
	}
	;

	inline void setNpe(const unsigned aNpe) {
		npe_ = aNpe;
	}
	;

	inline void setCrossTalk(const double &aCrossTalk) {
		crossTalk_ = aCrossTalk;
	}
	;

	inline void setIPCrossTalk(const double &aCrossTalk) {
		ipXtalk_ = aCrossTalk;
	}
	;

	inline void setNTotalPixels(const unsigned &aN) {
		nTotal_ = aN;
	}
	;

	inline void setSigmaPix(const unsigned aSigma) {
		sigmaPix_ = aSigma;
	}
	;

	inline void setNoise(const unsigned &alay, const double &aNoise) {
		noise_[alay] = aNoise;
	}
	;

	inline void setMipToADC(DetectorEnum adet, const double &aMipToADC) {
		mipToADC_[adet] = aMipToADC;
	}
	;

	inline void setMaxADC(DetectorEnum adet, const double &aMaxADC) {
		maxADC_[adet] = aMaxADC;
	}
	;

	inline void setGainSmearing(DetectorEnum adet, const double &aVal) {
		gainSmearing_[adet] = aVal;
	}
	;

	unsigned nRandomPhotoElec(const double &aMipE);

	unsigned nPixels(const double &aMipE);

	unsigned positiveRandomGaus(const unsigned &mean);

	double mipCor(const double &aMipE, const double &posx,
			const double &posy, const double &posz);

	double digiE(const double &aMipE, TH2F* &p_pixvspe, TH1F* &p_npixels,
			TH1F* &p_npixelssmeared, TH2F* &p_outvsnpix);

	double digiE(const double &aMipE);

	double ipXtalk(const std::vector<double> &aSimEvec);

	void addNoise(double &aDigiE, const unsigned &alay);

	unsigned adcConverter(double eMIP, DetectorEnum adet);

	double adcToMIP(const unsigned acdCounts, DetectorEnum adet, const bool smear = true);

	double MIPtoGeV(const HGCSSSubDetector &adet, const double &aMipE);

	double sumBins(const std::vector<TH2D*> &aHistVec, const double &aMipThresh);

	void Print(std::ostream &aOs) const;

private:

	unsigned seed_;
	unsigned npe_;
	double   crossTalk_;
	double   ipXtalk_;
	unsigned nTotal_;
	unsigned sigmaPix_;
	TRandom3 rndm_;
	std::map<DetectorEnum, unsigned> mipToADC_;
	std::map<DetectorEnum, unsigned> maxADC_;
	std::map<DetectorEnum, double>   gainSmearing_;
	std::map<unsigned, double>       noise_;

};

#endif
