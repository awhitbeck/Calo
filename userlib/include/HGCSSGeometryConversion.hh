#ifndef HGCSSGeometryConversion_h
#define HGCSSGeometryConversion_h

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "TH2D.h"
#include "TH2Poly.h"
#include "TMath.h"
#include "HGCSSDetector.hh"

class HGCSSGeometryConversion {

public:
	HGCSSGeometryConversion() {};

	HGCSSGeometryConversion(const unsigned model, const double cellsize);

	~HGCSSGeometryConversion();

	TH2Poly *hexagonMap() {
		static TH2Poly hc;
		return &hc;
	}
	;

        TH2Poly *squareMap() {
                static TH2Poly hsq;
                return &hsq;
        }
        ;

        std::map<int, std::pair<double, double> > hexaGeom;
        std::map<int, std::pair<double, double> > squareGeom;

        void initialiseSquareMap(const double xymin, const double side);

        void initialiseSquareMap(TH2Poly *map, const double xymin,
                        const double side, bool print);

	void initialiseHoneyComb(const double xymin, const double side);

	void initialiseHoneyComb(TH2Poly *map, const double xymin, const double side, bool print);

	void fillXY(TH2Poly* hist, std::map<int, std::pair<double, double> > &geom);

	void setGranularity(const std::vector<unsigned> &granul);

	unsigned getGranularity(const unsigned aLayer, const HGCSSSubDetector &adet);

	inline double getXYwidth() const {
		return width_;
	}
	;

	inline void setXYwidth(double width) {
		width_ = width;
	}
	;

	inline double cellSize() const {
		return cellSize_;
	}
	;

	double cellSize(const unsigned aLayer) const;

	double cellSizeInCm(const unsigned aLayer) const;

	void initialiseHistos(const bool recreate = false, std::string uniqStr = "",
			const bool print = true);

	void fill(std::vector<std::vector<int> > &filledBins, const DetectorEnum type, const unsigned newlayer, const double &weightedE, const double &posx, const double &posy);

	double sumBins(const std::vector<TH2Poly*> &aHistVec, const double &aMipThresh = 0.);

	void resetVector(std::vector<TH2Poly*> &aVec, std::string aVar,
			std::string aString, const HGCSSSubDetector &aDet,
			const unsigned nLayers, bool recreate = false, bool print = true);

	void deleteHistos(std::vector<TH2Poly*> &aVec);

	TH2Poly* get2DHist(const unsigned layer, std::string name);
        
	inline std::vector<TH2Poly*> & get2DEnergyVec(const DetectorEnum aDet) {
		return HistMapE_[aDet];
	}
	;

	inline std::vector<TH2Poly*> &get2DTimeVec(const DetectorEnum aDet) {
		return HistMapTime_[aDet];
	}
	;

	inline std::vector<TH2Poly*> &get2DNVec(const DetectorEnum aDet) {
		return HistMapN_[aDet];
	}
	;

private:

	void myHoneycomb(TH2Poly*map, Double_t xstart, Double_t ystart, Double_t a, // side length
			Int_t k,     // # hexagons in a column
			Int_t s);    // # columns

	double   width_;
	double   cellSize_;
	std::vector<unsigned> granularity_;
	unsigned model_;
	std::map<DetectorEnum, std::vector<TH2Poly*> > HistMapE_;
	std::map<DetectorEnum, std::vector<TH2Poly*> > HistMapTime_;
	std::map<DetectorEnum, std::vector<TH2Poly*> > HistMapN_;
	std::map<DetectorEnum, std::vector<double> >   avgMapN_;
	std::map<DetectorEnum, std::vector<double> >   avgMapE_;

};

#endif
