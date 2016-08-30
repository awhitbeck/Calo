#include "HGCSSGeometryConversion.hh"
#include <sstream>
#include <iostream>
#include <cmath>

HGCSSGeometryConversion::HGCSSGeometryConversion(const unsigned model,
		const double cellsize) {

	width_ = 200; //mm
	model_ = model;
	if (model == 0)
		width_ = 200;
	else if (model == 1)
		width_ = 500;
	else if (model == 2)
		width_ = 1500;
	else if (model == 3)
		width_ = 1000;

	cellSize_ = cellsize;

}

HGCSSGeometryConversion::~HGCSSGeometryConversion() {

	std::map<DetectorEnum, std::vector<TH2Poly*> >::iterator liter = HistMapE_.begin();
	for (; liter != HistMapE_.end(); ++liter) {

		deleteHistos(liter->second);
	}

	HistMapE_.clear();

        liter = HistMapN_.begin();
	for (; liter != HistMapN_.end(); ++liter) {

		deleteHistos(liter->second);
	}

	HistMapN_.clear();
	hexaGeom.clear();
}

void HGCSSGeometryConversion::initialiseSquareMap(const double xymin, const double side) {
        initialiseSquareMap(squareMap(), xymin, side, true);
        fillXY(squareMap(), squareGeom);
}

void HGCSSGeometryConversion::initialiseSquareMap(TH2Poly *map,
        const double xymin, const double side, bool print) {
        unsigned nx = static_cast<unsigned>(xymin * 2. / side);
        unsigned ny = nx;
        unsigned i, j;

        Double_t x1, y1, x2, y2;
        Double_t dx = side, dy = side;
        x1 = -1. * xymin;
        x2 = x1 + dx;

        for (i = 0; i < nx; i++) {
                y1 = -1. * xymin;
                y2 = y1 + dy;
                for (j = 0; j < ny; j++) {
                        map->AddBin(x1, y1, x2, y2);
                        y1 = y2;
                        y2 = y1 + dy;
                }
                x1 = x2;
                x2 = x1 + dx;
        }

        if (print) {
                std::cout << " -- Initialising squareMap with parameters: " << std::endl
                        << " ---- xymin = " << -1. * xymin << ", side = " << side
                        << ", nx = " << nx << ", ny=" << ny << std::endl;
        }

}

void HGCSSGeometryConversion::initialiseHoneyComb(const double width, const double side) {

	initialiseHoneyComb(hexagonMap(), width, side, true);
	fillXY(hexagonMap(), hexaGeom);
}

void HGCSSGeometryConversion::initialiseHoneyComb(TH2Poly *map,
		const double width, const double side, bool print) {
	//xstart,ystart,side length,

	// Center a cell at (x,y)=(0,0) and ensure coverage up to/past width/2 in all 4 directions,
	// assuming each cell is lying on a side.

	unsigned ncellwide = width / (2. * side);
	unsigned ny        = ncellwide + 1;
	unsigned nx        = ncellwide + 4;
	double   xstart    = -((double) ncellwide + 0.5) * side;
	double   ystart    = -((double) ncellwide + 1) * side * sqrt(3) / 2;

	if (print) {
		std::cout << " -- Initialising HoneyComb with parameters: " << std::endl
				<< " ---- (xstart,ystart) = (" << xstart << "," << ystart << ")"
				<< ", side = " << side << ", nx = " << nx << ", ny=" << ny
				<< std::endl;
	}
	//map->Honeycomb(-1.*xymin,-1.*xymin,side,nx,ny);
	myHoneycomb(map, xstart, ystart, side, ny, nx);
}

/// Bins the histogram using a honeycomb structure
/// 90 degree rotation, side up instead of corner up

void HGCSSGeometryConversion::myHoneycomb(TH2Poly* map, Double_t xstart,
		Double_t ystart, Double_t a,  // side length
		Int_t k,     // # hexagons in a column
		Int_t s)     // # columns
		{
	// Add the bins
	Double_t numberOfHexagonsInAColumn;
	Double_t x[6], y[6];
	Double_t xloop, yloop, ytemp;
	xloop = xstart;
	yloop = ystart + a * TMath::Sqrt(3) / 2.0;
	for (int sCounter = 0; sCounter < s; sCounter++) {

		ytemp = yloop; // Resets the temp variable

		// Determine the number of hexagons in that column
		if (sCounter % 2 == 0) {
			numberOfHexagonsInAColumn = k;
		} else {
			numberOfHexagonsInAColumn = k - 1;
		}

		for (int kCounter = 0; kCounter < numberOfHexagonsInAColumn;
				kCounter++) {

			// Go around the hexagon
			x[0] = xloop;
			y[0] = ytemp;
			x[1] = x[0] + a / 2.0;
			y[1] = y[0] + a * TMath::Sqrt(3) / 2.0;
			x[2] = x[1] + a;
			y[2] = y[1];
			x[3] = x[2] + a / 2.0;
			y[3] = y[1] - a * TMath::Sqrt(3) / 2.0;
			;
			x[4] = x[2];
			y[4] = y[3] - a * TMath::Sqrt(3) / 2.0;
			;
			x[5] = x[1];
			y[5] = y[4];

			map->AddBin(6, x, y);

			// Go up
			ytemp += a * TMath::Sqrt(3);
		}

		// Increment the starting position
		if (sCounter % 2 == 0)
			yloop += a * TMath::Sqrt(3) / 2.0;
		else
			yloop -= a * TMath::Sqrt(3) / 2.0;
		xloop += 1.5 * a;
	}
}

void HGCSSGeometryConversion::fillXY(TH2Poly* hist, std::map<int, std::pair<double, double> > &geom) {

	TIter next(hist->GetBins());
	TObject *obj = 0;
	TH2PolyBin *polyBin = 0;
	geom.clear();

	while ((obj = next())) {

		polyBin = (TH2PolyBin*) obj;
		int id = polyBin->GetBinNumber();
		std::pair<double, double> xy = std::pair<double, double>(
				(polyBin->GetXMax() + polyBin->GetXMin()) / 2.,
				(polyBin->GetYMax() + polyBin->GetYMin()) / 2.);
		geom.insert(std::pair<unsigned, std::pair<double, double> >(id, xy));
	}

	std::cout << " -- Check geomMap: size = " << geom.size() << std::endl;
	
}

void HGCSSGeometryConversion::deleteHistos(std::vector<TH2Poly*> &aVec) {

	if (aVec.size() != 0) {
		for (unsigned iL(0); iL < aVec.size(); ++iL) {
			//if (aVec[iL]) aVec[iL]->Delete();
			delete aVec[iL];
			aVec[iL] = 0;
		}
		aVec.clear();
	}
}

void HGCSSGeometryConversion::setGranularity(const std::vector<unsigned> &granul) {

	granularity_.reserve(granul.size());
	for (unsigned iL(0); iL < granul.size(); ++iL) {
		granularity_.push_back(granul[iL]);
	}
}

double HGCSSGeometryConversion::cellSize(const unsigned aLayer,	const double aR) const {

	if (theDetector().subDetectorByLayer(aLayer).isScint)
		return 10. * granularity_[aLayer];
	return cellSize_;
}

double HGCSSGeometryConversion::cellSizeInCm(const unsigned aLayer, const double aR) const {

	return cellSize(aLayer, aR) / 10.;
}

void HGCSSGeometryConversion::initialiseHistos(const bool recreate,
		std::string uniqStr, const bool print) {

	for (unsigned iS(0); iS < theDetector().nSections(); ++iS) {
		resetVector(HistMapE_[theDetector().detType(iS)], "EmipHits" + uniqStr,
				theDetector().detName(iS),
				theDetector().subDetectorBySection(iS),
				theDetector().nLayers(iS), recreate, print);

		std::vector<double> avgvecE;
		avgvecE.resize(theDetector().nLayers(iS), 0);
		avgMapE_[theDetector().detType(iS)] = avgvecE;

		resetVector(HistMapN_[theDetector().detType(iS)], "nHits" + uniqStr,
				theDetector().detName(iS),
				theDetector().subDetectorBySection(iS),
				theDetector().nLayers(iS), recreate, print);

		std::vector<double> avgvecN;
		avgvecN.resize(theDetector().nLayers(iS), 0);
		avgMapN_[theDetector().detType(iS)] = avgvecN;
	}
}

void HGCSSGeometryConversion::fill(std::vector<std::vector<int> > &filledBins, const DetectorEnum type, const unsigned newlayer, const double &weightedE, const double &aTime, const double &posx, const double &posy, const double &posz) {

        int currentBin = HistMapN_[type][newlayer]->FindBin(posx, posy);
        if (!(std::find(filledBins.at(newlayer).begin(), filledBins.at(newlayer).end(), currentBin) != filledBins.at(newlayer).end())) {

                HistMapN_[type][newlayer]->Fill(posx, posy, 1);
                avgMapN_[type][newlayer] += 1;
        }

        filledBins.at(newlayer).push_back(currentBin);

	HistMapE_[type][newlayer]->Fill(posx, posy, weightedE);
	avgMapE_[type][newlayer] += weightedE;

}

TH2Poly* HGCSSGeometryConversion::get2DHist(const unsigned layer, std::string name) {

	const HGCSSSubDetector &subdet = theDetector().subDetectorByLayer(layer);
	unsigned newlayer = layer - subdet.layerIdMin;
	if (name == "E")
		return HistMapE_[subdet.type][newlayer];
        else if (name == "N")
		return HistMapN_[subdet.type][newlayer];
	else {
		std::cerr << " ERROR !! Unknown histogram name. Exiting..."
				<< std::endl;
		exit(1);
	}
}

unsigned HGCSSGeometryConversion::getGranularity(const unsigned aLayer,
		const HGCSSSubDetector &adet) {

	unsigned idx = adet.layerIdMin + aLayer;
	return granularity_[idx];
}

void HGCSSGeometryConversion::resetVector(std::vector<TH2Poly*> &aVec,
		std::string aVar, std::string aString, const HGCSSSubDetector &aDet,
		const unsigned nLayers, bool recreate, bool print) {

	if (recreate) {

		for (unsigned iL(0); iL < aVec.size(); ++iL) {
			aVec[iL]->Delete();
		}
		aVec.clear();
	}
	if (aVec.size() != 0) {

		for (unsigned iL(0); iL < aVec.size(); ++iL) {
			aVec[iL]->Reset("");
		}

	} else {

		if (nLayers > 0) {
			aVec.resize(nLayers, 0);

			if (print)
				std::cout << " -- Creating " << nLayers << " 2D histograms for "
						<< aVar << " " << aString << std::endl;

			for (unsigned iL(0); iL < nLayers; ++iL) {

				std::ostringstream lname;
				lname << aVar << "_" << aString << "_" << iL;
				aVec[iL] = new TH2Poly();
				aVec[iL]->SetName(lname.str().c_str());

				if (print && aVar.find("EmipHits") != aVar.npos)
					std::cout << " ---- Layer " << iL;
                                if (aDet.isScint) {
                                        double newcellsize = 10. * getGranularity(iL, aDet);
                                        initialiseSquareMap(aVec[iL], width_, newcellsize, print && aVar.find("EmipHits") != aVar.npos);

				} else { 
                                        initialiseHoneyComb(aVec[iL], width_, cellSize_,
							print && aVar.find("EmipHits") != aVar.npos);
                                }
			}
		}
	}
}

