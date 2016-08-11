#include <string>
#include <set>
#include <iostream>
#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "Math/Vector4D.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

template<class T>

void extractParameterFromStr(std::string aStr, T &vec) {

	if (aStr == "")
		return;

	std::vector<std::string> layVec;
	boost::split(layVec, aStr, boost::is_any_of(","));

	for (unsigned iE(0); iE < layVec.size(); ++iE) { //loop on elements

		std::vector<std::string> lPair;
		boost::split(lPair, layVec[iE], boost::is_any_of(":"));

		if (lPair.size() != 2) {

			std::cout << " -- Wrong string for parameter given as input:"
					<< layVec[iE]
					<< " Try again, expecting exactly one symbol \":\" between two \",\" ..."
					<< std::endl;
			exit(1);
		}

		std::vector<std::string> lLay;
		boost::split(lLay, lPair[0], boost::is_any_of("-"));

		if (lLay.size() > 2) {

			std::cout << " -- Wrong string for granularities given as input:"
					<< lPair[0]
					<< " Try again, expecting at most one symbol \"-\"."
					<< std::endl;
			exit(1);
		}

		unsigned beginIdx = atoi(lLay[0].c_str());
		unsigned endIdx   = lLay.size() == 1 ? beginIdx : atoi(lLay[1].c_str());

		for (unsigned iL(beginIdx); iL < endIdx + 1; ++iL) {

			if (iL < vec.size())
				std::istringstream(lPair[1]) >> vec[iL];
			else {
				std::cout << " -- WARNING! Input parameter has more layers: "
						<< endIdx << " than detector : " << vec.size()
						<< ". Ignoring additional layer #" << iL
						<< "... PLEASE CHECK SETTINGS ARE CORRECT FOR EXISTING LAYERS!!"
						<< std::endl;
			}
		}
	} //loop on elements
}

void processHist(unsigned int layer, TH2Poly* histE, TH2Poly* histN, std::vector<TH2Poly*> histNTot, HGCSSDigitisation &myDigitiser, const HGCSSSubDetector &subdet, const std::vector<unsigned> &pThreshInADC) {

	DetectorEnum adet = subdet.type;
	bool isSi         = subdet.isSi;

        unsigned int xBins = histE->GetNbinsX();
        unsigned int yBins = histE->GetNbinsY();

        for (unsigned int x = 0; x < xBins; x++) {
                for (unsigned int y = 0; y < yBins; y++) {

                        int globalBin = histE->GetBin(x,y);
                        double digiE = histE->GetBinContent(globalBin);

                        if (digiE == 0) 
                                continue;
                        myDigitiser.addNoise(digiE, layer);
                        //for silicon-based Calo
                        unsigned adc = 0;

                        if (isSi) {

                                adc   = myDigitiser.adcConverter(digiE, adet);
                                digiE = myDigitiser.adcToMIP(adc, adet);
                        }

                        if (!(isSi && adc > pThreshInADC[layer])) {
                            
                                histN->SetBinContent(globalBin, 0);
                        }
                        else {
                                histE->SetBinContent(globalBin, digiE);
                                double temp = histNTot[layer]->GetBinContent(globalBin) + 1;

                                histNTot[layer]->SetBinContent(globalBin,temp);
                        }
                }
        }
}	//processHist

int main(int argc, char** argv) {	 

	/////////////////////////////////////////////////////////////
	//parameters
	/////////////////////////////////////////////////////////////

        int opt;
        int option_index = 0;
        static struct option long_options[] = {
            {"pNevts",      required_argument, 0, 'N'},
            {"inFilePath",  required_argument, 0, 'I'},
            {"outFilePath", required_argument, 0, 'O'},
            {"outFilename", required_argument, 0, 'F'},
            {"granulStr",   required_argument, 0, 'G'},
            {"noiseStr",    required_argument, 0, 'S'},
            {"threshStr",   required_argument, 0, 'T'},
            {"interCalib",  required_argument, 0, 'C'},
            {"nSiLayers",   required_argument, 0, 'L'}
        };

        int pNevts = -1;
        std::string inFilePath = "";
        std::string outFilePath = "";
        std::string outFilename = "";
        std::string granulStr = "";
        std::string noiseStr = "";
        std::string threshStr = "";
        double interCalib = 0;
        int nSiLayers = -1;

        while((opt = getopt_long(argc,argv,"N:I:O:F:G:S:T:C:L:", long_options, &option_index)) != -1) {
            switch(opt) {
                case 'N':
                    pNevts = int(atoi(optarg));
                    break;

                case 'I':
                    inFilePath = optarg;
                    break;

                case 'O':
                    outFilePath = optarg;
                    break;
                    
                case 'F':
                    outFilename = optarg;
                    break;

                case 'G':
                    granulStr = optarg;
                    break;

                case 'S':
                    noiseStr = optarg;
                    break;

                case 'T':
                    threshStr = optarg;
                    break;

                case 'C':
                    interCalib = double(atoi(optarg));
                    break;

                case 'L':
                    nSiLayers = int(atoi(optarg));
                    break;

            }
        }

        std::cout << inFilePath << std::endl;
	std::cout << " ----------------------------------------" << std::endl
			<< " -- Input parameters: " << std::endl << " -- Input file path: "
			<< inFilePath << std::endl << " -- Output file path: "
			<< outFilePath << std::endl << " -- Processing ";

	if (pNevts > 0)	std::cout << pNevts;
        
	else
		std::cout << "all";
	                std::cout << " events." << std::endl << " -- Granularities: " << granulStr
			<< std::endl << " -- noise: " << noiseStr << std::endl
			<< " -- thresholds: " << threshStr << std::endl
			<< " -- intercalibration factor (in %): " << interCalib << std::endl;

	unsigned pSeed = 0;
	/////////////////////////////////////////////////////////////
	//input
	/////////////////////////////////////////////////////////////

	std::string inputStr   = inFilePath;
	TFile       *inputFile = TFile::Open(inputStr.c_str());

	if (!inputFile) {

		std::cout << " -- Error, input file " << inputStr
				<< " cannot be opened. Exiting..." << std::endl;
		return 1;
	}

	TTree *inputTree = (TTree*) inputFile->Get("HGCSSTree");

	if (!inputTree) {

		std::cout << " -- Error, tree HGCSSTree  cannot be opened. Exiting..." << std::endl;
		return 1;
        }

	HGCSSInfo      *info         = (HGCSSInfo*) inputFile->Get("Info");
	const double   calorSizeXY   = info->calorSizeXY();
	const double   cellSize      = info->cellSize();
	const unsigned versionNumber = info->version();
	const unsigned model         = info->model();

	HGCSSEvent               *event  = 0;
	std::vector<HGCSSSimHit> *hitvec = 0;

	inputTree->SetBranchAddress("HGCSSEvent"    , &event);
	inputTree->SetBranchAddress("HGCSSSimHitVec", &hitvec);

	//initialise detector
	HGCSSDetector &myDetector = theDetector();

	std::cout << " -- Calor size XY = " << calorSizeXY << ", version number = "
			<< versionNumber << ", model = " << model << ", cellSize = "
			<< cellSize << std::endl;

        myDetector.buildDetector(versionNumber);
	//initialise calibration class
	HGCSSCalibration mycalib(nSiLayers);

	const unsigned nLayers = myDetector.nLayers();
	HGCSSGeometryConversion geomConv(model, cellSize, nSiLayers);
	geomConv.setXYwidth(calorSizeXY);
        geomConv.initialiseHoneyComb(calorSizeXY, cellSize);

	HGCSSDigitisation myDigitiser;

	std::vector<unsigned> granularity;
	granularity.resize(nLayers, 1);
	std::vector<double> pNoiseInMips;
	pNoiseInMips.resize(nLayers, 0.12);
	std::vector<unsigned> pThreshInADC;
	pThreshInADC.resize(nLayers, 5);

	extractParameterFromStr<std::vector<unsigned> >(granulStr, granularity);
	extractParameterFromStr<std::vector<double> >(noiseStr, pNoiseInMips);
	extractParameterFromStr<std::vector<unsigned> >(threshStr, pThreshInADC);

	std::cout << " -- Granularities and noise are setup like this:"
			<< std::endl;

	for (unsigned iL(0); iL < nLayers; ++iL) {

		std::cout << "Layer ";
		if (iL < 10)
			std::cout << " ";
		std::cout << iL << " : " << granularity[iL] << ", " << pNoiseInMips[iL]
				<< " mips, " << pThreshInADC[iL] << " adc - ";
		if (iL % 5 == 4)
			std::cout << std::endl;
		myDigitiser.setNoise(iL, pNoiseInMips[iL]);
	}

	std::cout << std::endl;

	geomConv.setGranularity(granularity);
	geomConv.initialiseHistos();
	TRandom3 *lRndm = new TRandom3();
	lRndm->SetSeed(pSeed);
	myDigitiser.setRandomSeed(pSeed);

	std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
			<< " ----------------------------------------" << std::endl;

	myDigitiser.setIntercalibrationFactor(interCalib);

	/////////////////////////////////////////////////////////////
	//output
	/////////////////////////////////////////////////////////////

	std::ostringstream outputStr;
        outputStr << outFilePath << "/" << outFilename << ".root";

	TFile *outputFile = TFile::Open(outputStr.str().c_str(), "RECREATE");

	if (!outputFile) {

		std::cout << " -- Error, output file " << outputStr.str()
				<< " cannot be opened. Exiting..." << std::endl;
		return 1;

	} else {

		std::cout << " -- File will be saved as " << outputStr.str()
				<< std::endl;
	}

	HGCSSInfo *lInfo = new HGCSSInfo();
	lInfo->calorSizeXY(calorSizeXY);
	lInfo->cellSize(cellSize);
	lInfo->version(versionNumber);
	lInfo->model(model);

	TTree *outputTree = new TTree("RecoTree", "HGC Standalone simulation reco tree");
	HGCSSRecoHitVec lRecoHits;
        HGCSSEvent      lEvent;
	unsigned maxRecHits = 0;
       	outputTree->Branch("HGCSSRecoHitVec", "std::vector<HGCSSRecoHit>", &lRecoHits);

	/////////////////////////////////////////////////////////////
	//Loop on events
	/////////////////////////////////////////////////////////////

	const unsigned nEvts =
			(pNevts > inputTree->GetEntries() || pNevts == 0) ?
					static_cast<unsigned>(inputTree->GetEntries()) : pNevts;
	std::cout << "- Processing = " << nEvts << " events out of "
			<< inputTree->GetEntries() << std::endl;

        // Initialize vector of histos for hits dist in layer
        std::vector<TH1F*> nHitsHist;
        std::vector<TH1F*> MIPsHist;
        std::vector<TH1F*> MIPsHistLayerSum;
        for (unsigned int firstLayer = 1; firstLayer < 7; firstLayer++) {
            for (unsigned int lastLayer = 6; lastLayer < 17; lastLayer++) {

                  char *histname = new char[20];
                  sprintf(histname, "MIPs_Sum_%d-%d", firstLayer, lastLayer);
                  MIPsHistLayerSum.push_back(new TH1F(histname, histname, 16000, -0.5, 3999.5));
            }
        }

        char *histname1 = new char[20];
        char *histname2 = new char[20];
        for (unsigned iL(0); iL < nLayers; ++iL) {
                sprintf(histname1, "nHits_Layer_%d", iL);
                nHitsHist.push_back(new TH1F(histname1, histname1, 1000, -0.5, 999.5));
                sprintf(histname2, "MIPs_Layer_%d", iL);
                MIPsHist.push_back(new TH1F(histname2, histname2, 16000, -0.5, 3999.5));

        }

        std::vector<int> eventSum(nLayers, 0);
        // Loop over all events in file
        std::vector<TH2Poly*> histNTot(100);
        for (unsigned ievt(0); ievt < nEvts; ++ievt) {

                inputTree->GetEntry(ievt);
		lEvent.vtx_x(event->vtx_x());
		lEvent.vtx_y(event->vtx_y());
		lEvent.vtx_z(event->vtx_z());
		mycalib.setVertex(lEvent.vtx_x(), lEvent.vtx_y(), lEvent.vtx_z());

	        if (ievt % 5 == 0)
			std::cout << "... Processing event: " << ievt << std::endl;

		unsigned prevLayer   = 10000;
		DetectorEnum type    = DetectorEnum::ECAL;
                unsigned subdetLayer = 0;
		bool     isScint     = false;
                std::vector<std::vector<int> > bins(100);  
                // Loop on simhits in a given event.
		for (unsigned iH(0); iH < (*hitvec).size(); ++iH) {

                        // Create SimHit object for iHth sim hit and acquire the simhits layer
			HGCSSSimHit lHit = (*hitvec)[iH];
                        unsigned layer = lHit.layer();
                        const HGCSSSubDetector &subdet = myDetector.subDetectorByLayer(layer);

			if (layer != prevLayer) {

				const HGCSSSubDetector &subdet = myDetector.subDetectorByLayer(layer);
				isScint = subdet.isScint;
				type = subdet.type;
                                subdetLayer = layer - subdet.layerIdMin;
				prevLayer = layer;
			}
                        
			std::pair<double, double> xy = lHit.get_xy(isScint, geomConv);
			double posx     = xy.first;  
			double posy     = xy.second;
			double radius   = sqrt(pow(posx, 2) + pow(posy, 2));
			double posz     = lHit.get_z();
                        double realtime = 0;
                        double energy   = lHit.energy() * mycalib.MeVToMip(layer, radius);

			if (energy > 0 && lHit.silayer() < 3) {
                                
		                geomConv.fill(bins, type, subdetLayer, energy, realtime, posx, posy, posz);
                        }
		}
                std::vector<double> layerSums(100,0.0);
                for (unsigned iL(0); iL < nLayers; ++iL) {

                        const HGCSSSubDetector &subdet = myDetector.subDetectorByLayer(iL);
                        TH2Poly *histN = geomConv.get2DHist(iL, "N");
                        TH2Poly *histE = geomConv.get2DHist(iL, "E");

                        if (ievt == 0) {
                                TH2Poly *histNClone = (TH2Poly*)histN->Clone();
                                histNTot[iL] = histNClone;
                        }

                        processHist(iL, histE, histN, histNTot, myDigitiser, subdet, pThreshInADC);
             
                        int Nintegral = histN->Integral();
                        double Eintegral = histE->Integral();

                        layerSums[iL] = Eintegral;
                        nHitsHist[iL]->Fill(Nintegral);
                        MIPsHist[iL]->Fill(Eintegral);

                        eventSum[iL] = Nintegral;

                }

                int counter = 0;
                for (int firstLayer = 1; firstLayer < 7; firstLayer++) {
                    for (int lastLayer = 6; lastLayer < 17; lastLayer++) {

                        double layerSum = 0;
                        for (int layer = firstLayer; layer <= lastLayer; layer++) {
                            layerSum += layerSums[layer];

                        }

                        MIPsHistLayerSum[counter]->Fill(layerSum);
                        counter += 1;

                    }
                }

                bins.clear();
                outputTree->Fill();
		//reserve necessary space and clear vectors.
		if (lRecoHits.size() > maxRecHits) {

			maxRecHits = 2 * lRecoHits.size();
			std::cout << " -- INFO: event " << ievt << " maxRecHits updated to "
					<< maxRecHits << std::endl;
		} 

		lRecoHits.clear();
		geomConv.initialiseHistos();
		lRecoHits.reserve(maxRecHits);

        }

        outputFile->cd();
        // Loop on every layer and acquire the filled TH2Poly histo summed over all events.
        for (unsigned iL(0); iL < nLayers; ++iL) {

                TH2Poly* histE = geomConv.get2DHist(iL, "E");
                //histNTot[iL]->Write();
                //histE->Write();
        }

	outputFile->cd();

        for (unsigned int i = 0; i < nHitsHist.size(); i++) {

                nHitsHist[i]->Write();
                MIPsHist[i]->Write();
        }

        for (unsigned int j = 0; j < MIPsHistLayerSum.size(); j++) {

                MIPsHistLayerSum[j]->Write();
        }

	outputFile->WriteObjectAny(lInfo, "HGCSSInfo", "Info");
	outputTree->Write();
	outputFile->Close();

	return 0;
}	
