#include <string>
#include <math.h>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include "TFile.h"
#include "TCutG.h"
#include "TRint.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TColor.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "Math/Vector4D.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSRecoJet.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSSamplingSection.hh"

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

int main(int argc, char** argv) {

        if (argc != 3) {
                std::cout << "Incorrect number (" << argc << ")" << " of arguments provided, 2 expected!" << std::endl;
        }
        std::string inputf1 = argv[1];
        int linecolor       = int(atoi(argv[2]));

        // Open CSV file ofstream and write top row
        /*ofstream csvFile1;
        csvFile1.open("4GeV_1M_AllLayers_NumHits.csv");
        for (int numHits = 0; numHits < 500; numHits++) {
                 
                 std::string header;

                 if (numHits == 0) {

                         csvFile1 << " ,NumHits 0,";
                 }
                 else if (numHits == 499) {
                         header = (boost::format("NumHits %s\n") %numHits).str();
                 }
                 else {
                     
                         header = (boost::format("NumHits %s,") %numHits).str();
                 }
                 csvFile1 << header;
        } 

        ofstream csvFile2;
        csvFile2.open("4GeV_1M_AllLayers_Energy.csv");
        for (int energy = 0; energy < 100; energy++) {
                 
                 std::string header;

                 if (energy == 0) {

                         csvFile2 << " ,energy bin 0,";
                 }
                 else if (energy == 99) {
                         header = (boost::format("energy bin %s\n") %energy).str();
                 }
                 else {
                     
                         header = (boost::format("energy bin %s,") %energy).str();
                 }
                 csvFile2 << header;
        } */

	TFile *infile1 = TFile::Open(inputf1.c_str());
	TTree *tree1 = (TTree*) infile1->Get("RecoTree");

	std::vector<HGCSSRecoHit> *hitVec1 = 0;
	tree1->SetBranchAddress("HGCSSRecoHitVec", &hitVec1);

        TFile *outfile = TFile::Open("analyzed_digi.root","RECREATE");
/*
        TCutG *polycut = new TCutG("bigHexagon", 6);

        for (int i = 0; i < 7; i++) {
              int x = 30*std::cos(2*M_PI*i/6);
              int y = 30*std::sin(2*M_PI*i/6);
              polycut->SetPoint(i, x, y);
        }*/
        for (int layer = 0; layer < 41; layer++) {

                // Write out CSV file row header with current layer
                /*std::string header = (boost::format("Layer %s,") %layer).str();
                csvFile1 << header;
                csvFile2 << header;*/
                
                /////////////////////////// String, Histo definitions ///////////////////////////

                std::string hLT = (boost::format("Hits in Layer %s") %layer).str();
                const char *hitsLayerTitle = hLT.c_str();

                //std::string nHT = (boost::format("Hits in Layer %s") %layer).str();
                //const char *nHitsTitle = nHT.c_str();
                
                //std::string nHHN = (boost::format("nHits_ECAL_%s") %layer).str();
                //const char *nHitsHistName = nHHN.c_str();

                std::string MHN = (boost::format("MIPs_Layer_%s") %layer).str();
                const char *MIPsHistName = MHN.c_str();
                
                std::string MHNW = (boost::format("MIPs_Layer_%s_Wide") %layer).str();
                const char *MIPsHistNameWide = MHNW.c_str();

                std::string MT = (boost::format("MIPs Deposit in Layer %s") %layer).str();
                const char *MIPsTitle = MT.c_str();

                //std::string nHF = (boost::format("nHits_ECAL_%s.png") %layer).str();
                //const char *nHitsFilename = nHF.c_str();

                //std::string MFW = (boost::format("MIPs_Layer_%s_wide.png") %layer).str();
                //const char *MIPsFilenameWide = MFW.c_str();

                //std::string MF = (boost::format("MIPs_Layer_%s.png") %layer).str();
                //const char *MIPsFilename = MF.c_str();

                //std::string nHLF = (boost::format("nHits_Layer_%s.png") %layer).str();
                //const char *nHitsLayerFilename = nHLF.c_str();

                std::string nHLHN = (boost::format("nHits_Layer_%s") %layer).str();
                const char *nHitsLayerHistName = nHLHN.c_str();  

                //TH2Poly *hits1   = (TH2Poly*)infile1->Get(nHitsHistName);
                //TH2Poly *MIPs1 = (TH2Poly*)infile1->Get(MIPsHistName);
                TH1F    *nHits1  = (TH1F*)infile1->Get(nHitsLayerHistName);
                TH1F    *MIPs1 = (TH1F*)infile1->Get(MIPsHistName);
                //TH1F    *MIPs2 = (TH1F*)MIPs1->Clone();

                ////////////////////////////////////////////////////////////////////////////////
               
                /////////////////////////// Fill CSV1 File /////////////////////////////////
                
                /*for (int bin = 0; bin < nHits1->GetSize()-2; bin++) {

                         std::string binContent;

                         if (bin == nHits1->GetSize()-3) {
                                  binContent = (boost::format("%s\n") %nHits1->GetBinContent(bin)).str();
                         }
                         else {
                                  binContent = (boost::format("%s,") %nHits1->GetBinContent(bin)).str();
                         }
                         
                         csvFile1 << binContent;
                }
                /////////////////////////// Fill CSV2 File /////////////////////////////////
                
                for (int bin = 0; bin < MIPs1->GetSize()-2; bin++) {

                         std::string binContent;

                         if (bin == MIPs1->GetSize()-3) {
                                  binContent = (boost::format("%s\n") %MIPs1->GetBinContent(bin)).str();
                         }
                         else {
                                  binContent = (boost::format("%s,") %MIPs1->GetBinContent(bin)).str();
                         }
                         
                         csvFile2 << binContent;
                }*/

                        
                /////////////////////////// Color Palette //////////////////////////////

                /*const Int_t NRGBs = 5;
                const Int_t NCont = 255;

                Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
                Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
                Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
                Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
                TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);*/

                ////////////////////////////////////////////////////////////////////////

                /////////////////////////// nHits TH2Poly //////////////////////////////
                
                /*TCanvas *myCanvas1 = new TCanvas(nHitsTitle, nHitsTitle, 1600, 1600);
                
                myCanvas1->SetLogz();

                hits1->SetTitle(nHitsTitle);
                hits1->SetStats(kFALSE);

                hits1->GetXaxis()->SetTitleSize(0.03);
                hits1->GetXaxis()->SetLabelSize(0.02);
                hits1->GetXaxis()->SetTitle("Position (cm)");
                hits1->GetXaxis()->SetRangeUser(-150,150);

                hits1->GetYaxis()->SetTitleSize(0.03);
                hits1->GetYaxis()->SetLabelSize(0.02);
                hits1->GetYaxis()->SetTitle("Position (cm)");
                hits1->GetYaxis()->SetRangeUser(-150,150);

                hits1->GetZaxis()->SetRangeUser(0,25000);
                hits1->GetZaxis()->SetLabelSize(0.02);

                hits1->SetMarkerSize(0.7);
                
                hits1->SetContour(NCont);

                outfile->cd();
                myCanvas1->Draw("colz[polycut]" "text");
          //      polycut->SetLineWidth(5);
            //    polycut->SetLineColor(1);
                hits1->Draw("colz text");
            //    polycut->Draw("same");
                myCanvas1->Update();
                myCanvas1->SaveAs(nHitsFilename);
                hits1->Write();*/

                //////////////////////////////////////////////////////////////////////

                ////////////////////////// Energy TH2Poly ////////////////////////////

                /*TCanvas *myCanvas2 = new TCanvas(energyTitle, energyTitle, 1200, 1200);

                myCanvas2->SetLogz();

                MIPs1->SetTitle(energyTitle);
                MIPs1->SetStats(kFALSE);

                MIPs1->GetXaxis()->SetRangeUser(-150,150);
                MIPs1->GetXaxis()->SetTitle("Position (cm)");
                MIPs1->GetXaxis()->SetLabelSize(0.02);
                MIPs1->GetXaxis()->SetTitleSize(0.03);

                MIPs1->GetYaxis()->SetRangeUser(-150,150);
                MIPs1->GetYaxis()->SetTitle("Position (cm)");
                MIPs1->GetYaxis()->SetLabelSize(0.02);
                MIPs1->GetYaxis()->SetTitleSize(0.03);

                MIPs1->GetZaxis()->SetRangeUser(0.0001,15000);
                MIPs1->GetZaxis()->SetLabelSize(0.02);

                MIPs1->SetMarkerSize(0.7);

                MIPs1->SetContour(NCont);
                outfile->cd();
                myCanvas2->Draw("colz[polycut]" "text");
              //  polycut->SetLineWidth(5);
              //  polycut->SetLineColor(1);
                MIPs1->Draw("colz text");
             //   polycut->Draw("same");
                myCanvas2->Update();
                myCanvas2->SaveAs(energyFilename);
                MIPs1->Write();*/

                //////////////////////////////////////////////////////////////////////

                //////////////////////////// nHits TH1F //////////////////////////////
                
                TCanvas *myCanvas3 = new TCanvas(nHitsLayerHistName, nHitsLayerHistName, 1200, 1200);
               
                myCanvas3->SetLogy();

                nHits1->SetTitle(hitsLayerTitle);
                nHits1->SetStats(kFALSE);
                nHits1->SetLineColor(6);

                nHits1->GetXaxis()->SetRangeUser(-0.5,29.5);
                nHits1->GetXaxis()->SetTitle("Number of Hits");
                nHits1->GetXaxis()->SetLabelSize(0.02);
                nHits1->GetXaxis()->SetTitleSize(0.03);

                nHits1->GetYaxis()->SetTitle("Number of Events");
                nHits1->GetYaxis()->SetLabelSize(0.02);
                nHits1->GetYaxis()->SetTitleSize(0.03);
                nHits1->SetLineWidth(3);
                nHits1->Draw();

                outfile->cd();
                //myCanvas3->SaveAs(nHitsLayerFilename);
                nHits1->Write();
                
                /////////////////////////////////////////////////////////////////////

                //////////////////////////// MIPs TH1F wide //////////////////////////////
                
                TCanvas *myCanvas5 = new TCanvas(MIPsHistName, MIPsHistName, 1200, 1200);
                
                myCanvas5->cd();
                myCanvas5->SetLogy();

                MIPs1->SetTitle(MIPsTitle);
                MIPs1->SetStats(kFALSE);
                MIPs1->SetLineColor(linecolor);

                MIPs1->GetXaxis()->SetRangeUser(-0.5,299.5);
                MIPs1->GetXaxis()->SetTitle("MIPs");
                MIPs1->GetXaxis()->SetLabelSize(0.02);
                MIPs1->GetXaxis()->SetTitleSize(0.03);

                MIPs1->GetYaxis()->SetTitle("Number of Events");
                MIPs1->GetYaxis()->SetLabelSize(0.02);
                MIPs1->GetYaxis()->SetTitleSize(0.03);
                MIPs1->SetLineWidth(3);
                MIPs1->Draw();

                outfile->cd();
                //myCanvas5->SaveAs(MIPsFilename);
                MIPs1->Write();

                
                //////////////////////////// MIPs TH1F wide //////////////////////////////
                //
                TCanvas *myCanvas4 = new TCanvas(MIPsHistNameWide, MIPsHistNameWide, 1200, 1200);
               
                myCanvas4->cd();
                myCanvas4->SetLogy();

                TH1 *MIPs2 = MIPs1->Rebin(8, MIPsHistNameWide);
                MIPs2->SetTitle(MIPsTitle);
                MIPs2->SetStats(kFALSE);
                MIPs2->SetLineColor(linecolor);

                MIPs2->GetXaxis()->SetRangeUser(-0.5,999.5);
                //MIPs2->Rebin(500);
                MIPs2->GetXaxis()->SetTitle("MIPs");
                MIPs2->GetXaxis()->SetLabelSize(0.02);
                MIPs2->GetXaxis()->SetTitleSize(0.03);

                MIPs2->GetYaxis()->SetTitle("Number of Events");
                MIPs2->GetYaxis()->SetLabelSize(0.02);
                MIPs2->GetYaxis()->SetTitleSize(0.03);
                MIPs2->SetLineWidth(3);
                MIPs2->Draw();

                outfile->cd();
                //myCanvas4->SaveAs(MIPsFilenameWide);
                MIPs2->Write();
                      
        }
        ///////////////////// Make Layer Sum Histos ////////////////////////////////

        for (unsigned int firstLayer = 1; firstLayer < 7; firstLayer++) {
            for (unsigned int lastLayer = 6; lastLayer < 17; lastLayer++) {
                 
                //char *oldhistnameWide      = new char[35];
                char *oldhistname          = new char[35];

                char *histnameWide      = new char[35];
                char *histname          = new char[35];


                char *histnameFile      = new char[35];
                char *histnameFileWide  = new char[35];

                char *histnameTitleWide = new char[35];
                char *histnameTitle     = new char[35];

                sprintf(histnameWide, "MIPs_Sum_Layers_%d-%d_Wide", firstLayer, lastLayer);
                sprintf(histname, "MIPs_Sum_Layers_%d-%d", firstLayer, lastLayer);
                sprintf(oldhistname, "MIPs_Sum_%d-%d", firstLayer, lastLayer);

                sprintf(histnameFile, "MIPs_Sum_Layers_%d-%d.png", firstLayer, lastLayer);
                sprintf(histnameFileWide, "MIPs_Sum_Layers_%d-%d_Wide.png", firstLayer, lastLayer);

                sprintf(histnameTitleWide, "MIPs Sum Layers %d-%d Wide", firstLayer, lastLayer);
                sprintf(histnameTitle, "MIPs Sum Layers %d-%d", firstLayer, lastLayer);

                TCanvas *myCanvas6 = new TCanvas(histname, histname, 1200, 1200);
                TCanvas *myCanvas7 = new TCanvas(histnameWide, histnameWide, 1200, 1200);

                //TH1F *MIPsLayerSumWide = (TH1F*)infile1->Get(histname);
                TH1F *oldMIPsLayerSum = (TH1F*)infile1->Get(oldhistname);
                TH1F *MIPsLayerSum     = (TH1F*)oldMIPsLayerSum->Clone(histname);
                TH1 *MIPsLayerSumWide = MIPsLayerSum->Rebin(8, histnameWide);

                myCanvas6->cd();
                myCanvas6->SetLogy();
                MIPsLayerSum->SetTitle(histnameTitle);
                //MIPsLayerSum->SetStats(kFALSE);
                MIPsLayerSum->SetLineColor(linecolor);
                
                MIPsLayerSum->GetXaxis()->SetRangeUser(-0.5,299.5);
                MIPsLayerSum->GetXaxis()->SetTitle("MIPs");
                MIPsLayerSum->GetXaxis()->SetLabelSize(0.02);
                MIPsLayerSum->GetXaxis()->SetTitleSize(0.03);
                
                MIPsLayerSum->GetYaxis()->SetTitle("Number of Events");
                MIPsLayerSum->GetYaxis()->SetLabelSize(0.02);
                MIPsLayerSum->GetYaxis()->SetTitleSize(0.03);
                MIPsLayerSum->SetLineWidth(3);
                MIPsLayerSum->Draw();
                 
                outfile->cd();
                //myCanvas6->SaveAs(histnameFile);
                MIPsLayerSum->Write();

                myCanvas7->cd();
                myCanvas7->SetLogy();
                MIPsLayerSumWide->SetTitle(histnameTitleWide);
                //MIPsLayerSumWide->SetStats(kFALSE);
                MIPsLayerSumWide->SetLineColor(linecolor);
                
                MIPsLayerSumWide->GetXaxis()->SetRangeUser(-0.5,999.5);
                MIPsLayerSumWide->GetXaxis()->SetTitle("MIPs");
                MIPsLayerSumWide->GetXaxis()->SetLabelSize(0.02);
                MIPsLayerSumWide->GetXaxis()->SetTitleSize(0.03);
                 
                MIPsLayerSumWide->GetYaxis()->SetTitle("Number of Events");
                MIPsLayerSumWide->GetYaxis()->SetLabelSize(0.02);
                MIPsLayerSumWide->GetYaxis()->SetTitleSize(0.03);
                MIPsLayerSumWide->SetLineWidth(3);
                MIPsLayerSumWide->Draw();
                 
                outfile->cd();
                //myCanvas7->SaveAs(histnameFileWide);
                MIPsLayerSumWide->Write();

            }
        }

        //csvFile1.close();
        //csvFile2.close();
        outfile->Close();
        return 0;

}
