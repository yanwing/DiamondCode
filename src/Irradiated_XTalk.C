#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sstream>
#include <string>
#include <vector> 
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TPaveLabel.h"
#include "TStyle.h"
#include "TTree.h"
#include "TMath.h"
#include "map"
#include "TROOT.h"
#include "TRint.h"

using namespace std;

void Irradiated_XTalk(int structure)

{

 gROOT->ForceStyle();

 TH1::SetDefaultSumw2();
 TH2::SetDefaultSumw2();
 TProfile::SetDefaultSumw2();
 TProfile2D::SetDefaultSumw2();
 
 //////////////////////////////////////////////////
 // Colouring                                    //
 //////////////////////////////////////////////////
 
 const Int_t NRGBs = 5;
 const Int_t NCont = 255;

 Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
 Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
 Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
 Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
 
 TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

 gStyle->SetNumberContours(NCont);
 gStyle->SetOptStat(0);
 
 //////////////////////////////////////////////////
 // Definition of structures                     //
 //////////////////////////////////////////////////
 
 std::map<int,std::string> struc_ID;
 struc_ID[1] = "Basic, 20 #mum";
 struc_ID[2] = "Rectangular B, 10 #mum";
 struc_ID[3] = "Rectangular A, 10 #mum";
 struc_ID[4] = "Varying, 10 #mum";
 struc_ID[5] = "Equalize, 10 #mum";
 struc_ID[6] = "Basic, 10 #mum";
 struc_ID[7] = "Varying, 20 #mum";
 struc_ID[8] = "Equalize, 20 #mum";
 struc_ID[9] = "Basic, 20 #mum";
 struc_ID[10] = "Rectangular B, 20 #mum";
 struc_ID[11] = "Rectangular A, 20 #mum";
 struc_ID[12] = "Basic, 10 #mum";
 
 
 //////////////////////////////////////////////////
 // Hit-Strip Identification                     //
 //////////////////////////////////////////////////
 
 int hit_strip_3[56]  = {7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62};
 
 int hit_strip_8[59] = {66,67,68,69,70,71,72,73,74,75,76,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123};
   
 int hit_strip_10[61] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61};
 
 int hit_strip_12[61] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62};

 int histo_number = 0;
 
 // set number of histograms
 if (structure == 3)
    histo_number = 56;
 else if (structure == 8)
    histo_number = 59;
 else if (structure == 10 || structure == 12)
    histo_number = 61;
 
 // Create File
 TFile* File = new TFile(Form("../data/XTalk_Irradiated_Structure_%d.root",structure),"RECREATE");
 
 // open data file
 TFile* DataFile = 0;
 DataFile = TFile::Open(Form("../data/Irradiated_Structure_Noise_Cut_%d_off.root",structure));
 
 // open noise file
 TFile* NoiseFile = 0;
 NoiseFile = TFile::Open("../data/Noise.root");
 
 int upper_strip = 0;
 int lower_strip = 0;
 int lower_diff = 0;
 int upper_diff = 0;
 
 // set range of channels with beam
 if(structure ==  3 || structure == 10 || structure == 12)
 {
 	lower_strip = 1;
 	upper_strip = 62;
 	lower_diff  = -127;
 	upper_diff  = 62;
 }
 else if(structure == 8)
 {
 	lower_strip = 62;
 	upper_strip = 127;
 	lower_diff  = -62;
 	upper_diff  = 127;
 }

 double Signal = 0;
 double SignalErr = 0;
 double sigma = 0;
 
 // create pointers for histograms
 TH1D* CrossTalkHisto;
 TH1F* DeConvHisto;
 TH1D* DeConvCrossTalk;
 TH2D* XTalkMap;
 TH2D* NoiseErrorMap;
 TH2F* dummy_histo = 0;
 
 // get the noise-histogram of the corresponding structure
 NoiseFile->cd();
 TH2F* NoiseMap = (TH2F*)gDirectory->Get(Form("../data/Noise_Histo_%d", structure));
 
 File->cd();
 
 //////////////////////////////////////////////////
 // Create Histograms                            //
 //////////////////////////////////////////////////
 
 CrossTalkHisto = new TH1D(Form("cross_talk_structure_%d",structure),Form("Crosstalk, embedded structure: %s; channel difference; arb. units",struc_ID[structure].c_str()),189, lower_diff-0.5, upper_diff-0.5);
  
 DeConvHisto = new TH1F(Form("deconvhisto_structure_%d",structure),"Deconvolution Histogram; channel difference; counts",189, lower_diff-0.5, upper_diff-0.5); 
 
 DeConvCrossTalk = new TH1D(Form("deconv_cross_talk_struc_%d",structure),Form("Crosstalk (normalised), embedded structure: %s; channel difference; arb. units",struc_ID[structure].c_str()),189, lower_diff-0.5, upper_diff-0.5);
 
 XTalkMap = new TH2D(Form("xtalkmap_struc_%d",structure),Form("Crosstalk Map, embedded structure: %s; channel difference; channel number of strip with beam; arb. units",struc_ID[structure].c_str()),256, -128.5, 127.5,histo_number,lower_strip-0.5,upper_strip-0.5);
 
 int c = 0;
 
 // create histogram for normalisation
 for (int a = 0; a < 128; a++)
 {
 	for (int b = 0; b < histo_number; b++)
  	{
   		if(structure == 8)   
   		c = -61 + a + b;
             	else if (structure == 3 || structure == 10|| structure == 12) 
             	c = -126 + a + b;
   		   		
   		DeConvHisto->Fill(c);
   		int d = DeConvHisto->GetXaxis()->FindBin(c);
   		DeConvHisto->SetBinError(d,0.0);
   		
  	}
 }

  int hitstrip = 0;
 
 double SignalRatio = 0;
 double NoiseContent = 0;
 double NoiseErr = 0;
 double threshold_step = 0;
 
 int channeldiff = 0;
 double xErr = 0;
 
 for (int i = 0; i < histo_number; i++) // data loop
 {	
 
 	threshold_step = 48.8; // first threshold
 
 	for (int j = 0; j < 16; j++) // threshold loop
 	{
 	
 		for (int k = 0; k <  128; k++) // channel loop
 		{
 		 
 		 // get channel with beam
 		 if(structure == 3)
 		 hitstrip = hit_strip_3[i];
 		 else if(structure == 8)
 		 hitstrip = hit_strip_8[i];
 		 else if(structure == 10)
 		 hitstrip = hit_strip_10[i];
 		 else if(structure == 12)
 		 hitstrip = hit_strip_12[i];

 		 DataFile->cd();
 		 
 		 // get hitmap
 		 dummy_histo = (TH2F*)gDirectory->Get(Form("Hitmap_Noise_Cut_Structure_%d_Channel_%d",structure,hitstrip));
 		 
 		 File->cd(); 
 	 	 
 	 	 // get signal and background data
                 Signal    = dummy_histo->GetBinContent((k+1),(j+1));                 
                 SignalErr = dummy_histo->GetBinError((k+1),(j+1));
                 
                 NoiseContent = NoiseMap->GetBinContent((k+1),(j+1));
                 NoiseErr = NoiseMap->GetBinError((k+1),(j+1));
                
                 
                 if (k != hitstrip)
                 {
                 	  // apply cut: signal has to be larger than noise
		          if (NoiseContent > 0  && Signal > NoiseContent && (Signal-NoiseContent)>(SignalErr + NoiseErr)) 
		          {
		          	   // calculate error
				   sigma = sqrt(((Signal*NoiseErr)/(NoiseContent*NoiseContent))*((Signal*NoiseErr)/(NoiseContent*NoiseContent)) + (SignalErr/NoiseContent)*(SignalErr/NoiseContent));
				   
				   // calculate siganl ratio
				   SignalRatio = Signal/NoiseContent; 
				  
				   // calculate distance between signal-channel and channel where signal > noise
				   channeldiff = hitstrip - k;
				    
				   // fill histograms
				   CrossTalkHisto->Fill(channeldiff,SignalRatio);				   
				   XTalkMap->Fill(channeldiff,hitstrip,SignalRatio);
						  
		          }
                 }      
                            
 		 } // closing channel loop
 		
 		threshold_step += 3.2; // increase threshold
 		
 	} // closing threshold loop
 	
 	threshold_step = 48.8; // set threshold back to intial value
 	
 } // closing file loop
 
 // normalise crosstalk-histogram
 DeConvCrossTalk->Divide(CrossTalkHisto,DeConvHisto,1.0,1.0);
 
 File->cd();
 
 // write file
 File->Write();

}
 		 
