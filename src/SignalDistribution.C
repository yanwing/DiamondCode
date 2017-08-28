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

void SignalDistribution(int structure)

{

 gROOT->ForceStyle();

 //TH1::SetDefaultSumw2();
 //TH2::SetDefaultSumw2();
 //TProfile::SetDefaultSumw2();
 //TProfile2D::SetDefaultSumw2();

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
 //gStyle->SetOptStat(0);
 
 //////////////////////////////////////////////////
 // Definition of structures                     //
 //////////////////////////////////////////////////
 
 std::map<int,std::string> struc_ID;
 struc_ID[1]  = "Basic, 20 #mum";
 struc_ID[2]  = "Rectangular B, 10 #mum";
 struc_ID[3]  = "Rectangular A, 10 #mum";
 struc_ID[4]  = "Varying, 10 #mum";
 struc_ID[5]  = "Equalize, 10 #mum";
 struc_ID[6]  = "Basic, 10 #mum";
 struc_ID[7]  = "Varying, 20 #mum";
 struc_ID[8]  = "Equalize, 20 #mum";
 struc_ID[9]  = "Basic, 20 #mum";
 struc_ID[10] = "Rectangular B, 20 #mum";
 struc_ID[11] = "Rectangular A, 20 #mum";
 struc_ID[12] = "Basic, 10 #mum";
 
 //////////////////////////////////////////////////
 // Strips hit by beam                           //
 //////////////////////////////////////////////////
 
 int hitstrip_upper[16] = {64,68,72,76,80, 84,88,92,96,100,104,108,112,116,120,124};
 int hitstrip_low[16] = {1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61};
 int hitstrip = 0;
 int upper_strip = 0;
 int lower_strip = 0;
 int upper_diff = 0;
 int lower_diff = 0;
 
 if      (structure == 9 ||
          structure == 8 ||
          structure == 7 ||
          structure == 6 ||
          structure == 5 ||
          structure == 4) {lower_strip = 64;
          		   upper_strip = 124;
          		   lower_diff = -68;
          		   upper_diff = 127;}
             		      
 else if (structure == 1 ||
          structure == 2 ||
          structure == 3 ||
          structure == 10||
          structure == 11||
          structure == 12) {lower_strip = 1;
          		    upper_strip = 61;
          		    lower_diff = -127;
          		    upper_diff = 68;}
          		    
 //////////////////////////////////////////////////
 // Hit-Strip Identification                     //
 //////////////////////////////////////////////////
 
 int hit_strip_3[28]   = {8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62}; 
 
 int hit_strip_8[29] = {66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122};
 
 int hit_strip_10[30]  = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60};  
 
 int hit_strip_12[30]  = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60};  
 
 int histo_number = 0;
 
 if(structure == 3)
    histo_number = 28;
 else if (structure == 8)
    histo_number = 29;
 else if (structure == 10 || structure == 12)
    histo_number = 30;
 
 upper_strip = 0; //Laura: removed "int"s here
 lower_strip = 0;
 lower_diff = 0;
 upper_diff = 0;
 
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
          		    
 //////////////////////////////////////////////////
 // Files and Histograms                         //
 //////////////////////////////////////////////////
 
 TFile* File = new TFile(Form("../data/Signal_Distribution_Structure_%d.root",structure),"RECREATE");
 TFile* OnFile = 0;
 TFile* OffFile = 0;
  
 TH2F* OnHisto[1];//Laura
 TH2F* OffHisto[1];//Laura

 OnFile = TFile::Open(Form("../data/Irradiated_Structure_Noise_Cut_%d_on.root",structure));
 if(structure == 3)
   TH2F* OnHisto[28];
 else if(structure == 8)
   TH2F* OnHisto[29];
 else if (structure == 10 || structure == 12)
   TH2F* OnHisto[30];
   
 OffFile = TFile::Open(Form("../data/Irradiated_Structure_Noise_Cut_%d_off.root",structure));
 if(structure == 3)
   TH2F* OffHisto[28];
 else if(structure == 8)
   TH2F* OffHisto[29];
 else if (structure == 10 || structure == 12)
   TH2F* OffHisto[30];

 //////////////////////////////////////////////////
 // Creating Histograms                          //
 //////////////////////////////////////////////////
 
 TH1F* DistributionOn = new TH1F(Form("distribution_on_struc_%d",structure),Form("Signal Distribution (on embedded), %s;signal - noise; number of entries",struc_ID[structure].c_str()), 250, -125.5, 124.5);
 
 TH1F* DistributionOff = new TH1F(Form("distribution_off_struc_%d",structure),Form("Signal Distribution (off embedded), %s;signal - noise; number of entries",struc_ID[structure].c_str()), 250, -125.5, 124.5);
 
 for (int i = 0; i < histo_number; i++)
 {
 	if(structure == 3)
	  hitstrip = hit_strip_3[i];
	else if(structure == 8)
	  hitstrip = hit_strip_8[i];
	else if(structure == 10)
	  hitstrip = hit_strip_10[i];
	else if(structure == 12)
	  hitstrip = hit_strip_12[i];
	  
 	OnFile->cd();
	OnHisto[i]   = (TH2F*)gDirectory->Get(Form("Hitmap_Noise_Cut_Structure_%d_Channel_%d", structure, hitstrip)); 
	
	OffFile->cd();
	OffHisto[i]  = (TH2F*)gDirectory->Get(Form("Hitmap_Noise_Cut_Structure_%d_Channel_%d", structure, hitstrip));
 }
	
	File->cd();
 
 double On            = 0;
 double Off           = 0;
 int    bin           = 0;
 int    binContentOn  = 0;
 int    binContentOff = 0;      
 
 
 //////////////////////////////////////////////////
 // Signal Distribution                          //
 //////////////////////////////////////////////////
  
 	for (int j = 0; j < histo_number; j++)// scann loop
 	{
 	
 		if(structure == 3)
		  hitstrip = hit_strip_3[j];
		else if(structure == 8)
		  hitstrip = hit_strip_8[j];
		else if(structure == 10)
		  hitstrip = hit_strip_10[j];
		else if(structure == 12)
		  hitstrip = hit_strip_12[j];
 		
 		for (int k = 0; k < 128; k++) // channel loop
 		{
 			for (int l = 0; l < 16; l++) // threshold loop 
 			{
	 			if (abs(hitstrip - k) > 5)
 				{
	 				
	 				On  = OnHisto[j]->GetBinContent(k+1,l+1);
	 				Off = OffHisto[j]->GetBinContent(k+1,l+1);
	 				
	 				DistributionOn->Fill(On);
	 				DistributionOff->Fill(Off);
	 				
	 			}
 				
 			} // threshold loop
 			
 		} // channel loop
 		
 	} // scann loop
  
 
 
 /////////////////////////////////////////////////
 
 TCanvas* c1 = new TCanvas("c1","c1", 1000,500);
 DistributionOn->Draw();
 DistributionOff->Draw("same");
 
 File->Write();

 }
 
