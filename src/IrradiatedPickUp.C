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

void IrradiatedPickUp(int structure)

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
 
 int upper_strip = 0;
 int lower_strip = 0;
 int lower_diff = 0;
 int upper_diff = 0;
 
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
 
 TFile* File = new TFile(Form("../data/Irradiated_PickUp_Structure_%d.root",structure),"RECREATE");
 TFile* OnFile = 0;
 TFile* OffFile = 0;
 TFile* NoiseFile = 0;
 TFile* DeConvFile = 0;
  
 TH2F* OnHisto[1];//Laura
 TH2F* OffHisto[1];//Laura
 TH2F* HistoDiff[1];//Laura
 
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
 
 NoiseFile = TFile::Open("../data/Noise.root");
 TH2F* NoiseHisto = (TH2F*)gDirectory->Get(Form("../data/Noise_Histo_%d", structure));
 
 DeConvFile = TFile::Open(Form("../data/XTalk_Irradiated_Structure_%d.root",structure));
 TH1F* DeConvHisto = (TH1F*)gDirectory->Get(Form("../data/deconv_cross_talk_struc_%d", structure));
 
 File->cd();
 
 if(structure == 3)
   TH2F* HistoDiff[28];
 else if(structure == 8)
   TH2F* HistoDiff[29];
 else if (structure == 10 || structure == 12)
   TH2F* HistoDiff[30];
 
 TH2F* PickUpRatioHisto = new TH2F(Form("PickUp_Ratio_Histo_structure_%d",structure), Form("Ratio of negative and positive signal rows (%s); #sigma; channel diff.; ratio (off-on)/(on-off)",struc_ID[structure].c_str()), 64, -128.5, 127.5, 15, -0.05, 1.45);
 
 TH2F* AntiPickUpHisto = new TH2F(Form("Anti_PickUp_Histo_structure_%d",structure), Form("Number of negative signal rows (%s); #sigma; channel diff.; Number of negative signal rows",struc_ID[structure].c_str()), 64, -128.5, 127.5, 15, -0.05, 1.45);
 
 TH2F* PickUpHisto = new TH2F(Form("PickUp_Histo_structure_%d",structure), Form("Number of positive signal rows (%s); #sigma; channel diff.; Number of positive signal rows",struc_ID[structure].c_str()), 64, -128.5, 127.5, 15, -0.05, 1.45);
 
 //////////////////////////////////////////////////
 // Creating Histograms                          //
 //////////////////////////////////////////////////
 
 int hitstrip = 0;
 
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
	
	File->cd();
	
	HistoDiff[i] = new TH2F(Form("HistoDiff_structure_%d_channel_%d",structure,hitstrip), Form("Signaldifference On/Off-Embedded %s",struc_ID[structure].c_str()), 128, -0.5, 127.5, 16,48.8,99.2);
 }
 
 
 File->cd();
 
 for (int i = 0; i < histo_number; i++)
 {	
 	HistoDiff[i]->Add(OnHisto[i],OffHisto[i],1,-1);
 }
 
 double SignalDiff = 0;
 double on = 0;
 double off = 0;
 double onErr = 0;
 double offErr = 0;
 double sigma = 0;
 double noise = 0; 
 double PickUpWeight = 0;
 int ChannelDiff = 0;
 int Counter = 0;
 int AntiCounter = 0;
 int binY = 0;
 
 //////////////////////////////////////////////////
 // Pick Up                                      //
 //////////////////////////////////////////////////
 
 for (int i = 0; i< histo_number; i++) // loop over scanns
 {
	if(structure == 3)
	  hitstrip = hit_strip_3[i];
	else if(structure == 8)
	  hitstrip = hit_strip_8[i];
	else if(structure == 10)
	  hitstrip = hit_strip_10[i];
	else if(structure == 12)
	  hitstrip = hit_strip_12[i];
		 
 	for (int l = 0; l<=15; l++) // loop over parts of sigma
	{
		 
		for (int k = 0; k<128; k++)  // loop over channels
		{
			for (int j = 0; j<16; j++)  // loop over thresholds
			{
				// channel number difference of strip with possible pick up and strip with beam
				ChannelDiff = hitstrip - k; 
			
				// signal difference of on and off embedded		
				SignalDiff = HistoDiff[i]->GetBinContent(k+1,j+1);
			
				// signal and error with beam on embedded
				on = OnHisto[i]->GetBinContent(k+1,j+1);
				onErr = OnHisto[i]->GetBinError(k+1,j+1);
			
				// signal and error with beam off embedded
				off = OffHisto[i]->GetBinContent(k+1,j+1);
				offErr = OffHisto[i]->GetBinError(k+1,j+1);
			
				// significance
				noise = NoiseHisto->GetBinContent(k+1,j+1); //0.1*l*sqrt(onErr*onErr + offErr*offErr);
				
				sigma = 0.1*l*sqrt(2)*noise;
				
				// condition if the signal difference for one channel is higher than zero 
				// for three or more thresholds in a row
				if (abs(SignalDiff) > sigma && abs(ChannelDiff) > 1)
				{
		
					if (Counter == 0 && SignalDiff <= 0)
					{
						continue;
					}
					else if (Counter == 0 && SignalDiff > 0)
					{
						Counter++;
						continue;
					}
					else if (Counter > 0 && Counter < 3 && SignalDiff <= 0)
					{
						Counter = 0;
						continue;
					}
					else if (Counter > 0 && Counter < 3 && SignalDiff > 0)
					{	
						Counter++;
						continue;
					}
					else if (Counter >= 3 && SignalDiff > 0)
					{
						Counter++;
						continue;
					}
					else if (Counter >= 3 && SignalDiff <= 0)
					{
						// weight to measure the significance of the posiible pick up
						int d = DeConvHisto->GetXaxis()->FindBin(i);
						if (d!=0)
						PickUpWeight = 1.0/d;//fabs((on-off)/noise)
						// PickUpWeight =threshold_steps[j]*(threshold_steps[j] - threshold_steps[j-Counter]);
						PickUpHisto->Fill(ChannelDiff,0.1*l,PickUpWeight);//(ChannelDiff,0.1*l,1)
						Counter = 0;
					}
			
				} // sigma condition
				
			} // threshold loop
		
			Counter     = 0;
		
		} // channel loop
		
	} // sigma loop
	
 } // scan loop
 
 //////////////////////////////////////////////////
 // Anti Pick Up                                 //
 //////////////////////////////////////////////////
 
 for (int i = 0; i< 16; i++) // loop over scanns
 {
	if(structure == 3)
	  hitstrip = hit_strip_3[i];
	else if(structure == 8)
	  hitstrip = hit_strip_8[i];
	else if(structure == 10)
	  hitstrip = hit_strip_10[i];
	else if(structure == 12)
	  hitstrip = hit_strip_12[i];
		 
 	for (int l = 0; l<=15; l++) // loop over parts of sigma
	{
		 
		for (int k = 0; k<128; k++)  // loop over channels
		{
			for (int j = 0; j<16; j++)  // loop over thresholds
			{
				ChannelDiff = hitstrip - k;
						
				SignalDiff = HistoDiff[i]->GetBinContent(k+1,j+1);
			
				on = OnHisto[i]->GetBinContent(k+1,j+1);
				onErr = OnHisto[i]->GetBinError(k+1,j+1);
			
				off = OffHisto[i]->GetBinContent(k+1,j+1);
				offErr = OffHisto[i]->GetBinError(k+1,j+1);
			
				noise = NoiseHisto->GetBinContent(k+1,j+1); //0.1*l*sqrt(onErr*onErr + offErr*offErr);
				
				sigma = 0.1*l*sqrt(2)*noise;
				
				if (abs(SignalDiff) > sigma && abs(ChannelDiff) > 1)
				{
				
					if (AntiCounter == 0 && SignalDiff > 0)
					{
						continue;
					}
					else if (AntiCounter == 0 && SignalDiff <= 0)
					{
						AntiCounter++;
						continue;
					}
					else if (AntiCounter > 0 && AntiCounter < 3 && SignalDiff > 0)
					{
						AntiCounter = 0;
						continue;
					}
					else if (AntiCounter > 0 && AntiCounter < 3 && SignalDiff <= 0)
					{	
						AntiCounter++;
						continue;
					}
					else if (AntiCounter >= 3 && SignalDiff <= 0)
					{
						AntiCounter++;
						continue;
					}
					else if (AntiCounter >= 3 && SignalDiff > 0)
					{
						//PickUpWeight = abs((on-off)/noise);
						//AntiPickUpHisto->Fill(ChannelDiff,0.1*l,1);
						int d = DeConvHisto->GetXaxis()->FindBin(i);
						if (d!=0)
						PickUpWeight = 1.0/d;//fabs((on-off)/noise)
						AntiPickUpHisto->Fill(ChannelDiff,0.1*l,PickUpWeight);//(ChannelDiff,0.1*l,1)
						AntiCounter = 0;
					}
			
				} // sigma condition
				
			} // threshold loop
		
			AntiCounter = 0;
		
		} // channel loop
		
	} // sigma loop
	
 } // scan loop
 
 //////////////////////////////////////////////////
 // Storing File                                 //
 //////////////////////////////////////////////////
 
 File->cd();
 File->Write();
 
 }
  
