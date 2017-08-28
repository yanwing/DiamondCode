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

using namespace std;

#ifndef SIGNAL_NOISE_RATIO_IRRADIATED_C
#define SIGNAL_NOISE_RATIO_IRRADIATED_C

void SignalNoiseRatioIrradiated(int structure, bool on)
{

 std::string onoffstring = "on";
 if (!on) onoffstring = "off";

 TH1::SetDefaultSumw2();
 TH2::SetDefaultSumw2();
 TProfile::SetDefaultSumw2();
 TProfile2D::SetDefaultSumw2();
 
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
 
 int hit_strip_on_3[29]   = {6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62}; 
 int hit_strip_off_3[56]  = {7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62};
 
 int hit_strip_on_8[30] = {64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122};
 int hit_strip_off_8[59] = {66,67,68,69,70,71,72,73,74,75,76,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123};
 
 int hit_strip_on_10[31]  = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62};  
 int hit_strip_off_10[61] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61};

 int hit_strip_on_12[31]  = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62};  
 int hit_strip_off_12[61] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62};

//////////////////////////////////////////////////
 // Create File & Histograms                     //
 //////////////////////////////////////////////////
 
 int histo_number = 0;
 
 if(structure == 3 && !on)
 {
 TH2F* Hitmap[56];
 histo_number = 56;
 }
 else if(structure == 3 && on)
 {
 TH2F* Hitmap[29];
 histo_number = 29;
 }
 else if((structure == 10 || structure == 12) && !on)
 {
 TH2F* Hitmap[61];
 histo_number = 61;
 }
 else if(structure == 10 && on)
 {
 TH2F* Hitmap[31];
 histo_number = 31;
 }
 else if (structure == 12 && on)
 {
  TH2F* Hitmap[31];
  histo_number = 31;
 }
  else if (structure == 8 && !on)
 {
  TH2F* Hitmap[59];
  histo_number = 59;
 }
 else if (structure == 8 && on)
 {
  TH2F* Hitmap[30];
  histo_number = 30;
 }
 
 TFile* File = new TFile(Form("../data/SignalNoiseRatio_Irradiated_Structure_%d_%s.root",structure,onoffstring.c_str()),"RECREATE");
 
 TFile* DataFile = 0; 
 DataFile = TFile::Open(Form("../data/Irradiated_Structure_%d_%s.root",structure, onoffstring.c_str()));//Laura: added onoffstring.c_str()
                                                                                                        //which ones were used before?
 //DataFile = TFile::Open(Form("../data/Irradiated_Structure_Noise_Cut_%d_%s.root",structure,onoffstring.c_str()));
 TFile* NoiseFile = 0;
 NoiseFile = TFile::Open("../data/Noise.root");
 
 File->cd();
 
 TH2F* DataHisto = 0;
 TProfile* NoiseProfile = new TProfile(Form("noise_profile_structure_%d",structure), Form("Noise Profile; threshold in [mV]"), 16, 48.8, 99.2);
 TProfile* SignalProfile = new TProfile(Form("sinal_profile_structure_%d",structure), Form("Signal Profile; threshold in [mV]"), 16, 48.8, 99.2);
 TProfile* SignalNoiseRatio = new TProfile(Form("signal_noise_ratio_structure_%d",structure), Form("Signal-To-Noise-Ratio (%s); threshold in [mV]; SNR", struc_ID[structure].c_str()), 16, 48.8, 99.2);
 
 NoiseFile->cd();
 TH2F* NoiseHisto = (TH2F*)gDirectory->Get(Form("../data/Noise_Histo_%d", structure));
 
 double weight = 0;
 double weight1= 0;
 double weight2= 0;
 double noise  = 0;
 double noise1 = 0;
 double noise2 = 0;
 double data   = 0;
 double data1  = 0;
 double data2  = 0;
 double ratio  = 0;
 double threshold_step = 48.8;
 int hit_strip = 0;
 
 for(int i = 0; i < 16; i++)
 {
 	for(int k = 0; k < histo_number; k++)
	{
	  
 	 if(structure == 3 && on)
   	 hit_strip = hit_strip_on_3[k];
   	 else if(structure == 3 && !on)
   	 hit_strip = hit_strip_off_3[k];
   	 
   	 else if(structure == 8 && on)
   	 hit_strip = hit_strip_on_8[k];
   	 else if(structure == 8 && !on)
   	 hit_strip = hit_strip_off_8[k];
   	    	 
   	 else if(structure == 10 && on)
   	 hit_strip = hit_strip_on_10[k];
   	 else if(structure == 10 && !on)
   	 hit_strip = hit_strip_off_10[k];
   	 
   	 else if(structure == 12 && on)
   	 hit_strip = hit_strip_on_12[k];
   	 else if(structure == 12 && !on)
   	 hit_strip = hit_strip_off_12[k];
   	 
   	 noise  = NoiseHisto->GetBinContent(hit_strip+1,i+1);
   	 noise1 = NoiseHisto->GetBinContent(hit_strip,i+1);
   	 noise2 = NoiseHisto->GetBinContent(hit_strip+2,i+1); 		 
 	 NoiseProfile->Fill(threshold_step,noise);
	 
	 DataFile->cd(Form("Irradiated_Structure_%d_%s",structure, onoffstring.c_str()));
	 //DataFile->cd();
	 	 
	 DataHisto = (TH2F*)gDirectory->Get(Form("Hitmap_Structure_%d_Channel_%d",structure,hit_strip));
	 //DataHisto = (TH2F*)gDirectory->Get(Form("Hitmap_Noise_Cut_Structure_%d_Channel_%d",structure,hit_strip));
	 	
	 data  = DataHisto->GetBinContent(hit_strip+1,i+1);
	 data1 = DataHisto->GetBinContent(hit_strip,i+1);
	 data2 = DataHisto->GetBinContent(hit_strip+2,i+1);
	 
	 if(noise1 != 0 && i == 0)
	 weight1 = data1/noise1;
	 
	 if(noise2 != 0 && i == 0)
	 weight2 = data2/noise2;
	 
	 weight = (weight1 + weight2)/2;
	 
	 if(noise != 0)
	 ratio = data/(noise*weight);
	 //ratio = data/(noise);
	 
	 SignalProfile->Fill(threshold_step,data);
	 
	 SignalNoiseRatio->Fill(threshold_step,ratio);
	 }
	 
	 threshold_step += 3.2;
 }
 
 TCanvas* c1 = new TCanvas("c1","c1", 1000,500);
 c1->SetLogy();
 SignalNoiseRatio->Draw("");
 //SignalNoiseRatio->GetYaxis()->SetRangeUser(0.001,101);
 
 File->Write();
 
}

#endif /* SINGAL_NOISE_RATIO_IRRADIATED_C */
