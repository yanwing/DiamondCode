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

#ifndef IRRADIATED_TEST_BEAM_ANALYSIS_C
#define IRRADIATED_TEST_BEAM_ANALYSIS_C

void irradiated_test_beam_analysis_03(int structure, bool on)
  
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
  
  // no statistic box
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
  // File Identification                          //
  //////////////////////////////////////////////////
  
  std::map<int,int> strcture_firstfile;
  std::map<int,int> strcture_lastfile;
  
  if(!on)
    {
      strcture_firstfile[3]  = 31;
      strcture_lastfile[3]   = 86;
      strcture_firstfile[10] = 116;
      strcture_lastfile[10]  = 176;
      strcture_firstfile[12] = 208;
      strcture_lastfile[12]  = 268;
      strcture_firstfile[8]  = 301;
      strcture_lastfile[8]   = 359;
    }
  
  if(on)
    {
      strcture_firstfile[3]  = 87;
      strcture_lastfile[3]   = 115;
      strcture_firstfile[10] = 177;
      strcture_lastfile[10]  = 207;
      strcture_firstfile[12] = 269;
      strcture_lastfile[12]  = 299;
      strcture_firstfile[8]  = 360;
      strcture_lastfile[8]   = 389;
    }
  
  //////////////////////////////////////////////////
  // Strip Identification                         //
  //////////////////////////////////////////////////
  
  std::map<int,int> strip_ID;
  
  strip_ID[1]  = 256;
  strip_ID[2]  = 128;
  strip_ID[3]  = 0;
  strip_ID[4]  = 256;
  strip_ID[5]  = 128;
  strip_ID[6]  = 0;
  strip_ID[7]  = 384;
  strip_ID[8]  = 512;
  strip_ID[9]  = 640;
  strip_ID[10] = 384;
  strip_ID[11] = 512;
  strip_ID[12] = 640;
  
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
  
  TH2F* Hitmap[61];//Laura
  
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
  //else if((structure == 10 || structure == 12) && !on) //Laura: separated 10, 12 cases
  else if(structure == 10 && !on)
    {
      TH2F* Hitmap[61];
      histo_number = 61;
    }
  else if(structure == 12 && !on) //Laura: separated 10, 12 cases
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
  
  std::string onoffstring = "on";
  if (!on) onoffstring = "off";
  
  // Create File
  TFile* File = new TFile(Form("../data/Irradiated_Structure_%d_%s.root",structure, onoffstring.c_str()),"RECREATE");
  TFile* RawData = 0;
  
  // Create Raw Data Histogram
  TH2F* raw_data_histo = 0;
  TH2F* raw_noise_histo = 0;
  
  // Create Noise Histogram
  TProfile2D* Average_Noise_Histo[12];
  
  for(int i = 0; i < 12; i++)
    {
      int j = i+1;	
      Average_Noise_Histo[i] = new TProfile2D(Form("Average_Noise_Histo_%d",j), Form("Noise (%s); Channel ; Threshold in [mV]",struc_ID[j].c_str()),128,-0.5,127.5, 16,48.8,99.2);
    }
  
  // Create Hitmaps
  
  //Laura (with some help from Root forum):
  TDirectory *dir = File->GetDirectory(Form("Irradiated_Structure_%d_%s",structure, onoffstring.c_str()));
  if (dir) {
    dir->cd();
  } else {
    //fprintf(stderr,"Missing directory analyzeHiMassTau\n");
    //}
    //if(! File->cd(Form("Irradiated_Structure_%d_%s",structure, onoffstring.c_str()))){
    //cout << "Directory " << Form("Irradiated_Structure_%d_%s",structure, onoffstring.c_str()) << " could not be accessed (yet)" << endl;
    File->mkdir(Form("Irradiated_Structure_%d_%s",structure, onoffstring.c_str()));
    File->cd(Form("Irradiated_Structure_%d_%s",structure, onoffstring.c_str()));
  }
  
  if(structure == 3 || structure == 10 || structure == 12 || structure == 8)
    { 
      for(int k = 0; k < histo_number; k++)
	{
	  int hit_strip = 0;
	  
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
	  
	  //Hitmap[k]= new TH2F(Form("Hitmap_Structure_%d_Channel_%d",structure,hit_strip), Form("Hitmap (%s), beam on channel %d; Channel ; Threshold in [mV]",struc_ID[i].c_str(),hit_strip),128,-0.5,127.5, 16,48.8,99.2);
	  Hitmap[k]= new TH2F(Form("Hitmap_Structure_%d_Channel_%d",structure,hit_strip), Form("Hitmap (%s), beam on channel %d; Channel ; Threshold in [mV]",struc_ID[structure].c_str(),hit_strip),128,-0.5,127.5, 16,48.8,99.2); //Laura
	  
	}
      
    }
  //else
  //  continue;
  //Laura: this continue statement was outside a loop and thus didn't work. Not needed at all?
  
  
  //////////////////////////////////////////////////
  // Read in Data                                 //
  //////////////////////////////////////////////////
  
  // Filling the Hitmaps
  for(int x = 1; x <= 12; x++)
    {
      
      int q = x - 1;
      
      for(int l = strcture_firstfile[structure]; l <= strcture_lastfile[structure]; l++) // loop over files
	{
	  int m = l - strcture_firstfile[structure]; // m = number of histogram
	  int channel = 0; // the 128 channels of a strucure	 
	  double threshold_step = 48.8;
	  
	 RawData = TFile::Open(Form("../raw_data/ABCN250_at_150V/strun287_%d.root",l));
	 
 	 // get the right data stream
	 if(x == 1 || x == 2 || x == 3 || x == 10 || x == 11 || x == 12)
	   raw_data_histo = (TH2F*)gDirectory->Get("h_scan0;1");
	 else if(x == 4 || x == 5 || x == 6 || x == 7 || x == 8 || x == 9)
	   raw_data_histo = (TH2F*)gDirectory->Get("h_scan1;1");
	 
	 int last_channel = strip_ID[x] + 128; // last channel of structure
	 double bin_content = 0;
	 double bin_error = 0;
	 
	 for(int n = 0; n < 768; n++) // loop over channels
	   { 	  	
	     for(int p = 0; p < 16; p++) // loop over thresholds
	       {
		 if(n >= strip_ID[x] && n < last_channel) // the 128 channels of a structure
		   {
		     channel = 128 - (last_channel - n);
		     
		     bin_content = raw_data_histo->GetBinContent(n,p);
		     if(bin_content != 0)
		       {
			 bin_error = 1/TMath::Sqrt(bin_content);
		       }
		     else
		       {
			 bin_error = 0.0;
		       }
		     // fill hitmaps
		     if(x == structure)
		       {
			 Hitmap[m]->Fill(channel,threshold_step,bin_content);
			 Hitmap[m]->SetBinError(channel,p,bin_error);
		       }
		     
		     // fill noise-maps
		     else
		       {
			 Average_Noise_Histo[q]->Fill(channel,threshold_step,bin_content);
			 Average_Noise_Histo[q]->SetBinError(channel,p,bin_error);
		       }
		     //}
		     //else continue;
		   }
		 else continue;
		 
		 threshold_step += 3.2; 
		 
	       } // close threshold loop
	     
	     threshold_step = 48.8;
	     
	   } // close channel loop
	 
	 RawData->Close();
	 
	} // close file loop
      
    }
  
  // write file
  File->Write();
  cout << "written to output file for structure " << structure << " and on-/off-embedded (1/0): " << on << endl;
  // close file
  File->Close();
  cout << "output file closed; end of script" << endl;
  
}

#endif /* IRRADIATED_TEST_BEAM_ANALYSIS_C */
