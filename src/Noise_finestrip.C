#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
//#include <typeinfo.h>
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

using namespace std;

#ifndef NOISE_finestrip_C
#define NOISE_finestrip_C

void Noise_finestrip()
  
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
  //gStyle->SetOptStat(0);
  
  // Numbering the Measurements 
  std::map<int,std::string> Data_ID;
  
  Data_ID[1] = "3_on";
  Data_ID[2] = "3_off";
  Data_ID[3] = "8_on";
  Data_ID[4] = "8_off";
  Data_ID[5] = "10_on";
  Data_ID[6] = "10_off";
  Data_ID[7] = "12_on";
  Data_ID[8] = "12_off";
  
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
  
  // Create File
  TFile* File = new TFile("../data/Noise_finestrip.root","RECREATE");
  TFile* RawData = 0;
  
  // Create Histograms
  
  TProfile2D* Noise_Histo[12];
  for(int i = 0; i < 12; i++)
    {
      int j = i+1;	
      Noise_Histo[i] = new TProfile2D(Form("Noise_Histo_%d",j), Form("Noise (%s); Channel ; Threshold in [mV]",struc_ID[j].c_str()),128,-0.5,127.5, 16,48.8,99.2);
    }
  TH2F* Raw_Histo = 0;
  
  /* for(int i = 0; i < 8; i++)
    {
      int j = i + 1;
      
      //cout<<"j = "<<j<<endl;
      */
      RawData = TFile::Open("../data/Irradiated_finestrip_Structure_3_on.root");
      
      for(int k = 0; k < 12; k++)
 	{
	  int l = k + 1;
	  
	  //cout<<"l = "<<l<<endl;
	  
	  Raw_Histo = (TH2F*)gDirectory->Get(Form("Average_Noise_Histo_%d",l));
	  Noise_Histo[k]->Add(Noise_Histo[k],Raw_Histo,1,1);
 	}
      
      RawData->Close();
      
  // }
 
  
  File->Write();
  
  File->Close();
  
}

#endif /* NOISE_finestrip_C */
