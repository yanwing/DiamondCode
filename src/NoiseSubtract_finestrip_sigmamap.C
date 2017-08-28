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

#ifndef NOISE_SUBTRACT_finestrip_C
#define NOISE_SUBTRACT_finestrip_C

void NoiseSubtract_finestrip(int structure, bool on)

{

  std::cout << "entering NoiseSubtract_finestrip.C" << std::endl;

 gROOT->ForceStyle();

 TH1::SetDefaultSumw2();
 TH2::SetDefaultSumw2();
 TProfile::SetDefaultSumw2();
 TProfile2D::SetDefaultSumw2();

 Double_t sum_hit;
 Double_t sum_noise;
 Double_t scale_factor;

 //////////////////////////////////////////////////
 // Colouring                                    //
 //////////////////////////////////////////////////
 const Int_t NRGBs = 5;
 const Int_t NCont = 255;

 Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
 Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
 Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
 Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
 
 Double_t Red[3] = {0.0, 1.0, 1.0};
 Double_t Green[3] = {0.0, 0.0, 1.0};
 Double_t Blue[3] = {1.0, 0.0, 1.0};
 Double_t Stops [3] = {0.0, 0.4, 1.0 };
 Double_t White [3] = {1.0, 1.0, 1.0};

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
 /*
 int hit_strip_on_3[29]   = {6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62}; 
 int hit_strip_off_3[56]  = {7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62};
 
 int hit_strip_on_8[30] = {64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122};
 int hit_strip_off_8[59] = {66,67,68,69,70,71,72,73,74,75,76,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123};
 
 int hit_strip_on_10[31]  = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62};  
 int hit_strip_off_10[61] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61};
 
 int hit_strip_on_12[31]  = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62};  
 int hit_strip_off_12[61] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62};
 */
 

 //////////////////////////////////////////////////
 // Create File & Histograms                     //
 //////////////////////////////////////////////////
 
 int histo_number = 0;//Laura
 TH2F* SNRmap[135];
 TH2F* Hitmap_Noise_Cut[135];//Laura (1 gives memory issues with array)
 /* 
 if(structure == 3 && !on)
   {
     TH2F* Hitmap_Noise_Cut[56];
     histo_number = 56;
   }
 else if(structure == 3 && on)
   {
     TH2F* Hitmap_Noise_Cut[29];
     histo_number = 29;
   }
 //else if((structure == 10 || structure == 12) && !on)
 else if(structure == 10 && !on) 
   {
     TH2F* Hitmap_Noise_Cut[61];
     histo_number = 61;
   }
 else if(structure == 12 && !on) 
   {
     TH2F* Hitmap_Noise_Cut[61];
     histo_number = 61;
   }//(split in 10 and 12 by Laura)
 else if(structure == 10 && on)
   {
     TH2F* Hitmap_Noise_Cut[31];
     histo_number = 31;
   }
 else if(structure == 12 && on)
   {
     TH2F* Hitmap_Noise_Cut[31];
     histo_number = 31;
   }
 else if (structure == 8 && !on)
   {
     TH2F* Hitmap_Noise_Cut[59];
     histo_number = 59;
   }
 else if (structure == 8 && on)
   {
     TH2F* Hitmap_Noise_Cut[30];
     histo_number = 30;
   }
 else
   {
     std::cout <<" . . . . . . . . - - - SOMETHING IS WEIRD! NONE OF THE DEFINED CASES! - - - . . . . . . . . "<< std::endl;
   }
 */
 std::string onoffstring = "on";
 if (!on) onoffstring = "off";
 
 // Create File
 TFile* File = new TFile(Form("../data/Irradiated_finestrip_Structure_Noise_Cut_%d_%s.root",structure, onoffstring.c_str()),"RECREATE");
 TFile* RawData = 0;
 TFile* NoiseData = 0;
 
 // Create Raw Data Histogram
 TH2F* raw_data_histo = 0;
 TH2F* raw_noise_histo = 0;
 int hit_strip = 0;

 	std::cout << "________________________Structure: " << structure << ", " << onoffstring.c_str() << std::endl;
	std::cout << "____________________histo_number = " << 135 << std::endl;


  	for(int k = 0; k < 134; k++)
   	{
	  /* hit_strip = 0;
	  
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
	    {
	      //std::cout << "hit_strip (before defining it via hit_strip_off_10[k]): " << hit_strip << std::endl;
	      //std::cout << "hit_strip_off_10[k] for k="<<k<< ": " << hit_strip_off_10[k] << std::endl;
	      hit_strip = hit_strip_off_10[k];
	      //std::cout << "hit_strip (after  defining it via hit_strip_off_10[k]): " << hit_strip << std::endl;
	      //std::cout << "hit_strip_off_10[1] = " << hit_strip_off_10[1] << std::endl;
	    }
	  else if(structure == 12 && on)
	    hit_strip = hit_strip_on_12[k];
	  else if(structure == 12 && !on)
	    hit_strip = hit_strip_off_12[k];
	  */

	  //defining Hitmap_Noise_Cut (incl. naming etc.) (Laura)

	  // Hitmap_Noise_Cut[k] = new TH2F(Form("Hitmap_Noise_Cut_Structure_%d_Channel_%d",structure,hit_strip), Form("Hitmap (%s), beam on channel %d, background subtracted; Channel ; Threshold in [mV]",struc_ID[structure].c_str(),hit_st\rip),128,-0.5,127.5, 16,48.8,99.2);
	  Hitmap_Noise_Cut[k] = new TH2F(Form("Hitmap_Noise_Cut_Structure_%d_Beamposition_%d",structure,k), Form("Hitmap (%s), beam position  %d, background subtracted; Channel ; Threshold in [mV]",struc_ID[structure].c_str(),k),128,-0.5,127.5, 16,48.8,99.2);
	  SNRmap[k] = new TH2F(Form("SNRmap_Structure_%d_Beampositio_%d",structure,k), Form("Hitmap (%s), beam position  %d, background subtracted; Channel ; Sigma",struc_ID[structure].c_str(),k),128,-0.5,127.5, 20,0,2);

	  //cout << "defining Hitmap Noise Cut" << k << endl;

  	}//end of loop for defining Hitmap_Noise_Cut
  	
	//std::cout << "IN BETWEEN THE LOOPS: hit_strip_off_10[0] = " << hit_strip_off_10[0] << std::endl;
	//std::cout << "IN BETWEEN THE LOOPS: hit_strip_on_10[0]  = " << hit_strip_on_10[0] << std::endl;

	for(int m = 0; m < 134  ; m++)
	  {
	    //std::cout << "IN THE 2nd LOOP: hit_strip_on_10[0]  = " << hit_strip_on_10[0] << std::endl;//ok
	    //std::cout << "IN THE 2nd LOOP: hit_strip_off_10[0] = " << hit_strip_off_10[0] << std::endl;//weird
	    /*
	    hit_strip = 0; //(Laura)

	    //	    std::cout << "loop (for filling Hitmap_Noise_Cut) defining 'hit_strip'..... counter: " << m << std::endl;
	   
	    if(structure == 3 && on)
	      hit_strip = hit_strip_on_3[m];
	    else if(structure == 3 && !on)
	      hit_strip = hit_strip_off_3[m];
	    
	    else if(structure == 8 && on)
	      hit_strip = hit_strip_on_8[m];
	    else if(structure == 8 && !on)
	      hit_strip = hit_strip_off_8[m];
	    
	    else if(structure == 10 && on){
	      hit_strip = hit_strip_on_10[m];
	    }
	    else if(structure == 10 && !on)
	      {
		//std::cout << "hit_strip (before defining it via hit_strip_off_10[m]): " << hit_strip << std::endl;
		//std::cout << "hit_strip_off_10[m] for m="<<m<< ": " << hit_strip_off_10[m] << std::endl;
		hit_strip = hit_strip_off_10[m];
		//std::cout << "hit_strip (after  defining it via hit_strip_off_10[m]): " << hit_strip << std::endl;
		//std::cout << "hit_strip_off_10[m+1] = " << hit_strip_off_10[m+1] << std::endl;
	      }
	    else if(structure == 12 && on)
	      hit_strip = hit_strip_on_12[m];
	    else if(structure == 12 && !on)
	      hit_strip = hit_strip_off_12[m];
	    */
	    //std::cout << "hit_strip = " << hit_strip << std::endl;
	    //cout << "open finestrip analysis file " << m << endl;
	     RawData = TFile::Open(Form("../data/Irradiated_finestrip_Structure_%d_%s.root",structure, onoffstring.c_str()));
	     // cout<< "write file " << m << endl;
	     RawData->cd(Form("Irradiated_finestrip_Structure_%d_%s",structure, onoffstring.c_str()));
	     //cout<< "open hitmap " << m <<endl;
	    raw_data_histo = (TH2F*)gDirectory->Get(Form("Hitmap_Structure_%d_Beamposition_%d",structure,m));
	    
	    // cout<< "open Noise file " << m << endl;
	    NoiseData = TFile::Open("../data/Noise.root");
	    NoiseData->cd();  	
	    raw_noise_histo = (TH2F*)gDirectory->Get(Form("Noise_Histo_%d",structure));
	    
	    
	    //std::cout << "entries in data histo for hit strip " << hit_strip << ": " << raw_data_histo->GetEntries() << std::endl;
	    //std::cout << "adding raw_data_histo, raw_noise_histo to obtain Hitmap_Noise_Cut " << m << std::endl;
	    
	    //Noramlization( Sam)
	    sum_hit = 0;
	    sum_noise = 0;
	    scale_factor = 1;
	    for (int ch = 1; ch < 129; ch++)
	      {
		sum_hit += raw_data_histo -> GetBinContent(ch,1);
		//cout<< "Bin Count of Hitmap : " << ch << ":"<< sum_hit << endl;
		sum_noise += raw_noise_histo -> GetBinContent(ch,1);
		//cout<< "Bin Count of Noise Map : " << ch << ":" << sum_noise << endl;
	      }
	    scale_factor = sum_hit/sum_noise;
	    // cout << "scale_factor = "<< scale_factor <<endl;
	    //Actually filling Hitmap_Noise_Cut
	    //cout << "Substract Noise from Hitmap " << m << endl;
	    Hitmap_Noise_Cut[m]->Add(raw_data_histo,raw_noise_histo,1,-scale_factor);
	    //cout << "Close file " << m << endl;
	    //signal Identification
 	    for (int cht = 1;cht <129; cht++)
	      {
		Double_t beam_hit[15];
		Double_t noise_error[15];
		//cout<< "define noise error" <<endl;
		for (int nq  = 1; nq < 11; nq++)
		  {
		    beam_hit[nq]=0;
		    beam_hit[nq] = Hitmap_Noise_Cut[m]->GetBinContent(cht,nq);
		    //cout<< nq <<endl;
		    noise_error[nq]=0;
		    noise_error[nq]= raw_noise_histo->GetBinContent(cht,nq);
		  }
		//cout<<"extracted noise bin for channel " <<cht<< "beam" <<m<<endl;
		for (int sig = 1 ; sig < 21; sig++)
		  {
		    Double_t count;
		    count=0;
		    for (int q = 1; q < 9; q++)
		      {
			if (beam_hit[q] > (sig*0.1)*noise_error[q])
			  {
			    if (beam_hit[q+1] > (sig*0.1)*noise_error[q+1])
			      {
				 if (beam_hit[q+2] > (sig*0.1)*noise_error[q+2])
				{
				    count = q;
				    //	    cout << "signal found in channel  " << cht-1 << " with  beam postion " << m << ", SNR = " <<sig*0.1 << ";" <<count<<endl;
				  }
			      else continue;
			       }
			    else continue;
			  }
			else continue;
		      }
		    SNRmap[m]->SetBinContent(cht,sig,(50+(count-1)*3.125));  
		  
		    if(sig == 2)
		      {
			if (count > 0) 
			  {
			    cout << "signal found in channel  " << cht-1 << " with  beam postion " << m << ", SNR = " <<sig*0.1 << ";" <<count<<endl;
			  }
		      }
		  }
	      }   
 
	    //adding a TCanvas that can be printed out as PDFs (Laura)
	    TCanvas *c = new TCanvas;
	    //Hitmap_Noise_Cut[m]->SetMinimum(0);
	    Hitmap_Noise_Cut[m]->SetMaximum(5000);
	    Hitmap_Noise_Cut[m]->Draw("COLZ");
	    c->Print(Form("../plots/Hitmap_Noise_Cut_%d.png",m));
	    delete c;
	    TCanvas *d =new TCanvas;
	    SNRmap[m]->SetMinimum(50);
	    SNRmap[m]->SetMaximum(80);
	    SNRmap[m]->GetZaxis()->SetTitle("Threshold");
	    SNRmap[m]->Draw("COLZ");
	    d->Print(Form("../plots/Sigma_%d.png",m));
	    delete d;

	    RawData->Close();
	    NoiseData->Close();
	  }//end of loop for filling Hitmap_Noise_Cut
	
	File->Write();
	File->Close();
	
}

#endif /* NOISE_SUBTRACT_finestrip_C */
