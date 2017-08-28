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
#include "THStack.h"
#include "TLatex.h"
using namespace std;

#ifndef NOISE_SUBTRACT_finestrip_C
#define NOISE_SUBTRACT_finestrip_C

void NoiseSubtract_finestrip(int structure, bool on, int bin, int PUcut,int pickup_crossing)

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
 
 ///////////////////////////////////////
 ///// Pick Up Channel elimaination/////
 ///////////////////////////////////////
 cout<<"set number of pick up channel"<<endl;
 int beam_channel[135] = {7,7,8,8,9,9,10,10,10,11,11,12,12,12,13,13,14,14,15,15,15,16,16,17,17,18,18,18,19,19,20,20,20,21,21,22,22,23,23,23,24,24,25,25,26,26,26,27,27,27,28,28,29,29,30,30,31,31,31,32,32,33,31,32,30,30,31,29,30,28,29,29,28,29,27,28,26,26,27,25,26,24,25,23,23,24,22,23,21,21,22,21,22,20,21,19,19,20,18,19,17,18,16,16,17,15,16,14,14,15,13,14,12,13,11,11,12,10,11,9,10,10,9,10,8,9,7,8,8,6,7,5,6,4,6};
 cout<<"set pick up channel number "<<endl;

 int hit_cht[135]={1,1,2,1,2,1,2,1,1,2,1,2,1,1,2,1,2,1,2,1,1,2,1,2,1,2,1,1,2,1,2,1,1,2,1,2,1,2,1,1,2,1,2,1,2,1,1,2,1,2,1,2,2,1,2,1,2,1,1,2,1,2,1,2,1,1,2,1,2,1,2,2,1,2,1,2,1,1,2,1,2,1,2,1,1,2,1,2,1,1,2,1,2,1,2,1,1,2,1,2,1,2,2,1,2,1,2,1,1,2,1,2,1,2,1,1,2,1,2,1,2,2,1,2,1,2,1,2,2,1,2,1,2,1,2};

 
 std::array<std::vector<int>, 135> beam;
 beam[0]={0,1,2,3,4,5,6};
 beam[1]={0,1,2,3,4,5,6};
 beam[2]={0,1,2,3,4,5,6,7};
 beam[3]={0,1,2,3,4,5,6,7};
 beam[4]={0,1,2,3,4,5,6,7,8};
 beam[5]={0,1,2,3,4,5,6,7,8};
 beam[6]={0,1,2,3,4,5,6,7,8,9};
 beam[7]={0,1,2,3,4,5,6,7,8,9};
 beam[8]={0,1,2,3,4,5,6,7,8,9};
 beam[9]={0,1,2,3,4,5,6,7,8,9,10};
 beam[10]={0,1,2,3,4,5,6,7,8,9,10};
 beam[11]={0,1,2,3,4,5,6,7,8,9,10,11};
 beam[12]={0,1,2,3,4,5,6,7,8,9,10,11};
 beam[13]={0,1,2,3,4,5,6,7,8,9,10,11};
 beam[14]={0,1,2,3,4,5,6,7,8,9,10,11,12};
 beam[15]={0,1,2,3,4,5,6,7,8,9,10,11,12};
 beam[16]={0,1,2,3,4,5,6,7,8,9,10,11,12,13};
 beam[17]={0,1,2,3,4,5,6,7,8,9,10,11,12,13};
 beam[18]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
 beam[19]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
 beam[20]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
 beam[21]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
 beam[22]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
 beam[23]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
 beam[24]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
 beam[25]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
 beam[26]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
 beam[27]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
 beam[28]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
 beam[29]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
 beam[30]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
 beam[31]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
 beam[32]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
 beam[33]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
 beam[34]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
 beam[35]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
 beam[36]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
 beam[37]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
 beam[38]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
 beam[39]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
 beam[40]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
 beam[41]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
 beam[42]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
 beam[43]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
 beam[44]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
 beam[45]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
 beam[46]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
 beam[47]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
 beam[48]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
 beam[49]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
 beam[50]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27};
 beam[51]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27};
 beam[52]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28};
 beam[53]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28};
 beam[54]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
 beam[55]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
 beam[56]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
 beam[57]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
 beam[58]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
 beam[59]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
 beam[60]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
 beam[61]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
 beam[62]={2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
 beam[63]={2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33};
 beam[64]={4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33};
 beam[65]={4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33};
 beam[66]={4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};
 beam[67]={6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};
 beam[68]={6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35};
 beam[69]={8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35};
 beam[70]={8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36};
 beam[71]={8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36};
 beam[72]={9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36};
 beam[73]={9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37};
 beam[74]={11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37};
 beam[75]={11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38};
 beam[76]={13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38};
 beam[77]={13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38};
 beam[78]={13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
 beam[79]={15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
 beam[80]={15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40};
 beam[81]={17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40};
 beam[82]={17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41};
 beam[83]={19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41};
 beam[84]={19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41};
 beam[85]={19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42};
 beam[86]={21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42};
 beam[87]={21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43};
 beam[88]={23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43};
 beam[89]={23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43};
 beam[90]={23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44};
 beam[91]={24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44};
 beam[92]={24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45};
 beam[93]={26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45};
 beam[94]={26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46};
 beam[95]={28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46};
 beam[96]={28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46};
 beam[97]={28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47};
 beam[98]={30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47};
 beam[99]={30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48};
 beam[100]={32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48};
 beam[101]={32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
 beam[102]={34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
 beam[103]={34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49};
 beam[104]={34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50};
 beam[105]={36,37,38,39,40,41,42,43,44,45,46,47,48,49,50};
 beam[106]={36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51};
 beam[107]={38,39,40,41,42,43,44,45,46,47,48,49,50,51};
 beam[108]={38,39,40,41,42,43,44,45,46,47,48,49,50,51};
 beam[109]={38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
 beam[110]={40,41,42,43,44,45,46,47,48,49,50,51,52};
 beam[111]={40,41,42,43,44,45,46,47,48,49,50,51,52,53};
 beam[112]={42,43,44,45,46,47,48,49,50,51,52,53};
 beam[113]={42,43,44,45,46,47,48,49,50,51,52,53,54};
 beam[114]={44,45,46,47,48,49,50,51,52,53,54};
 beam[115]={44,45,46,47,48,49,50,51,52,53,54};
 beam[116]={44,45,46,47,48,49,50,51,52,53,54,55};
 beam[117]={46,47,48,49,50,51,52,53,54,55};
 beam[118]={46,47,48,49,50,51,52,53,54,55,56};
 beam[119]={48,49,50,51,52,53,54,55,56};
 beam[120]={48,49,50,51,52,53,54,55,56,57};
 beam[121]={48,49,50,51,52,53,54,55,56,57};
 beam[122]={49,50,51,52,53,54,55,56,57};
 beam[123]={49,50,51,52,53,54,55,56,57,58};
 beam[124]={51,52,53,54,55,56,57,58};
 beam[125]={51,52,53,54,55,56,57,58,59};
 beam[126]={53,54,55,56,57,58,59};
 beam[127]={53,54,55,56,57,58,59,60};
 beam[128]={53,54,55,56,57,58,59,60};
 beam[129]={55,56,57,58,59,60};
 beam[130]={55,56,57,58,59,60,61};
 beam[131]={57,58,59,60,61};
 beam[132]={57,58,59,60,61,62};
 beam[133]={59,60,61,62};
 beam[134]={59,60,61,62,63,64};



 //
 //////////////////////////////////////////////////
 // Create File & Histograms                     //
 //////////////////////////////////////////////////
 
 // set scale                                                                                                                                                                                                                               
 int scale_xbin = 40;
 Double_t sigmulti = 0.5;
 Double_t scale_xmax = scale_xbin*sigmulti;
 int nchannel = 4;
 // int pickup_crossing = 5;

 int histo_number = 0;//Laura
 TH2F* SNRmap[135];
 // TH2F* SNRf = new TH2F(Form("Sigma_Fluctuation"),Form("Signal to Noise Fluctuation; Channel ; Sigma"),128,-0.5,127.5, 40,0,4);
 TH1F*SNR1HA = new TH1F (Form("Sigma_Fluctuation_1d_%d-bincut_all",bin),Form("Signal to Noise Error Significant : All, %d Bin cut;#frac{S}{#sqrt{B+S}};Counts",bin),scale_xbin,0,scale_xmax);
 TH1F*SNR1HN = new TH1F (Form("Sigma_Fluctuation_1d_%d-bincut_Noise",bin),Form("Signal to Noise Error Significant : Noise, %d Bin cut;#frac{S}{#sqrt{B+S}};Counts",bin),scale_xbin,0,scale_xmax);
 TH1F*SNR1HP = new TH1F (Form("Sigma_Fluctuation_1d_%d-bincut_Pickup",bin),Form("Signal to Noise Error Significant : Pick up Channel, %d Bin cut;#frac{S}{#sqrt{B+S}};Normalized Counts (AU)",bin),scale_xbin,0,scale_xmax);
 TH1F*SNR1HH = new TH1F (Form("Sigma_Fluctuation_1d_%d-bincut_Hit-channel",bin),Form("Signal to Noise Error : Hit Channel, %d Bin cut;#frac{S}{#sqrt{B+S}};Counts",bin),scale_xbin,0,scale_xmax);
 TH1F*SNR1HHon = new TH1F (Form("Sigma_Fluctuation_1d_%d-bincut_Hit-channel-on",bin),Form("Signal to Noise Error Significant : Hit Channel on strip, %d Bin cut;#frac{S}{#sqrt{B+S}};Counts",bin),scale_xbin,0,scale_xmax);
 TH1F*SNR1HHoff= new TH1F (Form("Sigma_Fluctuation_1d_%d-bincut_Hit-channel-off",bin),Form("Signal to Noise Error  : Hit Channel off strip, %d Bin cut;#frac{S}{#sqrt{B+S}};Counts",bin),scale_xbin,0,scale_xmax);
 TH1F*SNR1HHN = new TH1F (Form("Sigma_Fluctuation_1d_%d-bincut_Hit-Neighbor-channel",bin),Form("Signal to Noise Error Ratop : Neighbors of Hit Channel, %d Bin cut;#frac{S}{#sqrt{B+S}};Counts",bin),scale_xbin,0,scale_xmax);
 THStack *SNR1S = new THStack (Form("Sigma_Fluctuation_1d_%d-bincut_All-stack",bin),Form("Signal to Noise Error  : All Stacked, %d Bin cut;#frac{S}{#sqrt{B+S}};Counts",bin));
 TH2F*SNR2HN = new TH2F (Form("Sigma_Fluctuation_1d_%d-bincut_Noise",bin),Form("Signal to Noise Error Significant : Noise, %d Bin cut;Ratio;Threshold;counts",bin),scale_xbin,0,scale_xmax,8,50,75);
 //TH2F*SNR2HP = new TH2F (Form("Sigma_Fluctuation_1d_%d-bincut_Pickup",bin),Form("Signal to Noise Error Ratio : Pick up Channel, %d Bin cut;Ratio;Threshold;counts",bin),scale_xbin,0,scale_xmax,8,50,75);
 TH1F*SNR1HNL = new TH1F (Form("Sigma_Fluctuation_1d_%d-bincut_Noise_low",bin),Form("Signal to Noise Error Significant ( Metal Crossing < %d) : Noise, %d Bin cut;#frac{S}{#sqrt{B+S}};Counts",pickup_crossing,bin),scale_xbin,0,scale_xmax);
 TH1F*SNR1HPL = new TH1F (Form("Sigma_Fluctuation_1d_%d-bincut_Pickup_low",bin),Form("Signal to Noise Error Significant : Pick up Channel ( Metal Crossing < %d), %d Bin cut;#frac{S}{#sqrt{B+S}};Normailized Counts (AU)",pickup_crossing,bin),scale_xbin,0,scale_xmax);
 TH1F*SNR1HPN = new TH1F (Form("PickUp_Noise_Ratio_%d-bincut",bin),Form("Pick up to Noise , %d Bin Cut;#frac{S}{#sqrt{B+S}};AU",bin),scale_xbin,0,scale_xmax);
 TH1F*SNR1HNH = new TH1F (Form("Sigma_Fluctuation_1d_%d-bincut_Noise_high",bin),Form("Signal to Noise Error Significant : Noise ( Metal Crossing > %d), %d Bin cut;#frac{S}{#sqrt{B+S}};Counts",pickup_crossing,bin),scale_xbin,0,scale_xmax);
 TH1F*SNR1HPH = new TH1F (Form("Sigma_Fluctuation_1d_%d-bincut_Pickup_high",bin),Form("Signal to Noise Error Significant : Pick up Channel ( Metal Crossing > %d), %d Bin cut;#frac{S}{#sqrt{B+S}};Normailized Counts (AU)",pickup_crossing,bin),scale_xbin,0,scale_xmax);


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
	   Hitmap_Noise_Cut[k] = new TH2F(Form("Hitmap_Noise_Cut_Structure_%d_Scan_%d",structure,k), Form("Hitmap (%s), Scan %d, background subtracted; Channel ; Threshold in [mV]",struc_ID[structure].c_str(),k),128,-0.5,127.5, 16,48.8,99.2);
	   SNRmap[k] = new TH2F(Form("SNRmap_Structure_%d_%d-bincut_Scan_%d",structure,bin,k), Form("Signal to sigma ratio, Scan %d, %d_bin cut; Channel ; Ratio ",k,bin),128,-0.5,127.5, 40,0,4);
	   //cout << "defining Hitmap Noise Cut" << k << endl;

	 }//end of loop for defining Hitmap_Noise_Cut

	 //std::cout << "IN BETWEEN THE LOOPS: hit_strip_off_10[0] = " << hit_strip_off_10[0] << std::endl;
	 //std::cout << "IN BETWEEN THE LOOPS: hit_strip_on_10[0]  = " << hit_strip_on_10[0] << std::endl;

	 NoiseData = TFile::Open("../data/Noise.root");
	 NoiseData->cd();
	 raw_noise_histo = (TH2F*)gDirectory->Get(Form("Noise_Histo_%d",structure));


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
	     /*
	     NoiseData = TFile::Open("../data/Noise.root");
	     NoiseData->cd();  	
	     raw_noise_histo = (TH2F*)gDirectory->Get(Form("Noise_Histo_%d",structure));
	     */

	     //std::cout << "entries in data histo for hit strip " << hit_strip << ": " << raw_data_histo->GetEntries() << std::endl;
	     //std::cout << "adding raw_data_histo, raw_noise_histo to obtain Hitmap_Noise_Cut " << m << std::endl;

	     ////////////////////////////////////////////	    
	     ///////////Noramlization( Sam)//////////////
	     ///////////////////////////////////////////

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

	     ////////////////////////////////////////////////////
	     //////////////signal Identification/////////////////
	     ////////////////////////////////////////////////////
	     for (int cht = 1;cht <129; cht++)
	       {
		 Double_t id = 0;
		 Double_t count = 0;
		 Double_t threshold = 0;


		 //////////Pick up channel Idenification/////////
		 double_t pickup_ch = 0;
		 for (int pch = 0; pch < beam_channel[m]-hit_cht[m]-nchannel/*(neighbor channel)*/; pch++)
		   {
		     // cout<<"beam position " <<m << "number of channel"<< beam_channel[m]<<endl;
		     if (cht == beam[m].at(pch)+1)
		       {
			 pickup_ch = 1;
		       }
		     else continue; 
		   }
		 /////////Hit channel Idenification////////////
		 for (int hit_channel =0 ; hit_channel<hit_cht[m];hit_channel++) //loop option for 1 or 2 hit channel
		   {
		     if (cht == beam[m].at(beam_channel[m]-1)+1/*cht number correction*/ - hit_channel/*previous hit channel*/)
		       {
			 pickup_ch = 2;
		       }
		   }
		
		 //////////////Neighbor Channel Identification///////////////
		 for (int nch = 0 ; nch < nchannel ; nch++)
		   {
		     if (cht == beam[m].at(beam_channel[m]-1)-hit_cht[m]-nch/*last neighbor*/+1/*cht number correction*/)
		       {
			 pickup_ch=3;
		       }
		     
		     else if(cht == beam[m].at(beam_channel[m]-1)+1+nch/*next neighbor+neighbor channelS*/+1/*cht number correction*/)
		       {
			 pickup_ch=3;
		       }
		   }
		 /////////////////////////////////Fill SNR map////////////////////////////////////////                                                                                                                                

                 Double_t beam_hit[15];
                 Double_t noise_error[15];
                 //cout<< "define noise error" <<endl;                                                                                                                                                                                
                 for (int nq  = 1; nq < 11; nq++)
                   {
                     beam_hit[nq] = Hitmap_Noise_Cut[m]->GetBinContent(cht,nq);
		     noise_error[nq]= sqrt(raw_data_histo -> GetBinContent(cht,nq)); //use error from original Hitmap
		     //noise_error[nq]= raw_noise_histo->GetBinError(cht,nq); //use error from Noise Map
                   
		   }
                 //Sigma Map                                                                                                                                                                                                          

                 ///////////////////////   bin cut  ////////////////////////////                                                                                                                                                      
                 for (int sig = 0 ; sig < scale_xbin +1; sig++)
		   {
		     ////////Flucutation  Map with highest sigma entry only////////////                                                                                                                                           
		     for (int q = 1; q < 9; q++)
		       {
			 Double_t id_bin = 0;
			 if (noise_error[q] != 0)
			   {
			     for (int p = 0 ; p < bin ; p++)
			       {
				 if ( sig == 0 )
				   {
				     if (beam_hit[q+p] > (0.1)*noise_error[q+p])
				       {
					 id_bin =1;
				       }
				     else
				       {
					 id_bin=0;
					 break;
				       }
				   }
				 else
				   {
				     if (beam_hit[q+p] >sigmulti* (sig)*noise_error[q+p])
				       {
					 id_bin =1;
				       }
				     else
				       {
					 id_bin=0;
					 break;
				       }
				   }
			       }
			   }

			 if (id_bin == 1)
			   {
			     id=1;
			     threshold = q;

			     if(sig==0)
			       count = 1;
			     else
			       count = sig;

			   }
		       }
		     


			  SNRmap[m]->SetBinContent(cht,sig,50+(threshold-1)*3.125);
			 //SNRmap[m]->SetBinContent(cht,sig,counts); 
		       }      


		     /////////////////////////////////////////////////////////////////////////////////////

		     ///////////////////////// Fill 1D sigma flucation Histogram///////////////////////////
		 if (PUcut == 0)
		   {
		     if (id ==1)
		       {
			 //Option : Set SNR to highest sigma                                                                                                                                                                      
			 //SNRmap[m]->SetBinContent(cht,count,1);                                                                                                                                                                            
			 SNR1HA->AddBinContent(count);
			 if (count > 0)
			   {
			     cout<< "filled signal in channel" << cht-1 <<" with scan "<< m<< " with sigma = " << count*sigmulti <<endl;
			   }
		       }
		   }
		     ////////////////////////Option : pick up channel elimination/////////////////////////
		 else if (PUcut == 1)
		   {
		     
		     if (id ==1)
		       {
			 //Pick Up Channel
			 if (pickup_ch ==1)
			   {
			     // cout<<" signal in channel "<<cht-1<< " is elimiated" <<endl;
			     //beam_channel weighting not yet fixed
			     Double_t puweight = 1;
			       //19./(beam_channel[m]);
			     SNR1HP->AddBinContent(count,puweight);
			     //cout<< "count = "<<count <<endl;
			     cout<< "pick up channel weight = "<< puweight<<endl;
			     cout<< "Bin Content = " << SNR1HP->GetBinContent(count)<<endl;
			     // SNR2HP->SetBinContent(count,threshold-1, (SNR2HP->GetBinContent(count,threshold-1)+1));
			     if (beam_channel[m]<pickup_crossing)
			       {
				 //beam_channel weighting not yet fixed  
				 SNR1HPL->AddBinContent(count,puweight);
			       }
			     else if (beam_channel[m]>=pickup_crossing)
			       {
				 //beam_channel weighting not yet fixed  
				 SNR1HPH->AddBinContent(count,puweight);
			       }
			     if (count > 0)
                               {
				 cout<< "filled pick up channel" << cht-1 <<" with beam Position "<< m<< " with sigma = " << count*sigmulti <<endl;
			       }
			   }
			 //Hit Channel
			 else if (pickup_ch ==2)
			   {                           
			     SNR1HH->AddBinContent(count);
			     if(hit_cht[m]==1)
			       {
				 SNR1HHon->AddBinContent(count);
			       }
			     else if (hit_cht[m]==2)
			       {
				 SNR1HHoff->AddBinContent(count);
			       }
			     if (count > 0)
                               {
                                 cout<< "filled hit channel" << cht-1 <<" with beam Position "<< m<< " with sigma = " << count*sigmulti <<endl;
                               }
			   }
			 //Neighbor Channel
			 else if (pickup_ch ==3)
                           {
                             SNR1HHN->AddBinContent(count);
                             if (count > 0)
                               {
				 // cout<< "filled Neighbor channel" << cht-1 <<" with beam Position "<< m<< " with sigma = " << count*sigmulti <<endl;
                               }
                           }
			 //Noise 
			 else if (pickup_ch ==0)
			   {
			     //Option : Set SNR to highest sigma 
			     //SNRmap[m]->SetBinContent(cht,count,1);  
			     SNR1HN->AddBinContent(count);
			     SNR2HN->SetBinContent(count,threshold-1,(SNR2HN->GetBinContent(count,threshold-1)+1));
			     if (beam_channel[m]<pickup_crossing)
                               {
                                 SNR1HNL->AddBinContent(count);
                               }
                             else if (beam_channel[m]>=pickup_crossing)
			       {
                                 SNR1HNH->AddBinContent(count);
			       }

			     if (count > 0)
			       {
				 if( count == scale_xbin)
				   cout<< "filled Noise in channel" << cht-1 <<" with beam Position "<< m<< " with sigma = " << count*0.1 <<endl;
				 // cout<<"Bin content " <<"sigma = " << count*0.1 <<", Threshold "<< (50+(threshold-1)*3.125)<< " : "<<SNR2HN->GetBinContent(count,threshold-1)<<endl;
			       }
			   }
		       }
		   }
		 //////////////////////////////////////////////////////////////////////////////////////////
	       }
	     // SNRf->Add(SNRmap[m]);    
	     
	     ////////////Printing///////////
	     
	     ///Hitmap
	     //adding a TCanvas that can be printed out as PDFs (Laura)
	     /*
	     TCanvas *c = new TCanvas;
	     //Hitmap_Noise_Cut[m]->SetMinimum(0);
	     Hitmap_Noise_Cut[m]->SetMaximum(5000);
	     Hitmap_Noise_Cut[m]->Draw("COLZ");
	     c->Print(Form("../plots/Hitmap_Noise_Cut_%d.png",m));
	     delete c;


	     //Signal Map
	     TCanvas *d  =new TCanvas;
	     // SNRmap[m]->SetMinimum(50);
	     // SNRmap[m]->SetMaximum(80);
	     SNRmap[m]->GetZaxis()->SetTitle("Threshold");
	     SNRmap[m]->Draw("COLZ");
	     d->Print(Form("../plots/Sigma_%d-bincut_scan%d.png",bin,m));
	     delete d	    
	     */
	     //SNR fluctuation                                                                                                                                                                                                                    
	    RawData->Close();
	    //	    NoiseData->Close();
	 }
	//end of loop for filling Hitmap_Noise_Cut

	 SNR1HPN-> Divide(SNR1HP,SNR1HN,1,1);

	 TCanvas *g =new TCanvas;
	 raw_noise_histo -> Draw("COLZ");
	 g->Print(Form("../plots/raw_noise_histo.png"));
	 delete g;

	 NoiseData->Close();

	 //Print with Pick up option
	 if(PUcut ==0)
	   {
	     TCanvas *e = new TCanvas;
	     SNR1HA -> Draw();
	     e->Print(Form("../plots/Signal_%d-bincut_All.png",bin));
	     delete e;

	     TCanvas *f = new TCanvas;
	     gPad->SetLogy();
	     SNR1HA -> Draw();
	     f->Print(Form("../plots/Signal_%d-bincut_All_log.png",bin));
	     delete f;
	   }
	 else if (PUcut == 1)
	   {
	     SNR1HN -> SetFillColor(kYellow);
	     SNR1HNL -> SetFillColor(kYellow+2);
	     SNR1HNH -> SetFillColor(kYellow+1);
	     SNR1HP -> SetFillColor(kRed);
	     SNR1HPL -> SetFillColor(kRed+2);
	     SNR1HPH -> SetFillColor(kRed+1);
	     SNR1HH -> SetFillColor(kGreen);
	     SNR1HHon -> SetFillColor(kGreen+1);
	     SNR1HHoff -> SetFillColor(kGreen+2);
	     SNR1HHN -> SetFillColor(kBlue);
	     
	     //Noise histogram
	   TCanvas *n1 = new TCanvas;
	   SNR1HN -> Draw();
	   n1->Print(Form("../plots/Signal_%d-bincut_Noise.png",bin));
	   delete n1;

	   TCanvas *n1l = new TCanvas;
           SNR1HNL -> Draw();
           n1l->Print(Form("../plots/Signal_%d-bincut_Noise_low%d.png",bin,pickup_crossing));
           delete n1l;
	   TCanvas *n1h = new TCanvas;
           SNR1HNH -> Draw();
           n1h->Print(Form("../plots/Signal_%d-bincut_Noise_high_%d.png",bin,pickup_crossing));
           delete n1h;

	   
	   TCanvas *nlog = new TCanvas;
	   gPad->SetLogy();
	   SNR1HN -> Draw();
	   nlog->Print(Form("../plots/Signal_%d-bincut_Noise_log.png",bin));
	   delete nlog;
	  
	   TCanvas *n2 = new TCanvas;
           SNR2HN -> Draw("COLZ");
           n2->Print(Form("../plots/Signal_%d-bincut_2D_Noise.png",bin));
           delete n2;
	   /*
           TCanvas *n2log = new TCanvas;
           gPad->SetLogy();
           SNR2HN -> Draw("COLZ");
           n2log->Print(Form("../plots/Signal_%d-bincut_2D_Noise_log.png",bin));
           delete n2log;
	   */


	   //Pick up channel histogram
	   TCanvas *p1 = new TCanvas;
           SNR1HP -> Draw();
           p1->Print(Form("../plots/Signal_%d-bincut_Pickup.png",bin));
           delete p1;

	   TCanvas *p1l = new TCanvas;
           SNR1HPL -> Draw();
           p1l->Print(Form("../plots/Signal_%d-bincut_Pickup_low_%d.png",bin,pickup_crossing));
           delete p1l;
	   TCanvas *p1h = new TCanvas;
           SNR1HPH -> Draw();
           p1h->Print(Form("../plots/Signal_%d-bincut_Pickup_high_%d.png",bin,pickup_crossing));
           delete p1h;
	 
           TCanvas *plog = new TCanvas;
           gPad->SetLogy();
           SNR1HP -> Draw();
           plog->Print(Form("../plots/Signal_%d-bincut_Pickup_log.png",bin));
           delete plog;
	   /*
	   TCanvas *p2 = new TCanvas;
           SNR2HP -> Draw("COLZ");
           p2->Print(Form("../plots/Signal_%d-bincut_2D_Pickup.png",bin));
           delete p2;
	   
           TCanvas *p2log = new TCanvas;
           gPad->SetLogy();
           SNR2HP -> Draw("COLZ");
           p2log->Print(Form("../plots/Signal_%d-bincut_2D_Pickup_log.png",bin));
           delete p2log;
	   */

	   //Hit channel histogram
	   TCanvas *h1 = new TCanvas;
           SNR1HH -> Draw();
           h1->Print(Form("../plots/Signal_%d-bincut_Hit.png",bin));
           delete h1;

	   TCanvas *h1on = new TCanvas;
           SNR1HHon -> Draw();
           h1on->Print(Form("../plots/Signal_%d-bincut_Hit_on.png",bin));
           delete h1on;
	   TCanvas *h1off = new TCanvas;
           SNR1HHoff -> Draw();
           h1off->Print(Form("../plots/Signal_%d-bincut_Hit_off.png",bin));
           delete h1off;

	   /*
           TCanvas *hlog = new TCanvas;
           gPad->SetLogy();
           SNR1HH -> Draw();
           hlog->Print(Form("../plots/Signal_%d-bincut_Hit_log.png",bin));
           delete hlog;
	   */
	   //Neighbor channel histogram
	   TCanvas *hn1 = new TCanvas;
           SNR1HHN -> Draw();
           hn1->Print(Form("../plots/Signal_%d-bincut_Neighbor.png",bin));
           delete hn1;
	   /*
           TCanvas *hnlog = new TCanvas;
           gPad->SetLogy();
           SNR1HHN -> Draw();
           hnlog->Print(Form("../plots/Signal_%d-bincut_Neighbor_log.png",bin));
           delete hnlog;
	   */

	   TCanvas *hpn1 = new TCanvas;
	   SNR1HPN -> Draw();
	   hpn1 -> Print(Form("../plots/PicktoNoiseRatio_%d-bincut.png",bin));
	   delete hpn1;


	   //Stack histogram
	   SNR1S -> Add(SNR1HN);
	   SNR1S -> Add(SNR1HP);
	   SNR1S -> Add(SNR1HH);
	   SNR1S -> Add(SNR1HHN);
	   
	   TCanvas *s = new TCanvas;
	   SNR1S -> Draw();
	   TLegend *legs = new TLegend(0.68,0.58,0.85,0.86);
	   legs -> AddEntry(SNR1HN, "Noise","f");
	   legs -> AddEntry(SNR1HP, "Pick up","f");
	   legs -> AddEntry(SNR1HH, "Hit Channel","f");
	   legs -> AddEntry(SNR1HHN, "Neighbor Channel","f");
	   legs -> Draw("same");
	   s ->Print(Form("../plots/Signal_%d-bincut_Stack.png",bin));
	   delete s;

           TCanvas *slog = new TCanvas;
	   gPad->SetLogy();
	   SNR1S -> Draw();
	   TLegend *legslog = new TLegend(0.68,0.58,0.85,0.86);
	   legslog -> AddEntry(SNR1HN, "Noise","f");
           legslog -> AddEntry(SNR1HP, "Pick up (Normalized Count with AU) ","f");
           legslog -> AddEntry(SNR1HH, "Hit Channel","f");
           legslog -> AddEntry(SNR1HHN, "Neighbor Channel","f");
           legslog -> Draw("same");
           slog ->Print(Form("../plots/Signal_%d-bincut_Stack_log.png",bin));
           delete slog;

	   }
	
	File->Write();
	File->Close();
	
}

#endif /* NOISE_SUBTRACT_finestrip_C */
