#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
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
#include "map"
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include "irradiated_test_beam_analysis_03.C"
#include "Noise.C"
#include "NoiseSubtract.C"
#include "TROOT.h"
#include "TRint.h"

using namespace std;

void runAnalysis()
{
  
  gROOT->ProcessLine(".L irradiated_test_beam_analysis_03.C++");
  gROOT->ProcessLine(".L Noise.C++");
  gROOT->ProcessLine(".L NoiseSubtract.C++");
  
  cout<<" Structure 3, OFF "<<endl;
  irradiated_test_beam_analysis_03(3,0);
  cout<<" Structure 3, ON "<<endl;
  irradiated_test_beam_analysis_03(3,1);
  
  cout<<" Structure 8, OFF "<<endl;
  irradiated_test_beam_analysis_03(8,0);
  cout<<" Structure 8, ON "<<endl;
  irradiated_test_beam_analysis_03(8,1);
  
  cout<<" Structure 10, OFF "<<endl;
  irradiated_test_beam_analysis_03(10,0);
  cout<<" Structure 10, ON "<<endl;
  irradiated_test_beam_analysis_03(10,1);
  
  cout<<" Structure 12, OFF "<<endl;
  irradiated_test_beam_analysis_03(12,0);
  cout<<" Structure 12, ON "<<endl;
  irradiated_test_beam_analysis_03(12,1);
  
  cout<<"NOISE"<<endl;
  Noise();
  
  cout<<"SUBTRACT BACKGROUND"<<endl;
  NoiseSubtract(3,1);
  NoiseSubtract(3,0);
  NoiseSubtract(8,1);
  NoiseSubtract(8,0);
  NoiseSubtract(10,1);
  NoiseSubtract(10,0);
  NoiseSubtract(12,1);
  NoiseSubtract(12,0);
  
}
