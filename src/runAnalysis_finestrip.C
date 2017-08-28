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
#include "irradiated_finestrip_analysis.C"
#include "Noise_finestrip.C"
#include "NoiseSubtract_finestrip.C"
#include "TROOT.h"
#include "TRint.h"

using namespace std;

void runAnalysis_finestrip(int ana, int noi, int pick)
{
  
  gROOT->ProcessLine(".L irradiated_finestrip_analysis.C++");
 gROOT->ProcessLine(".L Noise_finestrip.C++");
  gROOT -> ProcessLine(".L NoiseSubtract_finestrip.C++");

   if(ana == 1)
    {
      cout<<" Structure 3, ON "<<endl;
      irradiated_finestrip_analysis(3,1);
    }

  if(noi==1)
    {
      cout<<"NOISE"<<endl;
      Noise_finestrip();
    }

 
 if(pick==1)
    {
      cout<<"SUBTRACT BACKGROUND"<<endl;
      for ( int i = 1; i<6;i++)
	{  
	  for ( int j =5; j<31; j++)
	    {
	      NoiseSubtract_finestrip(3,1,i,1,j) ;
	    }
	}
    }
}
