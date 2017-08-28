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
#include "Irradiated_XTalk.C"
#include "TROOT.h"
#include "TRint.h"

using namespace std;

void runXTalk()
     {

	gROOT->ProcessLine(".L Irradiated_XTalk.C");
	
	cout<<" Structure 3 "<<endl;
	Irradiated_XTalk(3);
	
	cout<<" Structure 8 "<<endl;
	Irradiated_XTalk(8);

  	cout<<" Structure 10 "<<endl;
  	Irradiated_XTalk(10);
  	
	cout<<" Structure 12 "<<endl;
	Irradiated_XTalk(12);

     }
