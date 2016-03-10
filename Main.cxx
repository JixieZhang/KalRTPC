//RTPC Kalman filter main program
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <TROOT.h>
#include <TApplication.h>
using namespace std;

#include "EXKalRTPC.h"
 
int main (int argc, char **argv)
{
  gROOT->SetBatch();
  TApplication app("EXKalRTPC", &argc, argv, 0, 0);

  int job;
  int nevents;
  double pt_min=0.1, pt_max=0.1;
  double costh_min=-0.00001, costh_max=0.00001;

  if(argc<3) {
    cerr << "Usage: "<<argv[0] <<" <job=0|1|2> <nevent> [pt_min_gev=0.1] [pt_max_gev=0.1]" 
      <<"[costh_min=-0.00001] [costh_max=0.00001]"<< endl;
    cerr << "\t  job: 0 generate helix, 1 loadtrack from geant4 root file, 2 generate circle\n";
    cerr << "\t  nevents: number of events to generate \n";
    cerr << "\t  pt_min_gev and pt_max_gev: specifiy the range of pt in Gev \n";
    cerr << "\t  Note that if pt is negative then anti-clockwise track will be generated \n";
    cerr << "\t  costh_min and costh_max: specifiy the range of costh, only for job==0\n";
    abort();
  }

  job = atoi(argv[1]);
  nevents = atoi(argv[2]);
  if(argc>3) pt_min = atof(argv[3]);
  if(argc>4) pt_max = atof(argv[4]);
  if(argc>5) costh_min = atof(argv[5]);
  if(argc>6) costh_max = atof(argv[6]);
  
  
  EXKalRTPC aKalFilter;
  aKalFilter.KalRTPC(job,nevents,pt_min,pt_max,costh_min,costh_max);

  return 0;
}
