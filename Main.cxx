//By Jixie Zhang: RTPC Kalman filter main program
//This program is used to test the kalman filter module
//
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <TROOT.h>
#include <TApplication.h>
using namespace std;

#include "EXKalManager.h"


int main (int argc, char **argv)
{
  gROOT->SetBatch();
  TApplication app("EXKalRTPC", &argc, argv, 0, 0);

  int job;
  int nevents;
  double pt_min=0.1, pt_max=0.1;
  double costh_min=-0.00001, costh_max=0.00001;
  double z_min=0.0, z_max=0.0;
  double error=0.05;
  char   infile[100];

  if(argc<3) {
    cerr << "Usage: "<<argv[0] <<" <job=0|1|2> <nevent> [pt_min_gev=0.1] [pt_max_gev=0.1] \\"<<endl 
         <<"        [costh_min=-0.00001] [costh_max=0.00001] [error=0.05] [z_min=0.0] [z_max=0.0] \\"<<endl
	 <<"        [infile=infile.root]"<<endl<<endl;
    cerr << "\t  job: 0 generate helix, 1 loadtrack from geant4 root file, 2 generate circle\n";
    cerr << "\t  nevents: number of events to generate \n";
    cerr << "\t  pt_min_gev and pt_max_gev: specifiy the range of pt in Gev \n";
    cerr << "\t  Note that if pt is negative then anti-clockwise track will be generated \n";
    cerr << "\t  costh_min and costh_max: specifiy the range of costh, only for job==0\n";
    cerr << "\t  error is used to initialize the comvariant matrix before fitting. \n";
    abort();
  }

  job = atoi(argv[1]);
  nevents = atoi(argv[2]);
  if(argc>3) pt_min = atof(argv[3]);
  if(argc>4) pt_max = atof(argv[4]);
  if(argc>5) costh_min = atof(argv[5]);
  if(argc>6) costh_max = atof(argv[6]);
  if(argc>7) error = atof(argv[7]);
  if(argc>8) z_min = atof(argv[8]);
  if(argc>9) z_max = atof(argv[9]);
  if(argc>10)  sprintf(infile, "%s",argv[10]);
  else  sprintf(infile, "%s","infile.root");

  EXKalManager aKalFilter;
  aKalFilter.SetCovMElement(error);
  aKalFilter.SetG4InputFile(infile);
  aKalFilter.Run(job,nevents,pt_min,pt_max,costh_min,costh_max,z_min,z_max);

  return 0;
}
