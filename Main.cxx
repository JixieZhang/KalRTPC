//By Jixie Zhang: RTPC Kalman filter main program
//This program is used to test the kalman filter module
//
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <iostream>

#include <TROOT.h>
#include <TApplication.h>

#include "EXKalManager.h"

using namespace std;

//new usage:  exe [-h] [-i] [-j <jobtype=0>] [-n <nevents=10000>]  [-e <error=0.05> ] 
//                [-c <space_cm=1.1> <max_ang_deg=30> <min_ang_deg=23.3> <ang_sep_cm=0.4>] 
//                [-p <pt_min_gev=0.1> <pt_max_gev=0.1>]  
//                [-t <costh_min=0> <costh_max=0>] 
//                [-z <z_min_cm=0> <z_max_cm=0>]   
//                [-i <infile="infile.root">]   
//also compatible to old usage:
//  exe <job=0|1|2> <nevent> [pt_min_gev=0.1] [pt_max_gev=0.1] 
//      [costh_min=-0.00001] [costh_max=0.00001] [error=0.05] [z_min=0.0] [z_max=0.0] 
//      [infile=infile.root] [g4type=0]

void usage()
{
  printf("/------------------------------usage start----------------------------------------\\ \n");
  printf("exe [-h] [-i] [-j <jobtype=0>] [-n <nevents=10000>]  [-e <error=0.05> ] \\\n");
  printf("    [-c <ntracks=1> [space_cm=1.9] [min_ang_deg=30] [max_ang_deg=40.0] [ang_sep_cm=0.4>]] \\\n");
  printf("    [-p <pt_min_gev=0.1> <pt_max_gev=0.1>]  \\\n");
  printf("    [-t <costh_min=-0.00001> <costh_max=0.00001>] \\\n");
  printf("    [-z <z_min_cm=0.0> <z_max_cm=0.0>] \\\n");  
  printf("    [-f <infile=\"infile.root\">] \\\n");   
  printf("    [-g <g4type=0> \n\n");  

  cerr << "-------------------The following old usage is also supported----------------------"<<endl;
  cerr <<"exe <job=0|1|2|3|4|5|6> <nevent> [pt_min_gev=0.1] [pt_max_gev=0.1] \\"<<endl 
    <<"   [costh_min=-0.00001] [costh_max=0.00001] [error=0.05] [z_min=0.0] [z_max=0.0] \\"<<endl
    <<"   [infile=infile.root] [g4type=0]"<<endl<<endl;

  cerr << "\t  -h show this help menu;\n";
  cerr << "\t  -i interactive mode;\n";
  cerr << "\t  job: 0: run KalmanFilter with generated helix;\n";
  cerr <<" \t       1: run KalmanFilter with geant4 track;\n";
  cerr << "\t       2: run KalmanFilter with generated circle;\n";
  cerr << "\t       3: run ChainFinder but no fitting with geant4 tracks, need to couple with -c option; \n";
  cerr << "\t       4: run ChainFinder + Global Helix Fitter with geant4 tracks\n";
  cerr << "\t       5: run ChainFinder + KalmanFilter with geant4 tracks\n";
  cerr << "\t  nevents: number of events to generate \n";
  cerr << "\t  pt_min_gev and pt_max_gev: specifiy the range of pt in Gev \n";
  cerr << "\t  Note that if pt is negative then anti-clockwise track will be generated \n";
  cerr << "\t  costh_min and costh_max: specifiy the range of costh, only for job==0\n";
  cerr << "\t  error is used to initialize the comvariant matrix before fitting. \n";
  cerr << "\t  z_min and z_max are used to tell the range of thrown z. \n";
  cerr << "\t  infile is used to specifiy the input root file. \n";
  cerr << "\t  g4type tells the type of the input file: 0, Jixie's Geant4 RTPC12 tree, 1 GEMC tree \n\n";
  cerr << "example1: exe 0 1000 ' -0.05' ' -0.07' ' -0.8' 0.8 0.05 \n";
  cerr << "example2: exe 0 1000 ' -0.05' ' -0.07' ' -0.8' 0.8 0.05 ' -20.0' 20.0 infile.root 1\n";
  cerr << "example3: exe -j 0 -n 1000 -p -0.05 -0.07 -c 25 1.9 30.0 40.0 0.4 -z -20 20 -i \n";
  printf("\\------------------------------usage  end------------------------------------------/ \n");
  abort();
}

int main(int argc, char **argv)
{
  const double rad2deg = 180./(4.*atan(1.));

  //by default run in batch mode, will switch into interactive mode by -i option
  gROOT->SetBatch(true);

  int    job=0;
  int    nevents=10000;
  double error=0.05;
  int    ntracks=1;  //how many track to dump into hitpool in an event
  double space=1.9, min_ang=30.0/rad2deg, max_ang=40.0/rad2deg, ang_sep=0.4;
  double pt_min=0.1, pt_max=0.1;
  double costh_min=-0.00001, costh_max=0.00001;
  double z_min=0.0, z_max=0.0;
  char   infile[100]; sprintf(infile, "%s","infile.root");
  int    g4type=0;

  if(argc<3) usage();

  int index, c;
  int optmarker = 0;
  while ((c = getopt (argc, argv, "hij:n:e:c:p:t:z:f:x:g:")) != -1) {
    switch (c) {
    case 'h':
      usage();
      break;      
    case 'i':
      gROOT->SetBatch(false);
      break;
    case 'j':
      job = atoi(optarg);
      break;
    case 'n':
      nevents = atoi(optarg);
      break;
    case 'e':
      error = atof(optarg);
      break;
    case 'x':
      optind--;
      for(int i=0;i<4;i++){
        if ( optind < argc && (argv[optind][0]!='-' || (argv[optind][0]=='-' && isdigit(argv[optind][1]))) ) {
          //fprintf(stderr, "\n-x option: argument_%d ='%s'\n",i,argv[optind]);
          if(i==0) space=atof(argv[optind]);
          if(i==1) max_ang=atof(argv[optind])/rad2deg;
          if(i==2) min_ang=atof(argv[optind])/rad2deg;
          if(i==3) ang_sep=atof(argv[optind]);
          optind++;
        } else {  
          fprintf(stderr, "\n-x option require 4 arguments \n\n");
          usage();
        }
      }
      break;
    case 'c':
      //this is another way to do the same thing as -c
      ntracks = atoi(argv[optind-1]);
      optmarker=optind-1;
      for( ;optind < argc && (argv[optind][0]!='-' || (argv[optind][0]=='-' && isdigit(argv[optind][1]))); optind++) {
        //printf("argu[0]=%c  argu[1]=%c\n",argv[optind][0],argv[optind][1]);
        //fprintf(stderr, "\n-c option argument_%d='%s'\n",optind-optmarker,argv[optind]);
        if(optind-optmarker==1) space=atof(argv[optind]);
        if(optind-optmarker==2) min_ang=atof(argv[optind])/rad2deg;
        if(optind-optmarker==3) max_ang=atof(argv[optind])/rad2deg;
        if(optind-optmarker==4) ang_sep=atof(argv[optind]);
        if(optind-optmarker==4) {optind++; break;}
      }
      break;
    case 'p':
      pt_min = atof(argv[optind-1]);
      if ( optind < argc && (argv[optind][0]!='-' || (argv[optind][0]=='-' && isdigit(argv[optind][1]))) ) {
        pt_max = atof(argv[optind]);
        optind++;
      } else {
        fprintf(stderr, "\n-p option require 2 arguments \n\n");
        usage();
      }
      break;
    case 't':
      costh_min = atof(argv[optind-1]);
      if ( optind < argc && (argv[optind][0]!='-' || (argv[optind][0]=='-' && isdigit(argv[optind][1]))) ) {
        costh_max = atof(argv[optind]);
        optind++;
      } else {
        fprintf(stderr, "\n-t option require 2 arguments \n\n");
        usage();
      }
      break;
    case 'z':
      z_min = atof(argv[optind-1]);
      if ( optind < argc && (argv[optind][0]!='-' || (argv[optind][0]=='-' && isdigit(argv[optind][1]))) ) {
        z_max = atof(argv[optind]);
        optind++;
      } else {
        fprintf(stderr, "\n-z option require 2 arguments \n\n");
        usage();
      }
      break;
    case 'f':
      sprintf(infile, "%s",argv[optind-1]);
      break;
    case 'g':
      g4type = atoi(optarg);
      break;
    case '?':
      if (strchr("j:n:e:c:p:t:z:f:x:g:",optopt))
        fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint (optopt))
        fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      else
        fprintf (stderr,
        "Unknown option character `\\x%x'.\n",
        optopt);
      return -1;
    default:
      abort();
    }
  }

  //taking care of those no option arguments
  //backward compatible
  index = 0;  
  for (; optind < argc; optind++) {
    index++;
    //printf ("Non-option argument argv[%d]=%s\n", index, argv[optind]);
    if(index==1) job = atoi(argv[optind]);
    if(index==2) nevents = atoi(argv[optind]);
    if(index==3) pt_min = atof(argv[optind]);
    if(index==4) pt_max = atof(argv[optind]);
    if(index==5) costh_min = atof(argv[optind]);
    if(index==6) costh_max = atof(argv[optind]);
    if(index==7) error = atof(argv[optind]);
    if(index==8) z_min = atof(argv[optind]);
    if(index==9) z_max = atof(argv[optind]);
    if(index==10)  sprintf(infile, "%s",argv[optind]);
    if(index==11) g4type = atoi(argv[optind]);
  }

  cout<<"arguments: job="<<job<<", nevents="<<nevents<<", error="<<error<<endl 
    <<"\t pt_min="<<pt_min<<", pt_max="<<pt_max<<", costh_min="<<costh_min<<", costh_max="<<costh_max<<endl
    <<"\t z_min="<<z_min<<", z_max="<<z_max<<", infile=\'"<<infile<<"\', g4type="<<g4type<<endl
    <<"\t ntracks="<<ntracks<<", space="<<space<<" min_ang="<<min_ang*rad2deg<<", max_ang="<<max_ang*rad2deg<<", ang_sep="<<ang_sep<<endl
    <<endl;


  ////////////////////////////////////////////////////////////////////////////////
  EXKalManager aKalMan;
  aKalMan.SetCovMElement(error);
  aKalMan.SetG4InputFile(infile);
  if(job>=3 && job<=5 ) {
    aKalMan.RunCFNFit(g4type, job, nevents,ntracks,space,min_ang,max_ang,ang_sep);
  } 
  else{
    aKalMan.RunKF(job,nevents,pt_min,pt_max,costh_min,costh_max,z_min,z_max);
  }

  return 0;
}
