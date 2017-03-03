#ifndef EXKalManager_H
#define EXKalManager_H 1

#include "TROOT.h"
#include "TApplication.h"
#include "TRint.h"
#include "TRandom.h" 
#include "TFile.h"      
#include "TTree.h"  

#include "TKalDetCradle.h"    // from KalTrackLib
#include "TKalTrackState.h"   // from KalTrackLib
#include "TKalTrackSite.h"    // from KalTrackLib
#include "TKalTrack.h"        // from KalTrackLib

#include "EXKalRTPC.h"
#include "EXKalDetector.h"
#include "EXEventGen.h"
#include "EXHit.h"

#include "ChainFinder.hh"
#include "NtReader.h"
/////////////////////////////////////////////////////////////////
//Maximum Number of Hit in a track, has been defined in "EXEventGen.h" 
#ifndef MaxHit
#define MaxHit 200
#define MinHit 5
#endif
/////////////////////////////////////////////////////////////////


class EXKalManager {

public:
  EXKalManager();
  virtual ~EXKalManager();

  void BeginOfRun();
  void EndOfRun();

  //read ntracks from G4MC_RTPC12 output tree and fill it into ChainFinder's hit pool.
  //return false if reach the end of file or fail
  //please note that the G4 root tree use unit of mm 
  bool FillCFHitPoolByG4Track(int ntracks);

  //read a track from G4MC_RTPC12 output tree and fill into fKalHits
  bool LoadAG4Track(bool bIncludeCurveBackHits=false);
 
  //run ChianFinder and global helix fitter
  int  RunCFNGHF(int nevents, int ntracks, double space, double min_ang, 
                 double max_ang, double ang_sep);
  
  //run ChianFinder + KalmanFilter
  int RunCFNKF(int nevents, int ntracks, double space, double min_ang, 
               double max_ang, double ang_sep);
  
  //run KalmanFilter only    
  int RunKF(int job, int nevents, double pt_min, double pt_max, double costh_min, 
            double costh_max, double z_min=0.0, double z_max=0.0);

  void SetCovMElement(double val) {fKalRTPC->SetCovMElement(val);};
  void SetG4InputFile(const char* val) {sprintf(fG4Inputfile, "%s","infile.root");};

private:

  void  Tree_Init();
  void  Tree_Fill(TKalTrack &kaltrack);
  void  Tree_Reset();

public:

  TFile* fFile;

  ChainFinder   *fChainFinder;// RTPC ChainFinder
  EXKalRTPC     *fKalRTPC;    // RTPC KalmanFilter system

private:

  TObjArray     *fKalHits;    // hit buffer
  EXEventGen    *fEventGen;   // enevt generator

  char          fG4Inputfile[255];
  NtReader      *fNtReader;           

private:
   
  int fChainFinderNTracks;
  //This part is used to fill the root tree
  //When load G4 track, we also need to load the following thrown parameters and store
  //them into the output tree. in unit of cm.
  //Fix me:
  //the chain finder might not put the track it found in the same order
  //as the original G4 tree, therefore I do not know how to incert these values
  double X0[MAX_CHAINS_PER_EVENT],Y0[MAX_CHAINS_PER_EVENT], Z0[MAX_CHAINS_PER_EVENT];
  double Theta0[MAX_CHAINS_PER_EVENT], Phi0[MAX_CHAINS_PER_EVENT], P0[MAX_CHAINS_PER_EVENT];
  
  //root variables
  //a lot of vaiable here are redundent, they have been defined in KalRTPC or EventGen
  //I just copy their values in Tree_fill
  //To run this program fsat, do not copy the value, just simply point
  //the tree leaves address to the original variables in  KalRTPC or EventGen

  TTree* fTree;
  int _eventtype_;  //to dientify what type of event this is, 0 helix, 1 g4track, 2 circle
  int _index_;      //store tree index 
  int eventid, trackid, ntrack;  //store real event id, trackid in current event, number of track

  //thrown info, from EXEventGen
  double p0,pt0,pz0,th0,ph0,_x0_,_y0_,_z0_;

  //KF result
  double p_rec,pt_rec,pz_rec,th_rec,ph_rec,x_rec,y_rec,z_rec;
  double r_rec,a_rec,b_rec;
  int ndf;
  double chi2,cl;

  //original hits info
  int npt0;   	   //npt0 is to store number of original hits 
  double step_x[MaxHit],step_y[MaxHit],step_z[MaxHit],step_phi[MaxHit],step_s[MaxHit];  
  
  //smeared or reconstructed hits info
  int npt;         //npt is to store number of used hits by kalman filter 
  int step_status[MaxHit];
  double step_x_rec[MaxHit],step_y_rec[MaxHit],step_z_rec[MaxHit],step_phi_rec[MaxHit],step_s_rec[MaxHit];
  double step_x_exp[MaxHit],step_y_exp[MaxHit],step_z_exp[MaxHit],step_phi_exp[MaxHit],step_s_exp[MaxHit];
  double step_x_fil[MaxHit],step_y_fil[MaxHit],step_z_fil[MaxHit],step_phi_fil[MaxHit],step_s_fil[MaxHit];

  //From EXKalRTPC
  double p_3pt,pt_3pt,th_3pt,r_3pt,a_3pt,b_3pt;

  double p_hel,pt_hel,pz_hel,th_hel,ph_hel,x_hel,y_hel,z_hel;
  double r_hel,a_hel,b_hel;
  double dca_hel,chi2_hel;

  double r_hel_raw,th_hel_raw,ph_hel_raw,a_hel_raw,b_hel_raw,z_hel_raw;

  //I want to study the Kalman Filter fitted resolution as a function of
  //initial parameter: Rho, Phi0, Theta
  //the true variable at 1st and last site
  double rho_1st, tnl_1st, phi0_1st;
  double rho_last, tnl_last, phi0_last;
  double rho_kal_ini, tnl_kal_ini, phi0_kal_ini;
};

#endif
