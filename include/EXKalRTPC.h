#ifndef EXKALRTPC_H
#define EXKALRTPC_H 1

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

#include "EXHYBTrack.h"        
#include "EXKalDetector.h"
#include "EXEventGen.h"
#include "EXHit.h"

/////////////////////////////////////////////////////////////////
//Maximum Number of Hits in a track, has been defined in "EXEventGen.h"
#ifndef MAX_HITS_PER_TRACK
#define MAX_HITS_PER_TRACK 200
#define MIN_HITS_PER_TRACK 5
#endif
/////////////////////////////////////////////////////////////////

static const double kMpr = 0.938272013;
static const double kPi = atan(1.0)*4;

class EXKalRTPC {

 public:
  EXKalRTPC();
  virtual ~EXKalRTPC();

  //create a track, store everything into fKalHits
  //prepare a track from xyz array, make sure in time increasing order    
  bool PrepareATrack_mm(double *x_mm, double *y_mm, double *z_mm, int npt,
                        bool smearing=false, bool bIncludeCurveBackHits=true);
  bool PrepareATrack(double *x, double *y,double *z, int npt, bool smearing=false,
                     bool bIncludeCurveBackHits=true);
  bool PrepareATrack(int job, double pt_min, double pt_max, double costh_min,
                     double costh_max, double z_min=0.0, double z_max=0.0, 
                     bool bIncludeCurveBackHits=true);

  //Provide suggestion if need to apply 2nd iteration kalman filter
  //if bRemoveBackwardHits==true,  it will remove backward hits, otherwise just 
  //copy all hits pointer into fKalHits_Forward
  //It will also fill the smeared hit position array
  //Note that the hit buffer must be sorted by time in increasing order
  bool JudgeFor2ndIteration(bool bRemoveBackwardHits=true);

  //this routine only apply 1 iteration KF fit
  int  DoFitAndFilter(double *x_cm, double *y_cm, double *z_cm, int n, 
                      bool bIncludeCurveBackHits=false);
                      
  //this routine allow to apply 2 iteration KF fit
  int  DoFitAndFilter(bool bApply2Iter=false);

  //Use the last site to swim back to the beam line
  void ReconVertex(TVKalState &state, double &p, double &pt, double &pz, 
                   double &th, double &ph, double &x, double &y, double &z, 
                   double &r_rec, double &a_rec, double &b_rec);

  void SetCovMElement(double val) {fCovMElement=val;};
  
  void  Reset();

  //Do global helix fit and apply my corrections
  //return chi2
  double DoGlobalHelixFit(double *x, double *y,double *z, int npt,bool bIncludeCurveBackHits=true); 

  void Example(int job, int nevents, double pt_min, double pt_max, double costh_min, 
               double costh_max, double z_min, double z_max);


  void DrawRawHits(Int_t color, const Char_t *opt=""); //Draw raw hits
  
 private:

  //Get vertex by finding TCylinder crossing point
  int GetVextex(THelicalTrack &hel, Double_t x_bpm, Double_t y_bpm, 
                TVector3 &xx,  Double_t &dfi, double &r_rec, double &a_rec, double &b_rec);

  //get vertex by finding dca to bpm point
  int GetVextex2(THelicalTrack &hel, Double_t x_bpm, Double_t y_bpm, 
                 TVector3 &xx,  Double_t &dfi, double &r_rec, double &a_rec, double &b_rec);

  //Create a helix from 3 points to get initial parameter for Kalman Filter
  //IterDirection=true is farward, otherwise backward
  THelicalTrack GetIniHelixBy3Pts(bool IterDirection=true);

  //Do global helix fit to get initial parameter for Kalman Filter
  //IterDirection=true is farward, otherwise backward
  THelicalTrack GetIniHelixByGHF(bool IterDirection=false);
  //Apply linear regression to "Rho*dPhi vs dZ" to determine theta and z of a helix
  void CorrHelixThetaZ(int npt,double szPos[][3], double Rho, double A, double B,
                       double& Theta0, double& Z0);

  void FitForward4InitHelix(THelicalTrack &Hel_last,TKalMatrix &C_last);
  void FitBackward4InitHelix(THelicalTrack &Hel_1st,TKalMatrix &C_1st);

 public:

  TFile* fFile;

  EXHYBTrack    *fKalTrack;   // The buffer to hold the fitted result
  TObjArray     *fKalHits;    // hit buffer to hold original hits, include the backward hits
  TObjArray     *fKalHits_Forward;    // hit buffer to hold only the forward part of hits
  TKalDetCradle *fCradle;     // detctor system
  EXKalDetector *fDetector;   // detector
  EXEventGen    *fEventGen;   // enevt generator


 private:
  double fCovMElement;

 public:
  //root variables, will be stored into root tree by manager

  //store kalman filter resonstriction result
  double P_rec, Pt_rec, Pz_rec, Theta_rec, Phi_rec, X_rec, Y_rec, Z_rec;
  double R_rec, A_rec, B_rec;
  int    NDF;
  double Chi2;

  //store the initial values providing to KF
  double rho_kal_ini, tnl_kal_ini, phi0_kal_ini;


  //store reconstructed hits or detector smeared thrown hits, in unit of cm
  //will be filled in RemoveBackwardHits(), plan to use these array to fill tree 
  int HitNum;  
  double StepX_rec[MAX_HITS_PER_TRACK],StepY_rec[MAX_HITS_PER_TRACK],StepZ_rec[MAX_HITS_PER_TRACK];
  double StepPhi_rec[MAX_HITS_PER_TRACK],StepS_rec[MAX_HITS_PER_TRACK];

  //step_status is used tell if this site has been used by KF, will be updated by DoFitAndFilter()
  int step_status[MAX_HITS_PER_TRACK];

  //store 3-point helix
  double P_3pt,Pt_3pt,Theta_3pt,R_3pt,A_3pt,B_3pt;
  //store raw_global helix result before Jixie's corrections
  double Phi_hel_raw,Theta_hel_raw,R_hel_raw,A_hel_raw,B_hel_raw,Z_hel_raw;
  double X_hel_raw, Y_hel_raw, DCA_hel_raw, Chi2_hel_raw;
  //store final Global Helix Fit result
  double P_hel, Phi_hel,Theta_hel,R_hel,A_hel,B_hel,Z_hel;
  double X_hel, Y_hel, DCA_hel, Chi2_hel;

  ClassDef(EXKalRTPC,1)   // KalRTPC kernel module
};

#endif
