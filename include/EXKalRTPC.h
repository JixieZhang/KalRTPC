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

#include "EXKalDetector.h"
#include "EXEventGen.h"
#include "EXHit.h"

#define MaxHit 200

class EXKalRTPC {

public:
  EXKalRTPC();
  virtual ~EXKalRTPC();

  void BeginOfRun();
  void EndOfRun();

  //create a track, store everything into fKalHits
  //prepare a track from xyz array, make sure radii are in increasing order  
  bool PrepareATrack(double *x_mm, double *y_mm,double *z_mm, int npt);
  bool PrepareATrack(int job, double pt_min, double pt_max, 
    double costh_min, double costh_max);
  
  void FitForward4InitHelix(THelicalTrack &Hel_last,TKalMatrix &C_last);
  void FitBackward4InitHelix(THelicalTrack &Hel_1st,TKalMatrix &C_1st);
  int DoFitAndFilter(double *x_mm, double *y_mm, double *z_mm, int n);

  int KalRTPC(int job, int nevents, double pt_min, double pt_max, 
    double costh_min, double costh_max);

  //reconstruct to vertex using helix at last hit
  void ReconVertex(TVKalState &state, double &p, double &pt, double &pz, 
    double &th, double &ph, double &x, double &y, double &z, 
    double &r_rec, double &a_rec, double &b_rec);
  
  void SetCovMElement(double val) {fCovMElement=val;};

private:

  void  Tree_Init();
  void  Tree_Fill(TKalTrack &kaltrack);
  void  Reset();

  //Get vertex by finding TCylinder crossing point
  int GetVextex(THelicalTrack &hel, Double_t x_bpm, Double_t y_bpm, 
    TVector3 &xx,  Double_t &dfi, double &r_rec, double &a_rec, double &b_rec);

  //get vertex by finding dca to bpm point
  int GetVextex2(THelicalTrack &hel, Double_t x_bpm, Double_t y_bpm, 
    TVector3 &xx,  Double_t &dfi, double &r_rec, double &a_rec, double &b_rec);
  

public:

  TFile* fFile;

  TObjArray     *fKalHits;    // hit buffer
  TKalDetCradle *fCradle;     // detctor system
  EXKalDetector *fDetector;   // detector
  EXEventGen    *fEventGen;   // enevt generator
  
private:
  double fCovMElement;

private:
  //root variables
  TTree* fTree;
  int _index_;
  double p_rec,pt_rec,pz_rec,th_rec,ph_rec,x_rec,y_rec,z_rec;
  double r_rec,a_rec,b_rec;
  int npt;
  double step_x[MaxHit],step_y[MaxHit],step_z[MaxHit];
  double step_px[MaxHit],step_py[MaxHit],step_pz[MaxHit];
  double step_bx[MaxHit],step_by[MaxHit],step_bz[MaxHit];
  int step_status[MaxHit];
  double step_x_exp[MaxHit],step_y_exp[MaxHit],step_z_exp[MaxHit];
  double step_x_fil[MaxHit],step_y_fil[MaxHit],step_z_fil[MaxHit];

  //varible from global fitter, come from g4 root tree 
  double p0,pt0,pz0,th0,ph0,_x0_,_y0_,_z0_;
  double p_hel,pt_hel,pz_hel,th_hel,ph_hel,x_hel,y_hel,z_hel;
  double r_hel,a_hel,b_hel;
  double p_3pt,pt_3pt,th_3pt,r_3pt,a_3pt,b_3pt;

  int ndf;
  double chi2,cl;

};

#endif
