#ifndef EXEVENTGEN_H
#define EXEVENTGEN_H

#include "TKalDetCradle.h"
#include "THelicalTrack.h"

/////////////////////////////////////////////////////////////////
//-----------------------------------
// RTPC Parameters
//-----------------------------------
static const double kRTPC_R_GEM1 = 7.0;
static const double kRTPC_R_Cathode = 3.0;

//Maximum Number of Hit in a track 
#define MaxHit 200
#define MinHit 5
//-----------------------------------
// Track Parameters
//-----------------------------------

/////////////////////////////////////////////////////////////////

class EXEventGen {
 public:
 EXEventGen(TKalDetCradle &cradle, TObjArray &kalhits)
   : fCradlePtr(&cradle), fHitBufPtr(&kalhits) {}
  virtual ~EXEventGen() {}

  THelicalTrack GenerateHelix(double pt_min, double pt_max, double cosmin=0.0, 
			      double cosmax=0.0, double z_min=0.0, double z_max=0.0);

  void Swim(THelicalTrack &heltrk, Bool_t bIncludeCurveBackHits, double mass);

  int  GenerateCircle(double pt_min, double pt_max, double cosmin=0.0, double cosmax=0.0, 
		      double z_min=0.0, double z_max=0.0, bool bIncludeCurveBackHits=true);

  //input: x y z in cm and in increasing time order 
  void MakeHitsFromTraj(double *x_cm, double *y_cm, double *z_cm, int npt, bool smearing=false,
			bool bIncludeCurveBackHits=true);
  void MakeHitsFromTraj_mm(double *x_mm, double *y_mm, double *z_mm, int npt, bool smearing=false,
			   bool bIncludeCurveBackHits=true);

  static void     SetT0(Double_t t0) { fgT0 = t0;   }
  static Double_t GetT0()            { return fgT0; }

  static void PrintHelix(THelicalTrack *aTrack, const char *title="helix");

 private:
  TKalDetCradle *fCradlePtr;     // pointer to detector system
  TObjArray     *fHitBufPtr;     // pointer to hit array

  static Double_t  fgT0;         // t0

 public:

  //some buff to hold some thrown variables, they will be used to fill root tree
  double X0,Y0,Z0;               //at vertex
  double P0,Theta0,Phi0;         //at vertex  
  
  //store the oringinal hit positions before smearing by detector resolution
  int StepNum;
  double StepX[MaxHit],StepY[MaxHit],StepZ[MaxHit],StepPhi[MaxHit],StepS[MaxHit];

  //record the generated helix at the 1st and last hit for studying Kalman filter
  //note that these values are not available for circle or geant4 track 
  double Rho_1st, TanLambda_1st, Phi0_1st;
  double Rho_last, TanLambda_last, Phi0_last;

  ClassDef(EXEventGen,1)   // Event Generator
    
};

#endif
