#ifndef EXEVENTGEN_H
#define EXEVENTGEN_H

#include "TKalDetCradle.h"
#include "THelicalTrack.h"

/////////////////////////////////////////////////////////////////
//-----------------------------------
// RTPC Track Parameters
//-----------------------------------
//Maximum Number of Hits in a track, has been defined in "ChainFinder.hh" 
#ifndef MAX_HITS_PER_TRACK
#define MAX_HITS_PER_TRACK 200
#define MIN_HITS_PER_TRACK 5
#endif

/////////////////////////////////////////////////////////////////

class EXEventGen {
public:
  EXEventGen(TKalDetCradle &cradle, TObjArray &kalhits);
  virtual ~EXEventGen();

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
  
  //binary searches for a value in an decreasing sorted array
  //   arr is an array to search in, in decreasing order
  // value is searched value
  //  left is an index of left boundary
  // right is an index of right boundary
  //returns position of searched value, if it presents in the array
  //        or -1*position_it_should_be_incerted, if it is absent
  static int BinarySearch(const double* arr, double value, int left, int right);

private:
  TKalDetCradle *fCradlePtr;     // pointer to detector system
  TObjArray     *fHitBufPtr;     // pointer to hit array

  static Double_t  fgT0;         // t0
    
  //use to store the detector layer boundaries
  //fDetLayerRBoundary[0]=kRTPC_R_GEM1, 
  //fDetLayerRBoundary[kNDetLayer]=kRTPC_R_Cathode, 
  // then fDetLayerRBoundary[i=1...kNDetLayer-1] := (kDetLayerRList[i-1]+kDetLayerRList[i])/2
  double *fDetLayerRBoundary;   //[kNDetLayer+1];

public:

  //some buff to hold some thrown variables, they will be used to fill root tree
  double X0,Y0,Z0;               //at vertex
  double P0,Theta0,Phi0;         //at vertex  
  
  //store the oringinal hit positions before smearing by detector resolution
  int StepNum;
  double StepX[MAX_HITS_PER_TRACK],StepY[MAX_HITS_PER_TRACK],StepZ[MAX_HITS_PER_TRACK];
  double StepPhi[MAX_HITS_PER_TRACK],StepS[MAX_HITS_PER_TRACK];

  //record the generated helix at the 1st and last hit for studying Kalman filter
  //note that these values are not available for circle or geant4 track 
  double Rho_1st, TanLambda_1st, Phi0_1st;
  double Rho_last, TanLambda_last, Phi0_last;

  ClassDef(EXEventGen,1)   // Event Generator
    
};

#endif
