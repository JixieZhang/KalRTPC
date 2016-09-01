#ifndef EXEVENTGEN_H
#define EXEVENTGEN_H

#include "TKalDetCradle.h"
#include "THelicalTrack.h"
#include "NtReader.h"

class EXEventGen : public NtReader {
public:
   EXEventGen(TKalDetCradle &cradle, TObjArray &kalhits)
             : fCradlePtr(&cradle), fHitBufPtr(&kalhits) {}
   virtual ~EXEventGen() {}

   THelicalTrack GenerateHelix(double pt_min, double pt_max,
                               double cosmin, double cosmax);
   

   //Create a helix from 3 points to get initial parameter for Kalman Filter
   //IterDirection=true is farward, otherwise backward
   THelicalTrack CreateInitialHelix(bool IterDirection=true);
 
   //Apply linear regression to "Rho*dPhi vs dZ" to determine theta and z of a helix
   void FitHelixThetaZ(int npt,double szPos[][3], double Rho, double A, double B,
                       double& Theta0, double& Z0);

   //Do global helix fit to get initial parameter for Kalman Filter
   //IterDirection=true is farward, otherwise backward
   THelicalTrack DoHelixFit(bool IterDirection=false);

   void          Swim(THelicalTrack &heltrk, Double_t mass = 0.13957018);

   int           GenCircle(double pt_min, double pt_max,
			  double cosmin=0, double cosmax=0);

   int           LoadOneTrack();
   //input: x y z in mm and in increasing order 
   void          MakeHitsFromTraj(double *x, double *y, double *z, int npt,
				  bool smearing=false);

   static void     SetT0(Double_t t0) { fgT0 = t0;   }
   static Double_t GetT0()            { return fgT0; }
   
   static void PrintHelix(THelicalTrack *aTrack, const char *title="helix");

private:
   TKalDetCradle *fCradlePtr;     // pointer to detector system
   TObjArray     *fHitBufPtr;     // pointer to hit array

   static Double_t  fgT0;         // t0
   
public:
  //debug 3-point helix
  double P_3pt,Pt_3pt,Theta_3pt,R_3pt,A_3pt,B_3pt;

   ClassDef(EXEventGen,1)   // Event Generator
};

#endif
