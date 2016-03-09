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
   
   //DO a global helix fit to get initial parameter for Kalman Filter
   THelicalTrack DoHelixFit();

   void          Swim(THelicalTrack &heltrk, Double_t mass = 0.13957018);

   int           GenCircle(double pt_min, double pt_max);

   int           LoadOneTrack();
   //input: x y z in mm and in increasing order 
   void          MakeHitsFromTraj(double *x, double *y, double *z, int npt,
				  bool smearing=false);

   static void     SetT0(Double_t t0) { fgT0 = t0;   }
   static Double_t GetT0()            { return fgT0; }

private:
   TKalDetCradle *fCradlePtr;     // pointer to detector system
   TObjArray     *fHitBufPtr;     // pointer to hit array

   static Double_t  fgT0;         // t0

   ClassDef(EXEventGen,1)   // Event Generator
};

#endif
