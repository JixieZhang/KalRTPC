#ifndef EXDETECTOR_H
#define EXDETECTOR_H

#include "TVector3.h"         // from ROOT
#include "EXVKalDetector.h"   // from KalTrackLib
#include "EXMeasLayer.h"

extern   const int    kNDetDummyLayer ;  // = 6;
extern   const int    kNDetLayer;        // = 35;
extern   const double kDetLayerRList[];  //in decreasing order

class EXKalDetector : public EXVKalDetector {
public:
   // Ctor and Dtor
   EXKalDetector(Int_t m = 100); //m is the maximum number of measurement layers
   ~EXKalDetector();

private:

   ClassDef(EXKalDetector,1)   // Sample hit class
};

#endif
