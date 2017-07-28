#ifndef EXDETECTOR_H
#define EXDETECTOR_H

#include "TVector3.h"         // from ROOT
#include "EXVKalDetector.h"   // from KalTrackLib
#include "EXMeasLayer.h"

//----------------------------------------------------------------------
// RTPC Parameters Start
//----------------------------------------------------------------------
static const double kRTPC_R_Target = 0.3;
static const double kRTPC_R_Cathode = 3.0;
static const double kRTPC_R_GEM1 = 7.0;
static const double kRTPC_Length = 40.0;
static const int    kNDetDummyLayer = 7;

#define RTPC_TIC 120    //TIC width is 200 ns

#if (RTPC_TIC == 200)
  static const int    kNDetLayer = 35;
  static const double kDetLayerRList[] = {
  6.916, 6.849, 6.779, 6.707, 6.632, 6.555, 6.474, 6.392, 6.306, 6.218, 
  6.127, 6.033, 5.937, 5.838, 5.736, 5.632, 5.525, 5.415, 5.303, 5.187, 
  5.070, 4.949, 4.826, 4.700, 4.572, 4.441, 4.307, 4.170, 4.031, 3.889, 
  3.745, 3.597, 3.448, 3.295, 3.140
};  // this list in unit of cm

#elif (RTPC_TIC == 120) 
  static const int    kNDetLayer = 42;
  static const double kDetLayerRList[] = {
  6.952, 6.857, 6.762, 6.667, 6.571, 6.476, 6.381, 6.286, 6.190, 6.095,
  6.000, 5.905, 5.810, 5.714, 5.619, 5.524, 5.429, 5.333, 5.238, 5.143,
  5.048, 4.952, 4.857, 4.762, 4.667, 4.571, 4.476, 4.381, 4.286, 4.190,
  4.095, 4.000, 3.905, 3.810, 3.714, 3.619, 3.524, 3.429, 3.333, 3.238,
  3.143, 3.048
};  // this list in unit of cm

#else 
  static const int    kNDetLayer = 55;
  static const double kDetLayerRList[] = {
  7.558, 7.493, 7.428, 7.363, 7.297, 7.230, 7.162, 7.094, 7.025, 6.956,
  6.886, 6.815, 6.743, 6.671, 6.598, 6.524, 6.449, 6.373, 6.296, 6.219,
  6.140, 6.060, 5.980, 5.898, 5.815, 5.730, 5.645, 5.558, 5.470, 5.380,
  5.289, 5.196, 5.101, 5.005, 4.906, 4.806, 4.703, 4.598, 4.491, 4.381,
  4.267, 4.151, 4.032, 3.908, 3.781, 3.649, 3.512, 3.369, 3.220, 3.064,
  2.898, 2.723, 2.534, 2.331, 2.107
};  // this list in unit of cm
#endif

//----------------------------------------------------------------------
// RTPC Parameters End
//----------------------------------------------------------------------

class EXKalDetector : public EXVKalDetector {
public:
   // Ctor and Dtor
   EXKalDetector(Int_t m = 100); //m is the maximum number of measurement layers
   ~EXKalDetector();

private:

   ClassDef(EXKalDetector,1)   // Sample hit class
};

#endif
