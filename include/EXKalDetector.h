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
static const int    kNDetDummyLayer = 6;

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
  static const int    kNDetLayer = 70;
  static const double kDetLayerRList[] = {
  6.944, 6.887, 6.831, 6.775, 6.718, 6.662, 6.606, 6.549, 6.493, 6.437,
  6.380, 6.324, 6.268, 6.211, 6.155, 6.099, 6.042, 5.986, 5.930, 5.873,
  5.817, 5.761, 5.704, 5.648, 5.592, 5.535, 5.479, 5.423, 5.366, 5.310,
  5.254, 5.197, 5.141, 5.085, 5.028, 4.972, 4.915, 4.859, 4.803, 4.746,
  4.690, 4.634, 4.577, 4.521, 4.465, 4.408, 4.352, 4.296, 4.239, 4.183,
  4.127, 4.070, 4.014, 3.958, 3.901, 3.845, 3.789, 3.732, 3.676, 3.620,
  3.563, 3.507, 3.451, 3.394, 3.338, 3.282, 3.225, 3.169, 3.113, 3.056
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
