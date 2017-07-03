#include "EXKalDetector.h"
#include "EXMeasLayer.h"
#include "EXHit.h"
#include "math.h"
#include <iostream>

#include "TRandom.h" // from ROOT

const int    kNDetDummyLayer = 6;
#define RTPC_TIC 120    //TIC width is 200 ns

#if (RTPC_TIC == 200)
const int    kNDetLayer = 35;
const double kDetLayerRList[] = {
  6.916, 6.849, 6.779, 6.707, 6.632, 6.555, 6.474, 6.392, 6.306, 6.218, 
  6.127, 6.033, 5.937, 5.838, 5.736, 5.632, 5.525, 5.415, 5.303, 5.187, 
  5.070, 4.949, 4.826, 4.700, 4.572, 4.441, 4.307, 4.170, 4.031, 3.889, 
  3.745, 3.597, 3.448, 3.295, 3.140
};  // this list in unit of cm

#elif (RTPC_TIC == 120) 
const int    kNDetLayer = 42;
const double kDetLayerRList[] = {
  6.952, 6.857, 6.762, 6.667, 6.571, 6.476, 6.381, 6.286, 6.190, 6.095,
  6.000, 5.905, 5.810, 5.714, 5.619, 5.524, 5.429, 5.333, 5.238, 5.143,
  5.048, 4.952, 4.857, 4.762, 4.667, 4.571, 4.476, 4.381, 4.286, 4.190,
  4.095, 4.000, 3.905, 3.810, 3.714, 3.619, 3.524, 3.429, 3.333, 3.238,
  3.143, 3.048
};  // this list in unit of cm

#else 
const int    kNDetLayer = 70;
const double kDetLayerRList[] = {
  6.944, 6.887, 6.831, 6.775, 6.718, 6.662, 6.606, 6.549, 6.493, 6.437,
  6.380, 6.324, 6.268, 6.211, 6.155, 6.099, 6.042, 5.986, 5.930, 5.873,
  5.817, 5.761, 5.704, 5.648, 5.592, 5.535, 5.479, 5.423, 5.366, 5.310,
  5.254, 5.197, 5.141, 5.085, 5.028, 4.972, 4.915, 4.859, 4.803, 4.746,
  4.690, 4.634, 4.577, 4.521, 4.465, 4.408, 4.352, 4.296, 4.239, 4.183,
  4.127, 4.070, 4.014, 3.958, 3.901, 3.845, 3.789, 3.732, 3.676, 3.620,
  3.563, 3.507, 3.451, 3.394, 3.338, 3.282, 3.225, 3.169, 3.113, 3.056
};  // this list in unit of cm
#endif

ClassImp(EXKalDetector)


EXKalDetector::EXKalDetector(Int_t m)
  : EXVKalDetector(m)
{
  //RTPC material table
  //Name        density(mg/cm3)     effective_radlen(mm)   Rin-Rout (mm)
  //7ATM_D2      1.17464            1.07246E6              0 - 3
  //Kapton       1.42E3             285.758                3 - 3.050
  //1ATM_He4     0.1665             5.665E6                3.050 - 20.0
  //Al-Mylay-Al  1.412              289.52                 20.0 - 20.00407
  //1ATM_He4DME  3.33504            1.63207E6              20.00407 - 30
  //Al-Mylay-Al  1.412              289.52                 30.0 - 30.00407
  //1ATM_He4DME  3.33504            1.63207E6              30.00407 - 70

  // effective atomic number, Zeff,  using Murty's paper: Nature 207, 398-399 (24 July 1965)
  //R. C. Murty, "Effective atomic numbers of heterogeneous materials", 
  //also see https://en.wikipedia.org/wiki/Effective_atomic_number
  //    Z_eff = pow(Sum{ f_i * Zi^2.94 },1/2.94), where f_i = n_i*Z_i/Z_total
  //==> Z_eff = pow(Sum{n_i*Z_i^3.94}/Z_total,1/2.94 )
  double A,Z,density,radlen,A_eff,Z_eff,Sum_A,Sum_Z;
  
  /*
  //vacuum 
  A       =  4.0026;                           // mass number
  Z       =  2.;                               // atomic number
  density = 0.1665E-13;                        // [g/cm^3]
  radlen  = 5.665E15;                          // [cm]
  TMaterial &Vacuum = *new TMaterial("Vacuum", "Vacuum", 
    A, Z, density, radlen, 0.);
  */

  //BONUS material table,
  //D2Gas at 7 atm 300k
  A       =  2.;                               // mass number
  Z       =  1.;                               // atomic number
  density = 1.17464E-3;                        // [g/cm^3]
  radlen  = 1.07246E5;                         // [cm]
  TMaterial &D2Gas = *new TMaterial("D2Gas", "D2Gas_7ATM_300k", 
    A, Z, density, radlen, 0.);

  //kapton, A and Z are effective numebers
  //(C22H10N2O5)n
  Sum_A = 12.*22+1.*10+14.*2+16.*5;
  Sum_Z =  6.*22+1.*10+ 7.*2+ 8.*5;
  Z_eff = pow((22*pow(6.,3.94)+10.+2*pow(7.,3.94)+5*pow(8.,3.94))/Sum_Z,1/2.94);
  A_eff = Z_eff * Sum_A / Sum_Z;
  A       = A_eff;                             // mass number
  Z       = Z_eff;                             // atomic number
  density = 1.42;                              // [g/cm^3]
  radlen  = 28.5758;                           // [cm]
  TMaterial &Kapton = *new TMaterial("Kapton", "Kapton", 
    A, Z, density, radlen, 0.);

  //He4 gas at 1 atm 300k
  A       =  4.0026;                           // mass number
  Z       =  2.;                               // atomic number
  density = 0.1665E-3;                         // [g/cm^3]
  radlen  = 5.665E5;                           // [cm]
  TMaterial &He4Gas = *new TMaterial("He4Gas", "He4_1ATM_300k", 
    A, Z, density, radlen, 0.);

  //aluminized mylar (C10H8O4)n
  Sum_A = 12.*10+1.*8+16.*4;
  Sum_Z =  6.*10+1.*8+ 8.*4;
  Z_eff = pow((10*pow(6.,3.94)+8.+4*pow(8.,3.94))/Sum_Z,1/2.94);
  A_eff = Z_eff * Sum_A / Sum_Z;
  A       = A_eff;                             // mass number
  Z       = Z_eff;                             // atomic number
  density = 1.412;                             // [g/cm^3]
  radlen  = 28.952;                            // [cm]
  TMaterial &AlMylar = *new TMaterial("AlMylar", "Al-Mylay-Al", 
    A, Z, density, radlen, 0.);

  /////////////////////////////////////////////////////////////////
  //He4+DME (C2H6O,density=2.11mg/cm^3 at 273k) gas mixture (8:2=4:1) 
  //at 1 atm 300k
  Sum_A = 12.*2+1.*6+16.*1+4.*4;
  Sum_Z =  6.*2+1.*6+ 8.*1+2.*4;
  Z_eff = pow((2*pow(6.,3.94)+6.+1*pow(8.,3.94)+4*pow(2.,3.94))/Sum_Z,1/2.94);
  A_eff = Z_eff * Sum_A / Sum_Z;
  A       =  A_eff;                             // mass number
  Z       =  Z_eff;                             // atomic number
  density = 0.51722E-3;                         // [g/cm^3]
  radlen  = 1.63207E5;                          // [cm]
  TMaterial &He4DME = *new TMaterial("He4DME", "He4_Mix_DME_80to20_1ATM_300k", 
    A, Z, density, radlen, 0.);

  
  //He4+CO2 (CO2,density=1.977mg/cm^3 at 273k) gas mixture (8:2=4:1) 
  //at 1 atm 300k
  Sum_A = 12.*1+16.*2+4.*4;
  Sum_Z =  6.*1+ 8.*2+2.*4;
  Z_eff = pow((1*pow(6.,3.94)+2*pow(8.,3.94)+4*pow(2.,3.94))/Sum_Z,1/2.94);
  A_eff = Z_eff * Sum_A / Sum_Z;

  //calculate mixtrue density
  double mixden;
  double den_he4=0.1786E-3*273/300.;
  double den_co2=2.11E-3*273/300.;
  mixden = den_he4*0.8+den_co2*0.2;

  //calculate radlength
  double mixradlen;   //using this formula: 1/X0 = sum (Wj/Xj)
  double wr_he4 = den_he4*0.8/mixden;  //weigth ratio
  double wr_co2 = den_co2*0.2/mixden;
  double radlen_he4 = 94.32/den_he4;
  double radlen_co2 = 36.20/den_co2;
  mixradlen = 1.0 / (wr_he4/radlen_he4 + wr_co2/radlen_co2);
  A       =  A_eff;                             // mass number
  Z       =  Z_eff;                             // atomic number
  density = mixden;                             // [g/cm^3]
  radlen  = mixradlen;                          // [cm]
  TMaterial &He4CO2 = *new TMaterial("He4CO2", "He4_Mix_CO2_80to20_1ATM_300k", 
    A, Z, density, radlen, 0.);
  std::cout<<"He4CO2: A_eff="<<A_eff<<",  Z_eff="<<Z_eff<<", radlen="<<mixradlen<<" cm\n";
  ////////////////////////////////////////////////////////////////////////////
  //start to build  detectors

  Int_t    nlayers = kNDetLayer; // # sampling layers 
  Double_t lhalf   = 20.;        // half length

  // Create dummy layers of the inner cylinder of the central tracker
  double r;   //in cm
 
   Bool_t active = EXMeasLayer::kActive;
   Bool_t dummy  = EXMeasLayer::kDummy;

  //add target gas and target wall and He4 gas
  Add(new EXMeasLayer(D2Gas, Kapton, r=0.3, lhalf, dummy));
  //Add(new EXMeasLayer(Kapton, He4Gas, r=0.3028, lhalf, EXVMeasLayer::kDummy));
  Add(new EXMeasLayer(Kapton, He4Gas, r=0.3050, lhalf, dummy));

  //add first aluminized mylar foil
  Add(new EXMeasLayer(He4Gas, AlMylar, r=2.0, lhalf, dummy));
  Add(new EXMeasLayer(AlMylar, He4DME, r=2.000407, lhalf, dummy));

  //add 2nd aluminized mylar foil
  Add(new EXMeasLayer(He4DME, AlMylar, r=3.0, lhalf, dummy));
  Add(new EXMeasLayer(AlMylar, He4DME, r=3.000407, lhalf, dummy));

  // Create measurement layers in the drift region
  for (Int_t layer = nlayers-1; layer >= 0; layer--) {
    EXMeasLayer* ml = new EXMeasLayer(He4DME, He4DME, kDetLayerRList[layer], lhalf, 
      active);
    ml->SetIndex(kNDetDummyLayer+nlayers-1-layer);
    Add(ml);
  }
  SetOwner(); 
}

EXKalDetector::~EXKalDetector()
{
  ;
}
