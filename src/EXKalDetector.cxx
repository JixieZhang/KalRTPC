#include "EXKalDetector.h"
#include "EXMeasLayer.h"
#include "EXHit.h"

#include "TRandom.h" // from ROOT
//By jixie: the system think the particle carry negative charge by default
//It doesnot allow me to give it negative field
Double_t EXKalDetector::fgBfield = 50.;

const int    kNDetDummyLayer = 6;
const int    kNDetLayer = 35;
const double kDetLayerRList[] = {
  6.916, 6.849, 6.779, 6.707, 6.632, 6.555, 6.474, 6.392, 6.306, 6.218, 
  6.127, 6.033, 5.937, 5.838, 5.736, 5.632, 5.525, 5.415, 5.303, 5.187, 
  5.070, 4.949, 4.826, 4.700, 4.572, 4.441, 4.307, 4.170, 4.031, 3.889, 
  3.745, 3.597, 3.448, 3.295, 3.140
};  // this list in unit of cm


ClassImp(EXKalDetector)


EXKalDetector::EXKalDetector(Int_t m)
  : TVKalDetector(m)
{
  //RTPC material table
  //Name        density(mg/cm3)     effective_radlen(mm)   Rin-Rout (mm)
  //7ATM_D2      1.17464            1.07246E6              0 - 3
  //Kapton       1.42E3             285.758                3 - 0.028
  //1ATM_He4     0.1665             5.665E6                3.028 - 20.0
  //Al-Mylay-Al  1.412              289.52                 20.0 - 20.00407
  //1ATM_He4DME  3.33504            1.63207E6              20.00407 - 30
  //Al-Mylay-Al  1.412              289.52                 30.0 - 30.00407
  //1ATM_He4DME  3.33504            1.63207E6              30.00407 - 70

  double A,Z,density,radlen;
  //BONUS material table, A and Z are fake numbers
  //D2Gas at 7 atm 300k
  A       =  2.;                               // mass number
  Z       =  1.;                               // atomic number
  density = 1.17464E-3;                        // [g/cm^3]
  radlen  = 1.07246E5;                         // [cm]
  TMaterial &D2Gas = *new TMaterial("D2Gas", "D2Gas_7ATM_300k", 
    A, Z, density, radlen, 0.);

  //kapton, A and Z is fake numebr
  A       = 16.0107;                           // mass number
  Z       =  7.0;                              // atomic number
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

  //aluminized mylar
  A       = 28.0107;                           // mass number
  Z       =  14.;                              // atomic number
  density = 1.412;                             // [g/cm^3]
  radlen  = 28.952;                            // [cm]
  TMaterial &AlMylar = *new TMaterial("AlMylar", "Al-Mylay-Al", 
    A, Z, density, radlen, 0.);

  //He4+DME gas mixture (1:9) at 1 atm 300k
  A       =  14.0026;                           // mass number
  Z       =  7.5;                               // atomic number
  density = 3.33504E-3;                         // [g/cm^3]
  radlen  = 1.63207E5;                          // [cm]
  TMaterial &He4DME = *new TMaterial("He4DME", "He4_Mix_DME_1to9_1ATM_300k", 
    A, Z, density, radlen, 0.);



  Int_t    nlayers = kNDetLayer; // # sampling layers 
  Double_t lhalf   = 20.;        // half length

  // Create dummy layers of the inner cylinder of the central tracker
  double r;   //in cm

  //add target gas and target wall and He4 gas
  Add(new EXMeasLayer(D2Gas, Kapton, r=0.3, lhalf, EXMeasLayer::kDummy));
  Add(new EXMeasLayer(Kapton, He4Gas, r=0.3028, lhalf, EXMeasLayer::kDummy));

  //add first aluminized mylar foil
  Add(new EXMeasLayer(He4Gas, AlMylar, r=2.0, lhalf, EXMeasLayer::kDummy));
  Add(new EXMeasLayer(AlMylar, He4DME, r=2.000407, lhalf, EXMeasLayer::kDummy));

  //add 2nd aluminized mylar foil
  Add(new EXMeasLayer(He4DME, AlMylar, r=3.0, lhalf, EXMeasLayer::kDummy));
  Add(new EXMeasLayer(AlMylar, He4DME, r=3.000407, lhalf, EXMeasLayer::kDummy));

  // Create measurement layers in the drift region
  for (Int_t layer = nlayers-1; layer >= 0; layer--) {
    EXMeasLayer* ml = new EXMeasLayer(He4DME, He4DME, kDetLayerRList[layer], lhalf, 
      EXMeasLayer::kActive);
    ml->SetIndex(kNDetDummyLayer+nlayers-1-layer);
    Add(ml);
  }
  SetOwner(); 
}

EXKalDetector::~EXKalDetector()
{
  ;
}
