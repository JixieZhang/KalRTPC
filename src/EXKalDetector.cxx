#include "EXKalDetector.h"
#include "EXMeasLayer.h"
#include "EXHit.h"
#include "math.h"

#include "TRandom.h" // from ROOT
//By jixie: the system think the particle carry negative charge by default
//It doesnot allow me to give it negative field
//Need to spend time to figure out why and solve it!
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

  // effective atomic number, Zeff,  using Murty's paper: Nature 207, 398-399 (24 July 1965)
  //R. C. Murty, "Effective atomic numbers of heterogeneous materials", 
  //also see https://en.wikipedia.org/wiki/Effective_atomic_number
  //    Z_eff = pow(Sum{ f_i * Zi^2.94 },1/2.94), where f_i = n_i*Z_i/Z_total
  //==> Z_eff = pow(Sum{n_i*Z_i^3.94}/Z_total,1/2.94 )
  double A,Z,density,radlen,A_eff,Z_eff,Sum_A,Sum_Z;
  
  //vacuum 
  A       =  4.0026;                           // mass number
  Z       =  2.;                               // atomic number
  density = 0.1665E-13;                        // [g/cm^3]
  radlen  = 5.665E15;                          // [cm]
  TMaterial &Vacuum = *new TMaterial("Vacuum", "Vacuum", 
    A, Z, density, radlen, 0.);


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
  TMaterial &He4DME = *new TMaterial("He4DME", "He4_Mix_DME_1to9_1ATM_300k", 
    A, Z, density, radlen, 0.);



  Int_t    nlayers = kNDetLayer; // # sampling layers 
  Double_t lhalf   = 20.;        // half length

  // Create dummy layers of the inner cylinder of the central tracker
  double r;   //in cm
 
  //add target gas and target wall and He4 gas
  Add(new EXMeasLayer(D2Gas, Kapton, r=0.3, lhalf, EXMeasLayer::kDummy));
  Add(new EXMeasLayer(Kapton, He4Gas, r=0.3028, lhalf, EXMeasLayer::kDummy));
  //Add(new EXMeasLayer(Kapton, He4Gas, r=0.3050, lhalf, EXMeasLayer::kDummy));

  //add first aluminized mylar foil
  Add(new EXMeasLayer(He4Gas, AlMylar, r=2.0, lhalf, EXMeasLayer::kDummy));
  Add(new EXMeasLayer(AlMylar, He4DME, r=2.000407, lhalf, EXMeasLayer::kDummy));

  //add 2nd aluminized mylar foil
  Add(new EXMeasLayer(He4DME, AlMylar, r=3.0, lhalf, EXMeasLayer::kDummy));
  Add(new EXMeasLayer(AlMylar, He4DME, r=3.000407, lhalf, EXMeasLayer::kDummy));
   /*
  //add target gas and target wall and He4 gas
  Add(new EXMeasLayer(Vacuum, Vacuum, r=0.3, lhalf, EXMeasLayer::kDummy));
  Add(new EXMeasLayer(Vacuum, Vacuum, r=0.3028, lhalf, EXMeasLayer::kDummy));

  //add first aluminized mylar foil
  Add(new EXMeasLayer(Vacuum, Vacuum, r=2.0, lhalf, EXMeasLayer::kDummy));
  Add(new EXMeasLayer(Vacuum, Vacuum, r=2.000407, lhalf, EXMeasLayer::kDummy));

  //add 2nd aluminized mylar foil
  Add(new EXMeasLayer(Vacuum, Vacuum, r=3.0, lhalf, EXMeasLayer::kDummy));
  Add(new EXMeasLayer(Vacuum, Vacuum, r=3.000407, lhalf, EXMeasLayer::kDummy));
  */

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
