//*************************************************************************
//* ===================
//*  EXMeasLayer Class
//* ===================
//*
//* (Description)
//*   Sample measurement layer class used by EXHit.
//* (Requires)
//* (Provides)
//*     class EXMeasLayer
//* (Update Recored)
//*   2016/02/23  Jixie Zhang, Modify for RTPC
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//

#include "EXMeasLayer.h"
#include "EXHit.h"
#include "EXKalDetector.h"
#include "TRandom.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
using namespace std;

#define _EXMeasLayerDebug_ 4

Bool_t   EXMeasLayer::kActive = kTRUE;
Bool_t   EXMeasLayer::kDummy = kFALSE;

ClassImp(EXMeasLayer)

EXMeasLayer::EXMeasLayer(
  TMaterial &min,      // material inside the layer
  TMaterial &mout,     // material outside the layer
  Double_t   r0,       // radius of the cylindrial layer
  Double_t   lhalf,    // half length of the cylindrical layer
  Bool_t     isactive) // flag to tell the layer is active
  : TVMeasLayer(min, mout, isactive),
  TCylinder(r0, lhalf)
{
  //Lorentz angle:  Phi(S) = 0.0278*S^2 - 0.5334*S + 2.3475  //in rad and cm
  //Ionization Location: S(t) = -0.0355*t^2 - 0.3208*t + 6.9481    //in us and cm
  //Drifting Time:  t(S) = -0.1642*S^2 -0.0947*S + 8.8001

  //==>  dPhi/dS = 0.0556*S - 0.5334
  //     dS/dt = -0.071*t - 0.3208
  //-->  dPhi/dt = (0.0556*S - 0.5334) (-0.071*t - 0.3208)
  //

  double S = r0;
  double t = -0.1642*(S*S) -0.0947*S + 8.8001;
  double dt = 0.1 ;  //100 ns
  double dS = fabs((-0.071*t - 0.3208)) * dt;
  double dPhi = fabs((0.0556*S - 0.5334)) * dS;

  fgSigmaX = S*dPhi;    //r*dphi 
  //fgSigmaX = dPhi;    //dphi 
  fgSigmaZ = 0.2;       //dz = 2 mm

#ifdef  _EXMeasLayerDebug_ 
  if(_EXMeasLayerDebug_>=2) {
    if(isactive) {
      cout<<"Construct active measurement layer:"
	<<"  TDC="<<setw(2)<<int(t/0.2)<<fixed
	<<"  r="<<setw(6)<<setprecision(4)<<S 
	<<"  dr="<<setw(6)<<setprecision(4)<<dS//<<endl
	<<"  dPhi="<<setw(6)<<setprecision(4)<<dPhi*57.3<<"(deg)"
	<<"  d(rphi)="<<setw(7)<<setprecision(5)<<fgSigmaX
	<<"  dz="<<setw(5)<<setprecision(3)<<fgSigmaZ
	<< setprecision(7)<< resetiosflags(ios::showpoint)<<endl;
    }
  }
#endif
}

EXMeasLayer::~EXMeasLayer()
{
}

TKalMatrix EXMeasLayer::XvToMv(const TVector3 &xv) const
{
  // Calculate hit coordinate information:
  //	mv(0,0) = r * phi 
  //     (1,0) = z

  TKalMatrix mv(kMdim,1);
  mv(0,0)  = GetR() * TMath::ATan2(xv.Y(), xv.X());
  mv(1,0)  = xv.Z();
  return mv;
}

TKalMatrix EXMeasLayer::XvToMv(const TVTrackHit &,
  const TVector3   &xv) const
{
  return XvToMv(xv);
}

TVector3 EXMeasLayer::HitToXv(const TVTrackHit &vht) const
{
  const EXHit &ht = dynamic_cast<const EXHit &>(vht);

  Double_t phi = ht(0,0) / GetR();
  Double_t z   = ht(1,0);
  Double_t x   = GetR() * TMath::Cos(phi);
  Double_t y   = GetR() * TMath::Sin(phi);

  return TVector3(x,y,z);
}

void EXMeasLayer::CalcDhDa(
  const TVTrackHit &,          // Hit: not used in this example
  const TVector3   &xxv,       // hit position vector
  const TKalMatrix &dxphiada,  // @x(\phi(a),a)/@a
  TKalMatrix &H)  const  // projector matrix = @h/@a
{
  // Calculate
  //    H = (@h/@a) = (@phi/@a, @z/@a)^t
  // where
  //        h(a) = (phi, z)^t: expected meas vector
  //        a = (drho, phi0, kappa, dz, tanl, t0)
  //

  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5,sdim-1);

  Double_t xv = xxv.X();
  Double_t yv = xxv.Y();
  Double_t xxyy = xv * xv + yv * yv;

  // Set H = (@h/@a) = (@d/@a, @z/@a)^t

  for (Int_t i=0; i<hdim; i++) {
    H(0,i) = - (yv / xxyy) * dxphiada(0,i) 
      + (xv / xxyy) * dxphiada(1,i);
    H(0,i) *= GetR();
    H(1,i) =  dxphiada(2,i);
  }
  if (sdim == 6) {
    H(0,sdim-1) = 0.;
    H(1,sdim-1) = 0.;
  }
}

void EXMeasLayer::ProcessHit(const TVector3    &xx,
  TObjArray   &hits, bool smearing)
{
  //By Jixie:
  //XvToMv() will use phi angle from vector xx and radius r from 
  //the mesuarement layer to calaulate expected hit xx', then it
  //will use xx' to construct h. Therefore one can not 
  //get the original xx when you use MvToXv(h), what it return is 
  //the 'expected hit' vector on the helix, xx'
  //to solve this issue, I add "TVector3 fRawXv;" into class EXHit
  //to store the raw (x,y,z)

  //Jixie: shift the measurement layer such that the hit is on the surface
  //try to solve the site discard problem by this
  this->SetR(xx.Perp());

  TKalMatrix h    = XvToMv(xx);
  Double_t   rphi = h(0, 0);
  Double_t   z    = h(1, 0);

  Double_t dx = fgSigmaX;
  Double_t dz = fgSigmaZ;
  if(smearing) {
    rphi += gRandom->Gaus(0., dx);   // smearing rphi
    z    += gRandom->Gaus(0., dz);   // smearing z
  }

  Double_t meas [2];
  Double_t dmeas[2];
  meas [0] = rphi;
  meas [1] = z;
  dmeas[0] = dx;
  dmeas[1] = dz;

  Double_t b = EXKalDetector::GetBfield();
  EXHit *aHit=new EXHit(*this, meas, dmeas, b);
  aHit->SetRawXv(xx);
  hits.Add(aHit);
}

