#include "EXEventGen.h"
#include "EXKalDetector.h"
#include "EXMeasLayer.h"
#include "TPlane.h"
#include "TRandom.h"

#include "EXHit.h"
#include "BonusHelixFit.hh"

//-----------------------------------
// Track Parameters
//-----------------------------------

#define __DR__     0.
#define __FI0__    0.
#define __DZ__     0.
#define __X0__     0.
#define __Y0__     0.
#define __Z0__     0.

ClassImp(EXEventGen)

#define _ExEventGenDebug_ 1

  Double_t EXEventGen::fgT0 = 14.; // [nsec]


THelicalTrack EXEventGen::GenerateHelix(double pt_min, double pt_max,
  double cosmin, double cosmax)
{
  const double PI=acos(0.0)*2;
  // ---------------------------
  //  Generate a helical track
  // ---------------------------
  //Jixie: pivot point is (0,0,0) this helix happens to go throught (0,0,0)
  Double_t dr  = __DR__;
  Double_t fi0 = __FI0__ + 2*PI*(gRandom->Uniform()-0.5);
  Double_t pt  = gRandom->Uniform(pt_min, pt_max);
  Double_t cpa = 1. / pt;
  Double_t dz  = __DZ__;
  Double_t cs  = gRandom->Uniform(cosmin, cosmax);
  Double_t tnl = cs / TMath::Sqrt((1-cs)*(1+cs)); 
  Double_t x0  = __X0__;
  Double_t y0  = __Y0__;
  Double_t z0  = __Z0__;

  Double_t b   = dynamic_cast<const EXKalDetector &>
    (dynamic_cast<EXMeasLayer *>
    (fCradlePtr->At(0))->GetParent(kFALSE)).GetBfield();

  THelicalTrack aTrack(dr,fi0,cpa,dz,tnl,x0,y0,z0,b);

  //pass these information to NtReader buffer so KalRTPC can load them 
  //into the output tree

  //double phi_c=atan2(this->B_rec,this->A_rec);
  double phi_c=(cpa>0)?fi0:fi0+PI;
  Phi0_p=(cpa>0.) ? phi_c+PI/2 : phi_c-PI/2;
  if(Phi0_p> PI) Phi0_p-=2*PI;
  if(Phi0_p<-PI) Phi0_p+=2*PI;

  this->X0=x0*10.;
  this->Y0=y0*10.;
  this->Z0=z0*10.;
  this->Theta0_p=acos(cs);
  this->P0_p=fabs(pt)/TMath::Sqrt((1-cs)*(1+cs));


#ifdef _ExEventGenDebug_
  //just for debug
  if(_ExEventGenDebug_>=1) {
    cout<<"\n Helix Event:  pt="<<pt<<"  Rho="<<aTrack.GetRho()
      <<", A="<<aTrack.GetXc()<<", B="<<aTrack.GetYc()
      <<", phi_c="<<phi_c*57.3<<"deg, fi0="<<fi0*57.3<<"deg "<<endl;
    cout<<"  P0_p="<<P0_p<<", Phi0_p="<<Phi0_p*57.3
      <<"deg, Theta0_p="<<Theta0_p*57.3<<"deg "<<endl;
  }
#endif

  return aTrack;
}

void EXEventGen::Swim(THelicalTrack &heltrk, Double_t mass)
{
  // ---------------------------
  //  Swim track and Make hits
  // ---------------------------

  Double_t dfi       = -dynamic_cast<TVSurface *>(fCradlePtr->At(0))
    ->GetSortingPolicy()
    / heltrk.GetRho();

  Int_t    nlayers   = fCradlePtr->GetEntries();
  Int_t    dlyr      = 1;
  Double_t dfisum    = 0.;

  for (Int_t lyr = 0; lyr >= 0; lyr += dlyr) { // loop over layers
    // change direction if it starts looping back
    if (lyr == nlayers - 1) dlyr = -1;

    EXMeasLayer &ml = *dynamic_cast<EXMeasLayer *>(fCradlePtr->At(lyr));
    TVSurface   &ms = *dynamic_cast<TVSurface *>(fCradlePtr->At(lyr));
    TVector3 xx;
    Double_t dfis = dfi;
    /*  
    //By Jixie:
    Int_t TCylinder::CalcXingPointWith(const TVTrack  &hel,
    TVector3 &xx, Double_t &phi, Int_t mode, Double_t eps) const
    */
    //TCylinder::CalcXingPointWith(), input is hel, output are xx,phi
    //It will calculate if helix hel have any crossing points with current
    //TCylinder, if it has, it will store the most closed crossing point into xx
    //and also the deflection angle dfi into phi
    if (!ms.CalcXingPointWith(heltrk,xx,dfi,1)
      || TMath::Abs(dfi) > TMath::Pi()
      || TMath::Abs(dfi + dfisum) > TMath::TwoPi()) {
	dfi = dfis;
	continue;
    }
    // should use the material behind the surface since dfi is measured 
    // from the last point to the current surface
    Bool_t   dir = dlyr < 0 ? kTRUE : kFALSE;

    const TMaterial &mat = ml.GetMaterial(dir);
    if (fCradlePtr->IsMSOn()) {
      TKalMatrix Qms(5,5);
      fCradlePtr->CalcQms(mat, heltrk, dfi, Qms, mass);
      Double_t sgphi  = TMath::Sqrt(Qms(1,1));
      Double_t sgtnl  = TMath::Sqrt(Qms(4,4));
      Double_t delphi = gRandom->Gaus(0.,sgphi);
      Double_t deltnl = gRandom->Gaus(0.,sgtnl);
#if 0
      dfi *= 0.5;
      TVector3 x0ms = heltrk.CalcXAt(dfi);
      heltrk.MoveTo(x0ms,dfi);     // M.S. at mid point

      heltrk.ScatterBy(delphi,deltnl);
      dfis = dfi;
#else
      heltrk.ScatterBy(delphi,deltnl); // multiple scattering
      dfis = 0.;
#endif
      // recalculate crossing point
      if (!ms.CalcXingPointWith(heltrk,xx,dfi,1)
	|| TMath::Abs(dfi) > TMath::Pi()
	|| TMath::Abs(dfi + dfisum) > TMath::TwoPi()) {
	  dfi = dfis;
	  continue;
      }
      dfis += dfi;
    }
    dfisum += dfi;

    heltrk.MoveTo(xx,dfi);	// move pivot to current hit

    if (fCradlePtr->IsDEDXOn()) {
      TKalMatrix av(5,1);
      heltrk.PutInto(av);
      av(2,0) += fCradlePtr->GetEnergyLoss(mat, heltrk, dfis, mass); // energy loss
      heltrk.SetTo(av, heltrk.GetPivot());
    }
    if (ml.IsActive()) {
      ml.ProcessHit(xx, *fHitBufPtr,true); // create hit point
    }
    if (lyr == nlayers - 1) break;
  }
}

int  EXEventGen::LoadOneTrack()
{
  int nhits = 0;
  double xx[200],yy[200],zz[200];
  while (nhits<5) {
    int n = NtReader::LoadATrack();
    if( n == -1 ) {
      cout<<"Reach the end of input root file \n";
      return -1;
    }
    nhits=0;
    for(int i=0;i<HitNum_m;i++) {
      if(StepID_m[i]>0) {
	xx[nhits]=StepX_rec_m[i];
	yy[nhits]=StepY_rec_m[i];
	zz[nhits]=StepZ_rec_m[i];
	nhits++;
	if(nhits>=200) break;
      }
    }
  }

#ifdef _ExEventGenDebug_
  //just for debug
  if(_ExEventGenDebug_>=1) {
    cout<<"\nNtuple Event "<<setw(5)<<Index<<":  HitNum_m="<<setw(2)<<HitNum_m
      <<",  Smax="<<setw(8)<<Smax<<",  Smin="<<setw(8)<<Smin<<endl
      <<"  P0="<<P0_p<<", Theta0="<<Theta0_p*57.3<<", Phi0="<<Phi0_p*57.3<<endl;
  }
  if(_ExEventGenDebug_>=4) {
    for(int i=0;i<HitNum_m;i++) {
      cout<<"Hit "<<setw(2)<<i<<"("<<setw(8)<<StepX_rec_m[i]<<", "
	<<setw(8)<<StepY_rec_m[i]<<", "
	<<setw(8)<<StepZ_rec_m[i]<<") ==>  S="<<setw(8)<<StepS_rec_m[i]
      <<" mm  Phi="<<setw(8)<<StepPhi_rec_m[i]*57.3<<" deg"<<endl;
    }
  }
#endif

  MakeHitsFromTraj(xx,yy,zz,nhits,false);
  return HitNum_m;
}

//generate a circle center at (a,b) and go through (0,0)
int  EXEventGen::GenCircle(double pt_min, double pt_max)
{
  const double PI = atan(1.)*4;
  double bfield_tesla = dynamic_cast<const EXKalDetector &>
    (dynamic_cast<EXMeasLayer *>
    (fCradlePtr->At(0))->GetParent(kFALSE)).GetBfield()/10;

  // r_m = pt_gev/(0.3*Bz_tesla);
  double pt = gRandom->Uniform(pt_min, pt_max);  //in gev

  double r = fabs(pt)/(0.3*bfield_tesla) * 100;  //in cm

  double cs=0.0;  //only for theta=90deg
  //phi_c is circle center phi angle in hall coordinate system
  double phi_c = 2*TMath::Pi()*gRandom->Uniform(); 
  double a = r * cos(phi_c);
  double b = r * sin(phi_c);

  double phi0_p = (pt>0) ? phi_c+PI/2 : phi_c-PI/2;
  if(phi0_p> PI) phi0_p-=2*PI;
  if(phi0_p<-PI) phi0_p+=2*PI;

  //base on phi angle in circle center coordinate to calculate lab x,y,z
  double phi_cir, phi_cir_0 = atan2(-b,-a);

  for(int i=0;i<kNDetLayer;i++) {
    if(pt<0) phi_cir = phi_cir_0+asin(kDetLayerRList[kNDetLayer-1-i]/2./r)*2;
    else phi_cir = phi_cir_0-asin(kDetLayerRList[kNDetLayer-1-i]/2./r)*2;

    StepX_rec_m[i]=(r*cos(phi_cir)+a)*10.;  //in mm
    StepY_rec_m[i]=(r*sin(phi_cir)+b)*10.;  //in mm
    StepZ_rec_m[i]=0;
  }
  HitNum_m=kNDetLayer;

  //pass these information to NtReader buffer so KalRTPC can load them 
  //into the output tree
  this->X0=0*10.;
  this->Y0=0*10.;
  this->Z0=0*10.;
  this->Phi0_p=phi0_p;
  this->Theta0_p=acos(cs);
  this->P0_p=fabs(pt)/sqrt((1-cs)*(1+cs));

  this->R_rec=r*10.;
  this->A_rec=a*10.;
  this->B_rec=b*10.;
  this->Theta_rec=Theta0_p;
  this->Phi_rec=phi0_p;
  this->Z_rec=0*10.;


#ifdef _ExEventGenDebug_
  //just for debug
  if(_ExEventGenDebug_>=1) {
    cout<<"\nCircle Event:  pt="<<pt<<"  Rho="<<r
      <<", A="<<a<<", B="<<b
      <<", phi_c="<<phi_c*57.3<<"deg "<<endl;
    cout<<"  P0_p="<<P0_p<<", Phi0_p="<<Phi0_p*57.3
      <<"deg, Theta0_p="<<Theta0_p*57.3<<"deg "<<endl;  
  }
  if(_ExEventGenDebug_>=4) {
    for(int i=0;i<HitNum_m;i++) {
      cout<<"Hit "<<setw(2)<<i<<"("<<setw(8)<<StepX_rec_m[i]<<", "
	<<setw(8)<<StepY_rec_m[i]<<", "
	<<setw(8)<<StepZ_rec_m[i]<<") ==>  S="<<setw(8)<<StepS_rec_m[i]
      <<" mm  Phi="<<setw(8)<<StepPhi_rec_m[i]*57.3<<" deg"<<endl;
    }
  }
#endif

  MakeHitsFromTraj(StepX_rec_m,StepY_rec_m,StepZ_rec_m,HitNum_m,true);

  return HitNum_m;
}

//x y z in mm and in increasing order 
void EXEventGen::MakeHitsFromTraj(double *x, double *y, double *z, int npt, bool smearing)
{
  // ---------------------------
  //  use given track to make hits
  // ---------------------------
  TVector3 xx; 
  const double RTPC_R_GEM1 = 7.0;
  const double RTPC_R_Cathode = 3.0;

  for (int i = 0; i < npt; i++) { 
    ///////////////////////////////////////////////////
    //determine which measurement layer this hit belongs to
    ///////////////////////////////////////////////////

    //Note that x,y,z from ntuple are in unit of mm
    xx.SetXYZ(x[i]/10.,y[i]/10.,z[i]/10.);
    int pLyrIndex=-1;    
    //pLyrIndex = Get Layer Index From R
    double r=xx.Perp();
    if(r>RTPC_R_GEM1 || r<RTPC_R_Cathode) continue;

    double pRmin, pRmax=RTPC_R_Cathode;
    for(int j=kNDetLayer-1;j>=0;j--){
      pRmin=pRmax;
      if(j>0) pRmax=(kDetLayerRList[j]+kDetLayerRList[j-1])/2;
      else pRmax=RTPC_R_GEM1;
      if(r>=pRmin && r<pRmax) {pLyrIndex=kNDetLayer-1-j; break;}
    }
#ifdef _ExEventGenDebug_
    if( _ExEventGenDebug_ >= 3) {
      cout<<"Hit_R="<<r<<" --> pLyrIndex="<<pLyrIndex+kNDetDummyLayer<<endl;
    }
#endif
    if(pLyrIndex<0) {
      cout<<"Warning: can not find Measurement Layer for s="<<r<<" cm"<<endl;
      continue;
    }
    pLyrIndex += kNDetDummyLayer;

    EXMeasLayer &ml = *dynamic_cast<EXMeasLayer *>(fCradlePtr->At(pLyrIndex));
    if (ml.IsActive()) {
      ml.ProcessHit(xx, *fHitBufPtr, smearing); // create hit point     

#ifdef _ExEventGenDebug_
      if( _ExEventGenDebug_ >= 3) {
	EXHit *hitp = dynamic_cast<EXHit *> (fHitBufPtr->Last());       
	TVector3 xv=ml.HitToXv((TVTrackHit&)(*hitp));   
	TVector3 xraw=hitp->GetRawXv();
	cerr << "MeasLayer "<<setw(2)<< ml.GetIndex()
	  <<": R="<<setw(6)<< ml.GetR()<<": ";
	cerr << "Xv  =("<<setw(8)<< xv.X()<<",  "
	  <<setw(8)<<xv.Y()<<", "<<setw(8)<<xv.Z()<<"); ";
	cerr << "Xraw=("<<setw(8)<< xraw.X()<<",  "
	  <<setw(8)<<xraw.Y()<<", "<<setw(8)<<xraw.Z()<<"): \n";
      }
#endif
    }

  }
}


//Create a helix from 3 points to get initial parameter for Kalman Filter
//IterDirection=true is farward, otherwise backward
THelicalTrack EXEventGen::CreateInitialHelix(bool IterDirection) 
{
  Int_t i1, i2, i3;
  if (IterDirection == kIterBackward) {
    i3 = 0;
    i1 = fHitBufPtr->GetEntries() - 1;
    i2 = i1 / 2;
  } else {
    i1 = 0;
    i3 = fHitBufPtr->GetEntries() - 1;
    i2 = i3 / 2;
  }

  EXHit   &h1 = *dynamic_cast<EXHit *>(fHitBufPtr->At(i1));   // first hit
  EXHit   &h2 = *dynamic_cast<EXHit *>(fHitBufPtr->At(i2));   // middle hit
  EXHit   &h3 = *dynamic_cast<EXHit *>(fHitBufPtr->At(i3));   // last hit
  TVector3 x1 = h1.GetMeasLayer().HitToXv(h1);
  TVector3 x2 = h2.GetMeasLayer().HitToXv(h2);
  TVector3 x3 = h3.GetMeasLayer().HitToXv(h3);
  double bfield = h1.GetBfield();  //in kGauss
  THelicalTrack aTrack(x1, x2, x3, bfield, IterDirection); // initial helix 


#ifdef _ExEventGenDebug_
  //just for debug
  double pRho = aTrack.GetRho();
  double pPt = pRho/aTrack.GetPtoR();
  double pA = aTrack.GetXc();
  double pB = aTrack.GetYc();
  double pPhi_c = atan2(pB,pA);
  double pFi0 = aTrack.GetPhi0(); 
  if(_ExEventGenDebug_>=1) {
    cout<<" 3-point Helix:  pt="<<pPt
      <<"  Rho="<<pRho<<", A="<<pA<<", B="<<pB
      <<", phi_c="<<pPhi_c*57.3<<"deg, fi0="<<pFi0*57.3<<"deg "<<endl;   
  }
#endif

  return aTrack;
}

//DO a global helix fit to get initial parameter for Kalman Filter
THelicalTrack EXEventGen::DoHelixFit() 
{
  const double PI=acos(0.0)*2;
  TIter next(fHitBufPtr, true);  

  // ---------------------------
  //  Start to extract all hits
  // ---------------------------
  //Here I assume all hits are in unit of cm
  double szPos[200][3];
  int npt = 0;
  EXHit *hitp = dynamic_cast<EXHit *>(next());
  while (hitp) {     // loop over hits    
    //fill the global variables for root tree
    //TVector3 xraw = hitp->GetRawXv();
    const EXMeasLayer &ml = dynamic_cast<const EXMeasLayer &>(hitp->GetMeasLayer());
    TVector3 xv = ml.HitToXv(*hitp);
    StepX_rec_m[npt]=xv.X();StepY_rec_m[npt]=xv.Y();StepZ_rec_m[npt]=xv.Z();
    szPos[npt][0]=xv.X();szPos[npt][1]=xv.Y();szPos[npt][2]=xv.Z();

    npt++;
    if(npt>=200) break;
    hitp = dynamic_cast<EXHit *>(next());
  }
  HitNum_m=npt;

  //do the helix fit and store all result into the tree leaves buffer
  double pRho, pA, pB, pPhi, pTheta, pX0, pY0, pZ0, pDCA, pChi2;
  int fit_to_beamline=1;
  helix_fit(npt,szPos,pRho,pA,pB,pPhi,pTheta,pX0,pY0,pZ0,pDCA,pChi2,
    fit_to_beamline);
  R_rec=fabs(pRho)*10;
  A_rec=pA*10;
  B_rec=pB*10;
  Z_rec=pZ0*10;
  Theta_rec=pTheta;
  Phi_rec=pPhi;
  DCA_rec=pDCA*10;

  //helix center phi angle in the hall
  double Phi_c = (pRho>0.) ? Phi_rec-PI/2 : Phi_rec+PI/2;
  //do not use the next line, 
  //if A_rec or B_rec are close to zero, it will give a wrong sign
  //Phi_c = atan2(B_rec, A_rec); 

  //////////////////////////////////////
  //now create a THelixTrack
  //////////////////////////////////////

  const double kGev=1.0e9;
  const double kLightVelocity=2.99792458e8;
  Double_t b   = dynamic_cast<const EXKalDetector &>
    (dynamic_cast<EXMeasLayer *>
    (fCradlePtr->At(0))->GetParent(kFALSE)).GetBfield();

  Double_t dr  = pDCA;
  Double_t pt  = (pRho/100.)*kLightVelocity*(b/10.) / kGev;
  Double_t cpa = 1. / pt;
  Double_t fi0 = (pRho>0)?Phi_c:Phi_c-PI;
  if(fi0> PI) fi0-=2*PI;
  if(fi0<-PI) fi0+=2*PI;
  Double_t dz  = 0;
  Double_t cs  = cos(Theta_rec);
  Double_t tnl = cs / TMath::Sqrt((1-cs)*(1+cs)); 
  Double_t x0  = 0;
  Double_t y0  = 0;
  Double_t z0  = Z_rec;

  THelicalTrack aTrack(dr,fi0,cpa,dz,tnl,x0,y0,z0,b);

  this->P0_rec_p=fabs(pt)/sin(Theta_rec);
  this->X0_rec_p=aTrack.GetXc()*10.;
  this->Y0_rec_p=aTrack.GetYc()*10.;
  this->Z0_rec_p=pZ0*10;
  this->Theta0_rec_p=Theta_rec;
  this->Phi0_rec_p=Phi_rec;


#ifdef _ExEventGenDebug_
  //just for debug
  if(_ExEventGenDebug_>=1) {
    cout<<"  global Helix:  pt="<<pt
      <<"  Rho="<<pRho<<", A="<<pA<<", B="<<pB
      <<", phi_c="<<Phi_c*57.3<<"deg, fi0="<<fi0*57.3<<"deg "<<endl;
    cout<<"  P_hel="<<P0_rec_p<<", Phi_hel="<<Phi0_rec_p*57.3
      <<"deg, Theta_hel="<<Theta0_rec_p*57.3<<"deg "<<endl;
  }
#endif

  return aTrack;
}
