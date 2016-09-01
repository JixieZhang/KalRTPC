#include "EXEventGen.h"
#include "EXKalDetector.h"
#include "EXMeasLayer.h"
#include "TPlane.h"
#include "TRandom.h"

#include "EXHit.h"
#include "BonusHelixFit.hh"
#include "GlobalDebuger.hh"

//-----------------------------------
// RTPC Parameters
//-----------------------------------
static const double kRTPC_R_GEM1 = 7.0;
static const double kRTPC_R_Cathode = 3.0;
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

//#define _ExEventGenDebug_ 1

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
  if(phi_c> PI) phi_c-=2*PI;
  if(phi_c<-PI) phi_c+=2*PI;
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
      <<"deg, Theta0_p="<<Theta0_p*57.3<<"deg  Z0="<<z0<<endl;
  }
#endif

  return aTrack;
}

//print the helix information for a given helix track point
void EXEventGen::PrintHelix(THelicalTrack *aTrack, const char *title)
{
  const double PI=acos(0.0)*2;
  // ---------------------------
  //  print a helical track
  // ---------------------------
  double tanLambda = aTrack->GetTanLambda(); 
  double cpa = aTrack->GetKappa(); 
  double fi0 = aTrack->GetPhi0(); 
  double rho = aTrack->GetRho();
  double A = aTrack->GetXc();
  double B = aTrack->GetYc();
  //double phi_c=atan2(B,A);
  double phi_c=(cpa>0)?fi0:fi0+PI;
  if(phi_c> PI) phi_c-=2*PI;
  if(phi_c<-PI) phi_c+=2*PI;
  double phi_p=(cpa>0.) ? phi_c+PI/2 : phi_c-PI/2;
  if(phi_p> PI) phi_p-=2*PI;
  if(phi_p<-PI) phi_p+=2*PI;

  double pt = fabs(1.0/cpa);
  //double pz = pt * tanLambda;
  double p  = pt * sqrt(1+tanLambda*tanLambda);   //p = pt / sinTheta 
  double th = (tanLambda>0) ? asin(pt/p) : PI-asin(pt/p);

  cout<<"\n+---------------------------------------------------------------------+\n";
  cout<<"  title="<<title<<endl;
  cout<<"  Kappa="<<cpa<<", tanLambda="<<tanLambda<<"\n";
  cout<<"  Rho="<<rho<<", A="<<A<<", B="<<B<<", phi_c="<<phi_c*57.3<<"deg, fi0="<<fi0*57.3<<"deg "<<endl;
  cout<<"  pt="<<pt<<"  P_p="<<p<<", Phi_p="<<phi_p*57.3<<"deg, Theta_p="<<th*57.3<<endl;
  cout<<"+---------------------------------------------------------------------+\n";
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
      <<"  P0="<<P0_p<<",  Pt="<<P0_p*sin(Theta0_p)<<", Theta0="
      <<Theta0_p*57.3<<", Phi0="<<Phi0_p*57.3<<"  Z0="<<Z0/10.<<"cm"<<endl;
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
//based on the helix definition: tanLambda = ctanTheta
//when rho>0, fi0 definition is different by PI, and dfi
//is in diff sign
//This routine has been fully debuged, it generates helix
//without energy loss or MSC, from (0,0,z0=0)
//One could use random z0 too 
int  EXEventGen::GenCircle(double pt_min, double pt_max, 
  double costh_min, double costh_max)
{
  const double PI = atan(1.)*4;
  double bfield_tesla = dynamic_cast<const EXKalDetector &>
    (dynamic_cast<EXMeasLayer *>
    (fCradlePtr->At(0))->GetParent(kFALSE)).GetBfield()/10;

  // r_m = pt_gev/(0.3*Bz_tesla);
  double pt = gRandom->Uniform(pt_min, pt_max);  //in gev
  double costh = gRandom->Uniform(costh_min, costh_max); 
  double sinth = sqrt(1.0-costh*costh);
  double tanlambda = costh/sinth;
  double rho = pt/(0.3*bfield_tesla) * 100;  //in cm
  double r = fabs(rho);
  double z0 = 0;

  //phi_c is circle center phi angle in hall coordinate system
  double phi_c = 2*TMath::Pi()*gRandom->Uniform(); 
  double a = r * cos(phi_c);
  double b = r * sin(phi_c);
  double fi0 = phi_c-PI;
  if(pt>0) fi0+=PI;
  if(fi0> PI) fi0-=2*PI;
  if(fi0<-PI) fi0+=2*PI;

  double phi0_p = (pt>0) ? phi_c+PI/2 : phi_c-PI/2;
  if(phi0_p> PI) phi0_p-=2*PI;
  if(phi0_p<-PI) phi0_p+=2*PI;

  //base on phi angle in circle center coordinate to calculate lab x,y,z
  //note that 
  //if (2r<S_cathode), no hit
  //if (2r>S_gem1) then only loop one iteration
  //if (2r<=S_gem1) then will curve back, 2 hits each layer
  double phi_cir;
  HitNum=0;
  //swim forward
  for(int i=0;i<kNDetLayer;i++) {
    if (kDetLayerRList[kNDetLayer-1-i]>2.*r) continue;   //no hit
    double dfi = asin(kDetLayerRList[kNDetLayer-1-i]/2./r)*2;
    if(pt>0) dfi*=-1;   //from definition
    phi_cir = fi0+dfi;

    StepX[HitNum]=(-rho*cos(phi_cir)+a)*10.;	    //in mm
    StepY[HitNum]=(-rho*sin(phi_cir)+b)*10.;	    //in mm
    StepZ[HitNum]=(z0-rho*tanlambda*dfi)*10.;	    //in mm
    StepS[HitNum]=sqrt(StepX[HitNum]*StepX[HitNum]+StepY[HitNum]*StepY[HitNum]);
    StepPhi[HitNum]=atan2(StepY[HitNum],StepX[HitNum]);
    HitNum++;
  }
  
  //swim backward
  if(kRTPC_R_GEM1>2.*r)
  {
    for(int i=0;i<kNDetLayer;i++) {
      if (kDetLayerRList[i]>2.*r) continue;   //no hit
      double dfi = 2*PI-asin(kDetLayerRList[i]/2./r)*2;
      if(pt>0) dfi*=-1; //from definition
      phi_cir = fi0+dfi;

      StepX[HitNum]=(-rho*cos(phi_cir)+a)*10.;	//in mm
      StepY[HitNum]=(-rho*sin(phi_cir)+b)*10.;	//in mm
      StepZ[HitNum]=(z0-rho*tanlambda*dfi)*10.;		//in mm
      StepS[HitNum]=sqrt(StepX[HitNum]*StepX[HitNum]+StepY[HitNum]*StepY[HitNum]);
      StepPhi[HitNum]=atan2(StepY[HitNum],StepX[HitNum]);
      HitNum++;
    }
  }
  //pass these information to NtReader buffer so KalRTPC can load them 
  //into the output tree
  this->X0=0*10.;
  this->Y0=0*10.;
  this->Z0=0*10.;
  this->Phi0_p=phi0_p;
  this->Theta0_p=acos(costh);
  this->P0_p=fabs(pt)/sinth;

#ifdef _ExEventGenDebug_
  //just for debug
  if(_ExEventGenDebug_>=1) {
    cout<<"\nCircle Event:  pt="<<pt<<"  Rho="<<rho
      <<", A="<<a<<", B="<<b
      <<", phi_c="<<phi_c*57.3<<"deg, fi0="<<fi0*57.3<<"deg "<<endl;
    cout<<"  P0_p="<<P0_p<<", Phi0_p="<<Phi0_p*57.3
      <<"deg, Theta0_p="<<Theta0_p*57.3<<"deg  Z0="<<0<<endl;  
  }
  if(_ExEventGenDebug_>=4) {
    for(int i=0;i<HitNum;i++) {
      cout<<"Hit "<<setw(2)<<i<<"("<<setw(8)<<StepX[i]<<", "
	<<setw(8)<<StepY[i]<<", "
	<<setw(8)<<StepZ[i]<<") ==>  S="<<setw(8)<<StepS[i]
      <<" mm  Phi="<<setw(8)<<StepPhi[i]*57.3<<" deg"<<endl;
    }
  }
#endif

  MakeHitsFromTraj(StepX,StepY,StepZ,HitNum,true);

  return HitNum;
}

//x y z in mm and in increasing order 
void EXEventGen::MakeHitsFromTraj(double *x, double *y, double *z, int npt, bool smearing)
{
  // ---------------------------
  //  use given track to make hits
  // ---------------------------
  TVector3 xx; 

  for (int i = 0; i < npt; i++) { 
    ///////////////////////////////////////////////////
    //determine which measurement layer this hit belongs to
    ///////////////////////////////////////////////////

    //Note that x,y,z from ntuple are in unit of mm
    xx.SetXYZ(x[i]/10.,y[i]/10.,z[i]/10.);
    int pLyrIndex=-1;    
    //pLyrIndex = Get Layer Index From R
    double r=xx.Perp();
    if(r>kRTPC_R_GEM1 || r<kRTPC_R_Cathode) continue;

    double pRmin, pRmax=kRTPC_R_Cathode;
    for(int j=kNDetLayer-1;j>=0;j--){
      pRmin=pRmax;
      if(j>0) pRmax=(kDetLayerRList[j]+kDetLayerRList[j-1])/2;
      else pRmax=kRTPC_R_GEM1;
      if(r>=pRmin && r<pRmax) {pLyrIndex=kNDetLayer-1-j; break;}
    }
#ifdef _ExEventGenDebug_
    if( _ExEventGenDebug_ >= 4) {
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
      if( _ExEventGenDebug_ >= 4) {
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

  
  const double PI=acos(0.0)*2;
  double pRho = aTrack.GetRho();
  Pt_3pt = pRho/aTrack.GetPtoR();
  R_3pt = pRho;
  A_3pt = aTrack.GetXc();
  B_3pt = aTrack.GetYc();
  double tanLambda = aTrack.GetTanLambda();
  P_3pt = fabs(Pt_3pt) * sqrt(1+tanLambda*tanLambda); 
  Theta_3pt = atan(1./tanLambda);
  if(Theta_3pt<0) Theta_3pt+=PI;

#ifdef _ExEventGenDebug_
  //just for debug
  double pPhi_c = atan2(B_3pt,A_3pt);
  double pFi0 = aTrack.GetPhi0(); 

  if(_ExEventGenDebug_>=2) {
    cout<<" 3-point Helix:  pt="<<Pt_3pt
      <<"  Rho="<<pRho<<", A="<<A_3pt<<", B="<<B_3pt
      <<", phi_c="<<pPhi_c*57.3<<"deg, fi0_last="<<pFi0*57.3<<"deg "<<endl;   
    cout<<"  P_3pt="<<P_3pt<<", Theta_3pt="<<Theta_3pt*57.3<<"deg"<<endl;
  }
#endif

  return aTrack;
}
   


//Apply linear regression to "-Rho*dPhi vs dZ" to determine theta and z of a helix
//according to definition, -Rho * tanLambda = dZdPhi
//tanTheta = 1/tanLambda = -Rho/dZdPhi = -Rho*dPhi/dZ
//it means that tanTheta is the slope of "-Rho*dPhi vs dZ"
//so we can do linear fit to get the slope
//if dPhi_vx is given, one can calcuate dZ_vx, then z_vx
//linear regression formula can be found here
//http://www.datagenetics.com/blog/august12013/index.html
//
void EXEventGen::FitHelixThetaZ(int npt,double szPos[][3], double Rho, double A, double B,
                                double& Theta0, double& Z0)
{
    const double PI=acos(0.0)*2;
    ///////////////////////////////////////////////////////////////////////    
    double rhodfi[200],dfi[200],dz[200];
    //angle in circle coordinate system
    double phi_c_1st=atan2(szPos[0][1]-B,szPos[0][0]-A);
    for(int i=1;i<npt;i++) {
      dz[i-1] = szPos[i][2]-szPos[0][2];
      dfi[i-1] = atan2(szPos[i][1]-B,szPos[i][0]-A) - phi_c_1st;
      if(dfi[i-1]<0) dfi[i-1] += 2*PI;
      rhodfi[i-1] = -Rho * dfi[i-1]; 
      if(i+1>=200) break;
    }

    //now do the linear regression on  "Rho*dPhi vs dZ"
    int n=npt-1;
    double M,C;    //the function is y = M*x + C
    double sumX=0,sumY=0,sumXY=0,sumX2=0;
    for(int j=0;j<n;j++) {
      sumX  += dz[j];
      sumY  += rhodfi[j];
      sumXY += dz[j]*rhodfi[j];
      sumX2 += dz[j]*dz[j];
    }
    M = (n*sumXY-sumX*sumY)/(n*sumX2-sumX*sumX);
    C = (sumY/n) - M*(sumX/n);

    //this method does not work well for theta0=90deg
#ifdef _ExEventGenDebug_
    if(_ExEventGenDebug_>=2) {
      cout<<"FitHelixThetaZ():  dfi_span="<<dfi[n-1]<<"rad, z_span="<<dz[n-1]<<"cm"
	<<",  <rhodfi>="<<sumY/n<<",  <dz>="<<sumX/n<<endl; 
    }
#endif
    if(fabs(dz[n-1])<0.4 && fabs(sumX/n)<0.2) {
#ifdef _ExEventGenDebug_
      if(_ExEventGenDebug_>=2) {
	cout<<"**FitHelixThetaZ():  dz_span too small, is within uncertianty, do nothing***\n";      
      }
#endif
      //Z0 = sumX/n;
      //Theta0 = PI/2;
      return;
    }
    
    //now convert M and C into theta and z according to definition
    //tanTheta = 1/tanLambda = Rho/dZdPhi = Rho*dPhi/dZ
    //it means that tanTheta is the slope of "Rho*dPhi vs dZ"
    double pTheta0 = atan(M);
    if(pTheta0<0) pTheta0 += PI;
    
    if((sumX/n<-0.2 && pTheta0*57.3<87) || (sumX/n>0.2 && pTheta0*57.3>93)) {
#ifdef _ExEventGenDebug_
      if(_ExEventGenDebug_>=1) {
	cout<<"**FitHelixThetaZ(): wrong fitted theta="<<pTheta0*57.3<<"deg, do nothing***"<<endl;
      }
      if(_ExEventGenDebug_>=2) {
	for(int j=0;j<n;j++) {
	  cout<<"point "<<setw(3)<<j<<":  rhodfi="<<setw(10)<<rhodfi[j]<<",  dfi="
	    <<setw(10)<<dfi[j]<<",  dz="<<setw(10)<<dz[j]<<endl; 
	}
      }
#endif
      return;
    }

    Theta0 = pTheta0;
    double phi_c_vx = atan2(0-B,0-A);
    double dfi_vx = phi_c_vx - phi_c_1st;
    double dz_vx = ( dfi_vx - C) / M;
    Z0 = dz_vx + szPos[0][2];
}

//Do global helix fit to get initial parameter for Kalman Filter
//IterDirection=true is farward, otherwise backward
THelicalTrack EXEventGen::DoHelixFit(bool IterDirection) 
{
  const double PI=acos(0.0)*2;
  //the buffer should be in  increasing order
  TIter next(fHitBufPtr, true);  //forward direction

  // ---------------------------
  //  Start to extract all hits
  // ---------------------------
  //in unit of cm
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
 
  ////////////////////////////////////////////////////////////////////////

  //do the helix fit and store all result into the tree leaves buffer
  double pRho, pA, pB, pPhi, pTheta, pX0, pY0, pZ0, pDCA, pChi2;
  int fit_to_beamline=1;
  helix_fit(npt,szPos,pRho,pA,pB,pPhi,pTheta,pX0,pY0,pZ0,pDCA,pChi2,
    fit_to_beamline);

   /////////////////////////////////////////////////////////////////////
  //TODO: 
  //check global helix fit why it return a wrong sign of rho for large curve track
  
  //Global helix fit might return the wrong sign, expecially for large curve back tracks
  //Once it happens, its theta and z are totally wrong, phi is off by PI according to definition 
  //here I determine the sign  
  //using dfi_vx2first, since it always less than PI. 
  //For clock-wise track, dfi_vx2first<0  
  //the next few lines will get the phi angle on circle system, then do a subtraction
  double phi_cir_vx = atan2(0-pB,0-pA);
  double phi_cir_first = atan2(StepY_rec_m[0]-pB,StepX_rec_m[0]-pA);
  double phi_cir_last = atan2(StepY_rec_m[npt-1]-pB,StepX_rec_m[npt-1]-pA);

  //double dfi_first2last = phi_cir_last - phi_cir_first;
  //if(dfi_first2last> PI) dfi_first2last-=2*PI;
  //if(dfi_first2last<-PI) dfi_first2last+=2*PI;
  
  double dfi_vx2first = phi_cir_first - phi_cir_vx;
  if(dfi_vx2first> PI) dfi_vx2first-=2*PI;
  if(dfi_vx2first<-PI) dfi_vx2first+=2*PI;
  
  double dfi_vx2last = phi_cir_last - phi_cir_vx;
  if(dfi_vx2last> PI) dfi_vx2last-=2*PI;
  if(dfi_vx2last<-PI) dfi_vx2last+=2*PI;
  
  int sign = (dfi_vx2first<0) ? 1 : -1;

#ifdef _ExEventGenDebug_
  //just for debug
  if(_ExEventGenDebug_>=5) {
    cout<<" First_hit=("<<StepX_rec_m[0]<<", "<<StepY_rec_m[0]<<", "<<StepZ_rec_m[0]<<"), ";
    cout<<"  Last_hit=("<<StepX_rec_m[npt-1]<<", "<<StepY_rec_m[npt-1]<<", "<<StepZ_rec_m[npt-1]<<") \n";
    cout<<"  phi_cir_vx="<<phi_cir_vx*57.3<<"  phi_cir_first="<<phi_cir_first*57.3
      <<"  phi_cir_last="<<phi_cir_last*57.3<<"  dfi_vx2first="<<dfi_vx2first*57.3
      <<"  dfi_vx2last="<<dfi_vx2last*57.3<<endl;
  } 
  if(_ExEventGenDebug_>=4) {
    cout<<"  npt="<<npt<<",  dz_span="<<StepZ_rec_m[npt-1]-StepZ_rec_m[0]
      <<"cm,  chi2="<<pChi2<<endl;
  }
#endif
  
  ////////////////////
  //just for debug, store r,theta,z to make figure 
  R_3pt = pRho;
  A_3pt = pTheta;
  B_3pt = pZ0;
  
  ////////////////////

  //make correction for global helix fit result
  if (sign*pRho<0)  {
#ifdef _ExEventGenDebug_
    cout<<"***Warning: global helix fit return wrong sign! Correct it back! \n";
    if(_ExEventGenDebug_>=2) {
      cout<<"***Before correction: Rho="<<setw(8)<<pRho<<", Phi="<<setw(8)<<pPhi*57.3
	<<"deg, Theta="<<setw(8)<<pTheta*57.3<<"deg, Z="<<pZ0<<"cm \n";
    }
#endif
    pRho *= -1.0;
    pPhi+=PI;
    if(pPhi> PI) pPhi-=2*PI;
    if(pPhi<-PI) pPhi+=2*PI;
    //todo:  need to get resonable theta and z
  
    FitHelixThetaZ(npt,szPos,pRho,pA,pB,pTheta,pZ0);
#ifdef _ExEventGenDebug_
    if(_ExEventGenDebug_>=2) {
      cout<<"*** After correction: Rho="<<setw(8)<<pRho<<", Phi="<<setw(8)<<pPhi*57.3
	<<"deg, Theta="<<setw(8)<<pTheta*57.3<<"deg, Z="<<pZ0<<"cm \n";
    }
#endif
  }
  else
  {
#ifdef _ExEventGenDebug_
    //sometimes the global helix return wrong theta, 
    //want to find out and correct it back
    double ppTheta=pTheta, ppZ0=pZ0;
    FitHelixThetaZ(npt,szPos,pRho,pA,pB,ppTheta,ppZ0);
    if((ppTheta-PI/2)*(pTheta-PI/2)<0) {
    if(_ExEventGenDebug_>=2) 
      cout<<"***Before FitHelixThetaZ: Theta="<<setw(8)<<pTheta*57.3<<"deg, Z="<<pZ0<<"cm \n";
      pTheta=ppTheta; pZ0=ppZ0;
    if(_ExEventGenDebug_>=2) 
      cout<<"*** After FitHelixThetaZ: Theta="<<setw(8)<<ppTheta*57.3<<"deg, Z="<<ppZ0<<"cm \n";
    }
#endif
  }
  
#ifdef _ExEventGenDebug_
  //just for debug
  if(_ExEventGenDebug_>=2) {
    if(fabs(Theta0_p-pTheta)*57.3>5 )
    {
      Pause4Debug();
    }
  }
#endif
  /////////////////////////////////////////////////////////////////////
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
  if(Phi_c> PI) Phi_c-=2*PI;
  if(Phi_c<-PI) Phi_c+=2*PI;

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
  Double_t fi0 = (pRho>0)?Phi_c:Phi_c-PI;   //the fi0 of vertex
  if(fi0> PI) fi0-=2*PI;
  if(fi0<-PI) fi0+=2*PI;
  Double_t dz  = 0;
  Double_t cs  = cos(Theta_rec);
  Double_t tnl = cs / TMath::Sqrt((1-cs)*(1+cs)); 
  Double_t x0  = 0;  // can also use pX0;
  Double_t y0  = 0;  // can also use pY0;
  Double_t z0  = pZ0;

  
  //IterDirection=true is farward, return helix state at 1st point
  //otherwise backward, return helix state at last point
  double fi0_last = fi0+dfi_vx2last;
  if(fi0_last> PI) fi0_last-=2*PI;
  if(fi0_last<-PI) fi0_last+=2*PI;

  double fi0_1st = fi0+dfi_vx2first;
  if(fi0_1st> PI) fi0_1st-=2*PI;
  if(fi0_1st<-PI) fi0_1st+=2*PI;

  //this helix track is not ideal, but only fi0,tnl and cpa will be used
  //by the caller
  double ret_fi0 = (IterDirection==kIterBackward)?fi0_last:fi0_1st;
  THelicalTrack aTrack(dr,ret_fi0,cpa,dz,tnl,x0,y0,z0,b);

  this->P0_rec_p=fabs(pt)/sin(pTheta);
  this->X0_rec_p=pX0*10.;
  this->Y0_rec_p=pY0*10.;
  this->Z0_rec_p=pZ0*10;
  this->Theta0_rec_p=pTheta;
  this->Phi0_rec_p=pPhi;


#ifdef _ExEventGenDebug_
  //just for debug
  if(_ExEventGenDebug_>=2) {
    cout<<"  global Helix:  pt="<<pt
      <<"  Rho="<<pRho<<", A="<<pA<<", B="<<pB
      <<", phi_c="<<Phi_c*57.3<<"deg,"
      <<" fi0_1st="<<fi0_1st*57.3<<"deg "
      <<" fi0_last="<<fi0_last*57.3<<"deg "<<endl;
    cout<<"  P_hel="<<P0_rec_p<<", Phi_hel="<<Phi0_rec_p*57.3
      <<"deg, Theta_hel="<<Theta0_rec_p*57.3<<"deg,  Z_hel="<<pZ0
      <<", fi0="<<fi0*57.3<<endl;
  }
#endif

  return aTrack;
}

