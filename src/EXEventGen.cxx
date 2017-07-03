#include "EXEventGen.h"
#include "EXKalDetector.h"
#include "EXMeasLayer.h"
#include "TPlane.h"
#include "TRandom.h"

#include "EXHit.h"
#include "GlobalDebuger.hh"

#include <stdio.h>
#include <iomanip>
#include <iostream>

using namespace std;

ClassImp(EXEventGen)

//#define _ExEventGenDebug_ 9

#ifdef _ExEventGenDebug_
#include "TBenchmark.h"
#endif

  Double_t EXEventGen::fgT0 = 14.; // [nsec]

EXEventGen::EXEventGen(TKalDetCradle &cradle, TObjArray &kalhits)
  : fCradlePtr(&cradle), fHitBufPtr(&kalhits) 
{
  //double fDetLayerRBoundary[kNDetLayer+1];
  //fDetLayerRBoundary[0]=kRTPC_R_GEM1, 
  //fDetLayerRBoundary[kNDetLayer]=kRTPC_R_Cathode, 
  // then fDetLayerRBoundary[i=1...kNDetLayer-1] := (kDetLayerRList[i-1]+kDetLayerRList[i])/2
  fDetLayerRBoundary = new double [kNDetLayer+1];
  fDetLayerRBoundary[0] = kRTPC_R_GEM1;
  fDetLayerRBoundary[kNDetLayer]=kRTPC_R_Cathode;
  for(int i=1;i<kNDetLayer-1; i++) {
    fDetLayerRBoundary[i] = (kDetLayerRList[i-1]+kDetLayerRList[i])/2;
  }
}

EXEventGen::~EXEventGen() 
{
  delete fDetLayerRBoundary;
}

//This routine will create parameters for a THelicalTrack,
//There are no hits inside, just a abstract track 
THelicalTrack EXEventGen::GenerateHelix(double pt_min, double pt_max,
  double cosmin, double cosmax, 
  double z_min, double z_max)
{
  const double kPi=acos(0.0)*2;
  // ---------------------------
  //  Generate a helical track
  // ---------------------------
  //Jixie: pivot point is (0,0,z0) this helix happens to go throught (0,0,z0)
  Double_t dr  = 0.;
  Double_t fi0 = 2*kPi*(gRandom->Uniform()-0.5);
  Double_t pt  = gRandom->Uniform(pt_min, pt_max);
  Double_t cpa = 1. / pt;
  Double_t dz  = 0.;
  Double_t cs  = gRandom->Uniform(cosmin, cosmax);
  Double_t tnl = cs / TMath::Sqrt((1-cs)*(1+cs)); 
  Double_t x0  = 0.;
  Double_t y0  = 0.;
  //Double_t z0  = 0.;
  Double_t z0  = gRandom->Uniform(z_min, z_max);

  Double_t b   = dynamic_cast<const EXKalDetector &>
    (dynamic_cast<EXMeasLayer *>
    (fCradlePtr->At(0))->GetParent(kFALSE)).GetBfield();

  THelicalTrack aTrack(dr,fi0,cpa,dz,tnl,x0,y0,z0,b);

  //pass these information into buffer so KalRTPC can load them 
  //into the output tree

  //double phi_c=atan2(this->B_rec,this->A_rec);
  double phi_c=(cpa>0)?fi0:fi0+kPi;
  if(phi_c> kPi) phi_c-=2*kPi;
  if(phi_c<-kPi) phi_c+=2*kPi;
  Phi0=(cpa>0.) ? phi_c+kPi/2 : phi_c-kPi/2;
  if(Phi0> kPi) Phi0-=2*kPi;
  if(Phi0<-kPi) Phi0+=2*kPi;

  this->X0=x0;
  this->Y0=y0;
  this->Z0=z0;
  this->Theta0=acos(cs);
  this->P0=fabs(pt)/TMath::Sqrt((1-cs)*(1+cs));


#ifdef _ExEventGenDebug_
  //just for debug
  if(_ExEventGenDebug_>=1) {
    cout<<"\n Helix Event:  pt="<<pt<<"  Rho="<<aTrack.GetRho()
      <<", A="<<aTrack.GetXc()<<", B="<<aTrack.GetYc()
      <<", phi_c="<<phi_c*57.3<<"deg, fi0="<<fi0*57.3<<"deg "<<endl;
    cout<<"  P0_p="<<P0<<", Phi0_p="<<Phi0*57.3
      <<"deg, Theta0_p="<<Theta0*57.3<<"deg  Z0="<<z0<<endl;
  }
#endif

  return aTrack;
}

//print the helix information for a given helix track point
void EXEventGen::PrintHelix(THelicalTrack *aTrack, const char *title)
{
  const double kPi=acos(0.0)*2;
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
  double phi_c=(cpa>0)?fi0:fi0+kPi;
  if(phi_c> kPi) phi_c-=2*kPi;
  if(phi_c<-kPi) phi_c+=2*kPi;
  double phi_p=(cpa>0.) ? phi_c+kPi/2 : phi_c-kPi/2;
  if(phi_p> kPi) phi_p-=2*kPi;
  if(phi_p<-kPi) phi_p+=2*kPi;

  double pt = fabs(1.0/cpa);
  //double pz = pt * tanLambda;
  double p  = pt * sqrt(1+tanLambda*tanLambda);   //p = pt / sinTheta 
  double th = (tanLambda>0) ? asin(pt/p) : kPi-asin(pt/p);

  cout<<"\n+---------------------------------------------------------------------+\n";
  cout<<"  title="<<title<<endl;
  cout<<"  Kappa="<<cpa<<", tanLambda="<<tanLambda<<"\n";
  cout<<"  Rho="<<rho<<", A="<<A<<", B="<<B<<", phi_c="<<phi_c*57.3<<"deg, fi0="<<fi0*57.3<<"deg "<<endl;
  cout<<"  pt="<<pt<<"  P_p="<<p<<", Phi_p="<<phi_p*57.3<<"deg, Theta_p="<<th*57.3<<endl;
  cout<<"+---------------------------------------------------------------------+\n";
}


void EXEventGen::Swim(THelicalTrack &heltrk, Bool_t bIncludeCurveBackHits, Double_t mass)
{
  // ---------------------------
  //  Swim track and Make hits
  // ---------------------------
  //-----------------------------------------
  //reset the values
  Rho_1st=TanLambda_1st=Phi0_1st=0.0;
  Rho_last=TanLambda_last=Phi0_last=0.0;
  //-----------------------------------------
  Double_t dfi       = -dynamic_cast<TVSurface *>(fCradlePtr->At(0))
    ->GetSortingPolicy()
    / heltrk.GetRho();

  Int_t    nlayers   = fCradlePtr->GetEntries();
  Int_t    dlyr      = 1;
  Double_t dfisum    = 0.;
  StepNum=0;

  for (Int_t lyr = 0; lyr >= 0; lyr += dlyr) { // loop over layers
    // change direction if it starts looping back
    if(!bIncludeCurveBackHits) {
      //by Jixie: do not include curve back hits
      if (lyr > nlayers -1 ) break; 
    }
    else {
      if (lyr == nlayers - 1) dlyr = -1;
    }

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
      //Jixie: need to take information out from here and store into the root file
      StepX[StepNum]=xx.X();StepY[StepNum]=xx.Y();StepZ[StepNum]=xx.Z();
      StepPhi[StepNum]=xx.Phi();StepS[StepNum]=xx.Perp();
      StepNum++;
      if(StepNum>=MaxHit) break;

      //here I try to store rho, tanLambda and fi0 for the 1st and last hit
      //But I can not tell when the track reach the last step.
      //Note that when the track die in the drift region, heltrk.GetRho() return NAN
      //it is not a number any more. Therefore I have to store it at each step
      //I use 'fabs(tmpRho)>0.01' to tell tmpRho is NAN or not
      double tmpRho=heltrk.GetRho();
      if(fabs(Rho_1st)<0.01 && fabs(tmpRho)>0.01) {
        Rho_1st = tmpRho;
        TanLambda_1st = heltrk.GetTanLambda();
        Phi0_1st = heltrk.GetPhi0(); 
      }

      if(fabs(tmpRho)>0.01) {
        Rho_last = tmpRho;
        TanLambda_last = heltrk.GetTanLambda();
        Phi0_last = heltrk.GetPhi0(); 
      }
    }
    if (lyr == nlayers - 1) break;
  }

}

//generate a circle center at (a,b) and go through (0,0)
//based on the helix definition: tanLambda = ctanTheta
//when rho>0, fi0 definition is different by kPi, and dfi
//is in diff sign
//This routine has been fully debuged, it generates helix
//without energy loss or MSC, from (0,0,z0=0)
//One could use random z0 too 
int  EXEventGen::GenerateCircle(double pt_min, double pt_max, double costh_min, double costh_max,
  double z_min, double z_max, bool bIncludeCurveBackHits)
{
  const double kPi = atan(1.)*4;
  double bfield_tesla = dynamic_cast<const EXKalDetector &>
    (dynamic_cast<EXMeasLayer *>
    (fCradlePtr->At(0))->GetParent(kFALSE)).GetBfield()/10;

  // r_m = pt_gev/(0.3*Bz_tesla);
  double pt = gRandom->Uniform(pt_min, pt_max);  //in gev
  double costh = gRandom->Uniform(costh_min, costh_max); 
  double sinth = sqrt(1.0-costh*costh);
  double tanlambda = costh/sinth;
  double rho = pt/(0.3*bfield_tesla) * 100;      //in cm
  double r = fabs(rho);
  double z0 = gRandom->Uniform(z_min, z_max);

  //phi_c is circle center phi angle in hall coordinate system
  double phi_c = 2*kPi*gRandom->Uniform(); 
  double a = r * cos(phi_c);
  double b = r * sin(phi_c);
  double fi0 = phi_c-kPi;
  if(pt>0) fi0+=kPi;
  if(fi0> kPi) fi0-=2*kPi;
  if(fi0<-kPi) fi0+=2*kPi;

  //phi0_p is the phi angle of the momentum at vertex
  double phi0_p = (pt>0) ? phi_c+kPi/2 : phi_c-kPi/2;
  if(phi0_p> kPi) phi0_p-=2*kPi;
  if(phi0_p<-kPi) phi0_p+=2*kPi;


  //record these variables for root tree
  //since there is no MS or ELoss, rho and tanlambda should not change
  Rho_1st=Rho_last=rho;
  TanLambda_1st=TanLambda_last=tanlambda;

  //base on phi angle in circle center coordinate to calculate lab x,y,z
  //note that 
  //if (2r<S_cathode), no hit
  //if (2r>S_gem1) then only loop one iteration
  //if (2r<=S_gem1) then will curve back,6 2 hits each layer
  double phi_cir;
  StepNum=0;
  //swim forward
  for(int i=0;i<kNDetLayer;i++) {
    if (kDetLayerRList[kNDetLayer-1-i]>2.*r) continue;   //no hit
    double dfi = asin(kDetLayerRList[kNDetLayer-1-i]/2./r)*2;
    if(pt>0) dfi*=-1;   //from definition
    phi_cir = fi0+dfi;

    //record these variables for root tree
    if(StepNum==0) {Phi0_1st=phi_cir;}

    StepX[StepNum]=-rho*cos(phi_cir)+a;	    //in cm
    StepY[StepNum]=-rho*sin(phi_cir)+b;	    //in cm
    StepZ[StepNum]=z0-rho*tanlambda*dfi;	  //in cm
    StepS[StepNum]=sqrt(StepX[StepNum]*StepX[StepNum]+StepY[StepNum]*StepY[StepNum]);
    StepPhi[StepNum]=atan2(StepY[StepNum],StepX[StepNum]);
    StepNum++;
    if(StepNum>=MaxHit) break;
  }

  //swim backward
  if(bIncludeCurveBackHits) {
    if(kRTPC_R_GEM1>2.*r) {
      for(int i=0;i<kNDetLayer;i++) {
        if (kDetLayerRList[i]>2.*r) continue;   //no hit
        double dfi = 2*kPi-asin(kDetLayerRList[i]/2./r)*2;
        if(pt>0) dfi*=-1; //from definition
        phi_cir = fi0+dfi;

        StepX[StepNum]=-rho*cos(phi_cir)+a;	    //in cm
        StepY[StepNum]=-rho*sin(phi_cir)+b;	    //in cm
        StepZ[StepNum]=z0-rho*tanlambda*dfi;	  //in cm
        StepS[StepNum]=sqrt(StepX[StepNum]*StepX[StepNum]+StepY[StepNum]*StepY[StepNum]);
        StepPhi[StepNum]=atan2(StepY[StepNum],StepX[StepNum]);
        StepNum++;
        if(StepNum>=MaxHit) break;
      }
    }
  }
  //pass these information to NtReader buffer so KalRTPC can load them 
  //into the output tree
  this->X0=0;
  this->Y0=0;
  this->Z0=z0;
  this->Phi0=phi0_p;
  this->Theta0=acos(costh);
  this->P0=fabs(pt)/sinth;

#ifdef _ExEventGenDebug_
  //just for debug
  if(_ExEventGenDebug_>=1) {
    cout<<"\nCircle Event:  pt="<<pt<<"  Rho="<<rho<<", A="<<a<<", B="<<b
      <<", phi_c="<<phi_c*57.3<<"deg, fi0="<<fi0*57.3<<"deg "<<endl;
    cout<<"  P0_p="<<P0<<", Phi0_p="<<Phi0*57.3
      <<"deg, Theta0_p="<<Theta0*57.3<<"deg  Z0="<<Z0<<endl;  
  }
  if(_ExEventGenDebug_>=4) {
    for(int i=0;i<StepNum;i++) {
      cout<<"Hit "<<setw(2)<<i<<"("<<setw(8)<<StepX[i]<<", "
        <<setw(8)<<StepY[i]<<", "
        <<setw(8)<<StepZ[i]<<") ==>  S="<<setw(8)<<StepS[i]
      <<" mm  Phi="<<setw(8)<<StepPhi[i]*57.3<<" deg"<<endl;
    }
  }
#endif

  MakeHitsFromTraj(StepX,StepY,StepZ,StepNum,true);

  return StepNum;
}


//binary searches for a value in an decreasing sorted array
//   arr is an array to search in, in decreasing order
// value is searched value
//  left is an index of left boundary
// right is an index of right boundary
//returns position of searched value, if it presents in the array
//        or -1*position_it_should_be_incerted, if it is absent
int EXEventGen::BinarySearch(const double *arr, double value, int left, int right) 
{
  while (left <= right) {
    int middle = (left + right) / 2;
    if (arr[middle] == value) {
      return middle;
    } else if (arr[middle] < value) {
      right = middle - 1;
    } else {
      left = middle + 1;
    }
  }
  return -((left + right) / 2 + 1);
}


//x y z in cm and in increasing time order 
void EXEventGen::MakeHitsFromTraj(double *x, double *y, double *z, int _npt_,
  bool smearing, bool bIncludeCurveBackHits)
{
  // ------------------------------------------------------
  //  use given track to make hits
  // ------------------------------------------------------
  int npt=0;

  if(!bIncludeCurveBackHits) { 
    //Do this block if you only keep forward going points
    double tmpR=0.0,tmpRmax=0.0;
    for (int jj=0; jj<_npt_; jj++) { 
      tmpR = sqrt(pow(x[jj],2)+pow(y[jj],2));
      //use 1 mm margin to determine if the track curve back or not
      if (tmpR+0.1 < tmpRmax) break; 
      if (tmpR>tmpRmax) tmpRmax=tmpR;
      npt++;
    }
  } else {
    npt = _npt_;
  }

  StepNum=0;
  TVector3 xx; 
  for (int i = 0; i < npt; i++) {
    ///////////////////////////////////////////////////
    //determine which measurement layer this hit belongs to
    ///////////////////////////////////////////////////
    //Note that kDetLayerRList is in decreasing order
    //By Jixie @20170308: 
    //After benchmark test, I found that, for kNDetLayer=35, binary search will
    //run faster than bruit-force search. It is not faster than expected simply 
    //because that there are only 35 detector layers and these hits are already sorted.
    //Assuming a track have 35 total hits and one hit in each layer, the total time 
    //it takes for search 10^8 times for the whole track are ~70 seconds for  
    // loop search and ~35 seconds for binary search.
    //
    //I also test kNDetLayer=70, in this case,  binary search is much fast. 
    // kNDetLayer=35, in each 10^8 search, binary serach takes 1.01 seconds
    // kNDetLayer=70, in each 10^8 search, binary serach takes 1.25 seconds
    // loop search searches to step 10 for 10^8 times will take 0.95 seconds
    // loop search searches to step 20 for 10^8 times will take 1.78 seconds
    // loop search searches to step 31 for 10^8 times will take 3.60 seconds
    // loop search searches to step 35 for 10^8 times will take 3.90 seconds
    // loop search searches to step 40 for 10^8 times will take 4.42 seconds
    // loop search searches to step 51 for 10^8 times will take 5.43 seconds
    // loop search searches to step 61 for 10^8 times will take 6.30 seconds
    // loop search searches to step 69 for 10^8 times will take 7.31 seconds
    //
    //Assuming a track have 70 total hits and one hit in each layer, the total time 
    //it takes for search 10^8 times for the whole track are ~259.5 seconds for  
    // loop search and ~87.5 seconds for binary search.
    //If kNDetLayer==21,  binary search uses the same time as loop search does.
    //
    //since the hits are almost sorted, I optimized the loop-search. it 
    //is much faster than binary search. It takes only 30-40% of the time  
    //binary search does.

    xx.SetXYZ(x[i],y[i],z[i]);
    int pLyrIndex=-1;    
    //pLyrIndex = Get Layer Index From R
    double r=xx.Perp();
    if(r>kRTPC_R_GEM1 || r<kRTPC_R_Cathode) continue;

    //using optimzed loop search to determine the detector layer index
    int jj = 6+int((kRTPC_R_GEM1-r)/(kRTPC_R_GEM1-kRTPC_R_Cathode)*kNDetLayer);
    if(jj>kNDetLayer) jj=kNDetLayer;
    for(int j=jj;j>0;j--){
      if(r>=fDetLayerRBoundary[j] && r<fDetLayerRBoundary[j-1]) {
        pLyrIndex=kNDetLayer-j; break;
      }
    }

#if defined _ExEventGenDebug_ && (_ExEventGenDebug_>=9)
    //Do bench mark test on brute force search and binuary search
    TBenchmark pBenchmark;
    int _ret_=0;

    pBenchmark.Start("stress_binaryserach");
    for(int k=0;k<100000000;k++) {
      ///////////////////
      _ret_=BinarySearch(fDetLayerRBoundary,r,0,kNDetLayer);
      pLyrIndex=kNDetLayer-abs(_ret_);
      ///////////////////
    }
    pBenchmark.Stop("stress_binaryserach");
    pBenchmark.Print("stress_binaryserach");


    pBenchmark.Start("stress_loopserach");
    for(int k=0;k<100000000;k++) {
      ///////////////////
      for(int j=kNDetLayer;j>0;j--){
        if(r>=fDetLayerRBoundary[j] && r<fDetLayerRBoundary[j-1]) {
          pLyrIndex=kNDetLayer-j; break;
        }
      }
      ///////////////////
    }
    pBenchmark.Stop("stress_loopserach");
    pBenchmark.Print("stress_loopserach");

    int pLyrIndex_opt=-1, jstart=0;
    pBenchmark.Start("stress_loopserach_opt");
    for(int k=0;k<100000000;k++) {
      ///////////////////
      jstart = 6+int((kRTPC_R_GEM1-r)/(kRTPC_R_GEM1-kRTPC_R_Cathode)*kNDetLayer);
      if(jstart>kNDetLayer) jstart=kNDetLayer;
      for(int j=jstart;j>0;j--){
        if(r>=fDetLayerRBoundary[j] && r<fDetLayerRBoundary[j-1]) {
          pLyrIndex_opt=kNDetLayer-j; break;
        }
      }
      ///////////////////
    }
    pBenchmark.Stop("stress_loopserach_opt");
    pBenchmark.Print("stress_loopserach_opt");

    if (_ExEventGenDebug_>=10) {
      for(int j=0;j<=kNDetLayer;j++) { printf("%6d ",j);}
      printf("\n"); 
      for(int j=0;j<=kNDetLayer;j++) { printf("%6.2f ",fDetLayerRBoundary[j]);}
      printf("\n"); 
    }
    cout<<"BinarySearch(): r="<<r<<" return="<<_ret_<<"  jstart="<<jstart
      <<"  loop_pLyrIndex_opt="<<pLyrIndex_opt
      <<"  loop_pLyrIndex="<<pLyrIndex
      <<"  binary_pLyrIndex="<<kNDetLayer-abs(_ret_)<<endl;
#endif


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

      //Jixie: need to take information out from here and store into the root file
      StepX[StepNum]=xx.X();StepY[StepNum]=xx.Y();StepZ[StepNum]=xx.Z();
      StepPhi[StepNum]=xx.Phi();StepS[StepNum]=xx.Perp();
      StepNum++;
      if(StepNum>=MaxHit) break;

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

//x y z in mm and in increasing time order 
void EXEventGen::MakeHitsFromTraj_mm(double *x_mm, double *y_mm, double *z_mm, 
  int npt, bool smearing, bool bIncludeCurveBackHits)
{
  double *x=new double[npt]; 
  double *y=new double[npt]; 
  double *z=new double[npt]; 

  //convert from mm to cm
  for (int i = 0; i < npt; i++) { 
    x[i]=x_mm[i]/10.0;
    y[i]=y_mm[i]/10.0;
    z[i]=z_mm[i]/10.0;
  }

  MakeHitsFromTraj(x, y, z, npt, smearing, bIncludeCurveBackHits);	

  delete x;
  delete y;
  delete z;
} 
