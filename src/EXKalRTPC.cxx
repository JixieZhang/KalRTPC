//This is the kernel class to do Kalman Filter for RTPC12
//Note that only fKalHits_Forward will be used by KF
// 
// To use this code, user need to do the following
// 1) Reset the buffer using Reset(); 
// 2) Prepare a chain in time increasing order using PrepareATrack(xxx);
// 3) Determine if it needs 2 iteration KF fitting and copy fKalHits to fKalHits_Forward
//    using JudgeFor2ndIteration(bool bIncludeCurveBackHits);
// 4) Do the fit using DoFitAndFilter(), which will fit fKalHits_Forward;
// see EXKalRTPC::Example() for details.

//////////////////////////////////////////////////////////////////////
#include "EXKalRTPC.h"
#include "BonusHelixFit.hh"
#include "GlobalDebuger.hh"
#include "CircleFitter_LM.h"

#include <iomanip>
#include <iostream>
using namespace std;

ClassImp(EXKalRTPC)

///////////////////////////////////////////////////////////////////////
//#define _EXKalRTPCDebug_ 0

static const Bool_t kDir = kIterBackward;
//static const Bool_t kDir = kIterForward;

CircleFitter_LM gLMFitter;

///////////////////////////////////////////////////////////////////////

EXKalRTPC::EXKalRTPC()
{
#ifdef _EXKalRTPCDebug_
  Global_Debug_Level=_EXKalRTPCDebug_;
#endif

  // ===================================================================
  //  Prepare a fDetector
  // ===================================================================
  fKalTrack = new TKalTrack();       // track buffer to hold the fitted track
  fKalHits  = new TObjArray();       // hit buffer to hold all hits
  fKalHits_Forward = new TObjArray();// hit buffer only forward hits
  fCradle   = new TKalDetCradle();   // detctor system
  fDetector = new EXKalDetector();   // fDetector

  fCradle->Install(*fDetector);      // install fDetector into its fCradle

#ifdef __MS_OFF__
  fCradle.SwitchOffMS();            // switch off multiple scattering
#endif

  // ===================================================================
  //  Prepare a Event Generator
  // ===================================================================
  fEventGen = new EXEventGen(*fCradle, *fKalHits);

  //initial value for the covariant matrix element
  //0.05 is tested to be the best for RTPC12 with 5T field
  //You might need to scan this value if condition change
  fCovMElement=5.0e-2;
}

EXKalRTPC::~EXKalRTPC()
{
  delete fEventGen;
  delete fDetector;
  delete fCradle;
  delete fKalHits;
  delete fKalHits_Forward;
}

int EXKalRTPC::GetVextex(THelicalTrack &hel, Double_t x_bpm, Double_t y_bpm,
			 TVector3 &xx,  Double_t &dfi, double &r_rec, double &a_rec, double &b_rec)
{
  r_rec = fabs(hel.GetRho());
  a_rec = hel.GetXc();
  b_rec = hel.GetYc();

  Double_t Rho2BPM = sqrt((a_rec-x_bpm)*(a_rec-x_bpm)+(b_rec-y_bpm)*(b_rec-y_bpm));
  Double_t pRVx = fabs(Rho2BPM - r_rec) + 1.0E-6;

  TCylinder pCylVx(pRVx,50.0,x_bpm,y_bpm,0.0);

  //get the vertex point from here, xx and dfi are outputs
  const TVTrack &trk = dynamic_cast<const TVTrack &>(hel);
  int n=pCylVx.CalcXingPointWith(trk,xx,dfi,0);
#ifdef _EXKalRTPCDebug_
  cout<<"GetVextex():  dfi="<<dfi*57.3<<"deg,   xx=("<<xx.X()<<", "<<xx.y()<<", "<<xx.Z()<<")"<<endl;
#endif
  return n;
}


int EXKalRTPC::GetVextex2(THelicalTrack &hel, double x_bpm, double y_bpm,
			  TVector3 &xx,  double &dfi, double &r_rec, double &a_rec, double &b_rec)
{
  TVector3 X0  = hel.GetPivot();
  r_rec = fabs(hel.GetRho());
  a_rec = hel.GetXc();
  b_rec = hel.GetYc();

  dfi = atan2(y_bpm-b_rec,x_bpm-a_rec)-atan2(X0.Y()-b_rec,X0.X()-a_rec);
  if( fabs(dfi) > kPi ) {
    dfi = (dfi>0) ? dfi-2*kPi : dfi+2*kPi;
  }
  xx = hel.CalcXAt(dfi);
#ifdef _EXKalRTPCDebug_
  if(Global_Debug_Level >= 1)
    cout<<"GetVextex2(): dfi="<<dfi*57.3<<"deg,   xx=("<<xx.X()<<", "<<xx.y()<<", "<<xx.Z()<<")"<<endl;
#endif
  return 1;
}

//Use the last site to swim back to the beam line
void EXKalRTPC::ReconVertex(TVKalState &state, double &p, double &pt, double &pz,
			    double &th, double &ph, double &x, double &y, double &z,
			    double &r_rec, double &a_rec, double &b_rec )
{
  TVector3 xv(-99,-99,-99);
  Double_t dfi=0;
  THelicalTrack hel = (dynamic_cast<TKalTrackState *>(&state))->GetHelix();

#ifdef _EXKalRTPCDebug_
  if(Global_Debug_Level >= 4) {
    EXEventGen::PrintHelix(&hel, "final helix at 1st point:");
  }
#endif

  //by helix definition" Pz = Pt * tanLambda, tanLambda = ctan(theta)
  double tanLambda = state(4,0); //or tanLambda = hel.GetTanLambda();
  double cpa = state(2,0);       //or cpa = hel.GetKappa();
  double fi0 = state(1,0);       //or fi0 = hel.GetPhi0();
  double rho = hel.GetRho();

  //GetVextex(hel,0,0,xv,dfi,r_rec,a_rec,b_rec);
  GetVextex2(hel,0,0,xv,dfi,r_rec,a_rec,b_rec);
  //this is the vertex
  x=xv.X();
  y=xv.Y();
  z=xv.Z();

  pt = fabs(1.0/cpa);
  pz = pt * tanLambda;
  p  = pt * sqrt(1+tanLambda*tanLambda);   //p = pt / sinTheta
  th = (tanLambda>0) ? asin(pt/p) : kPi-asin(pt/p);
  double fi0_vx = fi0 + dfi;  //vertex point fi0 angle in helix coordinate system
  //phi_vx_hel is the phi angle of vertex point on the circle coordinate system
  //by definition, positive helix fi0 is differ by pi from circle
  double phi_vx_hel=(cpa>0)?fi0_vx+kPi:fi0_vx;
  double phi_c_hall= phi_vx_hel - kPi ;  //phi angle of the helix center at hall coordinate system
  ph = (rho>0) ? phi_c_hall+kPi/2. : phi_c_hall-kPi/2.;
  if(ph> kPi) ph-=2*kPi;
  if(ph<-kPi) ph+=2*kPi;

#ifdef _EXKalRTPCDebug_
  //debug the vertex reconstruction
  if(Global_Debug_Level>=1) {
    TVector3 xv=hel.GetPivot();
    TVector3 x_fil=hel.CalcXAt(0.0);
    cout<<"   LastHit=("<<xv.X()<<", "<<xv.y()<<", "<<xv.Z()<<")"
	<<"   LastHit_fil=("<<x_fil.X()<<", "<<x_fil.y()<<", "<<x_fil.Z()<<")"
	<<",   Xc="<<hel.GetXc()<<",  Yc="<<hel.GetYc()<<endl;

    cout<<"Rec. to Vextex: p="<<p<<" ph="<<ph*57.3<<", th="<<th*57.3
	<<", rho="<<rho<<", fi0="<<fi0*57.3<<"deg,  fi0_vx="<<fi0_vx*57.3<<"deg"<<endl;
  }
#endif
}


void EXKalRTPC::Reset()
{
  //clear will only remove the pointer if it does own the object
  //the obj will stay
  fKalHits_Forward->Clear();
  fKalHits->Delete();
  
  rho_kal_ini=tnl_kal_ini=phi0_kal_ini=0.0;
  
  for(int i=0;i<HitNum;i++) {
    step_status[i]=0;
    StepX_rec[i]=StepY_rec[i]=StepZ_rec[i]=0.0;
    StepPhi_rec[i]=StepS_rec[i]=0.0;
  }

  P_rec=Pt_rec=Pz_rec=Theta_rec=Phi_rec=X_rec=Y_rec=Z_rec=0.0;
  R_rec=A_rec=B_rec=0.0;
  NDF=0;
  Chi2=0.0;

  P_3pt=Pt_3pt=Theta_3pt=0.0;
  R_3pt=A_3pt=B_3pt=0.0;

  P_hel=Phi_hel=Theta_hel=Z_hel=X_hel=Y_hel=0.0;
  R_hel=A_hel=B_hel=DCA_hel=Chi2_hel=0.0;
  
  Phi_hel_raw=Theta_hel_raw=Z_hel_raw=X_hel_raw=Y_hel_raw=0.0;
  R_hel_raw=A_hel_raw=B_hel_raw=DCA_hel_raw=Chi2_hel_raw=0.0; 
  
}

//Provide suggestion if need to apply 2nd iteration kalman filter
//if bIncludeCurveBackHits==false,  it will remove backward hits, otherwise just 
//copy all hits pointer into fKalHits_Forward
//It will also fill the smeared hit position array
//Note that the hit buffer must be sorted by time in increasing order
bool EXKalRTPC::JudgeFor2ndIteration(bool bIncludeCurveBackHits)
{
  //since the hits obj already exist in fKalHits, I simplely add their
  //pointers into fKalHits_Forward
  //the next line just make a copy for all the pointer, but not change ownership 
  //*fKalHits_Forward = *fKalHits;     
  //return false;
	
  bool bNeed2ndIter=false;
  TIter next(fKalHits, kIterForward);   
  TVector3 xv;
  double tmpS=0.0, Smax=0.0;
  int npt=0,idx=0;
  EXHit *hitp = dynamic_cast<EXHit *>(next());
  while (hitp) {     // loop over hits    
    //fill the global variables for root tree
    //xv = hitp->GetRawXv();
    const EXMeasLayer &ml = dynamic_cast<const EXMeasLayer &>(hitp->GetMeasLayer());
    xv = ml.HitToXv(*hitp);
    
    StepX_rec[npt]=xv.X();StepY_rec[npt]=xv.Y();StepZ_rec[npt]=xv.Z();
    StepPhi_rec[npt]=xv.Phi();StepS_rec[npt]=xv.Perp();
    tmpS=xv.Perp();
    if(tmpS+0.1 < Smax) {
      bNeed2ndIter=true;
      if(bIncludeCurveBackHits) {			
	fKalHits_Forward->Add(hitp); npt++;
      } else {
#ifdef _EXKalRTPCDebug_	
	if(Global_Debug_Level>=7) 
	  cout<<"Hit "<<setw(3)<<idx<<" removed! tmpS="<<tmpS<<"  Smax="<<Smax<<"\n";
#endif		
      }
    } else {
      if(tmpS>Smax) Smax=tmpS;
      fKalHits_Forward->Add(hitp); npt++;
    }
    if(npt>=MaxHit) break;
    hitp = dynamic_cast<EXHit *>(next());
    idx++;
  }
  HitNum=npt;

  //this objarry does not own these hits
  fKalHits_Forward->SetOwner(false);
  
  return bNeed2ndIter;
}

//Create a helix from 3 points to provide initial parameter for Kalman Filter
//IterDirection=true is farward, otherwise backward
THelicalTrack EXKalRTPC::GetIniHelixBy3Pts(bool IterDirection) 
{
  Int_t i1, i2, i3;
  if (IterDirection == kIterBackward) {
    i3 = 0;
    i1 = fKalHits_Forward->GetEntriesFast() - 1;
    i2 = i1 / 2;
  } else {
    i1 = 0;
    i3 = fKalHits_Forward->GetEntriesFast() - 1;
    i2 = i3 / 2;
  }

  EXHit   &h1 = *dynamic_cast<EXHit *>(fKalHits_Forward->At(i1));   // first hit
  EXHit   &h2 = *dynamic_cast<EXHit *>(fKalHits_Forward->At(i2));   // middle hit
  EXHit   &h3 = *dynamic_cast<EXHit *>(fKalHits_Forward->At(i3));   // last hit
  TVector3 x1 = h1.GetMeasLayer().HitToXv(h1);
  TVector3 x2 = h2.GetMeasLayer().HitToXv(h2);
  TVector3 x3 = h3.GetMeasLayer().HitToXv(h3);
  double bfield = h1.GetBfield();  //in kGauss
  THelicalTrack aTrack(x1, x2, x3, bfield, IterDirection); // initial helix 

  //store _3pt variable for root tree
  double pRho = aTrack.GetRho();
  Pt_3pt = pRho/aTrack.GetPtoR();
  R_3pt = pRho;
  A_3pt = aTrack.GetXc();
  B_3pt = aTrack.GetYc();
  double tanLambda = aTrack.GetTanLambda();
  P_3pt = fabs(Pt_3pt) * sqrt(1+tanLambda*tanLambda); 
  Theta_3pt = atan(1./tanLambda);
  if(Theta_3pt<0) Theta_3pt+=kPi;

#ifdef _EXKalRTPCDebug_
  //just for debug
  double pPhi_c = atan2(B_3pt,A_3pt);
  double pFi0 = aTrack.GetPhi0(); 

  if(_EXKalRTPCDebug_>=2) {
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
void EXKalRTPC::CorrHelixThetaZ(int npt,double szPos[][3], double Rho, double A, double B,
				double& Theta0, double& Z0)
{
  ///////////////////////////////////////////////////////////////////////    
  double rhodfi[200],dfi[200],dz[200];
  //angle in circle coordinate system
  double phi_c_1st=atan2(szPos[0][1]-B,szPos[0][0]-A);
  for(int i=1;i<npt;i++) {
    dz[i-1] = szPos[i][2]-szPos[0][2];
    dfi[i-1] = atan2(szPos[i][1]-B,szPos[i][0]-A) - phi_c_1st;
    if(dfi[i-1]<0) dfi[i-1] += 2*kPi;
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
#ifdef _EXKalRTPCDebug_
  if(_EXKalRTPCDebug_>=2) {
    cout<<"CorrHelixThetaZ():  dfi_span="<<dfi[n-1]<<"rad, z_span="<<dz[n-1]<<"cm"
	<<",  <rhodfi>="<<sumY/n<<",  <dz>="<<sumX/n<<endl; 
  }
#endif
  if(fabs(dz[n-1])<0.4 && fabs(sumX/n)<0.2) {
#ifdef _EXKalRTPCDebug_
    if(_EXKalRTPCDebug_>=2) {
      cout<<"**CorrHelixThetaZ():  dz_span too small, is within uncertianty, do nothing***\n";      
    }
#endif
    //Z0 = sumX/n;
    //Theta0 = kPi/2;
    return;
  }

  //now convert M and C into theta and z according to definition
  //tanTheta = 1/tanLambda = Rho/dZdPhi = Rho*dPhi/dZ
  //it means that tanTheta is the slope of "Rho*dPhi vs dZ"
  double pTheta0 = atan(M);
  if(pTheta0<0) pTheta0 += kPi;

  if((sumX/n<-0.2 && pTheta0*57.3<87) || (sumX/n>0.2 && pTheta0*57.3>93)) {
#ifdef _EXKalRTPCDebug_
    if(_EXKalRTPCDebug_>=1) {
      cout<<"**CorrHelixThetaZ(): wrong fitted theta="<<pTheta0*57.3<<"deg, do nothing***"<<endl;
    }
    if(_EXKalRTPCDebug_>=2) {
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


//Do global helix fit and apply my corrections
//return chi2
double EXKalRTPC::DoGlobalHelixFit(double *x, double *y,double *z, int _npt_, 
                                   bool bIncludeCurveBackHits) 
{
  //////////////////////////////////////////////////////////////////////
  int npt=0;
  double szPos[200][3];
  double tmpR=0.0, tmpRmax=0.0;
  for(int j=0;j<_npt_;j++) {
    //in case you do not want to include curve back hits!!!
    if(!bIncludeCurveBackHits) {
      tmpR = sqrt(x[j]*x[j]+y[j]*y[j]);
      if(tmpR > tmpRmax) tmpRmax = tmpR;  
      else if( tmpR+0.1 < tmpRmax) continue;
      //due to resolution, s value of some hits might be a little bit smaller 
      //than previous hit..
      //I set the margin here to be 1mm
    }
    szPos[npt][0] = x[j];
    szPos[npt][1] = y[j];
    szPos[npt][2] = z[j];
    npt++;
    if(npt>=200) break; //GHF only fit 200 hits
  }
  
  //copy the buffer to EXEventgen
  fEventGen->StepNum=npt;
  for(int j=0;j<npt;j++) {
    fEventGen->StepX[j]=x[j];
    fEventGen->StepY[j]=y[j];
    fEventGen->StepZ[j]=z[j];
    fEventGen->StepS[j]=sqrt(x[j]*x[j]+y[j]*y[j]);
    fEventGen->StepPhi[j]=atan2(y[j],x[j]);
  }

  if(npt<MinHit) return -1.0;
  //////////////////////////////////////////////////////////////////////
  
  //do the helix fit and store all results into the tree leaves buffer
  double pRho, pA, pB, pPhi, pTheta, pX0, pY0, pZ0, pDCA, pChi2=0.0;
  int fit_to_beamline=1;
  helix_fit(npt,szPos,pRho,pA,pB,pPhi,pTheta,pX0,pY0,pZ0,pDCA,pChi2,
	    fit_to_beamline);

  //store global helix result before my correction 
  Phi_hel_raw = pPhi;
  Theta_hel_raw = pTheta;
  R_hel_raw = pRho;
  A_hel_raw = pA;
  B_hel_raw = pB;
  Z_hel_raw = pZ0;
  X_hel_raw = pX0;
  Y_hel_raw = pY0;
  DCA_hel_raw = pDCA;
  Chi2_hel_raw = pChi2;


  /////////////////////////////////////////////////////////////////////
  //check global helix fit if it return a wrong sign of rho for large curve track

  //Global helix fit might return the wrong sign, expecially for large curve back tracks
  //Once it happens, its theta and z are totally wrong, phi is off by kPi according to definition 
  //here I determine the sign  
  //using dfi_vx2first, since it always less than kPi. 
  //For clock-wise track, dfi_vx2first<0  
  //the next few lines will get the phi angle on circle system, then do a subtraction
  double phi_cir_vx = atan2(0-pB,0-pA);
  double phi_cir_first = atan2(szPos[0][1]-pB,szPos[0][0]-pA);
  double phi_cir_last = atan2(szPos[npt-1][1]-pB,szPos[npt-1][0]-pA);

  //double dfi_first2last = phi_cir_last - phi_cir_first;
  //if(dfi_first2last> kPi) dfi_first2last-=2*kPi;
  //if(dfi_first2last<-kPi) dfi_first2last+=2*kPi;

  double dfi_vx2first = phi_cir_first - phi_cir_vx;
  if(dfi_vx2first> kPi) dfi_vx2first-=2*kPi;
  if(dfi_vx2first<-kPi) dfi_vx2first+=2*kPi;

  double dfi_vx2last = phi_cir_last - phi_cir_vx;
  if(dfi_vx2last> kPi) dfi_vx2last-=2*kPi;
  if(dfi_vx2last<-kPi) dfi_vx2last+=2*kPi;

  int sign = (dfi_vx2first<0) ? 1 : -1;

#ifdef _EXKalRTPCDebug_
  //just for debug
  if(_EXKalRTPCDebug_>=5) {
    cout<<" First_hit=("<<szPos[0][0]<<", "<<szPos[0][1]<<", "<<szPos[0][2]<<"), ";
    cout<<"  Last_hit=("<<szPos[npt-1][0]<<", "<<szPos[npt-1][1]<<", "<<szPos[npt-1][2]<<") \n";
    cout<<"  phi_cir_vx="<<phi_cir_vx*57.3<<"  phi_cir_first="<<phi_cir_first*57.3
	<<"  phi_cir_last="<<phi_cir_last*57.3<<"  dfi_vx2first="<<dfi_vx2first*57.3
	<<"  dfi_vx2last="<<dfi_vx2last*57.3<<endl;
  } 
  if(_EXKalRTPCDebug_>=4) {
    cout<<"  npt="<<npt<<",  dz_span="<<szPos[npt-1][2]-szPos[0][2]
	<<"cm,  chi2="<<pChi2<<endl;
  }
#endif

  ////////////////////

  //make correction for global helix fit result if needed
  if (sign*pRho<0) {
#ifdef _EXKalRTPCDebug_
    cout<<"***Warning: global helix fit return wrong sign! Correct it back! \n";
    if(_EXKalRTPCDebug_>=2) {
      cout<<"***Before correction: Rho="<<setw(8)<<pRho<<", Phi="<<setw(8)<<pPhi*57.3
	  <<"deg, Theta="<<setw(8)<<pTheta*57.3<<"deg, Z="<<pZ0<<"cm \n";
    }
#endif
    pRho *= -1.0;
    pPhi+=kPi;
    if(pPhi> kPi) pPhi-=2*kPi;
    if(pPhi<-kPi) pPhi+=2*kPi;

    //Corr theta and z
    CorrHelixThetaZ(npt,szPos,pRho,pA,pB,pTheta,pZ0);
    
#ifdef _EXKalRTPCDebug_
    if(_EXKalRTPCDebug_>=2) {
      cout<<"*** After correction: Rho="<<setw(8)<<pRho<<", Phi="<<setw(8)<<pPhi*57.3
	  <<"deg, Theta="<<setw(8)<<pTheta*57.3<<"deg, Z="<<pZ0<<"cm \n";
    }
#endif

  } else {
    //For those events that the sign of rho is right but 
#ifdef _EXKalRTPCDebug_
    //sometimes the global helix return wrong theta, 
    //want to find out and correct it back
    double ppTheta=pTheta, ppZ0=pZ0;
    CorrHelixThetaZ(npt,szPos,pRho,pA,pB,ppTheta,ppZ0);
    if((ppTheta-kPi/2)*(pTheta-kPi/2)<0) {
      pTheta=ppTheta; pZ0=ppZ0;
      if(_EXKalRTPCDebug_>=2) {
	cout<<"***Before CorrHelixThetaZ: Theta="<<setw(8)<<ppTheta*57.3<<"deg, Z="<<ppZ0<<"cm \n";
	cout<<"*** After CorrHelixThetaZ: Theta="<<setw(8)<<ppTheta*57.3<<"deg, Z="<<ppZ0<<"cm \n";
      }
    }
#endif
  }


  //corr phi and r
  //This part is missing right now because I found that GHF is already the 
  //best one in hand, CircleFitter_LM is worse.
  //you are welcome to develop one!

  //////////////////////////////////////
  //now Fit the circle using Levenberg-Marquardt method
  //////////////////////////////////////
  //sometimes this routine return very large values
  //gLMFitter.DoFit(npt,szPos,pA,pB,pRho);
  

  /////////////////////////////////////////////////////////////////////
  //store global helix result after my correction 
  Phi_hel = pPhi;
  Theta_hel = pTheta;
  R_hel = pRho;
  A_hel = pA;
  B_hel = pB;
  Z_hel = pZ0;
  X_hel = pX0;
  Y_hel = pY0;
  DCA_hel = pDCA;
  Chi2_hel = pChi2;

  
  //calculate pt and P_hel
  //Fix me:  if the field is not uniform this block will not work
  
  const double kGev=1.0e9;
  const double kLightVelocity=2.99792458e8;
  Double_t b   = dynamic_cast<const EXKalDetector &>
    (dynamic_cast<EXMeasLayer *>
     (fCradle->At(0))->GetParent(kFALSE)).GetBfield();
  Double_t pt  = (pRho/100.)*kLightVelocity*(b/10.) / kGev;
  
  P_hel = fabs(pt)/sin(pTheta);

#ifdef _EXKalRTPCDebug_
  //just for debug
  if(_EXKalRTPCDebug_>=2) {
    cout<<"  global Helix:  pt="<<pt
	<<"  Rho="<<pRho<<", A="<<pA<<", B="<<pB<<endl;
    cout<<"  P_hel="<<fabs(pt)/sin(pTheta)<<", Phi_hel="<<Phi_hel*57.3
	<<"deg, Theta_hel="<<Theta_hel*57.3<<"deg,  Z_hel="<<Z_hel
	<<", fi0="<<fi0*57.3<<endl;
  }
#endif

  return pChi2;
}


//Do global helix fit to get initial parameter for Kalman Filter
//IterDirection=true is farward, otherwise backward
THelicalTrack EXKalRTPC::GetIniHelixByGHF(bool IterDirection) 
{
  //the buffer should be in increasing time order
  TIter next(fKalHits_Forward, true);  //forward direction

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
    szPos[npt][0]=xv.X();szPos[npt][1]=xv.Y();szPos[npt][2]=xv.Z();

    npt++;
    if(npt>=200) break;
    hitp = dynamic_cast<EXHit *>(next());
  }

  ////////////////////////////////////////////////////////////////////////

  //do the helix fit and store all results into the tree leaves buffer
  double pRho, pA, pB, pPhi, pTheta, pX0, pY0, pZ0, pDCA, pChi2;
  int fit_to_beamline=1;
  helix_fit(npt,szPos,pRho,pA,pB,pPhi,pTheta,pX0,pY0,pZ0,pDCA,pChi2,
	    fit_to_beamline);

  //store global helix result before my correction 
  Phi_hel_raw = pPhi;
  Theta_hel_raw = pTheta;
  R_hel_raw = pRho;
  A_hel_raw = pA;
  B_hel_raw = pB;
  Z_hel_raw = pZ0;
  X_hel_raw = pX0;
  Y_hel_raw = pY0;
  DCA_hel_raw = pDCA;
  Chi2_hel_raw = pChi2;


  /////////////////////////////////////////////////////////////////////
  //check global helix fit if it return a wrong sign of rho for large curve track

  //Global helix fit might return the wrong sign, expecially for large curve back tracks
  //Once it happens, its theta and z are totally wrong, phi is off by kPi according to definition 
  //here I determine the sign  
  //using dfi_vx2first, since it always less than kPi. 
  //For clock-wise track, dfi_vx2first<0  
  //the next few lines will get the phi angle on circle system, then do a subtraction
  double phi_cir_vx = atan2(0-pB,0-pA);
  double phi_cir_first = atan2(szPos[0][1]-pB,szPos[0][0]-pA);
  double phi_cir_last = atan2(szPos[npt-1][1]-pB,szPos[npt-1][0]-pA);

  //double dfi_first2last = phi_cir_last - phi_cir_first;
  //if(dfi_first2last> kPi) dfi_first2last-=2*kPi;
  //if(dfi_first2last<-kPi) dfi_first2last+=2*kPi;

  double dfi_vx2first = phi_cir_first - phi_cir_vx;
  if(dfi_vx2first> kPi) dfi_vx2first-=2*kPi;
  if(dfi_vx2first<-kPi) dfi_vx2first+=2*kPi;

  double dfi_vx2last = phi_cir_last - phi_cir_vx;
  if(dfi_vx2last> kPi) dfi_vx2last-=2*kPi;
  if(dfi_vx2last<-kPi) dfi_vx2last+=2*kPi;

  int sign = (dfi_vx2first<0) ? 1 : -1;

#ifdef _EXKalRTPCDebug_
  //just for debug
  if(_EXKalRTPCDebug_>=5) {
    cout<<" First_hit=("<<szPos[0][0]<<", "<<szPos[0][1]<<", "<<szPos[0][2]<<"), ";
    cout<<"  Last_hit=("<<szPos[npt-1][0]<<", "<<szPos[npt-1][1]<<", "<<szPos[npt-1][2]<<") \n";
    cout<<"  phi_cir_vx="<<phi_cir_vx*57.3<<"  phi_cir_first="<<phi_cir_first*57.3
	<<"  phi_cir_last="<<phi_cir_last*57.3<<"  dfi_vx2first="<<dfi_vx2first*57.3
	<<"  dfi_vx2last="<<dfi_vx2last*57.3<<endl;
  } 
  if(_EXKalRTPCDebug_>=4) {
    cout<<"  npt="<<npt<<",  dz_span="<<szPos[npt-1][2]-szPos[0][2]
	<<"cm,  chi2="<<pChi2<<endl;
  }
#endif

  ////////////////////

  //make correction for global helix fit result if needed
  if (sign*pRho<0) {
#ifdef _EXKalRTPCDebug_
    cout<<"***Warning: global helix fit return wrong sign! Correct it back! \n";
    if(_EXKalRTPCDebug_>=2) {
      cout<<"***Before correction: Rho="<<setw(8)<<pRho<<", Phi="<<setw(8)<<pPhi*57.3
	  <<"deg, Theta="<<setw(8)<<pTheta*57.3<<"deg, Z="<<pZ0<<"cm \n";
    }
#endif
    pRho *= -1.0;
    pPhi+=kPi;
    if(pPhi> kPi) pPhi-=2*kPi;
    if(pPhi<-kPi) pPhi+=2*kPi;

    //Corr theta and z
    CorrHelixThetaZ(npt,szPos,pRho,pA,pB,pTheta,pZ0);
    
#ifdef _EXKalRTPCDebug_
    if(_EXKalRTPCDebug_>=2) {
      cout<<"*** After correction: Rho="<<setw(8)<<pRho<<", Phi="<<setw(8)<<pPhi*57.3
	  <<"deg, Theta="<<setw(8)<<pTheta*57.3<<"deg, Z="<<pZ0<<"cm \n";
    }
#endif

  } else {
    //For those events that the sign of rho is right but 
#ifdef _EXKalRTPCDebug_
    //sometimes the global helix return wrong theta, 
    //want to find out and correct it back
    double ppTheta=pTheta, ppZ0=pZ0;
    CorrHelixThetaZ(npt,szPos,pRho,pA,pB,ppTheta,ppZ0);
    if((ppTheta-kPi/2)*(pTheta-kPi/2)<0) {
      pTheta=ppTheta; pZ0=ppZ0;
      if(_EXKalRTPCDebug_>=2) {
	cout<<"***Before CorrHelixThetaZ: Theta="<<setw(8)<<ppTheta*57.3<<"deg, Z="<<ppZ0<<"cm \n";
	cout<<"*** After CorrHelixThetaZ: Theta="<<setw(8)<<ppTheta*57.3<<"deg, Z="<<ppZ0<<"cm \n";
      }
    }
#endif
  }


  //corr phi and r
  //This part is missing right now because I found that GHF is already the 
  //best one in hand, CircleFitter_LM is worse.
  //you are welcome to develop one!

  //////////////////////////////////////
  //now Fit the circle using Levenberg-Marquardt method
  //////////////////////////////////////
  //sometimes this routine return very large values
  //gLMFitter.DoFit(npt,szPos,pA,pB,pRho);
  

  /////////////////////////////////////////////////////////////////////
  //store global helix result after my correction 
  Phi_hel = pPhi;
  Theta_hel = pTheta;
  R_hel = pRho;
  A_hel = pA;
  B_hel = pB;
  Z_hel = pZ0;
  X_hel = pX0;
  Y_hel = pY0;
  DCA_hel = pDCA;
  Chi2_hel = pChi2;

  
  //////////////////////////////////////
  //now create a THelixTrack
  //////////////////////////////////////

  //helix center phi angle in the hall
  double Phi_c = (pRho>0.) ? pPhi-kPi/2 : pPhi+kPi/2;
  //do not use the next line, 
  //if A_rec or B_rec are close to zero, it will give a wrong sign
  //Phi_c = atan2(B_rec, A_rec); 
  if(Phi_c> kPi) Phi_c-=2*kPi;
  if(Phi_c<-kPi) Phi_c+=2*kPi;

  const double kGev=1.0e9;
  const double kLightVelocity=2.99792458e8;
  Double_t b   = dynamic_cast<const EXKalDetector &>
    (dynamic_cast<EXMeasLayer *>
     (fCradle->At(0))->GetParent(kFALSE)).GetBfield();

  Double_t dr  = pDCA;
  Double_t pt  = (pRho/100.)*kLightVelocity*(b/10.) / kGev;
  Double_t cpa = 1. / pt;
  Double_t fi0 = (pRho>0)?Phi_c:Phi_c-kPi;   //the fi0 of vertex
  if(fi0> kPi) fi0-=2*kPi;
  if(fi0<-kPi) fi0+=2*kPi;
  Double_t dz  = 0;
  Double_t cs  = cos(pTheta);
  Double_t tnl = cs / TMath::Sqrt((1-cs)*(1+cs)); 
  Double_t x0  = 0;  // can also use pX0;
  Double_t y0  = 0;  // can also use pY0;
  Double_t z0  = pZ0;


  //IterDirection=true is farward, return helix state at 1st point
  //otherwise backward, return helix state at last point
  double fi0_last = fi0+dfi_vx2last;
  if(fi0_last> kPi) fi0_last-=2*kPi;
  if(fi0_last<-kPi) fi0_last+=2*kPi;

  double fi0_1st = fi0+dfi_vx2first;
  if(fi0_1st> kPi) fi0_1st-=2*kPi;
  if(fi0_1st<-kPi) fi0_1st+=2*kPi;

  //this helix track is not ideal, but only fi0, tnl and cpa will be used
  //by the caller, these 3 are ok. 
  //Fix me:
  //The helix center of this helix provide wrong A and B. (these 2 values are not used) 
  //Here is one example: 
  /*  
  //Circle Event:  pt=-0.1  Rho=-6.66667, A=6.30317, B=-2.17129, phi_c=341.018deg, fi0=161.004deg 
  //P0_p=0.11547, Phi0_p=-109.016deg, Theta0_p=120.009deg  Z0=0
  +---------------------------------------------------------------------+
  title=GHF(backward): aTrack:
  Kappa=-10.1061, tanLambda=-0.62661
  Rho=-6.60124, A=4.74389, B=4.59041, phi_c=44.0612deg, fi0=-135.952deg 
  pt=0.0989501  P_p=0.116771, Phi_p=-45.9454deg, Theta_p=122.081
  +---------------------------------------------------------------------+
  +---------------------------------------------------------------------+
  title=KalRTPC(backward): 3-point helix at last point:
  Kappa=-10.1552, tanLambda=-0.670794
  Rho=-6.56935, A=6.21079, B=-2.23221, phi_c=43.5931deg, fi0=-136.42deg 
  pt=0.098472  P_p=0.118575, Phi_p=-46.4135deg, Theta_p=123.863
  +---------------------------------------------------------------------+
  */

  double ret_fi0 = (IterDirection==kIterBackward)?fi0_last:fi0_1st;
  THelicalTrack aTrack(dr,ret_fi0,cpa,dz,tnl,x0,y0,z0,b);
  
  //calculate P_hel
  P_hel = fabs(pt)/sin(pTheta);

#ifdef _EXKalRTPCDebug_
  //just for debug
  if(_EXKalRTPCDebug_>=2) {
    cout<<"  global Helix:  pt="<<pt
	<<"  Rho="<<pRho<<", A="<<pA<<", B="<<pB
	<<", phi_c="<<Phi_c*57.3<<"deg,"
	<<" fi0_1st="<<fi0_1st*57.3<<"deg "
	<<" fi0_last="<<fi0_last*57.3<<"deg "<<endl;
    cout<<"  P_hel="<<fabs(pt)/sin(pTheta)<<", Phi_hel="<<Phi_hel*57.3
	<<"deg, Theta_hel="<<Theta_hel*57.3<<"deg,  Z_hel="<<Z_hel
	<<", fi0="<<fi0*57.3<<endl;
  }
#endif

  return aTrack;
}


//prepare a track from xyz array in mm, make sure hits are in increasing time order   
bool EXKalRTPC::PrepareATrack_mm(double *x_mm, double *y_mm,double *z_mm, 
				 int npt, bool smearing, bool bIncludeCurveBackHits)
{
  fEventGen->MakeHitsFromTraj_mm(x_mm,y_mm,z_mm,npt,smearing,bIncludeCurveBackHits);
  return true;
}

//prepare a track from xyz array, make sure radii are in increasing order
bool EXKalRTPC::PrepareATrack(double *x, double *y,double *z, int npt, bool smearing,
			      bool bIncludeCurveBackHits)
{
  fEventGen->MakeHitsFromTraj(x,y,z,npt,smearing,bIncludeCurveBackHits);
  return true;
}

//Use the event generator to generate a track
bool EXKalRTPC::PrepareATrack(int job, double pt_min, double pt_max, double costh_min,
			      double costh_max, double z_min, double z_max, bool bIncludeCurveBackHits)
{
  //sometimes the user will provide too small pt_min, which will cause some
  //events have no hits, I have to avoid this situation
  //Here I set a maximum number of throw as 1000, if it fails then quit
  int pCounter=0,pMaxThrow=1000;
  
  if(job==0) {
    while(pCounter<pMaxThrow) {
      THelicalTrack hel = fEventGen->GenerateHelix(pt_min,pt_max,costh_min,costh_max,z_min,z_max);
#ifdef _EXKalRTPCDebug_
      if(Global_Debug_Level >= 4) {
	EXEventGen::PrintHelix(&hel, "GenerateHelix(backward): thrown helix at vertex:");
      }
#endif
      fEventGen->Swim(hel,bIncludeCurveBackHits,kMpr);
#ifdef _EXKalRTPCDebug_
      if(Global_Debug_Level >= 4) {
	EXEventGen::PrintHelix(&hel, "GenerateHelix(backward): swim thrown helix to last point:");
      }
#endif
      if(fKalHits->GetEntriesFast()>=MinHit) break;
      else fKalHits->Delete();
      pCounter++;
    }
    if(pCounter>=pMaxThrow) return false;
  }
  else {
    while(pCounter<pMaxThrow) {
      fEventGen->GenerateCircle(pt_min,pt_max,costh_min,costh_max,z_min,z_max,bIncludeCurveBackHits);
      if(fKalHits->GetEntriesFast()>=MinHit) break;
      else fKalHits->Delete();
      pCounter++;
    }
    if(pCounter>=pMaxThrow) return false;
  }
  return true;
}


// ============================================================
//  Do Kalman Filter
//  The hit buffer will be filled outside
//  Apply 2nd iteration will only help those track that curve back 
// ============================================================
int EXKalRTPC::DoFitAndFilter(bool bApply2Iter)
{
  // ============================================================
  // Get the initial helix track
  // ============================================================
  //  Get initial helix track by 3-point calculation or global helix fit.
  //  If apply two iterations, do the 1st iteration in backward direction.
  //  Use 1st iteration result to start the 2nd iteration
  THelicalTrack helstart;
  TKalMatrix C_start(kSdim,kSdim);


  if(bApply2Iter) {
#ifdef _EXKalRTPCDebug_
    //this part just for debug, show the 1st iter forward fitting result
    //to compare with 2nd iter
    if(Global_Debug_Level >= 4) {
      FitForward4InitHelix(helstart,C_start);
      EXEventGen::PrintHelix(&helstart, "FitForward result: helix at last point:");
      //cout<<"\nFitForward4InitHelix covMat:"<<endl;
      //C_start.DebugPrint(25);
    }
#endif
    //Fit backward and smooth it back to the most outside point
    //This has been confirmed to be the best choice
    FitBackward4InitHelix(helstart,C_start);
#ifdef _EXKalRTPCDebug_
    if(Global_Debug_Level >= 4) {
      EXEventGen::PrintHelix(&helstart, "Fitbackward and smooth back result: helix at last point:");
      //cout<<"\nFitBackward4InitHelix covMat(smoothed back):"<<endl;
      //C_start.DebugPrint(25);
    }
#endif
  } else {
    THelicalTrack hel_3point = this->GetIniHelixBy3Pts(kDir);
    THelicalTrack hel_global = this->GetIniHelixByGHF(kDir);
#ifdef _EXKalRTPCDebug_
    if(Global_Debug_Level >= 4) {
      EXEventGen::PrintHelix(&hel_3point, "KalRTPC(backward): 3-point helix at last point:");
      EXEventGen::PrintHelix(&hel_global, "KalRTPC(backward): global helix at last point:");
    }
#endif
    //for some curve back tracks, global helix fit fails,
    //use 3-point result if it happens
    helstart = (hel_global.GetRho()*hel_3point.GetRho()>0) ? hel_global : hel_3point;
  }
  
  // ---------------------------
  //store the initial parameter: Rho, Phi0, Theta
  // ---------------------------
  //helstart.GetPtoR() will return alpha.   alpha/kappa=rho
  rho_kal_ini  = helstart.GetRho();
  tnl_kal_ini  = helstart.GetTanLambda();
  phi0_kal_ini = helstart.GetPhi0();

  //debug: I want to use exact last point
  //rho_kal_ini  = fEventGen->Rho_last;
  //tnl_kal_ini  = fEventGen->TanLambda_last;
  //phi0_kal_ini = fEventGen->Phi0_last;


  // ---------------------------
  //  Create a dummy site: sited
  // ---------------------------
  Int_t i1 = (kDir == kIterBackward) ? fKalHits_Forward->GetEntries()-1 : 0;

  EXHit hitd = *dynamic_cast<EXHit *>(fKalHits_Forward->At(i1));
  hitd(0,1) = 1.e6;   // give a huge error to d
  hitd(1,1) = 1.e6;   // give a huge error to z

  TKalTrackSite &sited = *new TKalTrackSite(hitd);
  sited.SetOwner();   // site owns states

  // ---------------------------
  //  Set dummy state to sited
  // ---------------------------

  static TKalMatrix svd(kSdim,1);
  svd(0,0) = 0.;
  svd(1,0) = phi0_kal_ini;                   //helstart.GetPhi0();
  svd(2,0) = helstart.GetPtoR()/rho_kal_ini; //helstart.GetKappa();
  svd(3,0) = 0.;
  svd(4,0) = tnl_kal_ini;                    //helstart.GetTanLambda();
  if (kSdim == 6) svd(5,0) = 0.;


  static TKalMatrix C(kSdim,kSdim);
  if(bApply2Iter){
    //for (Int_t i=0; i<kSdim; i++) C(i,i) = C_start(i,i);
    C = C_start;
  }
  else
    {
      for (Int_t i=0; i<kSdim; i++) {
	//C(i,i) = 0.05;         // dummy error matrix
	C(i,i) = fCovMElement;   // dummy error matrix
      }
    }

  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

  // ---------------------------
  //  Add sited to the kaltrack
  // ---------------------------
  if(fKalTrack) delete fKalTrack;
  fKalTrack = new TKalTrack(); 
  //fKalTrack->Delete();     // a track is a kal system, reset this buffer 
  //It turns out that Delete() will only remove the elements from this objarray
  //but not reset chi2 and other properties, so I have to create a new instance  
  fKalTrack->SetMass(kMpr);
  fKalTrack->SetOwner();   // kaltrack owns sites
  fKalTrack->Add(&sited);  // add the dummy site to the track

  // ---------------------------
  //  Prepare hit iterrator
  // ---------------------------

  TIter next(fKalHits_Forward, kDir);   // come in to IP
  
  // ============================================================
  // Start Kalman Filter
  // ============================================================
  int npt = 0;
  EXHit *hitp = dynamic_cast<EXHit *>(next());
  while (hitp) {     // loop over hits
    step_status[npt]=1;   // store the status into root tree

#ifdef _EXKalRTPCDebug_    
    //Just for debug
    const EXMeasLayer &ml = dynamic_cast<const EXMeasLayer &>(hitp->GetMeasLayer());
    TVector3 xraw = hitp->GetRawXv();
    TVector3 xv = ml.HitToXv(*hitp);

    if(Global_Debug_Level >= 5) {
      cerr << "MeasLayer "<<setw(2)<<ml.GetIndex()
	   <<": R="<<setw(6)<<ml.GetR()<<": ";
      cerr << "xraw=("<<setw(8)<<xraw.X()<<",  "<<setw(8)<<xraw.Y()
	   <<", "<<setw(8)<<xraw.Z()<<"): ";
      cerr << "xv=("<<setw(8)<<xv.X()<<",  "<<setw(8)<<xv.Y()
	   <<", "<<setw(8)<<xv.Z()<<"): \n";
    }
#endif

    TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
    if (!fKalTrack->AddAndFilter(site)) {             // add and filter this site
      delete &site;                                   // delete this site, if failed
      step_status[npt]=0;                             // store the status into root tree 
      
#ifdef _EXKalRTPCDebug_
      cerr << "Fitter: site "<<npt<<" discarded: "
	   <<" xv=("<< xv.X()<<", "<<xv.Y()<<", "<<xv.Z()<<"), Phi="
	   <<xv.Phi()*57.3 <<"deg, S="<<xv.Perp()<<"cm"<< endl;
      Pause4Debug();
#endif
    } else {
#ifdef _EXKalRTPCDebug_
      //get filtered state then convert it to THelicalTrack
      //both ways to get state will work
      //TVKalState &state_fil1 = fKalTrack->GetState(TVKalSite::kFiltered);
      //TVKalState *state_fil2 = &fKalTrack->GetState(TVKalSite::kFiltered);

      if(Global_Debug_Level >= 7) {
	TVKalState *state_fil = (TVKalState*) &(site.GetCurState());
	THelicalTrack hel_fil = (dynamic_cast<TKalTrackState *>(state_fil))->GetHelix();
	TVector3 x_fil=hel_fil.CalcXAt(0.0);

	TVKalState *state_exp = &site.GetState(TVKalSite::kPredicted);
	THelicalTrack hel_exp = (dynamic_cast<TKalTrackState *>(state_exp))->GetHelix();
	TVector3 x_exp=hel_exp.CalcXAt(0.0);

	//debug: comparing the expected and filtered points
	cerr << "Event "<<"XXXX"<<": site "<<npt<<endl;
	cerr << "\t xraw =("<< xraw.X()<<",  "<<xraw.Y()<<", "<<xraw.Z()<<"): \n";
	cerr << "\t xv   =("<< xv.X()<<",  "<<xv.Y()<<", "<<xv.Z()<<"): \n";
	cerr << "\t x_exp=("<< x_exp.X()<<",  "<<x_exp.Y()<<", "<<x_exp.Z()<<"): \n";
	cerr << "\t x_fil=("<< x_fil.X()<<",  "<<x_fil.Y()<<", "<<x_fil.Z()<<"): \n";
      }
#endif
    }
    npt++;
    hitp = dynamic_cast<EXHit *>(next());
  }
  //this line will smooth the cursor back to the 1st site, which means GetCurSite()
  //will return the 1st site, not the last site
  //Please also note that 1st point is located at the cathode, the last point at GEM1
  //fKalTrack->SmoothBackTo(1);                          // smooth back.


#ifdef _EXKalRTPCDebug_
  if(Global_Debug_Level >= 3) {
    cout<<"Number of sites in kaltrack = "<<fKalTrack->GetEntries()<<endl;
  }
#endif
  // ============================================================
  //  Monitor Fit Result
  // ============================================================

  TVKalState *theLastState = (TVKalState*) &(fKalTrack->GetCurSite().GetCurState());
  ReconVertex(*theLastState, P_rec, Pt_rec, Pz_rec, Theta_rec, Phi_rec,
	      X_rec, Y_rec, Z_rec, R_rec, A_rec, B_rec );

#ifdef _EXKalRTPCDebug_
  if(Global_Debug_Level >= 4) {
    cout<<"\nFinal covMat:"<<endl;
    theLastState->GetCovMat().DebugPrint(25);
  }
#endif
  
  // ============================================================
  //  Monitor Fit Result
  // ============================================================
  //TMath::Prob(chi2, ndf) compute the probability for a certain Chi-squared (chi2)
  //and number of degrees of freedom (ndf).
  //Calculations are based on the incomplete gamma function P(a,x),
  //where a=ndf/2 and x=chi2/2.
  //P(a,x) represents the probability that the observed Chi-squared
  //for a correct model should be less than the value chi2.
  NDF  = fKalTrack->GetNDF();
  Chi2 = fKalTrack->GetChi2();
  //cl   = TMath::Prob(chi2, ndf);

  return 1;
}

//Do Kalman Filter with given x,y,x array
//Note that these hits should be in increasing order of time
//this routine only apply 1 iteration KF fit
int EXKalRTPC::DoFitAndFilter(double *x_cm, double *y_cm, double *z_cm, int n,
			      bool bIncludeCurveBackHits)
{
  if(n<5) return 0;

  // ============================================================
  // Fill the hit buffer with given array
  // ============================================================
  Reset();

  PrepareATrack(x_cm, y_cm, z_cm, n,bIncludeCurveBackHits);
  JudgeFor2ndIteration(bIncludeCurveBackHits);

  // ============================================================
  //  Do Kalman Filter
  // ============================================================

  // ---------------------------
  // Create initial helix
  // ---------------------------
  THelicalTrack hel_3point = this->GetIniHelixBy3Pts(kDir);
  THelicalTrack helstart = this->GetIniHelixByGHF(kDir);
#ifdef _EXKalRTPCDebug_
  if(Global_Debug_Level >= 4) {
    EXEventGen::PrintHelix(&hel_3point, "3-point helix at last point:");
    EXEventGen::PrintHelix(&helstart, "global helix at last point:");
  }
#endif

  //for some curve back tracks, global helix fit fail
  if(helstart.GetRho()*hel_3point.GetRho()<0) helstart=hel_3point;

  // ---------------------------
  //  Create a dummy site: sited
  // ---------------------------

  Int_t i1 = (kDir == kIterBackward) ? fKalHits_Forward->GetEntries() - 1 : 0;
  EXHit hitd = *dynamic_cast<EXHit *>(fKalHits_Forward->At(i1));
  hitd(0,1) = 1.e6;   // give a huge error to d
  hitd(1,1) = 1.e6;   // give a huge error to z

  TKalTrackSite &sited = *new TKalTrackSite(hitd);
  sited.SetOwner();   // site owns states


  // ---------------------------
  //  Set dummy state to sited
  // ---------------------------

  static TKalMatrix svd(kSdim,1);
  svd(0,0) = 0.;
  svd(1,0) = helstart.GetPhi0();
  svd(2,0) = helstart.GetKappa();
  svd(3,0) = 0.;
  svd(4,0) = helstart.GetTanLambda();
  if (kSdim == 6) svd(5,0) = 0.;

  static TKalMatrix C(kSdim,kSdim);
  for (Int_t i=0; i<kSdim; i++) {
    C(i,i) = fCovMElement;   // dummy error matrix
  }

  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

  // ---------------------------
  //  Add sited to the kaltrack
  // ---------------------------

  // a track is a kal system, it is also an objarray
  if(fKalTrack) delete fKalTrack;
  fKalTrack = new TKalTrack(); 
  //fKalTrack->Delete();     // a track is a kal system, reset this buffer 
  //It turns out that Delete() will only remove the elements from this objarray
  //but not reset chi2 and other properties, so I have to create a new instance  
  fKalTrack->SetMass(kMpr);
  fKalTrack->SetOwner();   // fKalTrack-> owns sites
  fKalTrack->Add(&sited);  // add the dummy site to the track

  // ---------------------------
  //  Prepare hit iterrator
  // ---------------------------

  TIter next(fKalHits_Forward, kDir);   // come in to IP

  // ---------------------------
  //  Start Kalman Filter
  // ---------------------------

  int npt = 0;
  EXHit *hitp = dynamic_cast<EXHit *>(next());
  while (hitp) {     // loop over hits
    step_status[npt]=1;
    TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
    if (!fKalTrack->AddAndFilter(site)) {               // add and filter this site
      delete &site;                                   // delete this site, if failed
      step_status[npt]=0;

#ifdef _EXKalRTPCDebug_
      const EXMeasLayer &ml = dynamic_cast<const EXMeasLayer &>(hitp->GetMeasLayer());
      TVector3 xv = ml.HitToXv(*hitp);
      cerr << "Fitter: site "<<npt<<" discarded: "
	   << " xv=("<< xv.X()<<",  "<<xv.Y()<<", "<<xv.Z()<<")"<< endl;
#endif
    }
    npt++;
    hitp = dynamic_cast<EXHit *>(next());
  }
  //this line will smooth the cursor back to the 1st site, which means GetCurSite()
  //will return the 1st site, not the last site
  //Please also note that 1st point is located at the cathode, the last point at GEM1
  //fKalTrack->SmoothBackTo(1);                          // smooth back.

  // ============================================================
  //  reconstruct to the vertrex
  // ============================================================

  TVKalState *theLastState = (TVKalState*) &(fKalTrack->GetCurSite().GetCurState());
  ReconVertex(*theLastState, P_rec, Pt_rec, Pz_rec, Theta_rec, Phi_rec,
	      X_rec, Y_rec, Z_rec, R_rec, A_rec, B_rec );
  
  NDF  = fKalTrack->GetNDF();
  Chi2 = fKalTrack->GetChi2();

  return 1;
}

//this routine is to fit the existing hits in forward direction to get
//a helix, then use this helix as input to the 2nd iteration
//note: hits are all stored in fKalHits and fKalHits_Forward
void EXKalRTPC::FitForward4InitHelix(THelicalTrack &Hel_last,TKalMatrix &C_last)
{
  // ============================================================
  //  Do Kalman Filter in forward direction to get a helix
  // ============================================================

  // ---------------------------
  // Create initial helix
  // ---------------------------
  THelicalTrack hel_3point = this->GetIniHelixBy3Pts(kIterForward);
  THelicalTrack helstart = this->GetIniHelixByGHF(kIterForward);
#ifdef _EXKalRTPCDebug_
  if(Global_Debug_Level >= 4) {
    EXEventGen::PrintHelix(&hel_3point, "FitForward4InitHelix(): 3-point helix at 1st point:");
    EXEventGen::PrintHelix(&helstart, "FitForward4InitHelix(): global helix at 1st point:");
  }
#endif

  //for some curve back tracks, global helix fit fail
  if(helstart.GetRho()*hel_3point.GetRho()<0) helstart=hel_3point;

  // ---------------------------
  //  Create a dummy site: sited
  // ---------------------------
  //fitting forward
  Int_t i1 = 0;
  EXHit hitd = *dynamic_cast<EXHit *>(fKalHits_Forward->At(i1));
  hitd(0,1) = 1.e6;   // give a huge error to d
  hitd(1,1) = 1.e6;   // give a huge error to z

  TKalTrackSite &sited = *new TKalTrackSite(hitd);
  sited.SetOwner();   // site owns states


  // ---------------------------
  //  Set dummy state to sited
  // ---------------------------

  static TKalMatrix svd(kSdim,1);
  svd(0,0) = 0.;
  svd(1,0) = helstart.GetPhi0();
  svd(2,0) = helstart.GetKappa();
  svd(3,0) = 0.;
  svd(4,0) = helstart.GetTanLambda();
  if (kSdim == 6) svd(5,0) = 0.;

  static TKalMatrix C(kSdim,kSdim);
  for (Int_t i=0; i<kSdim; i++) {
    C(i,i) = fCovMElement;   // dummy error matrix
  }

  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

  // ---------------------------
  //  Add sited to the kaltrack
  // ---------------------------

  TKalTrack kaltrack;    // a track is a kal system
  kaltrack.SetMass(kMpr);
  kaltrack.SetOwner();   // kaltrack owns sites
  kaltrack.Add(&sited);  // add the dummy site to the track

  // ---------------------------
  //  Prepare hit iterrator
  // ---------------------------

  TIter nextforward(fKalHits_Forward, kIterForward);   // looping forward

  // ---------------------------
  //  Start Kalman Filter
  // ---------------------------

  EXHit *hitp = dynamic_cast<EXHit *>(nextforward());
  while (hitp) {     // loop over hits
    TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
    if (!kaltrack.AddAndFilter(site)) {               // add and filter this site
      delete &site;                                   // delete this site, if failed
    }
    hitp = dynamic_cast<EXHit *>(nextforward());
  }
  //this line will smooth the cursor back to the 1st site, which means GetCurSite()
  //will return the 1st site, not the last site
  //Please also note that 1st point is located at the GEM1, the last point at GEM1
  //kaltrack.SmoothBackTo(1);                           // smooth back.

  // ============================================================
  //  Get the helix at the last site then return it back
  // ============================================================

  TKalTrackState *theLastState = dynamic_cast<TKalTrackState*> (&(kaltrack.GetCurSite().GetCurState()));
  Hel_last = theLastState->GetHelix();
  C_last = theLastState->GetCovMat();
#ifdef _EXKalRTPCDebug_
  if(Global_Debug_Level >= 4) {
    cout<<"\nFitForward4InitHelix covMat:"<<endl;
    theLastState->GetCovMat().DebugPrint(25);
  }
#endif

}


//this routine is to fit the existing hits in backward direction 
//then smooth it back to the 1st site to get a helix.
//This helix will then be used as input to the 2nd iteration
//note: hits are all stored in fKalHits and fKalHits_Forward
void EXKalRTPC::FitBackward4InitHelix(THelicalTrack &Hel_1st,TKalMatrix &C_1st)
{
  // ============================================================
  //  Do Kalman Filter in forward direction to get a helix
  // ============================================================

  // ---------------------------
  // Create initial helix
  // ---------------------------
  THelicalTrack hel_3point = this->GetIniHelixBy3Pts(kIterBackward);
  THelicalTrack helstart = this->GetIniHelixByGHF(kIterBackward);
#ifdef _EXKalRTPCDebug_
  if(Global_Debug_Level >= 4) {
    EXEventGen::PrintHelix(&hel_3point, "FitBackward4InitHelix(): 3-point helix at last point:");
    EXEventGen::PrintHelix(&helstart, "FitBackward4InitHelix(): global helix at last point:");
  }
#endif

  //for some curve back tracks, global helix fit fail
  if(helstart.GetRho()*hel_3point.GetRho()<0) helstart=hel_3point;

  // ---------------------------
  //  Create a dummy site: sited
  // ---------------------------
  //fitting backward
  Int_t i1 = fKalHits_Forward->GetEntriesFast() - 1;
  EXHit hitd = *dynamic_cast<EXHit *>(fKalHits_Forward->At(i1));
  hitd(0,1) = 1.e6;   // give a huge error to d
  hitd(1,1) = 1.e6;   // give a huge error to z

  TKalTrackSite &sited = *new TKalTrackSite(hitd);
  sited.SetOwner();   // site owns states


  // ---------------------------
  //  Set dummy state to sited
  // ---------------------------

  static TKalMatrix svd(kSdim,1);
  svd(0,0) = 0.;
  svd(1,0) = helstart.GetPhi0();
  svd(2,0) = helstart.GetKappa();
  svd(3,0) = 0.;
  svd(4,0) = helstart.GetTanLambda();
  if (kSdim == 6) svd(5,0) = 0.;

  static TKalMatrix C(kSdim,kSdim);
  for (Int_t i=0; i<kSdim; i++) {
    C(i,i) = fCovMElement;   // dummy error matrix
  }

  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

  // ---------------------------
  //  Add sited to the kaltrack
  // ---------------------------

  TKalTrack kaltrack;    // a track is a kal system
  kaltrack.SetMass(kMpr);
  kaltrack.SetOwner();   // kaltrack owns sites
  kaltrack.Add(&sited);  // add the dummy site to the track

  // ---------------------------
  //  Prepare hit iterrator
  // ---------------------------

  TIter nextbackward(fKalHits_Forward, kIterBackward);   // looping backward

  // ---------------------------
  //  Start Kalman Filter
  // ---------------------------

  EXHit *hitp = dynamic_cast<EXHit *>(nextbackward());
  while (hitp) {     // loop over hits
    TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
    if (!kaltrack.AddAndFilter(site)) {               // add and filter this site
      delete &site;                                   // delete this site, if failed
    }
    hitp = dynamic_cast<EXHit *>(nextbackward());
  }
  //this line will smooth the cursor back to the 1st site, which means GetCurSite()
  //will return the 1st site, not the last site
  //Please also note that 1st point is located at the cathode, the last point at GEM1
  kaltrack.SmoothBackTo(1);                           // smooth back.

  // ============================================================
  //  Get the last site then return it back
  // ============================================================

  TKalTrackState *theLastState = dynamic_cast<TKalTrackState*> (&(kaltrack.GetCurSite().GetCurState()));
  Hel_1st = theLastState->GetHelix();
  C_1st = theLastState->GetCovMat();
#ifdef _EXKalRTPCDebug_
  if(Global_Debug_Level >= 4) {
    cout<<"\nFitBackward4InitHelix covMat:"<<endl;
    theLastState->GetCovMat().DebugPrint(25);
  }
#endif

}

//this is just an example to illustrate how to use this code
void EXKalRTPC::Example(int job, int nevents, double pt_min, double pt_max, double costh_min, 
			double costh_max, double z_min, double z_max)
{
	
  // setup input|output root tree here 	
  //BeginOfRun();
	
  // ===================================================================
  //  set covariant matrix element if needed
  // ===================================================================
  // if you think 0.05 is not good, change it as you like	
  //this->SetCovMElement(0.05);
  
  
  // ===================================================================
  //  Event loop
  // ===================================================================
  for (Int_t eventno = 0; eventno < nevents; eventno++) { 

    // ============================================================
    //  Reset the buffer
    // ============================================================
    this->Reset();

    // ============================================================
    //  Generate a partcle and Swim the particle in fDetector
    // ============================================================
    bool bIncludeCurveBackHits=true;
    if(job == 0 || job == 2)
      PrepareATrack(job, pt_min, pt_max, costh_min, costh_max, z_min, z_max, bIncludeCurveBackHits);
    else if(job == 1) {
      //note that you have to provide x,y,z array
      //PrepareATrack(x,y,z,npt,false);
    }
    else {
      cout<<"This job("<<job<<") is not yet supported! I quit..."<<endl;
      exit(-3);
    }

    //Remove backward hits also judge whether or not need 2nd iteration
    bool bRemoveBackwardHits=true;
    bool bNeed2Iter = this->JudgeFor2ndIteration(bRemoveBackwardHits);

    //make sure there are at leat 5 hits in the buffer
    if(fKalHits_Forward->GetEntriesFast()<MinHit)  continue;

    // ============================================================
    //  Do KalmanFilter
    // ============================================================
    
    this->DoFitAndFilter(bNeed2Iter);
    
    //Fill the output root tree
    //Fill_Output_Tree();
  }
  
  //save the output
  //EndOfRun();
}
