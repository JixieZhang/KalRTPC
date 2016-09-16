// ********************************************************************
//
// $Id: BonusHelixFit.cc,v 1.1 2016/03/08 Bonus12 Exp $
// By Jixie Zhang, Modify for RTPC12
// ********************************************************************

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "BonusHelixFit.hh"
//#define HELIXFIT_DEBUG 1
/*
SUBROUTINE RWFTHL(NPT,RF,PF,WFI,ZF,WZF,IOPT,
1                  VV0,EE0,CH2PH,CH2Z,DEL,DELZ)
implicit none
C
C-----------------------------------------------------------------------
C! Fast helix fit
C
C   A generalization of the TFTHEL routine to allow it to be called
C   from a routine that contains any list of x and y values XF,YF for a
C   set of NPT points to be fitted.
C
C   Input:  NPT    /I     Number of 3-D points to be fit
C           XF     /R     Array of X-values of points to be fit
C           YF     /R     Array of Y-values of points to be fit
C           RF     /R     Array of R-values of points to be fit
C           PF     /R     Array of PHI-values of points to be fit
C           WFI    /R     Array of 1/(sig(rphi))**2 for each point
C           ZF     /R     Array of Z-values of points to be fit
C           WZF    /R     Array of 1/(sig(z))**2 for each point
C           IOPT = 0 -> DISTANCE**2=X**2+Y**2 MINIMISED
C                  1 -> WEIGHTED WITH 1/SIMA(R*PHI)**2
C                  2 -> ERROR MATRIX CALCULATED
C                  3 -> 3-DIMENSIONAL ITERATION
C  OUTPUT:   VV0 = 1/R*CHARGE    POS. IF CLOCKWISE
C                  TAN(LAMBDA)  {=DZ/DS}TAN(ANGLE TO X,Y PLANE)
C                  PHI0         {0,2PI} ANGLE TO X-AXIS at R=D0
C                  D0*SIGN      [CM]    MINIMAL DIST. TO Z-AXIS,
C                                       +VE IF AXIS ENCIRCLED
C                  Z0           [CM]    Z POS AT R=D0
C          EE0 = INVERSE OF ERROR MATRIX IN TRIANG. FORM
C          CH2PH = CHI SQUARED = SUM (PHI DEVIATIONS/ERRORS)**2
C          CH2Z  = CHI SQUARED = SUM (Z DEVIATIONS/ERRORS)**2
C          DEL = ??
c          DELZ= ??
C  NOTE: DEGREES OF FREEDOM = 2*NPT-5
C----------------------------------------------------------------
C     BASED ON  SUBROUTINE CIRCLE
C     REFERENCE:  COMPUTER PHYSICS COMMUNICATIONS VOL 33,P329
C
C   AUTHORS:  N. CHERNOV, G. OSOSKOV & M. POPPE
C   Modified by:  Fred Weber, 8 Jun 1989
C
C-----------------------------------------------------------------
*/

#ifdef __cplusplus
extern "C" {
#endif 
  //file  rwfthl.f
  void rwfthl_(int* npt, float* rf, float* pf, float* wfi, float* zf, float* wzf, 
    int* iopt, float* vv0, float* ee0, 
    float* ch2ph, float* ch2z,  /* return values */ 
    float* del, float* delz);

  //file  rwfthc.cc
  void rwfthc(int npt, float* rf, float* pf, float* wfi, float* zf, float* wzf, 
    int iopt, float* vv0, float* ee0, 
    float* ch2ph, float* ch2z,  /* return values */ 
    float* del, float* delz);

  void rwfthl_c(int npt, float* rf, float* pf, float* wfi, float* zf, float* wzf, 
    int iopt, float* vv0, float* ee0, 
    float* ch2ph, float* ch2z,  /* return values */ 
    float* del, float* delz)
  {
    int Helix_Use_C_Code=1;
    if(Helix_Use_C_Code ) 
      rwfthc(  npt,rf,pf,wfi,zf,wzf, iopt,vv0,ee0, ch2ph, ch2z, del,delz);
    else
      rwfthl_(&npt,rf,pf,wfi,zf,wzf,&iopt,vv0,ee0, ch2ph, ch2z, del,delz);
  }

  //if running windows, don't link the fortran code------------------
#ifdef WIN32_NO_G77_SUPPORT

  void rwfthl_(int* npt, float* rf, float* pf, float* wfi, float* zf, float* wzf, 
    int* iopt, float* vv0, float* ee0, 
    float* ch2ph, float* ch2z,  /* return values */ 
    float* del, float* delz)
  {
    return;
  }
#endif 
  //-----------------------------------------------------------------

#ifdef __cplusplus
}
#endif 

//apply corrction to R and Phi
void CorrHelixRPhi(double &Rho, double &Phi)
{
  std::cout<<"\nBefore CorrHelixRPhi():  R="<<Rho<<"  Phi="<<Phi<<std::endl;
  
/*******************************************   
ph_hel-ph0:r_hel          
Minimizer is Linear                                                    
Chi2                      =       3325.1                               
NDf                       =           76                               
p0                        =      1.44928   +/-   0.0193378             
p1                        =     -0.93305   +/-   0.0176701             
p2                        =      0.24896   +/-   0.00640558            
p3                        =   -0.0355206   +/-   0.00119808            
p4                        =   0.00292522   +/-   0.000125733           
p5                        = -0.000139073   +/-   7.4615e-06            
p6                        =  3.53636e-06   +/-   2.33549e-07           
p7                        = -3.71634e-08   +/-   2.99358e-09  
*///****************************************
  
  double R=fabs(Rho);
  if(R<20.0 && R>2.6) {
    double Para_dPhiVsR[] = {1.44928, -0.93305, 0.24896, -0.0355206,
      0.00292522, -0.000139073, 3.53636e-06, -3.71634e-08};
    double dPhi = Para_dPhiVsR[0];
    for(int i=1;i<=7;i++) dPhi += Para_dPhiVsR[i]*pow(R,i);
    Phi = Phi - dPhi;
  }


/******************************************
r_hel-rho_1st:r_hel, only good for 2.6<r<4.2          
Minimizer is Linear
Chi2                      =      635.106
NDf                       =           74
p0                        =     -1083.82   +/-   43.3963
p1                        =         1858   +/-   77.8616
p2                        =     -1343.62   +/-   59.2763
p3                        =      531.657   +/-   24.8203
p4                        =     -124.419   +/-   6.17324
p5                        =      17.2351   +/-   0.912016
p6                        =     -1.30966   +/-   0.0741081
p7                        =     0.042147   +/-   0.00255522       
*///****************************************
  
  if(R<4.2 && R>2.6) {
  double Para_dRVsR[] = { -1083.82, 1858, -1343.62, 531.657, -124.419,
			  17.2351, -1.30966, 0.042147};
    double dR = Para_dRVsR[0];
    for(int i=1;i<=7;i++) dR += Para_dRVsR[i]*pow(R,i);
    R = R - dR;
    Rho = (Rho<0)? -R : R;
    //std::cout<<"\t\t dR="<<dR<<std::endl;
  }

  std::cout<<"After  CorrHelixRPhi():  R="<<Rho<<"  Phi="<<Phi<<std::endl;

}

/*------------------------------------------------------------------------\
//Function name: void helix_fit(int PointNum,double szPos[][3], 
//  double& Rho, double& A, double& B,double& Phi, double& Theta, 
//  double& X0, double& Y0,double& Z0, double &DCA, double &chi2);
//
//  Calculate the raidus of the trajectory
//
//Input parameters:
//  PointNum:       number of x-y points
//  szPos[any number>3][3]:  xyz array, 
//OutPut: 
//  Rho:       Radius, positive is clockwise, 
//  A,B:       helix center X,Y position
//  Theta,Phi: theta and phi angle at the initial step
//  X0,Y0,Z0:  vertex position (x,y,z) for inital step
//  DCA:       distance of closest approach to beamline
//  Chi2:      ch2ph+ch2z/(npt-5) 
\------------------------------------------------------------------------*/
void helix_fit(int PointNum,double szPos[][3], double& Rho, double& A, double& B,
  double& Phi, double& Theta, double& X0, double& Y0,double& Z0,  
  double &DCA, double& Chi2,int fit_track_to_beamline)
{
  const int MAX_HITS_ON_CHAIN = 255;
  const float PI = acos(0.0)*2.0;

  int jj;
  float my_phi;
  float rf[MAX_HITS_ON_CHAIN];
  float pf[MAX_HITS_ON_CHAIN];
  float wfi[MAX_HITS_ON_CHAIN];
  float zf[MAX_HITS_ON_CHAIN];
  float wzf[MAX_HITS_ON_CHAIN];
  int iopt;
  int npt;
  float vv0[5];
  float ee0[15];
  float ch2ph;
  float ch2z;
  float del[MAX_HITS_ON_CHAIN];
  float delz[MAX_HITS_ON_CHAIN];

  float phi0;
  if(PointNum>MAX_HITS_ON_CHAIN) PointNum=MAX_HITS_ON_CHAIN-1;

  ///////////////////////////////////////////
  //By Jixie:  this helix fit does not work if the track curve back
  //Here I simply remove those points when they start to curve back.
  //I only keep those points that all the way reach maxR
  
  npt=0;
  double tmpR=0.0,tmpRmax=0.0;
  for (jj=0; jj<PointNum; jj++)
  { 
    tmpR = sqrt(pow(szPos[jj][0],2)+pow(szPos[jj][1],2));
    if (tmpR>tmpRmax) {
      tmpRmax=tmpR;
      npt++;
    }
  }
//#ifdef HELIXFIT_DEBUG
  if(npt<PointNum)
    std::cout<<"helix_fit(): input PointNum="<<PointNum
      <<", used PointNum="<<npt<<std::endl;
//#endif

  for (jj=0; jj<npt; jj++)
  {// r,phi,z coordinate
#ifdef HELIXFIT_DEBUG
    printf("point %3d: X=%.2f  Y=%.2f  Z=%.2f\n",jj+1,szPos[jj][0],szPos[jj][1],szPos[jj][2]);
#endif
    rf[jj] = sqrt(pow(szPos[jj][0],2)+pow(szPos[jj][1],2));
    pf[jj] = atan(szPos[jj][1]/szPos[jj][0]); //phi angle
    if(szPos[jj][1]>0 && szPos[jj][0]<0) pf[jj] +=PI;
    if(szPos[jj][1]<0 && szPos[jj][0]<0) pf[jj] +=PI;
    if(szPos[jj][1]<0 && szPos[jj][0]>0) pf[jj] +=2*PI;
    if(pf[jj]>2*PI)    pf[jj] -=2*PI;
    zf[jj] = szPos[jj][2];
    wfi[jj]= 1.0;
    wzf[jj]= 1.0;
  }
  
  if(fit_track_to_beamline)
  {
    rf[npt]= 0.0;
    pf[npt]= 0.0;
    zf[npt]= 0.0; 
    wfi[npt]= 1.0;
    wzf[npt]= 0.0; /* zero weight for Z on the beamline point*/
    //This means that don't calculate the chi square for Z on the beamline point
    npt++;
  }
  iopt=1; /* tells rwfthl what kind of fit to do */

  rwfthl_c(npt,rf,pf,wfi,zf,wzf,iopt,vv0,ee0,&ch2ph,&ch2z,del,delz);
  /*
  OUTPUT:      VV0 = 1/R*CHARGE   POS. IF CLOCKWISE
  C                  TAN(LAMBDA)  {=DZ/DS}TAN(ANGLE TO X,Y PLANE)
  C                  PHI0         {0,2PI} ANGLE TO X-AXIS at R=D0
  C                  D0*SIGN      [CM]    MINIMAL DIST. TO Z-AXIS,
  C                                       +VE IF AXIS ENCIRCLED
  C                  Z0           [CM]    Z POS AT R=D0
  */
  //reconstruct the output
  Rho  = (double)(1.0/vv0[0]); /* minimum distance to z=0 */
  phi0 = vv0[2]; /* in xy plane, direction of track relative to x axis */

  //This is from Nate, it is based on the definition of VV0
  my_phi = phi0+PI;
  if (vv0[0]<0.0) my_phi += PI;
  //vv0[0] negtive means curve to anti-CLOCKWISE direction, This is true in BONUS
  if(my_phi >  PI) my_phi -= 2.0*PI;
  if(my_phi < -PI) my_phi += 2.0*PI;


  //center of the circle
  A = (double)(-sin(my_phi)*((-vv0[3])+fabs(1.0/vv0[0])));
  B = (double)(cos(my_phi)*((-vv0[3])+fabs(1.0/vv0[0])));
  //position of the initial step
  Phi = (double)vv0[2];
  if(Phi> PI) Phi-=2*PI;
  if(Phi<-PI) Phi+=2*PI;
  Theta = PI/2.-atan(vv0[1]);
  X0 = -sin(my_phi)*(-vv0[3]);
  Y0 =  cos(my_phi)*(-vv0[3]);
  Z0 = (double)vv0[4];
  DCA = fabs(vv0[3]); /* dca = distance of closest approach to beamline */
  Chi2 = (npt>5)? (double)(ch2ph+ch2z/(npt-5)) : 9999.9;

  //By Jixie: aaply correction, only useful for RTPC12 
  CorrHelixRPhi(Rho,Phi);

#ifdef HELIXFIT_DEBUG
  printf("\nhelix_fit: fitting %d hits then return (a,b,r)= (%6.4f %6.4f %6.4f)\n",
    PointNum,X,Y,R);
#endif
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void HelixFit(int PointNum,double szPos[][3], double& R, double& A, double& B,
  double& Phi_deg, double& Theta_deg, double& Z0, int fit_track_to_beamline )
{
  const double PI=acos(0.0)*2;
  double Rho,Phi,Theta,X0,Y0,DCA,Chi2;
  helix_fit(PointNum, szPos, Rho, A, B, Phi, Theta, 
    X0, Y0, Z0, DCA, Chi2, fit_track_to_beamline);
  R=fabs(Rho);
  Phi_deg=Phi*180./PI;
  Theta_deg=Theta*180./PI;
}

