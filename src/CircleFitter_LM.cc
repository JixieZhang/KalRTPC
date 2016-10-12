//This is a circle fitter using Levenberg-Marquardt method

#include "CircleFitter_LM.h"
#include <iostream>
#include <cmath>

using namespace std;

//#define CircleFitter_LM_Debug 6

CircleFitter_LM::CircleFitter_LM() //Constructor (initializes the private data members)
{
  xCenter=0.0; yCenter=0.0; rHat=0.0; xIniCenter=0.0; yIniCenter=0.0; rIniHat=0.0;
  J = 0; dJdx = 0; dJdy = 0;
}

//We can remove iEvent and verbal from all the routines here. They are all for testing/debugging purposes
//There are some cases where it returns NAN results for the final fit, but the intial estimate is 
//  still not-NAN, therefore, we may use the initial estimate as input to the Kalman Filter.
//   Krishna Adhikari: 9/9/16
//
//   verbal               - a verbal flag (-1 gives minimal printout, can change to no-printout later)
//   nHits - # of hits to be used in the fit (15 was found to give the best result)
//                         but, if the # of good hits 'nGoodHits' is less than 15, then this # is
//                         set equal to 'nGoodHits'. 
//                         The # of good hits is counted as follows:  
//                             for(int i=0; i<allHits; i++) {if(step_status[i]==1) goodHits++;} 
//   hitX, hitY          - arrays of type double and of size nFitPoints=15 carrying (X,Y) values 
//                         of nFitPoints. If goodHits (# of good hits) is less than nFitPoints, 
//                         then some of the latter cells which wont have hit-info will be ignored
//                         based on the value of 'goodHits'
//  
//  Returns              - xOfCenter, yOfCenter, rFinal, fitStatus
//                         'fitStatus'  -  an integer: -1, 0 , 1 is for bad, okay & good fits
//
int CircleFitter_LM::DoFit(int nHits, double hitX[], double hitY[], 
  double & xOfCenter, double & yOfCenter, double & rFinal)
{
  //int nHits = nFitPoints; if(goodHits<nFitPoints-1) nHits = goodHits+1;
  xList[0]=0.0; yList[0]=0.0;
  int npt=1;
  for(int i=0;i<nHits;i++) {
    if(npt>=MaxNumberOfHit) break;
    xList[npt]=hitX[i]; 
    yList[npt]=hitY[i];
    npt++;
  }

  Initialize(npt, xList, yList);
  //if(goodHits<nFitPoints-1) Initialize(goodHits+1, hitX, hitY);
  //else                      Initialize(nFitPoints, hitX, hitY);
  
  int iter = 0, fitStatus=0;
  iter = Minimize(100, 0.1, 1.0e-12, npt, xList, yList);
  //if(goodHits<nFitPoints-1) iter = Minimize(100, 0.1, 1.0e-12, goodHits+1, hitX, hitY);
  //else                      iter = Minimize(100, 0.1, 1.0e-12, nFitPoints, hitX, hitY);
  
#ifdef  CircleFitter_LM_Debug
  if(CircleFitter_LM_Debug==1 || CircleFitter_LM_Debug == 6) {
    cout<<"\nAfter Initialize(...) call the initial fit values of x,y, and r are: \n"
    <<xIniCenter<<"  "<<yIniCenter<<" "<<rIniHat<<endl;
  }

  if(CircleFitter_LM_Debug==2 || CircleFitter_LM_Debug==6){
    cout<<"Converged after "<<iter<<" iterations.\n";
  }
#endif

  double rFin = GetRadius();
  xOfCenter = xCenter;   yOfCenter = yCenter;  rFinal = rFin;
  if( !(xCenter==xCenter) || !(yCenter==yCenter)  || !(rFin==rFin) )
  {
    xOfCenter = xIniCenter;   yOfCenter = yIniCenter;  rFinal = rIniHat; fitStatus = 0; 
    if( !(xIniCenter==xIniCenter) || !(yIniCenter==yIniCenter) || !(rIniHat==rIniHat) ) 
      fitStatus = -1;
  }
  else {xOfCenter = xCenter;   yOfCenter = yCenter;  rFinal = rFin; fitStatus = 1;}
  
#ifdef  CircleFitter_LM_Debug
  if(CircleFitter_LM_Debug==1 || CircleFitter_LM_Debug == 6) {
    cout<<"After Minimize(...) x,y, and r are: \n"
    <<xOfCenter<<"  "<<yOfCenter<<" "<<rFinal<<endl;
  }
#endif

  return fitStatus; 
}

int CircleFitter_LM::DoFit(int nHits, double hitXYZ[][3], 
  double & xOfCenter, double & yOfCenter, double & rFinal)
{
  double XX[200],YY[200];
  for(int i=0;i<nHits;i++) {
  XX[i]=hitXYZ[i][0];
  YY[i]=hitXYZ[i][1];
  }
  return DoFit(nHits,XX,YY,xOfCenter,yOfCenter,rFinal);
}

//XYPoint p1, p2, p3, circumCenter2, origin; 
//double p1x, p1y, p2x, p2y, p3x,p3y, circumCenter2x, circumCenter2y, originX, originY; 
void CircleFitter_LM::GetCircumCenterSum(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y,  
  bool Collinearity, double circumCenter2x, double circumCenter2y)
{
  double kpRadius=0.0, kpX=0.0, kpY=0.0, kpP1X=0.0, kpP1Y=0.0;
  GetCircumCenterOfThisTriplet(p1x,p1y,p2x,p2y,p3x,p3y,Collinearity,circumCenter2x,circumCenter2y); 
  if(Collinearity==true) return;//continue; //artificially given large value for GetCircumCenterOfThisTriplet of aligned triplets
  kpX = circumCenter2x; kpY = circumCenter2y; kpP1X = p1x; kpP1Y = p1y;
  kpRadius = sqrt( pow(kpX-kpP1X,2) + pow(kpY-kpP1Y,2));
  /*
  kpX = circumCenter2x; kpY = circumCenter2y; kpP1X = p1x; kpP1Y = p1y;
  if(!(kpX==kpX) || !(kpY==kpY) || !(kpP1X==kpP1X) || !(kpP1Y==kpP1Y) ) continue; //Avoid cases with 'nan' values

  kpRadius = sqrt( pow(kpX-kpP1X,2) + pow(kpY-kpP1Y,2)); 
  if(!(kpRadius==kpRadius)) continue;	
  */
  if(fabs(kpX) > 1000.0 || fabs(kpY) > 1000.0 || fabs(kpRadius) >1000) return;//continue;

  xCenter += circumCenter2x;
  yCenter += circumCenter2y;
}

void CircleFitter_LM::Initialize(int arrSize, double hitX[], double hitY[])//(XYPoint point [])
{
  bool Collinearity = false;
  xCenter=yCenter=0.0; 
  xIniCenter=yIniCenter=rIniHat=0.0;

  double p1x, p1y, p2x, p2y, p3x, p3y, circumCenter2x, circumCenter2y, originX, originY; 
  originX = 0.0, originY = 0.0; //origin.SetX(0.0); origin.SetY(0.0); 
  int nTriplets = 0;            //Number of non-aligned triplets
  circumCenter2x=circumCenter2y=0.0;

  //I want to try only those triplets which always uses the origin as one of the points, instead
  //  of using all possible combinations to reduce the computation time.

  p1x = originX; p1y = originY; //p1 = origin;   
  if(arrSize<9) 
  {
    //Choosing only some specific points (skipping some) to get the initial estimate  
    p3x = hitX[1]; p3y = hitY[1]; ////First hit
    for (int i = 2; i < arrSize; i++)
    {
      p2x = hitX[i]; p2y = hitY[i]; //p2 = points[i]; //any hit except 0th & 1st
      GetCircumCenterSum(p1x,p1y,p2x,p2y,p3x,p3y,Collinearity,circumCenter2x,circumCenter2y); 
      nTriplets++;	  
    }
    p3x = hitX[arrSize-1]; p3y = hitY[arrSize-1]; //Last of the selected good hits for 3rd point, 1st is (0,0) as before
    for (int i = 1; i < arrSize-1; i++)
    {
      p2x = hitX[i]; p2y = hitY[i]; //Any hit between hit1 and 3rd from last selected hit.
      GetCircumCenterSum(p1x,p1y,p2x,p2y,p3x,p3y,Collinearity,circumCenter2x,circumCenter2y); 
      nTriplets++;	  
    }	
  }
  else
  {
    //kp: In this block, I want to pair hits above half-way point to the initial hit
    //    and the other-half to the last hit in the sample selected for the fit, to keep distances large
    //    during this initial estimation.
    //
    for (int i = 1; i < arrSize-1; i++) {
      p2x = hitX[i]; p2y = hitY[i]; 
      if(i<arrSize/2) {p3x = hitX[arrSize-1]; p3y = hitY[arrSize-1]; }  
      else { p3x = hitX[1]; p3y = hitY[1];  }        
      GetCircumCenterSum(p1x,p1y,p2x,p2y,p3x,p3y,Collinearity,circumCenter2x,circumCenter2y); 
      nTriplets++;	
    }
  }

  if(nTriplets>0)
  {
    xCenter /= nTriplets;
    yCenter /= nTriplets;
  }

  UpdateRadius(arrSize, hitX, hitY);
  xIniCenter = xCenter; yIniCenter = yCenter; rIniHat = rHat;
}



void CircleFitter_LM::GetCircumCenterOfThisTriplet(double pIx, double pIy, double pJx, double pJy, double pKx, double pKy,
  bool &Collinearity, double & circumCenterX, double & circumCenterY) 
{
  Collinearity = false; 

  // some temporary variables
  double dIJx = pJx - pIx,      dIJy = pJy - pIy;
  double dJKx = pKx - pJx,      dJKy = pKy - pJy;
  double dKIx = pIx - pKx,      dKIy = pIy - pKy;
  double sqI  = pIx * pIx + pIy * pIy;
  double sqJ  = pJx * pJx + pJy * pJy;
  double sqK  = pKx * pKx + pKy * pKy;

  // determinant of the linear system: 0 for aligned points
  double det = dJKx * dIJy - dIJx * dJKy; 

  //if(det==0) { //For aligned points, give an unphysically large number for x and y
  if(fabs(det)<1.0e-50) { //Sometimes I was getting non-zero but very small value like det = 1.37943e-320
    Collinearity = true; 
    circumCenterX = 1.0e10; circumCenterY = 1.0e10;
    return; //return circumCenter;
  }
  else Collinearity = false; 

  double x_center = (sqI * dJKy + sqJ * dKIy + sqK * dIJy) / (2 * det);
  double y_center = -(sqI * dJKx + sqJ * dKIx + sqK * dIJx) / (2 * det);

  circumCenterX = x_center; circumCenterY = y_center; 
  return; //return circumCenter;
}


void CircleFitter_LM::UpdateRadius(int arrSize, double hitX[], double hitY[])
{
  rHat = 0;
  for (int i = 0; i < arrSize; ++i) {
    double dx = hitX[i] - xCenter;
    double dy = hitY[i] - yCenter;
    rHat += sqrt(dx * dx + dy * dy);
  }

  rHat /= arrSize; //kp: average of all likely radii (dist. betn.points and the initial center)
}



/** Minimize the distance residuals between the points and the circle.
* <p>We use a non-linear conjugate gradient method with the Polak and
* Ribiere coefficient for the computation of the search direction. The
* inner minimization along the search direction is performed using a
* few Newton steps. It is worthless to spend too much time on this inner
* minimization, so the convergence threshold can be rather large.</p>
* @param maxIter maximal iterations number on the inner loop (cumulated
* across outer loop iterations)
* @param innerThreshold inner loop threshold, as a relative difference on
* the cost function value between the two last iterations
* @param outerThreshold outer loop threshold, as a relative difference on
* the cost function value between the two last iterations
* @return number of inner loop iterations performed (cumulated
* across outer loop iterations)
* @exception LocalException if we come accross a singularity or if
* we exceed the maximal number of iterations
*/
//void ComputeCost(int arrSize, double hitX[], double hitY[]);
//double NewtonStep( double u, double v, int arrSize, double hitX[], double hitY[]);
int CircleFitter_LM::Minimize(int iterMax, double innerThreshold, double outerThreshold,
  int arrSize, double hitX[], double hitY[])
{
  ComputeCost(arrSize, hitX, hitY);

  if ((J < 1.0e-10) || (sqrt(dJdx * dJdx + dJdy * dJdy) < 1.0e-10)) 
  {
    // we consider we are already at a local minimum
    return 0;
  }
  double previousJ = J; 
  double previousU = 0.0, previousV = 0.0;
  double previousDJdx = 0.0, previousDJdy = 0.0;
  for (int iterations = 0; iterations < iterMax;) 
  {
    // search direction
    double u = -dJdx,       v = -dJdy;

    if (iterations > 0) 
    {
      // Polak-Ribiere coefficient
      double beta =
	(dJdx * (dJdx - previousDJdx) + dJdy * (dJdy - previousDJdy))
	/ (previousDJdx * previousDJdx + previousDJdy * previousDJdy);
      u += beta * previousU;
      v += beta * previousV;
    }
    previousDJdx = dJdx;
    previousDJdy = dJdy;
    previousU    = u;
    previousV    = v;

    // rough minimization along the search direction
    double innerJ;
    do 
    {
      innerJ = J;
      double lambda = NewtonStep(u, v, arrSize, hitX, hitY);
      xCenter += lambda * u;  yCenter += lambda * v;

      UpdateRadius(arrSize, hitX, hitY);
      ComputeCost(arrSize, hitX, hitY);
    } 
    while ((++iterations < iterMax)
      && ((fabs(J - innerJ) / J) > innerThreshold));

    // global convergence test
    if ((fabs(J - previousJ) / J) < outerThreshold) 
    {
      return iterations;
    }
    previousJ = J;
  }    
  //throw new LocalException("unable to converge after " + iterMax + " iterations");
  return 0;
}





/** Compute the cost function and its gradient.
* <p>The results are stored as instance attributes.</p>
*/
void CircleFitter_LM::ComputeCost(int arrSize, double hitX[], double hitY[]) {

  double dx=0.0, dy = 0.0, di=0.0, dr = 0.0, ratio=0.0, xx = 0.0, yy = 0.0; //8/28/16
  J    = 0;    dJdx = 0;    dJdy = 0; 
  for (int i = 0; i < arrSize; ++i) {
    //if(!(xCenter==xCenter) || !(yCenter==yCenter)) continue;
    //if(rHat>1000.0) continue;
    //xx = points[i].x();         yy = points[i].y();
    xx = hitX[i];         yy = hitY[i];
    dx = xx - xCenter;          dy = yy - yCenter;
    di = sqrt(dx * dx + dy * dy); //kp: 9/4/16: radius w.r.t the current center
    /*
    if (di < 1.0e-10) {
    throw new LocalException("cost singularity:" + " point at the circle center"); }
    */
    dr    = di - rHat;
    ratio = dr / di; //kp: 9/4/16: dr/r i.e. relative dr

    //cout<<"dx, dy, di: "<<dx<<"  "<<dy<<"  "<<di<<endl;

    //kp: 9/4/16: Now the weighted sums
    J    += dr * (di + rHat); //kp: 9/4/16: radial diff. 'dr' weighted by (di+rHat) which is roughly 2*rHat
    dJdx += dx * ratio; //kp: 9/4/16: weighted sum of dx (with weight = ratio)
    dJdy += dy * ratio; //kp: 9/4/16: weighted sum of dy (with weight = ratio)
  }
  dJdx *= 2.0;
  dJdy *= 2.0;
}


/** Compute the length of the Newton step in the search direction.
* @param u abscissa of the search direction
* @param v ordinate of the search direction
* @return value of the step along the search direction
*/

double CircleFitter_LM::NewtonStep(double u, double v, int arrSize, double hitX[], double hitY[]) {
  // compute the first and second derivatives of the cost
  // along the specified search direction
  double sum1 = 0, sum2 = 0, sumFac = 0, sumFac2R = 0;
  for (int i = 0; i < arrSize; ++i) {
    double dx     = xCenter - hitX[i]; 
    double dy     = yCenter - hitY[i]; 
    double di     = sqrt(dx * dx + dy * dy);
    double coeff1 = (dx * u + dy * v) /  di;
    double coeff2 = di - rHat;
    sum1         += coeff1 * coeff2;
    sum2         += coeff2 / di;
    sumFac       += coeff1;
    sumFac2R     += coeff1 * coeff1 / di;
  }

  // step length attempting to nullify the first derivative
  //double NewtonStep = -sum1 / ((u * u + v * v) * sum2  - sumFac * sumFac / arrSize  + rHat * sumFac2R);//8/29/16
  double denom1 = (u * u + v * v) * sum2,    denom2 = sumFac * sumFac / arrSize,     denom3 = rHat * sumFac2R;
  //double denom32 =  denom3 - denom2, denomTot = denom1 + denom32; 
  double denomTot = denom1 + denom3 - denom2 ;
  //Following line to avoid cases like the one when the values were -5.05178e-23, 3.34879e-11, 3.34879e-11, 0
  //   for denom1, denom2, denom3 and denomTot
  if(denomTot==0 && denom3==denom2) denomTot= denom1; if(denomTot==0 && denom1==denom2) denomTot= denom3;

  double NewtonStep = -sum1 / denomTot;//8/29/16

  return NewtonStep;
}
