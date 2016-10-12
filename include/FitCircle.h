//Header file for circle fit, c++ code!!!
#ifndef _FIT_CIRCLE_
#define _FIT_CIRCLE_
///////////////////////////////////////////////////////////////////////
namespace FitCircle {

  double MatrixDeterminant(double **a,int n);
  void   MatrixTranspose(double **a,int n);
  void   MatrixReverse(double **m, int n, double **r);
  void   MatrixTimesVector(double **m,int n,double *v, double *r);

  //use 3-D linear regression to determine R,A,B
  //Fit y vs x for a circle
  //http://www.had2know.com/academics/best-fit-circle-least-squares.html
  void   FitCircle_Algebra(int npt,double szPos[][3], double &R, double &A, double &B);
  
  //Journal of Electronic Imaging 12(1), 179 ¨C 193 (January 2003).
  //Classical geometrical approach to circle fitting¡ªreview and new developments
  //
  void   FitCircle_Kasa(int npt,double szPos[][3], double &R, double &A, double &B);
  
}
////////////////////////////////////////////////////////////////////////
#endif //#ifdef _FIT_CIRCLE_
