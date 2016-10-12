// ********************************************************************
//
// $Id: FitCircle.cc,v 1.1 2016/04/18 Bonus12 Exp $
// By Jixie Zhang, Modify for RTPC12
// ********************************************************************

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

namespace FitCircle {
  /*
  Recursive definition of determinate using expansion by minors.
  */
  double MatrixDeterminant(double **a,int n)
  {
    int i,j,j1,j2;
    double det = 0;
    double **m = NULL;

    if (n < 1) { /* Error */

    } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
    } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
	m = (double**)malloc((n-1)*sizeof(double *));
	for (i=0;i<n-1;i++)
	  m[i] = (double*)malloc((n-1)*sizeof(double));
	for (i=1;i<n;i++) {
	  j2 = 0;
	  for (j=0;j<n;j++) {
	    if (j == j1)
	      continue;
	    m[i-1][j2] = a[i][j];
	    j2++;
	  }
	}
	det += pow(-1.0,j1+2.0) * a[0][j1] * MatrixDeterminant(m,n-1);
	for (i=0;i<n-1;i++)
	  free(m[i]);
	free(m);
      }
    }
    return(det);
  }

  /*
  Find the cofactor matrix of a square matrix
  */
  void MatrixCoFactor(double **a,int n,double **b)
  {
    int i,j,ii,jj,i1,j1;
    double det;
    double **c;

    c = (double**)malloc((n-1)*sizeof(double *));
    for (i=0;i<n-1;i++)
      c[i] = (double*)malloc((n-1)*sizeof(double));

    for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

	/* Form the adjoint a_ij */
	i1 = 0;
	for (ii=0;ii<n;ii++) {
	  if (ii == i)
	    continue;
	  j1 = 0;
	  for (jj=0;jj<n;jj++) {
	    if (jj == j)
	      continue;
	    c[i1][j1] = a[ii][jj];
	    j1++;
	  }
	  i1++;
	}

	/* Calculate the determinate */
	det = MatrixDeterminant(c,n-1);

	/* Fill in the elements of the cofactor */
	b[i][j] = pow(-1.0,i+j+2.0) * det;
      }
    }
    for (i=0;i<n-1;i++)
      free(c[i]);
    free(c);
  }

  /*
  Transpose of a square matrix, do it in place
  */
  void MatrixTranspose(double **a,int n)
  {
    int i,j;
    double tmp;

    for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
	tmp = a[i][j];
	a[i][j] = a[j][i];
	a[j][i] = tmp;
      }
    }
  }

  //nxn matrix m times vector v, and store the reulst in vector r
  void MatrixTimesVector(double **m,int n,double *v, double *r)
  {
    int i,j;
    for (i=0;i<n;i++) {
      r[i] = 0.0;
      for (j=0;j<n;j++) {
	r[i] += m[i][j]*v[j];
      }
    }
  }

  //nxn matrix m times scaler s, and store the reulst in m
  void MatrixTimesScaler(double **m,int n,double s)
  {
    int i,j;
    for (i=0;i<n;i++) {
      for (j=0;j<n;j++) {
	m[i][j] *= s;
      }
    }
  }

  //calculate the reverse matrix of m and store it in r
  void MatrixReverse(double **m, int n, double **r)
  {
    double det = MatrixDeterminant(m,n);
    //store the cofactor matrix in r
    MatrixCoFactor(m,n,r);
    MatrixTranspose(r,n);
    MatrixTimesScaler(r,n,1.0/det);
  }


  //use 3-D linear regression to determine R,A,B
  //Fit y vs x for a circle
  //http://www.had2know.com/academics/best-fit-circle-least-squares.html
  void FitCircle_Algebra(int npt,double szPos[][3], double &R, double &A, double &B)
  {
    ///////////////////////////////////////////////////////////////////////
    //Fit y vs x for a circle
    double **MM;	  //the matrix
    double *VPara,*VN;	  //the vector
    MM = new double* [3];
    for(int i=0;i<3;i++) {
      MM[i] = new double [3];
      MM[i][0]=MM[i][1]=MM[i][2]=0;
    }
    double **MMReverse;
    MMReverse = new double* [3];
    for(int i=0;i<3;i++) MMReverse[i] = new double [3];

    VN = new double [3]; 
    VPara = new double [3];
    VN[0]=VN[1]=VN[2];
    VPara[0]=VPara[1]=VPara[2];

    for(int i=0;i<npt;i++) {
      MM[0][0] += szPos[i][0]*szPos[i][0];
      MM[0][1] += szPos[i][0]*szPos[i][1];
      MM[0][2] += szPos[i][0];
      MM[1][1] += szPos[i][0]*szPos[i][1];
      MM[1][2] += szPos[i][1];
      double x2addy2 = szPos[i][0]*szPos[i][0]+szPos[i][1]*szPos[i][1];
      VN[0] += szPos[i][0] * x2addy2;
      VN[1] += szPos[i][1] * x2addy2;
      VN[2] += x2addy2;
    }
    MM[1][0] = MM[0][1];
    MM[2][0] = MM[0][2];
    MM[2][1] = MM[1][2];    
    MM[2][2] = npt;

    //now get the reverse of MM, then VPara = MMReverse x VN
    MatrixReverse(MM,3,MMReverse);
    MatrixTimesVector(MMReverse,3,VN,VPara);
    //extract R, A, B
    A = -VPara[0]/2;
    B = -VPara[1]/2;    
    R = sqrt(4*VPara[2]+VPara[0]*VPara[0]+VPara[1]*VPara[1])/2;  

    //free the memory
    delete VPara;
    delete VN;  
    for(int i=0;i<3;i++) {
      delete MM[i];
      delete MMReverse[i];
    }
  }
  
    
  //Journal of Electronic Imaging 12(1), 179 ¨C 193 (January 2003).
  //Classical geometrical approach to circle fitting¡ªreview and new developments
  //
  void   FitCircle_Kasa(int npt,double szPos[][3], double &R, double &A, double &B)
  {
    ///////////////////////////////////////////////////////////////////////
    //Fit y vs x for a circle
    double alpha=0,beta=0,gamma=0,delta=0,epsilon=0;
    double sumX=0,sumX2=0,sumXY=0,sumY=0,sumY2=0,sumX3=0,sumXY2=0,sumX2Y=0,sumY3=0;
    for(int i=0;i<npt;i++) {
      sumX  += szPos[i][0];
      sumX2 += szPos[i][0]*szPos[i][0];
      sumX3 += szPos[i][0]*szPos[i][0]*szPos[i][0];
      sumXY += szPos[i][0]*szPos[i][1];
      sumX2Y+= szPos[i][0]*szPos[i][0]*szPos[i][1];
      sumXY2+= szPos[i][0]*szPos[i][1]*szPos[i][1];
      sumY  += szPos[i][1];
      sumY2 += szPos[i][1]*szPos[i][1];
      sumY3 += szPos[i][1]*szPos[i][1]*szPos[i][1];
    }

    alpha  = 2*(sumX*sumX - npt*sumX2);
    beta   = 2*(sumX*sumY - npt*sumXY);
    gamma  = 2*(sumY*sumY - npt*sumY2);
    delta  = sumX2*sumX - npt*sumX3 + sumX*sumY2 - npt*sumXY2;
    epsilon= sumX2*sumY - npt*sumY3 + sumY*sumY2 - npt*sumX2Y;
    
    //now calculate R,A,B
    double R2=0;
    A = (delta*gamma-epsilon*beta)/(alpha*gamma-beta*beta); 
    B = (alpha*epsilon-beta*delta)/(alpha*gamma-beta*beta); 
    R2 = A*A + B*B + (sumX2-2*A*sumX+sumY2-2*B*sumY)/npt;
    R = sqrt(R2);
  }

}  //end of namespace
