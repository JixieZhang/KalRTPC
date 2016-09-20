//This is a circle fitter using Levenberg-Marquardt method
#ifndef CIRCLE_FITTER_CLASS_LM_H
#define CIRCLE_FITTER_CLASS_LM_H


#define MaxNumberOfHit 200

class CircleFitter_LM
{
public:
  CircleFitter_LM(); //Constructor
  //int DoFit(int, int, double*, double*, double&, double&, double&);
  int DoFit(int nHits, double hitX[], double hitY[], 
    double & xOfCenter, double & yOfCenter, double & rFinal);
  int DoFit(int nHits, double hitXYZ[][3], 
  double & xOfCenter, double & yOfCenter, double & rFinal);

private:

  void Initialize(int arrSize, double hitX[], double hitY[]);

  void UpdateRadius(int arrSize, double hitX[], double hitY[]);
  void GetCircumCenterOfThisTriplet(double pIx, double pIy, double pJx, double pJy, double pKx, double pKy,
    bool &Collinearity, double & circumCenterX, double & circumCenterY);
  void GetCircumCenterSum(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y,  
    bool Collinearity, double circumCenter2x, double circumCenter2y);
  int Minimize(int iterMax, double innerThreshold, double outerThreshold, int arrSize, 
    double hitX[], double hitY[]);
  void ComputeCost(int arrSize, double hitX[], double hitY[]);
  double NewtonStep(double u, double v, int arrSize, double hitX[], double hitY[]);
  double GetRadius() { return rHat; }



private:
  double J, dJdx, dJdy;
  double xList[MaxNumberOfHit],yList[MaxNumberOfHit];

public:
  double xCenter, yCenter, rHat, xIniCenter, yIniCenter, rIniHat;

};

#endif
