// BField_Helm.h: interface for the BField_Helm class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(BFIELD_HELM_H)
#define BFIELD_HELM_H

#include  <math.h>
#include <CLHEP/Vector/Rotation.h>
using namespace CLHEP;

class BField_Helm
{
public:
	//Static method which returns the singleton pointer of this class.
	static BField_Helm* GetInstance();
private:
	static BField_Helm* fInstance;

public:
	BField_Helm(const char *inifile="BField_CLAS12Sol.ini", const char *mapfile="CLAS12Sol_rz.dat");
	virtual ~BField_Helm();
	bool GetBField(const double Pos[3],double B[3]);
	bool GetBField(const float fPos[3],float fB[3]);

private:
	bool ReadIni(const char *filename);
	bool ReadMap(const char *filename);
	bool Interpolation(double Pos[3],double B[3],int n=2);

public:

	void   Rotate_Field2Lab(const double FieldP[3],double LabP[3]);
	void   Transform_Field2Lab(const double FieldP[3],double LabP[3]);
	void   Rotate_Lab2Field(const double LabP[3],double FieldP[3]);
	void   Transform_Lab2Field(const double LabP[3],double FieldP[3]);
	double GetCurrentRatio(){ return mCurrentRatio;}
	void   SetCurrentRatio(double val){ mCurrentRatio=val;}
	bool   IsUniformField(){ return (mUseUniformB==0)?false:true;}
	void   GetUniformField(double pB[]){pB[0]=mUniformB[0];pB[1]=mUniformB[1];pB[2]=mUniformB[2];}
	void   GetOrigin(double pX[]){pX[0]=mOrigin[0];pX[1]=mOrigin[1];pX[2]=mOrigin[2];}

	CLHEP::HepRotation *GetRotation_L2F() {return mRotL2F;}
	void GetEulerAngles_L2F(double &phi, double &theta, double &psi) 
	{		
		phi=mRotL2F->getPhi();theta=mRotL2F->getTheta();psi=mRotL2F->getPsi(); 
	}

	//Create root ntuple: if verbose != 0 will print steps out 
	void CreateNtuple(const char *rootfile="targetfield.root",int verbose=0,
		double xmin=-50.0,double xmax=50.0, double xstep=1.0,
		double ymin=-50.0,double ymax=50.0, double ystep=1.0,
		double zmin=-100.0,double zmax=100.0, double zstep=1.0);

private:

	//dynalic allocated the 3 dimensional array: double mBField[indexR][indexZ][Variables];
	//mBField[indexR][indexZ][0] z
	//mBField[indexR][indexZ][1] r
	//mBField[indexR][indexZ][2] B_z
	//mBField[indexR][indexZ][3] B_r
	double ***mBField;

	//parameters from ini file
	int    mUseUniformB;
	double mUniformB[3];
	double mCurrentRatio;
	double mStepR,mStepZ;
	double mRmin,mRmax;
	double mZmin,mZmax;
	double mOrigin[3];
	int    mNPara;
	int    mRotAxis[3];
	double mRotAngle[3];

	bool   mDoShift, mDoRotation;
	int    mZNum,mRNum;
	CLHEP::HepRotation *mRotL2F;
	CLHEP::HepRotation *mRotF2L;

};
typedef BField_Helm HRSTargetField;
typedef BField_Helm G2PTargetField;

#endif // !defined(BFIELD_HELM_H)
