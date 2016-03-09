// ********************************************************************
// $Id: BField_Helm.cc,v 3.0, 2011/1/19  G2P Exp $
// Implementation of the BField_Helm class.
//
//////////////////////////////////////////////////////////////////////

#include "BField_Helm.hh"
#include "ConfigManager.hh"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include  <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include <CLHEP/Vector/ThreeVector.h>

//#define CREATE_MAP_NTUPLE 1
//#define BFIELD_HELM_DEBUG 4
//#define CREATE_NTUPLE_BY_FORTRAN 1

#ifdef BFIELD_HELM_DEBUG
#include "GlobalDebuger.hh"
#endif

#ifdef CREATE_NTUPLE_BY_FORTRAN
//This is the fortran routine from John, it can calculate the field for a hemlhotz coil
extern "C" void hlmhltz_(float* x,float* y,float* z,float* r,float* amp,
						 float* bx,float* by,float* bz);

//a c++ wrraper to call the fortran routine
void GetHelmField_F(double pos[3],double b[3],double helmr=0.20, double amp=116640.2)
{
	float x=(float)pos[0], y=(float)pos[1], z=(float)pos[2];
	float bx,by,bz;
	float helmr_=(float)helmr;  //m
	float amp_=(float)amp;	 //current ratio
	hlmhltz_(&x,&y,&z,&helmr_,&amp_,&bx,&by,&bz);
	b[0]=(double)bx;
	b[1]=(double)by;
	b[2]=(double)bz;
}
#endif


//When UseHallBRoutine is defined, the field map will be overwritten by GetHallBField_F
//switch on the hallb coil routine, if UseHallBRoutine<2 will recalculate the buffer elements
//and will do interpolate the buffer when GetField() is called, 
//Otherwise will do the calculation for the given position in each call
//By Jixie: only use it for hall b coil, do not define it if you are using other coils 
//#define UseHallBRoutine 1

#ifdef UseHallBRoutine
//Here is the fortran routine from HallB to calculate the field for HallB coil, which is 
//not sysmatric
extern "C" void sda_ptf_(float* x,float* b);

//a c++ wrraper to call the fortran routine
void GetHallBField_F(double pos[3],double field[3])
{
	float x[3],b[3];
	for(int i=0;i<3;i++)
	{
		x[i]=(float)pos[i];
	}

	//this fortran routine will return the field in kG, maximum is 5.0938709T
	//Here I match it to the TOSCA map maximum 4.97884764722
	sda_ptf_(x,b);
	for(int i=0;i<3;i++)
	{
		field[i]=b[i]/10.0 * 4.97884764722/5.0938709;
		//field[i]=b[i]/10.0 * 5.10/5.0938709;
	}
}
#endif

using namespace std;

BField_Helm* BField_Helm::fInstance=0;
BField_Helm* BField_Helm::GetInstance()
{ 
	if(!fInstance)  
	{
		//new BField_Helm();
		cout<<"BField_Helm is not initialized yet...exit...\n";
		exit(-99);
	}
	return fInstance; 
}


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BField_Helm::BField_Helm(const char *inifile,const char *mapfile)
{
#ifdef BFIELD_HELM_DEBUG
	if(BFIELD_HELM_DEBUG > Global_Debug_Level)
		SetGlobalDebugLevel("BField_Helm::BField_Helm()", (int)BFIELD_HELM_DEBUG);
#endif

	fInstance=this;
	this->ReadIni(inifile);
	

	mDoShift=(mOrigin[0]*mOrigin[0]+mOrigin[1]*mOrigin[1]+mOrigin[2]*mOrigin[2]>0)?true:false;
	//mDoRotation=(fabs(mRotAngle[0])+fabs(mRotAngle[1])+fabs(mRotAngle[2])>0)?true:false;

	mDoRotation=false;
	// Default constructor. Gives a unit matrix
	CLHEP::HepRotation pRot[3];

	for(int i=0;i<3;i++)
	{
		if(mRotAxis[i]==0 || fabs(mRotAngle[i])<1.0E-10) continue;
		mDoRotation=true;
		if(mRotAxis[i]==1) pRot[i].rotateX(mRotAngle[i]);
		else if(mRotAxis[i]==2) pRot[i].rotateY(mRotAngle[i]);
		else if(mRotAxis[i]==3) pRot[i].rotateZ(mRotAngle[i]);
	}

	mRotL2F=new CLHEP::HepRotation();
	mRotF2L=new CLHEP::HepRotation();
	*mRotL2F=pRot[2]*pRot[1]*pRot[0];
	*mRotF2L=mRotL2F->inverse(); 
#ifdef BFIELD_HELM_DEBUG
	double rad2deg=180./acos(-1.0);
	if(Global_Debug_Level>=1)
	{
		cout<<"\nCLHEP Lab2Field EulerAngles: phi="<<mRotL2F->getPhi()*rad2deg
			<<"  theta="<<mRotL2F->getTheta()*rad2deg
			<<"  psi="<<mRotL2F->getPsi()*rad2deg<<endl;
	}
#endif

	mRNum=int((mRmax-mRmin)/mStepR)+1;
	if (mRNum<2)
	{
		mRNum=3; //set the minimum to 3
		mRmax=mRmin+3*mStepZ;
	}
	mZNum=int((mZmax-mZmin)/mStepZ)+1;
	if (mZNum<2)
	{
		mZNum=3; //set the minimum to 3
		mZmax=mZmin+3*mStepZ;
	}

	//dynamicly allocated the 3 dimensional array: double mBField[R][Z][Variables];  
	//stored in the order of z r Bz Br
	///////////allocate start//////////
	int i,j,k;
	mBField=new double **[mRNum];
	for(i=0;i<mRNum;i++)
	{
		mBField[i]=new double *[mZNum];
		for (j=0;j<mZNum;j++)
		{
			mBField[i][j]=new double [mNPara];
			//initial the array
			for (k=0;k<mNPara;k++)    mBField[i][j][k]=0.0;
		}
	}
	/////////////allocate end////////////

	if(mUseUniformB!=1) ReadMap(mapfile);

}

BField_Helm::~BField_Helm()
{
	int i,j;
	for(i=0;i<mRNum;i++)
	{
		for (j=0;j<mZNum;j++)
		{
			delete [] mBField[i][j];
		}
		delete [] mBField[i];
	}
	delete mRotL2F;
	delete mRotF2L;
}

/////////////////////////////////////////////////////////////////////
bool BField_Helm::ReadIni(const char *filename)
{
	double deg2rad=acos(-1.0)/180.; 
#ifdef _CONFIGMANAGER_H_
	//By Jixie: I am not use this routine to read ini file any longer
	//I prefer to use ConfigManager::ReadFile
	//
	ConfigManager *pConfig=new ConfigManager(filename);

	pConfig->GetParameter("Helm_UseUniformB",mUseUniformB);
	pConfig->GetParameter("Helm_UniformB_x",mUniformB[0]);
	pConfig->GetParameter("Helm_UniformB_y",mUniformB[1]);
	pConfig->GetParameter("Helm_UniformB_z",mUniformB[2]);
	pConfig->GetParameter("Helm_Rmin",		mRmin);
	pConfig->GetParameter("Helm_Rmax",		mRmax);
	pConfig->GetParameter("Helm_Zmin",		mZmin);
	pConfig->GetParameter("Helm_Zmax",		mZmax);
	pConfig->GetParameter("Helm_NPara",		mNPara);
	pConfig->GetParameter("Helm_StepR",		mStepR);
	pConfig->GetParameter("Helm_StepZ",		mStepZ);
	pConfig->GetParameter("Helm_OriginX",	mOrigin[0]);
	pConfig->GetParameter("Helm_OriginY",	mOrigin[1]);
	pConfig->GetParameter("Helm_OriginZ",	mOrigin[2]);
	pConfig->GetParameter("Helm_RotAxis1",	mRotAxis[0]);
	pConfig->GetParameter("Helm_RotAxis2",	mRotAxis[1]);
	pConfig->GetParameter("Helm_RotAxis3",	mRotAxis[2]);
	pConfig->GetParameter("Helm_RotAngle1",	mRotAngle[0]); mRotAngle[0]*=deg2rad;
	pConfig->GetParameter("Helm_RotAngle2",	mRotAngle[1]); mRotAngle[1]*=deg2rad;
	pConfig->GetParameter("Helm_RotAngle3",	mRotAngle[2]); mRotAngle[2]*=deg2rad;
	pConfig->GetParameter("Helm_CurrentRatio",mCurrentRatio);

#ifdef BFIELD_HELM_DEBUG
	pConfig->PrintParamMap(); 
#endif
	return true;
#else
	//The following is the original code
	//please keep it this way
	FILE *ini;
	char ch[]="=;\n";
	char *name,*value,line[100];
	char *pDest;
	int iPos=-1;

	if((ini=fopen(filename,"r"))==NULL)
	{
		printf("***Error! Can not open configure file \"%s\"!",filename);
		return false;
	}

	printf("\nThe magnetic field configuration is:\n");
	while(!feof(ini))
	{
		fgets(line,100,ini);
		/* Search forward. for the '#' to skip the comment*/
		pDest = strchr( line, '#' );
		iPos = pDest - line + 1;
		//in Linux, if not found '#', iPos==1073747457//
		if(iPos>0 && iPos<100) continue;
		name=strtok(line,ch);
		value=strtok(0,ch);
		//show the confif info

		if(name&&value) printf("%15s = %s\n",name,value);
		else printf("read %s error\n",filename);

		if (strcmp(name,"UseUniformB")==0)			mUseUniformB=atoi(value);
		else if (strcmp(name,"Helm_UniformB_x")==0)		mUniformB[0]=atof(value);
		else if (strcmp(name,"Helm_UniformB_y")==0)		mUniformB[1]=atof(value);
		else if (strcmp(name,"Helm_UniformB_z")==0)		mUniformB[2]=atof(value);
		else if (strcmp(name,"Helm_Rmin")==0)			mRmin=atof(value);
		else if (strcmp(name,"Helm_Rmax")==0)			mRmax=atof(value);
		else if (strcmp(name,"Helm_Zmin")==0)			mZmin=atof(value);
		else if (strcmp(name,"Helm_Zmax")==0)			mZmax=atof(value);
		else if (strcmp(name,"Helm_Para")==0)			mNPara=atoi(value);
		else if (strcmp(name,"Helm_StepR")==0)			mStepR=atof(value);
		else if (strcmp(name,"Helm_StepZ")==0)			mStepZ=atof(value);
		else if (strcmp(name,"Helm_OriginX")==0)		mOrigin[0]=atof(value);
		else if (strcmp(name,"Helm_OriginY")==0)		mOrigin[1]=atof(value);
		else if (strcmp(name,"Helm_OriginZ")==0)		mOrigin[2]=atof(value);
		else if (strcmp(name,"Helm_RotAxis1")==0)		mRotAxis[0]=atoi(value);
		else if (strcmp(name,"Helm_RotAxis2")==0)		mRotAxis[1]=atoi(value);
		else if (strcmp(name,"Helm_RotAxis3")==0)		mRotAxis[2]=atoi(value);
		else if (strcmp(name,"Helm_RotAngle1")==0)		mRotAngle[0]=atof(value)*deg2rad;
		else if (strcmp(name,"Helm_RotAngle2")==0)		mRotAngle[1]=atof(value)*deg2rad;
		else if (strcmp(name,"Helm_RotAngle3")==0)		mRotAngle[2]=atof(value)*deg2rad;
		else if (strcmp(name,"Helm_CurrentRatio")==0)	mCurrentRatio=atof(value);
		else continue;
	}
	fclose(ini);

	if(mNPara<4) mNPara=4; //mNPara should not less than 4 columns
	return true;
#endif
}

/////////////////////////////////////////////////////////////////////
bool BField_Helm::ReadMap(const char *filename)
{
	if(mUseUniformB==1) return false;

	char strLog[1024];
	sprintf(strLog,"BField_Helm::ReadMap() is loading field map %s......\n",filename);
	cerr<<strLog;

	ifstream ins;
	int indexR=0,indexZ=0,col=0;
	double tempLine[7];
	ins.open(filename);
	if (ins.fail())
	{
		sprintf(strLog,"***ERROR! Can not open field map %s...exit!***\n",filename);
		cerr<<strLog;
		exit(-1);
		return false;
	}

	//eat the first line
	char tempname[256];
	ins.getline (tempname,256);

	while (!ins.eof())
	{
		ins.getline(tempname,256);				
		//check if it is an empty line	
		if(strlen(tempname) < size_t(2*mNPara-1)) continue;

		istringstream s1(tempname);
		for(col=0;col<mNPara;col++)    s1>>tempLine[col];

		//check for R and Z
		if (tempLine[1]>=mRmin && tempLine[1]<=mRmax && tempLine[0]>=mZmin && tempLine[0]<=mZmax)
		{//store the value

			//in case there is an empty line, r=z=Btot=0.0
			if( fabs(tempLine[0]) < 1.0E-8 && fabs(tempLine[1]) < 1.0E-8 && 
				sqrt(tempLine[2]*tempLine[2]+tempLine[3]*tempLine[3]) < 1.0E-8 ) 
			{
				cout<<"***Warning: ZERO field at  r="<<tempLine[1]<<"  z="<<tempLine[0]<<endl;
				cout<<"There could be an empty line in the map "<<filename<<endl;
			}
			else 
			{
				indexR=int((tempLine[1]-mRmin)/mStepR);
				indexZ=int((tempLine[0]-mZmin)/mStepZ);				

				//the map is already in Tesla
				for(col=0;col<mNPara;col++) mBField[indexR][indexZ][col]=tempLine[col];
			}
		}
	}
	ins.close();

	
#if defined UseHallBRoutine && (UseHallBRoutine<2)
	//For hall b coil, we have some fortran routine to calculate the field, Here 
	//I want to recalculate the field and fill the buffer
			
#if defined BFIELD_HELM_DEBUG && (BFIELD_HELM_DEBUG>=3)			
		printf("   r(cm)     z(cm) Bz_tosca(T) Bz_cal(T) dBz(gauss)     dBz/Bz\n");
#endif

	double pX[3],pB[3];
	for (indexR=0;indexR<mRNum;indexR++)
	{
		for (indexZ=0;indexZ<mZNum;indexZ++)
		{
			pX[2]=mBField[indexR][indexZ][0];	//z
			pX[1]=0;
			pX[0]=mBField[indexR][indexZ][1];	//r
			
			GetHallBField_F(pX,pB);
			
#if defined BFIELD_HELM_DEBUG && (BFIELD_HELM_DEBUG>=3)			
		printf("%8.2f  %8.2f  %9.6f  %9.6f  %8.3f  %10.2E \n", pX[0], pX[2], mBField[indexR][indexZ][2],
			pB[2],(mBField[indexR][indexZ][2]-pB[2])*10000.0,(mBField[indexR][indexZ][2]-pB[2])/pB[2]);
#endif
			mBField[indexR][indexZ][2]=pB[2];  //Bz
			mBField[indexR][indexZ][3]=sqrt(pB[0]*pB[0]+pB[1]*pB[1]); //Br
		}
	}
#endif
	

#ifdef CREATE_MAP_NTUPLE
	char rootfile[255];
	double x=0.0,y=0.0,z=0.0,Bx=0.0,By=0.0,Bz=0.0,r=0.0,Br=0.0,Btot=0.0;
	sprintf(rootfile,"%s.root",filename);
	TFile *file=new TFile(rootfile,"RECREATE");
	TTree *field=new TTree("field","field map");
	field->Branch("x",&x,"x/D");
	field->Branch("y",&y,"y/D");
	field->Branch("r",&r,"r/D");
	field->Branch("z",&z,"z/D");
	field->Branch("Bx",&Bx,"Bx/D");
	field->Branch("By",&By,"By/D");
	field->Branch("Br",&Br,"Br/D");
	field->Branch("Bz",&Bz,"Bz/D");
	field->Branch("Btot",&Btot,"Btot/D");


#ifdef BFIELD_HELM_DEBUG
	if(Global_Debug_Level>=4)
	{
		printf("The Magnetic field Map is:\n");		
		printf("       x         y        z       Bx       By       Bz        r        Br        Btot\n");
	}
#endif
	for (indexR=0;indexR<mRNum;indexR++)
	{
		for (indexZ=0;indexZ<mZNum;indexZ++)
		{
			z=mBField[indexR][indexZ][0];
			r=mBField[indexR][indexZ][1];	
			Bz=mBField[indexR][indexZ][2];				
			Br=mBField[indexR][indexZ][3];	
			x=r;
			y=0.0;
			Bx=Br;
			By=0.0;
			Btot=sqrt(Bx*Bx+By*By+Bz*Bz);

			if(fabs(r)<1.0E-8 && fabs(z)<1.0E-08 && fabs(Btot)<1.0E-8) 
			{
				cout<<"***Warning: ZERO field at  x="<<x<<"  y="<<y<<"  z="<<z<<endl;
				cout<<"***Warning: There is an empty line in the map buffer..."<<endl;
			}
			else field->Fill();

#ifdef BFIELD_HELM_DEBUG
			if(Global_Debug_Level>=4) 
			{					
				printf("%8.3f %8.3f %8.3f %8.6f %8.6f %8.6f %8.3f %8.6f %8.6f\n",
					x,y,z,Bx,By,Bz,r,Br,Btot);

			}
#endif
			//reset 
			x=y=z=Bx=By=Bz=r=Br=Btot=0.0;

		}
	}
	file->Write();
	file->Close();
	cout << "Close root file " << rootfile <<endl;  
	file->Delete();
#endif

#ifdef BFIELD_HELM_DEBUG
	if(Global_Debug_Level>=4)
	{
		printf("The Magnetic field Map is:\n");
		printf("      z         r        Bz        Br     B_tot\n");
		for (indexR=0;indexR<mRNum;indexR++)
		{
			for (indexZ=0;indexZ<mZNum;indexZ++)
			{
				for(col=0;col<2;col++) printf(" %8.2f ",mBField[indexR][indexZ][col]);
				for(col=2;col<mNPara;col++) printf(" %8.6f ",mBField[indexR][indexZ][col]);
				if(mNPara<5) printf(" %8.6f ",sqrt(pow(mBField[indexR][indexZ][2],2.0)+pow(mBField[indexR][indexZ][3],2.0)));
				printf("\n");
			}
		}
	}
#endif

	return true;
}

/////////////////////////////////////////////////////////////////////
bool BField_Helm::Interpolation(double Pos[3],double B[3],int n)
{/*////////////////////////////////////////
 //function: calculate the nth order interpolation
 //1)found out B[R0-Rn][Z0-Zn], (n+1)x(n+1) matrix
 //2)using this matrix to interpolate by Z first, found B[R0-Rn][z], 2x(n+1) matrix
 //3)then interpolate by Radius to get Br and Bz
 //Input:
 // Pos[3]: the position in x,y,z coordinate in cm,origin variable
 // n : calculate the Nth order,n could be 1 2 3
 //Output:
 // B[3]: the magnetic field in B_x,B_y,B_z in the same unit as the map
 // if the position is beyond the field map, return 0 field 
 *//////////////////////////////////////////

#ifdef BFIELD_HELM_DEBUG
	//I do not want to check these in release version, just make sure you provide the right parameters
	if(n<1) n=1; 
	if(n>4) n=4; 

	if (n>=mZNum || n>=mRNum)
	{
		printf("\n**Error! Too few points to finish %dth order interpolation!**",n);
		printf("**Please change the config file to load more data points or use lower order Number!**\n");
		return false;
	}
#endif
	double r=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]);
	double z=fabs(Pos[2]);
	int StartIndexZ=0,StartIndexR=0;  //the first point to do the interpolation
	StartIndexZ=int((z-mBField[0][0][0])/mStepZ);
	StartIndexR=int((r-mBField[0][0][1])/mStepR);

#ifdef BFIELD_HELM_DEBUG 
	if(BFIELD_HELM_DEBUG>=4) 
	  {
	    printf("StartIndexZ=%d  StartIndexR=%d\n",StartIndexZ,StartIndexR);
	  }
#endif
	if ((StartIndexR<0 || StartIndexR>=mRNum) || 
		(StartIndexZ<0 || StartIndexZ>=mZNum) )
	{
		//if the position is beyond the field map, return 0 field 
		B[0]=B[1]=B[2]=0.0;return false;
	}
#ifdef BFIELD_HELM_DEBUG
	if (StartIndexZ<0 || StartIndexZ>=mZNum)
	{
		printf("\n**Warning!!! Z=%.4f is out of range [%.4f,%.4f) !!!**\n",
			z,mBField[0][0][0],mBField[0][mZNum-1][0]);
		//B[0]=B[1]=B[2]=0.0;return false;
	}
	if (StartIndexR<0 || StartIndexR>=mRNum)
	{
		printf("\n**Warning!!! R=%.4f is out of range [%.4f,%.4f) !!!**\n",
			r,mBField[0][0][1],mBField[mRNum-1][0][1]);
		//B[0]=B[1]=B[2]=0.0;return false;
	}
#endif

	//choose a best position to interpolate, considering the lower and upper limits
	if (StartIndexR>mRNum-1-n) StartIndexR=mRNum-1-n;	
	else if (StartIndexR<0) StartIndexR=0;
	if (StartIndexZ>mZNum-1-n) StartIndexZ=mZNum-1-n;
	else if (StartIndexZ<0) StartIndexZ=0;

	// interpolate by Z first, 
	double temp=1.0;
	//put the n+1=5 to avoid new and delete, this can save a lot of time
	double tempB[2][5];
	//double *tempB[2];
	//tempB[0]=new double [n+1];
	//tempB[1]=new double [n+1];
	//tempB[2][n+1];
	//tempB[0][Point 1,2,3] :  B_R[n+1],
	//tempB[1][Point 1,2,3] :  B_Z[n+1],

	int iR,iZ,ii,jj,kk;
	for (ii=0;ii<=n;ii++)
	{
		iR=StartIndexR+ii;
		//initial
		tempB[0][ii]=0.;
		tempB[1][ii]=0.;
		for (kk=0;kk<=n;kk++)
		{
			temp=1.;
			for (jj=0;jj<=n;jj++)
			{
				iZ=StartIndexZ+jj;
				if(jj!=kk)
				{
					temp*=(z-mBField[iR][iZ][0])/(mBField[iR][StartIndexZ+kk][0]-mBField[iR][iZ][0]);
				}
			}
			//find out B[R_ii][Z]
			tempB[1][ii]+=temp*mBField[iR][StartIndexZ+kk][2];	//Bz
			tempB[0][ii]+=temp*mBField[iR][StartIndexZ+kk][3];  //Br
		}
	}

	//interpolation by Radius to get B[r][z]
	double B_polar[]={0.,0.};
	iZ=StartIndexZ;
	for (kk=0;kk<=n;kk++)
	{
		temp=1.;
		for (jj=0;jj<=n;jj++)
		{
			iR=StartIndexR+jj;
			if(jj!=kk)
				temp*=(r-mBField[iR][iZ][1])/(mBField[StartIndexR+kk][iZ][1]-mBField[iR][iZ][1]);
		}
		//found out B[R_ii][Z]
		B_polar[0]+=temp*tempB[0][kk];	
		B_polar[1]+=temp*tempB[1][kk];
	}

	//convert from Polar to Cartesian coordinate
	//in case r==0, move it to 1.0E-08 cm to avoid divided by zero
	//Note that I only consider the 1st section in this interpolation, 
	//One needs to flip the sign for other sections 
	if(r<1.0E-08) r=1.0E-08;
	B[0]=B_polar[0]*fabs(Pos[0])/r;
	B[1]=B_polar[0]*fabs(Pos[1])/r;
	B[2]=B_polar[1];

	//this is no need if not using dynamic array
	//delete [] tempB[0];
	//delete [] tempB[1];
	return  true;
}

void BField_Helm::Rotate_Lab2Field(const double LabP[3],double FieldP[3])
{
	Hep3Vector pFieldP(LabP[0],LabP[1],LabP[2]);
	pFieldP.transform(*mRotL2F);
	FieldP[0]=pFieldP.x();
	FieldP[1]=pFieldP.y();
	FieldP[2]=pFieldP.z();
}

void BField_Helm::Transform_Lab2Field(const double LabP[3],double FieldP[3])
{
	for(int i=0;i<3;i++) FieldP[i]=LabP[i]-mOrigin[i];
	if(mDoRotation) Rotate_Lab2Field(FieldP,FieldP);
}


void BField_Helm::Rotate_Field2Lab(const double FieldP[3],double LabP[3])
{
	//the following 3 lines do the same thing, but the 3rd line will also change pFieldP
	//pLabP=(*mRotF2L)(pFieldP);
	//pLabP=mRotF2L->operator ()(pFieldP);
	//pLabP=pFieldP.transform(*mRotF2L);
	Hep3Vector pLabP(FieldP[0],FieldP[1],FieldP[2]);
	pLabP.transform(*mRotF2L);
	LabP[0]=pLabP.x();
	LabP[1]=pLabP.y();
	LabP[2]=pLabP.z();
}

void BField_Helm::Transform_Field2Lab(const double FieldP[3],double LabP[3])
{
	Rotate_Field2Lab(FieldP,LabP);
	LabP[0]+=mOrigin[0];
	LabP[1]+=mOrigin[1];
	LabP[2]+=mOrigin[2];
}


/////////////////////////////////////////////////////////////////////
bool BField_Helm::GetBField(const double Pos[3],double B[3])
{//input x,y,z in centimeter, return B field in Tesla
	int i;
	if(mUseUniformB==1)
	{
		for (i=0;i<3;i++) B[i]=mUniformB[i];
		return true;
	}
	if(fabs(mCurrentRatio)<1.0E-10)
	{
		for (i=0;i<3;i++) B[i]=0.0;
		return true;
	}

	double pPos[3],pB[3]={0,0,0},flag[3]={1.0,1.0,1.0};
	//shift and rotate the origin to the field coordinate
	if(mDoShift || mDoRotation) Transform_Lab2Field(Pos,pPos);
	else
	{
		for (i=0;i<3;i++) pPos[i]=Pos[i];
	}
	//the map only provide fields at the 1st section (X>0,Y>0,and Z>0) 
	//this field goes along z, B_z always point to the same direction
	//need to deal with the other 7 sections
	//in X-Z plane, flip B_x in 2nd and fouth secions
	//in Y-Z plane, flip B_y in 2nd and fouth secions
	//By Jixie and Chao @20110629: No need to flip the sign if Bx=Br*x/r and By=Br*y/r
	//I changed the Interpolation such that it only interpolate in the 1st section
	//I still need to flip the sign here
	if (pPos[0]<0) {flag[0]*=-1.0;}
	if (pPos[1]<0) {flag[1]*=-1.0;}
	if (pPos[2]<0) {flag[0]*=-1.0;flag[1]*=-1.0;}

	//do the interpolation
#if defined UseHallBRoutine && (UseHallBRoutine>=2)
	GetHallBField_F(pPos,pB);
#else
	if(!Interpolation(pPos,pB,2)) {B[0]=B[1]=B[2]=0.0; return false;}
#endif

#ifdef BFIELD_HELM_DEBUG
	if(Global_Debug_Level>=3)
	{
		printf("Input position(x,y,z) in cm=(%f, %f, %f):==>\n",Pos[0],Pos[1],Pos[2]);
		printf("Map position(x,y,z) in cm=(%f, %f, %f);\n",pPos[0],pPos[1],pPos[2]);
		printf("The raw magnetic field in Tesla without apply Helm_CurrentRatio:\n");
		printf("(B_x=%7.5f, B_y=%7.5f, B_z=%7.5f);\n",
			pB[0]*=flag[0],pB[1]*=flag[1],pB[2]*=flag[2]);
	}
#endif

	//apply real current ratio and the flag
	for (i=0;i<3;i++) B[i]=pB[i]*mCurrentRatio*flag[i];

#ifdef BFIELD_HELM_DEBUG
	if(Global_Debug_Level>=2)
	{
		printf("The magnetic field in Tesla after apply Helm_CurrentRatio: \n");
		printf("(B_x=%7.5f, B_y=%7.5f, B_z=%7.5f);\n",B[0],B[1],B[2]);
	}
#endif

	//rotate back to the Hall coordinate
	if(mDoRotation) 
	{
		Rotate_Field2Lab(B,B);

#ifdef BFIELD_HELM_DEBUG
		if(Global_Debug_Level>=2)
		{
			printf("The magnetic field in Tesla after apply Rotation: \n");
			printf("(B_x=%7.5f, B_y=%7.5f, B_z=%7.5f);\n",B[0],B[1],B[2]);
		}
#endif
	}

	return true;
}

/////////////////////////////////////////////////////////////////////
bool BField_Helm::GetBField(const float fPos[3],float fB[3])
{//input x,y,z in centimeter
	bool status=false;
	double dPos[3],dB[3];
	int i;
	for ( i=0;i<3;i++)
	{
		dPos[i]=(double) fPos[i];
	}
	status=GetBField(dPos,dB);

	for ( i=0;i<3;i++)
	{
		fB[i]=(float) dB[i];
	}
	return status;
}

/////////////////////////////////////////////////////////////////////
//Create root ntuple: if verbose != 0 will print steps out 
//void CreateNtuple(const char *rootfile="targetfield.root",int verbose=0,
//				  double xmin=-50.0,double xmax=50.0, double xstep=1.0,
//				  double ymin=-50.0,double ymax=50.0, double ystep=1.0,
//				  double zmin=-100.0,double zmax=100.0, double zstep=1.0);

void BField_Helm::CreateNtuple(const char *rootfile,int verbose,double xmin,double xmax,double xstep,
							   double ymin,double ymax,double ystep,double zmin,double zmax,double zstep)
{
	double x,y,r,z,Bx,By,Br,Bz,Btot;

	TFile *file=new TFile(rootfile,"RECREATE");
	TTree *field=new TTree("field","field map");
	field->Branch("x",&x,"x/D");
	field->Branch("y",&y,"y/D");
	field->Branch("r",&r,"r/D");
	field->Branch("z",&z,"z/D");
	field->Branch("Bx",&Bx,"Bx/D");
	field->Branch("By",&By,"By/D");
	field->Branch("Br",&Br,"Br/D");
	field->Branch("Bz",&Bz,"Bz/D");
	field->Branch("Btot",&Btot,"Btot/D");

	if(verbose)
	{
		printf("The Magnetic field Map is:\n");
		printf("       x         y        z       Bx       By       Bz        r        Br        Btot\n");
	}

	int pXNum=(int)ceil((xmax-xmin)/xstep)+1;
	int pYNum=(int)ceil((ymax-ymin)/ystep)+1;
	int pZNum=(int)ceil((zmax-zmin)/zstep)+1;

	double pPos[3],pB[3];

#ifdef CREATE_NTUPLE_BY_FORTRAN
	//double pHelmR=0.234, pAmp=511580.0*2.59;
	double pHelmR=0.200, pAmp=511580.0*2.28;
#endif

	for (int ix=0;ix<pXNum;ix++)
	{
		pPos[0]=x=xmin+xstep*ix;
		for (int iy=0;iy<pYNum;iy++)
		{
			pPos[1]=y=ymin+ystep*iy;
			r=sqrt(x*x+y*y);
			for (int iz=0;iz<pZNum;iz++)
			{
				pPos[2]=z=zmin+zstep*iz;

				//Get The Field value
#ifndef CREATE_NTUPLE_BY_FORTRAN
				GetBField(pPos,pB);
#else	
				GetHelmField_F(pPos,pB,pHelmR,pAmp);
#endif			
				Bx=pB[0]; By=pB[1]; Bz=pB[2];  //in Tesla
				Br=sqrt(Bx*Bx+By*By);
				Btot=sqrt(Bx*Bx+By*By+Bz*Bz);

				if(fabs(Btot)<1.0E-13) 
				{
					cout<<"***Warning: Zero field at  x="<<x<<"  y="<<y<<"  z="<<z<<endl;
				}

				field->Fill();

				if(verbose)
				{
					printf("%8.3f %8.3f %8.3f %8.6f %8.6f %8.6f %8.3f %8.6f %8.6f\n",
						x,y,z,Bx,By,Bz,r,Br,Btot);
				}
			}
		}
	}

	file->Write();
	file->Close();
	file->Delete();

}
/////////////////////////////////////////////////////////////////////
