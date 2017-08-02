//Originally from Howard and Carlos
//Jixie change this code into a class
//
#include <iostream>
#include <math.h>
#include <stdio.h>
#include "TVector3.h"

#include <algorithm>
#include "ChainFinder.hh"

#include "TCanvas.h"
#include "TPad.h"
#include  "TPolyMarker3D.h"

//kRTPC_R_GEM1 and kRTPC_R_Cathode have been defined in EXKalDetector.h
#include "EXKalDetector.h"
//use '#include "EXKalDetector.h"' if you do not want to have the following 3 lines
//static const double kRTPC_R_GEM1 = 7.0;
//static const double kRTPC_R_Cathode = 3.0;
//static const double kRTPC_Length = 40.0;

///////////////////////////////////////////////////////////////////////
//how to define _ChainFinderDebug_?
//>=1: print hit pool
//>=3: add printing stored chain information of 'hit-id' and 'ThrownTID'
//>=4: add printing "Chain #: hit 'hit-id' is added to 'chain-position'" during searching a tree
//     also print the chain before sorting, only 'hit-id' and 'ThrownTID'
//>=5: add printing 2 more information of stored chain: 's' and 'phi'
//     add printing "seed=0   seed_pre=0   hit=1  : separation=  0.31 cm, angle=    0.0 deg    Chain 0: hit 1 is added to 1"
//>=6: add printing 1 more information of stored chain: 'TDC'  
//==9: will call BenchmarkSort(), printing the sorting result of each sort algorithm
//>=10:add printing chain buffer before and after sorting
//     also add printing the following: 
//"nhits=98 idxSmax=48  idxSmin=58  idxPhimax=58  idxPhimin=95
//         Smax=4.87793  Smin=2.99306  S[idxPhimax]=2.99306  S[idxPhimin]=2.996
//         Phimax=-16.8163  Phimin=-107.002  Phi[idxSmax]=-62.0518  Phi[idxSmin]=-16.8163"
//>=11: add printing details of processing backward chains
//>=12  add priting more information about insertion_sort
//>=13: add printing most details of insertion_sort
#define _ChainFinderDebug_ 3

#ifdef _ChainFinderDebug_
#include "GlobalDebuger.hh"
#include "TBenchmark.h"
#endif

///////////////////////////////////////////////////////////////////////

/*
//this is an example how to use this class 
void Example()
{
  ChainFinder pTF;  
  pTF->SetParameters(v1,v2,v3,v4,v5);
  for(int ii;ii<nevents;ii++) {
    pTF.Reset();
    pTF.PrepareHitPool(vvv);
    pTF.RemoveBadHitsFromPool();
    pTF.SearchChains();
    //call KF
    for(int i=0;i<fChainNum_Stored;i++) DoKalmanFilter();  
  }
}
*/

using namespace std;

static const double kPi = atan(1.0)*4;
static const double rad2deg = 180./(4.*atan(1.));


ChainFinder::ChainFinder() 
{
#if defined _ChainFinderDebug_ && defined _GLOBALDEBUGER_H_
  Global_Debug_Level=_ChainFinderDebug_;
#endif

  fHitNum = 0;
  fChainNum = 0;
  fChainNum_Stored = 0;
  for(int i=0;i<MAX_CHAINS_PER_EVENT;i++) {
    fHitNumInAChain[i] = 0;  
    for(int j=0;j<MAX_HITS_PER_CHAIN;j++) fHitIDInAChain[i][j]=-1;
  }
  //default values for these parameters
  Ini_Sep     = 1.0;             //cm  
  Max_Sep     = 1.1;             //cm
  Max_Sep_Ang = 23.3/rad2deg;    //rad
  Min_Sep     = 0.4;             //cm
  Min_Sep_Ang = 30.0/rad2deg;    //rad
}

ChainFinder::~ChainFinder() 
{
    //
}

void ChainFinder::Reset() 
{
  fHitNum = 0;
  for(int i=0;i<fChainNum;i++) {
    for(int j=0;j<fHitNumInAChain[i];j++) fHitIDInAChain[i][j]=-1;
    fHitNumInAChain[i] = 0;  
  }
  fChainNum = 0;
  fChainNum_Stored = 0;
}

//provide x,x,z in mm
void ChainFinder::PrepareHitPool_mm(int *id, int *tdc, int *adc, double *x_mm, double *y_mm, 
                                    double *z_mm, int n, int *throwntid, int append)
{
  double *x=new double [n];
  double *y=new double [n];
  double *z=new double [n];
  
  for(int i=0;i<n;i++) {
    x[i]=x_mm[i]/10.;
    y[i]=y_mm[i]/10.;
    z[i]=z_mm[i]/10.;
  }
  PrepareHitPool(id, tdc, adc, x, y, z, n, throwntid, append);
  
  delete x;
  delete y;
  delete z;
}

//provide x,x,z in cm
void ChainFinder::PrepareHitPool(int *id, int *tdc, int *adc, double *x, double *y,  
                                 double *z, int n, int *throwntid, int append)
{
  TVector3 pV3;
  if(append==0) fHitNum=0; 
  for(int i=0;i<n;i++) {
    if(fHitNum < MAX_HITS_PER_EVENT) {
      pV3.SetXYZ(x[i],y[i],z[i]);
      if (pV3.Perp()>kRTPC_R_GEM1+1.0 || pV3.Perp()<kRTPC_R_Cathode-1.0) continue;  
      
      fHitPool[fHitNum].ID=id[i];
      fHitPool[fHitNum].TDC=tdc[i];
      fHitPool[fHitNum].ADC=adc[i];
      fHitPool[fHitNum].X=x[i];
      fHitPool[fHitNum].Y=y[i];
      fHitPool[fHitNum].Z=z[i];
      fHitPool[fHitNum].Status=HUNTCHD;
      //The following is added by Jixie for sorting
      fHitPool[fHitNum].S=pV3.Perp();
      fHitPool[fHitNum].Phi=pV3.Phi();
      if(throwntid) fHitPool[fHitNum].ThrownTID=throwntid[i];
      else fHitPool[fHitNum].ThrownTID=-1;
      fHitPool[fHitNum].ChainInfo=-1;
      fHitNum++;
    } else {
      cout<<"MAX_HITS_PER_EVENT("<<MAX_HITS_PER_EVENT<<") is too small, tracks are potentially lost!\n";
    }
  }
}

//By Carlos: VECTOR CONTAINER ROUTINE
//provide x,x,z in cm
void ChainFinder::PrepareHitPool(vector<int> *id, vector<int> *tdc, vector<int> *adc, 
				   vector<double> *x, vector<double> *y, vector<double> *z, 
				   int n, vector<int> *throwntid, int append)
{
  TVector3 pV3;
  if(append==0) fHitNum=0; 
  for(int i=0;i<n;i++) {
    if(fHitNum < MAX_HITS_PER_EVENT) {
      pV3.SetXYZ((*x)[i],(*y)[i],(*z)[i]);
      if (pV3.Perp()>kRTPC_R_GEM1+1.0 || pV3.Perp()<kRTPC_R_Cathode-1.0) continue;  
      
      fHitPool[fHitNum].ID     = (*id)[i];//THIS IS PAD NUMBER
      fHitPool[fHitNum].TDC    = (*tdc)[i];
      fHitPool[fHitNum].ADC    = (*adc)[i];
      fHitPool[fHitNum].X      = (*x)[i];
      fHitPool[fHitNum].Y      = (*y)[i];
      fHitPool[fHitNum].Z      = (*z)[i];
      fHitPool[fHitNum].Status = HUNTCHD;
      //The following is added by Jixie for sorting
      fHitPool[fHitNum].S=pV3.Perp();
      fHitPool[fHitNum].Phi=pV3.Phi();

      if(throwntid) fHitPool[fHitNum].ThrownTID=(*throwntid)[i];
      else fHitPool[fHitNum].ThrownTID=-1;

      fHitPool[fHitNum].ChainInfo=-1;
      fHitNum++;
    } else {
      cout<<"MAX_HITS_PER_EVENT("<<MAX_HITS_PER_EVENT<<") is too small, tracks are potentially lost!\n";
    }
  }
}

void ChainFinder::InsertAHitToPool(int hitid, int id, int tdc, int adc, double x, double y, double z,
  int ThrownTID, int ChainInfo)
{
  for(int i=fHitNum;i>hitid;i--) {
      fHitPool[i].ID=fHitPool[i-1].ID;
      fHitPool[i].TDC=fHitPool[i-1].TDC;
      fHitPool[i].ADC=fHitPool[i-1].ADC;
      fHitPool[i].X=fHitPool[i-1].X;
      fHitPool[i].Y=fHitPool[i-1].Y;
      fHitPool[i].Z=fHitPool[i-1].Z;
      fHitPool[i].Status=fHitPool[i-1].Status;
      fHitPool[i].S=fHitPool[i-1].S;
      fHitPool[i].Phi=fHitPool[i-1].Phi;
      fHitPool[i].ThrownTID=fHitPool[i-1].ThrownTID;
      fHitPool[i].ChainInfo=fHitPool[i-1].ChainInfo;
  }

  fHitPool[hitid].ID=id;
  fHitPool[hitid].TDC=tdc;
  fHitPool[hitid].ADC=adc;
  fHitPool[hitid].X=x;
  fHitPool[hitid].Y=y;
  fHitPool[hitid].Z=z;
  fHitPool[hitid].Status=HUNTCHD;
  fHitPool[hitid].S=sqrt(x*x+y*y);
  fHitPool[hitid].Phi=atan2(y,x);
  fHitPool[hitid].ThrownTID=ThrownTID;
  fHitPool[hitid].ChainInfo=ChainInfo;
  fHitNum += 1;
}


int  ChainFinder::RemoveAHitFromPool(int hitid)
{
  int found = (hitid>=fHitNum) ? 0 : 1;
  if(!found) return 0;

  for(int i=hitid;i<fHitNum;i++) {
      fHitPool[i].X=fHitPool[i+1].X;
      fHitPool[i].Y=fHitPool[i+1].Y;
      fHitPool[i].Z=fHitPool[i+1].Z;
      fHitPool[i].Status=fHitPool[i+1].Status;
      fHitPool[i].S=fHitPool[i+1].S;
      fHitPool[i].Phi=fHitPool[i+1].Phi;
      fHitPool[i].ThrownTID=fHitPool[i+1].ThrownTID;
  }
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=2) {
    printf("\t*****Hit %2d is removed from pool!*****\n",hitid);
  }
#endif
  fHitNum -= 1;
  return found;
}

//This is to change status for hits before searching
//For example, if we want to remove bad hits or redundate hits
//to gain fast speed, we should remove these hits
int  ChainFinder::RemoveBadHitsFromPool()
{
  TVector3 pV3(0,0,0);
  TVector3 pV3_pre(0,0,0);
    // removes the same hits (perhaps redundant)
  for(int i = 0; i<fHitNum; i++) {
    if(fHitPool[i].Status & HITUNAV) continue;
    pV3.SetXYZ(fHitPool[i].X, fHitPool[i].Y, fHitPool[i].Z);
    
    if (pV3.Perp()>kRTPC_R_GEM1+1.0 || pV3.Perp()<kRTPC_R_Cathode-1.0) {
      fHitPool[i].Status |= HISUSED;
      //fHitNum -= RemoveAHitFromPool(i);
      continue;  
    }
    
    // removes the same hits (perhaps redundant)
    if(i>0) pV3_pre.SetXYZ(fHitPool[i-1].X, fHitPool[i-1].Y, fHitPool[i-1].Z);
    if ( (pV3-pV3_pre).Mag() <= 1.0E-5 ) {
      fHitPool[i].Status |= HISUSED;
      //fHitNum -= RemoveAHitFromPool(i);
      continue;
    } 
  }
  return 0;
}

//By Jixie: ----do not call this routine----
//RTPC12 will not use the same channel map as BoNuS6, therefore the order
//of HitPool is not identical to BoNuS6
//In RTPC12, we can only sort by connecter row increasing order. For the same 
//connector row, chan id increases along z
//the real data in bonus6 it in channel ID increasing order, for the same ID,
//TDC in decreasing order. This routine is only for simulation with BoNuS6 data
//
//sort hitpool to match this order
void ChainFinder::SortHitPoolByIDTDC()
{
  //make a copy of the whole public hit poll
  HitStruct *pHitPool = new HitStruct [fHitNum];
  int *pOrder = new int [fHitNum]; 
  for(int t=0;t<fHitNum;t++) {
    pOrder[t] = t;
    pHitPool[t].ID = fHitPool[t].ID;
    pHitPool[t].TDC = fHitPool[t].TDC;
    pHitPool[t].ADC = fHitPool[t].ADC;
    pHitPool[t].X = fHitPool[t].X;
    pHitPool[t].Y = fHitPool[t].Y;
    pHitPool[t].Z = fHitPool[t].Z;
    pHitPool[t].Status = fHitPool[t].Status;
    pHitPool[t].S = fHitPool[t].S;
    pHitPool[t].Phi = fHitPool[t].Phi;
    pHitPool[t].ThrownTID = fHitPool[t].ThrownTID;
    pHitPool[t].ChainInfo = fHitPool[t].ChainInfo;
  }
  
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=9) {  
    PrintHitPool("before sorting:");
  }
#endif

  QuickSort_IDTDC(pOrder, 0, fHitNum-1);

  //copy the value back to the public hit pool
  for(int i=0;i<fHitNum;i++) {
    int t = pOrder[i];
    fHitPool[i].ID =  pHitPool[t].ID;
    fHitPool[i].TDC =  pHitPool[t].TDC;
    fHitPool[i].ADC =  pHitPool[t].ADC;
    fHitPool[i].X =  pHitPool[t].X;
    fHitPool[i].Y =  pHitPool[t].Y;
    fHitPool[i].Z =  pHitPool[t].Z;
    fHitPool[i].Status =  pHitPool[t].Status;
    fHitPool[i].S =  pHitPool[t].S;
    fHitPool[i].Phi =  pHitPool[t].Phi;
    fHitPool[i].ThrownTID =  pHitPool[t].ThrownTID;
    fHitPool[i].ChainInfo =  pHitPool[t].ChainInfo;
  }
  
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=9) {  
    PrintHitPool("After sorted by ID_increasing then TDC_decreasing");
  }
#endif

  delete pHitPool;
  delete pOrder;
}

  
//sort by ID increasing order, for the same ID, sort TDC in decreasing order. 
void ChainFinder::PrintHitPool(const char *keywords)
{
  //print the keywords then content
  printf("\nHitPool, %d hits: %s\n",fHitNum, keywords);

  printf("%9s: ","PAD-ID");
  for(int i=0;i<fHitNum;i++) { printf("%7d",fHitPool[i].ID);}
  printf("\n");
  printf("%9s: ","TDC");
  for(int i=0;i<fHitNum;i++) { printf("%7d",fHitPool[i].TDC);}
  printf("\n");
  printf("%9s: ","ThrownTID");
  for(int i=0;i<fHitNum;i++) { printf("%7d",fHitPool[i].ThrownTID);}
  printf("\n");

#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=3) {    
    printf("%9s: ","S");
    for(int i=0;i<fHitNum;i++) { printf("%7.2f",fHitPool[i].S);}
    printf("\n");    
    printf("%9s: ","Phi");
    for(int i=0;i<fHitNum;i++) { printf("%7.1f",fHitPool[i].Phi*rad2deg);}
    printf("\n");
  }
#endif

  printf("\n");

}
  
//sort by ID increasing order, for the same ID, sort TDC in decreasing order. 
void ChainFinder::QuickSort_IDTDC(int *arr, int left, int right)
{
  int i = left, j = right;
  int tmp;
  int pivot = arr[(left + right) / 2];
  int ID_pivot = fHitPool[pivot].ID;
  int TDC_pivot = fHitPool[pivot].TDC;

  // partition 
  while (i <= j) {
    while (fHitPool[arr[i]].ID < ID_pivot || (fHitPool[arr[i]].ID == ID_pivot && fHitPool[arr[i]].TDC > TDC_pivot) ) i++;
    while (fHitPool[arr[j]].ID > ID_pivot || (fHitPool[arr[j]].ID == ID_pivot && fHitPool[arr[j]].TDC < TDC_pivot) ) j--;
    //while (fHitPool[arr[i]].ID < ID_pivot ) i++;
    //while (fHitPool[arr[j]].ID > ID_pivot ) j--;

    //swap i and j
    if (i <= j) {
      //cout<<"  QuickSort_IDTDC()  swap("<<i<<", "<<j<<") "<<endl;
      tmp = arr[i];
      arr[i] = arr[j];
      arr[j] = tmp;
      i++;
      j--;
    }
  };

  // recursion 
  if (left < j) QuickSort_IDTDC(arr, left, j);
  if (i < right) QuickSort_IDTDC(arr, i, right);
}


void ChainFinder::SetParameters(double max_sep, double max_sep_ang, double min_sep, 
                                double min_sep_ang, double ini_sep)
{
  // USER SHOULD PROVIDE THESE VARIABLES
  Max_Sep     = max_sep;
  Max_Sep_Ang = max_sep_ang;
  Min_Sep     = min_sep;
  Min_Sep_Ang = min_sep_ang;
  Ini_Sep     = ini_sep;

#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=1) {
  cout<<"\nChainFinder Parameters:"
      <<" Max_Sep="<<Max_Sep<<" cm, "
      <<" Max_Sep_Ang="<<Max_Sep_Ang*rad2deg<<" deg,"
      <<" Min_Sep="<<Min_Sep<<" cm, "
      <<" Min_Sep_Ang="<<Min_Sep_Ang*rad2deg<<" deg,"
      <<" Ini_Sep="<<Ini_Sep<<" cm.\n"<<endl;
  }
#endif
}
  

//I need to determine the location that this hit should be inserted
//I should compare the distances of all hits behine this seed in the seed_list
//By Jixie: this will slow down the search and not sure it can improve the result 
//if the chain is full return 0, otherwise 1
int  ChainFinder::AddAHitToChain(int seed, int chainid, int hitid, double space)
{
  if (fHitNumInAChain[chainid] >= MAX_HITS_PER_CHAIN) {
    cout<<"warning: MAX_HITS_PER_CHAIN("<<MAX_HITS_PER_CHAIN
        <<") is too small, hits are potentially lost!\n";
    return 0;
  }
  
  //determine the location this hit should be incerted
  int n=fHitNumInAChain[chainid];  //just make it short
  int t=n-1;  
  //require t>1 (not 0),  because this hit can not be inserted in front of its seed
  while(t>1 && fParentSeed[t]==seed && fDist2Seed[t]>space) t--;
  
  int position = (t==n-1) ? n : t;
  //if you want to always add this seed to the end of the list, set position = n;
  position = n;
  InsertAHitToChain(chainid,hitid,position,seed,space);
  return 1;
}

//operate to fHitIDInAChain[MAX_CHAINS_PER_EVENT][MAX_HITS_PER_CHAIN];
//will not touch fChainBuf yet
void ChainFinder::InsertAHitToChain(int chainid, int hitid, int position, int seed, double separation)
{ 
  int n = fHitNumInAChain[chainid];
  if(position>n) position=n;

  //shift the hits behind the given position 
  for(int i=n;i>position;i--) { 
    fHitIDInAChain[chainid][i] = fHitIDInAChain[chainid][i-1];
    //also need to shift parent seed and space
    fParentSeed[i] = fParentSeed[i-1];
    fDist2Seed[i] = fDist2Seed[i-1];
  }
  
  fHitPool[hitid].Status |= HISUSED; //mark this hit as used
  fHitIDInAChain[chainid][position] = hitid;
  fHitNumInAChain[chainid] += 1;
  
  fParentSeed[position] = seed;
  fDist2Seed[position] = separation;

#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=4) {
    cout<<"  Chain "<<chainid<<": hit "<<hitid<<" is "<< ((position==n) ? "added":"inserted")
      <<" to "<<position<<"\n";
  }
#endif 

}

//remove the first hit with id==hitid from fHitIDInAChain[MAX_CHAINS_PER_EVENT][MAX_HITS_PER_CHAIN];
//will not touch fChainBuf yet
//return number of hit that removed
int  ChainFinder::RemoveAHitFromChain(int chainid, int hitid)
{
  //get the position
  int found = 0;
  int idx=0;  //store the found position

  int n = fHitNumInAChain[chainid];
  for(int i=0;i<n;i++) { 
    if(fHitIDInAChain[chainid][i] == hitid) {
       found = 1;
       idx = i;
       break;
    }
  }
  if(!found) return 0;
  //shift all hits in the back by 1 position
  for(int i=idx;i<n-1;i++) {
      fHitIDInAChain[chainid][i] = fHitIDInAChain[chainid][i+1];
  }
  //shrink the chain buffer
  fHitIDInAChain[chainid][n-1] = -1;
  fHitNumInAChain[chainid] -= 1;
  return 1;
}

//remove the hit at given position  from fHitIDInAChain[MAX_CHAINS_PER_EVENT][MAX_HITS_PER_CHAIN];
//will not touch fChainBuf yet
//return number of hit that removed
int  ChainFinder::RemoveAHitFromChain_At(int chainid, int position)
{
  //check the position  
  int n = fHitNumInAChain[chainid];
  if (position>=n) return 0;

  //shift all hits in the back by 1 position
  for(int i=position;i<n-1;i++) {
      fHitIDInAChain[chainid][i] = fHitIDInAChain[chainid][i+1];
  }
  //shrink the chain buffer
  fHitIDInAChain[chainid][n-1] = -1;
  fHitNumInAChain[chainid] -= 1;
  return 1;
}



//return how many hit has been found
//fid hits for a given seed(pivot)
//nhits: number of hits inside the event
int ChainFinder::SearchHitsForASeed(int seed, int seed_pre) 
{
  //declare the 3 vector here to avoid construction and deconstruction frequently
  TVector3 pV3_seed_pre(fHitPool[seed_pre].X, fHitPool[seed_pre].Y, fHitPool[seed_pre].Z);
  TVector3 pV3_seed(fHitPool[seed].X, fHitPool[seed].Y, fHitPool[seed].Z);
  TVector3 pV3(0,0,0);
  TVector3 pV3Diff, pV3Diff_pre;
     
  pV3Diff_pre = pV3_seed-pV3_seed_pre;
  
  int found = 0;
  for(int i = 0; i<fHitNum; i++) {
    if(fHitPool[i].Status & HITUNAV) continue;
    if(seed == i) continue;
    
    pV3.SetXYZ(fHitPool[i].X, fHitPool[i].Y, fHitPool[i].Z);
    pV3Diff = pV3 - pV3_seed;
    
    //check the distance    
    //for the first seed of a chain, do not require angle in range
    //but require distance < 1.0 cm
    double separation = pV3Diff.Mag(); 
    if((seed == seed_pre && separation > Ini_Sep) || separation > Max_Sep ) continue; 
    
    //check the angle between pV3_diff and pV3_diff_pre
    double angle = pV3Diff.Angle(pV3Diff_pre);
    //not very sure we need this line    
    if(angle>90./rad2deg) angle = 180./rad2deg - angle;
    
#ifdef _ChainFinderDebug_
    if(_ChainFinderDebug_>=5) {
      printf("\t seed=%-3d seed_pre=%-3d hit=%-3d: separation=%6.2f cm, angle=%7.1f deg  ",
              seed, seed_pre, i, separation, angle*rad2deg);
    }
#endif

    if((separation >  Min_Sep && angle < Max_Sep_Ang) || 
       (separation <= Min_Sep && angle < Min_Sep_Ang) ) {
      int ret = AddAHitToChain(seed,fChainNum,i,separation);
      if(ret == 0) break; // this chain is full
      found += 1;
      continue;
    } 
    
#ifdef _ChainFinderDebug_
    if(_ChainFinderDebug_>=5) {
      printf("  ...... skipped \n");
    }
#endif
  }
  return found;
}


void ChainFinder::SearchChains(int do_sort) //HitStruct *fHitPool, int nhits)
{
  //reset fChainNum = 0;  
  fChainNum = 0;
  fHitNumInAChain[fChainNum] = 0;
  
  //loop over the seed pool, currently every point could be the seed
  int seed = 0, seed_pre = 0;
  for (int anchor_hit=0; anchor_hit < fHitNum; anchor_hit++) {  //anchor_hit
    //do nothing if this hit has been marked as used or unavailable
    if ( fHitPool[anchor_hit].Status & HITUNAV )  continue;
      
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=3) {
    cout<<"\n*****anchor_hit = "<<anchor_hit<<" *****"<<endl;
  }
#endif

    //add the initial seed into this chain  
    if(fHitNumInAChain[fChainNum] == 0) {  
      seed = anchor_hit;
      AddAHitToChain(seed,fChainNum,seed,0.0);
      //fHitNumInAChain[fChainNum] will be added by 1 inside 
    }
        
    //SEARCH ALGORITHM----->
    //search a chain: looping over all found seeds
    //Note that fHitNumInAChain[fChainNum] will self increasing if a hit is added.
    //currently MAX_HITS_PER_CHAIN=200, which is the maximum number of hits that the global 
    //helix fitter will take.  
    //If necessary, MAX_HITS_PER_CHAIN could be larger during chain search such that it will
    //loop over all hits belong to this chain. However, store only 200 hits into the final 
    //chain buffer (fChainBuf). 
    for (int seed_idx=0; seed_idx<fHitNumInAChain[fChainNum]; seed_idx++) {
      //thos sanity check is not necessary any more because I have check that in SearchHitsForASeed()
      //I want to gain speed by removing this sanity check
      //if(seed_idx >= MAX_HITS_PER_CHAIN) {
      //  cout<<"warning: MAX_HITS_PER_CHAIN("<<MAX_HITS_PER_CHAIN<<") is too small, hits are potentially lost!\n";
      //  break;
      //}
      seed = fHitIDInAChain[fChainNum][seed_idx]; 
      seed_pre = fParentSeed[seed_idx];
      int found = 0;
      found = SearchHitsForASeed(seed, seed_pre);
      
      //FIX me
      //By Jixie: this part is for cross chain, I do not know how to deal with it yet.
      //if(!found) {/**/;}
    }
    
    
    //judge if this chain is valid or not
    //do not deal with those chain if any of this is true: 
    //  a) number of hits < 5, or b) Smin > kRTPC_R_GEM1-2cm, or 
    //  c) Smax < kRTPC_R_Cathode+2cm, or d) Smax-Smin<2cm
    
    bool IsValidChain = false;
    
    //determine Smin Smax,
    double Smin=9999.0, Smax=0.0;    
    int *buf = fHitIDInAChain[fChainNum];   //just to make it short    
    for(int i=0;i<fHitNumInAChain[fChainNum];i++) {
      if (Smin > fHitPool[buf[i]].S) {Smin = fHitPool[buf[i]].S;}
      if (Smax < fHitPool[buf[i]].S) {Smax = fHitPool[buf[i]].S;}  
    }
#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=3
        if(_ChainFinderDebug_>=3) {
          cout<<"  Chain "<<fChainNum<<": "<<fHitNumInAChain[fChainNum]<<" hits,  Smax="
              <<Smax<<",  Smin="<<Smin<<endl;
        }
#endif
    
    if (fHitNumInAChain[fChainNum] < Min_HITS_PER_CHAIN) {
#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=3
        if(_ChainFinderDebug_>=3) {
          cout<<"\n***Warning: ignore this chain, it contains too few hits, HitNum="
              <<fHitNumInAChain[fChainNum]<<endl;
        }
#endif
    } else if (Smin > kRTPC_R_GEM1-2.0) {
#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=3
      if(_ChainFinderDebug_>=3) {
        cout<<"\n ***Warning: ignore this chain, Smin("<<Smin<<") > kRTPC_R_GEM1-2.0"<<endl;
      }
#endif
    } else if (Smax < kRTPC_R_Cathode+2.0) {
#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=3
      if(_ChainFinderDebug_>=3) {
        cout<<"\n ***Warning: ignore this chain, Smax("<<Smax<<") < kRTPC_R_Cathode+2.0"<<endl;
      }
#endif
    } else if (Smax - Smin < 2.0) {
#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=3
      if(_ChainFinderDebug_>=3) {
        cout<<"\n ***Warning: ignore this chain, Smax - Smin = "<<Smax - Smin<<endl;
      }
#endif
    }
    else IsValidChain = true;
    
    
    //judge if this chain is valid or not, do not store it if less than 5 hits
    if(IsValidChain) {      
#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=4
      if(_ChainFinderDebug_>=4) {
        PrintAChain(fChainNum, "Before sorting");
      }
#endif
      //now, sort the chain
#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=9
      if(_ChainFinderDebug_==9) BenchmarkSort(fChainNum);
#endif
      if(do_sort) SortAChain(fChainNum);
      StoreAChain(fChainNum);
      fChainNum++;
      if (fChainNum >= MAX_CHAINS_PER_EVENT) {
        cout<<"***Warning: MAX_CHAINS_PER_EVENT("<<MAX_CHAINS_PER_EVENT<<") is too small, stop search...\n";
        break;
      }
    } else {
      //If you do not want to save this segment, clear the buffer for this chain
      bool keep_all_segment = false;
      if (keep_all_segment) {
        fChainNum++;
        if (fChainNum >= MAX_CHAINS_PER_EVENT) {
          cout<<"***Warning: MAX_CHAINS_PER_EVENT("<<MAX_CHAINS_PER_EVENT<<") is too small, stop search...\n";
          break;
        }
      } else {
        for(int jj=0;jj<fHitNumInAChain[fChainNum];jj++) {
          fHitIDInAChain[fChainNum][jj]=-1;
        }
        fHitNumInAChain[fChainNum]=0;
      }
    }
    
  }// for (int anchor_hit = 0 ......

}

  
void ChainFinder::PrintAChain(int *buf, int size, const char *keywords)
{
  int n = size;
  printf("\nChain buffer##: %4d hits: %s\n", n, keywords);

  printf("%-8s: ","hit-id");
  for (int jj=0; jj<n; jj++) {
    printf("%7d", buf[jj]);
    if (!((jj+1)%15) && jj+1!=n) printf("\n%10s","");
  }
  printf("\n");

#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=3
  if(_ChainFinderDebug_>=5) {
    printf("%-8s: ","s(cm)");
    for (int jj=0; jj<n; jj++) {
      printf("%7.2f", fHitPool[buf[jj]].S);
      if (!((jj+1)%15) && jj+1!=n) printf("\n%10s","");
    }
    printf("\n");

    printf("%-8s: ","phi(deg)");
    for (int jj=0; jj<n; jj++) {
      printf("%7.1f", fHitPool[buf[jj]].Phi*rad2deg);
      if (!((jj+1)%15) && jj+1!=n) printf("\n%10s","");
    }
    printf("\n");
  }
  if(_ChainFinderDebug_>=6) {
    printf("%-8s: ","TDC");
    for (int jj=0; jj<n; jj++) {
      printf("%7d", fHitPool[buf[jj]].TDC);
      if (!((jj+1)%15) && jj+1!=n) printf("\n%10s","");
    }
    printf("\n");
  }

  if(_ChainFinderDebug_>=3) {
    printf("%-8s:","ThrownTID");
    for (int jj=0; jj<n; jj++) {
      printf("%7d", fHitPool[buf[jj]].ThrownTID);
      if (!((jj+1)%15) && jj+1!=n) printf("\n%10s","");
    }
    printf("\n");
  }
#endif

  printf("\n");
}

//print the chain buffer
void ChainFinder::PrintAChain(int chainid, const char *keywords)
{
  int n = fHitNumInAChain[chainid];
  printf("\nChain #%2d: %4d hits: %s\n", chainid, n, keywords);

  printf("%-8s: ","hit-id");
  for (int jj=0; jj<n; jj++) {
    printf("%7d", fHitIDInAChain[chainid][jj]);
    if (!((jj+1)%15) && jj+1!=n) printf("\n%10s","");
  }
  printf("\n");

#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=3
  if(_ChainFinderDebug_>=5) {
    printf("%-8s: ","s(cm)");
    for (int jj=0; jj<n; jj++) {
      printf("%7.2f", fHitPool[fHitIDInAChain[chainid][jj]].S);
      if (!((jj+1)%15) && jj+1!=n) printf("\n%10s","");
    }
    printf("\n");

    printf("%-8s: ","phi(deg)");
    for (int jj=0; jj<n; jj++) {
      printf("%7.1f", fHitPool[fHitIDInAChain[chainid][jj]].Phi*rad2deg);
      if (!((jj+1)%15) && jj+1!=n) printf("\n%10s","");
    }
    printf("\n");
  }

  if(_ChainFinderDebug_>=6) {
    printf("%-8s: ","TDC");
    for (int jj=0; jj<n; jj++) {
      printf("%7d", fHitPool[fHitIDInAChain[chainid][jj]].TDC);
      if (!((jj+1)%15) && jj+1!=n) printf("\n%10s","");
    }
    printf("\n");
  }
 
  if(_ChainFinderDebug_>=3) {
    printf("%-8s:","ThrownTID");
    for (int jj=0; jj<n; jj++) {
      printf("%7d", fHitPool[fHitIDInAChain[chainid][jj]].ThrownTID);
      if (!((jj+1)%15) && jj+1!=n) printf("\n%10s","");
    }
    printf("\n");
  }
#endif

  printf("\n");
}

//Bubble sort is simple but slow, the number of step is in the ordre of O(n^2)
//sort fHitIDInAChain[][], by S increaseing order
void ChainFinder::BubbleSort_S(int *arr, int size)  
{   
  int i, j, tmp;

  for(i=0; i<size; i++) { // Make a pass through the array for each element
    for(j=1; j<(size-i); j++) { // Go through the array beginning to end
      if(fHitPool[arr[j-1]].S > fHitPool[arr[j]].S) {
	// If the the previous number is greater, swap it 
	tmp = arr[j-1];
	arr[j-1] = arr[j];
	arr[j] = tmp;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//Quick Sort Algorithm
//The divide-and-conquer strategy is used in quicksort. Below the recursion step is described:
//
//1. Choose a pivot value. We take the value of the middle element as pivot value, but it 
//   can be any value, which is in range of sorted values, even if it doesn't present in the array.
//2. Partition. Rearrange elements in such a way, that all elements which are lesser than 
//   the pivot go to the left part of the array and all elements greater than the pivot, go 
//   to the right part of the array. Values equal to the pivot can stay in any part of the array. 
//   Notice, that array may be divided in non-equal parts.
//3. Sort both parts. Apply quicksort algorithm recursively to the left and the right parts.
//
//Partition algorithm in detail
//
//  There are two indices i and j and at the very beginning of the partition algorithm i points 
//  to the first element in the array and j points to the last one. Then algorithm moves i 
//  forward, until an element with value greater or equal to the pivot is found. Index j is moved 
//  backward, until an element with value lesser or equal to the pivot is found. If i ¡Ü j then 
//  they are swapped and i steps to the next position (i + 1), j steps to the previous one (j - 1). 
//  Algorithm stops, when i becomes greater than j.
//
//  After partition, all values before i-th element are less or equal than the pivot and all 
//  values after j-th element are greater or equal to the pivot.
/////////////////////////////////////////////////////////////////////////////////////////////////////
//quick sort faster than bubble sort, the number of step is in the ordre of O(n*log(n))
//sort fHitIDInAChain[][], by S increaseing order
void ChainFinder::QuickSort_S(int *arr, int left, int right) 
{
  int i = left, j = right;
  int tmp;
  int pivot = arr[(left + right) / 2];
  double S_pivot = fHitPool[pivot].S;

  // partition 
  while (i <= j) {
    while (fHitPool[arr[i]].S < S_pivot) i++;
    while (fHitPool[arr[j]].S > S_pivot) j--;

    //swap i and j
    if (i <= j) {
      tmp = arr[i];
      arr[i] = arr[j];
      arr[j] = tmp;
      i++;
      j--;
    }
  };

  // recursion 
  if (left < j) QuickSort_S(arr, left, j);
  if (i < right) QuickSort_S(arr, i, right);

}

//combile sorting S and Phi together is very slow!!!
//sort fHitIDInAChain[][], by S increaseing order, if S equal, by phi increasing order
//I could also use TDC instead of S, but I worry about the TDC offsets various from chanel to channel
void ChainFinder::QuickSort_SPhi(int *arr, int left, int right) 
{
  int i = left, j = right;
  int tmp;
  int pivot = (left + right) / 2;
  //double S_pivot = fHitPool[arr[pivot]].S;
  double tmpphi;

  
  //make a copy of the phi into array 
  //if cross the 180deg line, need to add 2Pi to negative phi 

  //determine Phimin Phimax 
  double Phimin,Phimax;
  int idxPhimin,idxPhimax;
  double *phi = new double [right+1];
  idxPhimin=idxPhimax=left;
  Phimin=Phimax=fHitPool[arr[left]].Phi;
  for(int i=left;i<=right;i++) {
    phi[i] = fHitPool[arr[i]].Phi;
    if(Phimin >= fHitPool[arr[i]].Phi) {Phimin = fHitPool[arr[i]].Phi; idxPhimin=i;}
    else if(Phimax <= fHitPool[arr[i]].Phi) {Phimax = fHitPool[arr[i]].Phi; idxPhimax=i;}
  }
  
  //this track aross the 180 deg line, need to check for phimin phimax again
  if(Phimax-Phimin > kPi) {   
  for(int i=left;i<=right;i++) {
      if (phi[i]<0.0) phi[i] += 2*kPi;
      if(Phimax <= phi[i]) {Phimax = phi[i]; idxPhimax=i;}   
      else if(Phimin >= phi[i]) {Phimin = phi[i]; idxPhimin=i;}
    }
  }

  //double phi_pivot = phi[pivot];

  //for (int t = left; t <= right; t++) {
  // printf("phi[t]=%f  fHitPool[arr[t]].Phi=%f \n",phi[t]*rad2deg, fHitPool[arr[t]].Phi*rad2deg);
  //}

  // partition 
  while (i <= j) {
    //while ( fHitPool[arr[i]].S < S_pivot || 
    //  (fabs(fHitPool[arr[i]].S-S_pivot)<1.0E-5 && phi[i]<phi_pivot) ) i++;
    //while ( fHitPool[arr[j]].S > S_pivot || 
    //  (fabs(fHitPool[arr[j]].S-S_pivot)<1.0E-5 && phi[j]>phi_pivot) ) j--;
    //Note that I have to use floating phi_pivot and S_pivot, the above will crash

    while ( fHitPool[arr[i]].S < fHitPool[arr[pivot]].S || 
      (fabs(fHitPool[arr[i]].S-fHitPool[arr[pivot]].S)<1.0E-5 && phi[i]<phi[pivot]) ) i++;
    while ( fHitPool[arr[j]].S > fHitPool[arr[pivot]].S || 
      (fabs(fHitPool[arr[j]].S-fHitPool[arr[pivot]].S)<1.0E-5 && phi[j]>phi[pivot]) ) j--;
    

    //swap i and j
    if (i <= j) {
      //cout<<"\QuickSort_SPhi(): swap("<<i<<", "<<j<<") ..."<<endl;
      //need to swap phi array too
      tmpphi = phi[i];
      phi[i] = phi[j];
      phi[j] = tmpphi;
      tmp = arr[i];
      arr[i] = arr[j];
      arr[j] = tmp;

      i++;
      j--;
    }
  };

  // recursion 
  if (left < j) QuickSort_SPhi(arr, left, j);
  if (i < right) QuickSort_SPhi(arr, i, right);

  delete phi;
  //Stop4Debug();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//Selection Sort Algorithm
//
// Array is imaginary divided into two parts - sorted one and unsorted one. At the beginning, 
// sorted part is empty, while unsorted one contains whole array. At every step, algorithm finds 
// minimal element in the unsorted part and adds it to the end of the sorted one. When unsorted 
// part becomes empty, algorithm stops.
//
// When algorithm sorts an array, it swaps first element of unsorted part with minimal element 
// and then it is included to the sorted part. This implementation of selection sort in not stable. 
// In case of linked list is sorted, and, instead of swaps, minimal element is linked to the 
// unsorted part, selection sort is stable.
/////////////////////////////////////////////////////////////////////////////////////////////////////
//selection sort is slow, the number of step is in the ordre of O(n^2)
//sort fHitIDInAChain[][], by S increaseing order
void ChainFinder::SelectSort_S(int *arr, int size) 
{
  int i, j, minIndex, tmp;    

  for (i = 0; i < size - 1; i++) {
    minIndex = i;

    for (j = i + 1; j < size; j++) {
      //find the minimum for the rest of the array 
      if (fHitPool[arr[j]].S < fHitPool[arr[minIndex]].S) minIndex = j;
    }

    if (minIndex != i) {
      //swap
      tmp = arr[i];
      arr[i] = arr[minIndex];
      arr[minIndex] = tmp;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//Insertion Sort Algorithm
//
// Insertion sort algorithm somewhat resembles selection sort. Array is imaginary divided into two 
// parts - sorted one and unsorted one. At the beginning, sorted part contains first element of the 
// array and unsorted one contains the rest. At every step, algorithm takes first element in the 
// unsorted part and inserts it to the right place of the sorted one. When unsorted part becomes 
// empty, algorithm stops.
/////////////////////////////////////////////////////////////////////////////////////////////////////
//insertion sort is slow, the number of step is in the ordre of O(n^2)
//sort fHitIDInAChain[][], by S increaseing order
void ChainFinder::InsertSort_S(int *arr, int size)
{
  int i, j, tmp;

  for (i = 1; i < size; i++) {
    j = i;
    while (j > 0 && fHitPool[arr[j]].S < fHitPool[arr[j-1]].S) {
      tmp = arr[j];
      arr[j] = arr[j-1];
      arr[j-1] = tmp;
      j--;
    }
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//Shell Sort Algorithm
//
// Insertion sort is quick if size is small and the array is nearly sorted. It is slow when the array
// is sorted in the opposite way.
// Shell sort is based on insertion sort. It repeatedly compares elements that are a certain distance 
// away from each other (gap represents this distance)
//
/////////////////////////////////////////////////////////////////////////////////////////////////////
//Shell Sort using half division
//Sort array fHitIDInAChain[][], by S increaseing order
void ChainFinder::ShellSort2_S(int *arr, int size)
{
  //Narrow the array by 2 everytime
  for (int gap = size / 2; gap > 0; gap /= 2) {
    for (int i = gap; i < size; ++i) {
      for (int j = i; j >= gap && fHitPool[arr[j]].S < fHitPool[arr[j-gap]].S; j -= gap) {
      //for (int j = i; j >= gap && arr[j] < arr[j-gap]; j -= gap) {
	//swap j and j-gap
	int tmp = arr[j];
	arr[j] = arr[j - gap];
	arr[j - gap] = tmp;
      }
    } 
  }
}


//Shell Sort using gap_sequence of 13,9,5,2,1
//Sort array fHitIDInAChain[][], by S increaseing order
void ChainFinder::ShellSort_Seq_S(int *arr, int size)
{
  //this is even worse than the next 
  static const int gap_sequence[] = { 13, 9, 5, 2, 1 }; 
  for(int gg=0;gg<5;gg++) {
    int gap = gap_sequence[gg];
    if( gap < size )
    {
      for(int i = gap ; i < size; ++i )
	for (int j = i-gap; j >= 0 && fHitPool[arr[j]].S > fHitPool[arr[j+gap]].S; j -= gap) {
	  int tmp = arr[j];
	  arr[j] = arr[j + gap];
	  arr[j + gap] = tmp;
	}
    }
  }
 
}


//sort fHitIDInAChain[][], by Phi increaseing order
void  ChainFinder::InsertSort_Phi(int *arr, int size) 
{
  int i, j, tmp;
  double tmpphi;
  //need to check if this check arossing phi=180deg line
  //Here I assumes that this chain is already sorted by S
  //therefore its [0] and [size-1] represent the two ends of phi
  double dphi = fHitPool[arr[0]].Phi - fHitPool[arr[size-1]].Phi;
  bool arossline = (fabs(dphi) > 3*kPi/4.0) ? true : false;

  //make a copy of the phi into array 
  //if cross the line, need to add 2Pi to negative phi 
  double *phi = new double [size];
  for (i = 0; i < size; i++) {
    phi[i] = fHitPool[arr[i]].Phi; 
    if (arossline && phi[i] < 0.0) phi[i] += 2*kPi; 
  }   

#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=11
  //print the array before shorting
      if(_ChainFinderDebug_>=13)
        cout<<"\tInsertSort_Phi(): before sorting ...**************************"<<endl;
      
      if(_ChainFinderDebug_>=13) {
	printf("%9s:","TDC");
	for(int t=0;t<size;t++) { 
	  printf("%7d",fHitPool[arr[t]].TDC);
	  if(!((t+1)%15)) printf("\n%10s","");
	}
	printf("\n");
	printf("%9s:","phi_deg");
	for(int t=0;t<size;t++) { 
	    printf("%7.1f",phi[t]*57.3); 
	  if(!((t+1)%15)) printf("\n%10s","");
	}
	printf("\n");
	printf("%9s:","ThrownTID");
	for(int t=0;t<size;t++) { 
	    printf("%7d",fHitPool[arr[t]].ThrownTID); 
	  if(!((t+1)%15)) printf("\n%10s","");
	}
	printf("\n");
      }
#endif

  //do sorting now 
  for (i = 1; i < size; i++) {
    j = i;
    while (j > 0 && fHitPool[arr[j]].TDC == fHitPool[arr[j-1]].TDC &&
           phi[j] < phi[j-1]) {
      //swap
      tmp = arr[j];
      arr[j] = arr[j-1];
      arr[j-1] = tmp;
      //need to swap phi also
      tmpphi = phi[j];
      phi[j] = phi[j-1];
      phi[j-1] = tmpphi;
	     
#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=10
      if(_ChainFinderDebug_>=11)
        cout<<"\tInsertSort_Phi(): swap("<<j<<", "<<j-1<<") ..."<<endl;
      
      if(_ChainFinderDebug_>=12) {
	printf("%9s:","TDC");
	for(int t=0;t<=j;t++) {
	  //if(t==j-1 || t==j) printf("%7d",fHitPool[arr[t]].TDC);
	  //else printf("%7s",".");
	  printf("%7d",fHitPool[arr[t]].TDC);
	  if(!((t+1)%15)) printf("\n%10s","");
	}
	printf("\n");
	printf("%9s:","phi_deg");
	for(int t=0;t<=j;t++) {
	  //if(t==j-1 || t==j) printf("%7.1f",phi[t]*57.3);
	  //else printf("%7s","."); 
	  printf("%7.1f",phi[t]*57.3);
	  if(!((t+1)%15)) printf("\n%10s","");
	}
	printf("\n");
      }
#endif

      j--;
    }
  }
  //free memory
  delete phi;
}

//sort fHitIDInAChain[][], will use it to store chain
//First sort by S increasing order, if their S is equal, sort by phi.
//in order of phi could be increasing or decreasing. But it should be 
//consistant with the whole track
void ChainFinder::BenchmarkSort(int chainid)
{
  //now sort array buf using its S and Phi 
#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=9
  //By Jixie; I want to do a bench mark test which sort algrithm works better 
  //Conclusion of sorting:
  //1) Bubble and selection sort are most slow in any case
  //2) Bubble, selection and incertion sort are stable (stable means will not change 
  //order if two values are equal to each other).  Quick and shell sort are not stable.
  //3) Incertion sort is the fastest only if the array already sorted. It will be as 
  //slow as bubble if the array is sorted in opposite way. For random array it is still 
  //faster than bubble and selection sort. 
  //4) Shell_2 and Shell_seq are about the same speed. 
  //5) For random array, if array size>=30, quick sort is the fastest one. 
  //6) For the chain from this chainfinder, (they are almost sorted only a few location 
  //need to adjust), if array size>50, quick sort is the fastest, otherwise incertion sort 
  //is the best. 

  int buf0[2000];
  int n=0; 
  bool test_by_one_chain = true;  // false will use the whole hit pool to do this test
  if(test_by_one_chain) {
    n=fHitNumInAChain[chainid]; 
    for(int i=0;i<n;i++) buf0[i]=fHitIDInAChain[chainid][i];
  }else{
    n=fHitNum; 
    for(int i=0;i<n;i++) buf0[i]=i;
  }
  //std::random_shuffle( buf0, buf0+n ) ;

  TBenchmark pBenchmark;
  int buf1[2000];

  /////////////////////////////////////////////////
  pBenchmark.Start("stress_copybuf");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf1[i]=buf0[i];
  }
  pBenchmark.Stop("stress_copybuf");
  pBenchmark.Print("stress_copybuf");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");

  /////////////////////////////////////////////////
  pBenchmark.Start("stress_bubble");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf1[i]=buf0[i];
    BubbleSort_S(buf1,n);
    InsertSort_Phi(buf1,n);
  }
  pBenchmark.Stop("stress_bubble");
  pBenchmark.Print("stress_bubble");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");

  /////////////////////////////////////////////////
  pBenchmark.Start("stress_select");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf1[i]=buf0[i];
    SelectSort_S(buf1,n);
    InsertSort_Phi(buf1,n);
  }
  pBenchmark.Stop("stress_select");
  pBenchmark.Print("stress_select");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");

  /////////////////////////////////////////////////
  pBenchmark.Start("stress_insert");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf1[i]=buf0[i];
    InsertSort_S(buf1,n);
    InsertSort_Phi(buf1,n);
  }
  pBenchmark.Stop("stress_insert");
  pBenchmark.Print("stress_insert");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");

  /////////////////////////////////////////////////
  pBenchmark.Start("stress_quick");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf1[i]=buf0[i];
    QuickSort_S(buf1,0,n-1);
    InsertSort_Phi(buf1,n);
  }
  pBenchmark.Stop("stress_quick");
  pBenchmark.Print("stress_quick");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");

  /////////////////////////////////////////////////
  pBenchmark.Start("stress_shell2");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf1[i]=buf0[i];
    ShellSort2_S(buf1,n);
    InsertSort_Phi(buf1,n);
  }
  pBenchmark.Stop("stress_shell2");
  pBenchmark.Print("stress_shell2");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");
  
  /////////////////////////////////////////////////
  pBenchmark.Start("stress_shell_seq");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf1[i]=buf0[i];
    ShellSort_Seq_S(buf1,n);
    InsertSort_Phi(buf1,n);
  }
  pBenchmark.Stop("stress_shell_seq");
  pBenchmark.Print("stress_shell_seq");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");
  
  /////////////////////////////////////////////////
  pBenchmark.Start("stress_quick_sphi");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf1[i]=buf0[i];
    QuickSort_SPhi(buf1,0,n-1);
  }
  pBenchmark.Stop("stress_quick_sphi");
  pBenchmark.Print("stress_quick_sphi");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");
  
  
  /////////////////////////////////////////////////
  int buf2[2000];
  for(int i=0;i<n;i++) buf1[i]=buf0[i];
  InsertSort_S(buf1,n);
  pBenchmark.Start("stress_InsertSort_Phi");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf2[i]=buf1[i];
    InsertSort_Phi(&buf2[0],n);
  }
  pBenchmark.Stop("stress_InsertSort_Phi");
  pBenchmark.Print("stress_InsertSort_Phi");
  for(int i=0;i<n;i++) printf("%3d ",buf2[i]);
  printf("\n\n");

  /////////////////////////////////////////////////
  pBenchmark.Start("stress_my_choice");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf1[i]=buf0[i];
    if(n<50) InsertSort_S(buf1,n);
    else QuickSort_S(buf1,0,n-1);
    InsertSort_Phi(buf1,n);
  }
  pBenchmark.Stop("stress_my_choice");
  pBenchmark.Print("stress_my_choice");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");

  
  /////////////////////////////////////////////////
  for(int i=0;i<n;i++) buf1[i]=buf0[i];
  pBenchmark.Start("stress_SortAChain");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf0[i]=buf1[i];
    SortAChain(chainid);
  }
  pBenchmark.Stop("stress_SortAChain");
  pBenchmark.Print("stress_SortAChain");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");
  for(int i=0;i<n;i++) buf0[i]=buf1[i];

#endif
}


//For froward tracks, sorting by either S then Phi will work.
//For curve back tracks, sorting by S then phi will not work.
//Here is my solution:  
//1) spilt the whole chain into two parts: forward part + backward part
//2) sort both parts by S then phi. 
//   Forward part: S increasing, if S equal then phi increasing;
//   Backward part: S decreasing, if S equal then phi increasing;
//3) Finally merge these 2 parts together 
void ChainFinder::SortAChain(int chainid) 
{  
  int nhits = fHitNumInAChain[chainid];
  int *buf = fHitIDInAChain[chainid];
  
  double *phi = new double [nhits];
  for(int i=0;i<nhits;i++) {phi[i]=fHitPool[buf[i]].Phi;}

  //first, we judge if it is a curve back track
  //if Smax and phi_min|phi_max is the same hit or neighbor hits
  //it is a forward track

  //determine Smin Smax Phimin Phimax
  double Smin,Smax,Phimin,Phimax;
  int idxSmin,idxSmax,idxPhimin,idxPhimax;
  
  idxSmin=idxSmax=idxPhimin=idxPhimax=0;
  Smin=9999.0; Smax=0.0;
  Phimin=9999.0;Phimax=-9999.0;
  for(int i=0;i<nhits;i++) {
    if(Smin >= fHitPool[buf[i]].S) {Smin = fHitPool[buf[i]].S; idxSmin=i;}
    if(Smax - fHitPool[buf[i]].S < 1.0E-5) {Smax = fHitPool[buf[i]].S; idxSmax=i;}   
    if(Phimin >= fHitPool[buf[i]].Phi) {Phimin = fHitPool[buf[i]].Phi; idxPhimin=i;}
    if(Phimax - fHitPool[buf[i]].Phi < 1.0E-5) {Phimax = fHitPool[buf[i]].Phi; idxPhimax=i;}
  }
  
  //this track aross the 180 deg line, need to check for phimin phimax again
  if(Phimax-Phimin > kPi) {
    Phimin=9999.0;
    Phimax=-9999.0;
    for(int i=0;i<nhits;i++) { 
      if (phi[i]<0.0) phi[i] += 2*kPi;
      if(Phimax <= phi[i]) {Phimax = phi[i]; idxPhimax=i;}   
      if(Phimin >= phi[i]) {Phimin = phi[i]; idxPhimin=i;}
    }
  }

  //determine track property: forward or backward
  //for positive track phi increasing
  //for backward track, Smax locate in the middle of the chain
  //but this chain buf is not sorted yet. so we have to sort it first ......
  bool bIsBackwardTrack=false;  
  //first of all, check coverage:
  double Sspan=Smax-Smin;
  double Phispan=Phimax-Phimin;
  //Theoretically, Phispan = Phi_at_cathode - Phi_at_GEM1, for 5T field and 4cm driftdisctance 
  //Phi_at_S = kPi/2-atan(S/2./R), where R = P_GeV/(0.3B) 
  //For Pt=0.30 GeV, Phispan = 5.7 deg
  //For Pt=0.20 GeV, Phispan = 8.3 deg
  //For Pt=0.15 GeV, Phispan = 10.8 deg
  //For Pt=0.053 GeV, Phispan = 66.8 deg
  //If Phispan is too small, it is a high Pt chain, which never curve back
  //If Phispan is very large, it is a small Pt chain, which is very likely to be a curve back.
  //I use 80 deg here to set a curve back track. Please note that this value depends on the
  //drift-path of the electron, we need to verify with garfield. 
  //Note that the phi uncertainty is about 2 deg
  
  if(Phispan < 12.0/rad2deg ) bIsBackwardTrack=false;
  //Smax and phimax stay together, it is a positive and forward track 
  else if (Smax-fHitPool[buf[idxPhimax]].S < 0.5) bIsBackwardTrack=false;
  //Smax and phimin stay together, it is a negative and forward track 
  else if (Smax-fHitPool[buf[idxPhimin]].S < 0.5) bIsBackwardTrack=false; 
  //take the differece of S between two ends of phi, if it is >0.5 cm less than Sspan, 
  //it must be a curve-back track. For a curve-back track, it is hard to determine which
  //end is head or tail
  else if (Sspan - fabs(fHitPool[buf[idxPhimin]].S-fHitPool[buf[idxPhimax]].S) > 0.5 ) {
    bIsBackwardTrack=true;
  }
  else if (Phispan > 80/rad2deg) bIsBackwardTrack=true;

#ifdef _ChainFinderDebug_
    if(_ChainFinderDebug_>=10) {
      cout<<"SortAChain():: overall information: nhits="<<nhits<<" idxSmax="<<idxSmax<<"  idxSmin="<<idxSmin
          <<"  idxPhimax="<<idxPhimax<<"  idxPhimin="<<idxPhimin <<endl
          <<"\t Smax="<<Smax<<"  Smin="<<Smin<<"  S[idxPhimax]="<<fHitPool[buf[idxPhimax]].S
          <<"  S[idxPhimin]="<<fHitPool[buf[idxPhimin]].S<<endl
          <<"\t Phimax="<<Phimax*rad2deg<<"  Phimin="<<Phimin*rad2deg<<"  Phi[idxSmax]="<<phi[idxSmax]*rad2deg
          <<"  Phi[idxSmin]="<<phi[idxSmin]*rad2deg<<endl
          <<endl;
      printf("***************  Chain #%02d is a %s chain  ***************\n\n",
              fChainNum, (bIsBackwardTrack ? "backward":"forward") );
    }
#endif

  //now do sorting
  if (!bIsBackwardTrack) { 
    //forward track
    //1st: sort by S
    if(nhits<50) InsertSort_S(buf,nhits);
    else QuickSort_S(buf,0,nhits-1);

#ifdef _ChainFinderDebug_
    if(_ChainFinderDebug_>=11) {
      printf("*****  Forward chain: after sortting by S: *****\n");
      PrintAChain(fChainNum, "forward chain: after sorting by S");
    }
#endif

    //2nd: sort by phi increasing for equal S hits
    InsertSort_Phi(buf,nhits);

#ifdef _ChainFinderDebug_
    if(_ChainFinderDebug_>=10) {
      printf("*****  Forward Chain: after sortting by S then Phi: *****\n");
      PrintAChain(fChainNum, "Forward Chain: after sorting by S then Phi");
    }
#endif

  } else {
    //backward tracks
    //1st: split the array into two part
    //forward: [0-idxSmax], backward: (idxSmax+1,nhits-1]

    //due to the fact that the original chain is in random order
    //the backward could contain some forward hits, I should split them out
    
    int count=0;
    int ii = 0, jj = nhits-1;
    int tmp;
    double tmpphi;
    double phi_pivot = phi[idxSmax];

    // move all hits that have larger phi to the backward part
    // if it is negative chain, we should do in the opposite way, but 
    // we do not know the charge yet ...
    while (ii <= jj) {
      while ( phi[ii] < phi_pivot) ii++;
      while ( phi[jj] > phi_pivot) jj--;

      if (ii <= jj) {
        if(ii!=jj) {
#ifdef _ChainFinderDebug_
          if(_ChainFinderDebug_>=11) {
            cout<<"\tSplit backward chain: swap("<<ii<<", "<< jj<<")"<<endl;
          }
#endif
          tmp = buf[ii];  buf[ii] = buf[jj];  buf[jj] = tmp;
          tmpphi = phi[ii]; phi[ii] = phi[jj]; phi[jj] = tmpphi;
          count++;
        }
        ii++;
        jj--;
      }
    }
    idxSmax = jj;
#ifdef _ChainFinderDebug_
    if(_ChainFinderDebug_>=11 && count) {
      cout<<"\n  Backward chain: after spliting, idxSmax="<<idxSmax<<endl;
      PrintAChain(buf, nhits, "Backward chain: after splitting and before sorting");
    }
#endif

    //take care of the forward part 
    //sorting by S inceasing order
    if(idxSmax<50) InsertSort_S(buf,idxSmax+1);
    else QuickSort_S(buf,0,idxSmax);   
    //sorting by phi increasing order
    InsertSort_Phi(buf,idxSmax+1);

#ifdef _ChainFinderDebug_
    if(_ChainFinderDebug_>=11) {
      printf("*****Forward part of the backward chain: after sortting by S then Phi: *****\n");
      PrintAChain(buf, idxSmax+1, "Backward chain: after farward part sorted");
    }
#endif

    /////////////////////////////////////////////////////////////////
    //take care of the backward part
    int n1 = nhits-idxSmax-1; 
    if(n1>0) {
      int *buf1 = new int [n1];
      //copy the backward part in oppsite order since they are nearly in S decreasing order
      //now buf1 is in S nearlly increasing order
      for(int i=0;i<n1;i++) {buf1[i]=buf[nhits-1-i];}

      //sorting by S inceasing order
      if(n1<50) InsertSort_S(buf1,n1);
      else QuickSort_S(buf1,0,n1-1);  

      //now put it back as S decreasing, noting that buf1 is in S increasing order now 
      for(int i=0;i<n1;i++) buf[idxSmax+1+i]=buf1[n1-1-i];
      //sorting by phi increasing order, add Smax point into it
      InsertSort_Phi(&buf[idxSmax],n1+1);

      delete buf1;

  #ifdef _ChainFinderDebug_
      if(_ChainFinderDebug_>=11) {
	printf("*****Backward part of the backward chain: after sortting by S then Phi: *****\n");
	PrintAChain(&buf[idxSmax], n1+1, "Backward chain: after backward part sorted");
      }
  #endif

  #ifdef _ChainFinderDebug_
      if(_ChainFinderDebug_>=10) {
	printf("*****Backward chain: after sortting by S then Phi: *****\n");
	PrintAChain(&buf[0], nhits, "Backward chain: after forward+backward sorted");
      }
  #endif
    }//end of if(n1>0) {
  }
  delete phi;
}

void ChainFinder::StoreAChain(int chainid)
{
  //store only MAX_HITS_PER_TRACK of hits into the primary chain buffer
  int n = fHitNumInAChain[chainid];
  if(n>MAX_HITS_PER_TRACK) n=MAX_HITS_PER_TRACK;
  
  for (int jj=0; jj<n; jj++){
    //Fill the pointer of found chain
    int hitid=fHitIDInAChain[chainid][jj];
    fChainBuf[chainid].Hits[jj] = &(fHitPool[hitid]);
   
    //By Jixie:  I add ChainInfo to tell ChainIndex cc and HitIndex jj
    //ChainInfo  = ccccjjjj,  where ccc is ChainIndex and jj is HitIndex 
    fHitPool[hitid].ChainInfo = chainid*1.0E4 + jj;
  }

  fChainBuf[chainid].HitNum = n; //number of hits in the chain
  fChainBuf[chainid].ID = chainid;  //this chain index 
  fChainNum_Stored++;
  
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=3) {
    printf("\n ***** Store chain #%d: %d/%d hits ***** \n", chainid,n,fHitNumInAChain[chainid]);
    char strmsg[255];
    sprintf(strmsg,"stored into fChainBuf[%d]!!!",fChainNum_Stored-1);
    PrintAChain(chainid,strmsg);
  }
#endif

}

void ChainFinder::DrawPool()
{
  if (!gPad || !fHitNum) return;
  gPad->cd();

  TPolyMarker3D *pm3dp = new TPolyMarker3D(fHitNum);
  pm3dp->SetBit(kCanDelete);
  pm3dp->SetMarkerColor(40);
  pm3dp->SetMarkerStyle(1);

  for (int hh=0;hh<fHitNum;hh++) { 
    pm3dp->SetPoint(hh, fHitPool[hh].X, fHitPool[hh].Y, fHitPool[hh].Z);
  }
  pm3dp->Draw("same");
  gPad->Update();
}

void ChainFinder::DrawChain()
{
  if (!gPad || !fChainNum_Stored) return;
  gPad->cd();

  TPolyMarker3D *pm3dp = 0;
  for (int cc=0;cc<fChainNum_Stored;cc++) {    
    //FixMe:  should there be memory leak here?  
    //TPolyMarker3D is created for every chain, but nowhere to clean it up
    pm3dp = new TPolyMarker3D(fChainBuf[cc].HitNum);
    pm3dp->SetBit(kCanDelete);
    pm3dp->SetMarkerColor(cc+3);
    pm3dp->SetMarkerStyle(6);

    //identify the true track whose tid==0 then change its color and marker style
    double tid=0; 
    for(int hh=0;hh<fChainBuf[cc].HitNum;hh++) {
      pm3dp->SetPoint(hh, fChainBuf[cc].Hits[hh]->X, fChainBuf[cc].Hits[hh]->Y, 
        fChainBuf[cc].Hits[hh]->Z);
      tid += int(fChainBuf[cc].Hits[hh]->ThrownTID/10000.0);  
    }
    tid /= fChainBuf[cc].HitNum;
    if(tid<0.01) {
      pm3dp->SetMarkerColor(1);
      pm3dp->SetMarkerStyle(2);
    }
    
    pm3dp->Draw("same");
    gPad->Update();
  }
}

