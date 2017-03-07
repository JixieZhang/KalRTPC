//Originally from Howard and Carlos
//Jixie change this code into a class
//
#include <iostream>
#include <math.h>
#include <stdio.h>
#include "TVector3.h"

#include "ChainFinder.hh"

extern const double kRTPC_R_GEM1 = 7.0;
extern const double kRTPC_R_Cathode = 3.0;

///////////////////////////////////////////////////////////////////////
#define _ChainFinderDebug_ 4

//#ifdef _ChainFinderDebug_
//#include "GlobalDebuger.hh"
//#endif

///////////////////////////////////////////////////////////////////////

/*
//this is an example how to use this class 
void EXample()
{
  ChainFinder pTF;  
  pTF->SetParameters(v1,v2,v3,v4);
  for(int ii;ii<nevents;ii++) {
    pTF.Reset();
    pTF.PrepareHitPool(vvv);
    pTF.RemoveBadHits();
    pTF.SearchChains();
    //call KF
    for(int i=0;i<fChainNum;i++) DoKalmanFilter();  
  }
}
*/

using namespace std;

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
  Max_Link_Sep = 1.1;             //cm
  Max_Ang      = 30.0/rad2deg;    //rad
  Min_Ang      = 23.3/rad2deg;    //rad
  Ang_Sep      = 0.4;             //cm
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
  fHitNum -= 1;
  return found;
}

//This is for changing status for hits before search
//For example, if we want to remove bad hits or 
int  ChainFinder::RemoveBadHitsFromPool()
{
  return 0;
}

void ChainFinder::SetParameters(double space, double min_ang, double max_ang, double ang_sep)
{
  // USER SHOULD PROVIDE THESE VARIABLES
  Max_Link_Sep = space;
  Min_Ang      = min_ang;
  Max_Ang      = max_ang;
  Ang_Sep      = ang_sep;

#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=1) {
  cout<<"\nChainFinder Parameters: space="<<space<<" cm, min_ang="<<min_ang*rad2deg
      <<" deg,  max_ang="<<max_ang*rad2deg<<" deg; ang_sep="<<ang_sep<<" cm\n"<<endl;
  }
#endif
}
  
void ChainFinder::AddAHitToChain(int chainid, int hitid)
{
  if (fHitNumInAChain[chainid] >= MAX_HITS_PER_CHAIN) {        
    printf("Too many hits for the chain list. Skip...\n"); 
  }
    
  //THE HIT INDEX WHICH ACOMPLISH THE CONDITION, IS STORED HERE-->
  fHitIDInAChain[chainid][fHitNumInAChain[chainid]] = hitid;
  
  /* mark it as used */
  fHitPool[hitid].Status |= HISUSED; //this kind of assignment is for 'historical' reasons
  fHitNumInAChain[chainid]++;      //ADD HIT TO THE CHAIN, INCREASE INDEX
        
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=4) {
    cout<<"  Chain "<<chainid<<": Hit "<<hitid<<" is added\n";
  }
#endif 
}

//operate to fHitIDInAChain[MAX_CHAINS_PER_EVENT][MAX_HITS_PER_CHAIN];
//will not touch fChainBuf yet
void ChainFinder::InsertAHitToChain(int chainid, int hitid, int position)
{
  int n = fHitNumInAChain[chainid];
  if(position>n) position=n;

  for(int i=n;i>position;i--) { 
    fHitIDInAChain[chainid][i] = fHitIDInAChain[chainid][i-1];
  }
  fHitIDInAChain[chainid][position] = hitid;
  fHitNumInAChain[chainid] += 1;
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


//remove redundate hit from fHitIDInAChain[MAX_CHAINS_PER_EVENT][MAX_HITS_PER_CHAIN];
void ChainFinder::RemoveRedundantFromChain(int chainid)
{
  
  
}

//return how many hit has been found
//fid hits for a given seed(pivot)
//nhits: number of hits inside the event
int ChainFinder::SearchHitsForASeed(int seed, int seed_pre) 
{
  //declare the 3 vector here to avoid construction and deconstruction frequently
  TVector3 pV3_seed_pre(fHitPool[seed_pre].X, fHitPool[seed_pre].Y, fHitPool[seed_pre].Z);
  TVector3 pV3_seed(fHitPool[seed].X, fHitPool[seed].Y, fHitPool[seed].Z);
  TVector3 pV3_pre(0,0,0);
  TVector3 pV3(0,0,0);
  TVector3 pV3Diff, pV3Diff_pre;
     
  pV3Diff_pre = pV3_seed-pV3_seed_pre;
  
  int found = 0;
  for(int i = 0; i<fHitNum; i++) {
    if(fHitPool[i].Status & HITUNAV) continue;
    if(seed == i) continue;
    
    pV3.SetXYZ(fHitPool[i].X, fHitPool[i].Y, fHitPool[i].Z);
    
    //By Jixie: this block has been taken care when filling the hit pool
    //for some reason, if the hit is invalid, it could be set to (0,0,0) or (9999.,9999.;9999.)
    //It will be better not to include these hit when filling the pool
    //This hard coded part should only be used by RTPC12
    //if (pV3.Perp()<kRTPC_R_GEM1+1.0 || pV3.Perp()>kRTPC_R_Cathode-1.0) continue;  
  
    if(i>0) pV3_pre.SetXYZ(fHitPool[i-1].X, fHitPool[i-1].Y, fHitPool[i-1].Z);
    
    // removes the same hits (perhaps redundant)
    //if (!(pV3-pV3_pre).Mag() == 0)  never use ==0 to judge a floating number
    if ( (pV3-pV3_pre).Mag() <= 1.0E-5 ) {
      fHitPool[i].Status |= HISUSED;
      continue;
    } 
    
    pV3Diff = pV3 - pV3_seed;
    
    //check the distance
    double separation = pV3Diff.Mag(), acceptance=0.0;
    if( separation > Max_Link_Sep ) continue; 
    
    //for the first seed of a chain, do not require angle in range
    //because this angle has a very large range
    //Fix me:  try to find a good cut for this
    if(seed == seed_pre) {
#ifdef _ChainFinderDebug_
      if(_ChainFinderDebug_>=5) {
	printf("\t seed=%-3d seed_pre=%-3d hit=%-3d: separation=%6.2f cm, angle=%7.1f deg\n",
	  seed, seed_pre, i, separation, acceptance*rad2deg);
      }
#endif
      AddAHitToChain(fChainNum,i);
      found++;
      continue;
    }
    
    
    //check the angle between pV3_diff and pV3_diff_pre
    acceptance = pV3Diff.Angle(pV3Diff_pre);
    
    //not very sure we need this line    
    if(acceptance>90./rad2deg) acceptance = 180./rad2deg - acceptance;
    
#ifdef _ChainFinderDebug_
    if(_ChainFinderDebug_>=5) {
      printf("\t seed=%-3d seed_pre=%-3d hit=%-3d: separation=%6.2f cm, angle=%7.1f deg\n",
	seed, seed_pre, i, separation, acceptance*rad2deg);
    }
#endif

    //in order to run fast, we should separate this condition judgement
    //we should always put large probability terms in the front    
    if (separation <= Ang_Sep && acceptance < Max_Ang) {
      AddAHitToChain(fChainNum,i);
      found++;
      continue;
    }
    if (separation >  Ang_Sep && acceptance < Min_Ang) {					
      AddAHitToChain(fChainNum,i);
      found++;
      continue;
    } 
    
  }
  return found;
}


void ChainFinder::SearchChains() //HitStruct *fHitPool, int nhits)
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
      AddAHitToChain(fChainNum,seed);
      //fHitNumInAChain[fChainNum] will be added by 1 inside 
    }
        
    //SEARCH ALGORITHM----->
    //search a chain: looping over all founded seeds
    //Note that fHitNumInAChain[fChainNum] will self increasing if hits added into the chain
    for (int seed_idx=0; seed_idx<fHitNumInAChain[fChainNum]; seed_idx++) {
    
      //Fix me: currently the chain result are not sorted yet
      //need to do this before passing to Kalman Filter
      seed = fHitIDInAChain[fChainNum][seed_idx]; 
      if (seed_idx==0) seed_pre=seed;
      int found = SearchHitsForASeed(seed, seed_pre);
      
      //check any hits are found based on this seed, if no more hits found and not reach 
      //the last hit of current chain yet, release this hit back to the pool
      if(found) seed_pre = seed;
      else {
	//Remove this seed from seed_list. This seed might be added back by other seeds, 
	//which in consequence A) mess up the order b)can have multiple redundate copies 
	//Therefore we have to remove redumdant seeds at the end of the search
	//after that, do not forget to sort the chain
        if(seed_idx+1<fHitNumInAChain[fChainNum]) {
          //fHitPool[seed].Status &= ~HISUSED;
          //RemoveAHitFromChain_At(fChainNum,seed_idx);
          //seed_idx--;
#ifdef _ChainFinderDebug_
	      if(_ChainFinderDebug_>=4) {
	       // cout<<"    seed="<<seed<<" seed_pre="<<seed_pre<<", hit "
	       //     <<seed<<" is released back to the pool\n";
	      }
#endif
        }
      }
    }
    
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=3) {
    PrintAChain(fChainNum);
  }
#endif

    //now, remove redundant seeds from chain buffer
    RemoveRedundantFromChain(fChainNum);

    //now, sort the chain
    SortAChain(fChainNum);
    
    //judge if this chain is valid or not, do not store it if less than 5 hits
    if( fHitNumInAChain[fChainNum] >= Min_HITS_PER_CHAIN) {
        StoreAChain(fChainNum);
      fChainNum++;
    } else {
      //If you do not want to save this segment, clear the buffer for this chain
      //not sure if we need to release hits back to the pool!
      bool keep_all_segment = true;
      if (keep_all_segment) {
	fChainNum++;
      }
      else {
	for(int jj=0;jj<fHitNumInAChain[fChainNum];jj++) {
	  fHitIDInAChain[fChainNum][jj]=-1;
	}
	fHitNumInAChain[fChainNum]=0;
      }
    }
  
  }// for (int anchor_hit = 0 ......

}
  
void ChainFinder::PrintAChain(int chainid)
{
  int n = fHitNumInAChain[chainid];
  printf("\nChain #%2d: %4d hits: \n", chainid, n);

  printf("%-8s: ","hit-id");
  for (int jj=0; jj<n; jj++) {
    printf("%7d", fHitIDInAChain[chainid][jj]);
    if (!((jj+1)%15) && jj+1!=n) printf("\n");
  }
  printf("\n");
  printf("%-8s: ","s(cm)");
  for (int jj=0; jj<n; jj++) {
    printf("%7.2f", fHitPool[fHitIDInAChain[chainid][jj]].S);
    if (!((jj+1)%15) && jj+1!=n) printf("\n");
  }
  printf("\n");
  printf("%-8s: ","phi(deg)");
  for (int jj=0; jj<n; jj++) {
    printf("%7.1f", fHitPool[fHitIDInAChain[chainid][jj]].Phi*rad2deg);
    if (!((jj+1)%15) && jj+1!=n) printf("\n");
  }
  printf("\n\n");
}


void ChainFinder::SortAChain(int chainid)
{
  ;
}

void ChainFinder::StoreAChain(int chainid)
{
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=1) {
    printf("  Store chain #%d: %d hits \n", chainid,fHitNumInAChain[chainid]);
  }
#endif
  for (int jj=0; jj<fHitNumInAChain[chainid]; jj++){
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=1) {
      printf(" %d", fHitIDInAChain[chainid][jj]);
  }
#endif
      //Fill the pointer of found chain
      int hitid=fHitIDInAChain[chainid][jj];
      fChainBuf[chainid].Hits[jj] = &(fHitPool[hitid]);
     
      //By Jixie:  I add ChainInfo to tell ChainIndex cc and HitIndex jj
      //ChainInfo  = cccjjj,  where ccc is ChainIndex and jj is HitIndex 
      fHitPool[hitid].ChainInfo = chainid*1.0E3 + jj;
    }

  fChainBuf[chainid].HitNum = fHitNumInAChain[chainid]; //number of hits in the chain
  fChainBuf[chainid].ID = chainid;  //this chain index 
  fChainNum_Stored++;

#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=1) {
  printf("\n\n");
  }
#endif

}




