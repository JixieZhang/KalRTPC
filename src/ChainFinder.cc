//Originally from Howard and Carlos
//Jixie change this code into a class
//
#include <iostream>
#include <math.h>
#include <stdio.h>
#include "TVector3.h"

#include <algorithm>
#include "ChainFinder.hh"

extern const double kRTPC_R_GEM1 = 7.0;
extern const double kRTPC_R_Cathode = 3.0;

///////////////////////////////////////////////////////////////////////
#define _ChainFinderDebug_ 10

#ifdef _ChainFinderDebug_
//#include "GlobalDebuger.hh"
#include "TBenchmark.h"
#endif

///////////////////////////////////////////////////////////////////////

/*
//this is an example how to use this class 
void Example()
{
  ChainFinder pTF;  
  pTF->SetParameters(v1,v2,v3,v4);
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
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=2) {
    printf("\t*****Hit %2d is removed from pool!*****\n",hitid);
  }
#endif
  fHitNum -= 1;
  return found;
}

//This is for changing status for hits before search
//For example, if we want to remove bad hits or redundate hit
//if we want fast speed, we should remove these hits
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
  
void ChainFinder::AddAHitToChain(int seed, int chainid, int hitid)
{
  if (fHitNumInAChain[chainid] >= MAX_HITS_PER_CHAIN) {        
    printf("Too many hits for the chain list. Skip...\n"); 
  }
    
  //THE HIT INDEX WHICH ACOMPLISH THE CONDITION, IS STORED HERE-->
  fHitIDInAChain[chainid][fHitNumInAChain[chainid]] = hitid;
  fParentSeed[fHitNumInAChain[chainid]] = seed;
  
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
  ;
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
    double separation = pV3Diff.Mag(); 
    if( separation > Max_Link_Sep ) continue; 
    
    //check the angle between pV3_diff and pV3_diff_pre
    double acceptance = pV3Diff.Angle(pV3Diff_pre);
    //not very sure we need this line    
    if(acceptance>90./rad2deg) acceptance = 180./rad2deg - acceptance;
    
#ifdef _ChainFinderDebug_
    if(_ChainFinderDebug_>=5) {
      printf("\t seed=%-3d seed_pre=%-3d hit=%-3d: separation=%6.2f cm, angle=%7.1f deg  ",
	seed, seed_pre, i, separation, acceptance*rad2deg);
    }
#endif

    //for the first seed of a chain, do not require angle in range
    //because this angle has a very large range
    //Fix me:  try to find a good cut for this
    if(seed == seed_pre) {
      AddAHitToChain(seed,fChainNum,i);
      found++;
      continue;
    }

    //in order to run fast, we should separate this condition judgement
    //we should always put large probability terms in the front    
    if (separation <= Ang_Sep && acceptance < Max_Ang) {
      AddAHitToChain(seed,fChainNum,i);
      found++;
      continue;
    }
    if (separation >  Ang_Sep && acceptance < Min_Ang) {					
      AddAHitToChain(seed,fChainNum,i);
      found++;
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
      AddAHitToChain(seed,fChainNum,seed);
      //fHitNumInAChain[fChainNum] will be added by 1 inside 
    }
        
    //SEARCH ALGORITHM----->
    //search a chain: looping over all founded seeds
    //Note that fHitNumInAChain[fChainNum] will self increasing if hits added into the chain
    for (int seed_idx=0; seed_idx<fHitNumInAChain[fChainNum]; seed_idx++) {
    
      //Fix me: currently the chain result are not sorted yet
      //need to do this before passing to Kalman Filter
      seed = fHitIDInAChain[fChainNum][seed_idx]; 
      seed_pre=fParentSeed[seed_idx];
      int found = SearchHitsForASeed(seed, seed_pre);
      
      //check any hits are found based on this seed, if no more hits found and not reach 
      //the last hit of current chain yet, release this hit back to the pool
      if(!found) {
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

    
    //judge if this chain is valid or not, do not store it if less than 5 hits
    if( fHitNumInAChain[fChainNum] >= Min_HITS_PER_CHAIN) {
      //now, sort the chain
      if(do_sort) SortAChain(fChainNum);
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
  printf("\n");

#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=9
  if(_ChainFinderDebug_>=10) {
    printf("%-8s: ","TDC");
    for (int jj=0; jj<n; jj++) {
      printf("%7d", fHitPool[fHitIDInAChain[chainid][jj]].TDC);
      if (!((jj+1)%15) && jj+1!=n) printf("\n");
    }
    printf("\n");
  }
    
  if(_ChainFinderDebug_>=11) {
    printf("%-8s:","ThrownTID");
    for (int jj=0; jj<n; jj++) {
      printf("%7d", fHitPool[fHitIDInAChain[chainid][jj]].ThrownTID);
      if (!((jj+1)%15) && jj+1!=n) printf("\n");
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

//shell_3 has bugs, it skip some hits 
void ChainFinder::ShellSort3_S(int *arr, int size)
{
  //Narrow the array by 2 everytime
  for (int gap = size / 3; gap > 0; gap /= 3) {
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
void ChainFinder::InsertSort_Phi(int *arr, int size) 
{
  int i, j, tmp;
  //need to check if this check arossing phi=180deg line
  double dphi = fHitPool[arr[0]].Phi - fHitPool[arr[size-1]].Phi;
  bool arossline = (fabs(dphi) > kPi-0.001) ? true : false;

  //make a copy of the phi into array 
  //if cross the line, need to add 2Pi to negative phi 
  double *phi = new double [size];
  for (i = 0; i < size; i++) {
    phi[i] = fHitPool[arr[i]].Phi; 
    if (arossline && phi[i] < 0.0) phi[i] += 2*kPi; 
  }

  //do sorting now 
  for (i = 1; i < size; i++) {
    j = i;
    while (j > 0 && fHitPool[arr[j]].TDC == fHitPool[arr[j-1]].TDC &&
      phi[j] < phi[j-1]) {
#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=9
	if(_ChainFinderDebug_>=10)
	  cout<<"\tInsertSort_Phi(): swap("<<j<<", "<<j-1<<") ..."<<endl;
#endif
	tmp = arr[j];
	arr[j] = arr[j-1];
	arr[j-1] = tmp;
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
void ChainFinder::SortAChain(int chainid)
{
  //now sort array buf using its S and Phi 
#if defined _ChainFinderDebug_ && _ChainFinderDebug_>=9
  //By Jixie; I want to do a bench mark test which sort algrithm works better 
  //Conclusion of sorting:
  //1) Bubble and selection sort are most slow in any case
  //2) Bubble, selection and incertion sort are stable (stable means will not change 
  //order if two values are equal to each other).  Quick and shell sort are not stable.
  //3) Incertion sort is the fatest only if the array already sorted. It will be as 
  //slow as bubble if the array is sorted in opposite way. For random array it is still 
  //faster than bubble and selection sort. 
  //4) Shell_3 is faster than shell_2 and shell_seq.
  //5) For random array, if array size>=30, quick sort is the fastest one. If array 	
  //size<30, shell_3 is faster than quick sort.
  //6) For the chain from this chainfinder, (they are almost sorted only a few location 
  //need to adjust), if array size<80, shell_3 is faster than quicksort.  
  //If size<40 incertion sort is the best. 

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
  }
  pBenchmark.Stop("stress_shell2");
  pBenchmark.Print("stress_shell2");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");
  
  /////////////////////////////////////////////////
  pBenchmark.Start("stress_shell3");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf1[i]=buf0[i];
    ShellSort3_S(buf1,n);
  }
  pBenchmark.Stop("stress_shell3");
  pBenchmark.Print("stress_shell3");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");

  /////////////////////////////////////////////////
  pBenchmark.Start("stress_shell_seq");
  for(int k=0;k<1000000;k++) {
    for(int i=0;i<n;i++) buf1[i]=buf0[i];
    ShellSort_Seq_S(buf1,n);
  }
  pBenchmark.Stop("stress_shell_seq");
  pBenchmark.Print("stress_shell_seq");
  for(int i=0;i<n;i++) printf("%3d ",buf1[i]);
  printf("\n\n");

#endif
   
  /////////////////////////////////////////////////
  //do the sort
  //get number of hits and make buf point to this array
  int nhits = fHitNumInAChain[chainid];
  int *buf = fHitIDInAChain[chainid];
  
  if(nhits<40) InsertSort_S(buf,nhits);
  //else if(nhits<80) ShellSort3_S(buf,nhits); //shell_3 has bugs
  else QuickSort_S(buf,0,nhits-1);

#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=3) {
    printf("*****  After sortting by S: *****\n");
    PrintAChain(fChainNum);
  }
#endif
  
  //sorting by phi
  InsertSort_Phi(buf,nhits);
    
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=3) {
    printf("*****  After sortting by Phi 1st: *****\n");
    PrintAChain(fChainNum);
  }
#endif

  //sorting by phi again
  InsertSort_Phi(buf,nhits);
    
#ifdef _ChainFinderDebug_
  if(_ChainFinderDebug_>=3) {
    printf("*****  After sortting by Phi 2nd: *****\n");
    PrintAChain(fChainNum);
  }
#endif

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

