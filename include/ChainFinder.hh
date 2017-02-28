//Modify by Jixie: 
//Change this code to a class
//The whole class is in unit of cm

#ifndef _ChainFinder_H_
#define _ChainFinder_H_ 1

#define MAX_HITS_PER_EVENT  5000   
#define MAX_CHAINS_PER_EVENT 100        
#define MAX_HITS_PER_CHAIN   200   

//This is the maximum separation to be included in the chain, in cm
//static double MAX_LINK_SEP = 11 ;  


typedef struct {
  double X;
  double Y;
  double Z;
  int Status;
  //the following is added by Jixie in order to do sorting
  double S;       // distance to beam line
  double Phi;     //from -pi to pi
  int ThrownTID;  //To tell which track it was originally from
}HitStruct;


typedef struct {
  int ID;     //Chain ID
  int HitNum; //Number of hits in this chain
  HitStruct* Hits[MAX_HITS_PER_CHAIN];
}ChainStruct;


class ChainFinder {
public:
  ChainFinder();
  virtual ~ChainFinder();
  
  void Reset();
  void PrepareHitPool(double *x, double *y, double *z, int n, int append=0);
  void AddHitToChain(int hitid);
  void SearchHitsForASeed(int seed, int seed_pre=0); //, HitStruct* fHitPoolPtr, int nhits);
  void SearchChains();   //HitStruct *fHitPoolPtr, int nhits);
  void StoreAChain();
  void SetHitStatus(HitStruct *hit);
  
  void InsertAHit(int hitid, double x, double y, double z, double ThrownTID);
  void RemoveAHit(int hitid);
  void RemoveBadHits();

  void SetParameters(double space, double max_ang, double min_ang, double ang_sep);


public:
  int       fHitNum;     //Number of Hits in the pool
  HitStruct fHitPool[MAX_HITS_PER_EVENT]; //Keep all hits in one event

  ChainStruct fChainBuf[MAX_CHAINS_PER_EVENT];  //Keep all chains
  int         fHitNumInAChain[MAX_CHAINS_PER_EVENT];  //keep number of hits on each chain
  int         fChainNum;

  //int    anchor_hit, seed_hit, next_hit, seed_index;

  //only keep the hit ID such that we can find the chains
  int    fHitIDInAChain[MAX_CHAINS_PER_EVENT][MAX_HITS_PER_CHAIN];

private:

  //TVector3 fV3Diff;      //store 3 vector of (current_hit - seed) 
  //TVector3 fV3Diff_pre;  //store previous fV3_deff


  double Max_Link_Sep;
  double Max_Ang;
  double Min_Ang;
  double Ang_Sep;

};

#endif 
