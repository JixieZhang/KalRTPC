//////////////////////////////////////////////////////////
// This class is used to read root tree 'rec_for_CF.root' by bonus_gemc.
// rec_for_CF.root contains two tree, 'Rec' and 'Gen', all the tree leaves 
// are vectors.
//////////////////////////////////////////////////////////

#ifndef ReadGEMC_h
#define ReadGEMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class ReadGEMC {
 public :
  ReadGEMC(const char *InFileName="infile.root");
  virtual ~ReadGEMC();
  virtual void     Init();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Int_t    LoadATrack();
  
 public :
  TTree          *fChain;    //!pointer to the analyzed TTree or TChain
  TTree          *fChainGen; //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent;  //!current event index in a TChain

  // Declaration of leaf types for Rec tree
  Int_t           EventID;
  vector<int>     *ChanID;
  vector<int>     *TrackId;
  vector<double>  *X;
  vector<double>  *Y;
  vector<double>  *Z;
  vector<double>  *ADC;
  vector<double>  *TDC;
  vector<double>  *TimeShift;

  // List of branches  for Rec tree
  TBranch        *b_EventID;   //!
  TBranch        *b_ChanID;   //!
  TBranch        *b_TrackId;   //!
  TBranch        *b_X;   //!
  TBranch        *b_Y;   //!
  TBranch        *b_Z;   //!
  TBranch        *b_ADC;   //!
  TBranch        *b_TDC;   //!
  TBranch        *b_TimeShift;   //!
  
  // Declaration of leaf types for Gen tree
  Int_t           event_v;
  vector<int>     *trackid_v;
  vector<double>  *z_v;
  vector<double>  *p_v;
  vector<double>  *pt_v;
  vector<double>  *th_v;
  vector<double>  *phi_v;

  // List of branches  for Gen tree
  TBranch        *b_event_v;   //!
  TBranch        *b_trackid_v;   //!
  TBranch        *b_z_v;   //!
  TBranch        *b_p_v;   //!
  TBranch        *b_pt_v;   //!
  TBranch        *b_theta_v;   //!
  TBranch        *b_phi_v;   //!

};

typedef  ReadGEMC GEMCReader;

#endif
