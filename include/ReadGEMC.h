//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 20 15:07:06 2017 by ROOT version 5.28/00
// from TTree eventtree/Reconstructed Data to be used only for the TrackFinder
// found on file: rec_for_CF.root
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
  vector<double>  *ChanID;
  vector<double>  *X;
  vector<double>  *Y;
  vector<double>  *Z;
  vector<double>  *ADC;
  vector<double>  *TDC;
  vector<int>     *TrackId;
  vector<int>     *TimeShift;

  // List of branches  for Rec tree
  TBranch        *b_ChanID;   //!
  TBranch        *b_X;   //!
  TBranch        *b_Y;   //!
  TBranch        *b_Z;   //!
  TBranch        *b_ADC;   //!
  TBranch        *b_TDC;   //!
  TBranch        *b_TrackId;   //!
  TBranch        *b_TimeShift;   //!
  
   // Declaration of leaf types for Gen tree
   Int_t           event_v;
   Double_t        z_v;
   Double_t        p_v;
   Double_t        pt_v;
   Double_t        th_v;
   Double_t        phi_v;

   // List of branches  for Gen tree
   TBranch        *b_EventID;   //!
   TBranch        *b_z_v;   //!
   TBranch        *b_p_v;   //!
   TBranch        *b_pt_v;   //!
   TBranch        *b_theta_v;   //!
   TBranch        *b_phi_v;   //!

};

typedef  ReadGEMC GEMCReader;

#endif
