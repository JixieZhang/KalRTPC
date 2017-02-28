//////////////////////////////////////////////////////////
// By Jixie:  ntuple reader to read ep tree, which is currently
// the stage1 tree of RTPC12 G4 sim. progrom
// Please note that the unit in this reader are
// mm for length, MeV for Edep,  GeV for momentum  
//////////////////////////////////////////////////////////

#ifndef NtReader_h
#define NtReader_h

#include <iostream>
#include <iomanip>
using namespace std;

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//print some information for debugging 
//#define _NtReaderDebug_  2

#define MaxHitInATrack  200

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class NtReader {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current event index in a TChain

  // Declaration of leaf types
  Int_t           ThrownIndex;
  Int_t           Index;
  Int_t           Pid;
  Double_t        Beam;
  Double_t        Ei;
  Double_t        X0;
  Double_t        Y0;
  Double_t        Z0;
  Double_t        P0_e;
  Double_t        Theta0_e;
  Double_t        Phi0_e;
  Double_t        Xvb_e;
  Double_t        Yvb_e;
  Double_t        Zvb_e;
  Double_t        Pvb_e;
  Double_t        Thetavb_e;
  Double_t        Phivb_e;
  Double_t        XS_e;
  Double_t        Theta0_rec_e;
  Double_t        Phi0_rec_e;
  Double_t        P0_rec_e;
  Double_t        X0_rec_e;
  Double_t        Y0_rec_e;
  Double_t        Z0_rec_e;
  Double_t        P0_p;
  Double_t        Theta0_p;
  Double_t        Phi0_p;
  Double_t        Xvb_p;
  Double_t        Yvb_p;
  Double_t        Zvb_p;
  Double_t        Pvb_p;
  Double_t        Thetavb_p;
  Double_t        Phivb_p;
  Double_t        Theta0_rec_p;
  Double_t        Phi0_rec_p;
  Double_t        P0_rec_p;
  Double_t        X0_rec_p;
  Double_t        Y0_rec_p;
  Double_t        Z0_rec_p;
  Double_t        Smin;
  Double_t        Smax;
  Double_t        MeanBz;
  Double_t        dEdX;
  Int_t           HitNum;
  Double_t        StepX[MaxHitInATrack];   //[HitNum]
  Double_t        StepY[MaxHitInATrack];   //[HitNum]
  Double_t        StepZ[MaxHitInATrack];   //[HitNum]
  Double_t        StepS[MaxHitInATrack];   //[HitNum]
  Double_t        StepPhi[MaxHitInATrack];   //[HitNum]
  Double_t        StepdE[MaxHitInATrack];   //[HitNum]
  Double_t        StepL[MaxHitInATrack];   //[HitNum]
  Int_t           StepID[MaxHitInATrack];   //[HitNum]
  Int_t           StepTDC[MaxHitInATrack];   //[HitNum]
  Int_t           StepADC[MaxHitInATrack];   //[HitNum]
  Double_t        StepX_rec[MaxHitInATrack];   //[HitNum]
  Double_t        StepY_rec[MaxHitInATrack];   //[HitNum]
  Double_t        StepZ_rec[MaxHitInATrack];   //[HitNum]
  Double_t        StepS_rec[MaxHitInATrack];   //[HitNum]
  Double_t        StepPhi_rec[MaxHitInATrack];   //[HitNum]
  Double_t        R_sim;
  Double_t        A_sim;
  Double_t        B_sim;
  Double_t        Theta_sim;
  Double_t        Phi_sim;
  Double_t        Z_sim;
  Double_t        DCA_sim;
  Int_t           HitNum_m;
  Int_t           StepID_m[MaxHitInATrack];   //[HitNum_m]
  Int_t           StepTDC_m[MaxHitInATrack];   //[HitNum_m]
  Int_t           StepADC_m[MaxHitInATrack];   //[HitNum_m]
  Double_t        StepX_rec_m[MaxHitInATrack];   //[HitNum_m]
  Double_t        StepY_rec_m[MaxHitInATrack];   //[HitNum_m]
  Double_t        StepZ_rec_m[MaxHitInATrack];   //[HitNum_m]
  Double_t        StepS_rec_m[MaxHitInATrack];   //[HitNum_m]
  Double_t        StepPhi_rec_m[MaxHitInATrack];   //[HitNum_m]
  Double_t        R_rec;
  Double_t        A_rec;
  Double_t        B_rec;
  Double_t        Theta_rec;
  Double_t        Phi_rec;
  Double_t        Z_rec;
  Double_t        DCA_rec;

  // List of branches
  TBranch        *b_ThrownIndex;   //!
  TBranch        *b_Index;   //!
  TBranch        *b_Pid;   //!
  TBranch        *b_Beam;   //!
  TBranch        *b_Ei;   //!
  TBranch        *b_X0;   //!
  TBranch        *b_Y0;   //!
  TBranch        *b_Z0;   //!
  TBranch        *b_P0_e;   //!
  TBranch        *b_Theta0_e;   //!
  TBranch        *b_Phi0_e;   //!
  TBranch        *b_Xvb_e;   //!
  TBranch        *b_Yvb_e;   //!
  TBranch        *b_Zvb_e;   //!
  TBranch        *b_Pvb_e;   //!
  TBranch        *b_Thetavb_e;   //!
  TBranch        *b_Phivb_e;   //!
  TBranch        *b_XS_e;   //!
  TBranch        *b_Theta0_rec_e;   //!
  TBranch        *b_Phi0_rec_e;   //!
  TBranch        *b_P0_rec_e;   //!
  TBranch        *b_X0_rec_e;   //!
  TBranch        *b_Y0_rec_e;   //!
  TBranch        *b_Z0_rec_e;   //!
  TBranch        *b_P0_p;   //!
  TBranch        *b_Theta0_p;   //!
  TBranch        *b_Phi0_p;   //!
  TBranch        *b_Xvb_p;   //!
  TBranch        *b_Yvb_p;   //!
  TBranch        *b_Zvb_p;   //!
  TBranch        *b_Pvb_p;   //!
  TBranch        *b_Thetavb_p;   //!
  TBranch        *b_Phivb_p;   //!
  TBranch        *b_Theta0_rec_p;   //!
  TBranch        *b_Phi0_rec_p;   //!
  TBranch        *b_P0_rec_p;   //!
  TBranch        *b_X0_rec_p;   //!
  TBranch        *b_Y0_rec_p;   //!
  TBranch        *b_Z0_rec_p;   //!
  TBranch        *b_Smin;   //!
  TBranch        *b_Smax;   //!
  TBranch        *b_MeanBz;   //!
  TBranch        *b_dEdX;   //!
  TBranch        *b_HitNum;   //!
  TBranch        *b_StepX;   //!
  TBranch        *b_StepY;   //!
  TBranch        *b_StepZ;   //!
  TBranch        *b_StepS;   //!
  TBranch        *b_StepPhi;   //!
  TBranch        *b_StepdE;   //!
  TBranch        *b_StepL;   //!
  TBranch        *b_StepID;   //!
  TBranch        *b_StepTDC;   //!
  TBranch        *b_StepADC;   //!
  TBranch        *b_StepX_rec;   //!
  TBranch        *b_StepY_rec;   //!
  TBranch        *b_StepZ_rec;   //!
  TBranch        *b_StepS_rec;   //!
  TBranch        *b_StepPhi_rec;   //!
  TBranch        *b_R_sim;   //!
  TBranch        *b_A_sim;   //!
  TBranch        *b_B_sim;   //!
  TBranch        *b_Theta_sim;   //!
  TBranch        *b_Phi_sim;   //!
  TBranch        *b_Z_sim;   //!
  TBranch        *b_DCA_sim;   //!
  TBranch        *b_HitNum_m;   //!
  TBranch        *b_StepID_m;   //!
  TBranch        *b_StepTDC_m;   //!
  TBranch        *b_StepADC_m;   //!
  TBranch        *b_StepX_rec_m;   //!
  TBranch        *b_StepY_rec_m;   //!
  TBranch        *b_StepZ_rec_m;   //!
  TBranch        *b_StepS_rec_m;   //!
  TBranch        *b_StepPhi_rec_m;   //!
  TBranch        *b_R_rec;   //!
  TBranch        *b_A_rec;   //!
  TBranch        *b_B_rec;   //!
  TBranch        *b_Theta_rec;   //!
  TBranch        *b_Phi_rec;   //!
  TBranch        *b_Z_rec;   //!
  TBranch        *b_DCA_rec;   //!

  NtReader(const char *InFileName="infile.root");
  virtual ~NtReader();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual void     Init();
  virtual void     Loop();
  virtual void     Show(Long64_t entry = -1);

  //Load a valid track (with NHitNum>5 and Smax>50mm), returm -1 if end of tree 
  virtual Int_t    LoadATrack();
};

#endif

#ifdef NtReader_cxx
NtReader::NtReader(const char *InFileName) : fChain(0), fCurrent(-1)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.

  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(InFileName);
  if (!f || !f->IsOpen()) f = TFile::Open(InFileName);
  if(!f) {
    cout<<"Error: can not open input root file '"<<InFileName<<"' ... "<<endl;
    return;  
  }
  f->GetObject("ep",fChain);

  Init();
}

NtReader::~NtReader()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t NtReader::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

void NtReader::Init()
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  if(!fChain) return;

  // Set branch addresses and branch pointers
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("ThrownIndex", &ThrownIndex, &b_ThrownIndex);
  fChain->SetBranchAddress("Index", &Index, &b_Index);
  fChain->SetBranchAddress("Pid", &Pid, &b_Pid);
  fChain->SetBranchAddress("Beam", &Beam, &b_Beam);
  fChain->SetBranchAddress("Ei", &Ei, &b_Ei);
  fChain->SetBranchAddress("X0", &X0, &b_X0);
  fChain->SetBranchAddress("Y0", &Y0, &b_Y0);
  fChain->SetBranchAddress("Z0", &Z0, &b_Z0);
  fChain->SetBranchAddress("P0_e", &P0_e, &b_P0_e);
  fChain->SetBranchAddress("Theta0_e", &Theta0_e, &b_Theta0_e);
  fChain->SetBranchAddress("Phi0_e", &Phi0_e, &b_Phi0_e);
  fChain->SetBranchAddress("Xvb_e", &Xvb_e, &b_Xvb_e);
  fChain->SetBranchAddress("Yvb_e", &Yvb_e, &b_Yvb_e);
  fChain->SetBranchAddress("Zvb_e", &Zvb_e, &b_Zvb_e);
  fChain->SetBranchAddress("Pvb_e", &Pvb_e, &b_Pvb_e);
  fChain->SetBranchAddress("Thetavb_e", &Thetavb_e, &b_Thetavb_e);
  fChain->SetBranchAddress("Phivb_e", &Phivb_e, &b_Phivb_e);
  fChain->SetBranchAddress("XS_e", &XS_e, &b_XS_e);
  fChain->SetBranchAddress("Theta0_rec_e", &Theta0_rec_e, &b_Theta0_rec_e);
  fChain->SetBranchAddress("Phi0_rec_e", &Phi0_rec_e, &b_Phi0_rec_e);
  fChain->SetBranchAddress("P0_rec_e", &P0_rec_e, &b_P0_rec_e);
  fChain->SetBranchAddress("X0_rec_e", &X0_rec_e, &b_X0_rec_e);
  fChain->SetBranchAddress("Y0_rec_e", &Y0_rec_e, &b_Y0_rec_e);
  fChain->SetBranchAddress("Z0_rec_e", &Z0_rec_e, &b_Z0_rec_e);
  fChain->SetBranchAddress("P0_p", &P0_p, &b_P0_p);
  fChain->SetBranchAddress("Theta0_p", &Theta0_p, &b_Theta0_p);
  fChain->SetBranchAddress("Phi0_p", &Phi0_p, &b_Phi0_p);
  fChain->SetBranchAddress("Xvb_p", &Xvb_p, &b_Xvb_p);
  fChain->SetBranchAddress("Yvb_p", &Yvb_p, &b_Yvb_p);
  fChain->SetBranchAddress("Zvb_p", &Zvb_p, &b_Zvb_p);
  fChain->SetBranchAddress("Pvb_p", &Pvb_p, &b_Pvb_p);
  fChain->SetBranchAddress("Thetavb_p", &Thetavb_p, &b_Thetavb_p);
  fChain->SetBranchAddress("Phivb_p", &Phivb_p, &b_Phivb_p);
  fChain->SetBranchAddress("Theta0_rec_p", &Theta0_rec_p, &b_Theta0_rec_p);
  fChain->SetBranchAddress("Phi0_rec_p", &Phi0_rec_p, &b_Phi0_rec_p);
  fChain->SetBranchAddress("P0_rec_p", &P0_rec_p, &b_P0_rec_p);
  fChain->SetBranchAddress("X0_rec_p", &X0_rec_p, &b_X0_rec_p);
  fChain->SetBranchAddress("Y0_rec_p", &Y0_rec_p, &b_Y0_rec_p);
  fChain->SetBranchAddress("Z0_rec_p", &Z0_rec_p, &b_Z0_rec_p);
  fChain->SetBranchAddress("Smin", &Smin, &b_Smin);
  fChain->SetBranchAddress("Smax", &Smax, &b_Smax);
  fChain->SetBranchAddress("MeanBz", &MeanBz, &b_MeanBz);
  fChain->SetBranchAddress("dEdX", &dEdX, &b_dEdX);
  fChain->SetBranchAddress("HitNum", &HitNum, &b_HitNum);
  fChain->SetBranchAddress("StepX", StepX, &b_StepX);
  fChain->SetBranchAddress("StepY", StepY, &b_StepY);
  fChain->SetBranchAddress("StepZ", StepZ, &b_StepZ);
  fChain->SetBranchAddress("StepS", StepS, &b_StepS);
  fChain->SetBranchAddress("StepPhi", StepPhi, &b_StepPhi);
  fChain->SetBranchAddress("StepdE", StepdE, &b_StepdE);
  fChain->SetBranchAddress("StepL", StepL, &b_StepL);
  fChain->SetBranchAddress("StepID", StepID, &b_StepID);
  fChain->SetBranchAddress("StepTDC", StepTDC, &b_StepTDC);
  fChain->SetBranchAddress("StepADC", StepADC, &b_StepADC);
  fChain->SetBranchAddress("StepX_rec", StepX_rec, &b_StepX_rec);
  fChain->SetBranchAddress("StepY_rec", StepY_rec, &b_StepY_rec);
  fChain->SetBranchAddress("StepZ_rec", StepZ_rec, &b_StepZ_rec);
  fChain->SetBranchAddress("StepS_rec", StepS_rec, &b_StepS_rec);
  fChain->SetBranchAddress("StepPhi_rec", StepPhi_rec, &b_StepPhi_rec);
  fChain->SetBranchAddress("R_sim", &R_sim, &b_R_sim);
  fChain->SetBranchAddress("A_sim", &A_sim, &b_A_sim);
  fChain->SetBranchAddress("B_sim", &B_sim, &b_B_sim);
  fChain->SetBranchAddress("Theta_sim", &Theta_sim, &b_Theta_sim);
  fChain->SetBranchAddress("Phi_sim", &Phi_sim, &b_Phi_sim);
  fChain->SetBranchAddress("Z_sim", &Z_sim, &b_Z_sim);
  fChain->SetBranchAddress("DCA_sim", &DCA_sim, &b_DCA_sim);
  fChain->SetBranchAddress("HitNum_m", &HitNum_m, &b_HitNum_m);
  fChain->SetBranchAddress("StepID_m", StepID_m, &b_StepID_m);
  fChain->SetBranchAddress("StepTDC_m", StepTDC_m, &b_StepTDC_m);
  fChain->SetBranchAddress("StepADC_m", StepADC_m, &b_StepADC_m);
  fChain->SetBranchAddress("StepX_rec_m", StepX_rec_m, &b_StepX_rec_m);
  fChain->SetBranchAddress("StepY_rec_m", StepY_rec_m, &b_StepY_rec_m);
  fChain->SetBranchAddress("StepZ_rec_m", StepZ_rec_m, &b_StepZ_rec_m);
  fChain->SetBranchAddress("StepS_rec_m", StepS_rec_m, &b_StepS_rec_m);
  fChain->SetBranchAddress("StepPhi_rec_m", StepPhi_rec_m, &b_StepPhi_rec_m);
  fChain->SetBranchAddress("R_rec", &R_rec, &b_R_rec);
  fChain->SetBranchAddress("A_rec", &A_rec, &b_A_rec);
  fChain->SetBranchAddress("B_rec", &B_rec, &b_B_rec);
  fChain->SetBranchAddress("Theta_rec", &Theta_rec, &b_Theta_rec);
  fChain->SetBranchAddress("Phi_rec", &Phi_rec, &b_Phi_rec);
  fChain->SetBranchAddress("Z_rec", &Z_rec, &b_Z_rec);
  fChain->SetBranchAddress("DCA_rec", &DCA_rec, &b_DCA_rec);   
}

Int_t NtReader::LoadATrack()
{
  // get a valid track which HitNum_m>5,  return -1 at end of file
  // otherwise return  HitNum_m
  HitNum_m = 0;

  while (HitNum_m<5 || Smax<50.) {
    fChain->GetEntry(++fCurrent);
#ifdef _NtReaderDebug_
    if(_NtReaderDebug_>=2) {
      cout<<"Ntuple Event "<<setw(5)<<Index<<":  HitNum_m="<<setw(2)<<HitNum_m
	<<",  Smax="<<setw(8)<<Smax<<",  Smin="<<setw(8)<<Smin<<endl;
      if(HitNum_m>5) {
	cout<<"\t P0="<<setw(8)<<P0_p<<",  Theta0="<<setw(8)<<Theta0_p*57.3
	  <<",  Phi0="<<setw(8)<<Phi0_p*57.3<<",  Z0="<<setw(8)<<Z0<<endl;
      }
    }       
#endif
    if(fCurrent>=int(fChain->GetEntries()))  return -1;
  }
  return HitNum_m;
}

void NtReader::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t NtReader::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  if(HitNum < 5) return 0;
  return 1;
}
#endif // #ifdef NtReader_cxx
