//GEMC root file reader
//It contains 2 trees, Rec and Gen, using std::vector other than array
//The Gen Tree contains only single value leaves
#include <iostream>
#include <iomanip>
using namespace std;

#include "ReadGEMC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

ReadGEMC::ReadGEMC(const char *InFileName) : fChain(0), fCurrent(-1)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.

  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(InFileName);
  if (!f || !f->IsOpen()) f = TFile::Open(InFileName);
  if(!f) {
    cout<<"Error: can not open input root file '"<<InFileName<<"' ... "<<endl;
    return;  
  }
  cout<<"ReadGEMC open file "<<InFileName<<endl; 
  f->GetObject("Rec",fChain);
  f->GetObject("Gen",fChainGen);

  Init();
}

ReadGEMC::~ReadGEMC()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t ReadGEMC::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

void ReadGEMC::Init()
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  if(!fChain) return;
  
  // Set object pointer
  ChanID =0;
  X = 0;
  Y = 0;
  Z = 0;
  ADC = 0;
  TDC = 0;
  // Set branch addresses and branch pointers
  fCurrent = -1;
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("Event", &EventID, &b_EventID);
  fChain->SetBranchAddress("ChanID", &ChanID, &b_ChanID);
  fChain->SetBranchAddress("TrackID", &TrackId, &b_TrackId);
  fChain->SetBranchAddress("X_rec", &X, &b_X);
  fChain->SetBranchAddress("Y_rec", &Y, &b_Y);
  fChain->SetBranchAddress("Z_rec", &Z, &b_Z);
  fChain->SetBranchAddress("ADC", &ADC, &b_ADC);
  fChain->SetBranchAddress("Time", &TDC, &b_TDC);
  fChain->SetBranchAddress("TimeShift", &TimeShift, &b_TimeShift);
   
  fChainGen->SetBranchAddress("event_v", &event_v, &b_event_v);
  fChainGen->SetBranchAddress("TrackID", &trackid_v, &b_trackid_v);
  fChainGen->SetBranchAddress("z_v", &z_v, &b_z_v);
  fChainGen->SetBranchAddress("p_v", &p_v, &b_p_v);
  fChainGen->SetBranchAddress("pt_v", &pt_v, &b_pt_v);
  fChainGen->SetBranchAddress("th_v", &th_v, &b_theta_v);
  fChainGen->SetBranchAddress("phi_v", &phi_v, &b_phi_v);
}


Int_t ReadGEMC::LoadATrack()
{
  // get a valid track which HitNum_m>=5,  return -1 at end of file
  // otherwise return  HitNum_m
  int HitNum_m = 0;

  while (HitNum_m<5) {
   
    if(TDC) {
      ChanID->clear();
      X->clear();
      Y->clear();
      Z->clear();
      ADC->clear();
      TDC->clear();
      TrackId->clear();
      TimeShift->clear();
    }
    
    fChain->GetEntry(++fCurrent);
    fChainGen->GetEntry(fCurrent);
    if(fCurrent>=int(fChain->GetEntriesFast()))  return -1;
    if(EventID!=event_v) {
      cout<<"\n***Warning: GEMMC root tree \'Rec\' is not syncronized with Gen tree***\n";
    }
    HitNum_m = TDC->size();
  }
  
  //print out for debug
  cout<<"\nReadGEMC:: GEMC event "<<event_v<<": "<<th_v->size()<<" tracks, "
      <<TDC->size()<<" hits in total."<<endl;
  /*
  cout<<"Number of entries:\n"
      <<"ChanID "<<" \t "<<ChanID->size()<<endl
      <<"TrackId"<<" \t "<<TrackId->size()<<endl
      <<"ADC    "<<" \t "<<ADC->size()<<endl
      <<"TDC    "<<" \t "<<TDC->size()<<endl
      <<"X      "<<" \t "<<X->size()<<endl
      <<"Y      "<<" \t "<<Y->size()<<endl
      <<"Z      "<<" \t "<<Z->size()<<endl
      <<"TimeShift"<<" \t "<<TimeShift->size()<<endl
      <<"p_v      "<<" \t "<<p_v->size()<<endl
      <<"th_v      "<<" \t "<<th_v->size()<<endl;
  
  for(int i=0;i<HitNum_m;i++) {
    cout<<(int)(TrackId->at(i))<<"  ";
    if( !((i+1)%20) ) cout<<endl;
  } 
  cout<<endl;
  */ 

  return HitNum_m;
}
