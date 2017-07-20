//This file is used to build a frame of input|output to EXKarRTPC,
//so we can separate EXKalRTPC as a clean module without reading|writing root tree
// Here is an example to use this class,
//
//  EXKalManager aKalFilter;
//  aKalFilter.SetCovMElement(error);
//  aKalFilter.SetG4InputFile(infile);
//  aKalFilter.Run(job,nevents,pt_min,pt_max,costh_min,costh_max);
//
#include "EXKalManager.h"
#include "GlobalDebuger.hh"

#include "TCanvas.h"
#include "TView.h"
#include "TView3D.h"
#include "TRotMatrix.h"
#include "TTUBE.h"
#include "TNode.h"
#include "TString.h"

#include <iomanip>
#include <iostream>
using namespace std;

ClassImp(EXKalManager)

///////////////////////////////////////////////////////////////////////
#define _EXKalManDebug_ 3

///////////////////////////////////////////////////////////////////////


TApplication* EXKalManager::fApp = 0; 

EXKalManager::EXKalManager()
{
#ifdef _EXKalManDebug_
  Global_Debug_Level=_EXKalManDebug_;
#endif

  // ===================================================================
  //  Prepare a KalRTPC
  // ===================================================================
  fChainFinder = new ChainFinder();
  fKalRTPC = new EXKalRTPC();

  // ===================================================================
  //  Set the Pointers
  // ===================================================================
  fKalHits  = fKalRTPC->fKalHits;     // hit buffer
  fEventGen = fKalRTPC->fEventGen;    // Jixie's track generator   

  _eventtype_ = 0;
  fNtReader = 0;     // will be instanced in BeginOfRun() if (_eventtype_ == 10)
  fGEMCReader = 0;   // will be instanced in BeginOfRun() if (_eventtype_ == 11)
  _index_=0;
  eventid=trackid=ntrack=0;
  
  if(fApp) delete fApp;
  fApp = new TApplication("EXKalRTPC", 0, 0, 0, 0);
}

EXKalManager::~EXKalManager()
{
  delete fKalRTPC;
  if(fNtReader) delete fNtReader;
  if(fApp) delete fApp;
}


void EXKalManager::Tree_Init()
{
  fFile = new TFile("h.root","RECREATE","Chain Finder and Kalman Filter for RTPC12");
  fTree = new TTree("t", "Chain Finder and Kalman Filter for RTPC12");

  TTree *t=fTree;
  t->Branch("index",&_index_,"index/I");
  t->Branch("type",&_eventtype_,"type/I");

  //to identify multiple track in an event
  //fill event track per tree entry such that the structure is simple
  t->Branch("eventid",&eventid,"eventid/I");
  t->Branch("trackid",&trackid,"trackid/I");
  t->Branch("ntrack",&ntrack,"ntrack/I");

  // For Chain Finder
  t->Branch("CF_ntrack_read",&CF_ntrack_read,"CF_ntrack_read/I");
  t->Branch("CF_ntrack_good",&CF_ntrack_good,"CF_ntrack_good/I");
  t->Branch("CF_HitNum",&CF_HitNum,"CF_HitNum/I");
  t->Branch("CF_ChainNum",&CF_ChainNum,"CF_ChainNum/I");
  t->Branch("CF_ID",CF_ID,"CF_ID[CF_HitNum]/I");
  t->Branch("CF_TDC",CF_TDC,"CF_TDC[CF_HitNum]/I");
  t->Branch("CF_ADC",CF_ADC,"CF_ADC[CF_HitNum]/I");
  t->Branch("CF_X",CF_X,"CF_X[CF_HitNum]/D");
  t->Branch("CF_Y",CF_Y,"CF_Y[CF_HitNum]/D");
  t->Branch("CF_Z",CF_Z,"CF_Z[CF_HitNum]/D");
  t->Branch("CF_Status",CF_Status,"CF_Status[CF_HitNum]/I");
  t->Branch("CF_S",CF_S,"CF_S[CF_HitNum]/D");
  t->Branch("CF_Phi",CF_Phi,"CF_Phi[CF_HitNum]/D");
  t->Branch("CF_ThrownTID",CF_ThrownTID,"CF_ThrownTID[CF_HitNum]/I");
  t->Branch("CF_ChainInfo",CF_ChainInfo,"CF_ChainInfo[CF_HitNum]/I");

  t->Branch("CF_X0",CF_X0,"CF_X0[CF_ntrack_read]/D");
  t->Branch("CF_Y0",CF_Y0,"CF_Y0[CF_ntrack_read]/D");
  t->Branch("CF_Z0",CF_Z0,"CF_Z0[CF_ntrack_read]/D");
  t->Branch("CF_Theta0",CF_Theta0,"CF_Theta0[CF_ntrack_read]/D");
  t->Branch("CF_Phi0",CF_Phi0,"CF_Phi0[CF_ntrack_read]/D");
  t->Branch("CF_P0",CF_P0,"CF_P0[CF_ntrack_read]/D");

  t->Branch("CF_ThrownTID_like",&CF_ThrownTID_like,"CF_ThrownTID_like/D");

  // thrown parameters
  t->Branch("p0",&p0,"p0/D");
  t->Branch("pt0",&pt0,"pt0/D");
  t->Branch("pz0",&pz0,"pz0/D");
  t->Branch("th0",&th0,"th0/D");
  t->Branch("ph0",&ph0,"ph0/D");
  t->Branch("x0",&_x0_,"x0/D");
  t->Branch("y0",&_y0_,"y0/D");
  t->Branch("z0",&_z0_,"z0/D");

  //KF result
  t->Branch("p_rec",&p_rec,"p_rec/D");
  t->Branch("pt_rec",&pt_rec,"pt_rec/D");
  t->Branch("pz_rec",&pz_rec,"pz_rec/D");
  t->Branch("th_rec",&th_rec,"th_rec/D");
  t->Branch("ph_rec",&ph_rec,"ph_rec/D");
  t->Branch("x_rec",&x_rec,"x_rec/D");
  t->Branch("y_rec",&y_rec,"y_rec/D");
  t->Branch("z_rec",&z_rec,"z_rec/D");
  t->Branch("r_rec",&r_rec,"r_rec/D");
  t->Branch("a_rec",&a_rec,"a_rec/D");
  t->Branch("b_rec",&b_rec,"b_rec/D");
  t->Branch("ndf",&ndf,"ndf/I");
  t->Branch("chi2",&chi2,"chi2/D");
  t->Branch("cl",&cl,"cl/D");

  //original hits info, from 
  t->Branch("npt0",&npt0,"npt0/I");
  t->Branch("step_x",step_x,"step_x[npt0]/D");
  t->Branch("step_y",step_y,"step_y[npt0]/D");
  t->Branch("step_z",step_z,"step_z[npt0]/D");
  t->Branch("step_phi",step_phi,"step_phi[npt0]/D");
  t->Branch("step_s",step_s,"step_s[npt0]/D");
  t->Branch("step_s_min",&step_s_min,"step_s_min/D");
  t->Branch("step_s_max",&step_s_max,"step_s_max/D");

  //smeared or reconstructed hits info
  t->Branch("npt",&npt,"npt/I");
  t->Branch("step_status",step_status,"step_status[npt]/I");

  t->Branch("step_x_rec",step_x_rec,"step_x_rec[npt]/D");
  t->Branch("step_y_rec",step_y_rec,"step_y_rec[npt]/D");
  t->Branch("step_z_rec",step_z_rec,"step_z_rec[npt]/D");
  t->Branch("step_phi_rec",step_phi_rec,"step_phi_rec[npt]/D");
  t->Branch("step_s_rec",step_s_rec,"step_s_rec[npt]/D");


  t->Branch("step_x_exp",step_x_exp,"step_x_exp[npt]/D");
  t->Branch("step_y_exp",step_y_exp,"step_y_exp[npt]/D");
  t->Branch("step_z_exp",step_z_exp,"step_z_exp[npt]/D");
  t->Branch("step_phi_exp",step_phi_exp,"step_phi_exp[npt]/D");
  t->Branch("step_s_exp",step_s_exp,"step_s_exp[npt]/D");

  t->Branch("step_x_fil",step_x_fil,"step_x_fil[npt]/D");
  t->Branch("step_y_fil",step_y_fil,"step_y_fil[npt]/D");
  t->Branch("step_z_fil",step_z_fil,"step_z_fil[npt]/D");
  t->Branch("step_phi_fil",step_phi_fil,"step_phi_fil[npt]/D");
  t->Branch("step_s_fil",step_s_fil,"step_s_fil[npt]/D");

  t->Branch("p_3pt",&p_3pt,"p_3pt/D");
  t->Branch("pt_3pt",&pt_3pt,"pt_3pt/D");
  t->Branch("th_3pt",&th_3pt,"th_3pt/D");
  t->Branch("r_3pt",&r_3pt,"r_3pt/D");
  t->Branch("a_3pt",&a_3pt,"a_3pt/D");
  t->Branch("b_3pt",&b_3pt,"b_3pt/D");

  t->Branch("p_hel",&p_hel,"p_hel/D");
  t->Branch("pt_hel",&pt_hel,"pt_hel/D");
  t->Branch("pz_hel",&pz_hel,"pz_hel/D");
  t->Branch("th_hel",&th_hel,"th_hel/D");
  t->Branch("ph_hel",&ph_hel,"ph_hel/D");
  t->Branch("x_hel",&x_hel,"x_hel/D");
  t->Branch("y_hel",&y_hel,"y_hel/D");
  t->Branch("z_hel",&z_hel,"z_hel/D");
  t->Branch("r_hel",&r_hel,"r_hel/D");
  t->Branch("a_hel",&a_hel,"a_hel/D");
  t->Branch("b_hel",&b_hel,"b_hel/D");
  t->Branch("dca_hel",&dca_hel,"dca_hel/D");
  t->Branch("chi2_hel",&chi2_hel,"chi2_hel/D");

  t->Branch("r_hel_raw",&r_hel_raw,"r_hel_raw/D");
  t->Branch("th_hel_raw",&th_hel_raw,"th_hel_raw/D");
  t->Branch("ph_hel_raw",&ph_hel_raw,"ph_hel_raw/D");
  t->Branch("a_hel_raw",&a_hel_raw,"a_hel_raw/D");
  t->Branch("b_hel_raw",&b_hel_raw,"b_hel_raw/D");
  t->Branch("z_hel_raw",&z_hel_raw,"z_hel_raw/D");

  t->Branch("rho_1st",&rho_1st,"rho_1st/D");
  t->Branch("tnl_1st",&tnl_1st,"tnl_1st/D");
  t->Branch("phi0_1st",&phi0_1st,"phi0_1st/D");
  t->Branch("rho_last",&rho_last,"rho_last/D");
  t->Branch("tnl_last",&tnl_last,"tnl_last/D");
  t->Branch("phi0_last",&phi0_last,"phi0_last/D");
  t->Branch("rho_kal_ini",&rho_kal_ini,"rho_kal_ini/D");
  t->Branch("tnl_kal_ini",&tnl_kal_ini,"tnl_kal_ini/D");
  t->Branch("phi0_kal_ini",&phi0_kal_ini,"phi0_kal_ini/D");

}

void EXKalManager::Tree_Reset_CF()
{
  //ChainFinder tree buffer, number of found tracks store 
  //at 'ntrack' and also CF_ChainNum
  for(int i=0;i<CF_HitNum;i++) {
    CF_ID[i]=CF_TDC[i]=CF_ADC[i]=-1; 
    CF_X[i]=CF_Y[i]=CF_Z[i]=CF_S[i]=CF_Phi[i]=0.0; 
    CF_Status[i]=CF_ThrownTID[i]=CF_ChainInfo[i]=-1; 
  }

  for(int i=0;i<CF_ntrack_read;i++) {
    CF_X0[i]=CF_Y0[i]=CF_Z0[i]=CF_Theta0[i]=CF_Phi0[i]=CF_P0[i]=0.0;
  }

  CF_ntrack_read=CF_ntrack_good=0;  
  CF_HitNum=CF_ChainNum=0;
  CF_ThrownTID_like=0.0;
}

void EXKalManager::Tree_Reset()
{
  //The following are KF buffer
  p0=pt0=pz0=th0=ph0=_x0_=_y0_=_z0_=0.0;

  p_rec=pt_rec=pz_rec=th_rec=ph_rec=x_rec=y_rec=z_rec=0.0;
  r_rec=a_rec=b_rec=0.0;
  ndf=0;
  chi2=cl=0.0;

  for(int i=0;i<npt0;i++) {
    step_x[i]=step_y[i]=step_z[i]=step_phi[i]=step_s[i]=0.0;
  }
  step_s_min=99999.0;
  step_s_max=-99999.0;

  for(int i=0;i<npt;i++) {
    step_status[i]=0;
    step_x_rec[i]=step_y_rec[i]=step_z_rec[i]=step_phi_rec[i]=step_s_rec[i]=0.0;
    step_x_exp[i]=step_y_exp[i]=step_z_exp[i]=step_phi_exp[i]=step_s_exp[i]=0.0;
    step_x_fil[i]=step_y_fil[i]=step_z_fil[i]=step_phi_fil[i]=step_s_fil[i]=0.0;
  }

  p_3pt=pt_3pt=th_3pt=9999.0;
  r_3pt=a_3pt=b_3pt=0.0; 

  p_hel=pt_hel=pz_hel=th_hel=ph_hel=x_hel=y_hel=z_hel=9999.0;
  r_hel=a_hel=b_hel=0.0; 
  dca_hel=chi2_hel=0.0;

  r_hel_raw=th_hel_raw=ph_hel_raw=a_hel_raw=b_hel_raw=z_hel_raw=9999.0;

  rho_1st=tnl_1st=phi0_1st=0;
  rho_last=tnl_last=phi0_last=0;
  rho_kal_ini=tnl_kal_ini=phi0_kal_ini=0;
}

void EXKalManager::BeginOfRun(int eventtype)
{
  // ===================================================================
  //  Prepare a Root tree input
  // ===================================================================
  if(eventtype==10 || eventtype==1) {
    fNtReader = new NtReader(fG4Inputfile);
    if(! fNtReader->fChain ) {
      cout<<"Error: can not open input root file, I quit ..."<<endl;
      exit(-2);  
    }
  }
  else if(eventtype==11) {
    fGEMCReader = new GEMCReader(fG4Inputfile);
    if(! fGEMCReader->fChain ) {
      cout<<"Error: can not open input root file, I quit ..."<<endl;
      exit(-2);  
    }
  }

  // ===================================================================
  //  Prepare a Root tree output
  // ===================================================================
  Tree_Init();
}

void EXKalManager::EndOfRun()
{
  fFile->Write("",TObject::kOverwrite);
}

//read ntracks from G4MC_RTPC12 output tree and fill it into ChainFinder's hit pool.
//please note that the G4 root tree use unit of mm 
int EXKalManager::FillChainFinderHitPoolGEMC(int ntracks)
{
  int nhits = 0;

  int id[MAX_HITS_PER_EVENT],tdc[MAX_HITS_PER_EVENT],adc[MAX_HITS_PER_EVENT];
  int throwntid[MAX_HITS_PER_EVENT];
  double xx[MAX_HITS_PER_EVENT],yy[MAX_HITS_PER_EVENT],zz[MAX_HITS_PER_EVENT];

  CF_ntrack_read=0;
  for(int t=0;t<ntracks;t++) {
    int n = fGEMCReader->LoadATrack();
    if( n == -1 ) {
      cout<<"Reach the end of input root file \n";
      return -1;
    }

    nhits=0;
    int HitNum_m = fGEMCReader->TDC->size();
    for(int i=0;i<HitNum_m;i++) {
      if(fGEMCReader->ChanID->at(i)>0) {
        //if negative means this hit is invalid

        id[nhits] = (int)(fGEMCReader->ChanID->at(i));
        tdc[nhits]= (int)(fGEMCReader->TDC->at(i));
        adc[nhits]= (int)(fGEMCReader->ADC->at(i));
        xx[nhits] = fGEMCReader->X->at(i)/10.;
        yy[nhits] = fGEMCReader->Y->at(i)/10.;
        zz[nhits] = fGEMCReader->Z->at(i)/10.;
        throwntid[nhits] = fGEMCReader->TrackId->at(i)*1.0E4+i;
        nhits++;
        if(nhits >= MAX_HITS_PER_EVENT) break;
      }
    }
    if(nhits==0) continue;

    CF_ntrack_read++;
    if(nhits>=MIN_HITS_PER_TRACK) CF_ntrack_good++;

    //append this track into the hit pool
    int append=1;
    fChainFinder->PrepareHitPool(id, tdc, adc, xx, yy, zz, nhits, throwntid, append);


#ifdef _EXKalManDebug_
    //just for debug
    if(_EXKalManDebug_>=3) {
      cout<<"\nGEMC Ntuple Event "<<setw(5)<<fGEMCReader->event_v<<":  HitNum_m="<<setw(2)<<fGEMCReader->TDC->size()<<endl
        <<"  P0="<<fGEMCReader->p_v<<",  Pt="<<fGEMCReader->pt_v<<", Theta0="
        <<fGEMCReader->th_v*57.3<<", Phi0="<<fGEMCReader->phi_v*57.3<<"  Z0="<<fGEMCReader->z_v/10.<<"cm"<<endl;
    }
    if(_EXKalManDebug_>=4) {
      for(int i=0;i<(int)fGEMCReader->TDC->size();i++) {
        cout<<"Hit "<<setw(2)<<i<<"("<<setw(8)<<fGEMCReader->X->at(i)<<", "
          <<setw(8)<<fGEMCReader->Y->at(i)<<", "
          <<setw(8)<<fGEMCReader->Z->at(i)<<") ==>  S="<<setw(8)<<sqrt(pow(fGEMCReader->X->at(i),2.0)+pow(fGEMCReader->Y->at(i),2.0))
          <<" mm  Phi="<<setw(8)<<atan2(fGEMCReader->Y->at(i),fGEMCReader->X->at(i))*57.3<<" deg"<<endl;
      }
    }    
#endif

    //This part is used to fill the root tree
    //During loading a G4 track, we also need to load some thrown parameters and store
    //them into fEXEventGen. Note that G4 tree use unit of mm
    //Fix me:
    //the chain finder might not put the track it found in the same order
    //as the original G4 tree, therefore I do not know how to incert these values
    CF_X0[t]=0.;
    CF_Y0[t]=0.;
    CF_Z0[t]=fGEMCReader->z_v/10.;
    CF_Theta0[t]=3.14159/2-fGEMCReader->th_v;
    CF_Phi0[t]=fGEMCReader->phi_v ;
    CF_P0[t]=fGEMCReader->p_v/1000.;
  }

  return nhits;
}


//read a track from G4MC_RTPC12 output tree, require at least 5 hits
//please note that the G4 root tree use unit of mm 
bool EXKalManager::LoadAGEMCTrack(bool bIncludeCurveBackHits)
{
  int nhits = 0;
  double xx[MAX_HITS_PER_TRACK],yy[MAX_HITS_PER_TRACK],zz[MAX_HITS_PER_TRACK];
  while (nhits<5) {
    int n = fGEMCReader->LoadATrack();
    if( n == -1 ) {
      cout<<"Reach the end of input root file \n";
      return false;
    }
    double tmpR=0.0,tmpRmax=0.0;
    nhits=0;
    int HitNum_m = fGEMCReader->TDC->size();
    for(int i=0;i<HitNum_m;i++) {
      if(fGEMCReader->ChanID->at(i)>0) {
        //if negative means this hit is invalid

        //in case you do not want to include curve back hits!!!
        if(!bIncludeCurveBackHits) {
          tmpR = sqrt(pow(fGEMCReader->X->at(i),2.0)+pow(fGEMCReader->Y->at(i),2.0))/10.;
          if(tmpR > tmpRmax) tmpRmax = tmpR;  
          else if( tmpR+0.1 < tmpRmax) continue;
          //due to resolution, s value of some hits might be a little bit smaller 
          //than previous hit..
          //I set the margin here to be 1mm
        }
        xx[nhits] = fGEMCReader->X->at(i)/10.;
        yy[nhits] = fGEMCReader->Y->at(i)/10.;
        zz[nhits] = fGEMCReader->Z->at(i)/10.;
        nhits++;
        if(nhits >= MAX_HITS_PER_TRACK) break;
      }
    }
  }

#ifdef _EXKalManDebug_
  //just for debug
  if(_EXKalManDebug_>=3) {
    cout<<"\nGEMC Ntuple Event "<<setw(5)<<fGEMCReader->event_v<<":  HitNum_m="<<setw(2)<<fGEMCReader->TDC->size()<<endl
      <<"  P0="<<fGEMCReader->p_v<<",  Pt="<<fGEMCReader->pt_v<<", Theta0="
      <<fGEMCReader->th_v*57.3<<", Phi0="<<fGEMCReader->phi_v*57.3<<"  Z0="<<fGEMCReader->z_v/10.<<"cm"<<endl;
  }
  if(_EXKalManDebug_>=4) {
    for(int i=0;i<(int)fGEMCReader->TDC->size();i++) {
      cout<<"Hit "<<setw(2)<<i<<"("<<setw(8)<<fGEMCReader->X->at(i)<<", "
        <<setw(8)<<fGEMCReader->Y->at(i)<<", "
        <<setw(8)<<fGEMCReader->Z->at(i)<<") ==>  S="<<setw(8)<<sqrt(pow(fGEMCReader->X->at(i),2.0)+pow(fGEMCReader->Y->at(i),2.0))
        <<" mm  Phi="<<setw(8)<<atan2(fGEMCReader->Y->at(i),fGEMCReader->X->at(i))*57.3<<" deg"<<endl;
    }
  }
#endif

  bool smearing=false;
  fEventGen->MakeHitsFromTraj(xx,yy,zz,nhits,smearing,bIncludeCurveBackHits);

  //When load G4 track, we also need to load some thrown parameters and store
  //them into fEXEventGen. Note that G4 tree use unit of mm

  fEventGen->X0=0.;
  fEventGen->Y0=0.;
  fEventGen->Z0=fGEMCReader->z_v/10.;
  fEventGen->Theta0=3.14159/2-fGEMCReader->th_v;
  fEventGen->Phi0=fGEMCReader->phi_v ;
  fEventGen->P0=fGEMCReader->p_v/1000.;

  return true;
}

//read ntracks from G4MC_RTPC12 output tree and fill it into ChainFinder's hit pool.
//please note that the G4 root tree use unit of mm 
int EXKalManager::FillChainFinderHitPool(int ntracks)
{
  int nhits = 0;

  int id[MAX_HITS_PER_EVENT],tdc[MAX_HITS_PER_EVENT],adc[MAX_HITS_PER_EVENT];
  int throwntid[MAX_HITS_PER_EVENT];
  double xx[MAX_HITS_PER_EVENT],yy[MAX_HITS_PER_EVENT],zz[MAX_HITS_PER_EVENT];

  CF_ntrack_read=0;
  for(int t=0;t<ntracks;t++) {
    int n = fNtReader->LoadATrack();
    if( n == -1 ) {
      cout<<"Reach the end of input root file \n";
      return -1;
    }

    nhits=0;
    for(int i=0;i<fNtReader->HitNum_m;i++) {
      if(fNtReader->StepID_m[i]>0) {
        //if negative means this hit is invalid

        id[nhits] = fNtReader->StepID_m[i];
        tdc[nhits]= fNtReader->StepTDC_m[i];
        adc[nhits]= fNtReader->StepADC_m[i];
        xx[nhits] = fNtReader->StepX_rec_m[i]/10.;
        yy[nhits] = fNtReader->StepY_rec_m[i]/10.;
        zz[nhits] = fNtReader->StepZ_rec_m[i]/10.;
        throwntid[nhits] = t*1.0E4+i;
        throwntid[nhits] = fNtReader->StepADC_m[i]*1.0E4+i;  //just for debug
        nhits++;
        if(nhits >= MAX_HITS_PER_EVENT) break;
      }
    }
    if(nhits==0) continue;

    CF_ntrack_read++;
    if(nhits>=MIN_HITS_PER_TRACK) CF_ntrack_good++;

    //append this track into the hit pool
    int append=1;
    fChainFinder->PrepareHitPool(id, tdc, adc, xx, yy, zz, nhits, throwntid, append);


#ifdef _EXKalManDebug_
    //just for debug
    if(_EXKalManDebug_>=3) {
      cout<<"\nNtuple Event "<<setw(5)<<fNtReader->Index<<":  HitNum_m="<<setw(2)<<fNtReader->HitNum_m
        <<",  Smax="<<setw(8)<<fNtReader->Smax<<",  Smin="<<setw(8)<<fNtReader->Smin<<endl
        <<"  P0="<<fNtReader->P0_p<<",  Pt="<<fNtReader->P0_p*sin(fNtReader->Theta0_p)<<", Theta0="
        <<fNtReader->Theta0_p*57.3<<", Phi0="<<fNtReader->Phi0_p*57.3<<"  Z0="<<fNtReader->Z0/10.<<"cm"<<endl;
    }
    if(_EXKalManDebug_>=4) {
      for(int i=0;i<fNtReader->HitNum_m;i++) {
        cout<<"Hit "<<setw(2)<<i<<"("<<setw(8)<<fNtReader->StepX_rec_m[i]<<", "
          <<setw(8)<<fNtReader->StepY_rec_m[i]<<", "
          <<setw(8)<<fNtReader->StepZ_rec_m[i]<<") ==>  S="<<setw(8)<<fNtReader->StepS_rec_m[i]
        <<" mm  Phi="<<setw(8)<<fNtReader->StepPhi_rec_m[i]*57.3<<" deg"<<endl;
      }
    }
#endif

    //This part is used to fill the root tree
    //When load G4 track, we also need to load some thrown parameters and store
    //them into fEXEventGen. Note that G4 tree use unit of mm
    //Fix me:
    //the chain finder might not put the track it found in the same order
    //as the original G4 tree, therefore I do not know how to incert these values
    CF_X0[t]=fNtReader->X0/10.;
    CF_Y0[t]=fNtReader->Y0/10.;
    CF_Z0[t]=fNtReader->Z0/10.;
    CF_Theta0[t]=fNtReader->Theta0_p;
    CF_Phi0[t]=fNtReader->Phi0_p ;
    CF_P0[t]=fNtReader->P0_p;
  }

  return nhits;
}


//read a track from G4MC_RTPC12 output tree, require at least 5 hits
//please note that the G4 root tree use unit of mm 
bool EXKalManager::LoadAG4Track(bool bIncludeCurveBackHits)
{
  int nhits = 0;
  double xx[MAX_HITS_PER_TRACK],yy[MAX_HITS_PER_TRACK],zz[MAX_HITS_PER_TRACK];
  while (nhits<5) {
    int n = fNtReader->LoadATrack();
    if( n == -1 ) {
      cout<<"Reach the end of input root file \n";
      return false;
    }
    double tmpR=0.0,tmpRmax=0.0;
    nhits=0;
    for(int i=0;i<fNtReader->HitNum_m;i++) {
      if(fNtReader->StepID_m[i]>0) {
        //By Jixie @ 20160914:  do not include curve back hits!!!
        if(!bIncludeCurveBackHits) {
          tmpR = fNtReader->StepS_rec_m[i]/10.;
          if(tmpR > tmpRmax) tmpRmax = tmpR;  
          else if( tmpR+0.1 < tmpRmax) continue;
          //due to resolution, s value of some hits might be a little bit smaller 
          //than previous hit..
          //I set the margin here to be 1mm
        }
        xx[nhits] = fNtReader->StepX_rec_m[i]/10.;
        yy[nhits] = fNtReader->StepY_rec_m[i]/10.;
        zz[nhits] = fNtReader->StepZ_rec_m[i]/10.;
        nhits++;
        if(nhits >= MAX_HITS_PER_TRACK) break;
      }
    }
  }

#ifdef _EXKalManDebug_
  //just for debug
  if(_EXKalManDebug_>=3) {
    cout<<"\nNtuple Event "<<setw(5)<<fNtReader->Index<<":  HitNum_m="<<setw(2)<<fNtReader->HitNum_m
      <<",  Smax="<<setw(8)<<fNtReader->Smax<<",  Smin="<<setw(8)<<fNtReader->Smin<<endl
      <<"  P0="<<fNtReader->P0_p<<",  Pt="<<fNtReader->P0_p*sin(fNtReader->Theta0_p)<<", Theta0="
      <<fNtReader->Theta0_p*57.3<<", Phi0="<<fNtReader->Phi0_p*57.3<<"  Z0="<<fNtReader->Z0/10.<<"cm"<<endl;
  }
  if(_EXKalManDebug_>=4) {
    for(int i=0;i<fNtReader->HitNum_m;i++) {
      cout<<"Hit "<<setw(2)<<i<<"("<<setw(8)<<fNtReader->StepX_rec_m[i]<<", "
        <<setw(8)<<fNtReader->StepY_rec_m[i]<<", "
        <<setw(8)<<fNtReader->StepZ_rec_m[i]<<") ==>  S="<<setw(8)<<fNtReader->StepS_rec_m[i]
      <<" mm  Phi="<<setw(8)<<fNtReader->StepPhi_rec_m[i]*57.3<<" deg"<<endl;
    }
  }
#endif

  bool smearing=false;
  fEventGen->MakeHitsFromTraj(xx,yy,zz,nhits,smearing,bIncludeCurveBackHits);

  //When load G4 track, we also need to load some thrown parameters and store
  //them into fEXEventGen. Note that G4 tree use unit of mm

  fEventGen->X0=fNtReader->X0/10.;
  fEventGen->Y0=fNtReader->Y0/10.;
  fEventGen->Z0=fNtReader->Z0/10.;
  fEventGen->Theta0=fNtReader->Theta0_p;
  fEventGen->Phi0=fNtReader->Phi0_p ;
  fEventGen->P0=fNtReader->P0_p;

  return true;
}

//To identify which thrown track this chain corresponding to
//return the ThrownTID, also return the likelyhood, which is
//defined as occurance/total-hits
int EXKalManager::IdentifyThrownTID(int chainid, double &likelyhood) 
{
  int n = fChainFinder->fChainBuf[chainid].HitNum;
  int *thrownTID_list = new int [n];
  int *thrownTID_occur = new int [n];

  int nSource = 0;
  thrownTID_list[nSource]=fChainFinder->fChainBuf[chainid].Hits[0]->ThrownTID/1000;
  thrownTID_occur[nSource++]=1;

  //first need to check how many thrownTID exist
  for(int j=1;j<n;j++) {
    int thrownTID = fChainFinder->fChainBuf[chainid].Hits[j]->ThrownTID/1000;
    int found = 0;
    for(int t=0;t<nSource;t++) {
      if(thrownTID == thrownTID_list[t]) {
        thrownTID_occur[t]++;
        found = 1;
        break;
      }
    }
    if(!found) {
      thrownTID_list[nSource] = thrownTID;
      thrownTID_occur[nSource++]=1;
    }
  }

  //check which thrownTID dominate, store it
  int theThrownTID = thrownTID_list[0];
  int occurMax = thrownTID_occur[0];
  for(int t=1;t<nSource;t++) {
    if(occurMax < thrownTID_occur[t]) {
      occurMax = thrownTID_occur[t];
      theThrownTID = thrownTID_list[t];
    }
  }

  likelyhood = double(occurMax)/double(n);
  delete thrownTID_list;
  delete thrownTID_occur;
  return theThrownTID;
}


//run ChainFinder to search for chains then do fitting
//in each event, read multiple tracks from G4 root tree and store them into hit pool
//job := 3, no fit; 4 call GHF; 5 call KF;
int EXKalManager::RunCFNFit(int treetype, int job, int nevents, int ntracks, double max_sep, 
  double max_sep_ang, double min_sep, double min_sep_ang, double ini_sep)
{
  double xx[MAX_HITS_PER_TRACK],yy[MAX_HITS_PER_TRACK],zz[MAX_HITS_PER_TRACK];
  //note that chain finder does not sort the found chain yet, 
  //we should always include all hits
  bool bIncludeCurveBackHits=true;

  //BeginOfRun() will base on _eventtype_ to determine if or not to start NtReader
  //therefore we have to set _eventtype_ value here
  _eventtype_= 10+treetype;
  BeginOfRun(_eventtype_);

  //void SetParameters(double max_sep, double max_sep_ang, double min_sep, double min_sep_ang,
  //                   double ini_sep=1.5);
  fChainFinder->SetParameters(max_sep, max_sep_ang, min_sep, min_sep_ang, ini_sep);

  // ===================================================================
  //  Event loop
  // ===================================================================

  for (Int_t eventno = 0; eventno < nevents; eventno++) { 

#ifdef _EXKalManDebug_
    cerr << "\n------ Event " << eventno << " ------" << endl;
#endif
    eventid = eventno;  
    trackid = 0; //will be updated if _eventtype_>3

    // ============================================================
    //  Reset the buffer
    Tree_Reset_CF();
    Tree_Reset();
    fKalRTPC->Reset();

    // ===================================================================
    //execute Chain Finder
    fChainFinder->Reset();
    int ret=0;
    if(treetype==0)  ret=FillChainFinderHitPool(ntracks);   
    else  ret=FillChainFinderHitPoolGEMC(ntracks);  
    if(ret == -1) break;
    else if(ret<MIN_HITS_PER_TRACK) continue;
    fChainFinder->RemoveBadHitsFromPool();
    if(fChainFinder->fHitNum < MIN_HITS_PER_TRACK) continue;
    fChainFinder->SortHitPoolByIDTDC();
    fChainFinder->SearchChains();

    // ===================================================================
    //judge if to do fitting to the found chains
    if(job==3) {
      Tree_Fill(*(fKalRTPC->fKalTrack));
    } else {
      // ===================================================================
      //now call GHF or KF to do fitting
      int CF_npt = 0;  //number of hit in this chain
      for(int i=0;i<fChainFinder->fChainNum_Stored;i++) {
        //extract the hits in the found chain
        CF_npt = fChainFinder->fChainBuf[i].HitNum;  
        if(CF_npt>MAX_HITS_PER_TRACK) CF_npt=MAX_HITS_PER_TRACK;     
        for(int j=0;j<CF_npt;j++) {
          xx[j] = fChainFinder->fChainBuf[i].Hits[j]->X;
          yy[j] = fChainFinder->fChainBuf[i].Hits[j]->Y;
          zz[j] = fChainFinder->fChainBuf[i].Hits[j]->Z;
        }

        //only for simulation
        //to tell which thrown track this chain is
        //It is slow, do not call if only one thrown track per event 
        if(CF_ntrack_read<2) {trackid=0;CF_ThrownTID_like=1.0;}
        else trackid = IdentifyThrownTID(i, CF_ThrownTID_like); 

        if(job==4) {
          // ============================================================
          //  call global helix fit
          // ============================================================
          fKalRTPC->DoGlobalHelixFit(xx,yy,zz,CF_npt,bIncludeCurveBackHits);       
        } else if(job==5) {
          // ============================================================
          //  call KalmanFilter
          // ============================================================

          //  Generate a partcle and Swim the particle in fDetector
          bool bTrackReady = false;
          bool smearing = false; 
          //bTrackReady = fKalRTPC->PrepareATrack(xx,yy,zz,CF_npt,smearing,bIncludeCurveBackHits);
          bTrackReady = fKalRTPC->PrepareATrack(xx,yy,zz,CF_npt,smearing,false);
          if(!bTrackReady) continue;

          //Remove backward hits also judge whether or not need 2nd iteration
          bool bRemoveBackwardHits=true;
          bool bNeed2Iter = fKalRTPC->JudgeFor2ndIteration(bRemoveBackwardHits);
          if(fKalRTPC->fKalHits_Forward->GetEntriesFast()<MIN_HITS_PER_TRACK) continue;

          //Do KalmanFilter
          fKalRTPC->DoFitAndFilter(bNeed2Iter);
        }

        // ============================================================
        //fill output tree
        // ============================================================
        Tree_Fill(*(fKalRTPC->fKalTrack));
        
        
        // ============================================================
        //  Very Primitive Event Display
        // ============================================================
        EventVisulization();
        
        //  Reset the buffer
        fKalRTPC->Reset();
        Tree_Reset();
      } //end of for(int i=0;i<fChainFinder->fChainNum;i++)
    } //end of if(job==3)

#ifdef _EXKalManDebug_
    Stop4Debug(_index_);
#endif
  } //end of envet loop

  //ROOT Tree will be written and closed inside EndOfRun();
  EndOfRun();
  return 1;
}

void EXKalManager::EventVisulization2()
{
  // ============================================================
  //  Very Primitive Event Display
  // ============================================================

  static TCanvas *cvp = 0;
  if (!gROOT->IsBatch()) {
    if (!cvp) {
      cvp = new TCanvas("OED", "Event Display", 10, 10, 610, 610);
    } else {
      cvp->cd();cvp->Clear();
    }

    //TView   *vwp = new TView(1);
    TView *vwp = TView::CreateView(1,0,0);
    double xmin,ymin,zmin,xmax,ymax,zmax;
    Int_t ierr;
#if 0
    xmin=ymin=zmin=-kRTPC_Length/2.0;
    xmax=ymax=zmax=kRTPC_Length/2.0;
    vwp->SetRange(xmin,ymin,zmin,xmax,ymax,zmax);
    vwp->SetView(10.,80.,80.,ierr);   
#else
    xmin=ymin=zmin=-(kRTPC_R_GEM1+1.0);
    xmax=ymax=zmax=kRTPC_R_GEM1+1.0;
    vwp->SetRange(xmin,ymin,zmin,xmax,ymax,zmax);
    //vwp->SetView(0.,0.,0.,ierr);  //side view      
    vwp->SetView(0.,0.,90.,ierr);   //beam view
#endif

    fKalRTPC->fDetector->Draw(40);
    fKalRTPC->fKalTrack->Draw(2,"");         

    cout << "Next? [y(yes)/n(no, abort)/e(edit)/q(quit)] " << flush;
    static const Int_t kMaxLen = 255;
    Char_t temp[kMaxLen];
    cin.getline(temp,kMaxLen);
    TString opts(temp);
    opts.ToLower();
    if (!opts.Length()) { /*continue*/;}
    else if (opts[0] == 'n' || opts[0] == 'q') { abort(); } 
    else if (opts[0] == 'e') {
      cout << "Select \"Quit ROOT\" from \"File\" to display next" << endl;
      cout << "\"CTRL+C\" to really quit" << endl;
      fApp->Run(kTRUE);
    }
  } // end of event display
}

void EXKalManager::EventVisulization()
{
  // ============================================================
  //  Very Primitive Event Display
  // ============================================================

  static TCanvas *cvp = 0;
  if (!gROOT->IsBatch()) {
    if (!cvp) {
      cvp = new TCanvas("OED", "Event Display", 10, 10, 410, 810);
    } else {
      cvp->cd();cvp->Clear();
      cvp->Divide(1,2);
    }
    
    cvp->cd();
    //TView   *vwp = new TView(1);
    double xmin,ymin,zmin,xmax,ymax,zmax;
    Int_t ierr;
    
    cvp->cd(1);
    TView *vwp1 = TView::CreateView(1,0,0);
    //TView3D   *vwp1 = new TView3D();
    xmin=ymin=zmin=-kRTPC_Length/2.0;
    xmax=ymax=zmax=kRTPC_Length/2.0;
    vwp1->SetRange(xmin,ymin,zmin,xmax,ymax,zmax);
    vwp1->SetView(10.,80.,80.,ierr);   
    //vwp1->SetView(0.,90.,90.,ierr);  //side view  
    fKalRTPC->fDetector->Draw(40);
    fKalRTPC->fKalTrack->Draw(2,"");   

    cvp->cd(2);
    TView *vwp2 = TView::CreateView(1,0,0);
    //TView3D   *vwp2 = new TView3D();
    xmin=ymin=zmin=-(kRTPC_R_GEM1+1.0);
    xmax=ymax=zmax=kRTPC_R_GEM1+1.0;
    vwp2->SetRange(xmin,ymin,zmin,xmax,ymax,zmax);
    vwp2->SetView(0.,0.,90.,ierr);   //beam view
    fKalRTPC->fDetector->Draw(40);
    fChainFinder->DrawPool();
    fChainFinder->DrawChain();       
    fKalRTPC->fKalTrack->Draw(2,"same");  
  
    cvp->cd(0);
    cout << "Next? [y(yes)/n(no, abort)/e(edit)/q(quit)] " << flush;
    static const Int_t kMaxLen = 255;
    Char_t temp[kMaxLen];
    cin.getline(temp,kMaxLen);
    TString opts(temp);
    opts.ToLower();
    if (!opts.Length()) { /*continue*/;}
    else if (opts[0] == 'n' || opts[0] == 'q') { abort(); } 
    else if (opts[0] == 'e') {
      cout << "Select \"Quit ROOT\" from \"File\" to display next" << endl;
      cout << "\"CTRL+C\" to really quit" << endl;
      fApp->Run(kTRUE);
    }
  } // end of event display
}

//test the KF only
//cerr << "Usage: "<<argv[0] <<" <job=0|1|2> <nevent> [pt_min_gev=0.1] [pt_max_gev=0.1]" << endl;
//cerr << "\t  job: 0 generate helix, 1 loadtrack from geant4 root file, 2 generate circle\n";
//cerr << "\t  nevents: number of events to generate \n";
//cerr << "\t  pt_min_gev and pt_max_gev: specifiy the range of pt in Gev \n";
//cerr << "\t  Note that if pt is negative then anti-clockwise track will be generated \n";
int EXKalManager::RunKF(int job, int nevents, double pt_min, double pt_max, double costh_min, 
  double costh_max, double z_min, double z_max)
{
  //User will fill BeginOfRun() and EndOfRun();
  //ROOT tree will be opened and initilized inside BeginOfRun()

  //BeginOfRun() will base on _eventtype_ to determine if or not to start NtReader
  //therefore we have to set _eventtype_ value here
  _eventtype_= job;
  BeginOfRun(_eventtype_);

  // ===================================================================
  //  Event loop
  // ===================================================================

  for (Int_t eventno = 0; eventno < nevents; eventno++) { 
#ifdef _EXKalManDebug_
    cerr << "\n------ Event " << eventno << " ------" << endl;
#endif
    //I try to mimic one event contains multiple track,
    //currently these value is set in this way. 
    //The user can manipulate them later


    eventid=eventno;
    trackid=0;
    ntrack=1;

    // ============================================================
    //  Reset the buffer
    // ============================================================
    fKalRTPC->Reset();
    Tree_Reset();

    // ============================================================
    //  Generate a partcle and Swim the particle in fDetector
    // ============================================================
    bool bTrackReady = false;
    bool bIncludeCurveBackHits = true;  //this only work for EXEventGen

    if(_eventtype_ == 0 || _eventtype_ == 2)
      bTrackReady = fKalRTPC->PrepareATrack(_eventtype_, pt_min, pt_max, costh_min, 
      costh_max, z_min, z_max, bIncludeCurveBackHits);
    else if(_eventtype_ == 1)
      bTrackReady=LoadAG4Track(bIncludeCurveBackHits);
    else {
      cout<<"This _eventtype_("<<_eventtype_<<") is not yet supported! I quit..."<<endl;
      exit(-3);
    }
    if(!bTrackReady) break;

    //Remove backward hits also judge whether or not need 2nd iteration
    bool bRemoveBackwardHits=true;
    bool bNeed2Iter = fKalRTPC->JudgeFor2ndIteration(bRemoveBackwardHits);
    if(fKalRTPC->fKalHits_Forward->GetEntriesFast()<MIN_HITS_PER_TRACK) continue;

    // ============================================================
    //  Do KalmanFilter
    // ============================================================
    fKalRTPC->DoFitAndFilter(bNeed2Iter);


    // ============================================================
    //  Fill the KalmanFilter result into the root tree
    // ============================================================
    Tree_Fill(*(fKalRTPC->fKalTrack));

    // ============================================================
    //  Very Primitive Event Display
    // ============================================================
    EventVisulization();
   

#ifdef _EXKalManDebug_
    Stop4Debug(_index_);
#endif
  } // end of event loop

  //ROOT Tree will be written and closed inside EndOfRun();
  EndOfRun();
  return 1;
}


//I separate filling the tree into an individual subroutine since 
//filling the tree will not be necessary when provided to CLAS12 software
void EXKalManager::Tree_Fill(EXHYBTrack &kaltrack)
{
  //Note that all variable in NtReader are in unit of mm
  //all come from KalRTPC are in unit of cm

  //This block is for chain-finder
  //It will be called CF_ChainNum times in a given event
  //I do not care the speed for this class, since it is just for illustration
  CF_HitNum = fChainFinder->fHitNum;
  ntrack = fChainFinder->fChainNum_Stored;
  CF_ChainNum = fChainFinder->fChainNum_Stored;
  for(int i=0;i<CF_HitNum;i++) {
    CF_ID[i]=fChainFinder->fHitPool[i].ID;
    CF_TDC[i]=fChainFinder->fHitPool[i].TDC;
    CF_ADC[i]=fChainFinder->fHitPool[i].ADC;
    CF_X[i]=fChainFinder->fHitPool[i].X;
    CF_Y[i]=fChainFinder->fHitPool[i].Y;
    CF_Z[i]=fChainFinder->fHitPool[i].Z;
    CF_S[i]=fChainFinder->fHitPool[i].S;
    CF_Phi[i]=fChainFinder->fHitPool[i].Phi;
    CF_Status[i]=fChainFinder->fHitPool[i].Status;
    CF_ThrownTID[i]=fChainFinder->fHitPool[i].ThrownTID;
    CF_ChainInfo[i]=fChainFinder->fHitPool[i].ChainInfo;
  }


  //Load the thrown parameters, most of them available from fEvenGen
  if(_eventtype_>=3) {
    _x0_=CF_X0[trackid];
    _y0_=CF_Y0[trackid];
    _z0_=CF_Z0[trackid];
    p0=CF_P0[trackid];
    th0=CF_Theta0[trackid];
    ph0=CF_Phi0[trackid]; 
  }else {
    _x0_=fEventGen->X0;
    _y0_=fEventGen->Y0;
    _z0_=fEventGen->Z0;
    p0=fEventGen->P0;
    th0=fEventGen->Theta0;
    ph0=fEventGen->Phi0; 
  }
  if(ph0> kPi) ph0-=2*kPi;
  if(ph0<-kPi) ph0+=2*kPi;
  pt0=p0*sin(th0);
  pz0=p0*cos(th0);

  npt0=fEventGen->StepNum;

  step_s_min = 99999.0;
  step_s_max =-99999.0;
  for(int i=0;i<npt0;i++) {
    step_x[i]=fEventGen->StepX[i];
    step_y[i]=fEventGen->StepY[i];
    step_z[i]=fEventGen->StepZ[i];
    step_phi[i]=fEventGen->StepPhi[i];
    step_s[i]=fEventGen->StepS[i];
    if(step_s_min>step_s[i]) step_s_min=step_s[i];
    if(step_s_max<step_s[i]) step_s_max=step_s[i];
  }

  //Load kalman filter result from EXKalRTPC
  p_rec =fKalRTPC->P_rec;
  th_rec=fKalRTPC->Theta_rec;
  ph_rec=fKalRTPC->Phi_rec; 
  pt_rec=p_rec*sin(th_rec);
  pz_rec=p_rec*cos(th_rec);
  x_rec =fKalRTPC->X_rec;
  y_rec =fKalRTPC->Y_rec;
  z_rec =fKalRTPC->Z_rec;
  r_rec =fKalRTPC->R_rec;
  a_rec =fKalRTPC->A_rec;
  b_rec =fKalRTPC->B_rec;
  ndf   =fKalRTPC->NDF;
  chi2  =fKalRTPC->Chi2;
  cl    =TMath::Prob(chi2, ndf);

  //copy the status array
  npt = fKalRTPC->fKalHits_Forward->GetEntriesFast(); 
  for(int i=0;i<npt;i++) step_status[i]=fKalRTPC->step_status[i];

  //load various sites
  TIter next(fKalRTPC->fKalHits_Forward, kIterBackward);   
  //Note that the first site of kaltrack is a dummy site
  //do not include it into the output root tree
  int nn=0;
  int iSite=kaltrack.GetEntriesFast()-1;  //site 0 is dummy site, start from 1
  EXHit *hitp = dynamic_cast<EXHit *>(next());  
  while (hitp) {     // loop over hits    
    //if(nn<0) {cout<<"something wrong!!! npt="<<npt<<", nsite="<<kaltrack.GetEntriesFast()-1<<"\n";break;}
    //fill the global variables for root tree
    //TVector3 x_rec = hitp->GetRawXv();  //I currently use raw info to debug step
    const EXMeasLayer &ml = dynamic_cast<const EXMeasLayer &>(hitp->GetMeasLayer());
    TVector3 x_rec = ml.HitToXv(*hitp); 
    step_x_rec[nn]=x_rec.X();step_y_rec[nn]=x_rec.Y();step_z_rec[nn]=x_rec.Z();
    step_phi_rec[nn]=x_rec.Phi();step_s_rec[nn]=x_rec.Perp();

    //filtered state vector exist only if it is not removed
    if(step_status[nn]) {
      TKalTrackSite *site = (TKalTrackSite *)kaltrack.At(iSite--); // load the site

      TVKalState *state_exp = (TVKalState*) &(site->GetState(TVKalSite::kPredicted));
      THelicalTrack hel_exp = (dynamic_cast<TKalTrackState *>(state_exp))->GetHelix();
      TVector3 x_exp=hel_exp.CalcXAt(0.0);
      step_x_exp[nn]=x_exp.X();step_y_exp[nn]=x_exp.Y();step_z_exp[nn]=x_exp.Z();
      step_phi_exp[nn]=x_exp.Phi();step_s_exp[nn]=x_exp.Perp();

      //TVKalState *state_fil = (TVKalState*) &(site->GetCurState());
      TVKalState *state_fil = (TVKalState*) &(site->GetState(TVKalSite::kFiltered));
      THelicalTrack hel_fil = (dynamic_cast<TKalTrackState *>(state_fil))->GetHelix();
      TVector3 x_fil=hel_fil.CalcXAt(0.0);
      step_x_fil[nn]=x_fil.X();step_y_fil[nn]=x_fil.Y();step_z_fil[nn]=x_fil.Z();
      step_phi_fil[nn]=x_fil.Phi();step_s_fil[nn]=x_fil.Perp();
      //cout<<"good hit "<<nn<<"/"<<npt<<endl;
    }
    else {
      step_x_exp[nn]=step_y_exp[nn]=step_z_exp[nn]=step_phi_exp[nn]=step_s_exp[nn]=0.0;
      step_x_fil[nn]=step_y_fil[nn]=step_z_fil[nn]=step_phi_fil[nn]=step_s_fil[nn]=0.0;
      //cout<<"*bad* hit "<<nn<<"/"<<npt<<endl;
    }
    hitp = dynamic_cast<EXHit *>(next());
    nn++;
  }

  //store global helix result before my correction 
  r_hel_raw = fKalRTPC->R_hel_raw;
  th_hel_raw= fKalRTPC->Theta_hel_raw;
  ph_hel_raw= fKalRTPC->Phi_hel_raw; 
  a_hel_raw = fKalRTPC->A_hel_raw;
  b_hel_raw = fKalRTPC->B_hel_raw;
  z_hel_raw = fKalRTPC->Z_hel_raw;

  //store global helix result after my correction 
  p_hel =fKalRTPC->P_hel;
  th_hel=fKalRTPC->Theta_hel;
  ph_hel=fKalRTPC->Phi_hel; 
  pt_hel=p_hel*sin(th_hel);
  pz_hel=p_hel*cos(th_hel);
  x_hel =fKalRTPC->X_hel;
  y_hel =fKalRTPC->Y_hel;
  z_hel =fKalRTPC->Z_hel;
  r_hel =fKalRTPC->R_hel;
  a_hel =fKalRTPC->A_hel;
  b_hel =fKalRTPC->B_hel;
  dca_hel =fKalRTPC->DCA_hel;
  chi2_hel =fKalRTPC->Chi2_hel;

  p_3pt =fKalRTPC->P_3pt;
  th_3pt=fKalRTPC->Theta_3pt;
  pt_3pt=fKalRTPC->Pt_3pt;
  r_3pt =fKalRTPC->R_3pt;
  a_3pt =fKalRTPC->A_3pt;
  b_3pt =fKalRTPC->B_3pt;


  rho_kal_ini  =fKalRTPC->rho_kal_ini;
  tnl_kal_ini  =fKalRTPC->tnl_kal_ini;
  phi0_kal_ini =fKalRTPC->phi0_kal_ini;

  rho_1st  =fEventGen->Rho_1st;
  tnl_1st  =fEventGen->TanLambda_1st;
  phi0_1st =fEventGen->Phi0_1st;
  rho_last =fEventGen->Rho_last;
  tnl_last =fEventGen->TanLambda_last;
  phi0_last=fEventGen->Phi0_last;

  fTree->Fill();
  if(!((_index_+1)%100)) fFile->Write("",TObject::kOverwrite);
  _index_++;
}

