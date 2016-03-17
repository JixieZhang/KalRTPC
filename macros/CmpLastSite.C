//this script is used to compare kalman filter result with global
//helix fit result
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <iomanip>

#include "math.h"
using namespace std;

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TQObject.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TString.h"
#include "TCut.h"
#include "TCutG.h"
#include "TPaveText.h"
#include "TText.h"
#include "TPad.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

void FitGauss(TH1F *h1, double &mean, double &sigma)
{
  if(!h1) return;
  h1->Fit("gaus","","Q");
  TF1 *f = (TF1 *)h1->GetListOfFunctions()->FindObject("gaus");
  mean=f->GetParameter(1);
  sigma=f->GetParameter(2);
  
  h1->Fit("gaus","","R+",mean-1.5*sigma,mean+1.5*sigma);
  f = (TF1 *)h1->GetListOfFunctions()->FindObject("gaus");
  mean=f->GetParameter(1);
  sigma=f->GetParameter(2);
  h1->Draw();
}

void FitGauss(TH1F *h1)
{
  if(!h1) return;
  h1->Fit("gaus","","Q");
  TF1 *f = (TF1 *)h1->GetListOfFunctions()->FindObject("gaus");
  double mean=f->GetParameter(1);
  double sigma=f->GetParameter(2);
  
  h1->Fit("gaus","","R+",mean-1.5*sigma,mean+1.5*sigma);
  h1->Draw();
}

void CmpLastSite(){
  double mean_dZ,sigma_dZ;
  double mean_dT,sigma_dT;
  double mean_dF,sigma_dF;
  double mean_dP,sigma_dP;

  TTree *t=(TTree*) gROOT->FindObject("t"); 
  //TTree *t = (TTree*)gDirectory->Get("t");
  //////////////////////////////////////////////////////////////
  TCanvas *c12 = new TCanvas("c12","",900,900);
  c12->Divide(3,3);
  c12->cd(1);
  t->Draw("step_x-step_x_exp>>hd0x_exp(50,-1,1)","","");
  TH1F *hd0x_exp = (TH1F*) gROOT->FindObject("hd0x_exp");
  FitGauss(hd0x_exp);

  c12->cd(2);
  t->Draw("step_x-step_x_fil>>hd0x_fil(50,-1,1)","","");
  TH1F *hd0x_fil = (TH1F*) gROOT->FindObject("hd0x_fil");
  FitGauss(hd0x_fil);

  c12->cd(3);
  t->Draw("step_x_fil-step_x_exp>>hdx_exp_fil(50,-1,1)","","");
  TH1F *hdx_exp_fil = (TH1F*) gROOT->FindObject("hdx_exp_fil");
  FitGauss(hdx_exp_fil);

 
  c12->cd(4);
  t->Draw("step_y-step_y_exp>>hd0y_exp(50,-1,1)","","");
  TH1F *hd0y_exp = (TH1F*) gROOT->FindObject("hd0y_exp");
  FitGauss(hd0y_exp);

  c12->cd(5);
  t->Draw("step_y-step_y_fil>>hd0y_fil(50,-1,1)","","");
  TH1F *hd0x_fil = (TH1F*) gROOT->FindObject("hd0y_fil");
  FitGauss(hd0y_fil);

  c12->cd(6);
  t->Draw("step_y_fil-step_y_exp>>hdy_exp_fil(50,-1,1)","","");
  TH1F *hdy_exp_fil = (TH1F*) gROOT->FindObject("hdy_exp_fil");
  FitGauss(hdx_exp_fil);


  c12->cd(7);
  t->Draw("step_z-step_z_exp>>hd0z_exp(50,-1,1)","","");
  TH1F *hd0z_exp = (TH1F*) gROOT->FindObject("hd0z_exp");
  FitGauss(hd0z_exp);

  c12->cd(8);
  t->Draw("step_z-step_z_fil>>hd0z_fil(50,-1,1)","","");
  TH1F *hd0z_fil = (TH1F*) gROOT->FindObject("hd0z_fil");
  FitGauss(hd0z_fil);

  c12->cd(9);
  t->Draw("step_z_fil-step_z_exp>>hdz_exp_fil(50,-1,1)","","");
  TH1F *hdz_exp_fil = (TH1F*) gROOT->FindObject("hdz_exp_fil");
  FitGauss(hdz_exp_fil);

  c12->SaveAs("delta_lastsite.png");

  //t->Scan("step_x:step_x_fil:step_y:step_y_fil:step_z:step_z_fil:step_status")

}
