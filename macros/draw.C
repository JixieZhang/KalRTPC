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

void draw(int log_hel=0,int log_rec=1, const char *key="")
{
  double mean_dZ,sigma_dZ;
  double mean_dT,sigma_dT;
  double mean_dF,sigma_dF;
  double mean_dP,sigma_dP;

  ofstream fout;
  fout.open("table.txt",ios_base::app);

  TTree *t=(TTree*) gROOT->FindObject("t"); 
  //TTree *t = (TTree*)gDirectory->Get("t");
  //////////////////////////////////////////////////////////////
  TCanvas *c13 = new TCanvas("c13","",900,700);
  c13->Divide(2,2);
  c13->cd(1);
  t->Draw("z0-z_hel>>hd0hz(80,-1,1)","","");
  TH1F *hd0hz = (TH1F*) gROOT->FindObject("hd0hz");
  FitGauss(hd0hz,mean_dZ,sigma_dZ);

  c13->cd(2);
  t->Draw("th0-th_hel>>hd0ht(40,-0.2,0.2)","","");
  TH1F *hd0ht = (TH1F*) gROOT->FindObject("hd0ht");
  FitGauss(hd0ht,mean_dT,sigma_dT);

  c13->cd(3);
  t->Draw("ph0-ph_hel>>hd0hph(100,-0.5,0.5)","","");
  TH1F *hd0hph = (TH1F*) gROOT->FindObject("hd0hph");
  FitGauss(hd0hph,mean_dF,sigma_dF);

  c13->cd(4);
  t->Draw("pt0-pt_hel>>hd0hp(100,-0.05,0.05)","","");
  TH1F *hd0hp = (TH1F*) gROOT->FindObject("hd0hp");
  FitGauss(hd0hp,mean_dP,sigma_dP);

  if(log_hel) {
    fout<<"\n"<<setw(40)<<"FileName"<<"  "
	<<setw(18)<<"dZ(cm)"<<"  "
	<<setw(18)<<"dTheta(rad)"<<"  "
	<<setw(18)<<"dPhi(rad)"<<"  "
	<<setw(19)<<"dPt(GeV/c)"<<endl;
    fout<<setw(40)<<"global_helix_fit"<<fixed<<setprecision(3)<<"  "
	<<setw(8)<<mean_dZ<<" +/- " <<setw(5)<<sigma_dZ<<"  "
	<<setw(8)<<mean_dT<<" +/- " <<setw(5)<<sigma_dT<<"  "
	<<setw(8)<<mean_dF<<" +/- " <<setw(5)<<sigma_dF<<"  "
	<<setprecision(4)
	<<setw(8)<<mean_dP<<" +/- " <<setw(6)<<sigma_dP<<endl;
  }
  c13->SaveAs(Form("thrown_vs_hel_%s.png",key));

  //////////////////////////////////////////////////////////////
  TCanvas *c12 = new TCanvas("c12","",900,700);
  c12->Divide(2,2);
  c12->cd(1);
  t->Draw("z0-z_rec>>hd0z(80,-1,1)","","");
  TH1F *hd0z = (TH1F*) gROOT->FindObject("hd0z");
  FitGauss(hd0z,mean_dZ,sigma_dZ);

  c12->cd(2);
  t->Draw("th0-th_rec>>hd0t(40,-0.2,0.2)","","");
  TH1F *hd0t = (TH1F*) gROOT->FindObject("hd0t");
  FitGauss(hd0t,mean_dT,sigma_dT);

  c12->cd(3);
  t->Draw("ph0-ph_rec>>hd0ph(100,-0.5,0.5)","","");
  TH1F *hd0ph = (TH1F*) gROOT->FindObject("hd0ph");
  FitGauss(hd0ph,mean_dF,sigma_dF);

  c12->cd(4);
  t->Draw("pt0-pt_rec>>hd0p(100,-0.05,0.05)","","");
  TH1F *hd0p = (TH1F*) gROOT->FindObject("hd0p");
  FitGauss(hd0p,mean_dP,sigma_dP);
 
  c12->SaveAs(Form("thrown_vs_rec_%s.png",key));

  if(log_rec) {
    fout<<setw(40)<<gFile->GetName()<<fixed<<setprecision(3)<<"  "
	<<setw(8)<<mean_dZ<<" +/- " <<setw(5)<<sigma_dZ<<"  "
	<<setw(8)<<mean_dT<<" +/- " <<setw(5)<<sigma_dT<<"  "
	<<setw(8)<<mean_dF<<" +/- " <<setw(5)<<sigma_dF<<"  "
	<<setprecision(4)
	<<setw(8)<<mean_dP<<" +/- " <<setw(6)<<sigma_dP<<endl;
  }
  fout.close();
  
  //////////////////////////////////////////////////////////////
  TCanvas *c11 = new TCanvas("c11","",900,700);
  c11->Divide(2,2);
  c11->cd(1);
  t->Draw("z_hel-z_rec>>hdz(80,-1,1)","","");
  TH1F *hdz = (TH1F*) gROOT->FindObject("hdz");
  FitGauss(hdz);

  c11->cd(2);
  t->Draw("th_hel-th_rec>>hdt(40,-0.2,0.2)","","");
  TH1F *hdt = (TH1F*) gROOT->FindObject("hdt");
  FitGauss(hdt);

  c11->cd(3);
  t->Draw("ph_hel-ph_rec>>hdph(100,-0.5,0.5)","","");
  TH1F *hdph = (TH1F*) gROOT->FindObject("hdph");
  FitGauss(hdph);

  c11->cd(4);
  t->Draw("pt_hel-pt_rec>>hdp(100,-0.05,0.05)","","");
  TH1F *hdp = (TH1F*) gROOT->FindObject("hdp");
  FitGauss(hdp);

  c11->SaveAs(Form("hel_vs_rec_%s.png",key));

}
