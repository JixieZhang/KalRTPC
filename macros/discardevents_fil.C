//this script is used to study events that discarded by Kalman Filter
#include "stdlib.h"
#include <iostream>
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

void discardevents_old()
{
  TTree *t=(TTree*) gROOT->FindObject("t");
  t->Draw(">>elist1","step_status==0");
  TEventList *elist1 = (TEventList*)gDirectory->Get("elist1");
  
  t->SetEventList(elist1);
 
  TCanvas *c20=new TCanvas("c20","discard events",800,800);
  TH2F *h2Frame=new TH2F("h2Frame","Hits: Accepted(blue), Discarded(red); x(cm); y(cm)",80,-8,8,80,-8,8);
  h2Frame->Draw();
  t->Draw("step_y:step_x","step_status==0","same");
  TGraph *gr0=(TGraph*) (gROOT->FindObject("Graph")->Clone("gr0"));
  gr0->SetTitle("discraded hits(red)");
  gr0->SetMarkerColor(2);
  gr0->SetMarkerStyle(2);

  c20->Clear();
  h2Frame->Draw();
  t->Draw("step_y:step_x","step_status==1","same");
  TGraph *gr1=(TGraph*) (gROOT->FindObject("Graph")->Clone("gr1"));
  gr1->SetTitle("accepted hits(blue)");
  gr1->SetMarkerColor(4);
  gr1->SetMarkerStyle(4);

  c20->Clear();
  h2Frame->Draw();
  gr0->Draw("psame");
  gr1->Draw("psame");
  c20->SaveAs("_discard.png");

}

void discardevents_fil()
{
  TTree *t=(TTree*) gROOT->FindObject("t");
  t->Draw(">>elist1","step_status==0");
  TEventList *elist1 = (TEventList*)gDirectory->Get("elist1");
  
  t->SetEventList(elist1);
 
  TCanvas *c20=new TCanvas("c20","discard events",800,800);
  TH2F *h2Frame=new TH2F("h2Frame","Hits: Accepted(blue), Discarded(red); x(cm); y(cm)",80,-8,8,80,-8,8);
  h2Frame->Draw();
  t->Draw("step_y:step_x","step_status==0","same");
  TGraph *gr0=(TGraph*) (gROOT->FindObject("Graph")->Clone("gr0"));
  gr0->SetTitle("discraded hits(red)");
  gr0->SetMarkerColor(2);
  gr0->SetMarkerStyle(2);

  c20->Clear();
  h2Frame->Draw();
  t->Draw("step_y:step_x","step_status==1","same");
  TGraph *gr1=(TGraph*) (gROOT->FindObject("Graph")->Clone("gr1"));
  gr1->SetTitle("accepted hits(blue)");
  gr1->SetMarkerColor(4);
  gr1->SetMarkerStyle(4);

  c20->Clear();
  h2Frame->Draw();
  t->Draw("step_y_fil:step_x_fil","step_status==1","same");
  TGraph *gr_fil=(TGraph*) (gROOT->FindObject("Graph")->Clone("gr_fil"));
  gr_fil->SetTitle("optimized hits(purple)");
  gr_fil->SetMarkerColor(6);
  gr_fil->SetMarkerStyle(22);

  c20->Clear();
  h2Frame->Draw();
  gr0->Draw("psame");
  gr1->Draw("psame");
  gr_fil->Draw("psame");
  c20->SaveAs("_discard.png");

}
