//Draw tracking efficiency
#include <iostream>
#include "string.h"
#include "math.h"

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"

using namespace std;

void TrackingEff()
{
  TChain *ep=new TChain("ep");
  ep->Add("nt_P-10T-10_Step1mm_SingleTDC_51tracks.root");
  TChain *t=new TChain("t");
  t->Add("h.root");
  
  TH1F *h1p_0=new TH1F("h1p_0"," ;P (GeV/c) ; ",30,0.05,0.35);
  TH1F *h1p_1=new TH1F("h1p_1"," ;P (GeV/c) ; ",30,0.05,0.35);
  h1p_0->SetLineColor(1);
  h1p_1->SetLineColor(2);
  
  TH1F *h1th_0=new TH1F("h1th_0"," ;#theta (deg) ; ",45,0.,180.);
  TH1F *h1th_1=new TH1F("h1th_1"," ;#theta (deg) ; ",45,0.,180.);
  h1th_0->SetLineColor(1);
  h1th_1->SetLineColor(2);

  TH2F *h2thp_0=new TH2F("h1thp_0"," ;P (GeV/c) ;#theta (deg) ",30,0.05,0.35,45,0.,180.);
  TH2F *h2thp_1=new TH2F("h1thp_1"," ;P (GeV/c) ;#theta (deg) ",30,0.05,0.35,45,0.,180.);
  //h2thp_0->SetMarkerColor(1);
  //h2thp_1->SetMarkerColor(2);

  TCut cut0 = "ShiftTDC==0 && Smin<50. && Smax>50. && HitNum_m>=5 && Smax-Smin>20.";
  TCut cut1 = "shifttime==0";


  ep->Project("h1p_0","P0_p",cut0);
  t->Project("h1p_1","p0",cut1);
  TH1F* h1p_eff = (TH1F*)h1p_1->Clone("h1p_eff");
  h1p_eff->Divide(h1p_0);
  h1p_eff->Scale(100);
  h1p_eff->SetTitle("Chain Finder efficiency (%)");
  h1p_eff->GetXaxis()->SetNdivisions(10,5,0);
  h1p_eff->GetYaxis()->SetNdivisions(10,5,0);

  ep->Project("h1th_0","Theta0_p*57.3",cut0);
  t->Project("h1th_1","th0*57.3",cut1);
  TH1F* h1th_eff = (TH1F*)h1th_1->Clone("h1th_eff");
  h1th_eff->Divide(h1th_0);
  h1th_eff->Scale(100);
  h1th_eff->SetTitle("Chain Finder efficiency (%)");
  h1th_eff->GetXaxis()->SetNdivisions(10,5,0);
  h1th_eff->GetYaxis()->SetNdivisions(10,5,0);
  
  ep->Project("h1thp_0","Theta0_p*57.3:P0_p",cut0);
  t->Project("h1thp_1","th0*57.3:p0",cut1);
  TH2F* h2thp_eff = (TH2F*)h2thp_1->Clone("h2thp_eff");
  h2thp_eff->Divide(h2thp_0);
  h2thp_eff->Scale(100);
  h2thp_eff->SetTitle("Chain Finder efficiency (%)");

  
  TCanvas *c1=new TCanvas("c1","",900,700);
  c1->Divide(1,2);
  
  c1->cd(1);gPad->SetGrid(1,1);
  h1p_eff->Draw();
  
  c1->cd(2);gPad->SetGrid(1,1);
  h1th_eff->Draw();

  c1->SaveAs("Tracking_Eff_1D.png");
  
  TCanvas *c2=new TCanvas("c2","",800,600);
  c2->cd();
  c2->SetRightMargin(0.12);
  h2thp_eff->SetMinimum(50);
  h2thp_eff->Draw("TEXT colz");
  c2->SaveAs("Tracking_Eff_2D.png");
}
