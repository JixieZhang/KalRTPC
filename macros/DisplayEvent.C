//This script is to plot the trajectory of RTPC events
//
#include <iostream>
using namespace std;

#include <TQObject.h>
#include "TROOT.h"
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TVirtualFitter.h>
#include <TMinuit.h>
#include <TMath.h>

#include "TLegend.h"
#include "TString.h"
#include "TCut.h"
#include "TCutG.h"
#include "TPaveText.h"
#include "TText.h"
#include "TSystem.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include <TChain.h>
#include "TLine.h"

TH2F *h2yxframe = 0;
TH2F *h2xzframe = 0;
TH2F *h2yzframe = 0;
TH2F *h2rzframe = 0;
TCanvas *c1,*c2;
int  thisevent = 0;


void Init()
{
  c1 = new TCanvas("c1","RTPC track display",0,20,800,800);
  c2 = new TCanvas("c2","RTPC track display",900,20,600,600);
  
  h2yxframe=new TH2F("yxframe","y-x",
		     200,-10,10,200,-10,10);

  h2yxframe->SetXTitle("X (mm)  ");
  h2yxframe->SetYTitle("Y (mm)  ");
  h2yxframe->SetMarkerStyle(1);
  h2yxframe->SetMarkerColor(3);

  for(int i=0;i<3600;i++)
    {
      double phi = i*2*3.14159/3600;
      h2yxframe->Fill(2*cos(phi),2*sin(phi));
      h2yxframe->Fill(3*cos(phi),3*sin(phi));
      h2yxframe->Fill(7*cos(phi),7*sin(phi));
      h2yxframe->Fill(8*cos(phi),8*sin(phi));
    }

  h2yzframe=new TH2F("yzframe","y-z",
		     400,-20,20,200,-10,10);

  h2yzframe->SetXTitle("Z (mm)  ");
  h2yzframe->SetYTitle("Y (mm)  ");
  h2yzframe->SetMarkerStyle(1);
  h2yzframe->SetMarkerColor(3);
  
  
  h2xzframe=new TH2F("xzframe","x-z",
		     400,-20,20,200,-10,10);

  h2xzframe->SetXTitle("Z (mm)  ");
  h2xzframe->SetYTitle("X (mm)  ");
  h2xzframe->SetMarkerStyle(1);
  h2xzframe->SetMarkerColor(3);

  h2rzframe=new TH2F("rzframe","r-z",
		     400,-20,20,100,0,10);

  h2rzframe->SetXTitle("Z (mm)  ");
  h2rzframe->SetYTitle("r (mm)  ");
  h2rzframe->SetMarkerStyle(1);
  h2rzframe->SetMarkerColor(3);

  for(int i=0;i<400;i++)
    {
      double z=-20.+i*40./400;
      h2yzframe->Fill(z,-2);
      h2yzframe->Fill(z,-3);
      h2yzframe->Fill(z,-7);
      h2yzframe->Fill(z,-8);
      h2yzframe->Fill(z,2);
      h2yzframe->Fill(z,3);
      h2yzframe->Fill(z,7);
      h2yzframe->Fill(z,8);

      h2xzframe->Fill(z,-2);
      h2xzframe->Fill(z,-3);
      h2xzframe->Fill(z,-7);
      h2xzframe->Fill(z,-8);
      h2xzframe->Fill(z,2);
      h2xzframe->Fill(z,3);
      h2xzframe->Fill(z,7);
      h2xzframe->Fill(z,8);

      h2rzframe->Fill(z,2);
      h2rzframe->Fill(z,3);
      h2rzframe->Fill(z,7);
      h2rzframe->Fill(z,8);
    }

}

double GetVariable(const char *var,TCut cut)
{
  TTree *t=(TTree*) gROOT->FindObject("t");

  TH1F *htemp = (TH1F*) gROOT->FindObject("htemp");
  if(htemp) delete htemp;
  t->Draw(Form("%s>>htemp",var),cut);
  htemp = (TH1F*) gROOT->FindObject("htemp");
  return double(htemp->GetMean());
}

int GetEntries(const char *var,TCut cut)
{
  TTree *t=(TTree*) gROOT->FindObject("t");

  TH1F *htemp = (TH1F*) gROOT->FindObject("htemp");
  if(htemp) delete htemp;
  t->Draw(Form("%s>>htemp",var),cut);
  htemp = (TH1F*) gROOT->FindObject("htemp");
  //cout<<"Number of entries = "<<htemp->GetEntries()<<endl;
  return htemp->GetEntries();
}

//display the trajectory of an event
//if entry>0, display the specified event
//otherwise, display the current event
int  DisplayEvent( int entry=0, TCut extracut="")
{
  if(thisevent<=0) Init();

  TTree *t=(TTree*) gROOT->FindObject("t");
  t->SetAlias("step_s","sqrt(step_x*step_x+step_y*step_y)");
  t->SetAlias("step_s_exp","sqrt(step_x_exp*step_x_exp+step_y_exp*step_y_exp)");
  t->SetAlias("step_s_fil","sqrt(step_x_fil*step_x_fil+step_y_fil*step_y_fil)");


  int iEvent = (entry>0) ? entry : thisevent++;
  if(iEvent>=t->GetEntries()) {
    cout<<"Reach the end of root file ...\n";
    return -1;
  }

  char strcut[255]; sprintf(strcut,"index==%d",iEvent);
  char strAcccut[255]; sprintf(strAcccut,"index==%d && step_status==1",iEvent);
  char strDiscut[255]; sprintf(strDiscut,"index==%d && step_status==0",iEvent);

  TCut cut = strcut;
  TCut Acccut = strAcccut;
  TCut Discardcut = strDiscut;

  cut += extracut;
  Acccut += extracut;
  Discardcut += extracut;

  c1->Clear();
  c1->Divide(2,3);
  c1->cd(1);
  double Pt = GetVariable("pt0",cut);
  c1->cd(2);
  double Theta0 = GetVariable("th0*57.3",cut); 
  c1->cd(3);
  double Phi0 = GetVariable("ph0*57.3",cut);
  c1->cd(4);
  int NHits = GetVariable("npt",cut);
  c1->cd(5);
  double Z0 = GetVariable("z0",cut);

  c1->Clear();
  c1->Divide(2,3);
  c1->cd(1);
  double Pt_rec = GetVariable("pt_rec",cut);
  c1->cd(2);
  double Theta_rec = GetVariable("th_rec*57.3",cut); 
  c1->cd(3);
  double Phi_rec = GetVariable("ph_rec*57.3",cut);
  c1->cd(4);
  double Z_rec = GetVariable("z_rec",cut);

  if(NHits<1) return 0;

  TPaveText *pt = new TPaveText(0.65,0.65,1-gStyle->GetPadRightMargin(),0.89,"brNDC");
  pt->SetFillColor(0);
  if(gROOT->IsBatch()) pt->SetFillStyle(4000);
  pt->SetBorderSize(0);
  //pt->AddText(Form("Event %d", iEvent));
  pt->AddText(Form("Pt_{0}=%.4f",Pt));
  pt->AddText(Form("#theta_{0}=%.1f^{o}",Theta0));
  pt->AddText(Form("#phi_{0}=%.1f^{o}",Phi0));
  pt->AddText(Form("Z_{0}=%.2f",Z0));
  //pt->AddText(Form("NHits=%d",NHits));

  TPaveText *pt_rec = new TPaveText(0.65,0.65,1-gStyle->GetPadRightMargin(),0.89,"brNDC");
  pt_rec->SetFillColor(0);
  if(gROOT->IsBatch()) pt_rec->SetFillStyle(4000);
  pt_rec->SetBorderSize(0);
  //pt_rec->AddText(Form("Event %d", iEvent));
  pt_rec->AddText(Form("Pt_{rec}=%.4f",Pt_rec));
  pt_rec->AddText(Form("#theta_{rec}=%.1f^{o}",Theta_rec));
  pt_rec->AddText(Form("#phi_{rec}=%.1f^{o}",Phi_rec));
  pt_rec->AddText(Form("Z_{rec}=%.2f",Z_rec));

  h2yxframe->SetTitle(Form("Event %d: y-x, %d Hits", iEvent,NHits));

  TObject *obj=0;
  if(obj=gROOT->FindObject("gryx")) delete obj;
  if(obj=gROOT->FindObject("gryz")) delete obj;
  if(obj=gROOT->FindObject("grxz")) delete obj;
  if(obj=gROOT->FindObject("grrz")) delete obj;

  c1->Clear();
  //Get the TGraph object, in one canvas there can be only one 
  c1->cd(); h2yxframe->Draw();
  if(!gROOT->IsBatch()) h2yxframe->Draw("same");
  t->Draw("step_y:step_x",Acccut,"same");
  TGraph *gryx=(TGraph*)(gROOT->FindObject("Graph")->Clone("gryx"));
  gryx->SetName("gryx");
  gryx->SetMarkerColor(4);gryx->SetMarkerStyle(4);

  c1->Clear();
  c1->cd(); h2yxframe->Draw();
  t->Draw("step_y_exp:step_x_exp",Acccut,"same");
  TGraph *gryx_exp=(TGraph*)(gROOT->FindObject("Graph")->Clone("gryx_exp"));
  gryx_exp->SetName("gryx_exp");
  gryx_exp->SetMarkerColor(6);gryx_exp->SetMarkerStyle(20);

  c1->Clear();
  c1->cd(); h2yxframe->Draw();
  t->Draw("step_y:step_x",Discardcut,"same");
  TGraph *gryx_discard=0;
  if(gROOT->FindObject("Graph")) {
    TGraph *gryx_discard=(TGraph*)(gROOT->FindObject("Graph")->Clone("gryx_discard"));
    gryx_discard->SetName("gryx_discard");
    gryx_discard->SetMarkerColor(2);gryx_discard->SetMarkerStyle(2);
  }

  c1->Clear();
  c1->cd(); h2yzframe->Draw();
  t->Draw("step_y:step_z",Acccut,"same");
  TGraph *gryz=(TGraph*)(gROOT->FindObject("Graph")->Clone("gryz"));
  gryz->SetName("gryz");
  gryz->SetMarkerColor(4);gryz->SetMarkerStyle(4);

  c1->cd(); h2xzframe->Draw();
  t->Draw("step_x:step_z",Acccut,"same");
  TGraph *grxz=(TGraph*)(gROOT->FindObject("Graph")->Clone("grxz"));
  grxz->SetName("grxz");
  grxz->SetMarkerColor(4);grxz->SetMarkerStyle(4);
		
  c1->cd(); h2rzframe->Draw();
  t->Draw("step_s:step_z",Acccut,"same");
  TGraph *grrz=(TGraph*)(gROOT->FindObject("Graph")->Clone("grrz"));
  grrz->SetName("grrz");
  grrz->SetMarkerColor(4);grrz->SetMarkerStyle(4);



  c1->Clear();
  c1->Divide(2,2);

  c1->cd(1); h2yxframe->Draw();
  pt->Draw();
  if(!gROOT->IsBatch()) h2yxframe->Draw("same");
  gryx->Draw("same p");
  if(gryx_discard) gryx_discard->Draw("same p");
  t->Draw("step_y_fil:step_x_fil>>h2yx_fil",Acccut,"same*");

  c1->cd(2); h2yzframe->Draw();
  pt_rec->Draw();
  gryz->Draw("same p");
  t->Draw("step_y_fil:step_z_fil>>h2yz_fil",Acccut,"same*");

  c1->cd(3); h2xzframe->Draw();
  grxz->Draw("same p");
  t->Draw("step_x_fil:step_z_fil>>h2xz_fil",Acccut,"same*");

		
  c1->cd(4); h2rzframe->Draw();
  grrz->Draw("same p");
  t->Draw("step_s_fil:step_z_fil>>h2rz_fil",Acccut,"same*");

  c1->Update();
  c1->SaveAs(Form("Movie/RTPC_Event_%03d.png",iEvent));

  //////////////plot zoom in figure////////////
  c2->Clear();
  c2->cd();
  gryx->SetTitle(Form("Event %d: Raw(Blue),Expected(Purple), Filtered(Black)",iEvent));
  gryx->Draw("AP");
  gryx_exp->Draw("samep");
  if(gryx_discard) gryx_discard->Draw("samep");
  t->Draw("step_y_fil:step_x_fil",Acccut,"same*");
  c2->Update();
  c2->SaveAs(Form("Movie/RTPC_Event_Zoom_%03d.png",iEvent));


  t->Scan("p0:pt0:th0*57.3:ph0*57.3:z0:npt:pt_rec:th_rec*57.3:ph_rec*57.3:z_rec",cut);

  return NHits;
}

//slice show n events, each last for specified second
//if istart<=0, start from current events, otherwise from specified event
void go(int n=1,int second=2, int istart=0)
{ 
  if(istart>0)  thisevent=istart;
  for(int i=0;i<n;i++)
    {
      int nHits = 0;
      while (nHits==0) 
	{
	  nHits = DisplayEvent(thisevent++);
	  if(nHits == -1) break;
	}
      c1->Update(); 
      if(nHits == -1) break;
      gSystem->Sleep(second*1000);
    }

}

void DisplayMyEvent(TCut cut="")
{ 
  int nHits = 0;
  while (nHits==0) 
    {
      nHits = DisplayEvent(0,cut);
      if(nHits == -1) break;
    }
  c1->Update();
}

//slice show n events, each last for specified second
//if istart<=0, start from current events, otherwise from specified event
void SliceShowMyEvent(int n=1, int second=2, int istart=0, 
		      TCut cut="")
{ 
  if(istart>0)  thisevent=istart;
  for(int i=0;i<n;i++)
    {
      int nHits = 0;
      while (nHits==0) 
	{
	  nHits = DisplayEvent(thisevent++,cut);
	  if(nHits == -1) break;
	}
      c1->Update();
      if(nHits == -1) break;
      gSystem->Sleep(second*1000);

    }
}
