//this script is used to check kalman filter eficiency
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
#include "TLegend.h"
#include "TLegendEntry.h"


TF1* FitGauss(TH1* h1, double &mean, double &sigma, double range_in_sigma=1.5)
{
	//cout<<"h1->GetEntries()="<<h1->GetEntries()<<endl;
	if(h1->GetEntries()<50) return NULL;

	double xmin=h1->GetMean()-5.0*h1->GetRMS();
	double xmax=h1->GetMean()+5.0*h1->GetRMS();
	//h1->Fit("gaus","RQ","",xmin,xmax);
	h1->Fit("gaus","Q");
	TF1 *f=(TF1*)h1->GetListOfFunctions()->FindObject("gaus");
	mean=f->GetParameter(1);
	sigma=f->GetParameter(2);
	xmin=mean-range_in_sigma*sigma;
	xmax=mean+range_in_sigma*sigma;
	double chi2_1 = f->GetChisquare()/f->GetNDF();

	h1->Fit("gaus","RQ","",xmin,xmax);
	f=(TF1*)h1->GetListOfFunctions()->FindObject("gaus");	
	double chi2_2 = f->GetChisquare()/f->GetNDF();
	//cout<<"chi2_1="<<chi2_1<<"\t chi2_2="<<chi2_2<<endl;
	if(chi2_2>1.1*chi2_1) {
      //2nd iteration is worse than 1st iteratin, roll back to previous fit
	  h1->Fit("gaus","Q");
	  f=(TF1*)h1->GetListOfFunctions()->FindObject("gaus");	
	}
		
	mean=f->GetParameter(1);
	sigma=f->GetParameter(2);
	  

	if(gStyle->GetOptFit()==0)
	{
		char str[100];
		TText *text=0;

		double xx=gStyle->GetPadLeftMargin()+0.03; 
		TPaveText *pt = new TPaveText(xx,0.20,xx+0.45,0.45,"brNDC");
		pt->SetBorderSize(0);
		pt->SetFillColor(0);
		sprintf(str,"Mean = %.3G",mean);
		text=pt->AddText(str);
		text->SetTextColor(2);
		sprintf(str,"Sigma = %.3G",sigma);
		text=pt->AddText(str);
		text->SetTextColor(2);
		pt->Draw("same");
	}

	return f;
}

TF1* FitGauss(TH1* h1,double range_in_sigma=1.5)
{
	double mean,sigma;
	return FitGauss(h1,mean,sigma,range_in_sigma);
}


void CheckEff(int range=50, const char* cut="",const char* key="")
{
  //////////////////////////////////////////////////////////////
  //only show the entries in the opt pave text
	//gStyle->SetOptStat(10);	
	//gStyle->SetOptFit(0);
	//gStyle->SetStatH(0.12);  //If OptFit is 1 than this number will not work

	TCut cut_dPt_50 = "abs(pt0-pt_rec-0.0087)<0.01";
	TCut cut_dTh_50 = "abs(th0-th_rec)<0.075";
	TCut cut_dPh_50 = "abs(ph0-ph_rec-0.012)<0.06";

	TCut cut_dPt_70 = "abs(pt0-pt_rec-0.0027)<0.021";
	TCut cut_dTh_70 = "abs(th0-th_rec)<0.060";
	TCut cut_dPh_70 = "abs(ph0-ph_rec-0.001)<0.045";
/*
TFile *_file0 = TFile::Open("nt_50to70MeV_All_2Iter_Backward.root")
t->Draw("pt0-pt_rec","abs(pt0-pt_rec-0.0087)<0.01","")
t->Draw("pt0-pt_rec","abs(pt0-pt_rec-0.0087)<0.01 && abs(th0-th_rec)<0.075","same")
t->Draw("pt0-pt_rec","abs(pt0-pt_rec-0.0087)<0.01 && abs(th0-th_rec)<0.075 && abs(ph0-ph_rec-0.012)<0.06","same")
t->Draw("pt0-pt_rec","abs(pt0-pt_rec-0.0087)<0.01 && abs(th0-th_rec)<0.075 && abs(ph0-ph_rec-0.012)<0.06","same")
.q
TFile *_file0 = TFile::Open("nt_70to250MeV_All_1Iter.root")
t->Draw("pt0-pt_rec","abs(pt0-pt_rec-0.0027)<0.021","")
t->Draw("pt0-pt_rec","abs(pt0-pt_rec-0.0027)<0.021 && abs(th0-th_rec)<0.060","same")
t->Draw("pt0-pt_rec","abs(pt0-pt_rec-0.0027)<0.021 && abs(th0-th_rec)<0.060 && abs(ph0-ph_rec-0.001)<0.045","same")
t->Draw("pt0-pt_rec","abs(pt0-pt_rec-0.0027)<0.021 && abs(th0-th_rec)<0.060 && abs(ph0-ph_rec-0.001)<0.045","same")
.q
*/

    TCut cut_user = cut;
	TCut cut_dPt, cut_dTh, cut_dPh;
	if(range==50) {
	cut_dPt=cut_dPt_50;
	cut_dTh=cut_dTh_50;
	cut_dPh=cut_dPh_50;
	}
	else {
	cut_dPt=cut_dPt_70;
	cut_dTh=cut_dTh_70;
	cut_dPh=cut_dPh_70;
	}


	TTree *t = (TTree*)gDirectory->Get("t");

	TCanvas *c11 = new TCanvas("c11","",600,450);
	double N000,N001,N011,N111;
	t->Draw("pt0>>h000",cut_user,"");
	TH1F *h000 = (TH1F*) gROOT->FindObject("h000");
	h000->SetLineColor(1);
	N000=h000->GetEntries();
	h000->SetTitle(Form("All: 100%%"));

	t->Draw("pt0>>h001",cut_user && cut_dPt,"same");
	TH1F *h001 = (TH1F*) gROOT->FindObject("h001");
	h001->SetLineColor(2);
	h001->SetMarkerColor(2);
	h001->SetMarkerStyle(2);
	N001=h001->GetEntries();
	h001->Divide(h000);
	h001->SetTitle(Form("3-#sigma dPt: %.1f%%",N001/N000*100));

	t->Draw("pt0>>h011",cut_user && cut_dPt && cut_dTh,"same");
	TH1F *h011 = (TH1F*) gROOT->FindObject("h011");
	h011->SetLineColor(3);
	h011->SetMarkerColor(3);
	h011->SetMarkerStyle(4);
	N011=h011->GetEntries();
	h011->Divide(h000);
	h011->SetTitle(Form("3-#sigma dPt & d#theta: %.1f%%",N011/N000*100));

	t->Draw("pt0>>h111",cut_user && cut_dPt && cut_dTh && cut_dPh,"same");
	TH1F *h111 = (TH1F*) gROOT->FindObject("h111");
	h111->SetLineColor(4);
	h111->SetMarkerColor(4);
	h111->SetMarkerStyle(20);
	N111=h111->GetEntries();
	h111->Divide(h000);
	h111->SetTitle(Form("3-#sigma dPt & d#theta & d#phi: %.1f%%",N111/N000*100));

	//prepare the legend
	double xx=gStyle->GetPadLeftMargin()+0.03; 
	double yy=0.3+gStyle->GetPadBottomMargin();
	TLegend *leg = new TLegend(xx,yy-0.06*3,0.58,yy,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextFont(62);
	leg->SetFillColor(0);
	leg->SetFillStyle(4000);
	TLegendEntry *entry;
	entry=leg->AddEntry(h001,"","lf");
	entry->SetLineWidth(3);
	entry->SetLineColor(h001->GetLineColor());
	entry->SetTextColor(h001->GetLineColor());
	
	entry=leg->AddEntry(h011,"","lf");
	entry->SetLineWidth(3);
	entry->SetLineColor(h011->GetLineColor());
	entry->SetTextColor(h011->GetLineColor());
	
	entry=leg->AddEntry(h111,"","lf");
	entry->SetLineWidth(3);
	entry->SetLineColor(h111->GetLineColor());
	entry->SetTextColor(h111->GetLineColor());
	
	//make the horizontal line at y=1.0
	double x1=h000->GetXaxis()->GetXmin();
	double x2=h000->GetXaxis()->GetXmax();
	TLine hline(x1,1.0,x2,1.0);
	hline.SetLineWidth(2);
	hline.SetLineColor(1);
	hline.SetLineStyle(7);
	
	
	TH1F *hFrame = (TH1F*)h000->Clone("Frame");
	hFrame->Scale(1.2);
	hFrame->Divide(h000);
	hFrame->SetTitle("; pt (GeV/c); ratio_to_thrown");
	hFrame->SetLineColor(0);
	hFrame->SetMarkerColor(0);
	hFrame->Draw();
	leg->Draw("same");
	hline.Draw("plsame");
	h001->Draw("psame");
	h011->Draw("psame");
	h111->Draw("psame");

	if(strlen(key)>2) {
		c11->SaveAs(Form("Efficiency_%s.png",key));
	}
	else {	
		if(range==50) c11->SaveAs("Efficiency_Pt50to70.png");
		else c11->SaveAs("Efficiency_Pt70to250.png");
	}
}
