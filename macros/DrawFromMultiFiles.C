//this script is used to draw a histo using multiple root files then compare them
//user need to provide target string, binning string, cut string

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
#include "TChain.h"
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

// key is used to tell what is the feature for all input files 
// hname# is used to tell what is the special feature of file# among all input files 
void DrawFromMultiFiles(const char* tgStr, const char *binStr,
			const char *cutStr="", const char* treename="t", const char *key="",
			const char *file0=0, int color0=1, const char *hname0="0", 
			const char *file1=0, int color1=2, const char *hname1="1", 
			const char *file2=0, int color2=4, const char *hname2="2", 
			const char *file3=0, int color3=6, const char *hname3="3", 
			const char *file4=0, int color4=8, const char *hname4="4")
{

  system("mkdir -p Graph");
  ofstream fout;
  fout.open("table.txt",ios_base::app);

  const int NFILE = 5;
  TChain *tree[NFILE]; 
  const char *hname[NFILE]={hname0,hname1,hname2,hname3,hname4};
  const char *fname[NFILE]={file0,file1,file2,file3,file4};
  int color[NFILE]={color0,color1,color2,color3,color4};
  TH1F *h[NFILE];
  double mean[NFILE],sigma[NFILE];
  
  int n=0;  //number of valid files
  for(int i=0;i<NFILE;i++) {
    if(fname[i]) {
      tree[i]=new TChain(treename); 
      tree[i]->Add(fname[i]); 
      n=i+1;
    }
  }
  

  //////////////////////////////////////////////////////////////
  TCanvas *c12 = new TCanvas("c12","",900,700);
  c12->cd();
  //create the histo
  for(int i=0;i<n;i++) {
    char tg[200];
    sprintf(tg,"%s>>%s%s",tgStr,hname[i],binStr);
    tree[i]->Draw(tg,cutStr);
    h[i] = (TH1F*) gROOT->FindObject(hname[i]);
    h[i]->SetLineColor(color[i]);
    h[i]->SetTitle(hname[i]);
  }
  
  //now draw all histo into the same pad
 
  
  //now determine which histo has the GetMaximum height then use it as frame
  
  int pMaxIndex=0;
  double max= h[0]->GetMaximum(); 
  for(int i=1;i<n;i++) {
	if(max < h[i]->GetMaximum()) {max=h[i]->GetMaximum(); pMaxIndex=i;}
  }
	
  TH1F *hFrame=(TH1F*)h[pMaxIndex]->Clone("hFrame");
  hFrame->SetTitle(key);            //change the title
  hFrame->Draw();
  //now draw the legend so it is on the background
  bool drawLegend=true;
  if(drawLegend){

	double xx=gStyle->GetPadLeftMargin()+0.03; 
	double yy=1.0-gStyle->GetPadTopMargin();

	TLegend *leg = new TLegend(xx,yy-0.06*n,0.52,yy,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextFont(62);
	leg->SetFillColor(0);
	leg->SetFillStyle(4000);
	
	for(int i=0;i<n;i++) {
	TLegendEntry *entry=leg->AddEntry(hname[i],hname[i],"lf");
	entry->SetLineWidth(3);
	entry->SetLineColor(color[i]);
	entry->SetTextColor(color[i]);
	}
	leg->Draw();
  }
  //now draw the histo so they will all above the legend
  for(int i=0;i<n;i++) h[i]->Draw("same");  
    
  c12->SaveAs(Form("Graph/Cmp_%s.png",key));
 
  
  //////////////////////////////////////////////////////////////
  int height=900, width=700;
  if (n==2 || n==3) width=500;
  TCanvas *c13 = new TCanvas("c13","",width,height);
  
  if(n>=5) c13->Divide(2,3);
  else if(n==4) c13->Divide(2,2);
  else if(n>=2) c13->Divide(1,n);
  
  for(int i=0;i<n;i++) {
    c13->cd(i+1);    
    FitGauss(h[i],mean[i],sigma[i]);
    //h[i]->Draw();   
  }
  c13->SaveAs(Form("Graph/Cmp_MultiPad_%s.png",key));

  
  bool log=true;
  if(log) {
	fout<<"\n"<<key<<endl;
    for(int i=0;i<n;i++) {
      fout<<fixed<<setprecision(4)
	  <<setw(8)<<mean[i]<<" +/- " <<setw(5)<<sigma[i]<<"  "
	  <<setw(30)<<hname[i]<<endl;
    }
  }
  fout.close();

}
