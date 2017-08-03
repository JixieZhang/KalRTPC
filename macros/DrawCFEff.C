void DrawCFEff(double p1=1.9,double p2=30,double p3=40,double p4=0.4)
{
  int n = t->Draw("CF_ChainNum>>h11(4,-0.5,3.5)","CF_ChainNum<3.5","");
  t->Draw("CF_ChainNum>>h12(4,-0.5,3.5)","CF_ChainNum<3.5","");

  h11->Scale(100./n);
  h12->Scale(100./n);
  h11->SetTitle("ChainFinder Efficiency; Number of Found Chain ; Efficiency(%) ");
  h12->SetTitle("ChainFinder Efficiency; Number of Found Chain ; Efficiency(%) ");
  h11->SetMaximum(100);
  h11->Draw();
  h12->Draw("textsame");
  
  TLatex *tex = new TLatex(1.6,75.0,Form("para = {%.1f, %.1f, %.1f, %.1f}",p1,p2,p3,p4));
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  
  c1->SaveAs("CF_eff_shift.png");
  
  
}


