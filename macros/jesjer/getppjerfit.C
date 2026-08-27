#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

void getppjerfit()
{

  SetsPhenixStyle();
  TCanvas *c = new TCanvas();
  //SetsPhenixStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TH1F *h = new TH1F("h","",5,15,60);
  h->GetYaxis()->SetRangeUser(0.1,0.5);
  h->GetXaxis()->SetTitle("p_{T}^{truth} [GeV]");
  h->GetYaxis()->SetTitle("#sigma(p_{T}^{reco}/p_{T}^{truth})/<p_{T}^{reco}/p_{T}^{truth}>");
  TFile *f = new TFile("ppjesjer.root","READ");
  TGraphErrors *g_pp = (TGraphErrors*)f->Get("h_jer_run21_r04_z0_eta0123");
  TF1 *f_fitjer = new TF1("f_fitjer", "sqrt( [0]*[0] + ([1]*[1]/x) + ([2]*[2]/(x*x)) )", 0, 100);
  f_fitjer->SetParNames("Constant term", "Stochastic term", "Noise term");
  f_fitjer->SetLineColor(kRed);

  float p_a[5];
  float p_b[5];
  float p_c[5];
  float e_a[5];
  float e_b[5];
  float e_c[5];
  
  f_fitjer->SetParameters(0.1, 0.6, 2);
  double fit_min, fit_max;
  fit_min = 15;
  fit_max = 60;
  g_pp->Fit(f_fitjer, "", "", fit_min, fit_max);
  h->Draw();
  g_pp->Draw("PSAME");
  p_a[0] = f_fitjer->GetParameter(1);
  p_b[0] = f_fitjer->GetParameter(2);
  p_c[0] = f_fitjer->GetParameter(0);
  e_a[0] = f_fitjer->GetParError(1);
  e_b[0] = f_fitjer->GetParError(2);
  e_c[0] = f_fitjer->GetParError(0);

  TLegend *leg = new TLegend(.15,.7,.45,.9);
  leg->SetFillStyle(0);
  leg->AddEntry("","#it{#bf{sPHENIX}} Internal","");
  leg->AddEntry("","p+p #sqrt{s_{NN}}=200 GeV","");
  leg->AddEntry("","anti-#it{k}_{#it{t}} #it{R} = 0.3","");
  leg->AddEntry("","Pythia8","");
  
  TLegend *leg2 = new TLegend(.5,.4,.8,.9);
  leg2->SetFillStyle(0);
  leg2->AddEntry("","Fit: #frac{a}{#sqrt{p_{T}^{truth}}} #oplus #frac{b}{p_{T}^{truth}} #oplus c","");
  leg2->AddEntry("",Form("a = %1.3f #pm %1.3f",p_a[0],e_a[0]),"");
  leg2->AddEntry("",Form("b = %1.3f #pm %1.3f",p_b[0],e_b[0]),"");
  leg2->AddEntry("",Form("c = %1.3f #pm %1.3f",p_c[0],e_c[0]),"");
  
  leg->Draw("SAME");
  leg2->Draw("SAME");
  
  c->Print("jerfitpp.pdf");
  c->Clear();
  
  TFile *f2 = new TFile("hists_highcut_v14.root","READ");
  TGraphErrors *g_auau[5];
  string centbins[] = {"0-100%","0-10%","10-30%","30-50%","50-90%"};
  for(int i = 0; i < 5; i++)
    {
      TLegend *leg3 = new TLegend(.15,.7,.45,.9);
      leg3->SetFillStyle(0);
      leg3->AddEntry("","#it{#bf{sPHENIX}} Internal","");
      leg3->AddEntry("","Au+Au #sqrt{s_{NN}}=200 GeV","");
      leg3->AddEntry("","anti-#it{k}_{#it{t}} #it{R} = 0.3","");
      leg3->AddEntry("","Pythia8 + Hijing","");
      leg3->AddEntry("",centbins[i].c_str(),"");


      h->Draw();
      f_fitjer->SetParameters(0.1, 0.6, 5);
      f_fitjer->FixParameter(1,p_a[0]);
      f_fitjer->FixParameter(0,p_c[0]);
      g_auau[i] = (TGraphErrors*) f2->Get(Form("g_jer_cent%i",i));
      g_auau[i]->SetMarkerColor(kBlack);
      g_auau[i]->SetMarkerStyle(20);
      g_auau[i]->Fit(f_fitjer, "", "", fit_min, fit_max);
      h->Draw();
      g_auau[i]->Draw("PSAME");

      p_a[i+1] = f_fitjer->GetParameter(1);
      p_b[i+1] = f_fitjer->GetParameter(2);
      p_c[i+1] = f_fitjer->GetParameter(0);
      e_a[i+1] = f_fitjer->GetParError(1);
      e_b[i+1] = f_fitjer->GetParError(2);
      e_c[i+1] = f_fitjer->GetParError(0);
      
      TLegend *leg4 = new TLegend(.5,.4,.8,.9);
      leg4->SetFillStyle(0);
      leg4->AddEntry("","Fit: #frac{a}{#sqrt{p_{T}^{truth}}} #oplus #frac{b}{p_{T}^{truth}} #oplus c","");
      leg4->AddEntry("",Form("a = %1.3f #pm %1.3f",p_a[i+1],e_a[i+1]),"");
      leg4->AddEntry("",Form("b = %1.3f #pm %1.3f",p_b[i+1],e_b[i+1]),"");
      leg4->AddEntry("",Form("c = %1.3f #pm %1.3f",p_c[i+1],e_c[i+1]),"");

      leg3->Draw("SAME");
      leg4->Draw("SAME");
      
      c->Print(Form("jerfit%i_v14.pdf",i));
      c->Clear();
    }
  //c->Print("jerfits.pdf)");




   
   TH1F *h_par[3];
   TH1F *h_sys = new TH1F("h_sys","",6,0,6);
   h_sys->SetMarkerColor(kRed);
   h_sys->SetLineColor(kRed);
   float sysweight[] = {0.,0.,0.996711,1.05849,1.12146,0.974294};
   for(int p = 0; p < 3; p++)
     {
       h_par[p] = new TH1F(Form("h_par%i",p),"",6,0,6);
       h_par[p]->GetXaxis()->SetBinLabel(1,"pp");
       h_par[p]->GetXaxis()->SetBinLabel(2,"0-100\%");
       h_par[p]->GetXaxis()->SetBinLabel(3,"0-10\%");
       h_par[p]->GetXaxis()->SetBinLabel(4,"10-30\%");
       h_par[p]->GetXaxis()->SetBinLabel(5,"30-50\%");
       h_par[p]->GetXaxis()->SetBinLabel(6,"50-90\%");
       for(int i = 0; i < 6; i++)
	 {
	   if(p == 0)
	     {
	       h_par[p]->SetBinContent(i+1,p_a[i]);
	       h_par[p]->SetBinError(i+1,e_a[i]);
	       h_par[p]->GetYaxis()->SetTitle("a");
	     }
	   else if(p == 1)
	     {
	       h_par[p]->SetBinContent(i+1,p_b[i]);
               h_par[p]->SetBinError(i+1,e_b[i]);
	       h_par[p]->GetYaxis()->SetTitle("b");
	       
	       h_sys->SetBinContent(i+1,p_b[i]*sysweight[i]);
	       h_sys->SetBinError(i+1,e_b[i]*sysweight[i]);
             }
	   else	if(p == 2)
             {
               h_par[p]->SetBinContent(i+1,p_c[i]);
               h_par[p]->SetBinError(i+1,e_c[i]);
	       h_par[p]->GetYaxis()->SetTitle("c");
             }
	 }
       h_par[p]->Draw();
       if(p == 1)
	 {
	   h_sys->Draw("SAME");
	   TLegend *legb = new TLegend(.18,.7,.33,.85);
	   legb->SetFillStyle(0);
	   legb->AddEntry("","#it{#bf{sPHENIX}} Internal","");
	   //legb->AddEntry("","p+p #sqrt{s_{NN}}=200 GeV","");
	   legb->AddEntry(h_par[p],"b_{MC}","p");
	   legb->AddEntry(h_sys,"b_{Data}","p");
	   legb->Draw("SAME");
	 }
       c->Print(Form("par%i_v14.pdf",p));
     }

   TF1 *f_jer[6];
   TF1 *f_jersys[6];
   TF1 *f_jersysup[6];
   TF1 *f_jersysdown[6];
 
   for(int i = 0; i < 6; i++)
     {
       f_jer[i] = new TF1(Form("f_jer%i",i), "sqrt( [0]*[0] + ([1]*[1]/x) + ([2]*[2]/(x*x)) )", 0, 100);
       f_jer[i]->FixParameter(1,p_a[0]);
       f_jer[i]->FixParameter(2,p_b[i]);
       f_jer[i]->FixParameter(0,p_c[0]);
       f_jersys[i] = new TF1(Form("f_jersys%i",i), "sqrt( [0]*[0] + ([1]*[1]/x) + ([2]*[2]/(x*x)) )", 0, 100);
       f_jersys[i]->FixParameter(1,p_a[0]);
       f_jersys[i]->FixParameter(2,p_b[i]*sysweight[i]);
       f_jersys[i]->FixParameter(0,p_c[0]);

       f_jersysup[i] = new TF1(Form("f_jersysup%i",i), "sqrt( [0]*[0] + ([1]*[1]/x) + ([2]*[2]/(x*x)) )", 0, 100);
       f_jersysup[i]->FixParameter(1,p_a[0]);
       f_jersysup[i]->FixParameter(2,(p_b[i]+e_b[i])*sysweight[i]);
       f_jersysup[i]->FixParameter(0,p_c[0]);

       f_jersysdown[i] = new TF1(Form("f_jersysdown%i",i), "sqrt( [0]*[0] + ([1]*[1]/x) + ([2]*[2]/(x*x)) )", 0, 100);
       f_jersysdown[i]->FixParameter(1,p_a[0]);
       f_jersysdown[i]->FixParameter(2,(p_b[i]-e_b[i])*sysweight[i]);
       f_jersysdown[i]->FixParameter(0,p_c[0]);

     }
   TFile *f_out = new TFile("jerfits_v14.root","RECREATE");
   f_jer[2]->Write("f_jer_mc_0_10");
   f_jer[3]->Write("f_jer_mc_10_30");
   f_jer[4]->Write("f_jer_mc_30_50");
   f_jer[5]->Write("f_jer_mc_50_90");
   f_jersys[2]->Write("f_jer_data_0_10");
   f_jersys[3]->Write("f_jer_data_10_30");
   f_jersys[4]->Write("f_jer_data_30_50");
   f_jersys[5]->Write("f_jer_data_50_90");
   f_jersysup[2]->Write("f_jer_sysup_0_10");
   f_jersysup[3]->Write("f_jer_sysup_10_30");
   f_jersysup[4]->Write("f_jer_sysup_30_50");
   f_jersysup[5]->Write("f_jer_sysup_50_90");
   f_jersysdown[2]->Write("f_jer_sysdown_0_10");
   f_jersysdown[3]->Write("f_jer_sysdown_10_30");
   f_jersysdown[4]->Write("f_jer_sysdown_30_50");
   f_jersysdown[5]->Write("f_jer_sysdown_50_90");
   
   //c->Print("jerfunctions.pdf(");
   for(int i = 2; i < 6; i++)
     {
       h->GetXaxis()->SetTitle("p_{T}^{truth} [GeV]");
       h->GetYaxis()->SetTitle("JER fit");
       h->Draw();
       f_jer[i]->Draw("SAME");
       f_jersys[i]->SetLineColor(kRed);
       f_jersys[i]->Draw("SAME");
       f_jersysup[i]->SetLineColor(kBlue);
       f_jersysup[i]->Draw("SAME");
       f_jersysdown[i]->SetLineColor(kGreen+3);
       f_jersysdown[i]->Draw("SAME");


       TLegend *legg = new TLegend(.5,.6,.85,.9);
       legg->SetFillStyle(0);
       legg->AddEntry(f_jer[i],"Pythia+Hijing","l");
       legg->AddEntry(f_jersys[i],"b scaled to data","l");
       legg->AddEntry(f_jersysup[i],"b scaled to data + uncert.","l");
       legg->AddEntry(f_jersysdown[i],"b scaled to data - uncert.","l");
       legg->AddEntry("",centbins[i-1].c_str(),"");
       legg->Draw("SAME");
       
       c->Print(Form("jerfunctions%i_v14.pdf",i));
     }
   //c->Print("jerfunctions.pdf)");
}
