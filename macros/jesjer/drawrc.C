#include "sPhenixStyle.h"
#include "sPhenixStyle.C"
void drawrc()
{
  SetsPhenixStyle();
  TCanvas *c = new TCanvas();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  
  TFile *f_mc = new TFile("rc_sub1_r03_MB.root","READ");
  TFile *f_data = new TFile("rc_sub1_r03.root","READ");

  TGraphErrors *g_sigma_mc = (TGraphErrors*) f_mc->Get("rc_sub1_sigma_cent");
  g_sigma_mc->SetName("g_sigma_mc");
  TGraphErrors *g_sigma_data = (TGraphErrors*) f_data->Get("rc_sub1_sigma_cent");
  g_sigma_data->SetName("g_sigma_data");
  TGraphErrors *g_mean_mc = (TGraphErrors*) f_mc->Get("rc_sub1_mean_cent");
  g_mean_mc->SetName("g_mean_mc");
  TGraphErrors *g_mean_data = (TGraphErrors*) f_data->Get("rc_sub1_mean_cent");
  g_mean_data->SetName("g_mean_data");

  g_sigma_mc->Draw();
  g_sigma_data->Draw("SAME");

  int ncent = 4; //0-10,10-30,30-50, 50-90
  TH1D *h_mc[ncent];
  TH1D *h_data[ncent];

  h_mc[0] = (TH1D*) f_mc->Get("rc_sub1_cent_1d_cent0_10");
  h_mc[0]->SetName("h_mc0");
  
  h_data[0] = (TH1D*) f_data->Get("rc_sub1_cent_1d_cent0_5");
  h_data[0]->SetName("h_data0");
  h_data[0]->Add((TH1D*) f_data->Get("rc_sub1_cent_1d_cent5_10"));

  h_mc[1] = (TH1D*) f_mc->Get("rc_sub1_cent_1d_cent10_20");
  h_mc[1]->SetName("h_mc1");
  h_mc[1]->Add((TH1D*) f_mc->Get("rc_sub1_cent_1d_cent20_30"));

  h_data[1] = (TH1D*) f_data->Get("rc_sub1_cent_1d_cent10_20");
  h_data[1]->SetName("h_data1");
  h_data[1]->Add((TH1D*) f_data->Get("rc_sub1_cent_1d_cent20_30"));

  h_mc[2] = (TH1D*) f_mc->Get("rc_sub1_cent_1d_cent30_40");
  h_mc[2]->SetName("h_mc2");
  h_mc[2]->Add((TH1D*) f_mc->Get("rc_sub1_cent_1d_cent40_50"));

  h_data[2] = (TH1D*) f_data->Get("rc_sub1_cent_1d_cent30_40");
  h_data[2]->SetName("h_data2");
  h_data[2]->Add((TH1D*) f_data->Get("rc_sub1_cent_1d_cent40_50"));

  h_mc[3] = (TH1D*) f_mc->Get("rc_sub1_cent_1d_cent50_60");
  h_mc[3]->SetName("h_mc3");
  h_mc[3]->Add((TH1D*) f_mc->Get("rc_sub1_cent_1d_cent60_70"));
  h_mc[3]->Add((TH1D*) f_mc->Get("rc_sub1_cent_1d_cent70_80"));
  h_mc[3]->Add((TH1D*) f_mc->Get("rc_sub1_cent_1d_cent80_90"));

  h_data[3] = (TH1D*) f_data->Get("rc_sub1_cent_1d_cent50_60");
  h_data[3]->SetName("h_data3");
  h_data[3]->Add((TH1D*) f_data->Get("rc_sub1_cent_1d_cent60_70"));
  h_data[3]->Add((TH1D*) f_data->Get("rc_sub1_cent_1d_cent70_80"));
  h_data[3]->Add((TH1D*) f_data->Get("rc_sub1_cent_1d_cent80_90"));
  
  //c->Print("fluctuations.pdf(");
  string centbins[] = {"0-10%","10-30%","30-50%","50-90%"};  
  for(int i = 0; i < 4; i++)
    {
      h_mc[i]->Scale(1./h_mc[i]->Integral());
      h_data[i]->Scale(1./h_data[i]->Integral());
      h_mc[i]->GetYaxis()->SetRangeUser(0,0.6);
      h_mc[i]->Draw();
      h_data[i]->SetMarkerColor(kRed);
      h_data[i]->SetLineColor(kRed);
      h_data[i]->Draw("SAME");
      
      TLegend *leg = new TLegend(.12,.7,.42,.9);
      leg->SetFillStyle(0);
      leg->AddEntry("","#it{#bf{sPHENIX}} Internal","");
      leg->AddEntry("","Au+Au #sqrt{s_{NN}}=200 GeV","");
      leg->AddEntry("","#it{R} = 0.3 Random Cones","");
      leg->AddEntry("",centbins[i].c_str(),"");
      
      TLegend *leg2 = new TLegend(.5,.6,.8,.9);
      leg2->SetFillStyle(0);
      leg2->AddEntry(h_mc[i],"MB Hijing","p");
      leg2->AddEntry("",Form("Mean = %1.3f #pm %1.3f",h_mc[i]->GetMean(),h_mc[i]->GetMeanError()),"");
      leg2->AddEntry("",Form("Sigma = %1.3f #pm %1.3f",h_mc[i]->GetRMS(),h_mc[i]->GetRMSError()),"");
      leg2->AddEntry(h_data[i],"MB Data","p");
      leg2->AddEntry("",Form("Mean = %1.3f #pm %1.3f",h_data[i]->GetMean(),h_data[i]->GetMeanError()),"");
      leg2->AddEntry("",Form("Sigma = %1.3f #pm %1.3f",h_data[i]->GetRMS(),h_data[i]->GetRMSError()),"");

      std::cout<<centbins[i]<<" "<<h_data[i]->GetRMS()/h_mc[i]->GetRMS()<<std::endl;
      
      leg->Draw("SAME");
      leg2->Draw("SAME");
      c->Print(Form("fluctuations%i.pdf",i));
    }
  //c->Print("fluctuations.pdf)");
}
