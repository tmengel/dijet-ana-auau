#include "sPhenixStyle.h"
#include "sPhenixStyle.C"
void drawfractions()
{
  SetsPhenixStyle();
  TCanvas *c = new TCanvas();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  
  TFile *f_mc = new TFile("MCfractions.root","READ");
  TFile *f_data = new TFile("datafractions.root","READ");

  
  string centbins[] = {"0-10%","10-30%","30-60%","60-100%"};  

  TLegend *leg = new TLegend(.15,.7,.45,.9);
  leg->SetFillStyle(0);
  leg->AddEntry("","#it{#bf{sPHENIX}} Internal","");
  leg->AddEntry("","p+p #sqrt{s_{NN}}=200 GeV","");
  leg->AddEntry("","anti-#it{k}_{#it{t}} #it{R} = 0.3","");  

  TH3F *h_data = (TH3F*)f_data->Get("h_ET_EMfrac_reco");
  h_data->SetName("h_data");
  TH3F *h_mc = (TH3F*)f_mc->Get("h_ET_EMfrac_reco");
  h_mc->SetName("h_mc");

  int npt = h_data->GetNbinsX();
  int ncent = h_data->GetNbinsZ();
  c->Print("fractions.pdf(");
  for(int j = 0; j < ncent; j++)
    {
      for(int i = 0; i < npt-1; i++)
	{
	  if(h_data->GetXaxis()->GetBinLowEdge(i+2) < 20) continue;
	  h_data->GetXaxis()->SetRange(i+1,i+1);
	  h_mc->GetXaxis()->SetRange(i+1,i+1);
	  h_data->GetZaxis()->SetRange(j+1,j+1);
          h_mc->GetZaxis()->SetRange(j+1,j+1);
	  TH1F *hd = (TH1F*) h_data->Project3D("y");
	  hd->SetName(Form("hd%i",i));
	  if(hd->Integral() == 0) continue;
	  hd->Scale(1./hd->Integral());
	  hd->Draw();
	  TH1F *hm = (TH1F*) h_mc->Project3D("y");
	  hm->SetName(Form("hm%i",i));
	  if(hm->Integral() == 0) continue;
	  hm->Scale(1./hm->Integral());
	  hm->SetMarkerColor(kRed);
	  hm->SetLineColor(kRed);
	  hm->Draw("SAME");

	  leg->Draw("SAME");
	  TLegend *leg2 = new TLegend(.7,.6,.9,.9);
	  leg2->SetFillStyle(0);
	  leg2->AddEntry(hd,"Data","p");
	  leg2->AddEntry(hm,"MC","p");
	  leg2->AddEntry("",Form("%2.0f GeV < p_{T} < %2.0f GeV",h_data->GetXaxis()->GetBinLowEdge(i+1),h_data->GetXaxis()->GetBinLowEdge(i+2)),"");
	  leg2->AddEntry("",centbins[j].c_str(),"");
	  leg2->Draw("SAME");
	  
	  c->Print("fractions.pdf");
	}
      h_data->GetXaxis()->SetRangeUser(20,60);
      h_mc->GetXaxis()->SetRangeUser(20,60);
      h_data->GetZaxis()->SetRange(j+1,j+1);
      h_mc->GetZaxis()->SetRange(j+1,j+1);
      TH1F *hd2 = (TH1F*) h_data->Project3D("y");
      hd2->SetName("hdinc");
      if(hd2->Integral() == 0) continue;
      hd2->Scale(1./hd2->Integral());
      hd2->Draw();
      TH1F *hm2 = (TH1F*) h_mc->Project3D("y");
      hm2->SetName("hminc");
      if(hm2->Integral() == 0) continue;
      hm2->Scale(1./hm2->Integral());
      hm2->SetMarkerColor(kRed);
      hm2->SetLineColor(kRed);
      hm2->Draw("SAME");

      leg->Draw("SAME");
      TLegend *leg3 = new TLegend(.75,.6,.95,.9);
      leg3->SetFillStyle(0);
      leg3->AddEntry(hd2,"Data","p");
      leg3->AddEntry(hm2,"MC","p");
      leg3->AddEntry("",centbins[j].c_str(),"");
      leg3->AddEntry("","20<p_{T}<60GeV","");
      leg3->Draw("SAME");
      c->Print("fractions.pdf");

      hd2->Divide(hm2);
      hd2->GetXaxis()->SetRangeUser(0,1);
      hd2->GetYaxis()->SetTitle("Data/MC");
      hd2->Draw();
      leg->Draw("SAME");
      leg3->Draw("SAME");
      c->Print("fractions.pdf");
    }
  c->Print("fractions.pdf)");
}
