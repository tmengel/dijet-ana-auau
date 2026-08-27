#include "sPhenixStyle.h"
#include "sPhenixStyle.C"


void makejesjer() 
{

  SetsPhenixStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  
  //if you want to combine multiple files use this
  /*for(int i = 0; i < 100; i++){
    ct->Add(Form("../macro/condor/output_%i.root",i));
    }*/

  vector<float> *eta = 0;
  vector<float> *phi = 0;
  vector<float> *pt = 0;
  vector<float> *e = 0;
  vector<float> *truthEta = 0;
  vector<float> *truthPhi = 0;
  vector<float> *truthPt = 0;
  vector<float> *truthE = 0;
  int cent;
  vector<int> *truthflavor = 0;

  const Float_t pt_range[] = {13,16,21,25,30,35,40,45,50,55,60};
  /*Float_t pt_range[14];
  for(int i = 0; i < 14; i++)
    {
      pt_range[i] = 10+i*5;
      }*/
  const int pt_N = sizeof(pt_range)/sizeof(float) - 1;
  const int resp_N = 60;
  Float_t resp_bins[resp_N+1];
  for(int i = 0; i < resp_N+1; i++){
    resp_bins[i] = 1.2/resp_N * i;
  }
  const int eta_N = 40;
  Float_t eta_bins[eta_N+1];
  for(int i = 0; i < eta_N+1; i++){
    eta_bins[i] = -1.1 + 2.2/eta_N * i;
  }
  const int phi_N = 40;
  Float_t phi_bins[phi_N+1];
  for(int i = 0; i < phi_N+1; i++){
    phi_bins[i] = -TMath::Pi() + 2*TMath::Pi()/phi_N * i;
  }
  const int subet_N = 400;
  Float_t subet_bins[subet_N+1];
  for(int i = 0; i < subet_N+1; i++){
    subet_bins[i] = i*0.5;
  }
  const float cent_bins[] = {-1, 0, 10, 30, 50, 90}; //the first bin is for inclusive in centrality/pp events
  const int cent_N = sizeof(cent_bins)/sizeof(float) - 1;

  TH3F *h_MC_Reso_pt = new TH3F("h_MC_Reso", "" , pt_N, pt_range, resp_N, resp_bins, cent_N, cent_bins);;
  h_MC_Reso_pt->GetXaxis()->SetTitle("p_{T}^{truth} [GeV]");
  h_MC_Reso_pt->GetYaxis()->SetTitle("p_{T}^{reco}/p_{T}^{truth}");

  TH2F *h_pt_reco = new TH2F("h_pt_reco","",subet_N,subet_bins, cent_N, cent_bins);
  h_pt_reco->GetXaxis()->SetTitle("p_{T} [GeV]");
  TH2F *h_pt_true = new TH2F("h_pt_true","",pt_N,pt_range, cent_N, cent_bins);
  h_pt_true->GetXaxis()->SetTitle("p_{T} [GeV]");
  TH2F *h_pt_true_matched = new TH2F("h_pt_true_matched","",subet_N,subet_bins, cent_N, cent_bins);
  h_pt_true_matched->GetXaxis()->SetTitle("p_{T} [GeV]");
  TH3F *h_eta_phi_reco = new TH3F("h_eta_phi_reco","",eta_N, eta_bins, phi_N, phi_bins, cent_N, cent_bins);
  h_eta_phi_reco->GetXaxis()->SetTitle("#eta");
  h_eta_phi_reco->GetYaxis()->SetTitle("#phi");
  TH3F *h_eta_phi_true = new TH3F("h_eta_phi_true","",eta_N, eta_bins, phi_N, phi_bins, cent_N, cent_bins);
  TH3F *h_subEt_pt = new TH3F("h_subEt_pt","",pt_N,pt_range, subet_N, subet_bins, cent_N, cent_bins);


  TGraphErrors *g_jes[cent_N];
  TGraphErrors *g_jer[cent_N];
  for(int icent = 0; icent < cent_N; icent++){
    g_jes[icent] = new TGraphErrors(pt_N);
    g_jer[icent] = new TGraphErrors(pt_N);
  }

  float weight[] = {3.997e+06,6.218e+04,2.502e+03};
  int leadcut[] = {13,16,21};
  for(int ifile = 0; ifile < 3; ifile++)
    {
      TChain * ct = new TChain("ttree");
      ct->Add(Form("/sphenix/tg/tg01/jets/dlis/sim/hijing/v15/TREE_JET_SIM_v15_%i0_new_ProdA_2024-00000030.root",ifile+1));
  
      ct->SetBranchAddress("jet_eta_3_sub",&eta);
      ct->SetBranchAddress("jet_phi_3_sub",&phi);
      ct->SetBranchAddress("jet_pt_3_sub",&pt);
      ct->SetBranchAddress("jet_e_3_sub",&e);
      
      ct->SetBranchAddress("truth_jet_eta_3",&truthEta);
      ct->SetBranchAddress("truth_jet_phi_3",&truthPhi);
      ct->SetBranchAddress("truth_jet_pt_3",&truthPt);
      ct->SetBranchAddress("truth_jet_flavor_3",&truthflavor);
      ct->SetBranchAddress("centrality", &cent);
      int nentries = ct->GetEntries();
      for(int i = 0; i < nentries; i++){
	ct->GetEntry(i);
	
	int njets = truthPt->size();
	int nrecojets = pt->size();
	int isgood = 0;
	
	for(int tj = 0; tj < njets; tj++){
	  //if(truthPt->at(tj) > leadcut[ifile]) isgood = 1;
	  if(ifile == 2 && truthPt->at(tj) > leadcut[ifile]) isgood = 1;
	  else if(truthPt->at(tj) > leadcut[ifile] && truthPt->at(tj) < leadcut[ifile+1]) isgood = 1;
	}
	if(isgood == 0) continue;
	for(int tj = 0; tj < njets; tj++){
	  if(truthflavor->at(tj) != 21) continue;
	  
	  //fill truth hists
	  h_pt_true->Fill(truthPt->at(tj), cent, weight[ifile]);
	  h_pt_true->Fill(truthPt->at(tj), -1, weight[ifile]);
	  h_eta_phi_true->Fill(truthEta->at(tj),truthPhi->at(tj), cent, weight[ifile]);
	  h_eta_phi_true->Fill(truthEta->at(tj),truthPhi->at(tj), -1, weight[ifile]);
	  //do reco to truth jet matching
	  float matchEta, matchPhi, matchPt, matchE, matchsubtracted_et, dR;
	  float dRMax = 100;
	  for(int rj = 0; rj < nrecojets; rj++){
	    if(e->at(rj) < 0) continue;
	    float dEta = truthEta->at(tj) - eta->at(rj);
	    float dPhi = truthPhi->at(tj) - phi->at(rj);
	    while(dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
	    while(dPhi < -TMath::Pi()) dPhi += 2*TMath::Pi();
	    dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
	    if(dR < dRMax){
	      matchEta = eta->at(rj);
	      matchPhi = phi->at(rj);
	      matchPt = pt->at(rj);
	      matchE = e->at(rj);
	      dRMax = dR;
	    }
	  }
	  
	  if(dRMax > 0.3) continue;
	  h_MC_Reso_pt->Fill(truthPt->at(tj),matchPt/truthPt->at(tj), cent, weight[ifile]);
	  h_pt_true_matched->Fill(truthPt->at(tj), cent, weight[ifile]);
	  h_eta_phi_reco->Fill(matchEta,matchPhi, cent, weight[ifile]);
	  h_pt_reco->Fill(matchPt, cent, weight[ifile]);
	  h_MC_Reso_pt->Fill(truthPt->at(tj),matchPt/truthPt->at(tj), -1, weight[ifile]);
	  h_pt_true_matched->Fill(truthPt->at(tj), -1, weight[ifile]);
	  h_eta_phi_reco->Fill(matchEta,matchPhi, -1, weight[ifile]);
	  h_pt_reco->Fill(matchPt, -1, weight[ifile]);
	}
      }
    }

  TCanvas *c = new TCanvas("c","c");
  c->Print("fits04_gluon.pdf("); //uncomment these to get a pdf with all the fits
  TLegend *leg = new TLegend(.25,.2,.6,.5);
  leg->SetFillStyle(0);
  gStyle->SetOptFit(1);
  for(int icent = 0; icent < cent_N; ++icent){
    for (int i = 0; i < pt_N; ++i){
      TF1 *func = new TF1("func","gaus",0,1.2);
      h_MC_Reso_pt->GetXaxis()->SetRange(i+1,i+1);
      h_MC_Reso_pt->GetZaxis()->SetRange(icent+1,icent+1);
      TH1F *h_temp = (TH1F*) h_MC_Reso_pt->Project3D("y");
      h_temp->Fit(func,"","",0,1.2);
      h_temp->Fit(func,"","",func->GetParameter(1)-1.5*func->GetParameter(2),func->GetParameter(1)+1.5*func->GetParameter(2));
      func->SetLineColor(kRed);
      h_temp->Draw();
      leg->AddEntry("",Form("%2.0f%%-%2.0f%%",h_MC_Reso_pt->GetZaxis()->GetBinLowEdge(icent+1),h_MC_Reso_pt->GetZaxis()->GetBinLowEdge(icent+2)),"");
      leg->AddEntry("",Form("%2.0f < p_T < %2.0f GeV",h_MC_Reso_pt->GetXaxis()->GetBinLowEdge(i+1),h_MC_Reso_pt->GetXaxis()->GetBinLowEdge(i+2)),"");
      leg->Draw("SAME");
      /*-------for calculating the JER uncertainty----*/
      float dsigma = func -> GetParError(2);
      float denergy = func -> GetParError(1);
      float sigma = func -> GetParameter(2);
      float energy = func -> GetParameter(1);

      float djer = dsigma/energy + sigma*denergy/pow(energy,2);//correct way to calculate jer
                                                               //accounting for the fact that the 
                                                               //mean response and width are correlated
      c->Print("fits04_gluon.pdf");
      leg->Clear();
      g_jes[icent]->SetPoint(i,0.5*(pt_range[i]+pt_range[i+1]),func->GetParameter(1));
      g_jes[icent]->SetPointError(i,0.5*(pt_range[i+1]-pt_range[i]),func->GetParError(1));
      g_jer[icent]->SetPoint(i,0.5*(pt_range[i]+pt_range[i+1]),func->GetParameter(2)/func->GetParameter(1));
      
      //g_jer[icent]->SetPointError(i,0.5*(pt_range[i+1]-pt_range[i]),func->GetParError(2));
      g_jer[icent]->SetPointError(i,0.5*(pt_range[i+1]-pt_range[i]),djer);
    }
  }
  c->Print("fits04_gluon.pdf)");



  TFile *f_out = new TFile("hists_gluon.root","RECREATE");
  for(int icent = 0; icent < cent_N; ++icent){
    g_jes[icent]->Write(Form("g_jes_cent%i",icent));
    g_jer[icent]->Write(Form("g_jer_cent%i",icent));
  }

  h_pt_true->Write();
  h_eta_phi_true->Write();
  h_pt_true_matched->Write();
  h_eta_phi_reco->Write();
  h_pt_reco->Write();
}
