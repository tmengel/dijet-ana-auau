void fractions()
{

  float dphicut = 3*TMath::Pi()/4.;
  int isMC = 0;

  TFile *f = new TFile("/sphenix/tg/tg01/jets/dlis/data/v101/all/TREE_DIJET_v10_1_492_2024p020_v007_gl10-all.root","READ");
  TFile *f10 = new TFile("/sphenix/tg/tg01/jets/dlis/sim/hijing/v14/TREE_JET_SIM_v14_10_new_ProdA_2024-00000030.root","READ");
  TFile *f20 = new TFile("/sphenix/tg/tg01/jets/dlis/sim/hijing/v14/TREE_JET_SIM_v14_20_new_ProdA_2024-00000030.root","READ");
  TFile *f30 = new TFile("/sphenix/tg/tg01/jets/dlis/sim/hijing/v14/TREE_JET_SIM_v14_30_new_ProdA_2024-00000030.root","READ");


  TFile *corrFile = new TFile("/sphenix/user/hanpuj/JES_MC_Calibration/offline/JES_Calib_Default.root", "READ");
  TF1 *f_corr = (TF1*)corrFile->Get("JES_Calib_Default_Func");  

  TFile *frw = new TFile("reweight.root","READ");
  //TH2F *h_rw = (TH2F*) frw->Get("h_ET_EMfrac_recoData");
  //h_rw->SetName("h_rw");
  
  std::vector<float> *jet_et = 0;
  std::vector<float> *jet_e = 0;
  std::vector<float> *jet_pt = 0;
  std::vector<float> *jet_eta = 0;
  std::vector<float> *jet_phi = 0;
  
  std::vector<float> *truth_jet_pt = 0;
  std::vector<float> *truth_jet_eta = 0;
  std::vector<float> *truth_jet_phi = 0;
  std::vector<float> *truth_jet_et = 0;

  ULong64_t gl1_scaled;
  float mbd_vertex_z;
  float mbd_time_zero;
  int cent;
  
  std::vector<float> *jet_EM = 0;
  std::vector<float> *jet_IH = 0;
  std::vector<float> *jet_OH = 0;


  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();



  //set up hists
  TH3F *h_ET_EMfrac_truth = new TH3F("h_ET_EMfrac_truth",";E_{T}^{truth};E_{T}^{EM}/E_{T}",25,10,60,50,-0.5,1.5,4,0,4);
  TH3F *h_ET_EMfrac_match = new TH3F("h_ET_EMfrac_match",";E_{T}^{reco};E_{T}^{EM}/E_{T}",25,10,60,50,-0.5,1.5,4,0,4);
  TH3F *h_ET_EMfrac_reco = new TH3F("h_ET_EMfrac_reco",";E_{T}^{reco};E_{T}^{EM}/E_{T}",25,10,60,50,-0.5,1.5,4,0,4);
  //TH2F *h_ET_EMfrac_smeared = new TH2F("h_ET_EMfrac_smeared",";E_{T}^{smeared};E_{T}^{EM}/E_{T}",25,10,60,50,-0.5,1.5);

  TH1F *h_cent = new TH1F("h_cent","",100,0,100);
  
  
  //and get the response
  TH2F *h_MC_Reso_pt = new TH2F("h_MC_Reso", "" , 25,10,60,100,0,2.);;
  h_MC_Reso_pt->GetXaxis()->SetTitle("p_{T}^{truth} [GeV]");
  h_MC_Reso_pt->GetYaxis()->SetTitle("p_{T}^{reco}/p_{T}^{truth}");

  TH2F *h_MC_Reso_ptrw = new TH2F("h_MC_ResoRW", "" , 25,10,60,100,0,2.);;
  h_MC_Reso_ptrw->GetXaxis()->SetTitle("p_{T}^{truth} [GeV]");
  h_MC_Reso_ptrw->GetYaxis()->SetTitle("p_{T}^{reco}/p_{T}^{truth}");

  TH2F *h_MC_Reso_ptem = new TH2F("h_MC_ResoEM", "" , 25,10,60,100,0,2.);;
  h_MC_Reso_ptem->GetXaxis()->SetTitle("p_{T}^{truth} [GeV]");
  h_MC_Reso_ptem->GetYaxis()->SetTitle("p_{T}^{reco}/p_{T}^{truth}");

  int nfiles = 1;
  if(isMC) nfiles = 3;
  for(int ifile = 0; ifile < nfiles; ifile++)
    {
      TTree *t;
      float weight = 1;
      if(isMC)
	{
	  if(ifile == 0)
	    {
	      t = (TTree*) f10->Get("ttree");
	      weight = 3.997e+06;
	    }
	  if(ifile == 1)
            {
              t = (TTree*) f20->Get("ttree");
              weight = 6.218e+04;
            }
	  if(ifile == 2)
	    {
	      t = (TTree*) f30->Get("ttree");
	      weight = 2.502e+03;
	    }
	}
      else t = (TTree*) f->Get("ttree");
      t->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);
      //t->SetBranchAddress("time_zero", &mbd_time_zero);
      //t->SetBranchAddress("gl1_scaled", &gl1_scaled);
      t->SetBranchAddress("jet_pt_3_sub", &jet_pt);
      t->SetBranchAddress("jet_et_3_sub", &jet_et);
      t->SetBranchAddress("jet_e_3_sub", &jet_e);
      t->SetBranchAddress("jet_eta_3_sub", &jet_eta);
      t->SetBranchAddress("jet_phi_3_sub", &jet_phi);
      t->SetBranchAddress("jet_emcal_3_sub", &jet_EM);
      t->SetBranchAddress("jet_hcalin_3_sub", &jet_IH);
      t->SetBranchAddress("jet_hcalout_3_sub", &jet_OH);
      t->SetBranchAddress("centrality", &cent);
      
      if(isMC)
	{
	  t->SetBranchAddress("truth_jet_pt_3", &truth_jet_pt);
	  t->SetBranchAddress("truth_jet_eta_3", &truth_jet_eta);
	  t->SetBranchAddress("truth_jet_phi_3", &truth_jet_phi);
	}


      int nev = t->GetEntries();
      std::cout<<"Running over " << nev <<" events"<<std::endl;
      for (int i = 0; i < nev; i++)
	{
	  t->GetEntry(i);
	  if(i%100000 == 0) std::cout<<"Event "<<i<<std::endl;
	  if (fabs(mbd_vertex_z) > 60) continue;
	  int centbin = 0;
	  if(cent > 10) centbin = 1;
	  if(cent > 30) centbin = 2;
	  if(cent > 60) centbin = 3;
	  
	  if(isMC)//look at truth jets
	    {
	      int ntruthjets = truth_jet_pt->size();
	      float leadjetpttrue = 0;
	      for (int j = 0; j < ntruthjets;j++)
		{
		  if(truth_jet_pt->at(j) > leadjetpttrue)
                    {
                      leadjetpttrue = truth_jet_pt->at(j);
                    }
		}
	      if(ifile == 0 && leadjetpttrue > 21 || leadjetpttrue < 13) continue;
              if(ifile == 1 && leadjetpttrue > 33 || leadjetpttrue < 21) continue;
	      if(ifile == 2 && leadjetpttrue < 33) continue;
	      
	      for (int j = 0; j < ntruthjets;j++)
                {
		  if(truth_jet_eta->at(j) > 0.7) continue;
		  
		  //match to reco jet
		  int nrecojets = jet_et->size();
		  float matchEta, matchPhi, matchPt, matchE, dR, matchEM;
		  float dRMax = 100;
		  for(int rj = 0; rj < nrecojets; rj++){
		    if(abs(jet_eta->at(rj)) > 0.7 || jet_e->at(rj) < 0) continue;
		    float dEta = truth_jet_eta->at(j) - jet_eta->at(rj);
		    float dPhi = truth_jet_phi->at(j) - jet_phi->at(rj);
		    while(dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
		    while(dPhi < -TMath::Pi()) dPhi += 2*TMath::Pi();
		    dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
		    if(dR < dRMax){
		      matchEta = jet_eta->at(rj);
		      matchPhi = jet_phi->at(rj);
		      matchPt = jet_pt->at(rj);
		      matchE = jet_e->at(rj);
		      dRMax = dR;
		      matchEM = jet_EM->at(rj);
		    }
		  }

		  if(dRMax > 0.3) continue;
		  h_ET_EMfrac_truth->Fill(truth_jet_pt->at(j), matchEM, centbin, weight);
		  h_ET_EMfrac_match->Fill(matchPt, matchEM, centbin, weight);
		
		  //smear EM component by additional 10% and add scale shift
		  float emet = matchPt * matchEM;
		  float emetsmear = emet + gRandom->Gaus(0,truth_jet_pt->at(j)*0.10) + truth_jet_pt->at(j)*0.15;
		  //float emetsmear = emet + truth_jet_pt->at(j)*0.10;
		  float etsmear = matchPt + emetsmear - emet;
		  /*if(truth_jet_pt->at(j) > 10){
		    std::cout<<"Nominal is: "<<emet<<" plus "<<truth_jet_pt->at(j)*0.1<<" and smeared gives "<<emetsmear<<std::endl;
		    std::cout<<"Nominal frac: "<<matchEM<<" alt frac: "<<emetsmear/etsmear<<std::endl;
		    }*/		  
		  //h_ET_EMfrac_smeared->Fill(etsmear, emetsmear/etsmear, weight);
		
		  //and get the response
		  /*h_MC_Reso_pt->Fill(truth_jet_pt->at(j), matchPt/truth_jet_pt->at(j),weight);
		  int bin = h_rw->FindBin(matchPt,matchEM);
		  if (matchPt > 40) bin = h_rw->FindBin(40,matchEM);
		  if (matchPt < 10) bin = h_rw->FindBin(10,matchEM); 

		  if(matchEM > 0 && matchEM < 1 && h_rw->GetBinContent(bin) < 10) h_MC_Reso_ptrw->Fill(truth_jet_pt->at(j), matchPt/truth_jet_pt->at(j),weight*h_rw->GetBinContent(bin));
		  else if(h_rw->GetBinContent(bin) > 10) h_MC_Reso_ptrw->Fill(truth_jet_pt->at(j), matchPt/truth_jet_pt->at(j),weight*10);
		  else h_MC_Reso_ptrw->Fill(truth_jet_pt->at(j), matchPt/truth_jet_pt->at(j),weight);
		  if(matchEM > 0.8) h_MC_Reso_ptem->Fill(truth_jet_pt->at(j), matchPt/truth_jet_pt->at(j),weight);*/
		}
	    }
	  //look at just reco jets

	  int nrecojets = jet_et->size();
	  float leadjetpt = 0;
          float subleadjetpt = 0;
          float leadjetphi = -999;
          float subleadjetphi = -999;
	  float leadjetem = -999;
          float subleadjetem = -999;

	  for(int rj = 0; rj < nrecojets; rj++)
	    {
	      if(abs(jet_eta->at(rj)) > 0.7 || jet_e->at(rj) < 0) continue;

	      if(jet_pt->at(rj) > leadjetpt)
                {
                  subleadjetpt = leadjetpt;
                  subleadjetphi = leadjetphi;
		  subleadjetem = leadjetem;
		  leadjetpt = jet_pt->at(rj);
                  leadjetphi = jet_phi->at(rj);
		  leadjetem = jet_EM->at(rj);
                }
		else if(jet_pt->at(rj) > subleadjetpt)
                {
                  subleadjetpt = jet_pt->at(rj);
                  subleadjetphi = jet_phi->at(rj);
		  subleadjetem = jet_EM->at(rj);
                }
	    }
	  float dphi = fabs(leadjetphi - subleadjetphi);
          while(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
          while(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
          if (fabs(dphi) < 3*TMath::Pi()/4. || subleadjetpt/leadjetpt < 0.3) continue;
	  h_ET_EMfrac_reco->Fill(leadjetpt,leadjetem, centbin,weight);
	  h_ET_EMfrac_reco->Fill(subleadjetpt,subleadjetem, centbin,weight);
	}
    }
  TFile *fout;
  if(isMC) fout = new TFile("MCfractions.root","RECREATE");
  else fout = new TFile("datafractions.root","RECREATE");
  fout->cd();
  h_ET_EMfrac_reco->Write();
  h_ET_EMfrac_truth->Write();
  h_ET_EMfrac_match->Write();
  //h_ET_EMfrac_smeared->Write();
  h_MC_Reso_ptrw->Write();
  h_MC_Reso_pt->Write();
  h_MC_Reso_ptem->Write();
}
