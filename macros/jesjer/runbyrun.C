#include "sPhenixStyle.h"
#include "sPhenixStyle.C"
#include <vector>
#include <iostream>
using std::vector;

#include <odbc++/connection.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <odbc++/drivermanager.h>
#pragma GCC diagnostic pop
#include <odbc++/resultset.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>

int getRunTime(int runnumber = 0) {

  string runtime;
  odbc::Connection *m_OdbcConnection = nullptr;
  int icount = 0;
  TRandom3 rnd = TRandom3();
  do {
    try {
      m_OdbcConnection = odbc::DriverManager::getConnection("daq","","");
    } catch (odbc::SQLException &e) {
      std::cout << " Exception caught during DriverManager::getConnection" << std::endl;
      std::cout << "Message: " << e.getMessage() << std::endl;
    }
    icount++;
    int wait = (int) 20 + rnd.Uniform()*300;
    if (!m_OdbcConnection) { std::this_thread::sleep_for(std::chrono::seconds(wait)); }
    
  } while (!m_OdbcConnection && icount < 10);
  if (icount == 10) { return 1; }
  std::string sql = "SELECT brtimestamp FROM run WHERE runnumber = " + std::to_string(runnumber) + ";";
  
  odbc::Statement* stmt = m_OdbcConnection->createStatement();
  odbc::ResultSet* resultSet = stmt->executeQuery(sql);
  
  if (resultSet && resultSet->next()) {
    odbc::Timestamp brtimestamp = resultSet->getTimestamp("brtimestamp");
    runtime = brtimestamp.toString(); // brtimestamp is in 'America/New_York' time zone 
  }
  
  std::cout << "Runtime for Run " << runnumber << ": " << runtime << std::endl;

  string output_runtime = runtime;
  for (char& c : output_runtime) 
  {
    if (c == '-' || c == ':') 
    {
      c = ' ';
    }
  }

  int second = -9999;
  int minute = -9999;
  int hour = -9999;
  int day = -9999;
  int month = -9999;
  int year = -9999;

  
  std::istringstream iss(output_runtime);
  iss >> year >> month >> day >> hour >> minute >> second;

  //calculate seconds from 10/9 00:00:00
  int time = (day - 9)*86400;
  time += hour*3600;
  time += minute*60;
  time += second;

  std::cout<<"Time since 10/9 "<<time<<" seconds"<<std::endl;
  
  delete resultSet;
  delete stmt;
  delete m_OdbcConnection;
  return time;
  
}

void runbyrun()
{
  SetsPhenixStyle();
  TCanvas *c = new TCanvas();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();


  TH1F *h_jets = new TH1F("h_jets","",800,54200,55000);
  TH1F *h_gl1 = new TH1F("h_gl1","",800,54200,55000);
  
  
  std::vector<float> *jet_et = 0;
  std::vector<float> *jet_e = 0;
  std::vector<float> *jet_pt = 0;
  std::vector<float> *jet_eta = 0;
  std::vector<float> *jet_phi = 0;

  ULong64_t live_first[64];
  ULong64_t live_last[64];

  int runs[] = {54263,54264,54266,54268,54269,54360,54361,54362,54364,54373,54374,54377,54379,54404,54408,54412,54414,54415,54422,54423,54424,54425,54431,54432,54433,54435,54436,54437,54438,54451,54463,54465,54466,54467,54468,54469,54470,54471,54479,54481,54482,54483,54484,54485,54486,54487,54488,54489,54490,54491,54493,54494,54495,54496,54497,54498,54527,54528,54531,54532,54533,54534,54535,54543,54544,54545,54546,54547,54548,54549,54584,54587,54588,54589,54590,54592,54593,54594,54595,54596,54597,54598,54599,54600,54602,54603,54604,54676,54682,54707,54747,54748,54749,54803,54804,54806,54849,54850,54862,54863,54865,54866,54867,54871,54872,54873,54896,54897,54899,54900,54906,54907,54908,54909,54910,54911,54912,54913,54914,54915,54916,54917,54918,54919,54920,54921,54931,54933,54934,54935,54936,54937,54938,54944,54945,54948,54949,54951,54952,54953,54965,54966,54967,54968,54969,54971,54972,54973,54974};
  int nruns = sizeof(runs) / sizeof(runs[0]);

  TGraphErrors *g_time = new TGraphErrors(nruns);
  
  for(int irun = 0; irun < nruns; irun++)
    {
      TFile *f = new TFile(Form("/sphenix/tg/tg01/jets/dlis/data/v101/runs/TREE_DIJET_v10_1_492_2024p020_v007_gl10-000%i.root",runs[irun]),"READ");
      TTree *t = (TTree*) f->Get("ttree");
      
      t->SetBranchAddress("jet_pt_3_sub", &jet_pt);
      t->SetBranchAddress("jet_et_3_sub", &jet_et);
      t->SetBranchAddress("jet_e_3_sub", &jet_e);
      t->SetBranchAddress("jet_eta_3_sub", &jet_eta);
      t->SetBranchAddress("jet_phi_3_sub", &jet_phi);
      
      
      int nev = t->GetEntries();
      std::cout<<"Running over " << nev <<" events"<<std::endl;
      int njets = 0;
      for (int i = 0; i < nev; i++)
	{
	  t->GetEntry(i);
	  if(i%100000 == 0) std::cout<<"Event "<<i<<std::endl;
	  int nj = jet_pt->size();
	  for(int ij = 0; ij < nj; ij++)
	    {
	      if(jet_pt->at(ij) > 20) h_jets->Fill(runs[irun]);;
	    }
	}

      TTree *te = (TTree*)f->Get("ttree_end");
      te->SetBranchAddress("live_first",live_first);
      te->SetBranchAddress("live_last",live_last);
      //te->SetBranchAddress("live_scalers",&live_scalers);
      //te->SetBranchAddress("scaled_scalers",&scaled_scalers);
      Long64_t nseg = te->GetEntries();
      for(Long64_t i=0;i<nseg;++i){
	te->GetEntry(i);
	ULong64_t seg_live = live_last[10] - live_first[10];
	//std::cout<<live_last[10]<<std::endl;
	h_gl1->Fill(runs[irun],seg_live);
      }
      int time  = getRunTime(runs[irun]);
      int bin = h_jets->FindBin(runs[irun]);
      g_time->SetPoint(irun,time,h_jets->GetBinContent(bin)/h_gl1->GetBinContent(bin));
      g_time->SetPointError(irun,0,h_jets->GetBinContent(bin)/h_gl1->GetBinContent(bin)*sqrt(h_jets->GetBinError(bin)/h_jets->GetBinContent(bin)*h_jets->GetBinError(bin)/h_jets->GetBinContent(bin) + h_gl1->GetBinError(bin)/h_gl1->GetBinContent(bin)*h_gl1->GetBinError(bin)/h_gl1->GetBinContent(bin)));
    }

  TFile *fout = new TFile("runbyrun.root","RECREATE");
  h_jets->Write();
  h_gl1->Write();
  g_time->Write("g_time");
}
