#include "EventSelector.h"
#include "EventCut.h"
#include "EventCutReport.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>


EventSelector::~EventSelector()
{
  for (auto &cut : m_cuts){
    delete cut;
  }
  m_cuts.clear();

  if(m_report){
    delete m_report;
  }
}

void EventSelector::AddCut(EventCut * cut)
{
  if(!CutInVector(cut->Name())){
    m_cuts.push_back(cut);
  } else {
    std::cerr << Name() + "::AddCut(EventCut * cut) Cut " + cut->Name() + " already exists" << std::endl;
    if(Verbosity()){
      GetCut(cut->Name())->identify();
    }
  }
  return ;
}

void EventSelector::AddCuts(const std::vector<EventCut *> &cuts)
{
  for(auto cut: cuts){
    AddCut(cut);
  }
  return ;
}

EventCut* EventSelector::GetCut(const std::string &name)
{
  if(CutInVector(name)){
    return m_cuts.at(FindCutIdx(name));
  } else {
    std::cerr << Name() + "::GetCut(const std::string &name) Cut " + name + " not found" << std::endl;
    return nullptr;
  }
}

void EventSelector::RemoveCut(const std::string &name)
{
  if(CutInVector(name)){
    delete m_cuts.at(FindCutIdx(name));
    m_cuts.erase(m_cuts.begin() + FindCutIdx(name));
  } else {
    std::cerr << Name() + "::RemoveCut(const std::string &name) Cut " + name + " not found" << std::endl;
  }
  return ;
}

void EventSelector::PrintCuts(std::ostream &os) const
{
  os << Name() + "::PrintCuts(std::ostream &os) Printing cuts:" << std::endl;
  for(auto cut: m_cuts){
    cut->identify(os);
  }
  return ;
}

unsigned int EventSelector::FindCutIdx(const std::string &name)
{
  for(unsigned int i = 0; i < m_cuts.size(); i++){
    if(m_cuts.at(i)->Name() == name){
      return i;
    }
  }
  return m_cuts.size();
}

bool EventSelector::CutInVector(const std::string &name)
{
  return FindCutIdx(name) < m_cuts.size();
}

int EventSelector::InitRun(PHCompositeNode * /*topNode*/)
{
  m_report = new EventCutReport();
  for(auto cut: m_cuts){
    cut->Verbosity(Verbosity());
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int EventSelector::process_event(PHCompositeNode *topNode)
{

  bool passed = true;
  m_nevents_processed++;
  for(auto cut: m_cuts){
    if(!(*cut)(topNode)){
      passed = false;
    }
    m_report->addResult(cut);
  }
  if(passed){
    m_nevents_passed++;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  // if we get here, the event failed
  return Fun4AllReturnCodes::ABORTEVENT;

}

int EventSelector::ResetEvent(PHCompositeNode * /*topNode*/)
{
  for(auto cut: m_cuts){
    cut->clear();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int EventSelector::End(PHCompositeNode * /*topNode*/)
{
  // if(Verbosity()){
    std::cout << Name() + "::End(PHCompositeNode *topNode) Printing report" << std::endl;
    std::cout << "Events Summary: " << m_nevents_passed << "/" << m_nevents_processed << " (" << 100.0 * m_nevents_passed / m_nevents_processed << "%)" << std::endl;
    m_report->printReport();
  // }
  return Fun4AllReturnCodes::EVENT_OK;
}



