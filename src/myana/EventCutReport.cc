#include "EventCutReport.h"
#include "EventCut.h"

#include <sstream>


EventCutReport::~EventCutReport()
{
    m_cut_names.clear();
    m_cut_details.clear();
    m_number_of_events_passed.clear();
    m_number_of_events_total.clear();
}

void EventCutReport::addResult(EventCut* cut)
{
    if(CutIdx(cut) == k_unregistered_cut) {
        RegisterCut(cut);
        UpdateResults(cut);
    } else {
        UpdateResults(cut);
    }
    return ;
}

void EventCutReport::printReport(std::ostream &os) const
{
    for(unsigned int i = 0; i < m_cut_names.size(); i++){
        os << "Cut: " << m_cut_names.at(i) << std::endl;
        os << "Details: " << m_cut_details.at(i);
        os << "Events Passed: " << m_number_of_events_passed.at(i) <<
         "/" << m_number_of_events_total.at(i) << " (" <<
          100.0 * m_number_of_events_passed.at(i) / m_number_of_events_total.at(i) 
          << "%)" << std::endl;
    }
    return ;
}

unsigned int EventCutReport::CutIdx(EventCut* cut) const
{
    for(unsigned int i = 0; i < m_cut_names.size(); i++){
        if(m_cut_names.at(i) == cut->Name()){
            return i;
        }
    }
    return k_unregistered_cut;
}

void EventCutReport::RegisterCut(EventCut* cut)
{
    m_cut_names.push_back(cut->Name());
    std::stringstream ss;
    cut->identify(ss);
    m_cut_details.push_back(ss.str());
    m_number_of_events_total.push_back(0);
    m_number_of_events_passed.push_back(0);
    return ;
}

void EventCutReport::UpdateResults(EventCut* cut)
{
    unsigned int idx = CutIdx(cut);
    m_number_of_events_total.at(idx)++;
    if(cut->Passed()){
        m_number_of_events_passed.at(idx)++;
    }
    return ;
}

