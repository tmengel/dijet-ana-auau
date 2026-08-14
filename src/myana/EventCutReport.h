/*!
 * \file EventCutReport.h
 * \brief Helper class to store the results of the event cuts
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Verison: 1.0.1 $
 * \date $Date: 09/13/2024 $
 */

#ifndef _EVENTCUTREPORT_H_
#define _EVENTCUTREPORT_H_

#include <iostream>
#include <string>
#include <vector>
#include <climits>

class EventCut;

class EventCutReport
{
  public:

    EventCutReport() {}
    ~EventCutReport();

    void addResult(EventCut* cut);

    void printReport(std::ostream &os = std::cout) const;

  private:

    std::vector<std::string> m_cut_names{};
    std::vector<std::string> m_cut_details{};
    std::vector<unsigned int> m_number_of_events_total{};
    std::vector<unsigned int> m_number_of_events_passed{};

    unsigned int CutIdx(EventCut* cut) const;
    const unsigned int k_unregistered_cut = {UINT_MAX};

    void RegisterCut(EventCut* cut);
    void UpdateResults(EventCut* cut);

};

#endif // _EVENTCUTREPORT_H_
