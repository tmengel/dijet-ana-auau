/*!
 * \file EventSelector.h
 * \brief SubsysReco module to select events based on user input EventCuts
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Verison: 1.0.1 $
 * \date $Date: 09/10/2024 $
 */
 
#ifndef _EVENTSELECTOR_H_
#define _EVENTSELECTOR_H_

#include <fun4all/SubsysReco.h>

#include <iostream>
#include <string>
#include <vector>

class PHCompositeNode;
class EventCut;
class EventCutReport;

class EventSelector : public SubsysReco
{
 public:

    // Constructor and destructor
    EventSelector(const std::string &name = "EventSelector") : SubsysReco(name) {}
    ~EventSelector() override;

    // Interface for adding and removing cuts
    void AddCut(EventCut * cut);
    void AddCuts( const std::vector<EventCut* > &cuts );
    EventCut * GetCut( const std::string &name );
    void RemoveCut( const std::string &name );

    // Print the cuts
    void PrintCuts( std::ostream &os = std::cout ) const;

    // Standard Fun4All functions
    int InitRun(PHCompositeNode *topNode) override;
    int process_event(PHCompositeNode *topNode) override;
    int ResetEvent(PHCompositeNode *topNode) override;
    int End(PHCompositeNode *topNode) override;

 private:

    std::vector<EventCut *> m_cuts{};
    unsigned int m_nevents_processed{0};
    unsigned int m_nevents_passed{0};
    
    // helper functions
    unsigned int FindCutIdx(const std::string &name);
    bool CutInVector(const std::string &name);

    EventCutReport * m_report{nullptr};
};

#endif // EVENTSELECTION_EVENTSELECTOR_H