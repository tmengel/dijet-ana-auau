/*!
 * \file TriggerSelect.h
 * \brief Event selection cut for trigger selection.
 * \author Tanner Mengel <tmengel@bnl.gov>
 */

#ifndef EVENTSELECTION_TRIGGERSELECT_H
#define EVENTSELECTION_TRIGGERSELECT_H

#include "EventCut.h"

class PHCompositeNode;

class TriggerSelect  : public EventCut
{
  public:

    TriggerSelect();
    ~TriggerSelect() override {} 

    void identify(std::ostream &os = std::cout) const override;
    bool operator()(PHCompositeNode* topNode) override;

    void SelectTrigger(size_t trigger_id) { m_trigger_ids.push_back(trigger_id); }
    void SetPacket(const int packet_id) { SetNodeName(std::to_string(packet_id)); }
    int GetPacket() const { return std::stoi(GetNodeName()); }
  
  private:
    std::vector<size_t> m_trigger_ids; // list of trigger IDs to select
  };

#endif // EVENTSELECTION_TriggerSelect_H