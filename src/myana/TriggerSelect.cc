#include "TriggerSelect.h"

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1Packetv1.h>

#include <cstdlib>
#include <cstdint>

TriggerSelect::TriggerSelect() : EventCut("TriggerSelect")
{
    // Set the default node name
    SetNodeName("14001");
}

void TriggerSelect::identify(std::ostream &os) const
{
    os << Name() + "::identify: " << std::endl;
    os << "  Node name: " << GetPacket() << std::endl;
    return;
}

bool TriggerSelect::operator()(PHCompositeNode *topNode)
{
    // get minbias info
    auto * node = findNode::getClass<Gl1Packet>(topNode, GetPacket());
    if( !node ){
        std::cerr << Name() + "::operator(PHCompositeNode *topNode) Could not find " << GetPacket() << std::endl;
        exit(-1); // this is a fatal error
    }

    uint64_t gl1_scaled = node->getScaledVector(); 
    for (size_t i = 0; i < m_trigger_ids.size(); ++i) {
        Passed(((gl1_scaled >> m_trigger_ids[i]) & 0x1U) == 0x1U);
        if (!Passed()) {
            if (Verbosity()) {
                std::cout << Name() + "::operator(PHCompositeNode *topNode) Trigger " << m_trigger_ids[i] << " did not fire." << std::endl;
            }
            break; // exit loop if any trigger did not fire
        }
    }

    return Passed();
}
