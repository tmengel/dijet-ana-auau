#include "MinBiasCut.h"

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

// minbias
#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/MinimumBiasInfov1.h>

#include <cstdlib>

MinBiasCut::MinBiasCut() : EventCut("MinBiasCut")
{
    // Set the default node name
    SetNodeName("MinimumBiasInfo");
}

void MinBiasCut::identify(std::ostream &os) const
{
    os << Name() + "::identify: " << std::endl;
    os << "  Node name: " << GetNodeName() << std::endl;
    return;
}

bool MinBiasCut::operator()(PHCompositeNode *topNode)
{
    // get minbias info
    auto * node = findNode::getClass<MinimumBiasInfov1>(topNode, GetNodeName());
    if( !node ){
        std::cerr << Name() + "::operator(PHCompositeNode *topNode) Could not find " + GetNodeName() << std::endl;
        exit(-1); // this is a fatal error
    }

    Passed( node->isAuAuMinimumBias() );
    if(!Passed() && Verbosity()){
       std::cout << Name() + "::operator(PHCompositeNode *topNode) Event failed" + Name() + " event selection" << std::endl;
    }

    return Passed();
}
