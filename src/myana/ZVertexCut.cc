#include "ZVertexCut.h"

// phool includes
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

// global vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMapv1.h>

#include <cstdlib>


ZVertexCut::ZVertexCut(float high, float low) : EventCut("ZVertexCut") 
{
        SetNodeName("GlobalVertexMap"); // default node name
        if(low == -999.0){
            SetRange(high, -high);
        } else {
            SetRange(high, low);
        }
}

void ZVertexCut::identify(std::ostream &os) const
{
    os << Name() + "::identify: " << std::endl;
    os << "  Node name: " << GetNodeName() << std::endl;
    os << "  z-vertex range: [" << GetRange().second << ", " << GetRange().first << "]" << std::endl;
    return;
}

bool ZVertexCut::operator()(PHCompositeNode *topNode)
{

    // get global vertex map node
    auto *vrtxMap = findNode::getClass<GlobalVertexMapv1>(topNode, GetNodeName());
    if ( !vrtxMap ) {
        std::cerr << Name() + "::operator(PHCompositeNode *topNode) Could not find " + GetNodeName() << std::endl;
        // exit(-1); // this is a fatal error
        return false;
    } 

    // get primary vertex
    auto * vrtx = vrtxMap->get(0);
    if ( !vrtx ) {
        
        if(Verbosity()){
            std::cout << Name() + "::operator(PHCompositeNode *topNode) Could not find primary vertex in " + GetNodeName() << std::endl;
        }
        return false;

    } 

    Passed((vrtx->get_z() <= GetRange().second) && (vrtx->get_z() >= GetRange().first)); // check z-vertex
    if(!Passed() && Verbosity()){
        std::cout << Name() + "::operator(PHCompositeNode *topNode) Event failed" + Name() + " event selection. Vertex z: " + std::to_string(vrtx->get_z()) << std::endl;
    }

    return Passed();
}
