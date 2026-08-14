/*!
 * \file ZVertexCut.h
 * \brief Cuts events based on the z-vertex
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Verison: 1.0.1 $
 * \date $Date: 09/10/2024 $
 */

#ifndef _ZVERTEXCUT_H_
#define _ZVERTEXCUT_H_

#include "EventCut.h"

class PHCompositeNode;

class ZVertexCut  : public EventCut
{
  public:

    ZVertexCut(float high, float low = -999.0);
    ~ZVertexCut() override {}

    void identify(std::ostream &os = std::cout) const override;
    bool operator()(PHCompositeNode* topNode) override;
  
};

#endif // EVENTSELECTION_ZVERTEXCUT_H