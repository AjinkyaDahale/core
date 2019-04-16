/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_ERSHORTEDGEREMOVER_H
#define MA_ERSHORTEDGEREMOVER_H

/* short edge removal is related
   to edge collapsing but not identical.
   Here we try to eliminate the short edge
   by any means necessary, typically because
   collapsing it has failed.
   
   We begin with an approach suggested by
   Jie Wan: try to collapse adjacent edges. */

/* note - this algorithm is quite agressive, and requires two
   layers of elements around the target vertices.
   It should be run sparingly. */

#include "maVertRemover.h"
#include "maElemRemover.h"

namespace ma {

class ERShortEdgeRemover
{
  public:
    ERShortEdgeRemover(Adapt* a);
    void setEdge(Entity* e);
    bool requestLocality(apf::CavityOp* o);
    void findEdges();
    bool didImproveQuality();
    bool run();
  private:
    Adapt* adapter;
    Mesh* mesh;
    Entity* edge;
    ElemRemover elemRemover;
};

}

#endif // MA_ERSHORTEDGEREMOVER_H
