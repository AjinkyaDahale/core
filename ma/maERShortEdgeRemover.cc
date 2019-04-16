/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maERShortEdgeRemover.h"
#include "maShape.h"
#include "maAdapt.h"
#include <apfCavityOp.h>
#include <pcu_util.h>

namespace ma {

ERShortEdgeRemover::ERShortEdgeRemover(Adapt* a)
{
  adapter = a;
  elemRemover.Init(a);
  mesh = a->mesh;
  edge = 0;
}

void ERShortEdgeRemover::setEdge(Entity* e)
{
  edge = e;
}

bool ERShortEdgeRemover::requestLocality(apf::CavityOp* o)
{
  Entity* verts[2];
  mesh->getDownward(edge, 0, verts);
  return o->requestLocality(verts, 2);
}

bool ERShortEdgeRemover::run()
{
  // TODO: set old elements
  Upward adjacent;
  Entity* verts[2];
  mesh->getDownward(edge, 0, verts);
  mesh->getAdjacent(edge, mesh->getDimension(), adjacent);
  EntitySet mids(adjacent.begin(), adjacent.end());
  EntityArray oldCav;
  APF_ITERATE(EntitySet, mids, it)
    oldCav.append(*it);
  for (int i = 0; i < 2; i++) {
    mesh->getAdjacent(verts[i], mesh->getDimension(), adjacent);
    APF_ITERATE(Upward, adjacent, it)
      if ( ! mids.count(*it))
        oldCav.append(*it);
  }

  // TODO: any classifngroups and things to ignore
  elemRemover.addClassifnGroups(oldCav);
  if (mesh->getModelType(mesh->toModel(edge)) < mesh->getDimension()) {
    for (int d = mesh->getDimension()-1; d > 0; --d) {
      apf::Adjacent ups;
      mesh->getAdjacent(edge, d, ups);
      elemRemover.setIgnoredModelFaces(ups);
      elemRemover.addClassifnGroups(ups);
    }
  }

  if (mesh->getModelType(mesh->toModel(verts[0])) < mesh->getDimension()) {
    for (int d = mesh->getDimension()-1; d > 0; --d) {
      apf::Adjacent ups;
      mesh->getAdjacent(verts[0], d, ups);
      elemRemover.setIgnoredModelFaces(ups);
      elemRemover.addClassifnGroups(ups);
    }
  }

  if (mesh->getModelType(mesh->toModel(verts[1])) < mesh->getDimension()) {
    for (int d = mesh->getDimension()-1; d > 0; --d) {
      apf::Adjacent ups;
      mesh->getAdjacent(verts[1], d, ups);
      elemRemover.setIgnoredModelFaces(ups);
      elemRemover.addClassifnGroups(ups);
    }
  }

  // TODO: try the element removal
  bool isEROK = elemRemover.setCavity(oldCav);
  if (isEROK) {
    double qualityToBeat = getWorstQuality(adapter, oldCav);
    isEROK = elemRemover.makeNewElements(qualityToBeat);
  }

  // TODO: cleanup or cancel
  if (isEROK) {
    elemRemover.destroyOldElements();
    elemRemover.unmark();
    return true;
  } else {
    elemRemover.cancel();
    return false;
  }
  return false;
}

}
