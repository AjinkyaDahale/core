/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maVertRemover.h"
#include "maAdapt.h"

namespace ma {

void VertRemover::Init(Adapt* a)
{
  adapter = a;
  mesh = a->mesh;
  collapse.Init(a);
}

void VertRemover::setVert(Entity* v)
{
  vert = v;
}

void VertRemover::findEdges()
{
  mesh->getUp(vert,edges);
}

bool VertRemover::tryToCollapse(Entity* e)
{
  if (!setupCollapse(collapse, e, vert))
    return false;
  double oldQuality = collapse.getOldQuality();
  if ( ! collapse.tryThisDirection(oldQuality))
    return false;
  collapse.destroyOldElements();
  return true;
}

bool VertRemover::run()
{
  Entity* shortEdge[2] = {edges.e[0], edges.e[1]};
  double minLength, l;
  minLength = adapter->sizeField->measure(edges.e[0]);
  for (int i = 0; i < edges.n; ++i) {
    l = adapter->sizeField->measure(edges.e[i]);
    if (l < minLength)
    {
      minLength = l;
      shortEdge[1] = shortEdge[0];
      shortEdge[0] = edges.e[i];
    }
  }
  for (int i = 0; i < 2; ++i)
    if (tryToCollapse(shortEdge[i]))
      return true;
  return false;
}

}
