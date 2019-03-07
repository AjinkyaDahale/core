#include "maMesh.h"
#include "maSnap.h"
#include "maDBG.h"
#include "maElemRemCollapse.h"
#include "maSolutionTransfer.h"
#include "maShapeHandler.h"

namespace ma
{
void ElemRemCollapse::Init(Adapt* a)
{
  Collapse::Init(a);
  elemRemover.Init(a);
}

bool ElemRemCollapse::tryThisDirectionNoCancel(double qualityToBeat)
{
  Mesh* m = adapt->mesh;
  PCU_ALWAYS_ASSERT( ! m->isShared(vertToCollapse));
  // Edges on closure of region have some peculiarities not addressed yet
  // if ( m->getDimension() > m->getModelType(m->toModel(edge)))
  //   return false;
  apf::Adjacent oldCav;
  m->getAdjacent(vertToCollapse, m->getDimension(), oldCav);
  elemRemover.addClassifnGroups(oldCav);
  if (m->getModelType(m->toModel(vertToCollapse)) < m->getDimension()) {
    for (int d = m->getDimension()-1; d > 0; --d) {
      apf::Adjacent ups;
      m->getAdjacent(vertToCollapse, d, ups);
      elemRemover.setIgnoredModelFaces(ups);
      elemRemover.addClassifnGroups(ups);
    }
  }
  elemRemover.setCavity(oldCav);
  bool newCavOK = elemRemover.makeNewElements(qualityToBeat);
  if (!newCavOK)
    return false;
  // since they are okay in a linear sense, now fit and do a quality assessment
  fitElements();

  return true;
}

void ElemRemCollapse::cancel()
{
  cancel(false);
}

void ElemRemCollapse::cancel(bool cavOnly)
{
  if (!cavOnly)
    Collapse::unmark();

  elemRemover.cancel();
}

void ElemRemCollapse::transfer()
{
  // keeping as placeholder for whenever element insertion is used
}

void ElemRemCollapse::unmark(bool cavOnly)
{
  if (!cavOnly)
    Collapse::unmark();

  elemRemover.unmark();
}

void ElemRemCollapse::destroyOldElements()
{
  elemRemover.destroyOldElements();
}

void ElemRemCollapse::destroyNewElements()
{
  elemRemover.destroyNewElements();
}

}
