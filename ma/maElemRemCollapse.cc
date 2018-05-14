#include "maMesh.h"
#include "maSnap.h"
#include "maDBG.h"
#include "maElemRemCollapse.h"
#include "maSolutionTransfer.h"
#include "maShapeHandler.h"

namespace ma
{

ElemRemCollapse::ElemRemCollapse(Adapt* a):
adapter(a)
{
}

bool ElemRemCollapse::setCavity(apf::DynamicArray<Entity*> elems)
{
  int numElems = elems.getSize();
  PCU_ALWAYS_ASSERT(numElems);
  Mesh* m = adapter->mesh;
  // TODO: set boundary
  int d = m->getDimension(); 
  // int bcount = 0;
  // Entity* b;
  apf::Downward bs;

  for (int i = 0; i < numElems; ++i) {
    PCU_ALWAYS_ASSERT_VERBOSE(apf::Mesh::typeDimension[m->getType(elems[i])] == d,
                              "Desired cavity contains entities of different dimension than that of mesh.\n");
    int numBs = m->getDownward(elems[i], d-1, bs);
    // TODO: iter through bs, mark cav boundary entities and set bcount
    for (int j = 0; j < numBs; ++j) {
      setFlags(adapter, bs[j],
               getFlags(adapter, bs[j]) ^ MARKED);
    }
  }

  // Loop through the elems[i]s and then the bs and add bs to boundaryEnts
  for (int i = 0; i < numElems; ++i) {
    int numBs = m->getDownward(elems[i], d-1, bs);
    cavityEnts[i];
    for (int j = 0; j < numBs; ++j) {
      if (getFlag(adapter, bs[j], MARKED)) {
        // clearFlag(adapter, bs[j], MARKED);
        boundaryEnts.append(bs[j]);
        cavityEnts.append(elems[i]);
        cavEntPositive.append(elems[i]);
        Entity* edges[3];
        m->getDownward(bs[j], 1, edges);
        for (int k = 0; k < 3; k++) {
          if(!getFlag(adapter, edges[k], MARKED)) {
	    // TODO: Add the angle here
            bEdges.append(std::make_pair(edges[k], 0.0));
            setFlag(adapter, edges[k], MARKED);
          }
        }
        // This one's just during the debug phase, because
        // fields on faces aren't written to VTK
        Entity* vs[3];
        m->getDownward(bs[j], 0, vs);
        for (int k = 0; k < 3; k++) {
          setFlags(adapter, vs[k], MARKED);
        }
      }
    }
  }

  ma_dbg::createCavityMesh(adapter, elems, "the_cavity");
  ma_dbg::dumpMeshWithFlag(adapter, 0, 0, ma::MARKED, "MARKED", "marked_ents");

  return true;
}

bool ElemRemCollapse::addElement(Entity* // e
                                 )
{
  return false;
}

bool ElemRemCollapse::makeNewElements()
{
  return false;
}

void ElemRemCollapse::cancel()
{
}

void ElemRemCollapse::transfer()
{
}

void ElemRemCollapse::destroyOldElements()
{
}

}
