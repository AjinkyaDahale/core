#include "maMesh.h"
#include "maSnap.h"
#include "maDBG.h"
#include "maElemRemCollapse.h"
#include "maSolutionTransfer.h"
#include "maShapeHandler.h"

namespace ma
{

  static double getCosDihedral(Mesh* m, Entity* edge, Entity* face1, Entity* face2, const apf::Matrix3x3& Q = apf::Matrix3x3(1.,0.,0.,0.,1.,0.,0.,0.,1.))
{
  Entity* vs[2];
  m->getDownward(edge, 0, vs);
  // This should work for linear elements:
  // get normals at the same vertex but from different faces and compare
  apf::Vector3 norm1 = computeFaceNormalAtVertex(m, face1, vs[0], Q);
  apf::Vector3 norm2 = computeFaceNormalAtVertex(m, face2, vs[0], Q);
  // For now inverting one of the normals to get dihedral angle
  norm2 = apf::Vector3(0.0, 0.0, 0.0) - norm2;
  return norm1 * norm2;
}

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
	    // TODO: face pairs
	    double cda = getCosDihedral(m, edges[k], bs[j], bs[k]);
            bEdges.append(std::make_pair(edges[k], cda));
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
