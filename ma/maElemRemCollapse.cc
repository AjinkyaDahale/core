#include "maMesh.h"
#include "maSnap.h"
#include "maDBG.h"
#include "maElemRemCollapse.h"
#include "maSolutionTransfer.h"
#include "maShapeHandler.h"

namespace ma
{

static int rots[4][4] = {{-1,  1,  2,  0},
                         { 4, -1,  3,  5},
                         { 7,  8, -1,  6},
                         {10,  9, 11, -1}};

// Orients the four vertices of the element-to-be-made (3 downward verts
// of `face` and `opVert`) according to a base tet
static void orientForBuild(Mesh* m, Entity* opVert, Entity* face, Entity* tet,
                           bool dontInvert, Downward vs)
{
  Entity *tvi[4], *fvi[3];
  m->getDownward(tet, 0, tvi);
  m->getDownward(face, 0, fvi);
  // Position of vert 0 and vert 3 of element-to-be-made in
  // old tet
  int p0, p3;
  p0 = findIn(tvi, 4, fvi[0]);
  for (p3 = 0; p3 < 4; ++p3)
    if (-1==findIn(fvi, 3, tvi[p3]))
      break;
  rotateEntity(m->getType(tet), tvi, rots[p0][p3], vs);
  vs[3] = opVert;
  // Inversion logic. The `!=` acts as logical XOR.
  if ((fvi[1]==vs[1])!=dontInvert) {
    vs[1] = fvi[2];
    vs[2] = fvi[1];
  }
}

// TODO: Might as well use apf::getFaceFaceAngleInTet()
static double getCosDihedral(Mesh* m, Entity* edge, Entity* face1,
        Entity* face2, const apf::Matrix3x3& Q = apf::Matrix3x3(1.,0.,0.,0.,1.,0.,0.,0.,1.))
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
  int d = m->getDimension(); 
  // int bcount = 0;
  // Entity* b;
  apf::Downward bs;

  for (int i = 0; i < numElems; ++i) {
    PCU_ALWAYS_ASSERT_VERBOSE(apf::Mesh::typeDimension[m->getType(elems[i])] == d,
                              "Desired cavity contains entities of different dimension than that of mesh.\n");
    setFlag(adapter, elems[i], CAV_OLD);
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
    for (int j = 0; j < numBs; ++j) {
      if (getFlag(adapter, bs[j], MARKED)) {
        // clearFlag(adapter, bs[j], MARKED);
        PCU_ALWAYS_ASSERT(bFaceMap.count(bs[j]) == 0);
        bFaceMap[bs[j]] = std::make_pair(elems[i], true);

        Entity* edges[3];
        m->getDownward(bs[j], 1, edges);
        for (int k = 0; k < 3; k++) {
          if(!getFlag(adapter, edges[k], MARKED)) {
            // TODO: Add the angle here
            // TODO: face pairs
            PCU_ALWAYS_ASSERT(bEdgeMap.count(edges[k]) == 0);
            BEdge be1;
            be1.edge = edges[k];
            be1.face1 = bs[j];
            bEdgeMap[edges[k]] = be1;
            setFlag(adapter, edges[k], MARKED);
          } else {
            PCU_ALWAYS_ASSERT(bEdgeMap[edges[k]].face1);
            bEdgeMap[edges[k]].face2 = bs[j];
            bEdgeMap[edges[k]].cda = getCosDihedral(m, edges[k], bs[j],
                                                   bEdgeMap[edges[k]].face1);
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

bool ElemRemCollapse::removeEdge(Entity* e)
{
  Mesh* m = adapter->mesh;
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::EDGE);
  PCU_ALWAYS_ASSERT(getFlag(adapter, e, MARKED));
  BEdge& bedge = bEdgeMap[e];
  Entity *face1, *opVert;
  Entity *tet = bFaceMap[face1].first;
  bool dontInvert = bFaceMap[face1].second;
  face1 = bedge.face1;
  opVert = getTriVertOppositeEdge(m, bedge.face2, e);
  Entity* vs[4];
  orientForBuild(m, opVert, face1, tet, dontInvert, vs);
  apf::buildElement(m, NULL, apf::Mesh::TET, vs);
  // TODO: make changes in data structures
  // TODO: return appropriately
  return false;
}

bool ElemRemCollapse::addElement(Entity* e)
{
  Mesh* m = adapter->mesh;
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::TET);
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
