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
  // Inversion logic.
  PCU_ALWAYS_ASSERT((fvi[1]==vs[1])||(fvi[1]==vs[2]));
  PCU_ALWAYS_ASSERT((fvi[2]==vs[1])||(fvi[2]==vs[2]));
  if (!dontInvert) {
    vs[1] = fvi[2];
    vs[2] = fvi[1];
  }
}

// TODO: Might as well use apf::getFaceFaceAngleInTet()
static double getCosDihedral(Mesh* m, Entity* edge,
			     Entity* face1, Entity* tet1,
			     Entity* face2, Entity* tet2,
			     const apf::Matrix3x3& Q = apf::Matrix3x3(1.,0.,0.,0.,1.,0.,0.,0.,1.))
{
  apf::Vector3 norm1 = computeFaceNormalAtEdgeInTet(m, tet1, face1, edge, Q);
  apf::Vector3 norm2 = computeFaceNormalAtEdgeInTet(m, tet2, face2, edge, Q);
  // For now inverting one of the normals to get dihedral angle
  norm2 = apf::Vector3(0.0, 0.0, 0.0) - norm2;
  return norm1 * norm2;
}

static void showBFaces(Adapt* a, BFaceMap& bFaceMap, const char* name)
{
  EntityArray faces;
  for (BFaceMap::iterator it = bFaceMap.begin(); it != bFaceMap.end(); ++it) {
    faces.append(it->first);
  }
  ma_dbg::createCavityMesh(a, faces, name, apf::Mesh::TRIANGLE);
  
}

ElemRemCollapse::ElemRemCollapse(Adapt* a):
adapter(a)
{
}

void ElemRemCollapse::markEdges(Mesh* m, Entity* face)
{
  Entity* edges[3];
  m->getDownward(face, 1, edges);
  for (int k = 0; k < 3; k++) {
    if(!getFlag(adapter, edges[k], MARKED)) {
      // TODO: Add the angle here
      // TODO: face pairs
      PCU_ALWAYS_ASSERT(bEdgeMap.count(edges[k]) == 0);
      BEdge be1;
      be1.edge = edges[k];
      be1.face1 = face;
      bEdgeMap[edges[k]] = be1;
      setFlag(adapter, edges[k], MARKED);
    } else {
      PCU_ALWAYS_ASSERT(bEdgeMap[edges[k]].face1);
      PCU_ALWAYS_ASSERT(!bEdgeMap[edges[k]].face2);
  Entity* face1 = bEdgeMap[edges[k]].face1;
      bEdgeMap[edges[k]].face2 = face;
      bEdgeMap[edges[k]].cda = getCosDihedral(m, edges[k], face, bFaceMap[face].first,
                                              face1, bFaceMap[face1].first);
  // The `!=` acts as XOR. We invert if exactly one of the reference elements is negative
  if (bFaceMap[face1].second != bFaceMap[face].second)
    bEdgeMap[edges[k]].cda *= -1;
    }
  }
}

void ElemRemCollapse::unmarkEdges(Mesh* m, Entity* face)
{
  Entity* edges[3];
  m->getDownward(face, 1, edges);
  for (int k = 0; k < 3; k++) {
    PCU_ALWAYS_ASSERT(getFlag(adapter, edges[k], MARKED));
    PCU_ALWAYS_ASSERT(bEdgeMap.count(edges[k]));
    if (bEdgeMap[edges[k]].face2){
      // Ensuring that face1 survives
      if (bEdgeMap[edges[k]].face1 == face){
	bEdgeMap[edges[k]].face1 = bEdgeMap[edges[k]].face2;
      } else
	PCU_ALWAYS_ASSERT(bEdgeMap[edges[k]].face2 == face);
      bEdgeMap[edges[k]].face2 = NULL;
    } else {
      PCU_ALWAYS_ASSERT(bEdgeMap[edges[k]].face1 == face);
      clearFlag(adapter, edges[k], MARKED);
      bEdgeMap.erase(edges[k]);
    }
  }
}

bool ElemRemCollapse::setCavity(apf::DynamicArray<Entity*> elems)
{
  int numElems = elems.getSize();
  PCU_ALWAYS_ASSERT(numElems);
  oldEnts.setSize(numElems);
  Mesh* m = adapter->mesh;
  int d = m->getDimension();
  // int bcount = 0;
  // Entity* b;
  apf::Downward bs;

  // "Dry run" to mark all boundary faces
  for (int i = 0; i < numElems; ++i) {
    PCU_ALWAYS_ASSERT_VERBOSE(apf::Mesh::typeDimension[m->getType(elems[i])] == d,
                              "Desired cavity contains entities of different dimension than that of mesh.\n");
    setFlag(adapter, elems[i], CAV_OLD);
    oldEnts.append(elems[i]);
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

	markEdges(m, bs[j]);
	
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
  showBFaces(adapter, bFaceMap, "first_bfaces");
  ma_dbg::dumpMeshWithFlag(adapter, 0, 0, ma::MARKED, "MARKED", "marked_ents");

  return true;
}

Entity* ElemRemCollapse::removeEdge(Entity* e)
{
  Mesh* m = adapter->mesh;
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::EDGE);
  PCU_ALWAYS_ASSERT(getFlag(adapter, e, MARKED));
  BEdge& bedge = bEdgeMap[e];
  Entity *face1, *opVert;
  Entity *tet;
  bool dontInvert;
  face1 = bedge.face1;
  tet = bFaceMap[face1].first;
  dontInvert  = bFaceMap[face1].second;
  opVert = getTriVertOppositeEdge(m, bedge.face2, e);
  Entity* vs[4];
  orientForBuild(m, opVert, face1, tet, dontInvert, vs);
  Entity* newTet = buildElement(adapter, NULL, apf::Mesh::TET, vs);
  if(findTetRotation(m, newTet, vs) == -1) {
    // The vertex ordering of returned element is negative of desired.
    // This might be because an element as such already exists.
    newTet = NULL;
  }
  return newTet;
}

bool ElemRemCollapse::removeElement(Entity* e)
{
  Mesh* m = adapter->mesh;
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::TET);
  // TODO: mark this entity
  newEnts.append(e);
  ma_dbg::createCavityMesh(adapter, newEnts, "cavity_now");
  setFlag(adapter, e, CAV_NEW);
  // TODO: switch marks on faces
  Entity* fs[4];
  m->getDownward(e, 2, fs);
  for (int j = 0; j < 4; ++j) {
    setFlags(adapter, fs[j], getFlags(adapter, fs[j]) ^ MARKED);

    if (!getFlag(adapter, fs[j], MARKED)){
      bFaceMap.erase(fs[j]);
      unmarkEdges(m, fs[j]);
    }
  }
  for (int j = 0; j < 4; ++j) {
    if (getFlag(adapter, fs[j], MARKED)){
      bFaceMap[fs[j]] = std::make_pair(e, false);
      markEdges(m, fs[j]);
    }
  }
  // TODO: change marks on edges
  return true;
}

bool ElemRemCollapse::addElement(Entity* e)
{
  Mesh* m = adapter->mesh;
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::TET);
  return false;
}

bool ElemRemCollapse::makeNewElements()
{
  size_t count =0;
  for (BEdgeMap::iterator it = bEdgeMap.begin(); it != bEdgeMap.end(); ++it) {
    Entity* newTet = removeEdge(it->first);
    if (newTet && adapter->shape->getQuality(newTet) > 0) {
      removeElement(newTet);
      it = bEdgeMap.begin();
      std::stringstream ss;
      count++;
      ss << "bfaces_" << count;
      showBFaces(adapter, bFaceMap, ss.str().c_str());
    } else {
      adapter->mesh->destroy(newTet);
    }
  }

  ma_dbg::createCavityMesh(adapter, newEnts, "the_new_cavity");
  showBFaces(adapter, bFaceMap, "final_bfaces");

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
