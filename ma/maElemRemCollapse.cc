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
    Entity* temp_v = vs[1];
    vs[1] = vs[2];
    vs[2] = temp_v;
  }
}

// TODO: Might as well use apf::getFaceFaceAngleInTet()
static double getCosDihedral(Mesh* m, Entity* edge,
			     Entity* face1, Entity* tet1,
			     Entity* face2, Entity* tet2,
			     const apf::Matrix3x3& Q = apf::Matrix3x3(1.,0.,0.,0.,1.,0.,0.,0.,1.))
{
  PCU_ALWAYS_ASSERT(m->getType(face1) == apf::Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(face2) == apf::Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(tet1) == apf::Mesh::TET);
  PCU_ALWAYS_ASSERT(m->getType(tet2) == apf::Mesh::TET);
  apf::Vector3 norm1 = computeFaceNormalAtEdgeInTet(m, tet1, face1, edge, Q);
  apf::Vector3 norm2 = computeFaceNormalAtEdgeInTet(m, tet2, face2, edge, Q);
  // Inverting one of the normals to get dihedral angle
  norm2 = apf::Vector3(0.0, 0.0, 0.0) - norm2;
  return norm1 * norm2;
}

static void showBFaces(Adapt* a, const BFaceMap& bFaceMap, const char* name)
{
  EntityArray faces;
  faces.setSize(bFaceMap.size());
  size_t i = 0;
  for (BFaceMap::const_iterator it = bFaceMap.begin();
       it != bFaceMap.end(); ++it) {
    faces[i] = it->first;
    ++i;
  }
  ma_dbg::createCavityMesh(a, faces, name, apf::Mesh::TRIANGLE);
  
}

ElemRemCollapse::ElemRemCollapse(Adapt* a):
adapter(a)
{
}

bool ElemRemCollapse::markEdges(Mesh* m, Entity* face, bool dryRun)
{
  Entity* edges[3];
  m->getDownward(face, 1, edges);

  bool areNewAnglesGood = true;

  for (int k = 0; k < 3; k++) {
    if(!getFlag(adapter, edges[k], MARKED)) {
      PCU_ALWAYS_ASSERT(bEdgeMap.count(edges[k]) == 0);
      // Unmarked edges can always be marked
      if (dryRun) continue;
      BEdge be1;
      be1.edge = edges[k];
      be1.face1 = face;
      bEdgeMap[edges[k]] = be1;
      setFlag(adapter, edges[k], MARKED);
      // edgesInQueue.push_back(edges[k]);
    } else {
      PCU_DEBUG_ASSERT(bEdgeMap.count(edges[k]) > 0);
      PCU_ALWAYS_ASSERT(bEdgeMap[edges[k]].face1);
      if (bEdgeMap[edges[k]].face2) {
        if (dryRun) return false;
        else
          PCU_ALWAYS_ASSERT_VERBOSE(false, "Attempted to mark an edge with more than one faces already defined.\n");}
      if (dryRun) continue;
      Entity* face1 = bEdgeMap[edges[k]].face1;
      bEdgeMap[edges[k]].face2 = face;
      bEdgeMap[edges[k]].cda = getCosDihedral(m, edges[k],
                                              face, bFaceMap[face].first,
                                              face1, bFaceMap[face1].first);
      // The `!=` acts as XOR. We invert if exactly one of the reference elements is negative
      if (bFaceMap[face1].second != bFaceMap[face].second)
        bEdgeMap[edges[k]].cda *= -1;

      // TODO: is there a better way to compute angles than making the element
      // If element made is negative, the dihedral angle > pi
      bool tryEntMade;
      Entity* tryEnt = removeEdge(edges[k], &tryEntMade);
      double tryEntQual = -1;
      if (tryEnt)
        tryEntQual = adapter->shape->getQuality(tryEnt);
      if (tryEntQual < adapter->input->validQuality) {
        bEdgeMap[edges[k]].cda = -2 - bEdgeMap[edges[k]].cda;
      } else if (tryEntQual < adapter->input->goodQuality) {
        areNewAnglesGood = false;
      }
      if (tryEnt && tryEntMade) destroyElement(adapter, tryEnt);
    }
  }

  if (dryRun)
    return true;
  else
    return areNewAnglesGood;
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
  Mesh* m = adapter->mesh;
  int d = m->getDimension();
  // int bcount = 0;
  // Entity* b;
  apf::Downward bs;

  oldEnts.insert(elems.begin(), elems.end());
  // "Dry run" to mark all boundary faces
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

	markEdges(m, bs[j]);
      }
    }
  }

  ma_dbg::createCavityMesh(adapter, elems, "the_cavity");
  showBFaces(adapter, bFaceMap, "first_bfaces");

  return true;
}

Entity* ElemRemCollapse::removeEdge(Entity* e, bool* elemMade)
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
  // If an element already exists, it was either added to the cavity,
  // or is outside of it. In either case, removing it causes
  // complications.
  Entity* newTet = apf::findElement(m, apf::Mesh::TET, vs);
  if (newTet) {
    if (elemMade) *elemMade = false;
    if(findTetRotation(m, newTet, vs) == -1) {
      // The vertex ordering of returned element is negative of desired.
      // Scrapping the operation.
      return NULL;
    } else {
      return newTet;
    }
  }
  newTet = buildElement(adapter, NULL, apf::Mesh::TET, vs);
  if (elemMade) *elemMade = true;
  if(findTetRotation(m, newTet, vs) == -1) {
    // The vertex ordering of returned element is negative of desired.
    // Scrapping the operation.
    newTet = NULL;
  }
  return newTet;
}

Entity* ElemRemCollapse::removeFace(Entity* face, bool* elemMade)
{
  Mesh* m = adapter->mesh;

  PCU_ALWAYS_ASSERT(m->getType(face) == apf::Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(getFlag(adapter, face, MARKED));

  Entity* edges[3];
  m->getDownward(face, 1, edges);

  bool resultMade = false;
  Entity* result = NULL;
  double qualityToBeat = adapter->input->validQuality;

  for (int i = 0; i < 3; ++i) {
    bool newTetMade;
    Entity* newTet = removeEdge(edges[i], &newTetMade);
    if (!newTet) continue;
    double newQual = adapter->shape->getQuality(newTet);
    if (newQual > qualityToBeat){
      if (result && resultMade) destroyElement(adapter, result);
      result = newTet;
      resultMade = newTetMade;
      qualityToBeat = newQual;
    }
    else
      if (newTet && (result != newTet) && newTetMade)
        destroyElement(adapter, newTet);
  }

  if (elemMade) *elemMade = resultMade;

  return result;
}

bool ElemRemCollapse::removeElement(Entity* e)
{
  Mesh* m = adapter->mesh;
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::TET);
  PCU_ALWAYS_ASSERT(!getFlag(adapter, e, CAV_NEW));
  // TODO: mark this entity
  Entity* fs[4];
  m->getDownward(e, 2, fs);

  for (int j = 0; j < 4; ++j) {
    // switch marks on faces
    setFlags(adapter, fs[j], getFlags(adapter, fs[j]) ^ MARKED);

    if (!getFlag(adapter, fs[j], MARKED)){
      unmarkEdges(m, fs[j]);
    }
  }

  bool canMark = true;
  for (int j = 0; j < 4; ++j) {
    if (getFlag(adapter, fs[j], MARKED)){
      canMark = canMark && markEdges(m, fs[j], true);
    }
  }

  if (!canMark) {
    for (int j = 0; j < 4; ++j) {
      // Switch-back marks on faces
      setFlags(adapter, fs[j], getFlags(adapter, fs[j]) ^ MARKED);
    }
  }

  bool areNewAnglesGood = true;
  for (int j = 0; j < 4; ++j) {
    if (getFlag(adapter, fs[j], MARKED)){
      if (canMark)
        bFaceMap[fs[j]] = std::make_pair(e, !canMark);
      areNewAnglesGood = markEdges(m, fs[j]) && areNewAnglesGood;
    } else {
      bFaceMap.erase(fs[j]);
    }
  }

  newEnts.insert(e);
  if (canMark) {
    // if (!areNewAnglesGood) {
    //   // There are some bad angles introduced, try to undo
    //   if (addElement(e))
    //     return false;
    //   // If failed: what's done is done, go forward
    // }
    ma_dbg::createCavityMesh(adapter, newEnts, "cavity_now");
    setFlag(adapter, e, CAV_NEW);
  } else {
    oldEnts.insert(e);
    // destroyElement(adapter, e);
  }
  
  return canMark;
}

bool ElemRemCollapse::addElement(Entity* e, bool isOld)
{
  Mesh* m = adapter->mesh;
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::TET);

  Entity* fs[4];
  m->getDownward(e, 2, fs);

  for (int j = 0; j < 4; ++j) {
    // switch marks on faces
    setFlags(adapter, fs[j], getFlags(adapter, fs[j]) ^ MARKED);

    if (!getFlag(adapter, fs[j], MARKED)){
      unmarkEdges(m, fs[j]);
    }
  }

  bool canMark = true;
  for (int j = 0; j < 4; ++j) {
    if (getFlag(adapter, fs[j], MARKED)){
      canMark = canMark && markEdges(m, fs[j], true);
    }
  }

  if (!canMark) {
    for (int j = 0; j < 4; ++j) {
      // Switch-back marks on faces
      setFlags(adapter, fs[j], getFlags(adapter, fs[j]) ^ MARKED);
    }
  }

  for (int j = 0; j < 4; ++j) {
    if (getFlag(adapter, fs[j], MARKED)) {
      if (canMark)
        bFaceMap[fs[j]] = std::make_pair(e, canMark);
      markEdges(m, fs[j]);
    } else {
      bFaceMap.erase(fs[j]);
    }
  }
  
  oldEnts.insert(e);
  if (isOld) setFlag(adapter, e, CAV_OLD);
  if (canMark) {
    clearFlag(adapter, e, CAV_NEW);
  } else {
    newEnts.insert(e);
  }

  return canMark;
}

bool ElemRemCollapse::makeNewElements()
{
  size_t count = 0;
  compareEdgeByCosAngle comp(bEdgeMap);
  
  std::stringstream ss;
  ss << "bfaces_" << count;
  showBFaces(adapter, bFaceMap, ss.str().c_str());
  // at worst during marking, all the tet's edges will be added.
  edgesInQueue.reserve(edgesInQueue.size()+6);
  
  edgesInQueue.clear();
  for (BEdgeMap::iterator it = bEdgeMap.begin(); it != bEdgeMap.end(); ++it)
    edgesInQueue.push_back(it->first);

  std::vector<Entity*>::iterator first_it = edgesInQueue.begin();
  std::vector<Entity*>::iterator last_it = edgesInQueue.end();

  std::make_heap(first_it, last_it, comp);

  for (; (first_it != last_it);) {
    std::pop_heap(first_it, last_it, comp);
    Entity* edge = edgesInQueue.back();
    edgesInQueue.pop_back();
    first_it = edgesInQueue.begin();
    last_it = edgesInQueue.end();

    if (!getFlag(adapter, edge, MARKED)) continue;

    bool elementRemoved = false;
    Entity* face1 = bEdgeMap[edge].face1;
    Entity* face2 = bEdgeMap[edge].face2;
    bool newTetMade;
    Entity* newTet = removeFace(face1, &newTetMade);

    if (newTet &&
        (adapter->shape->getQuality(newTet) >
         adapter->input->goodQuality)) {
      elementRemoved = elementRemoved || removeElement(newTet);
    } else {
      if (newTet && newTetMade) destroyElement(adapter, newTet);
    }

    if (!elementRemoved) {
      newTet = removeFace(face2, &newTetMade);

      if (newTet &&
          (adapter->shape->getQuality(newTet) >
           adapter->input->goodQuality)) {
        elementRemoved = elementRemoved || removeElement(newTet);
      } else {
        if (newTet && newTetMade) destroyElement(adapter, newTet);
      }
    }

    if (elementRemoved) {
      edgesInQueue.clear();
      for (BEdgeMap::iterator it = bEdgeMap.begin(); it != bEdgeMap.end(); ++it)
        edgesInQueue.push_back(it->first);

      first_it = edgesInQueue.begin();
      last_it = edgesInQueue.end();
      std::make_heap(first_it, last_it, comp);
      count++;
      std::stringstream ss;
      ss << "bfaces_" << count;
      showBFaces(adapter, bFaceMap, ss.str().c_str());
    }
  }

  ma_dbg::createCavityMesh(adapter, newEnts, "the_new_cavity");
  showBFaces(adapter, bFaceMap, "final_bfaces");

  return false;
}

void ElemRemCollapse::cancel()
{
  destroyNewElements();
  unmark();
  bEdgeMap.clear();
  bFaceMap.clear();
  oldEnts.clear();
  newEnts.clear();
  edgesInQueue.clear();
}

void ElemRemCollapse::transfer()
{
  // keeping as placeholder for whenever element insertion is used
}

void ElemRemCollapse::unmark()
{
  APF_ITERATE(BEdgeMap,bEdgeMap,it) {
    clearFlag(adapter, it->first, MARKED);
  }

  APF_ITERATE(BFaceMap,bFaceMap,it) {
    clearFlag(adapter, it->first, MARKED);
  }

  APF_ITERATE(EntitySet,oldEnts,it) {
    clearFlag(adapter, *it, CAV_OLD);
    clearFlag(adapter, *it, CAV_NEW);
  }

  APF_ITERATE(EntitySet,newEnts,it) {
    clearFlag(adapter, *it, CAV_OLD);
    clearFlag(adapter, *it, CAV_NEW);
  }
}

void ElemRemCollapse::destroyOldElements()
{
  APF_ITERATE(ma::EntitySet,oldEnts,it)
    if (!getFlag(adapter, *it, CAV_NEW))
      destroyElement(adapter, *it);
}

void ElemRemCollapse::destroyNewElements()
{
  APF_ITERATE(ma::EntitySet,newEnts,it)
    if (!getFlag(adapter, *it, CAV_OLD))
      destroyElement(adapter, *it);
}

}
