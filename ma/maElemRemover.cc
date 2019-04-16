#include "maMesh.h"
#include "maSnap.h"
#include "maDBG.h"
#include "maElemRemover.h"
#include "maSolutionTransfer.h"
#include "maShapeHandler.h"

namespace ma
{

// Using rots[i][j] will result in the the v[i] to be v'[0] after rotation
// and v[j] to be v'[4], where v comes from getDownward(tet, 0, v).
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

// Returns the cos of the angle between normals on `face1` and `face2` at the
// mid-point of `edge`, pointing away from `tet1` and `tet2` respectively.
// NOTE: Might as well modify apf::getFaceFaceAngleInTet()
static double getCosDihedral(Adapt* a, Entity* edge,
                             Entity* face1, Entity* tet1,
                             Entity* face2, Entity* tet2)
{
  Mesh* m = a->mesh;

  PCU_ALWAYS_ASSERT(m->getType(face1) == apf::Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(face2) == apf::Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(tet1) == apf::Mesh::TET);
  PCU_ALWAYS_ASSERT(m->getType(tet2) == apf::Mesh::TET);

  Matrix Q;
  apf::MeshElement* me = apf::createMeshElement(m, edge);
  Vector center(0.0, 0.0, 0.0);
  a->sizeField->getTransform(me,center,Q);
  apf::destroyMeshElement(me);

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

static double getEstimatedTetQual(Mesh* m, Entity** verts)
{
  Vector v[4];
  for (size_t i = 0; i < 4; ++i) {
    m->getPoint(verts[i], 0, v[i]);
  }
  return measureLinearTetQuality(v);
}

static double getTetVol(Mesh* m, Entity** verts)
{
  Vector v[4];
  for (size_t i = 0; i < 4; ++i) {
    m->getPoint(verts[i], 0, v[i]);
  }
  return cross((v[1] - v[0]), (v[2] - v[0])) * (v[3] - v[0]);
}

// static void showNewEnts(Adapt* a, const EntitySet& newEnts, const char* name)
// {
//   EntityArray tets;
//   tets.setSize(newEnts.size());
//   size_t i = 0;
//   APF_ITERATE(EntitySet,newEnts,it) {
//     if (getFlag(a, *it, CAV_NEW)) {
//       tets[i] = *it;
//       ++i;
//     }
//   }
//   tets.setSize(i);
//   ma_dbg::createCavityMesh(a, tets, name, apf::Mesh::TET);
// }

Entity* ElemRemover::buildOrFind(Adapt* a, int type, Entity** vs,
                                     bool* elemMade, bool* elemInverted)
{
  Mesh* m = a->mesh;

  // If an element already exists, it was either added to the cavity,
  // or is outside of it. In either case, removing it causes
  // complications.
  Entity* newTet = apf::findElement(m, type, vs);
  if (newTet) {
    if (elemMade) *elemMade = false;
  } else {
    ElementBuilder2 b(m,type,classifnGroups,a->buildCallback,elemMade);
    newTet = b.run(type,vs);
    if (elemMade) *elemMade = true;
  }

  if (elemInverted) *elemInverted = false;
  if (findTetRotation(m, newTet, vs) == -1) {
    // The vertex ordering of returned element is negative of desired.
    // Scrapping the operation.
    if (elemInverted)
      *elemInverted = true;
    else
      return NULL;
  }

  return newTet;
}

bool ElemRemover::edgeIntersectsFace(Adapt* a, Entity* face, Entity* edge)
{
  Mesh* m = a->mesh;
  Entity* fv[3]; // face vertices
  Entity* ev[2]; // edge vertices
  m->getDownward(edge, 0, ev);
  m->getDownward(face, 0, fv);

  double tol = 1e-10;

  // If one vertex is common to both edge and face, only non-trivial
  // intersection is if they are coplanar. We let it be handled somewhere else.
  if (findIn(fv, 3, ev[0])!=-1 || findIn(fv, 3, ev[1])!=-1)
    return false;

  Vector fp[3]; // face point coordinates
  Vector ep[2]; // edge point coordinates
  for (size_t i = 0; i < 3; ++i) {
    m->getPoint(fv[i], 0, fp[i]);
  }
  for (size_t i = 0; i < 2; ++i) {
    m->getPoint(ev[i], 0, ep[i]);
  }

  // Following the Moller-Trumbore intersection check algorithm, without
  // back face culling. We ultimately estimate parametric coordinates of the
  // intersection point of the line of the edge and the plane of the face w.r.t.
  // the edge (t) and the face (u, v, 1-u-v). All intermediate entities
  // computed are for optimization purposes.
  Vector D = ep[1] - ep[0];
  Vector E0 = fp[1] - fp[0];
  Vector E1 = fp[2] - fp[0];
  Vector A = cross(D, E1);
  double det = A * E0;

  // To make comparison dimensionless, we compare the dot product with the
  // product of magnitudes.
  double detTol = std::fabs(tol * A.getLength() * E0.getLength());

  Vector T, B;
  if (det > detTol) {
    T = ep[0] - fp[0];
  } else if (det < -detTol) {
    // For further computations, we assume det > 0
    det = -det;
    T = fp[0] - ep[0];
  } else {
    return false;
  }

  // The parameters we compute are dimensionless, but we instead choose to
  // scale the tolerance and avoid some division.
  double paramTol = tol * det;

  // Once we have the parametric coordinates, we theoretically want them to be
  // between 0 and 1 to say intersection happens. However, we choose to err on
  // side of ensuring there is no intersection and change the limits to
  // -tol and 1 + tol.
  double u = A * T;

  B = cross(T, E0);
  double v = B * D;

  double t = B * E1;
  if ((t < -paramTol) || (t - det > paramTol))
    return false;

  if ((u > -0.2*det) && (v > -0.2*det) && (det - u - v > -0.2*det)) {
    Entity* qcheckvs[4];
    qcheckvs[2] = ev[0]; qcheckvs[3] = ev[1];
    if ((u > -0.2*det) && (u < -paramTol)) {
      qcheckvs[0] = fv[0]; qcheckvs[1] = fv[2];
      if (getEstimatedTetQual(m, qcheckvs) < qualToBeat)
	return true;
    }
    if ((v > -0.2*det) && (v < -paramTol)) {
      qcheckvs[0] = fv[0]; qcheckvs[1] = fv[1];
      if (getEstimatedTetQual(m, qcheckvs) < qualToBeat)
	return true;
    }
    if ((det - u - v > -0.2*det) && (det - u - v < -paramTol)) {
      qcheckvs[0] = fv[1]; qcheckvs[1] = fv[2];
      if (getEstimatedTetQual(m, qcheckvs) < qualToBeat)
	return true;
    }
  }

  if ((u < -paramTol) || (u - det > paramTol))
    return false;

  if ((v < -paramTol) || (v + u - det > paramTol)) {
    return false;
  }

  return true;
}

bool ElemRemover::newTetClear(Adapt* a, Entity* tet)
{
  Mesh* m = a->mesh;
  Entity* fs[4];
  m->getDownward(tet, 2, fs);

  for (size_t i = 0; i < 4; ++i) {
    if (getFlag(a, fs[i], MARKED)) continue;
    APF_ITERATE(BEdgeMap,bEdgeMap,it) {
      if (!ignoredFaceClassifn.empty() &&
          ignoredFaceClassifn.count(m->toModel(it->first)))
        continue;
      if (edgeIntersectsFace(a, fs[i], it->first))
        return false;
    }
  }

  return true;
}

bool ElemRemover::markEdges(Mesh* m, Entity* face, bool dryRun)
{
  Entity* edges[3];
  m->getDownward(face, 1, edges);

  bool areNewAnglesGood = true;

  for (int k = 0; k < 3; k++) {
    if(!getFlag(adapt, edges[k], MARKED)) {
      PCU_ALWAYS_ASSERT(bEdgeMap.count(edges[k]) == 0);
      // Unmarked edges can always be marked
      if (dryRun) continue;
      BEdge be1;
      be1.edge = edges[k];
      be1.face1 = face;
      bEdgeMap[edges[k]] = be1;
      setFlag(adapt, edges[k], MARKED);
      // edgesInQueue.push_back(edges[k]);
    } else {
      PCU_DEBUG_ASSERT(bEdgeMap.count(edges[k]) > 0);
      PCU_ALWAYS_ASSERT(bEdgeMap[edges[k]].face1);
      if (bEdgeMap[edges[k]].face2) {
        if (dryRun) return false;
        else {
          ma_dbg::createCavityMesh(adapt, oldEnts, "crossedge_old_ents");
          ma_dbg::createCavityMesh(adapt, newEntsArray, "crossedge_new_ents", apf::Mesh::TET);
          showBFaces(adapt, bFaceMap, "crossedge_bfaces");
          PCU_ALWAYS_ASSERT_VERBOSE(false, "Attempted to mark an edge with more"
                                    " than one faces already defined.\n");
        }
      }
      if (dryRun) continue;
      Entity* face1 = bEdgeMap[edges[k]].face1;
      bEdgeMap[edges[k]].face2 = face;
      bEdgeMap[edges[k]].cda = getCosDihedral(adapt, edges[k],
                                              face, bFaceMap[face].first,
                                              face1, bFaceMap[face1].first);
      // The `!=` acts as XOR. We invert if exactly one of the reference
      // elements is negative
      if (bFaceMap[face1].second != bFaceMap[face].second)
        bEdgeMap[edges[k]].cda *= -1;

      // A check to see if the dihedral angle is < pi or not
      // TODO: is there a better way to compute angles than making the element
      // If element made is negative, the dihedral angle > pi
      Entity* vs[4];
      orientForBuild(m, getTriVertOppositeEdge(m, face, edges[k]), face1,
                     bFaceMap[face1].first, bFaceMap[face1].second, vs);
      if (!(getTetVol(m, vs) > 0)) {
        bEdgeMap[edges[k]].cda = -2 - bEdgeMap[edges[k]].cda;
      } else if (bEdgeMap[edges[k]].cda > 0 && \
                 getEstimatedTetQual(m, vs) < this->qualToBeat) {
        // TODO: could we use a less arbitrary condition for too acute angles?
        // Just reply that the angle is too acute for our liking
        if (ignoredFaceClassifn.empty() ||
            ignoredFaceClassifn.count(m->toModel(edges[k]))==0)
          areNewAnglesGood = false;
      }
    }
  }

  if (dryRun)
    return true;
  else
    return areNewAnglesGood;
}

void ElemRemover::unmarkEdges(Mesh* m, Entity* face)
{
  Entity* edges[3];
  m->getDownward(face, 1, edges);
  for (int k = 0; k < 3; k++) {
    PCU_ALWAYS_ASSERT(getFlag(adapt, edges[k], MARKED));
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
      clearFlag(adapt, edges[k], MARKED);
      bEdgeMap.erase(edges[k]);
    }
  }
}

void ElemRemover::Init(Adapt* a)
{
  adapt = a;
}

bool ElemRemover::setCavity(apf::DynamicArray<Entity*>& elems)
{
  int numElems = elems.getSize();
  PCU_ALWAYS_ASSERT(numElems);
  Mesh* m = adapt->mesh;
  int d = m->getDimension();
  // int bcount = 0;
  // Entity* b;
  apf::Downward bs;

  // Expected that anything in these to be destroyed will have
  // been destroyed by the time we're here
  PCU_DEBUG_ASSERT(oldEnts.empty());
  PCU_DEBUG_ASSERT(newEnts.empty());
  PCU_DEBUG_ASSERT(newEntsArray.getSize() == 0);
  PCU_DEBUG_ASSERT(bFaceMap.empty());
  PCU_DEBUG_ASSERT(bEdgeMap.empty());

  modelEnt = m->toModel(elems[0]);

  oldEnts.insert(elems.begin(), elems.end());
  // "Dry run" to mark all boundary faces
  for (int i = 0; i < numElems; ++i) {
    PCU_ALWAYS_ASSERT_VERBOSE(apf::Mesh::typeDimension[m->getType(elems[i])] == d,
                              "Desired cavity contains entities of different "
                              "dimension than that of mesh.\n");
    // We do not want to cross model region boundaries
    if (m->toModel(elems[i]) != modelEnt) continue;
    setFlag(adapt, elems[i], CAV_OLD);
    int numBs = m->getDownward(elems[i], d-1, bs);
    // TODO: iter through bs, mark cav boundary entities and set bcount
    for (int j = 0; j < numBs; ++j) {
      setFlags(adapt, bs[j],
               getFlags(adapt, bs[j]) ^ MARKED);
    }
  }

  bool canMarkAllBEdges = true;

  // Loop through the elems[i]s and then the bs and add bs to boundaryEnts
  for (int i = 0; i < numElems; ++i) {
    int numBs = m->getDownward(elems[i], d-1, bs);
    for (int j = 0; j < numBs; ++j) {
      if (getFlag(adapt, bs[j], MARKED)) {
        // clearFlag(adapt, bs[j], MARKED);
        if (bFaceMap.count(bs[j]) != 0) {
          ma_dbg::createCavityMesh(adapt, elems, "badfaces_old_ents");
          ma_dbg::createCavityMesh(adapt, newEntsArray, "badfaces_new_ents", apf::Mesh::TET);
          showBFaces(adapt, bFaceMap, "badfaces_bfaces");
          print("Something bad happening at i = %d, j = %d", i, j);
          apf::fail("A face is being added twice to bfacemap");
        }
        bFaceMap[bs[j]] = std::make_pair(elems[i], true);

        if (markEdges(m, bs[j], true)) {
          markEdges(m, bs[j]);
        } else {
          canMarkAllBEdges = false;
        }
      }
    }
  }

  // ma_dbg::createCavityMesh(adapt, elems, "the_cavity");
  // showBFaces(adapt, bFaceMap, "first_bfaces");

  return canMarkAllBEdges;
}

Entity* ElemRemover::removeEdge(Entity* e, bool* elemMade)
{
  Mesh* m = adapt->mesh;

  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::EDGE);
  PCU_ALWAYS_ASSERT(getFlag(adapt, e, MARKED));

  BEdge& bedge = bEdgeMap[e];
  Entity *face1, *face2, *opVert;
  Entity *tet;
  bool dontInvert;
  face1 = bedge.face1;
  face2 = bedge.face2;

  if (!ignoredFaceClassifn.empty())
    if(ignoredFaceClassifn.count(m->toModel(face1)) ||
       ignoredFaceClassifn.count(m->toModel(face2))) {
      if (elemMade) *elemMade = false;
      return NULL;
    }

  tet = bFaceMap[face1].first;
  dontInvert  = bFaceMap[face1].second;
  opVert = getTriVertOppositeEdge(m, face2, e);

  Entity* vs[4];
  orientForBuild(m, opVert, face1, tet, dontInvert, vs);

  return buildOrFind(adapt, apf::Mesh::TET, vs, elemMade, NULL);
}

bool ElemRemover::removeElement(Entity* e)
{
  Mesh* m = adapt->mesh;
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::TET);
  if (getFlag(adapt, e, CAV_NEW)) {
    return false;
    // ma_dbg::createCavityMesh(adapt, newEntsArray, "new_ents_array", apf::Mesh::TET);
    // showBFaces(adapt, bFaceMap, "bfaces");
    // EntityArray offender;
    // offender.append(e);
    // ma_dbg::createCavityMesh(adapt, offender, "offending_tet", apf::Mesh::TET);
    // apf::fail("Some tet is being removed twice!\n");
  }
  // TODO: mark this entity
  Entity* fs[4];
  m->getDownward(e, 2, fs);

  newEnts.insert(e);

  if (!newTetClear(adapt, e)) {
    oldEnts.insert(e);
    return false;
  }

  for (int j = 0; j < 4; ++j) {
    // switch marks on faces
    setFlags(adapt, fs[j], getFlags(adapt, fs[j]) ^ MARKED);

    if (!getFlag(adapt, fs[j], MARKED)){
      unmarkEdges(m, fs[j]);
    }
  }

  bool canMark = true;
  for (int j = 0; j < 4; ++j) {
    if (getFlag(adapt, fs[j], MARKED)){
      canMark = canMark && markEdges(m, fs[j], true);
    }
  }

  if (!canMark) {
    for (int j = 0; j < 4; ++j) {
      // Switch-back marks on faces
      setFlags(adapt, fs[j], getFlags(adapt, fs[j]) ^ MARKED);
    }
  }

  bool areNewAnglesGood = true;
  for (int j = 0; j < 4; ++j) {
    if (getFlag(adapt, fs[j], MARKED)){
      if (canMark)
        bFaceMap[fs[j]] = std::make_pair(e, !canMark);
      areNewAnglesGood = markEdges(m, fs[j]) && areNewAnglesGood;
    } else {
      bFaceMap.erase(fs[j]);
    }
  }

  if (canMark) {
    if (!areNewAnglesGood) {
      // There are some bad angles introduced, try to undo
      if (addElement(e)) {
        return false;
      } else {
        ma_dbg::createCavityMesh(adapt, oldEnts, "oldEnts_cantadd", Mesh::TET);
        ma_dbg::createCavityMesh(adapt, newEnts, "newEnts_cantadd", Mesh::TET);
        ma_dbg::createCavityMesh(adapt, newEntsArray, "newEntsArray_cantadd", Mesh::TET);
        apf::fail("cannot add back!");
      }
      // If failed: what's done is done, go forward
    }
    // ma_dbg::createCavityMesh(adapt, newEnts, "cavity_now");
    setFlag(adapt, e, CAV_NEW);
    newEntsArray.append(e);
  } else {
    oldEnts.insert(e);
    // destroyElement(adapt, e);
  }

  return canMark;
}

bool ElemRemover::addElement(Entity* e, bool isOld)
{
  Mesh* m = adapt->mesh;
  PCU_ALWAYS_ASSERT(m->getType(e) == apf::Mesh::TET);

  Entity* fs[4];
  m->getDownward(e, 2, fs);

  for (int j = 0; j < 4; ++j) {
    // switch marks on faces
    setFlags(adapt, fs[j], getFlags(adapt, fs[j]) ^ MARKED);

    if (!getFlag(adapt, fs[j], MARKED)){
      unmarkEdges(m, fs[j]);
    }
  }

  bool canMark = true;
  for (int j = 0; j < 4; ++j) {
    if (getFlag(adapt, fs[j], MARKED)){
      canMark = canMark && markEdges(m, fs[j], true);
    }
  }

  if (!canMark) {
    for (int j = 0; j < 4; ++j) {
      // Switch-back marks on faces
      setFlags(adapt, fs[j], getFlags(adapt, fs[j]) ^ MARKED);
    }
  }

  for (int j = 0; j < 4; ++j) {
    if (getFlag(adapt, fs[j], MARKED)) {
      if (canMark)
        bFaceMap[fs[j]] = std::make_pair(e, canMark);
      markEdges(m, fs[j]);
    } else {
      bFaceMap.erase(fs[j]);
    }
  }

  oldEnts.insert(e);
  if (isOld) setFlag(adapt, e, CAV_OLD);
  if (canMark) {
    clearFlag(adapt, e, CAV_NEW);
  } else {
    newEnts.insert(e);
  }

  return canMark;
}

void ElemRemover::setIgnoredModelFaces(EntityArray& ups)
{
  Mesh* m = adapt->mesh;
  int dim = m->getDimension();
  APF_ITERATE(EntityArray,ups,it) {
    Model* c = m->toModel(*it);
    if (m->getModelType(c) < dim)
      ignoredFaceClassifn.insert(c);
  }
}

void ElemRemover::addClassifnGroups(EntityArray& ups)
{
  Mesh* m = adapt->mesh;
  APF_ITERATE(EntityArray,ups,it) {
    Model* c = m->toModel(*it);
    apf::Downward vs;
    int nv = m->getDownward(*it, 0, vs);
    for (int i = 0; i < nv; ++i)
      classifnGroups[c].insert(vs[i]);
  }
}

bool ElemRemover::makeNewElements(double qualityToBeat)
{
  this->qualToBeat = qualityToBeat;
  size_t count = 0;
  CompareEdgeByCosAngle comp(bEdgeMap);

  // std::stringstream ss;
  // ss << "bfaces_" << count;
  // showBFaces(adapt, bFaceMap, ss.str().c_str());
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

    if (!getFlag(adapt, edge, MARKED)) continue;

    bool elementRemoved = false;
    // Entity* face1 = bEdgeMap[edge].face1;
    // Entity* face2 = bEdgeMap[edge].face2;
    bool newTetMade;
    Entity* newTet = removeEdge(edge, &newTetMade);

    if (newTet && (adapt->shape->getQuality(newTet) > qualityToBeat)) {
    // if (newTet &&
    //     (adapt->shape->getQuality(newTet) > qualityToBeat) &&
    //     (adapt->sizeField->measure(
    //       getTetEdgeOppositeEdge(adapt->mesh, newTet, edge)) < 1.5)) {
      elementRemoved = elementRemoved || removeElement(newTet);
    } else {
      if (newTet && newTetMade) destroyElement(adapt, newTet);
      else PCU_ALWAYS_ASSERT((!newTet) || adapt->mesh->toModel(newTet) != NULL);
    }

    if (elementRemoved) {
      edgesInQueue.clear();
      for (BEdgeMap::iterator it = bEdgeMap.begin(); it != bEdgeMap.end(); ++it)
        edgesInQueue.push_back(it->first);

      first_it = edgesInQueue.begin();
      last_it = edgesInQueue.end();
      std::make_heap(first_it, last_it, comp);
      count++;
      // std::stringstream ss;
      // ss << "bfaces_" << count;
      // showBFaces(adapt, bFaceMap, ss.str().c_str());
    }
  }

  // showNewEnts(adapt, newEnts, "the_new_cavity");
  // showBFaces(adapt, bFaceMap, "final_bfaces");

  bool retVal = false;
  if (!ignoredFaceClassifn.empty()) {
    bool ignoreRest = true;
    APF_ITERATE(BFaceMap,bFaceMap,it) {
      ignoreRest = ignoreRest &&
        (ignoredFaceClassifn.count(adapt->mesh->toModel(it->first)));
    }
    retVal = ignoreRest;
  } else {
    retVal = (bEdgeMap.size() == 0);
  }

  if (adapt->printERFails && !retVal && (adapt->erfails++ % 20 == 0)) {
    reportState(adapt->erfails);
  }

  return retVal;
}

void ElemRemover::reportState(int count)
{
  std::stringstream ss;
  ss << "cavity_old_" << count;
  ma_dbg::createCavityMesh(adapt, oldEnts, ss.str().c_str());
  ss.str(std::string());
  ss << "new_ents_array_" << count;
  ma_dbg::createCavityMesh(adapt, newEntsArray, ss.str().c_str(), apf::Mesh::TET);
  ss.str(std::string());
  ss << "final_bfaces_" << count;
  showBFaces(adapt, bFaceMap, ss.str().c_str());
}

void ElemRemover::cancel()
{
  destroyNewElements();
  unmark();
  bEdgeMap.clear();
  bFaceMap.clear();
  oldEnts.clear();
  newEnts.clear();
  newEntsArray.setSize(0);
  edgesInQueue.clear();
}

void ElemRemover::transfer()
{
  // keeping as placeholder for whenever element insertion is used
}

void ElemRemover::unmark()
{
  APF_ITERATE(BEdgeMap,bEdgeMap,it) {
    clearFlag(adapt, it->first, MARKED);
  }
  bEdgeMap.clear();

  APF_ITERATE(BFaceMap,bFaceMap,it) {
    clearFlag(adapt, it->first, MARKED);
  }
  bFaceMap.clear();

  APF_ITERATE(EntitySet,oldEnts,it) {
    clearFlag(adapt, *it, CAV_OLD);
    clearFlag(adapt, *it, CAV_NEW);
  }
  oldEnts.clear();

  APF_ITERATE(EntitySet,newEnts,it) {
    clearFlag(adapt, *it, CAV_OLD);
    clearFlag(adapt, *it, CAV_NEW);
  }
  newEnts.clear();
  newEntsArray.setSize(0);

  classifnGroups.clear();
}

void ElemRemover::destroyOldElements()
{
  APF_ITERATE(ma::EntitySet,oldEnts,it)
    if (!getFlag(adapt, *it, CAV_NEW))
      destroyElement(adapt, *it);
    else PCU_ALWAYS_ASSERT((!*it) || adapt->mesh->toModel(*it) != NULL);
  oldEnts.clear();
}

void ElemRemover::destroyNewElements()
{
  APF_ITERATE(ma::EntitySet,newEnts,it)
    if (!getFlag(adapt, *it, CAV_OLD))
      destroyElement(adapt, *it);
    else PCU_ALWAYS_ASSERT((!*it) || adapt->mesh->toModel(*it) != NULL);
  newEnts.clear();
  newEntsArray.setSize(0);
}

}
