#ifndef MA_ELEM_REMOVER_H
#define MA_ELEM_REMOVER_H

#include "maAdapt.h"
#include <pcu_util.h>

#include <queue>
#include <map>

namespace ma {

class BEdge
{
public:
  BEdge()
  {
    edge = NULL;
    face1 = NULL;
    face2 = NULL;
    cda = 1.0;
  }
  Entity* edge;
  Entity* face1;
  Entity* face2;
  double cda;
};

typedef std::map<Entity*, BEdge> BEdgeMap;
// The values correspond to <reference element, is ref elem positive>
typedef std::map<Entity*, std::pair<Entity*, bool> > BFaceMap;
typedef std::map<Model*, EntitySet> ClassifnGroups;
class ElemRemover
{
public:
  void Init(Adapt* a);
  bool setCavity(apf::DynamicArray<Entity*>& elems);
  /** Add a single element to the cavity.
   * \param e: the entitiy to be added.
   * \param isOld: `true` if the element is from the unmodified mesh.
   * If the operation cancels, do not destroy this.
   */
  bool addElement(Entity* e, bool isOld = false);
  void setIgnoredModelFaces(EntityArray& ups);
  void addClassifnGroups(EntityArray& ups);

  /** Only creates an element with the edge and adjacent
      faces on cavity surface. May return NULL under certain circumstances.*/
  Entity* removeEdge(Entity* e, bool* elemMade);
  /** Remove the highest dimension entity */
  bool removeElement(Entity* e);

  void reportState();

  bool makeNewElements(double qualityToBeat);
  void cancel();
  void transfer();
  void destroyOldElements();
  void destroyNewElements();
  void unmark();

private:
  bool markEdges(Mesh* m, Entity* face, bool dryRun = false);
  void unmarkEdges(Mesh* m, Entity* face);

  Entity* buildOrFind(Adapt* a, int type, Entity** vs,
                      bool* elemMade, bool* elemInverted = NULL);

  bool edgeIntersectsFace(Adapt* a, Entity* face, Entity* edge);
  bool newTetClear(Adapt* a, Entity* tet);

  Adapt* adapt;

  ClassifnGroups classifnGroups;

  double qualToBeat;
  apf::ModelEntity* modelEnt;

  EntitySet oldEnts;
  EntitySet newEnts;
  EntityArray newEntsArray;
  BEdgeMap bEdgeMap;
  BFaceMap bFaceMap;
  std::vector<Entity*> edgesInQueue;

  std::set<Model*> ignoredFaceClassifn;

  class CompareEdgeByCosAngle {
  public:
    CompareEdgeByCosAngle(BEdgeMap& bEdgeMap) : _bEdgeMap(&bEdgeMap) {}
    bool operator()(Entity* a, Entity* b)
    {
      // Recall that this returns whether a compares _less_ that b, and
      // and std::make_heap returns a _max_ heap.

      // If a is not on cavity boundary, pop it out first
      if(_bEdgeMap->count(a) == 0) return false;
      // likewise for b
      if(_bEdgeMap->count(b) == 0) return true;

      return (*_bEdgeMap)[a].cda < (*_bEdgeMap)[b].cda;
    }
  private:
    BEdgeMap* const _bEdgeMap;
  };

  class ElementBuilder2 : public apf::ElementVertOp
  {
  public:
    ElementBuilder2(Mesh* m, int baseType, ClassifnGroups& classifnGroups,
      apf::BuildCallback* cb, bool* elemMade) :
      _classifnGroups(classifnGroups)
    {
      mesh = m;
      _elemMade = elemMade;
      _baseType = baseType;
      callback = cb;
    }
    virtual Entity* apply(int type, Entity** down)
    {
      Model* c = NULL;
      int c_dim = 4;
      int min_dim = apf::Mesh::typeDimension[type];
      int nv = apf::Mesh::adjacentCount[type][0];
      for (int i = 0; i < nv; ++i) {
        int dn_dim = mesh->getModelType(mesh->toModel(down[i]));
        min_dim = (min_dim < dn_dim) ? dn_dim : min_dim;
      }
      // Find the lowest group containing all
      APF_ITERATE(ClassifnGroups,_classifnGroups,it) {
        int nc_dim = mesh->getModelType(it->first);
        if (nc_dim < c_dim &&
            nc_dim >= min_dim){
          bool allIn = true;
          for (int i = 0; allIn && i < nv; ++i)
            allIn = allIn && it->second.count(down[i]);
          if (allIn){
            c = it->first;
            c_dim = nc_dim;
          }
        }
      }
      PCU_ALWAYS_ASSERT(c != NULL);
      Entity* ans = NULL;
      if (type == _baseType)
        ans = makeOrFind(mesh,c,type,down,callback,_elemMade);
      else
        ans = makeOrFind(mesh,c,type,down,callback);
      APF_ITERATE(ClassifnGroups,_classifnGroups,it) {
        int nc_dim = mesh->getModelType(it->first);
        bool allIn = true;
        for (int i = 0; allIn && i < nv; ++i)
          allIn = allIn && it->second.count(down[i]);
        if (allIn &&
            nc_dim >= min_dim) {
          _classifnGroups[it->first].insert(ans);
        }
      }
      return ans;
    }
  private:
    Mesh* mesh;
    ClassifnGroups& _classifnGroups;
    bool* _elemMade;
    int _baseType;
    apf::BuildCallback* callback;
  };

};

}
#endif // MA_ELEM_REMOVER_H
