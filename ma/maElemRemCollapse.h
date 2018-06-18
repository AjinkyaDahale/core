#ifndef MA_ELEM_REM_COLLAPSE
#define MA_ELEM_REM_COLLAPSE

#include "maAdapt.h"
#include "maCollapse.h"
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
 class ElemRemCollapse : public Collapse
{
  public:
    virtual void Init(Adapt* a);

    bool setCavity(apf::DynamicArray<Entity*>& elems);
    bool addElement(Entity* e, bool isOld = false);
  void setIgnoredModelFaces(EntityArray& ups);
  void addClassifnGroups(EntityArray& ups);

    /** Only creates an element with the edge and adjacent
	faces on cavity surface. May return NULL under certain circumstances.*/
    Entity* removeEdge(Entity* e, bool* elemMade);
    /** Removes an element determined by the given face on cavity boundary
        and another cavity boundary face sharing an edge with it. */
    Entity* removeFace(Entity* e, bool* elemMade);
    /** Remove the highest dimension entity */
    bool removeElement(Entity* e);

    bool makeNewElements(double qualityToBeat);
    void cancel(bool cavOnly);
    virtual void cancel();
    void transfer();
    virtual void destroyOldElements();
    virtual void destroyNewElements();
    virtual bool tryThisDirectionNoCancel(double qualityToBeat);
    void unmark(bool cavOnly = false);

  private:
    bool markEdges(Mesh* m, Entity* face, bool dryRun = false);
    void unmarkEdges(Mesh* m, Entity* face);
  bool edgeIntersectsFace(Adapt* a, Entity* face, Entity* edge);

  Entity* buildOrFind(Adapt* a, int type, Entity** vs,
                      bool* elemMade, bool* elemInverted=NULL);

    ClassifnGroups classifnGroups;
  
    bool newTetClear(Adapt* a, Entity* tet);

    double qualToBeat;
    /* Adapt* adapt; */
    apf::ModelEntity* modelEnt;
    
    EntitySet oldEnts;
    EntitySet newEnts;
    EntityArray newEntsArray;
    BEdgeMap bEdgeMap;
    BFaceMap bFaceMap;
    std::vector<Entity*> edgesInQueue;

  std::set<Model*> ignoredFaceClassifn;

    class compareEdgeByCosAngle {
    public:
    compareEdgeByCosAngle(BEdgeMap& bEdgeMap) : _bEdgeMap(&bEdgeMap) {}
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
    ElementBuilder2(
        Mesh* m,
        int baseType,
        ClassifnGroups& classifnGroups,
        apf::BuildCallback* cb,
        bool* elemMade) :
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
#endif // MA_ELEM_REM_COLLAPSE
