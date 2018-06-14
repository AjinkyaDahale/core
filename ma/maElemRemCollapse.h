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

    std::map< Model*, EntitySet> classifnGroups;
  
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
};

}
#endif // MA_ELEM_REM_COLLAPSE
