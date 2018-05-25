#ifndef MA_ELEM_REM_COLLAPSE
#define MA_ELEM_REM_COLLAPSE

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

class ElemRemCollapse
{
  public:
    ElemRemCollapse(Adapt* a);

    bool setCavity(apf::DynamicArray<Entity*> elems);
    bool addElement(Entity* e, bool isOld = false);

    /** Only creates an element with the edge and adjacent
	faces on cavity surface. May return NULL under certain circumstances.*/
    Entity* removeEdge(Entity* e, bool* elemMade);
    /** Removes an element determined by the given face on cavity boundary
        and another cavity boundary face sharing an edge with it. */ 
    Entity* removeFace(Entity* e, bool* elemMade);
    /** Remove the highest dimension entity */
    bool removeElement(Entity* e);

    bool makeNewElements();
    void cancel();
    void transfer();
    void destroyOldElements();

  private:
    bool markEdges(Mesh* m, Entity* face, bool dryRun = false);
    void unmarkEdges(Mesh* m, Entity* face);
    
    Adapt* adapter;
    EntitySet oldEnts;
    EntitySet newEnts;
    BEdgeMap bEdgeMap;
    BFaceMap bFaceMap;
    std::vector<Entity*> edgesInQueue;

    class compareEdgeByCosAngle {
    public:
    compareEdgeByCosAngle(BEdgeMap& bEdgeMap) : _bEdgeMap(&bEdgeMap) {}
      bool operator()(Entity* a, Entity* b)
      {
        // recall that this returns whether a compares _less_ that b, and
        // and std::make_heap returns a _max_ heap
        // if a is not on cavity boundary, pop it out first
        if(_bEdgeMap->count(a) == 0) return false;
        // likewise for b
        if(_bEdgeMap->count(b) == 0) return true;
        return (*_bEdgeMap)[a].cda > (*_bEdgeMap)[b].cda;
      }
    private:
      BEdgeMap* const _bEdgeMap;
    };
};

}
#endif // MA_ELEM_REM_COLLAPSE
