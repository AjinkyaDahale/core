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
typedef std::map<Entity*, std::pair<Entity*, bool> > BFaceMap;

class ElemRemCollapse
{
  public:
    ElemRemCollapse(Adapt* a);

    bool setCavity(apf::DynamicArray<Entity*> elems);
    bool addElement(Entity* e);

    Entity* removeEdge(Entity* e);

    bool makeNewElements();
    void cancel();
    void transfer();
    void destroyOldElements();
  private:
    Adapt* adapter;
    Entity* inPoint;
    EntityArray oldEnts;
    BFaceMap bFaceMap;
    EntityArray newEnts;
    BEdgeMap bEdgeMap;
};

}
#endif // MA_ELEM_REM_COLLAPSE
