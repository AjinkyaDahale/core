#ifndef MA_ELEM_REM_COLLAPSE
#define MA_ELEM_REM_COLLAPSE

#include "maAdapt.h"
#include <pcu_util.h>

#include <queue>
#include <map>

namespace ma {

typedef std::pair<Entity*, double> BEdge;

class BEdge1
{
public:
  BEdge1()
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

typedef std::map<Entity*, BEdge1> EdgeMap;

class ElemRemCollapse
{
  public:
    ElemRemCollapse(Adapt* a);

    bool setCavity(apf::DynamicArray<Entity*> elems);
    bool addElement(Entity* e);

    bool removeEdge(Entity* e);

    bool makeNewElements();
    void cancel();
    void transfer();
    void destroyOldElements();
  private:
    Adapt* adapter;
    Entity* inPoint;
    EntityArray cavityEnts;
    EntityArray boundaryEnts;
    apf::DynamicArray<bool> cavEntPositive;
    apf::DynamicArray<EntityArray> newEnts;
    EdgeMap edgeMap;
};

}
#endif // MA_ELEM_REM_COLLAPSE
