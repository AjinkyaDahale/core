#ifndef MA_ELEM_REM_COLLAPSE
#define MA_ELEM_REM_COLLAPSE

#include "maAdapt.h"
#include <pcu_util.h>

#include <queue>

namespace ma {

typedef std::pair<Entity*, double> BEdges;

class ElemRemCollapse
{
  public:
    ElemRemCollapse(Adapt* a);

    bool setCavity(apf::DynamicArray<Entity*> elems);
    bool addElement(Entity* e);

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
    apf::DynamicArray<BEdges> bEdges;
};

}
#endif // MA_ELEM_REM_COLLAPSE
