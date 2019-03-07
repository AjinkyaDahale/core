#ifndef MA_ELEM_REM_COLLAPSE_H
#define MA_ELEM_REM_COLLAPSE_H

#include "maAdapt.h"
#include "maCollapse.h"
#include "maElemRemover.h"
#include <pcu_util.h>

namespace ma {

class ElemRemCollapse : public Collapse
{
public:
  virtual void Init(Adapt* a);
  void cancel(bool cavOnly);
  virtual void cancel();
  void transfer();
  virtual void destroyOldElements();
  virtual void destroyNewElements();
  virtual bool tryThisDirectionNoCancel(double qualityToBeat);
  void unmark(bool cavOnly = false);

private:
  ElemRemover elemRemover;
};

}
#endif // MA_ELEM_REM_COLLAPSE_H
