#ifndef LIBDISTL_TYPES_H
#define LIBDISTL_TYPES_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/constants.h>

namespace Distl {

template<class T>
struct list_types
{
  typedef scitbx::af::shared<T> list_t;
};

}
#endif
