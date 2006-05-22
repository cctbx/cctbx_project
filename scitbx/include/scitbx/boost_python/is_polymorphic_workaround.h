#ifndef SCITBX_BOOST_PYTHON_IS_POLYMORPHIC_WORKAROUND_H
#define SCITBX_BOOST_PYTHON_IS_POLYMORPHIC_WORKAROUND_H

#if defined(__APPLE__) && defined(__MACH__) \
 && defined(__GNUC__) && __GNUC__ == 3 && __GNUC_MINOR__ == 3

#include <boost/type_traits/is_polymorphic.hpp>

#define SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(T) \
namespace boost \
{ \
  struct is_polymorphic< \
    T > \
      : mpl::false_ \
  {}; \
}

#else

#define SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(T)

#endif

#endif // SCITBX_BOOST_PYTHON_IS_POLYMORPHIC_WORKAROUND_H
