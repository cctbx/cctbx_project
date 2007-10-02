#ifndef BOOST_ADAPTBX_OPTIONAL_FWD_H
#define BOOST_ADAPTBX_OPTIONAL_FWD_H

#include <boost/optional.hpp>

namespace boost_adaptbx { namespace optional_conversions {
  namespace boost_python {

  struct optional_fwd
  {
    friend void f(boost::optional<double>);
  };

  }}}

#endif
