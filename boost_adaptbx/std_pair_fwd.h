#ifndef BOOST_ADAPTBX_STD_PAIR_FWD_H
#define BOOST_ADAPTBX_STD_PAIR_FWD_H

#include <utility>
#include <boost/optional.hpp>

namespace boost_adaptbx { namespace std_pair_conversions {
  namespace boost_python {

  struct std_pair_fwd
  {
    friend void f(std::pair< boost::optional<std::size_t>,
                             boost::optional<int>       >);
    friend void f(std::pair< boost::optional<std::size_t>,
                             boost::optional<unsigned>       >);
    friend void f(std::pair< boost::optional<std::size_t>,
                             boost::optional<float>       >);
    friend void f(std::pair< boost::optional<std::size_t>,
                             boost::optional<double>       >);
  };

}}}

#endif
