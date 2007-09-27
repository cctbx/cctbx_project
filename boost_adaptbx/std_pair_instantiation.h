#ifndef BOOST_ADAPTBX_STD_PAIR_FWD_H
#define BOOST_ADAPTBX_STD_PAIR_FWD_H

/* This header has a dual purpose:
- it is to be included by std_pair_exp.cpp, which then create an instance of
"instantiate" to actually create the conversions
- it is to be included by any Boost.Python code which wishes to use those
conversions with the same intended purpose as
scitbx/array_family/boost_python/flex_fwd.h
(the reader is referred to the comments in that header)
*/

#include <boost_adaptbx/std_pair_conversion.h>
#include <boost/optional.hpp>

namespace boost_adaptbx { namespace std_pair_conversions {

  struct instantiate {
    to_python<int, double> a;

    // next 4 needed for flex.find_partial_sum_xxx
    to_python< boost::optional<std::size_t>, boost::optional<int>      > b;
    to_python< boost::optional<std::size_t>, boost::optional<unsigned> > c;
    to_python< boost::optional<std::size_t>, boost::optional<float>    > d;
    to_python< boost::optional<std::size_t>, boost::optional<double>   > e;
  };

}}

#endif
