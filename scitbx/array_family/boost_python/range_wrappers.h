#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_RANGE_WRAPPERS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_RANGE_WRAPPERS_H

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <scitbx/array_family/range.h>

namespace scitbx { namespace af { namespace boost_python {

  template <typename ValueType,
            typename IntType,
            template<typename> class CheckType=range_args::no_check>
  struct range_wrappers
  {
    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      def(python_name,
        (shared<ValueType>(*)(IntType const&, IntType const&, IntType const&))
          range<ValueType, IntType, CheckType>::array, (
            arg_("start"), arg_("stop"), arg_("step")));
      def(python_name,
        (shared<ValueType>(*)(IntType const&, IntType const&))
          range<ValueType, IntType, CheckType>::array, (
            arg_("start"), arg_("stop")));
      def(python_name,
        (shared<ValueType>(*)(IntType const&))
          range<ValueType, IntType, CheckType>::array, (
            arg_("stop")));
    }
  };

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_RANGE_WRAPPERS_H
