/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/ref.h>
#include <scitbx/boost_python/pickle_double_buffered.h>

namespace scitbx { namespace af { namespace boost_python {
namespace pickle_double_buffered {

  struct to_string
  : scitbx::boost_python::pickle_double_buffered::to_string
  {
    using
      scitbx::boost_python::pickle_double_buffered::to_string::operator<<;

    template <typename ElementType,
              typename AccessorType>
    to_string& operator<<(af::const_ref<ElementType, AccessorType> const& val)
    {
      for(std::size_t i=0;i<val.size();i++) {
        *this << val[i];
      }
      return *this;
    }
  };

  struct from_string
  : scitbx::boost_python::pickle_double_buffered::from_string
  {
    from_string(PyObject* str_obj)
    : scitbx::boost_python::pickle_double_buffered::from_string(str_obj)
    {}

    using
      scitbx::boost_python::pickle_double_buffered::from_string::operator>>;

    template <typename ElementType>
    from_string& operator>>(ref<ElementType> val)
    {
      for(std::size_t i=0;i<val.size();i++) {
        *this >> val[i];
      }
      return *this;
    }
  };

}}}} // namespace scitbx::af::boost_python::pickle_double_buffered
