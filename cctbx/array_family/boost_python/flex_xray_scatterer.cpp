/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <cctbx/xray/scatterer.h>
#include <boost/python/return_internal_reference.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  struct to_string : pickle_double_buffered::to_string
  {
    using pickle_double_buffered::to_string::operator<<;

    to_string& operator<<(cctbx::xray::scatterer<> const& val)
    {
      *this << val.label
            << val.caasf.label()
            << val.fp_fdp;
      *this << val.site.const_ref()
            << val.occupancy
            << val.anisotropic_flag
            << val.u_iso;
      *this << val.u_star.const_ref()
            << val.multiplicity()
            << val.weight();
      return *this;
    }
  };

  struct from_string : pickle_double_buffered::from_string
  {
    from_string(PyObject* str_obj)
    : pickle_double_buffered::from_string(str_obj)
    {}

    using pickle_double_buffered::from_string::operator>>;

    from_string& operator>>(cctbx::xray::scatterer<>& val)
    {
      std::string caasf_label;
      int multiplicity;
      cctbx::xray::scatterer<>::float_type weight;
      *this >> val.label
            >> caasf_label
            >> val.fp_fdp;
      *this >> val.site.ref()
            >> val.occupancy
            >> val.anisotropic_flag
            >> val.u_iso;
      *this >> val.u_star.ref()
            >> multiplicity
            >> weight;
      val.setstate(caasf_label, multiplicity, weight);
      return *this;
    }
  };

} // namespace <anonymous>

  void wrap_flex_xray_scatterer()
  {
    using namespace cctbx;

    flex_wrapper<cctbx::xray::scatterer<>,
                 boost::python::return_internal_reference<>
                >::plain("xray_scatterer")
      .def_pickle(flex_pickle_double_buffered<
        cctbx::xray::scatterer<>, to_string, from_string>());
  }

}}} // namespace scitbx::af::boost_python
