// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: Using af::shared (rwgk)
     2001 Dec 21: Using iterator-based fftbx interface (rwgk)
     2001 Nov 03: Created (R.W. Grosse-Kunstleve)
 */

#include <boost/python/cross_module.hpp>
#include <cctbx/array_family/tiny_bpl.h>
#include <cctbx/array_family/small_bpl.h>
#include <cctbx/array_family/shared_bpl_.h>
#include <cctbx/array_family/flex_types.h>
#include <cctbx/fftbx/gridding.h>
#include <cctbx/fftbx/complex_to_complex_3d.h>
#include <cctbx/fftbx/real_to_complex_3d.h>

namespace cctbx { namespace af { namespace bpl { namespace {

  void import_flex()
  {
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(int, "int")
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(double, "double")
    CCTBX_ARRAY_FAMILY_FLEX_IMPORT(std::complex<double>, "complex_double");
  }

}}}} // namespace cctbx::af::bpl<anonymous>

CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(int)
CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(double)

using namespace cctbx;

namespace {

  void raise_size_error()
  {
    PyErr_SetString(PyExc_RuntimeError, "Array is too small.");
    boost::python::throw_error_already_set();
  }

  void assert_0_based_1d_size(
    af::flex_grid<> const& grid,
    std::size_t sz)
  {
    bpl_utils::assert_0_based_1d(grid);
    if (grid.size1d() < sz) raise_size_error();
  }

  void assert_0_based_3d_size(
    af::flex_grid<> const& grid,
    af::int3 const& fft_n)
  {
    bpl_utils::assert_0_based_3d(grid);
    for(std::size_t i=0;i<3;i++) {
      if (grid.grid()[i] != fft_n[i]) raise_size_error();
    }
  }

  int adjust_gridding_2(int min_grid,
                        int max_prime)
  {
    return fftbx::adjust_gridding(min_grid, max_prime);
  }

  int adjust_gridding_3(int min_grid,
                        int max_prime,
                        int mandatory_factor)
  {
    return fftbx::adjust_gridding(min_grid, max_prime, mandatory_factor);
  }

  af::flex_grid_default_index_type
  adjust_gridding_triple_2(
    af::flex_grid_default_index_type const& min_grid,
    int max_prime)
  {
    return fftbx::adjust_gridding_array_flex(min_grid, max_prime);
  }

  af::flex_grid_default_index_type
  adjust_gridding_triple_3(
    af::flex_grid_default_index_type const& min_grid,
    int max_prime,
    af::flex_grid_default_index_type const& mandatory_factors)
  {
    return fftbx::adjust_gridding_array_flex(min_grid, max_prime,
                                             mandatory_factors);
  }

  typedef af::flex_double flex_real_array;
  typedef af::ref<double, af::grid<3> > ref_3d_real_array;
  typedef af::flex_complex_double flex_complex_array;
  typedef af::ref<std::complex<double>, af::grid<3> > ref_3d_complex_array;

  flex_complex_array
  cc_forward_complex(fftbx::complex_to_complex<double>& fft,
                     flex_complex_array a)
  {
    assert_0_based_1d_size(a.accessor(), fft.N());
    fft.forward(a.begin());
    return flex_complex_array(a, af::flex_grid<>(fft.N())
      .set_layout(fft.N()));
  }

  flex_complex_array
  cc_backward_complex(fftbx::complex_to_complex<double>& fft,
                      flex_complex_array a)
  {
    assert_0_based_1d_size(a.accessor(), fft.N());
    fft.backward(a.begin());
    return flex_complex_array(a, af::flex_grid<>((fft.N()))
      .set_layout(fft.N()));
  }

  flex_complex_array
  cc_forward_real(fftbx::complex_to_complex<double>& fft,
                  flex_real_array a)
  {
    assert_0_based_1d_size(a.accessor(), 2 * fft.N());
    fft.forward(a.begin());
    return flex_complex_array(a.handle(), af::flex_grid<>((fft.N()))
      .set_layout(fft.N()));
  }

  flex_complex_array
  cc_backward_real(fftbx::complex_to_complex<double>& fft,
                   flex_real_array a)
  {
    assert_0_based_1d_size(a.accessor(), 2 * fft.N());
    fft.backward(a.begin());
    return flex_complex_array(a.handle(), af::flex_grid<>((fft.N()))
      .set_layout(fft.N()));
  }

  flex_complex_array
  rc_forward_complex(fftbx::real_to_complex<double>& fft,
                     flex_complex_array a)
  {
    assert_0_based_1d_size(a.accessor(), fft.Ncomplex());
    fft.forward(a.begin());
    return flex_complex_array(a, af::flex_grid<>((fft.Ncomplex()))
      .set_layout(fft.Ncomplex()));
  }

  flex_real_array
  rc_backward_complex(fftbx::real_to_complex<double>& fft,
                      flex_complex_array a)
  {
    assert_0_based_1d_size(a.accessor(), fft.Ncomplex());
    fft.backward(a.begin());
    return flex_real_array(a.handle(), af::flex_grid<>((fft.Mreal()))
      .set_layout(fft.Nreal()));
  }

  flex_complex_array
  rc_forward_real(fftbx::real_to_complex<double>& fft,
                  flex_real_array a)
  {
    assert_0_based_1d_size(a.accessor(), fft.Mreal());
    fft.forward(a.begin());
    return flex_complex_array(a.handle(), af::flex_grid<>((fft.Ncomplex()))
      .set_layout(fft.Ncomplex()));
  }

  flex_real_array
  rc_backward_real(fftbx::real_to_complex<double>& fft,
                   flex_real_array a)
  {
    assert_0_based_1d_size(a.accessor(), fft.Mreal());
    fft.backward(a.begin());
    return flex_real_array(a, af::flex_grid<>((fft.Mreal()))
      .set_layout(fft.Nreal()));
  }

  flex_complex_array
  cc_3d_forward_complex(fftbx::complex_to_complex_3d<double>& fft,
                        flex_complex_array a)
  {
    assert_0_based_3d_size(a.accessor(), fft.N());
    ref_3d_complex_array map(a.begin(), af::grid<3>(fft.N()));
    fft.forward(map);
    return flex_complex_array(a, af::flex_grid<>(af::adapt((fft.N())))
      .set_layout(af::adapt(fft.N())));
  }

  flex_complex_array
  cc_3d_backward_complex(fftbx::complex_to_complex_3d<double>& fft,
                         flex_complex_array a)
  {
    assert_0_based_3d_size(a.accessor(), fft.N());
    ref_3d_complex_array map(a.begin(), af::grid<3>(fft.N()));
    fft.backward(map);
    return flex_complex_array(a, af::flex_grid<>(af::adapt((fft.N())))
      .set_layout(af::adapt(fft.N())));
  }

  flex_complex_array
  cc_3d_forward_real(fftbx::complex_to_complex_3d<double>& fft,
                     flex_real_array a)
  {
    assert_0_based_3d_size(a.accessor(), fftbx::Nreal_from_Ncomplex(fft.N()));
    ref_3d_real_array map(
      a.begin(), af::grid<3>(fftbx::Nreal_from_Ncomplex(fft.N())));
    fft.forward(map);
    return flex_complex_array(a.handle(), af::flex_grid<>(af::adapt((fft.N())))
      .set_layout(af::adapt(fft.N())));
  }

  flex_complex_array
  cc_3d_backward_real(fftbx::complex_to_complex_3d<double>& fft,
                      flex_real_array a)
  {
    assert_0_based_3d_size(a.accessor(), fftbx::Nreal_from_Ncomplex(fft.N()));
    ref_3d_real_array map(
      a.begin(), af::grid<3>(fftbx::Nreal_from_Ncomplex(fft.N())));
    fft.backward(map);
    return flex_complex_array(a.handle(), af::flex_grid<>(af::adapt((fft.N())))
      .set_layout(af::adapt(fft.N())));
  }

  flex_complex_array
  rc_3d_forward_complex(fftbx::real_to_complex_3d<double>& fft,
                        flex_complex_array a)
  {
    assert_0_based_3d_size(a.accessor(), fft.Ncomplex());
    ref_3d_complex_array map(a.begin(), af::grid<3>(fft.Ncomplex()));
    fft.forward(map);
    return flex_complex_array(a, af::flex_grid<>(af::adapt((fft.Ncomplex())))
      .set_layout(af::adapt(fft.Ncomplex())));
  }

  flex_real_array
  rc_3d_backward_complex(fftbx::real_to_complex_3d<double>& fft,
                         flex_complex_array a)
  {
    assert_0_based_3d_size(a.accessor(), fft.Ncomplex());
    ref_3d_complex_array map(a.begin(), af::grid<3>(fft.Ncomplex()));
    fft.backward(map);
    return flex_real_array(a.handle(), af::flex_grid<>(af::adapt((fft.Mreal())))
      .set_layout(af::adapt(fft.Nreal())));
  }

  flex_complex_array
  rc_3d_forward_real(fftbx::real_to_complex_3d<double>& fft,
                     flex_real_array a)
  {
    assert_0_based_3d_size(a.accessor(), fft.Mreal());
    ref_3d_real_array map(a.begin(), af::grid<3>(fft.Mreal()));
    fft.forward(map);
    return flex_complex_array(a.handle(),
      af::flex_grid<>(af::adapt((fft.Ncomplex())))
      .set_layout(af::adapt(fft.Ncomplex())));
  }

  flex_real_array
  rc_3d_backward_real(fftbx::real_to_complex_3d<double>& fft,
                      flex_real_array a)
  {
    assert_0_based_3d_size(a.accessor(), fft.Mreal());
    ref_3d_real_array map(a.begin(), af::grid<3>(fft.Mreal()));
    fft.backward(map);
    return flex_real_array(a, af::flex_grid<>(af::adapt((fft.Mreal())))
      .set_layout(af::adapt(fft.Nreal())));
  }

#   include <cctbx/basic/from_bpl_import.h>

  void init_module(python::module_builder& this_module)
  {
    const std::string Revision = "$Revision$";
    this_module.add(ref(to_python(
        Revision.substr(11, Revision.size() - 11 - 2))), "__version__");

    af::bpl::import_flex();

    class_builder<fftbx::factorization>
    py_factorization(this_module, "factorization");
    class_builder<fftbx::complex_to_complex<double> >
    py_complex_to_complex(this_module, "complex_to_complex");
    class_builder<fftbx::real_to_complex<double> >
    py_real_to_complex(this_module, "real_to_complex");
    class_builder<fftbx::complex_to_complex_3d<double> >
    py_complex_to_complex_3d(this_module, "complex_to_complex_3d");
    class_builder<fftbx::real_to_complex_3d<double> >
    py_real_to_complex_3d(this_module, "real_to_complex_3d");

    py_complex_to_complex.declare_base(
      py_factorization, boost::python::without_downcast);
    py_real_to_complex.declare_base(
      py_factorization, boost::python::without_downcast);

    this_module.def(adjust_gridding_2, "adjust_gridding");
    this_module.def(adjust_gridding_3, "adjust_gridding");
    this_module.def(adjust_gridding_triple_2, "adjust_gridding_triple");
    this_module.def(adjust_gridding_triple_3, "adjust_gridding_triple");

    py_factorization.def(constructor<>());
    py_factorization.def(constructor<std::size_t, bool>());
    py_factorization.def(&fftbx::factorization::N, "N");
    py_factorization.def(&fftbx::factorization::Factors, "Factors");

    py_complex_to_complex.def(constructor<>());
    py_complex_to_complex.def(constructor<std::size_t>());
    py_complex_to_complex.def(
      &fftbx::complex_to_complex<double>::WA, "WA");
    py_complex_to_complex.def(cc_forward_complex, "forward");
    py_complex_to_complex.def(cc_backward_complex, "backward");
    py_complex_to_complex.def(cc_forward_real, "forward");
    py_complex_to_complex.def(cc_backward_real, "backward");

    py_real_to_complex.def(constructor<>());
    py_real_to_complex.def(constructor<std::size_t>());
    py_real_to_complex.def(&fftbx::real_to_complex<double>::Nreal, "Nreal");
    py_real_to_complex.def(&fftbx::real_to_complex<double>::Mreal, "Mreal");
    py_real_to_complex.def(
      &fftbx::real_to_complex<double>::Ncomplex, "Ncomplex");
    py_real_to_complex.def(&fftbx::real_to_complex<double>::WA, "WA");
    py_real_to_complex.def(rc_forward_complex, "forward");
    py_real_to_complex.def(rc_backward_complex, "backward");
    py_real_to_complex.def(rc_forward_real, "forward");
    py_real_to_complex.def(rc_backward_real, "backward");

    py_complex_to_complex_3d.def(constructor<>());
    py_complex_to_complex_3d.def(
      constructor<std::size_t, std::size_t, std::size_t>());
    py_complex_to_complex_3d.def(
      constructor<const af::int3&>());
    py_complex_to_complex_3d.def(
      &fftbx::complex_to_complex_3d<double>::N, "N");
    py_complex_to_complex_3d.def(cc_3d_forward_complex, "forward");
    py_complex_to_complex_3d.def(cc_3d_backward_complex, "backward");
    py_complex_to_complex_3d.def(cc_3d_forward_real, "forward");
    py_complex_to_complex_3d.def(cc_3d_backward_real, "backward");

    py_real_to_complex_3d.def(constructor<>());
    py_real_to_complex_3d.def(
      constructor<std::size_t, std::size_t, std::size_t>());
    py_real_to_complex_3d.def(
      constructor<const af::int3&>());
    py_real_to_complex_3d.def(
      &fftbx::real_to_complex_3d<double>::Nreal, "Nreal");
    py_real_to_complex_3d.def(
      &fftbx::real_to_complex_3d<double>::Mreal, "Mreal");
    py_real_to_complex_3d.def(
      &fftbx::real_to_complex_3d<double>::Ncomplex, "Ncomplex");
    py_real_to_complex_3d.def(rc_3d_forward_complex, "forward");
    py_real_to_complex_3d.def(rc_3d_backward_complex, "backward");
    py_real_to_complex_3d.def(rc_3d_forward_real, "forward");
    py_real_to_complex_3d.def(rc_3d_backward_real, "backward");
  }

}

BOOST_PYTHON_MODULE_INIT(fftbx)
{
  boost::python::module_builder this_module("fftbx");
  init_module(this_module);
}
