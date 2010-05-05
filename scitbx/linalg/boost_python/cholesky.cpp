#include <boost/python/class.hpp>
#include <boost/python/def.hpp>

#include <scitbx/matrix/cholesky.h>


namespace scitbx { namespace matrix { namespace boost_python {

  struct matrix_cholesky_gill_murray_wright_decomposition_in_place_wrappers
  {
    typedef matrix::cholesky::gill_murray_wright_decomposition_in_place<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("gill_murray_wright_cholesky_decomposition_in_place",
                  no_init)
        .def(init<af::shared<double> const &, optional<double> >(
             (arg("packed_u"), arg("epsilon"))))
        .def_readonly("epsilon", &w_t::epsilon)
        .add_property("packed_u", make_getter(&w_t::packed_u, rbv()))
        .add_property("e", make_getter(&w_t::e, rbv()))
        .add_property("pivots", make_getter(&w_t::pivots, rbv()))
        .def("solve", &w_t::solve, (arg("b")))
      ;
    }
  };

  struct cholesky_failure_info_wrapper
  {
    typedef cholesky::failure_info<double> wt;

    static bool nonzero(wt const &self) {
      return self.failed;
    }

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def_readonly("index", &wt::index)
        .def_readonly("value", &wt::value)
        .def("__nonzero__", nonzero)
        ;
    }
  };

  template <class DecompositionType>
  struct cholesky_decomposition_for_python : DecompositionType
  {
    typedef typename DecompositionType::scalar_t scalar_t;
    typedef typename DecompositionType::accessor_type accessor_type;
    typedef af::ref<scalar_t, accessor_type> packed_ref_t;

    af::shared<scalar_t> a_;

    cholesky_decomposition_for_python(af::shared<scalar_t> const &a)
      : DecompositionType(packed_ref_t(const_cast<scalar_t *>(a.begin()),
                                       af::dimension_from_packed_size(a.size()))),
        a_(a) // so as to make a lives longer than this
    {}
  };

  template <class DecompositionType>
  struct cholesky_decomposition_wrapper
  {
    typedef cholesky_decomposition_for_python<DecompositionType> wt;
    typedef typename wt::scalar_t scalar_t;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<af::shared<scalar_t> const &>())
        .def_readonly("failure", &wt::failure)
        .def("solve", &wt::solve)
        ;
    }
  };

  void wrap_cholesky() {
    using namespace boost::python;

    matrix_cholesky_gill_murray_wright_decomposition_in_place_wrappers::wrap();

    cholesky_decomposition_wrapper<
      cholesky::l_l_transpose_decomposition_in_place<double>
    >::wrap("l_l_transpose_cholesky_decomposition_in_place");
    cholesky_decomposition_wrapper<
      cholesky::u_transpose_u_decomposition_in_place<double>
    >::wrap("u_transpose_u_cholesky_decomposition_in_place");

    cholesky_failure_info_wrapper::wrap("cholesky_failure_info");
  }
}}}
