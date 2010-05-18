#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/random/boost_python/random.h>

#include <boost/random/uniform_real.hpp>
#include <scitbx/sparse/random.h>

namespace scitbx { namespace sparse { namespace boost_python {

  template <class ElementDistribution>
  struct matrix_distribution_traits
  {};

  template<typename FloatType, class ElementDistribution>
  struct matrix_distribution
  {
    typedef sparse::matrix_distribution<FloatType, ElementDistribution>  wt;
    typedef typename wt::index_type index_type;

    static std::string name() {
      return matrix_distribution_traits<ElementDistribution>::name();
    }

    static wt *make_1(index_type n_rows, index_type n_cols, double density,
                      ElementDistribution &e)
    {
      return new wt(n_rows, n_cols, density, e);
    }

    static wt *make_2(index_type n_rows, index_type n_cols, index_type n_zeroes,
                      ElementDistribution &e)
    {
      return new wt(n_rows, n_cols, n_zeroes, e);
    }

    static void wrap_specific(boost::python::class_<wt> &klass) {
      using namespace boost::python;
      klass
        .add_property("n_rows", &wt::n_rows)
        .add_property("n_cols", &wt::n_cols)
        .add_property("non_zeroes", &wt::non_zeroes)
        ;
      def("matrix_distribution", make_1,
          return_value_policy<manage_new_object>(),
           ((arg("n_rows"), arg("n_cols"), arg("density"),
             arg("elements"))));
      def("matrix_distribution", make_2,
          return_value_policy<manage_new_object>(),
           ((arg("n_rows"), arg("n_cols"), arg("non_zeroes"),
             arg("elements"))));
    }
  };

  template <>
  struct matrix_distribution_traits<boost::uniform_real<double> >
  {
    static std::string name() { return "matrix_with_uniform_element"; }
  };


  template <class ElementDistribution>
  struct vector_distribution_traits
  {};

  template<typename FloatType, class ElementDistribution>
  struct vector_distribution
  {
    typedef sparse::vector_distribution<FloatType, ElementDistribution>  wt;
    typedef typename wt::index_type index_type;

    static std::string name() {
      return vector_distribution_traits<ElementDistribution>::name();
    }

    static wt *make_1(index_type size, double density,
                      ElementDistribution &e)
    {
      return new wt(size, density, e);
    }

    static wt *make_2(index_type size, index_type n_zeroes,
                      ElementDistribution &e)
    {
      return new wt(size, n_zeroes, e);
    }

    static void wrap_specific(boost::python::class_<wt> &klass) {
      using namespace boost::python;
      klass
        .add_property("size", &wt::size)
        .add_property("non_zeroes", &wt::non_zeroes)
        ;
    def("vector_distribution", make_1,
        return_value_policy<manage_new_object>(),
         ((arg("size"), arg("density"), arg("elements"))));
    def("vector_distribution", make_2,
        return_value_policy<manage_new_object>(),
         ((arg("size"), arg("non_zeroes"), arg("elements"))));
    }
  };

  template <>
  struct vector_distribution_traits<boost::uniform_real<double> >
  {
    static std::string name() { return "vector_with_uniform_element"; }
  };

void wrap_random() {
  using scitbx::random::boost_python::wrap_distribution_and_variate;

  wrap_distribution_and_variate<
    matrix_distribution<double, boost::uniform_real<double> > >();

  wrap_distribution_and_variate<
    vector_distribution<double, boost::uniform_real<double> > >();
}

}}} // scitbx::sparse::boost_python
