#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/def.hpp>

#include <scitbx/math/distributions.h>

#if defined(__GNUC__) && __GNUC__ == 3 && __GNUC_MINOR__ == 2
# define SCITBX_MATH_STUDENTS_T_DISABLED // to avoid compilation errors
#endif

namespace scitbx { namespace math {

/*! Wrappers for boost::math statistical distributions.
    See also:
    http://www.boost.org/libs/math/doc/sf_and_dist/html/math_toolkit/dist.html
 */
namespace {

  template <typename FloatType>
  struct normal_distribution_wrappers
  {
    typedef boost::math::normal_distribution<FloatType> wt;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<wt>("normal_distribution", no_init)
        .def(init<FloatType, FloatType>((arg("mean")=0, arg("sd")=1)))
      ;
    }
  };

#ifndef SCITBX_MATH_STUDENTS_T_DISABLED
  template <typename FloatType>
  struct students_t_distribution_wrappers
  {
    typedef boost::math::students_t_distribution<FloatType> wt;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<wt>("students_t_distribution", no_init)
        .def(init<FloatType>(arg("v")))
        .def("degrees_of_freedom", &wt::degrees_of_freedom)
        .def("find_degrees_of_freedom", &wt::find_degrees_of_freedom)
      ;
    }
  };
#endif

  template <typename FloatType, class Distribution>
  struct non_member_function_wrappers
  {
    typedef Distribution wt;

    static void
    wrap()
    {
      using namespace boost::python;
      def("mean", (FloatType(*)(wt const&)) boost::math::mean);
      def("median", (FloatType(*)(wt const&)) boost::math::median);
      def("mode", (FloatType(*)(wt const&)) boost::math::mode);
      def("variance", (FloatType(*)(wt const&)) boost::math::variance);
      def("standard_deviation",
        (FloatType(*)(wt const&)) boost::math::standard_deviation);
      def("skewness", (FloatType(*)(wt const&)) boost::math::skewness);
      def("kurtosis", (FloatType(*)(wt const&)) boost::math::kurtosis);
      def("pdf", (FloatType(*)(wt const&, FloatType const&)) boost::math::pdf);
      def("cdf", (FloatType(*)(wt const&, FloatType const&)) boost::math::cdf);
      def("quantile", (FloatType(*)(wt const&, FloatType const&))
        boost::math::quantile);
/* this is changed due to the compilation error with MSVC 2005 x64 form:
	  def("quantiles",(
        scitbx::af::shared<FloatType>(*)(wt const&, std::size_t))
		quantiles);
*/
	  def("quantiles", quantiles<FloatType,wt>);
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_distributions()
  {
    normal_distribution_wrappers<double>::wrap();
    non_member_function_wrappers<
      double, boost::math::normal_distribution<double> >::wrap();
#ifndef SCITBX_MATH_STUDENTS_T_DISABLED
    students_t_distribution_wrappers<double>::wrap();
    non_member_function_wrappers<
      double, boost::math::students_t_distribution<double> >::wrap();
#endif
  }

}}} // namespace scitbx::math::boost_python
