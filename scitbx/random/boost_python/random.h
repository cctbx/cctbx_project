#ifndef SCITBX_RANDOM_BOOST_PYTHON_RANDOM_H
#define SCITBX_RANDOM_BOOST_PYTHON_RANDOM_H

#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/random/variate_generator.h>
#include <scitbx/random/mersenne_twister.h>
#include <boost/optional.hpp>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/return_value_policy.hpp>


namespace scitbx { namespace random { namespace boost_python {

  /// Boost.Python wrapper for variate_generator<Engine, Distribution>
  /** The Python binding is a proper iterator with 'next' and '__iter__'.
      It has no constructor: instead the factory function
      variate(engine, distribution), also created, here shall be used.
   */
  template<class Engine, class Distribution>
  struct variate_generator_wrappers
  {
    typedef scitbx::random::variate_generator<Engine, Distribution> wt;
    typedef typename wt::result_type result_type;
    typedef boost::optional<std::size_t> opt_size_t;

    static boost::python::object
    generate_one_or_many(wt &self, opt_size_t size) {
      using namespace boost::python;
      if (!size) return object(self());
      else return object(self(*size));
    }

    static result_type generate_one(wt &self) {
      return self();
    }

    static wt &identity(wt &self) {
      return self;
    }

    static wt *make(Engine e, Distribution d) {
      return new wt(e, d);
    }

    static void
    wrap(char const *name)
    {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def("__call__", generate_one_or_many, arg("size")=opt_size_t())
        .def("next", generate_one)
        .def("__iter__", identity, return_self<>())
        .def("__next__", generate_one)
        ;
      def("variate", make, return_value_policy<manage_new_object>(),
          (arg("engine"), arg("distribution")));
    }
  };


  /// Boost.Python wrapper for
  template <class Delegate>
  struct wrap_distribution_and_variate
  {
    typedef typename Delegate::wt wt;

    wrap_distribution_and_variate() {
      using namespace boost::python;
      std::string class_name = Delegate::name() + "_distribution";
      class_<wt> klass(class_name.c_str(), no_init);
      klass.def("reset", &wt::reset);
      Delegate::wrap_specific(klass);

      std::string name = Delegate::name() + std::string("_variate_generator");
      variate_generator_wrappers<
      boost_random::mt19937 &,
      typename Delegate::wt>::wrap(name.c_str());
    }
  };



}}}

#endif // GUARD
