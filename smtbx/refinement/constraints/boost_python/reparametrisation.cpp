#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/pure_virtual.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/docstring_options.hpp>

#include <boost/operators.hpp>
#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

  /*** To support polymorphism fully, some extra wrappers must be written
   *** (c.f. Boost.Python tutorial)
   ***/

  /* To avoid the traps of multiple inheritance armed by Boost.Python
     ways to deal with polymorphic calls, we must resort to macros
   */

  #define SMTBX_CONSTRAINTS_OVERRIDE_SIZE                                      \
  std::size_t size() const {                                                   \
    return this->get_override("size")();                                       \
  }                                                                            \

  #define SMTBX_CONSTRAINTS_OVERRIDE_LINEARISE                                 \
  void linearise(uctbx::unit_cell const &unit_cell,                            \
                         sparse_matrix_type *jacobian_transpose)               \
  {                                                                            \
    this->get_override("linearise")(unit_cell, jacobian_transpose);            \
  }

  #define SMTBX_CONSTRAINTS_OVERRIDE_STORE                                     \
  virtual void store(uctbx::unit_cell const &unit_cell) const {                \
    this->get_override("store")(unit_cell);                                    \
  }

  #define SMTBX_CONSTRAINTS_BEFRIEND_INIT_PARAM                                \
    template<class ParameterType>                                              \
    friend void initialise_parameter(ParameterType *p,                         \
                                     boost::python::tuple const &arguments,    \
                                     int first);

  #define SMTBX_CONSTRAINTS_PARAMETER_CONSTRUCTOR(klass)                       \
    SMTBX_CONSTRAINTS_BEFRIEND_INIT_PARAM                                      \
    py_##klass(boost::python::tuple const &arguments)                          \
    : klass(boost::python::len(arguments))                                     \
    {                                                                          \
      initialise_parameter(this, arguments);                                   \
    }

  template <class ParameterType>
  void initialise_parameter(ParameterType *p,
                            boost::python::tuple const &arguments,
                            int first = 0)
  {
    using namespace boost::python;
    for (int i=first; i<p->n_arguments(); ++i) {
      extract<parameter *> proxy(arguments[i]);
      if (!proxy.check()) {
        PyErr_Format(PyExc_ValueError,
                     "Argument #%i is not of type "
                     "'smtbx.refinement.constraints.parameter'",
                     i);
        throw_error_already_set();
      }
      p->set_argument(i, proxy());
    }
  }

  struct py_parameter : parameter,
                        boost::python::wrapper<parameter>
  {
    SMTBX_CONSTRAINTS_PARAMETER_CONSTRUCTOR(parameter)

    SMTBX_CONSTRAINTS_OVERRIDE_SIZE
    SMTBX_CONSTRAINTS_OVERRIDE_LINEARISE
  };

  struct py_crystallographic_parameter
  : public crystallographic_parameter,
    boost::python::wrapper<crystallographic_parameter>
  {
    SMTBX_CONSTRAINTS_PARAMETER_CONSTRUCTOR(crystallographic_parameter)

    SMTBX_CONSTRAINTS_OVERRIDE_SIZE
    SMTBX_CONSTRAINTS_OVERRIDE_LINEARISE
    SMTBX_CONSTRAINTS_OVERRIDE_STORE
  };

  struct py_site_parameter : site_parameter,
                             boost::python::wrapper<site_parameter>
  {

    SMTBX_CONSTRAINTS_OVERRIDE_LINEARISE

    static std::size_t checked(std::size_t n_arguments) {
      SMTBX_ASSERT(n_arguments >= 1)(n_arguments);
      return n_arguments;
    }

    SMTBX_CONSTRAINTS_BEFRIEND_INIT_PARAM

    py_site_parameter(boost::python::tuple const &arguments)
    : site_parameter(boost::python::extract<scatterer_type *>(arguments[0])(),
                     checked(boost::python::len(arguments)))
    {
      using namespace boost::python;
      initialise_parameter(this, arguments, 1);
    }

  };


  /*** The usual wrappers
   ***
   *** It should be noted that all non-virtual class wrappers needs
   *** (a) to have an auto_ptr<> storage, and
   *** (b) to instruct Boost.Python about the implicit conversion
   ***     from auto_ptr<wrapped_type> to auto_ptr<parameter>.
   ***/

  struct parameter_wrapper
  {
    typedef parameter wt;
    typedef py_parameter pywt;

    static void wrap() {
      using namespace boost::python;
      return_internal_reference<> rir;
      class_<pywt, boost::noncopyable>("parameter", no_init)
        .def(init<boost::python::tuple>(arg("arguments")))
        .add_property("n_arguments", &wt::n_arguments)
        .def("argument", &wt::argument, rir)
        .add_property("index", &wt::index)
        .add_property("is_root", &wt::is_root)
        .add_property("is_variable", &wt::is_variable)
        .def("size", pure_virtual(&wt::size))
        .def("evaluate", &wt::evaluate)
        .def("linearise", pure_virtual(&wt::linearise))
        ;
    }
  };


  struct crystallographic_parameter_wrapper
  {
    typedef crystallographic_parameter wt;
    typedef py_crystallographic_parameter pywt;

    static void wrap() {
      using namespace boost::python;
      class_<pywt, bases<parameter>,
             boost::noncopyable>("crystallographic_parameter", no_init)
        .def(init<boost::python::tuple>(arg("arguments")))
        .def("store", pure_virtual(&wt::store), arg("unit_cell"))
        ;
    }
  };

  struct independent_scalar_parameter_wrapper
  {
    typedef independent_scalar_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt, bases<parameter>,
             std::auto_ptr<wt>,
             boost::noncopyable>("independent_scalar_parameter", no_init)
        .def(init<double, bool>((arg("value"), arg("variable")=true)))
        .def_readwrite("value", &wt::value)
        .def("size", &wt::size)
        .def("linearise", &wt::linearise, (arg("unit_cell"),
                                           arg("jacobian_transpose")))
        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct site_parameter_wrapper
  {
    typedef site_parameter wt;
    typedef py_site_parameter pywt;

    static void wrap() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      class_<pywt, bases<crystallographic_parameter>,
             boost::noncopyable>("site_parameter", no_init)
        .def(init<boost::python::tuple>(arg("scatterer")))
        .add_property("value",
                      make_getter(&wt::value, rbv), make_setter(&wt::value))
        .def("size", &wt::size)
        .def("store", &wt::store, arg("unit_cell"))
        ;
    }
  };

  struct independent_site_parameter_wrapper
  {
    typedef independent_site_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt, bases<site_parameter>,
             std::auto_ptr<wt>,
             boost::noncopyable>("independent_site_parameter", no_init)
        .def(init<wt::scatterer_type *>(arg("scatterer")))
        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };


  struct reparametrisation_wrapper
  {
    typedef reparametrisation wt;

    // Tricky code ahead!
    static parameter *add(boost::python::tuple args,
                          boost::python::dict kwds)
    {
      using namespace boost::python;
      // We have been called as add(self, klass, *args, **kwds)
      object py_self = args[0], klass = args[1], param_args = args.slice(2, _);
      wt &self = extract<wt &>(py_self)();

      /* Construct a new parameter, taking advantage of Python dynamic
         nature to call the constructor of the proper class through its wrapper
       */
      object py_param = klass(*param_args, **kwds);

      /* All Python wrapper in the parameter hierarchy have std::auto_ptr<>
         storage and instruct Boost.Python of the conversion
         from std::auto_ptr<> to std::auto_ptr<parameter> (downcast to base
         class): thus the following 'extract' works.
       */
      std::auto_ptr<parameter>
      borrowed_param = extract<std::auto_ptr<parameter> >(py_param)();

      /* Apply the official recipe to transfer ownership from Boost.Python
         hands to the reparametrisation instance self
         (http://www.boost.org/doc/libs/1_43_0/libs/python/doc/v2/faq.html#ownership).
       */
      parameter *p = borrowed_param.get();
      self.add(p);
      borrowed_param.release();
      return p;
    }

    static void wrap() {
      using namespace boost::python;
      class_<wt> klass("reparametrisation", no_init);
      klass
        .def(init<uctbx::unit_cell const &>(arg("unit_cell")))
        .def("finalise", &wt::finalise)
        .def("linearise", &wt::linearise)
        .def("store", &wt::store)
        ;
      docstring_options no_signature(true, false);
      klass
        .def("add",
             raw_function(make_function(add, return_internal_reference<>())),
             "r.add(klass, *args, **kwds): create the parameter "
             "klass(*args, **kwds) and add it to the reparametrisation r. "
             "'klass' is therefore typically one of the class in the "
             "parameter hierarchy, e.g.\n"
             "r.add(independent_scalar_parameter, value=1., variable=True)\n"
             "but it may also be a factory function.")
        ;
    }
  };

  void debug(char const *msg, parameter const *p) {
    std::cout << msg << p << std::endl;
  }

  void wrap_reparametrisation() {
    using namespace boost::python;
    parameter_wrapper::wrap();
    independent_scalar_parameter_wrapper::wrap();
    crystallographic_parameter_wrapper::wrap();
    site_parameter_wrapper::wrap();
    independent_site_parameter_wrapper::wrap();
    reparametrisation_wrapper::wrap();
    def("debug", &debug);
  }

}}}}
