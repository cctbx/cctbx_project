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

#include <boost_adaptbx/iterator_range.h>
#include <scitbx/boost_python/container_conversions.h>

#include <boost/format.hpp>

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

  /* Three goals may be pursued in relation to inheritance:

      (1) Polymorphic calls from Python to C++
          e.g. given a Python object s of type independent_site_parameter,
          properly dispatch s.scatterers() to the appropriate ancestor in C++

      (2) Given a C++ function taking an argument of type site_parameter *,
          make sure that from the Python side one can pass to it an
          object of type special_position_site_parameter, e.g.

      (3) Make it possible to write new constraints in Python, which requires
          polymorphic calls from C++ to Python.

   Goals (1) and (2) are easily achieved with the usual flavour
   of Boost.Python wrapper. That is to say that virtual or pure virtual
   member functions are wrapped as any other member function, once and only
   once in the class where they are introduced in the inheritance hierarchy,
   e.g. linearise is wrapped as part of the wrapper of class parameter.
   Boost.Python then auto-magically care of all the details. The only necessary
   trick is to instruct Boost.Python to use boost::noncopyable storage for
   any class with non-overriden pure virtual member functions,
   as otherwise the code won't compile since Boost.Python will try
   to construct an object.

   Goal (3) requires to write an inordinate amount of boiler code. The
   Boost.Python tutorial looks easy enough on that subject but the featured
   example has only one base class and one derived class. For the non-trivial
   hiearchy of ours, it scales very badly. In an earlier version, we had
   attempted to reach that goal, but we have now decided against it,
   until the need for writing reparametrisations in Python arises.

   An orthogonal question is that of memory management for the instances
   of parameter classes (and heirs). Since class reparametrisation owns the
   pointers to the objects it refers too, we need to be very carefull that
   the Python side releases ownership and transfers it correctly. The first
   step to do that is to instruct Boost.Python to use std::auto_ptr storage
   for all instanciable classes (i.e. without any non-overriden pure
   virtual member functions). The second step is the Boost.Python binding
   for class reparametrisation (see below).
   */

  struct parameter_wrapper
  {
    typedef parameter wt;

    static void wrap() {
      using namespace boost::python;
      return_internal_reference<> rir;
      class_<wt,
             boost::noncopyable>("parameter", no_init)
        .add_property("index", &wt::index)
        .add_property("n_arguments", &wt::n_arguments)
        .def("argument", &wt::argument, rir)
        .add_property("is_independent", &wt::is_independent)
        .add_property("is_root", &wt::is_root)
        .def("size", &wt::size)
        .add_property("is_variable", &wt::is_variable)
        .def("evaluate", &wt::evaluate, arg("unit_cell"))
        .def("linearise", &wt::linearise,
             (arg("unit_cell"), arg("jacobian_transpose")))
        ;
    }
  };

  struct asu_parameter_wrapper
  {
    typedef asu_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<parameter>,
             boost::noncopyable>("asu_parameter", no_init)
        .def("component_indices_for",
             &wt::component_indices_for,
             arg("scatterer"))
        .def("store", &wt::store, arg("unit_cell"))
        .add_property("scatterers", &wt::scatterers)
        ;
    }
  };

  struct single_asu_scatterer_parameter_wrapper
  {
    typedef single_asu_scatterer_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<asu_parameter>,
             boost::noncopyable>("single_asu_scatterer_parameter", no_init)
        ;
    }
  };

  struct scalar_parameter_wrapper
  {
    typedef scalar_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<parameter>,
             boost::noncopyable>("scalar_parameter", no_init)
        .def_readonly("value", &wt::value);
    }
  };

  struct independent_scalar_parameter_wrapper
  {
    typedef independent_scalar_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<scalar_parameter>,
             std::auto_ptr<wt> >("independent_scalar_parameter", no_init)
        .def(init<double, bool>((arg("value"), arg("variable"))));
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct twin_fraction_parameter_wrapper
  {
    typedef twin_fraction_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<independent_scalar_parameter>,
             std::auto_ptr<wt> >("twin_fraction_parameter", no_init)
        .def(init<cctbx::xray::twin_fraction<double> *>((
          arg("twin_fraction"))));
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct extinction_parameter_wrapper
  {
    typedef extinction_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<independent_scalar_parameter>,
             std::auto_ptr<wt> >("extinction_parameter", no_init)
        .def(init<cctbx::xray::extinction_correction<double> *>((
          arg("extinction"))));
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  template <int N>
  struct small_vector_parameter_wrapper
  {
    typedef small_vector_parameter<N> wt;

    static void wrap() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      std::string name = (
        boost::format("small_%1%_vector_parameter") % N).str();
      class_<wt,
             bases<parameter>,
             boost::noncopyable>(name.c_str(), no_init)
        .add_property("value",
                      make_getter(&wt::value, rbv),
                      make_setter(&wt::value));
    }
  };

  template <int N>
  struct independent_small_vector_parameter_wrapper
  {
    typedef independent_small_vector_parameter<N> wt;

    static void wrap() {
      using namespace boost::python;
      std::string name = (
        boost::format("independent_small_%1%_vector_parameter") % N).str();
      class_<wt,
             bases<small_vector_parameter<N> >,
             std::auto_ptr<wt> >(name.c_str(), no_init)
        .def(init<af::small<double, N> const &, bool>
             ((arg("value"), arg("variable"))));
        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct site_parameter_wrapper
  {
    typedef site_parameter wt;

    static void wrap() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      class_<wt,
             bases<parameter>,
             boost::noncopyable>("site_parameter", no_init)
        .add_property("value",
                      make_getter(&wt::value, rbv),
                      make_setter(&wt::value));
        ;
    }
  };

  struct asu_site_parameter_wrapper
  {
    typedef asu_site_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<site_parameter, single_asu_scatterer_parameter>,
             boost::noncopyable>("asu_site_parameter", no_init)
        ;
    }
  };

  struct independent_site_parameter_wrapper
  {
    typedef independent_site_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<asu_site_parameter>,
             std::auto_ptr<wt> >("independent_site_parameter", no_init)
        .def(init<asu_parameter::scatterer_type *>(arg("scatterer")))
        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct u_star_parameter_wrapper
  {
    typedef u_star_parameter wt;

    static void wrap() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      class_<wt,
             bases<parameter>,
             boost::noncopyable>("u_star_parameter", no_init)
        .add_property("value",
                      make_getter(&wt::value, rbv),
                      make_setter(&wt::value));
        ;
    }
  };

  struct asu_u_star_parameter_wrapper
  {
    typedef asu_u_star_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<u_star_parameter, single_asu_scatterer_parameter>,
             boost::noncopyable>("asu_u_star_parameter", no_init)
        ;
    }
  };

  struct independent_u_star_parameter_wrapper
  {
    typedef independent_u_star_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<asu_u_star_parameter>,
             std::auto_ptr<wt> >("independent_u_star_parameter", no_init)
        .def(init<asu_parameter::scatterer_type *>(arg("scatterer")))
        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct asu_occupancy_parameter_wrapper
  {
    typedef asu_occupancy_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<scalar_parameter, single_asu_scatterer_parameter>,
             boost::noncopyable>("asu_occupancy_parameter", no_init)
        ;
    }
  };

  struct independent_occupancy_parameter_wrapper
  {
    typedef independent_occupancy_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<asu_occupancy_parameter>,
             std::auto_ptr<wt> >("independent_occupancy_parameter", no_init)
        .def(init<asu_parameter::scatterer_type *>(arg("scatterer")))
        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct asu_u_iso_parameter_wrapper
  {
    typedef asu_u_iso_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<scalar_parameter, single_asu_scatterer_parameter>,
             boost::noncopyable>("asu_u_iso_parameter", no_init)
        ;
    }
  };

  struct independent_u_iso_parameter_wrapper
  {
    typedef independent_u_iso_parameter wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<asu_u_iso_parameter>,
             std::auto_ptr<wt> >("independent_u_iso_parameter", no_init)
        .def(init<asu_parameter::scatterer_type *>(arg("scatterer")))
        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct index_range_to_tuple
  {
    static PyObject *convert(index_range const &ir) {
      using namespace boost::python;
      tuple t = ir.is_valid() ? make_tuple(ir.first(), ir.last()) : tuple();
      return incref(t.ptr());
    }

    static PyTypeObject const *get_pytype() { return &PyTuple_Type; }

    index_range_to_tuple() {
      using namespace boost::python;
      to_python_converter<index_range, index_range_to_tuple, true>();
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

      boost_adaptbx::iterator_range_wrapper<wt::range>
      ::wrap("parameter_iterator");

      class_<wt> klass("reparametrisation", no_init);
      klass
        .def(init<uctbx::unit_cell const &>(arg("unit_cell")))
        .def("finalise", &wt::finalise)
        .def("linearise", &wt::linearise)
        .def("store", &wt::store)
        .def("parameters", &wt::parameters)
        .add_property("jacobian_transpose",
                      make_getter(&wt::jacobian_transpose))
        .def("apply_shifts", &wt::apply_shifts)
        .add_property("norm_of_independent_parameter_vector",
                      &wt::norm_of_independent_parameter_vector)
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

    index_range_to_tuple();
    scitbx::boost_python::container_conversions::to_tuple_mapping<
      asu_parameter::scatterer_sequence_type>();

    parameter_wrapper::wrap();
    asu_parameter_wrapper::wrap();
    single_asu_scatterer_parameter_wrapper::wrap();

    scalar_parameter_wrapper::wrap();
    independent_scalar_parameter_wrapper::wrap();

    twin_fraction_parameter_wrapper::wrap();

    extinction_parameter_wrapper::wrap();

    small_vector_parameter_wrapper<3>::wrap();
    independent_small_vector_parameter_wrapper<3>::wrap();

    small_vector_parameter_wrapper<6>::wrap();
    independent_small_vector_parameter_wrapper<6>::wrap();

    site_parameter_wrapper::wrap();
    asu_site_parameter_wrapper::wrap();
    independent_site_parameter_wrapper::wrap();

    u_star_parameter_wrapper::wrap();
    asu_u_star_parameter_wrapper::wrap();
    independent_u_star_parameter_wrapper::wrap();

    asu_occupancy_parameter_wrapper::wrap();
    independent_occupancy_parameter_wrapper::wrap();

    asu_u_iso_parameter_wrapper::wrap();
    independent_u_iso_parameter_wrapper::wrap();

    reparametrisation_wrapper::wrap();

    def("debug", &debug);
  }

}}}}
