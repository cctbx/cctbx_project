#ifndef MMTBX_GEOMETRY_CONTAINMENT_PYTHON_H
#define MMTBX_GEOMETRY_CONTAINMENT_PYTHON_H

#include <string>

#include <boost/range/adaptor/filtered.hpp>
#include <boost/python/class.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <boost/python/stl_iterator.hpp>

#include <boost_adaptbx/boost_range_python.hpp>
#include <boost_adaptbx/exporting.hpp>

#include <mmtbx/geometry/containment.hpp>

namespace mmtbx
{

namespace geometry
{

namespace containment
{

namespace python
{

struct add_neighbours_from_range_export
{
  typedef void result_type;

  template< typename Checker, typename Range >
  void operator ()(
    boost::python::class_< Checker >& myclass,
    boost::mpl::identity< Range > myexport
    ) const
  {
    myclass.def(
      "add",
      &process< Checker, Range >,
      boost::python::arg( "neighbours" )
      );
  }

  template< typename Checker, typename Range >
  static void process(Checker& checker, const Range& neighbours)
  {
    checker.add( neighbours.begin(), neighbours.end() );
  }
};


template< typename InputRangeList, typename TransformedRange >
struct checker_export
{
  typedef TransformedRange transformed_points_range;

  template< typename Export >
  void operator ()(boost::mpl::identity< Export > myexport)
  {
    using namespace boost::python;
    typedef typename Export::first checker_type;
    typedef typename Export::second name_type;

    std::string prefix = std::string( boost::mpl::c_str< name_type >::value );
    boost_adaptbx::python::generic_range_wrapper< typename checker_type::storage_type >
      ::wrap( ( prefix + "spheres_range" ).c_str() );

    class_< checker_type > myclass( ( prefix + "_checker" ).c_str(), no_init );
    myclass.def( init<>() )
      .def(
        "add",
        add_neighbours_from_list< checker_type >,
        arg( "neighbours" )
        )
      .def(
        "neighbours",
        &checker_type::neighbours,
        return_internal_reference<>()
        )
      .def(
        "__call__",
        &checker_type::operator (),
        arg( "point" )
        )
      ;

    boost_adaptbx::exporting::method_list< InputRangeList >::process(
      myclass,
      add_neighbours_from_range_export()
      );

    typedef boost::filtered_range< checker_type, transformed_points_range >
      filtered_transformed_points_range;
    boost_adaptbx::python::generic_range_wrapper<
      filtered_transformed_points_range
      >
      ::wrap( "filtered_transformed_points_range" );

    filtered_transformed_points_range
      (*filterfunc)( transformed_points_range&, checker_type ) =
        &boost::adaptors::filter< transformed_points_range, checker_type >;

    def(
      "filter",
      filterfunc,
      with_custodian_and_ward_postcall< 0, 1 >(),
      ( arg( "range" ), arg( "predicate" ) )
      );
  }

  template< typename Checker >
  static void add_neighbours_from_list(Checker& checker, boost::python::object neighbours)
  {
    typedef typename Checker::neighbour_type neighbour_type;
    checker.add(
      boost::python::stl_input_iterator< neighbour_type >( neighbours ),
      boost::python::stl_input_iterator< neighbour_type >()
      );
  }
};


} // namespace python
} // namespace containment
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_CONTAINMENT_PYTHON_H
