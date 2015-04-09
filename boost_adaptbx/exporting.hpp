#ifndef BOOST_ADAPTBX_EXPORTING_H
#define BOOST_ADAPTBX_EXPORTING_H

#include <boost/mpl/identity.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/copy_if.hpp>

#include <boost/mpl/set.hpp>
#include <boost/mpl/string.hpp>
#include <boost/mpl/vector.hpp>

#include <boost/python/class.hpp>
#include <boost/python/operators.hpp>

#include <boost/bind.hpp>
#include <boost/type_traits.hpp>

#include <string>

namespace boost_adaptbx
{

namespace exporting
{

template< typename ExportList >
struct class_list
{
  template< typename ExportSpec >
  static void process(ExportSpec const& export_spec)
  {
    boost::mpl::for_each<
      ExportList,
      boost::mpl::make_identity< boost::mpl::placeholders::_ >
      >( export_spec );
  }
};

template< typename ExportList >
struct method_list
{
  template< typename ClassType, typename ExportSpec >
  static void process(ClassType& myclass, ExportSpec const& export_spec)
  {
    boost::mpl::for_each<
      ExportList,
      boost::mpl::make_identity< boost::mpl::placeholders::_ >
      >( boost::bind( export_spec, myclass, _1 ) );
  }
};

template< typename TypeList >
struct unique_type_set
{
  typedef typename boost::mpl::fold<
    TypeList,
    boost::mpl::set<>,
    boost::mpl::insert< boost::mpl::placeholders::_1, boost::mpl::placeholders::_2 >
    >::type type;
};

template< typename TypeList >
struct novel_python_type
{
  typedef typename boost::mpl::copy_if<
    TypeList,
    boost::mpl::or_<
      boost::is_class< boost::mpl::placeholders::_1 >,
      boost::is_union< boost::mpl::placeholders::_1 >
      >,
    boost::mpl::back_inserter< boost::mpl::vector<> >
    >::type type;
};

template< typename Type >
struct python_type_export_traits
{
  static void process(boost::python::class_< Type > myclass)
  {}
};

template< typename Prefix >
struct python_type_export
{
  bool m_enable_equality;
  bool m_enable_comparison;

  python_type_export()
    : m_enable_equality( false ), m_enable_comparison( false )
  {}

  python_type_export& enable_equality_operators()
  {
    m_enable_equality = true;
    return *this;
  }

  python_type_export& enable_comparison_operators()
  {
    m_enable_comparison = true;
    return *this;
  }

  template< typename Type >
  void operator ()(boost::mpl::identity< Type > mytype) const
  {
    std::string prefix( boost::mpl::c_str< Prefix >::value );
    boost::python::class_< Type > myclass( ( prefix + typeid( Type ).name() ).c_str() );

    if ( m_enable_equality )
    {
      enable_equality_operators( myclass );
    }

    if ( m_enable_comparison )
    {
      enable_comparison_operators( myclass );
    }

    python_type_export_traits< Type >::process( myclass );
  }

  template< typename Type >
  void enable_equality_operators(boost::python::class_< Type >& myclass) const
  {
    using namespace boost::python;
    myclass.def( self == self )
      .def( self != self )
      ;
  }

  template< typename Type >
  void enable_comparison_operators(boost::python::class_< Type >& myclass) const
  {
    using namespace boost::python;
    myclass.def( self < self )
      .def( self <= self )
      .def( self > self )
      .def( self >= self )
      ;
  }
};

} // namespace exporting
} // namespace boost_adaptbx

#endif // BOOST_ADAPTBX_EXPORTING_H
