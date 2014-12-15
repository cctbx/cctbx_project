#ifndef BOOST_ADAPTBX_EXPORTING_H
#define BOOST_ADAPTBX_EXPORTING_H

#include <boost/mpl/identity.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/placeholders.hpp>

#include <boost/bind.hpp>

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

} // namespace exporting
} // namespace boost_adaptbx

#endif // BOOST_ADAPTBX_EXPORTING_H
