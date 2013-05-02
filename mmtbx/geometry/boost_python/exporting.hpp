#ifndef MMTBX_GEOMETRY_EXPORTING_H
#define MMTBX_GEOMETRY_EXPORTING_H

#include <boost/mpl/identity.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/placeholders.hpp>

#include <boost/bind.hpp>


namespace mmtbx
{

namespace geometry
{

namespace exporting
{

template< typename ExportList, typename ExportSpec >
struct class_list
{
  static void process()
  {
    boost::mpl::for_each<
      ExportList,
      boost::mpl::make_identity< boost::mpl::placeholders::_ >
      >( ExportSpec() );
  }
};

template< typename ExportList, typename ExportSpec >
struct method_list
{
  static void process(typename ExportSpec::export_class_type& myclass)
  {
    boost::mpl::for_each<
      ExportList,
      boost::mpl::make_identity< boost::mpl::placeholders::_ >
      >( boost::bind( ExportSpec(), myclass, _1 ) );
  }
};

} // namespace exporting
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_EXPORTING_H
