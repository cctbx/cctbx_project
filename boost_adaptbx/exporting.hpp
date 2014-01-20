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

template< typename ExportList, typename ExportSpec >
struct class_list
{
  static void process(ExportSpec const& export_spec = ExportSpec())
  {
    boost::mpl::for_each<
      ExportList,
      boost::mpl::make_identity< boost::mpl::placeholders::_ >
      >( export_spec );
  }
};

template< typename ExportList, typename ExportSpec >
struct method_list
{
  static void process(
    typename ExportSpec::export_class_type& myclass,
    ExportSpec const& export_spec = ExportSpec()
    )
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
