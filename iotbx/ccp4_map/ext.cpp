#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include <iotbx/error.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <cmaplib.h>

namespace iotbx {

//! Interfaces to CCP4 cmaplib.
/*! Links:
      http://www.ccp4.ac.uk/dist/html/maplib.html
      http://www.ccp4.ac.uk/html/C_library/cmaplib_8h.html
 */
namespace ccp4_map {

  namespace af = scitbx::af;

  class map_reader
  {
    public:
      map_reader() {}

      map_reader(
        std::string const& file_name)
      {
        boost::shared_ptr<CMap_io::CMMFile> mfile(
          static_cast<CMap_io::CMMFile*>(
            CMap_io::ccp4_cmap_open(file_name.c_str(), O_RDONLY)),
          cmap_close_ptr_deleter);
        if (mfile.get() == 0) {
          throw std::runtime_error(
            "iotbx.ccp4_map: error opening file for reading: \""
            + file_name + "\"");
        }
        int datamode = CMap_io::ccp4_cmap_get_datamode(mfile.get());
        if (datamode != BYTE && datamode != FLOAT32) {
          throw std::runtime_error(
            "iotbx.ccp4_map: unsupported map data mode.");
        }
        CMap_io::ccp4_cmap_get_mapstats(
          mfile.get(), &header_min, &header_max, &header_mean, &header_rms);
        float cell_float[6];
        CMap_io::ccp4_cmap_get_cell(mfile.get(), cell_float);
        std::copy(cell_float, cell_float+6, unit_cell_parameters.begin());
        space_group_number = CMap_io::ccp4_cmap_get_spacegroup(mfile.get());
        int origin[3];
        CMap_io::ccp4_cmap_get_origin(mfile.get(), origin);
        int dim[3];
        CMap_io::ccp4_cmap_get_dim(mfile.get(), dim);
        for(unsigned i=0;i<3;i++) {
          IOTBX_ASSERT(dim[i] >= 1);
        }
        int order[3];
        CMap_io::ccp4_cmap_get_order(mfile.get(), order);
        for(unsigned i=0;i<3;i++) {
          IOTBX_ASSERT(order[i] >= 1);
          IOTBX_ASSERT(order[i] <= 3);
          order[i]--;
        }
        af::flex_grid<>::index_type fg_origin;
        for(unsigned i=0;i<3;i++) {
          fg_origin.push_back(origin[order[i]]);
        }
        af::flex_grid<>::index_type fg_last;
        for(unsigned i=0;i<3;i++) {
          fg_last.push_back(origin[order[i]]+dim[order[i]]);
        }
        unsigned n_crs[3]; // column-row-section = fast-medium-slow
        for(unsigned i=0;i<3;i++) {
          n_crs[i] = static_cast<unsigned>(dim[i]);
        }
        data = af::versa<float, af::flex_grid<> >(
          af::flex_grid<>(fg_origin, fg_last, true));
        af::ref<float, af::c_grid<3> > data_ref(
          data.begin(),
          af::c_grid<3>(
            n_crs[order[0]],
            n_crs[order[1]],
            n_crs[order[2]]));
        unsigned section_size = n_crs[0] * n_crs[1];
        boost::scoped_array<float> section(new float [section_size]);
        unsigned char* section_char = 0;
        if (datamode == BYTE) {
          section_char = reinterpret_cast<unsigned char*>(section.get());
        }
        unsigned i_crs[3];
        for(i_crs[2]=0;i_crs[2]<n_crs[2];i_crs[2]++) {
          if (CMap_io::ccp4_cmap_read_section(
                mfile.get(), section.get()) != 1) {
            throw std::runtime_error(
              "iotbx.ccp4_map: ccp4_cmap_read_section() call failed.");
          }
          unsigned index = 0;
          for(i_crs[1]=0;i_crs[1]<n_crs[1];i_crs[1]++) {
            for(i_crs[0]=0;i_crs[0]<n_crs[0];i_crs[0]++) {
              unsigned i = i_crs[order[0]];
              unsigned j = i_crs[order[1]];
              unsigned k = i_crs[order[2]];
              if (datamode == BYTE) {
                data_ref(i,j,k) = static_cast<float>(section_char[index++]);
              }
              else {
                data_ref(i,j,k) = section[index++];
              }
            }
          }
        }
      }

      float header_min;
      float header_max;
      double header_mean;
      double header_rms;
      af::tiny<double, 6> unit_cell_parameters;
      int space_group_number;
      af::versa<float, af::flex_grid<> > data;

      protected:

        static void
        cmap_close_ptr_deleter(
          CMap_io::CMMFile* ptr)
        {
          if (ptr != 0) {
            CMap_io::ccp4_cmap_close(ptr);
          }
        }
  };

  void
  init_module()
  {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;
    class_<map_reader>("map_reader", no_init)
      .def(init<std::string const&>((arg_("file_name"))))
      .def_readonly("header_min", &map_reader::header_min)
      .def_readonly("header_max", &map_reader::header_max)
      .def_readonly("header_mean", &map_reader::header_mean)
      .def_readonly("header_rms", &map_reader::header_rms)
      .add_property("unit_cell_parameters",
        make_getter(&map_reader::unit_cell_parameters, rbv()))
      .def_readonly("space_group_number", &map_reader::space_group_number)
      .def_readonly("data", &map_reader::data)
    ;
  }

}} // namespace iotbx::ccp4_map

BOOST_PYTHON_MODULE(iotbx_ccp4_map_ext)
{
  iotbx::ccp4_map::init_module();
}
