#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include <iotbx/error.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded_periodic.h>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <cmaplib.h>

#if !defined(CCP4_BYTE)
# define CCP4_BYTE BYTE
#endif
#if !defined(CCP4_FLOAT32)
# define CCP4_FLOAT32 FLOAT32
#endif

namespace iotbx {

//! Interfaces to CCP4 cmaplib.
/*! Links:
      http://www.ccp4.ac.uk/dist/html/maplib.html
      http://www.ccp4.ac.uk/html/C_library/cmaplib_8h.html
 */
namespace ccp4_map {

  namespace af = scitbx::af;

  void
  cmap_close_ptr_deleter(
    CMap_io::CMMFile* ptr)
  {
    if (ptr != 0) {
      CMap_io::ccp4_cmap_close(ptr);
    }
  }

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
        if (datamode != CCP4_BYTE && datamode != CCP4_FLOAT32) {
          throw std::runtime_error(
            "iotbx.ccp4_map: unsupported map data mode.");
        }
        CMap_io::ccp4_cmap_get_mapstats(
          mfile.get(), &header_min, &header_max, &header_mean, &header_rms);
        CMap_io::ccp4_cmap_get_grid(mfile.get(), unit_cell_grid.begin());
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
        int order_xyz[3];
        {
          int order_crs[3]; // column-row-section = fast-medium-slow
          CMap_io::ccp4_cmap_get_order(mfile.get(), order_crs);
          for(unsigned i=0;i<3;i++) {
            IOTBX_ASSERT(order_crs[i] >= 1);
            IOTBX_ASSERT(order_crs[i] <= 3);
            order_xyz[order_crs[i]-1] = i;
          }
        }
        af::flex_grid<>::index_type fg_origin;
        for(unsigned i=0;i<3;i++) {
          fg_origin.push_back(origin[order_xyz[i]]);
        }
        af::flex_grid<>::index_type fg_last;
        for(unsigned i=0;i<3;i++) {
          fg_last.push_back(origin[order_xyz[i]]+dim[order_xyz[i]]);
        }
        unsigned n_crs[3];
        for(unsigned i=0;i<3;i++) {
          n_crs[i] = static_cast<unsigned>(dim[i]);
        }
        data = af::versa<float, af::flex_grid<> >(
          af::flex_grid<>(fg_origin, fg_last, true));
        af::ref<float, af::c_grid<3> > data_ref(
          data.begin(),
          af::c_grid<3>(
            n_crs[order_xyz[0]],
            n_crs[order_xyz[1]],
            n_crs[order_xyz[2]]));
        unsigned section_size = n_crs[0] * n_crs[1];
        boost::scoped_array<float> section(new float [section_size]);
        unsigned char* section_char = 0;
        if (datamode == CCP4_BYTE) {
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
              unsigned i = i_crs[order_xyz[0]];
              unsigned j = i_crs[order_xyz[1]];
              unsigned k = i_crs[order_xyz[2]];
              if (datamode == CCP4_BYTE) {
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
      af::tiny<int, 3> unit_cell_grid;
      af::tiny<double, 6> unit_cell_parameters;
      int space_group_number;
      af::versa<float, af::flex_grid<> > data;

  };

  void
  write_ccp4_map_p1_cell(
    std::string const& file_name,
    cctbx::uctbx::unit_cell const& unit_cell,
    cctbx::sgtbx::space_group const& space_group,
    af::int3 const& gridding_first,
    af::int3 const& gridding_last,
    af::const_ref<double, af::c_grid_padded_periodic<3> > const& map_data,
    af::const_ref<std::string> const& labels)
  {
    IOTBX_ASSERT(labels.size() <= 10);
    // TODO write symmetry operators
    boost::shared_ptr<CMap_io::CMMFile> mfile(
      static_cast<CMap_io::CMMFile*>(
        CMap_io::ccp4_cmap_open(file_name.c_str(), O_WRONLY)),
      cmap_close_ptr_deleter);
    if (mfile.get() == 0) {
      throw std::runtime_error(
        "iotbx.ccp4_map: error opening file for writing: \""
        + file_name + "\"");
    }
    CMap_io::ccp4_cmap_set_datamode(mfile.get(), CCP4_FLOAT32);
    for (int i = 0; i < labels.size(); i++) {
      CMap_io::ccp4_cmap_set_label(mfile.get(), labels[i].c_str(), i);
    }
    // symmetry
    af::double6 const& unit_cell_parameters = unit_cell.parameters();
    float cell_float[6];
    for(unsigned i=0;i<6;i++) {
      cell_float[i] = static_cast<float>(unit_cell_parameters[i]);
    }
    CMap_io::ccp4_cmap_set_cell(mfile.get(), cell_float);
    int space_group_number = space_group.type().number();
    CMap_io::ccp4_cmap_set_spacegroup(mfile.get(), space_group_number);
    int grid[3];
    af::tiny<int, 3> n_real(af::adapt(map_data.accessor().focus()));
    std::copy(n_real.begin(), n_real.end(), grid);
    CMap_io::ccp4_cmap_set_grid(mfile.get(), grid);
    int order[3] = {3, 2, 1};
    CMap_io::ccp4_cmap_set_order(mfile.get(), order);
    int dim[3];
    dim[0] = gridding_last[2] - gridding_first[2] + 1;
    dim[1] = gridding_last[1] - gridding_first[1] + 1;
    dim[2] = gridding_last[0] - gridding_first[0] + 1;
    CMap_io::ccp4_cmap_set_dim(mfile.get(), dim);
    int origin[3];
    origin[0] = gridding_first[2];
    origin[1] = gridding_first[1];
    origin[2] = gridding_first[0];
    CMap_io::ccp4_cmap_set_origin(mfile.get(), origin);
    unsigned section_size = dim[0] * dim[1];
    boost::scoped_array<float> section(new float [section_size]);
    for (int i = gridding_first[0]; i <= gridding_last[0]; i++) {
      unsigned index = 0;
      for (int j = gridding_first[1]; j <= gridding_last[1]; j++) {
        for (int k = gridding_first[2]; k <= gridding_last[2]; k++) {
          section[index++] = static_cast<float>(map_data(i,j,k));
        }
      }
      CMap_io::ccp4_cmap_write_section(mfile.get(), section.get());
    }
  }

  void
  init_module()
  {
    using namespace boost::python;
    def("write_ccp4_map", write_ccp4_map_p1_cell, (
        arg("file_name"),
        arg("unit_cell"),
        arg("space_group"),
        arg("gridding_first"),
        arg("gridding_last"),
        arg("map_data"),
        arg("labels")));
    typedef return_value_policy<return_by_value> rbv;
    class_<map_reader>("map_reader", no_init)
      .def(init<std::string const&>((arg("file_name"))))
      .def_readonly("header_min", &map_reader::header_min)
      .def_readonly("header_max", &map_reader::header_max)
      .def_readonly("header_mean", &map_reader::header_mean)
      .def_readonly("header_rms", &map_reader::header_rms)
      .add_property("unit_cell_grid",
        make_getter(&map_reader::unit_cell_grid, rbv()))
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
