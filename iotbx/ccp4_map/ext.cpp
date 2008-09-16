#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include <iotbx/error.h>
#include <cctbx/maptbx/accessors/c_grid_padded_p1.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <boost/scoped_array.hpp>

#include <cmaplib.h>

namespace iotbx { namespace ccp4_map {

  namespace af = scitbx::af;

  class map_reader
  {
    public:
      map_reader() {}

      map_reader(
        std::string const& file_name)
      {
        CMap_io::CMMFile* mfile = static_cast<CMap_io::CMMFile*>(
          CMap_io::ccp4_cmap_open(file_name.c_str(), O_RDONLY));
        IOTBX_ASSERT(mfile != 0);
        IOTBX_ASSERT(CMap_io::ccp4_cmap_get_datamode(mfile) == FLOAT32);

        float cell_float[6];
        CMap_io::ccp4_cmap_get_cell(mfile, cell_float);
        std::copy(cell_float, cell_float+6, cell.begin());
        CMap_io::ccp4_cmap_get_grid(mfile, grid.begin());
        CMap_io::ccp4_cmap_get_origin(mfile, origin.begin());
        CMap_io::ccp4_cmap_get_order(mfile, axes_order.begin());
        CMap_io::ccp4_cmap_get_dim(mfile, map_dim.begin());
        sg = CMap_io::ccp4_cmap_get_spacegroup(mfile);
        CMap_io::ccp4_cmap_get_mapstats(
          mfile, &min, &max, &average, &standard_deviation);

        data = af::versa<float, af::flex_grid<> >(
          af::flex_grid<>(map_dim[0], map_dim[1], map_dim[2]));
        af::ref<float, af::c_grid<3> > data_ref(
          data.begin(),
          af::c_grid<3>(af::adapt(data.accessor().all())));

#ifdef THIS_DOES_NOT_SEEM_TO_BE_CORRECT
        int num_cols = map_dim[0];
        int num_rows = map_dim[1];
        boost::scoped_array<float> row(new float [num_rows]);
        int k = 0;
        int sec = -1;
        int nsec = num_rows * num_cols;
        IOTBX_ASSERT(nsec != 0);
        for (int i=0; i<num_rows; i++) {
          if (k % nsec == 0) sec++;
          IOTBX_ASSERT(CMap_io::ccp4_cmap_read_row(mfile, row.get()) == 1);
          for (int j=0; j<num_cols; j++) {
            data_ref(j,i,sec) = row[j];
            k++;
          }
        }
#endif
        CMap_io::ccp4_cmap_close(mfile);
      }

      af::versa<float, af::flex_grid<> > data;
      int sg;
      af::tiny<double, 6> cell;
      float min;
      float max;
      af::tiny<int, 3> grid;
      af::tiny<int, 3> origin;
      af::tiny<int, 3> axes_order;
      af::tiny<int, 3> map_dim;
      double average;
      double standard_deviation;
  };

  void
  init_module()
  {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;
    class_<map_reader>("map_reader", no_init)
      .def(init<std::string const&>((arg_("file_name"))))
      .def_readonly("data", &map_reader::data)
      .def_readonly("space_group", &map_reader::sg)
      .add_property("unit_cell", make_getter(&map_reader::cell, rbv()))
      .def_readonly("map_min", &map_reader::min)
      .def_readonly("map_max", &map_reader::max)
      .add_property("map_grid", make_getter(&map_reader::grid, rbv()))
      .add_property("map_origin", make_getter(&map_reader::origin, rbv()))
      .add_property("axes_order", make_getter(&map_reader::axes_order, rbv()))
      .add_property("map_dim", make_getter(&map_reader::map_dim, rbv()))
      .def_readonly("average", &map_reader::average)
      .def_readonly("standard_deviation", &map_reader::standard_deviation)
    ;
  }

}} // namespace iotbx::ccp4_map

BOOST_PYTHON_MODULE(iotbx_ccp4_map_ext)
{
  iotbx::ccp4_map::init_module();
}
