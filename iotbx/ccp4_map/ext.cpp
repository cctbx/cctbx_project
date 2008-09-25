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

        CMap_io::ccp4_cmap_get_mapstats(
          mfile, &header_min, &header_max, &header_mean, &header_rms);

        float cell_float[6];
        CMap_io::ccp4_cmap_get_cell(mfile, cell_float);
        std::copy(cell_float, cell_float+6, cell.begin());
        CMap_io::ccp4_cmap_get_grid(mfile, grid.begin());
        CMap_io::ccp4_cmap_get_origin(mfile, origin.begin());
        CMap_io::ccp4_cmap_get_order(mfile, axes_order.begin());
        CMap_io::ccp4_cmap_get_dim(mfile, map_dim.begin());
        sg = CMap_io::ccp4_cmap_get_spacegroup(mfile);

        // based on clipper/ccp4/ccp4_map_io.cpp CCP4MAPfile::import_xmap
        int orderfms[3], orderxyz[3], dim[3], gfms0[3], gfms1[3];
        CMap_io::ccp4_cmap_get_order( mfile, orderfms );
        CMap_io::ccp4_cmap_get_dim( mfile, dim );
        CMap_io::ccp4_cmap_get_origin( mfile, gfms0 );
        int dmode = CMap_io::ccp4_cmap_get_datamode( mfile );
        if ( dmode != 0 && dmode != 2 ) {
          throw std::runtime_error("CCP4CCP4MAPfile: unsupported data mode");
        }

        for ( int i = 0; i < 3; i++ ) gfms1[i] = gfms0[i] + dim[i] - 1;
        for ( int i = 0; i < 3; i++ ) orderxyz[orderfms[i]-1] = i;

        data = af::versa<float, af::flex_grid<> >(
          af::flex_grid<>(
            map_dim[orderxyz[0]],
            map_dim[orderxyz[1]],
            map_dim[orderxyz[2]]));
        af::ref<float, af::c_grid<3> > data_ref(
          data.begin(),
          af::c_grid<3>(af::adapt(data.accessor().all())));

        int n0 = (gfms1[0]-gfms0[0]+1);
        int n1 = n0 * (gfms1[1]-gfms0[1]+1);
        IOTBX_ASSERT(dim[0] == n0);
        IOTBX_ASSERT(dim[0]*dim[1] == n1);
        boost::scoped_array<float> section(new float [n1]);
        int g[3];
        for (g[2]=gfms0[2];g[2]<=gfms1[2];g[2]++) {
          CMap_io::ccp4_cmap_read_section(mfile, section.get());
          if (dmode == 0) { // deal with byte maps
            for (int i=n1-1; i>=0; i--) {
              section[i] = static_cast<float>(
                reinterpret_cast<unsigned char*>(&section[0])[i]);
            }
          }
          int index = 0;
          for (g[1]=gfms0[1];g[1]<=gfms1[1];g[1]++) {
            for (g[0]=gfms0[0];g[0]<=gfms1[0];g[0]++) {
              int i = g[orderxyz[0]]-gfms0[orderxyz[0]];
              int j = g[orderxyz[1]]-gfms0[orderxyz[1]];
              int k = g[orderxyz[2]]-gfms0[orderxyz[2]];
              data_ref(i,j,k) = section[index++];
            }
          }
        }
        CMap_io::ccp4_cmap_close(mfile);
      }

      af::versa<float, af::flex_grid<> > data;
      float header_min;
      float header_max;
      double header_mean;
      double header_rms;
      int sg;
      af::tiny<double, 6> cell;
      af::tiny<int, 3> grid;
      af::tiny<int, 3> origin;
      af::tiny<int, 3> axes_order;
      af::tiny<int, 3> map_dim;
  };

  void
  init_module()
  {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;
    class_<map_reader>("map_reader", no_init)
      .def(init<std::string const&>((arg_("file_name"))))
      .def_readonly("data", &map_reader::data)
      .def_readonly("header_min", &map_reader::header_min)
      .def_readonly("header_max", &map_reader::header_max)
      .def_readonly("header_mean", &map_reader::header_mean)
      .def_readonly("header_rms", &map_reader::header_rms)
      .def_readonly("space_group", &map_reader::sg)
      .add_property("unit_cell", make_getter(&map_reader::cell, rbv()))
      .add_property("map_grid", make_getter(&map_reader::grid, rbv()))
      .add_property("map_origin", make_getter(&map_reader::origin, rbv()))
      .add_property("axes_order", make_getter(&map_reader::axes_order, rbv()))
      .add_property("map_dim", make_getter(&map_reader::map_dim, rbv()))
    ;
  }

}} // namespace iotbx::ccp4_map

BOOST_PYTHON_MODULE(iotbx_ccp4_map_ext)
{
  iotbx::ccp4_map::init_module();
}
