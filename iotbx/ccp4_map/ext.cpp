#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include <cctbx/maptbx/accessors/c_grid_padded_p1.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/accessors/c_grid.h>

#include <cmaplib.h>
#include <iostream>

namespace af = scitbx::af;

namespace iotbx { namespace ccp4_map {

// ===

  void
  test_read(const char* file_name)
  {
    int ii,sg,num_rows,num_cols,points;
//    const void *row,*data;
    float *row,*data;
    float cell[6],min,max;
    int grid[3],origin[3],axes_order[3],map_dim[3];
    double mean,rms;

    CMap_io::CMMFile* mfile = static_cast<CMap_io::CMMFile*>(
      CMap_io::ccp4_cmap_open(file_name, O_RDONLY));

      CMap_io::ccp4_cmap_get_cell(mfile, cell);
      CMap_io::ccp4_cmap_get_grid(mfile, grid);
      CMap_io::ccp4_cmap_get_origin(mfile, origin);
      CMap_io::ccp4_cmap_get_order(mfile, axes_order);
      CMap_io::ccp4_cmap_get_dim(mfile, map_dim);
      sg = CMap_io::ccp4_cmap_get_spacegroup(mfile);
      CMap_io::ccp4_cmap_get_mapstats(mfile, &min, &max, &mean, &rms);

      num_cols=map_dim[0];
      num_rows=map_dim[1];
      points=map_dim[0]*map_dim[1]*map_dim[2];

        std::cout << "nr:" << num_rows << "np:" << points << std::endl;
        std::cout << "dx:" << min << std::endl;
        std::cout << "dy:" << max << std::endl;
        std::cout << "dz:" << mean << std::endl;
        std::cout << "ox:" << rms << std::endl;
        std::cout << "oy:" << origin[1] << std::endl;
        std::cout << "oz:" << origin[2] << std::endl;
        std::cout << "gx: " << grid[0] << "gy: " << grid[1] << "gz: " << grid[2] << std::endl;
        std::cout << "f: " << axes_order[0] << "m: " << axes_order[1] << "s: " << axes_order[2] << std::endl;

      row = new float [num_rows];
      data = new float [points];
      int k=0;
      for (int i=0; i<num_rows; i++)
       {
         ii=CMap_io::ccp4_cmap_read_row(mfile, row);
         for (int j=0; j<num_cols; j++)
            {*(data+k)=*(row+j);
             k++;
                }
        }
        std::cout << "data:" << *(data+10) << *(data+20) << std::endl;

    CMap_io::ccp4_cmap_close(mfile);
  }

// ======= map read as of xplor

  class map_reader
  {
    public:
      map_reader() {}

      map_reader(
        std::string const& file_name)
      {
    int ii,num_rows,num_cols,points;
    float *row;

    CMap_io::CMMFile* mfile = static_cast<CMap_io::CMMFile*>(
      CMap_io::ccp4_cmap_open(file_name.c_str(), O_RDONLY));

      float cell_float[6];
      CMap_io::ccp4_cmap_get_cell(mfile, cell_float);
      std::copy(cell_float, cell_float+6, cell.begin());
      CMap_io::ccp4_cmap_get_grid(mfile, grid.begin());
      CMap_io::ccp4_cmap_get_origin(mfile, origin.begin());
      CMap_io::ccp4_cmap_get_order(mfile, axes_order.begin());
      CMap_io::ccp4_cmap_get_dim(mfile, map_dim.begin());
      sg = CMap_io::ccp4_cmap_get_spacegroup(mfile);
      CMap_io::ccp4_cmap_get_mapstats(mfile, &min, &max, &mean, &rms);

      num_cols=map_dim[0];
      num_rows=map_dim[1];
      points=map_dim[0]*map_dim[1]*map_dim[2];

      data = af::versa<float, af::flex_grid<> >(
        af::flex_grid<>(map_dim[0], map_dim[1], map_dim[2]));
      af::ref<float, af::c_grid<3> > data_ref(
        data.begin(),
        af::c_grid<3>(af::adapt(data.accessor().all())));
//...
      row = new float [num_rows];

      int k=0;
      int sec=-1;
      int nsec=num_rows*num_cols;
      for (int i=0; i<num_rows; i++)
       {
        if (k % nsec == 0) sec++;
         ii=CMap_io::ccp4_cmap_read_row(mfile, row);
         for (int j=0; j<num_cols; j++)
            {
              data_ref(j,i,sec) = *(row+j);
              k++;
                }
        }

    CMap_io::ccp4_cmap_close(mfile);

//...
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
    double mean;
    double rms;
  };

// ======= end

  void
  init_module()
  {
    using namespace boost::python;
    def("test_read", test_read, (arg_("file_name")));

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
    ;
  }

}} // namespace iotbx::ccp4_map

BOOST_PYTHON_MODULE(iotbx_ccp4_map_ext)
{
  iotbx::ccp4_map::init_module();
}
