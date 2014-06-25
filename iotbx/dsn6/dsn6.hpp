
// DSN6 map output
//
// FIXME I am not sure about the scaling here, need to re-examine

#include <iotbx/error.h>
#include <cctbx/uctbx.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded_periodic.h>

#include <iostream>
#include <fstream>

// http://stackoverflow.com/questions/3916097/integer-byte-swapping-in-c
#define SWAP_BYTES_INTEGER(n) \
  short hibyte = (n & 0xff00) >> 8; \
  short lobyte = (n & 0xff); \
  n = lobyte << 8 | hibyte;

namespace iotbx { namespace dsn6 {

  namespace af = scitbx::af;

  // http://www.uoxray.uoregon.edu/tnt/manual/node104.html
  // http://lists.bioxray.dk/pipermail/o-info/2002/001921.html
  void write_dsn6_map(
    std::string const& file_name,
    cctbx::uctbx::unit_cell const& unit_cell,
    af::int3 const& gridding_first,
    af::int3 const& gridding_last,
    af::const_ref<double, af::c_grid_padded_periodic<3> > const& map_data)
  {
    af::double6 const& unit_cell_parameters = unit_cell.parameters();
    af::tiny<int, 3> n_real(af::adapt(map_data.accessor().focus()));
    int dim[3];
    for (unsigned i = 0; i < 3; i++) {
      dim[i] = (gridding_last[i] - gridding_first[i] + 1);
    }
    double rho_max = 0, rho_min = 2.e16;
    for (unsigned i = 0; i < map_data.size(); i++) {
      double rho = map_data[i];
      if (rho > rho_max) rho_max = rho;
      if (rho < rho_min) rho_min = rho;
    }
    int i1 = 80; //, i2 = 100;
    for (unsigned i = 0; i < 6; i++) {
      double param = unit_cell_parameters[i];
      if (i1 * param > 32760) {
        i1 = std::min(i1, int (32760 / (float) ((int) (0.999 + param))));
        i1 = std::max(1, i1);
      }
    }
    double rho_range = rho_max - rho_min;
    IOTBX_ASSERT(rho_range > 0);
    std::ofstream mapfile(file_name.c_str(), std::ios::out | std::ios::binary);
    // origin
    for (unsigned i = 0; i < 3; i++) {
      short n = static_cast<short>(gridding_first[i]);
      SWAP_BYTES_INTEGER(n);
      mapfile.write((char *) &n, 2);
    }
    // extent
    for (unsigned i = 0; i < 3; i++) {
      short n = static_cast<short>(dim[i]);
      SWAP_BYTES_INTEGER(n);
      mapfile.write((char *) &n, 2);
    }
    // unit cell grid dimensions
    for (unsigned i = 0; i < 3; i++) {
      short n = static_cast<short>(n_real[i]);
      SWAP_BYTES_INTEGER(n);
      mapfile.write((char *) &n, 2);
    }
    // unit cell parameters
    for (unsigned i = 0; i < 6; i++) {
      short n = static_cast<short>(i1 *
        unit_cell_parameters[i]);
      SWAP_BYTES_INTEGER(n);
      mapfile.write((char *) &n, 2);
    }
    // header(16)
    double header_16 = 100 * 250 / rho_range;
    {
      short n = static_cast<short>(header_16);
      SWAP_BYTES_INTEGER(n);
      mapfile.write((char *) &n, 2);
    }
    // header(17)
    double header_17 = (3*rho_max - 253*rho_min) / rho_range;
    {
      short n = static_cast<short>(header_17);
      SWAP_BYTES_INTEGER(n);
      mapfile.write((char *) &n, 2);
    }
    // header(18) - unit cell scaling factor
    {
      short n = static_cast<short>(i1);
      SWAP_BYTES_INTEGER(n);
      mapfile.write((char *) &n, 2);
    }
    // header(19) - just 100
    {
      short n = 100;
      SWAP_BYTES_INTEGER(n);
      mapfile.write((char *) &n, 2);
    }
    for (unsigned i = 20; i <= 256; i++) {
      char n = 0;
      mapfile.write(&n, 1); mapfile.write(&n, 1);
    }
    // END of header
    int n_blocks[3];
    for (unsigned i = 0; i < 3; i++) {
      n_blocks[i] = (int) std::ceil(((float) dim[i]) / 8.0);
    }
    char brick[512];
    char brick_swapped[512];
    double scale_factor = header_16 / 100;
    for (int zz = 0; zz < n_blocks[2]; zz++) {
      for (int yy = 0; yy < n_blocks[1]; yy++) {
        for (int xx = 0; xx < n_blocks[0]; xx++) {
          unsigned i_byte = 0;
          for (int z = 0; z < 8; z++) {
            for (int y = 0; y < 8; y++) {
              for (int x = 0; x < 8; x++) {
                int i = (xx * 8) + x + gridding_first[0];
                int j = (yy * 8) + y + gridding_first[1];
                int k = (zz * 8) + z + gridding_first[2];
                double rho = map_data(i,j,k);
                double rho_scaled = (rho - rho_min) * scale_factor;
                char rhoc = static_cast<char>(rho_scaled);
                brick[i_byte++] = rhoc;
          }}}
          IOTBX_ASSERT(i_byte == 512);
          // XXX not sure about the byte swapping here, but this was
          // necessary to be able to read DSN6 maps from the EDS...
          for (unsigned i = 0; i < 256; i++) {
            brick_swapped[i*2] = brick[(i*2)+1];
            brick_swapped[(i*2)+1] = brick[i*2];
          }
          mapfile.write(brick_swapped, 512);
    }}}
    mapfile.close();
  }

}} // namespace iotbx::dsn6
