#include "cctbx/maptbx/skeletons.h"
#include "iotbx/xplor/map_reader.h"

#include <iostream>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>

using namespace cctbx::maptbx;
using namespace cctbx::sgtbx;
using namespace cctbx::uctbx;
using namespace iotbx::xplor;
using namespace std;
const char nl='\n';

static string symbol(const list<string> &title_lines)
{
  for(const auto &s : title_lines)
  {
    size_t pos = s.find("SPACE GROUP HALL");
    if( pos!=string::npos )
    {
      pos = s.find(':', pos+14);
      if( pos!=string::npos )
        return s.substr(++pos);
    }
  }
  return "";
}

static string symbol(int argc, char *argv[])
{
  for(int i=1; i<argc; ++i)
  {
    if( string(argv[i]).find("-space_group")!=string::npos
      && 1+i<argc)
    {
      return argv[1+i];
    }
  }
  return "";
}

asymmetric_map make_map(int argc, char *argv[], unit_cell &cell)
{
  if( argc<2 )
    throw cctbx::error("Provide map file on command line");
  string file_name(argv[1]);
  map_reader rdr(file_name);
  cell = unit_cell(rdr.unit_cell_parameters);
  string spgrs = symbol(rdr.title_lines);
  if( spgrs.empty() )
    spgrs = symbol(argc, argv);
  if( spgrs.empty() )
    throw cctbx::error("Space group symbol is missing, use: -space_group "
      "command line option");
  space_group_symbols smb("Hall: "+spgrs);
  space_group group(smb);
  scitbx::af::c_interval_grid<3> grid(rdr.grid_first, rdr.grid_last);
  long sz = boost::accumulate(rdr.grid_size, 1L, multiplies<long>());
  if( grid.size_1d()*2UL <= sz && group.order_z()>1U )
    return asymmetric_map(group.type(), move(rdr.data), rdr.grid_size);
  CCTBX_ASSERT( grid.size_1d() == sz );
  cout << "Creating asymmetric map from whole unit cell map\n";
  return asymmetric_map(group.type(), rdr.data.const_ref());
}

void run(int argc, char *argv[])
{
  unit_cell cell;
  asymmetric_map amap = make_map(argc, argv, cell);
  CCTBX_ASSERT( cell.volume()>1E-3 );
  cout
    << "Crystal: " << cell.parameters() << "  (" << cell.volume() << ")  "
    << amap.space_group().type().lookup_symbol()
    << "\nMap grid: [" << amap.box_begin() << ", " << amap.box_end() << ") on "
    << amap.unit_cell_grid_size()
    << ";  size: " << amap.data().accessor().size_1d()
    << nl;
  double sigma = 3.;
  skeleton skelet = swanson(amap, sigma);
  auto components=skeleton_components(skelet.joins, skelet.maximums.size()+1U);
  cout << "Blob count: " << components.first << nl;
  auto sizes = mask_components(skelet.marks, components.second);
  {
    asymmetric_map tmap = amap.explicit_copy();
    const char fn[] = "higest_blob-xplor.map";
    cout << "Most intense blob size? " << sizes[1] << "  DEBUG: " << sizes[0]
      << "\nSaving blob in: " << fn << nl;
    mask_density_map(tmap.data_ref(), skelet.marks, 1);
    tmap.save(fn, cell);
  }
  {
    CCTBX_ASSERT( !sizes.empty() );
    auto i = boost::max_element(sizes);
    long m = i-sizes.begin();
    const char fn[] = "largest_blob-xplor.map";
    cout << "Largest blob: " << m << " size: " << *i
      << "\nSaving in: " << fn << nl;
    CCTBX_ASSERT( m!=0 );
    // boost::partial_sort(sizes, sizes.begin()+3);
    mask_density_map(amap.data_ref(), skelet.marks, m);
    amap.save(fn, cell);
  }
}

int main(int argc, char *argv[])
{
  try {
    run(argc, argv);
    return 0;
  }
  catch(std::exception &err)
  {
    std::cerr << "ERROR! " << err.what() << std::endl;
  }
  return -1;
}
