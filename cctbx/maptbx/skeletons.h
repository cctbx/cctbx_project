//// DO NOT USE. UNDER CONSTRUCTION.

#include <array>
#include <set>

#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/accessors/c_grid_padded_periodic.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <cctbx/xray/scatterer.h>

namespace cctbx { namespace maptbx
{


typedef scitbx::af::versa<double,scitbx::af::c_grid_padded<3> > map_t;

typedef scitbx::af::const_ref<double,scitbx::af::c_grid_padded<3> > const_map_t;

typedef scitbx::af::versa<std::size_t,scitbx::af::c_grid_padded<3> > marks_t;

using std::tuple;
using std::array;
using std::get;

typedef scitbx::vec3<int> int3_t;
typedef tuple< int3_t, double > xyzm_t;

struct join_t
{
  std::size_t ilt, igt;
};

inline bool operator< (const join_t &a, const join_t &b)
{
  return
    (a.ilt < b.ilt) ||
    (a.ilt == b.ilt && a.igt < b.igt )
  ;
}

typedef std::set<join_t> joins_t;

class skeleton
{
public:
  std::vector<int3_t> maximums;
  joins_t joins;
};

skeleton swanson(const_map_t &map, double sigma=3.);

std::vector<std::size_t> find_clusters(const skeleton &skelet);

class shortest_paths
{
public:
  std::vector<std::size_t> predecessors, distances;
};


}}
