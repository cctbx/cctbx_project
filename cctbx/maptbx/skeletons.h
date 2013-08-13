//// DO NOT USE. UNDER CONSTRUCTION.

#include <array>
#include <set>
#include <vector>

#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/accessors/c_grid_padded_periodic.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/vec3.h>

#include <cctbx/maptbx/asymmetric_map.h>
#include <cctbx/sgtbx/space_group_type.h>

namespace cctbx { namespace maptbx
{


typedef scitbx::af::versa<double,scitbx::af::c_grid_padded<3> > map_t;

typedef scitbx::af::const_ref<double,scitbx::af::c_grid_padded<3> > const_map_t;

typedef scitbx::af::versa<unsigned,cctbx::maptbx::asymmetric_map::asu_grid_t >
  marks_t;

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
  std::size_t min_count, grows_count, join_count;
  marks_t marks;
};

skeleton swanson(const const_map_t &map,
  const cctbx::sgtbx::space_group_type &spgr, double sigma);

skeleton swanson(const cctbx::maptbx::asymmetric_map &map, double sigma);

std::vector<std::size_t> find_clusters(const skeleton &skelet);

class shortest_paths
{
public:
  std::vector<std::size_t> predecessors, distances;
};


}}
