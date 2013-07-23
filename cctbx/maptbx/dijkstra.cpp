#include "skeletons.h"

namespace cctbx { namespace maptbx
{

shortest_paths dijkstra(const skeleton &skelet,
    const std::vector<std::size_t> &molids, std::size_t iatom)
{
  shortest_paths result;
  if( iatom> skelet.maximums.size() )
    throw std::runtime_error("wrong vertex");
  CCTBX_ASSERT( skelet.maximums.size() == molids.size() );
  std::size_t molid = molids[iatom];
  std::size_t nv = std::count(molids.begin(), molids.end(), molid);
  result.distances.resize(nv,0);
  result.predecessors.resize(nv,0);
  std::vector<std::size_t> Q(nv);
  for(std::size_t i=0; i<Q.size(); ++i)
    Q[i] = i;
  std::size_t j=0;
  std::size_t inv = std::numeric_limits<std::size_t>::max();
  do
  {
    std::size_t mloc = std::min_element(result.distances.begin(),
       result.distances.end()) - result.distances.begin();
    mloc = Q[mloc];
    if( mloc == inv )
      break;
    if( j>nv )
      throw std::logic_error("index is out of range");
    ++j;
    for(joins_t::const_iterator i=skelet.joins.begin();
        i!=skelet.joins.end(); ++i)
    {
      const join_t &join = *i; // kelet.joins[i];
      std::size_t adj = 0;
      if( join.ilt == mloc )
        adj = join.igt;
      if( join.igt == mloc )
        adj = join.ilt;
      if( adj!=0 )
      {
        if( result.distances[adj] > result.distances[mloc]+1 )
        {
          result.distances[adj] = result.distances[mloc] + 1;
          result.predecessors[adj] = mloc;
        }
      }
    }
    Q[mloc] = inv;
  } while(true);
  return result;
}

}}
