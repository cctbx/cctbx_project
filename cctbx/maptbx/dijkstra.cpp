// UNDER CONSTRUCTION
#include <utility>

#include "skeletons.h"
#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/config.hpp>

#include <cctbx/error.h>

namespace cctbx { namespace maptbx
{

typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS>
  graph_t;

std::pair<size_t,std::vector<int> > skeleton_components(const joins_t &joins,
  std::size_t nv)
{
  namespace b = boost;
  graph_t g(nv);
  for( const auto & j : joins )
    b::add_edge(j.ilt, j.igt, g);
  size_t nvj = b::num_vertices(g);
  CCTBX_ASSERT( nv == nvj );
  std::vector<int> cm(nv);
  size_t n = b::connected_components(g,&cm[0]);
  return std::make_pair(n,std::move(cm));
}

std::vector<std::size_t> mask_components(marks_t &mask,
  const std::vector<int> &components)
{
  std::vector<std::size_t> sizes(components.size(),0);
  for( auto & m : mask )
  {
    if( m!=0 )
    {
      m = components.at(m);
      ++(sizes.at(m));
    }
  }
  return sizes;
}

void mask_density_map(asymmetric_map::data_ref_t map, const marks_t &mask,
  unsigned val)
{
  // CCTBX_ASSERT( map.accessor() == mask.accessor() );

  CCTBX_ASSERT( map.size() == mask.size() );
  for(size_t i=0; i<map.size(); ++i)
  {
    if( mask[i]!=val )
      map[i] = 0.;
  }
}


// consider useing boost::graph
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
