//// DO NOT USE. UNDER CONSTRUCTION.

#include <cstdio>
#include <string>
#include <limits>
#include <set>

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <cctbx/error.h>
// iostream is included by this one
#include <cctbx/sgtbx/direct_space_asu/proto/asymmetric_unit.h>

#include "cctbx/maptbx/skeletons.h"

namespace cctbx { namespace maptbx
{

inline bool any_lt(const int3_t &a, const int3_t &b)
{
  for(short j=0; j<3; ++j)
    if( a[j]<b[j] )
      return true;
  return false;
}

inline bool any_lt(const int3_t &a, int b)
{
  for(short j=0; j<3; ++j)
    if( a[j]<b )
      return true;
  return false;
}

inline bool any_ge(const int3_t &a, const int3_t &b)
{
  for(short j=0; j<3; ++j)
    if( a[j]>=b[j] )
      return true;
  return false;
}

// faster version of small size std::set ?
template<typename T, unsigned short N=1000> class array_as_set :
  private scitbx::af::small<T,N>
{
public:
  typedef scitbx::af::small<T,N> buf_t;
  using buf_t::begin;
  using buf_t::end;
  using buf_t::empty;
  using buf_t::size;
  using typename buf_t::iterator;
  using typename buf_t::const_iterator;
  void insert(const T &v)
  {
    iterator l=std::lower_bound(this->begin(), this->end(),v);
    // *l >= v;
    if( l!=this->end() )
    {
      if( *l != v )
        this->buf_t::insert(l,v);
    }
    else
        this->buf_t::insert(l,v);
    //CCTBX_ASSERT( !this->empty() );
    //CCTBX_ASSERT( std::is_sorted(this->begin(), this->end()) );
    //CCTBX_ASSERT(std::adjacent_find(this->begin(),this->end())==this->end());
  }
};

skeleton swanson(const cctbx::maptbx::asymmetric_map &amap, double sigma)
{
  double mean=0., esd=0., mx=-9.E200;
  const auto &map_data = amap.data().const_ref();
  for(std::size_t ii=0; ii<map_data.size(); ++ii)
  {
    //! @todo discard points outside asu
    double m = map_data[ii];
    mean += m;
    esd += m*m;
    if( m > mx )
      mx = m;
  }
  mean /= map_data.size();
  esd = esd/map_data.size() - mean*mean;
  esd = std::sqrt(esd);
  double mapcutoff = mean + esd * sigma;

  std::vector<xyzm_t> xyzm;
  xyzm.reserve(10000);
  //! @todo code duplication: mmtbx::masks::atom_mask::mask_asu
  std::size_t inside=0, tot=0;
  const auto &opt_asu = amap.optimized_asu();
  const double *md = map_data.begin();
  for(auto i3=amap.grid_begin(); !i3.over(); i3.incr(), ++md)
  {
    int3_t p=i3();
    double m = *md;
    if( m>mapcutoff )
    {
      ++tot;
      if( opt_asu.where_is(p) != 0 ) // inside or on the face
      {
        ++inside;
        xyzm_t x(p,m);
        xyzm.push_back( x );
      }
    }
  }
  std::sort(xyzm.begin(), xyzm.end(),
    [](const xyzm_t &a, const xyzm_t &b) {return (get<1>(a)) > (get<1>(b));}
  );

  CCTBX_ASSERT( get<1>(xyzm.front()) == mx );
  marks_t marks(map_data.accessor(), 0UL); // 0 means no mark
  std::size_t nmarks=0, nbonds=0, min_count = 0, grows_count=0, join_count=0;
  const unsigned short cube_size = 26; // 3*3 + (3*3-1) + (3*3) == 3^3 - 1
  typedef int3_t i3t;
  int3_t cube[cube_size] = {
    i3t(-1,-1,-1),  i3t(0,-1,-1),  i3t(1,-1,-1),  i3t(-1,0,-1),  i3t(0,0,-1),
    i3t(1,0,-1),    i3t(-1,1,-1),  i3t(0,1,-1),   i3t(1,1,-1),   i3t(-1,-1,0),
    i3t(0,-1,0),    i3t(1,-1,0),   i3t(-1,0,0),   i3t(1,0,0),    i3t(-1,1,0),
    i3t(0,1,0),     i3t(1,1,0),    i3t(-1,-1,1),  i3t(0,-1,1),   i3t(1,-1,1),
    i3t(-1,0,1),    i3t(0,0,1),    i3t(1,0,1),    i3t(-1,1,1),   i3t(0,1,1),
    i3t(1,1,1)
  };

  auto ibox_min = amap.box_begin();
  auto ibox_max = amap.box_end();
  skeleton result;
  result.maximums.reserve(1000);
  for(std::size_t jj=0; jj<xyzm.size(); ++jj)
  {
    // typedef std::set<std::size_t> feature_set_t;
    typedef array_as_set<std::size_t> feature_set_t;
    feature_set_t featureset;
    int3_t x = get<0>(xyzm[jj]);
    for(unsigned short ic=0; ic<cube_size; ++ic)
    {
      int3_t neighbor = x + cube[ic]; //! @todo check if in cell ?
      if( any_lt(neighbor,ibox_min) || any_ge(neighbor,ibox_max) )
        continue;
      unsigned mark = marks(neighbor);
      if( mark != 0U ) // 0 is not a mark
        featureset.insert(mark);
    }

    const unsigned fssize = featureset.size();
    if( fssize==0U ) // all_eq(featureset, 0) )
    {
      // 2a: new maximum
      ++nmarks;
      marks(x) = nmarks;
      result.maximums.push_back(x);
    }
    else if( fssize == 1 )
    {
      // 2b: part of the growing nodule
      std::size_t mark = *featureset.begin();
      CCTBX_ASSERT( mark != 0 );
      marks(x) = mark;
      ++grows_count;
    }
    else if( fssize == cube_size )
    {
      // 2d: local minimum will not be in the neighborhood of any point
      ++min_count;
      marks(x) = *featureset.begin(); // highest nodule
    }
    else if( fssize>1 )
    {
      // 2c: two or more nodules merging
      ++join_count;
      marks(x) = *featureset.begin(); // highest nodule
      for(feature_set_t::const_iterator i=featureset.begin();
          i!=featureset.end(); ++i)
      {
        std::size_t fi = *i;
        CCTBX_ASSERT( fi>0U && fi<=nmarks );
        feature_set_t::const_iterator j=i;
        ++j;
        for( ; j!=featureset.end(); ++j)
        {
          std::size_t fj = *j;
          CCTBX_ASSERT( fj>0U && fj<=nmarks );
          CCTBX_ASSERT( fi<fj );
          join_t join;
          join.ilt = fi;
          join.igt = fj;
          result.joins.insert(join); // inserts only new joins
        }
      }
    }
  }
  CCTBX_ASSERT( nmarks + grows_count + join_count + min_count == xyzm.size() );
  result.min_count = min_count;
  result.grows_count = grows_count;
  result.join_count = join_count;
  result.marks = std::move(marks);
  result.mapcutoff = mapcutoff;
  return result;
}

std::vector<std::size_t> find_clusters(const skeleton &skelet)
{
  std::size_t nmaxs = skelet.maximums.size(), mol_id=0;
  std::vector<std::size_t> atoms(nmaxs,0);
  for(joins_t::const_iterator i=skelet.joins.begin(); i!=skelet.joins.end();
    ++i)
  {
    const join_t &join = *i;
    long i_cat = join.ilt;
    long i_an = join.igt;
    CCTBX_ASSERT( i_cat < i_an && i_cat>0 && i_an>0 && i_cat<=nmaxs
      && i_an<=nmaxs );
    if( atoms[i_cat-1U]==0 && atoms[i_an-1U]==0 )
    {
      ++mol_id;
      atoms[i_cat] = mol_id;
      atoms[i_an] = mol_id;
    }
    else
    {
      if(atoms[i_cat]!=0 && atoms[i_an]!=0 && atoms[i_cat]!=atoms[i_an])
      {
        std::size_t imn = std::min(atoms[i_cat],atoms[i_an]);
        atoms[i_cat] = imn;
        atoms[i_an] = imn;
      }
      else
      {
        if( atoms[i_cat]!=0 )
          atoms[i_an] = atoms[i_cat];
        if( atoms[i_an]!=0 )
          atoms[i_cat] = atoms[i_an];
      }
    }
  }
  std::size_t imn = 0, imx=0;
  do
  {
    imx = *std::max_element(atoms.begin(), atoms.end());
    for(joins_t::const_iterator i=skelet.joins.begin(); i!=skelet.joins.end();
      ++i)
    {
      const join_t &join = *i;
      long i_cat = join.ilt;
      long i_an = join.igt;
      if( atoms[i_cat]==0 || atoms[i_an]==0 )
        throw std::logic_error("disaster");
      if( atoms[i_cat] != atoms[i_an] )
      {
        imn = std::min(atoms[i_cat],atoms[i_an]);
        atoms[i_cat] = imn;
        atoms[i_an] = imn;
      }
    }
  }while( imn!=0 );
  imx = *std::max_element(atoms.begin(), atoms.end());

  std::vector< array<std::size_t,2> > molecules;
  molecules.reserve(imx);
  std::size_t imxsz = 0;
  for(std::size_t i=0; i<imx; ++i)
  {
    imn = std::count(atoms.begin(), atoms.end(), i);
    array< std::size_t,2 > m = {imn,0};
    molecules.push_back(m);
    if( imn!=0 )
      ++mol_id;
    if( imn>imxsz )
    {
      imx = i;
      imxsz = imn;
    }
  }
  if( imx!=0 )
  {
    mol_id = 1;
    for(std::size_t i=0; i<molecules.size(); ++i)
    {
      get<1>(molecules[i]) = mol_id;
      ++mol_id;
    }
  }
  std::vector<std::size_t> molids;
  for(std::size_t i=0; i<atoms.size(); ++i)
  {
    if( atoms[i] == 0 )
      continue; // ?
    molids.push_back( get<1>(molecules[atoms[i]]) );
  }
  return molids;
}

}} // cctbx::maptbx
