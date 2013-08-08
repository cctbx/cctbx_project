#include <cstdio>
#include <string>
#include <limits>
#include <set>


#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>

#include <cctbx/error.h>

//// DO NOT USE. UNDER CONSTRUCTION.

#include "cctbx/maptbx/skeletons.h"

namespace cctbx { namespace maptbx
{

inline bool gr(const xyzm_t &a, const xyzm_t &b)
{
  return (get<1>(a)) > (get<1>(b));
}

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

inline bool all_eq(const std::set<std::size_t> &s, std::size_t a)
{
#if 0
  if( s.empty() )
    return false;
  std::set<std::size_t>::const_iterator b=s.begin(), e = s.end();
  --e;
  return (*b == a) && (*e == a);
#endif
  return boost::algorithm::all_of_equal(s,a);
}

inline bool all_ne(const std::set<std::size_t> &s, std::size_t a)
{
#if 0
  for(std::set<std::size_t>::const_iterator i=s.begin(); i != s.end(); ++i)
  {
    if( *i == a )
      return false;
  }
  return true;
#endif
  return !boost::algorithm::any_of_equal(s,a);
}


skeleton swanson(const_map_t &map, double sigma)
{
  double mean=0., esd=0., mx=-9.E200;
  for(std::size_t ii=0; ii<map.size(); ++ii)
  {
    double m = map[ii];
    mean += m;
    esd += m*m;
    if( m > mx )
      mx = m;
  }
  mean /= map.size();
  esd = esd/map.size() - mean*mean;
  esd = std::sqrt(esd);
  double mapcutoff = mean + esd * sigma;

  std::vector<xyzm_t> xyzm;
  xyzm.reserve(10000);
  const scitbx::af::tiny<std::size_t,3> ndim( map.accessor().focus() );
  int3_t indim(ndim);
  for(std::size_t i=0; i<ndim[0]; ++i)
  {
    for(std::size_t j=0; j<ndim[1]; ++j)
    {
      for(std::size_t k=0; k<ndim[2]; ++k)
      {
        double m = map(i,j,k);
        if( m>mapcutoff )
        {
          int3_t p(i,j,k);
          xyzm_t x(p,m);
          xyzm.push_back( x );
        }
      }
    }
  }
  std::sort(xyzm.begin(), xyzm.end(), gr);
  //[](const xyzm_t &a, const xyzm_t &b) {return (get<1>(a)) > (get<1>(b));}

  CCTBX_ASSERT( get<1>(xyzm.front()) == mx );
  marks_t marks(ndim, 0UL); // 0 means no mark
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

  skeleton result;
  for(std::size_t jj=0; jj<xyzm.size(); ++jj)
  {
    //! @todo very inefficient, consider boost::pool_allocator
    std::set<std::size_t> featureset;
    int3_t x = get<0>(xyzm[jj]);
    for(unsigned short ic=0; ic<cube_size; ++ic)
    {
      int3_t neighbor = x + cube[ic];
      if( any_lt(neighbor,0) || any_ge(neighbor,indim) )
        continue;
      std::size_t mark = marks(neighbor);
      if( mark != 0U ) // 0 is not a mark
        featureset.insert(mark);
    }

    const unsigned fssize =featureset.size();
    if( featureset.empty() ) // all_eq(featureset, 0) )
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
    else if( fssize>1 )
    {
      ++join_count;
      // 2c: two or more nodules merging
      for(std::set<std::size_t>::const_iterator i=featureset.begin();
          i!=featureset.end(); ++i)
      {
        std::size_t fi = *i;
        CCTBX_ASSERT( fi>0U && fi<=nmarks );
        std::set<std::size_t>::const_iterator j=i;
        ++j;
        for( ; j!=featureset.end(); ++j)
        {
          std::size_t fj = *j;
          CCTBX_ASSERT( fj>0U );
          if( fj>nmarks )
            throw std::logic_error("wrong featureset j");
          join_t join;
          if( fi<fj )
          {
            join.ilt = fi;
            join.igt = fj;
          }
          else
          {
            join.ilt = fj;
            join.igt = fi;
          }
          //! @todo search existing joins
          bool bnewjoing = true;
          result.joins.insert(join);
        }
      }
    }
    else if( fssize == cube_size ) // all_ne(featureset, 0) )
    {
      // 2d: local minimum will not be in the neighborhood of any point
      ++min_count;
    }
    else
      throw std::logic_error("impossible");
  }
  CCTBX_ASSERT( nmarks + grows_count + join_count + min_count == xyzm.size() );
  result.min_count = min_count;
  result.grows_count = grows_count;
  result.join_count = join_count;
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
    //CCTBX_ASSERT( i_cat < i_an && i_cat>=0 && i_an>=0 && i_cat<nmaxs
    //  && i_an<nmaxs );
    if( atoms[i_cat]==0 && atoms[i_an]==0 )
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
