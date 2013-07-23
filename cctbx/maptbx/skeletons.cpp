#include <cstdio>
#include <string>
#include <limits>
#include <set>


// #include "util/push_disable_warnings.hxx"
#include <scitbx/error.h>
#include <cctbx/error.h>
#include <cctbx/uctbx.h>
//#include "util/pop_enable_warnings.hxx"


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

inline bool all_eq(const std::set<int> &s, int a)
{
  if( s.empty() )
    return false;
  std::set<int>::const_iterator b=s.begin(), e = s.end();
  --e;
  return (*b == a) && (*e == a);
}

inline bool all_ne(const std::set<int> &s, int a)
{
  for(std::set<int>::const_iterator i=s.begin(); i != s.end(); ++i)
  {
    if( *i == a )
      return false;
  }
  return true;
}


skeleton swanson(const_map_t &map, double sigma)
{
  double mean=0., esd=0.;
  for(std::size_t ii=0; ii<map.size(); ++ii)
  {
    double m = map[ii];
    mean += m;
    esd += m*m;
  }
  mean /= map.size();
  esd = esd/map.size() - mean*mean;
  double mapcutoff = mean + esd * sigma;

  std::vector<xyzm_t> xyzm;
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
          xyzm_t x({i,j,k},m);
          xyzm.push_back( x );
        }
      }
    }
  }
  std::sort(xyzm.begin(), xyzm.end(), gr);
  //[](const xyzm_t &a, const xyzm_t &b) {return (get<1>(a)) > (get<1>(b));}

  marks_t marks;
  marks.resize( ndim );
  marks.fill(0);
  std::size_t nmarks=0, nbonds=0;
  const unsigned short cube_size = 26;
  int3_t cube[cube_size] = {
    {-1,-1,-1},  {0,-1,-1},  {1,-1,-1},  {-1,0,-1},  {0,0,-1},  {1,0,-1},
        {-1,1,-1},  {0,1,-1},  {1,1,-1},
    {-1,-1,0},   {0,-1,0},   {1,-1,0},   {-1,0,0},
        {1,0,0},   {-1,1,0},   {0,1,0},   {1,1,0},
    {-1,-1,1},   {0,-1,1},   {1,-1,1},   {-1,0,1},   {0,0,1},   {1,0,1},
        {-1,1,1},   {0,1,1},   {1,1,1}
  };

  skeleton result;
  for(std::size_t jj=0; jj<xyzm.size(); ++jj)
  {
    std::set<int> featureset;
    //featureset.fill(-1);
    xyz_t x = get<0>(xyzm[jj]);
    for(unsigned short ic=0; ic<cube_size; ++ic)
    {
      xyz_t neighbor = x + cube[ic];
      if( any_lt(neighbor,0) || any_ge(neighbor,indim) )
        continue;
      int mark = marks(neighbor);
      featureset.insert(mark);
    }

    const unsigned fssize =featureset.size();
    CCTBX_ASSERT( fssize!=0U );
    if( all_eq(featureset, 0) )
    {
      // new maximum
      ++nmarks;
      marks(x) = nmarks;
      result.maximums.push_back(x);
    }
    else if( all_ne(featureset, 0) )
    {
      ; // local minimum will not be in the neighborhood of any point
    }
    else if( fssize == 2 )
    {
      // part of the growing nodule
      std::set<int>::const_iterator b = featureset.begin();
      if( *b != 0 )
        marks(x) = *b;
      else
        marks(x) = *(++b);
    }
    else if( fssize>2 )
    {
      for(std::set<int>::const_iterator i=featureset.begin();
          i!=featureset.end(); ++i)
      {
        int fi = *i;
        if( fi<=0 )
          continue;
        if( fi>nmarks )
          throw std::logic_error("wrong feautreset i");
        std::set<int>::const_iterator j=i;
        ++j;
        for( ; j!=featureset.end(); ++j)
        {
          int fj = *j;
          if( fj<=0 )
            continue;
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
          bool bnewjoing = true;
          result.joins.insert(join);
        }
      }
    }
    else
      throw std::logic_error("impossible");
  }
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
    molecules.push_back({imn,0});
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
