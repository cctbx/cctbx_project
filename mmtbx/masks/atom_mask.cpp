#include <numeric>
#include <algorithm>

#include <scitbx/array_family/flex_types.h>
#include <cctbx/maptbx/structure_factors.h>
#include <scitbx/fftpack/real_to_complex_3d.h>
#include <cctbx/maptbx/gridding.h>
#include <scitbx/array_family/tiny_types.h>
#include <cctbx/sgtbx/direct_space_asu/proto/small_vec_math.h>

#include "atom_mask.h"

namespace mmtbx { namespace masks {

  using namespace cctbx::sgtbx::asu;
  typedef double f_t;

  namespace {

    inline int ifloor(f_t const& x)
    {
      return scitbx::math::float_int_conversions<f_t, int>::ifloor(x);
    }

    inline int iceil(f_t const& x)
    {
      return scitbx::math::float_int_conversions<f_t, int>::iceil(x);
    }

    inline double approx_surface_fraction_under_symmetry(
      std::size_t n,
      std::size_t n_solvent,
      std::size_t space_group_order_z)
    {
      std::size_t n_non_solvent = (n - n_solvent)
                                * space_group_order_z;
      if (n_non_solvent >= n) return 0.0;
      return static_cast<double>(n - n_non_solvent) / n;
    }

  }


  unsigned short site_symmetry_order(const cctbx::sgtbx::space_group &group, const rvector3_t &v, scitbx::double3 &delta)
  {
    unsigned short nops = 0;
    scitbx::double3 dv(boost::rational_cast<double,int>(v[0]),
      boost::rational_cast<double,int>(v[1]), boost::rational_cast<double,int>(v[2]));
    dv -= floor(dv);
    for(size_t i=0; i<group.order_z(); ++i)
    {
      const cctbx::sgtbx::rt_mx op = group(i);
      scitbx::double3 sv0 = op * dv;
      scitbx::double3 sv1 = sv0 - scitbx::floor(sv0);
      // scitbx::double3 sv = scitbx::abs(dv - sv1);
      scitbx::double3 sv( std::fabs(dv[0]-sv1[0]), std::fabs(dv[1]-sv1[1]), std::fabs(dv[2]-sv1[2]) );
      for(unsigned char j=0; j<3; ++j)
      {
        MMTBX_ASSERT( sv[j]<=1.0 );
        MMTBX_ASSERT( sv[j]>=0.0 );
        MMTBX_ASSERT( delta[j]>=0.0 );
        if( sv[j] > 1.0-delta[j] )
          sv[j] = std::fabs( sv[j]-1.0 );
      }
      if( scitbx::le_all(sv , delta) )
        ++nops;
    }
    MMTBX_ASSERT( nops>0U );
    return nops;
  }

  void atom_mask::mask_asu()
  {
    unsigned short order = group.order_z();
    MMTBX_ASSERT( order>0 );
    const af::ref<data_type, grid_t > data_ref = data.ref();
    const scitbx::vec3<size_t> n_ = data_ref.accessor();
    const scitbx::vec3<int> n(n_[0], n_[1], n_[2]);
    cctbx::sgtbx::asu::rvector3_t mn, mx;
    MMTBX_ASSERT( n[0]>0 && n[1]>0 && n[2] >0 );
    scitbx::double3 tolerance( 0.05/n[0], 0.05/n[1], 0.05/n[2]);
    asu.box_corners(mn,mx);
    mul(mn, n);
    mul(mx, n);
    cctbx::sg_vec3 imx = ceil(mx), imn = floor(mn);
    MMTBX_ASSERT( imx[0]>imn[0] && imx[1] > imn[1] && imx[2] > imn[2] );
    MMTBX_ASSERT( imx[0] <= n[0]+1 && imx[1] <= n[1] + 1 && imx[2] <= n[2] + 1 ); // this should go away

    for(long i=imn[0]-1; i<=imx[0]+1 && i<=n[0]+1; ++i)
    {
      const rational_t x(i,n[0]);
      long i_c = i % n[0];
      if( i_c<0 )
        i_c += n[0];
      for(long j=imn[1]-1; j<=imx[1]+1 && j<=n[1]+1; ++j)
      {
        const rational_t y(j,n[1]);
        long  j_c = j % n[1];
        if( j_c<0 )
          j_c += n[1];
        for(long k=imn[2]-1; k<=imx[2]+1 && k<=n[2]+1; ++k)
        {
          const rational_t z(k,n[2]);
          long  k_c = k % n[2];
          if( k_c<0 )
            k_c += n[2];
          MMTBX_ASSERT( i_c<n[0] && j_c<n[1] && k_c<n[2] );
          const rvector3_t pos(x,y,z);

          if( asu.is_inside(pos) )
          {
            unsigned short nops = site_symmetry_order(group, pos, tolerance);
            MMTBX_ASSERT( nops>0 );
            MMTBX_ASSERT( order%nops == 0);
            nops = order / nops;
            MMTBX_ASSERT( data_ref(i_c,j_c,k_c)==0 || nops == data_ref(i_c,j_c,k_c) );
            data_ref(i_c,j_c,k_c) = nops;
          }
          //else
            //data_ref(i_c,j_c,k_c ) = 0;
        }
      }
    }
  }

  void atom_mask::compute(
    const coord_array_t & sites_frac,
    const double_array_t & atom_radii)
  {
    this->atoms_to_asu(sites_frac, atom_radii);
    this->mask_asu();
    size_t n1 = this->get_mask().size() - std::count( this->get_mask().begin(), this->get_mask().end(), 0 );
    MMTBX_ASSERT( n1>0 );
    size_t nn = 0;
    nn = std::accumulate( this->get_mask().begin(), this->get_mask().end(), nn );
    MMTBX_ASSERT( nn > 0 );
    MMTBX_ASSERT( nn == this->get_mask().size() );
    size_t n_solvent = this->compute_accessible_surface(asu_atoms, asu_radii);
    size_t tmp=0;
    tmp = std::count( get_mask().begin(), get_mask().end(), 0);
    tmp = std::count_if( get_mask().begin(), get_mask().end(), std::bind2nd(std::less<data_type>(),0) );
    size_t nn_solv = 0;
    nn_solv = std::accumulate( this->get_mask().begin(), this->get_mask().end(), nn_solv );
    this->compute_contact_surface(n_solvent);
    nn_solv = 0U;
    nn_solv = std::accumulate( this->get_mask().begin(), this->get_mask().end(), nn_solv );
    contact_surface_fraction = accessible_surface_fraction = double(nn_solv)/nn;
    tmp = std::count( get_mask().begin(), get_mask().end(), 0);
    tmp = std::count_if( get_mask().begin(), get_mask().end(), std::bind2nd(std::less<data_type>(),0) );
    MMTBX_ASSERT( tmp==0 );
  }

  void atom_mask::atoms_to_asu(
    const coord_array_t & sites_frac,
    const double_array_t & atom_radii)
  {
    MMTBX_ASSERT(sites_frac.size() == atom_radii.size());
    asu_atoms.clear();
    asu_radii.clear();
    const cctbx::sgtbx::space_group &grp = group;
    const scitbx::af::tiny<double,6> rcell = cell.reciprocal_parameters();
    const scitbx::double3 rp(rcell[0], rcell[1], rcell[2]);

    std::vector< scitbx::tiny3 > cells;
    this->asu.get_adjacent_cells(cells);

    for(size_t iat=0; iat<sites_frac.size(); ++iat)
    {
      scitbx::double3 at = sites_frac[iat];
      const double at_r =  atom_radii[iat];
      const double radius = at_r + solvent_radius;
      MMTBX_ASSERT( radius >= 0.0 );
      scitbx::double3 box = rp * radius;
      box *= 1.05;

      for(size_t isym=0; isym<grp.order_z(); ++isym)
      {
        const cctbx::sgtbx::rt_mx symop = grp(isym);
        scitbx::double3 sym_at = symop*at;
        sym_at -= scitbx::floor(sym_at);
        for(size_t icell=0; icell<cells.size(); ++icell)
        {
          const scitbx::double3 sym_at_cell = sym_at + cells[icell];
          const cctbx::sgtbx::asu::intersection_kind  intersection = asu.does_intersect(sym_at_cell, box);
          if( intersection != cctbx::sgtbx::asu::none )
          {
            this->asu_atoms.push_back(sym_at);
            this->asu_radii.push_back(at_r);
            if( intersection == cctbx::sgtbx::asu::fully )
              goto end_sym_loop; // need to break out of the symmetry loop here
            break;
          }
        }
      }
      end_sym_loop:
        ;
    }
  }


  typedef scitbx::af::c_grid_padded<3> padded_grid_t;
  typedef scitbx::af::versa<double, padded_grid_t > versa_3d_padded_real_array;
  typedef scitbx::af::versa<std::complex<double>, padded_grid_t > versa_3d_padded_complex_array;


  scitbx::af::shared< std::complex<double> > atom_mask::structure_factors(
      const scitbx::af::const_ref< cctbx::miller::index<> > &indices ) const
  {
    const mask_array_t &msk = this->get_mask();
    const mmtbx::masks::grid_t grid = msk.accessor();
    scitbx::fftpack::real_to_complex_3d<double> fft(grid);

    // m_real : physical dims, n_real - focus dims
    // assert( n_real <= m_real )
    const cctbx::sg_vec3 mdim = fft.m_real(), ndim = fft.n_real();
    const padded_grid_t pad( mdim, ndim );
    versa_3d_padded_real_array padded_real(pad, 0.0);
    MMTBX_ASSERT( padded_real.size() >= msk.size() );
    // convert non-padded integer mask to padded real array
    for(size_t i=0; i< ndim[0]; ++i)
      for( size_t j=0; j<ndim[1]; ++j)
        for( size_t k=0; k<ndim[2]; ++k)
          padded_real(i,j,k) = msk(i,j,k);

    fft.forward(padded_real); // in-place forward FFT
    const padded_grid_t pad_complex( fft.n_complex(), fft.n_complex() );
    versa_3d_padded_complex_array result(padded_real.handle(), pad_complex );

    double scale = cell.volume() / ( ndim.product() * group.order_z() ) ;
    // result *= scale;
    for(size_t i=0; i<result.size(); ++i)
      result[i] *= scale;

    const cctbx::maptbx::structure_factors::from_map<double>  the_from_map (
      group,
      false, // anomalous flag
      indices,
      result.const_ref(),
      true); // conjugate_flag
    return the_from_map.data();
  }


  void atom_mask::determine_gridding(cctbx::sg_vec3 &grid, double resolution, double factor, int use_symmetry ) const
  {
    MMTBX_ASSERT( factor > 0.0 );
    double step = resolution/factor;
    if(step < 0.15)
      step = 0.15;
    step = std::min(0.8, step);
    const double d_min = 2.0*step, resolution_factor = 0.5;
    if( use_symmetry == 2 )
    {
      const cctbx::sgtbx::search_symmetry_flags use_all(true,
          0,
          true,
          true,
          true);
      grid = cctbx::maptbx::determine_gridding<int>(cell, d_min, resolution_factor, use_all, group.type());
      grid = group.refine_gridding( grid );
    }
    else if( use_symmetry == 1 )
    {
      const cctbx::sgtbx::search_symmetry_flags use_all(true,
          0,
          true,
          true,
          true);
      grid = cctbx::maptbx::determine_gridding<int>(cell, d_min, resolution_factor, use_all, group.type());
    }
    else if( use_symmetry == 3 )
    {
      grid = cctbx::maptbx::determine_gridding<int>(cell, d_min, resolution_factor);
      grid = group.refine_gridding( grid );
    }
    else
    {
      grid = cctbx::maptbx::determine_gridding<int>(cell, d_min, resolution_factor);
    }
  }

  size_t atom_mask::compute_accessible_surface(
    const coord_array_t & sites_frac,
    const double_array_t & atom_radii)
  {
    MMTBX_ASSERT( sites_frac.size() == atom_radii.size() );
    cctbx::uctbx::unit_cell const& unit_cell = this->cell;
    const size_t space_group_order_z = this->group_order_z;

    // Severe code duplication: cctbx/maptbx/average_densities.h
    af::ref<data_type, af::c_grid<3> > data_ref = data.ref();
    std::size_t n_solvent = data_ref.size();
    int nx = static_cast<int>(data_ref.accessor()[0]);
    int ny = static_cast<int>(data_ref.accessor()[1]);
    int nz = static_cast<int>(data_ref.accessor()[2]);
    MMTBX_ASSERT( nx>0 && ny>0 && nz>0 );
    f_t mr1= static_cast<f_t>(unit_cell.metrical_matrix()[0]); // a*a
    f_t mr5= static_cast<f_t>(unit_cell.metrical_matrix()[1]); // b*b
    f_t mr9= static_cast<f_t>(unit_cell.metrical_matrix()[2]); // c*c
    // a*b*cos(gamma)
    f_t mr2= static_cast<f_t>(unit_cell.metrical_matrix()[3]);
    // a*c*cos(beta)
    f_t mr3= static_cast<f_t>(unit_cell.metrical_matrix()[4]);
    // c*b*cos(alpha)
    f_t mr6= static_cast<f_t>(unit_cell.metrical_matrix()[5]);
    f_t tmr2 = mr2*2; //2*a*b*cos(gamma);
    f_t tmr3 = mr3*2; //2*a*c*cos(beta);
    f_t tmr6 = mr6*2; //2*b*c*cos(alpha);
    f_t sx = 1/static_cast<f_t>(nx); f_t tsx= sx*2; f_t sxsq=mr1*sx*sx;
    f_t sy = 1/static_cast<f_t>(ny); f_t tsy= sy*2; f_t sysq=mr5*sy*sy;
    f_t sz = 1/static_cast<f_t>(nz); f_t tsz= sz*2; f_t szsq=mr9*sz*sz;
    f_t w1=mr1*sx*tsx; f_t w4=mr5*sy*tsy;
    f_t w2=mr2*sx*tsy; f_t w5=mr6*sy*tsz;
    f_t w3=mr3*sx*tsz; f_t w6=mr9*sz*tsz;
    f_t tsxg1=tsx*mr1; f_t tsyg4=tsy*mr2; f_t tszg3=tsz*mr3;
    f_t tsxg4=tsx*mr2; f_t tsyg5=tsy*mr5; f_t tszg8=tsz*mr6;
    f_t tsxg7=tsx*mr3; f_t tsyg8=tsy*mr6; f_t tszg9=tsz*mr9;
    f_t rp[3];
    for(unsigned i=0;i<3;i++) {
      rp[i] = static_cast<f_t>(unit_cell.reciprocal_parameters()[i]);
    }
    std::vector<int> mys;
    std::vector<int> mzs;
    for(std::size_t i_site=0;i_site<sites_frac.size();i_site++) {
      cctbx::fractional<> const& site = sites_frac[i_site];
      f_t xfi=static_cast<f_t>(site[0]);
      f_t yfi=static_cast<f_t>(site[1]);
      f_t zfi=static_cast<f_t>(site[2]);
      const f_t atmrad = atom_radii[i_site];
      MMTBX_ASSERT( atmrad >= 0.0 );
      f_t cutoff=static_cast<f_t>(atmrad+solvent_radius);
      f_t radsq=static_cast<f_t>(atmrad*atmrad);
      f_t cutoffsq=cutoff*cutoff;
      f_t coas = cutoff*rp[0];
      int x1box=ifloor(nx*(xfi-coas));
      int x2box=iceil(nx*(xfi+coas));
      f_t cobs = cutoff*rp[1];
      int y1box=ifloor(ny*(yfi-cobs));
      int y2box=iceil(ny*(yfi+cobs));
      f_t cocs = cutoff*rp[2];
      int z1box=ifloor(nz*(zfi-cocs));
      int z2box=iceil(nz*(zfi+cocs));
      f_t sxbcen=xfi-x1box*sx;
      f_t sybcen=yfi-y1box*sy;
      f_t szbcen=zfi-z1box*sz;
      f_t distsm=mr1*sxbcen*sxbcen+mr5*sybcen*sybcen+mr9*szbcen*szbcen
            +tmr2*sxbcen*sybcen+tmr3*sxbcen*szbcen+tmr6*sybcen*szbcen;
      f_t w7=tsxg1*sxbcen+tsxg4*sybcen+tsxg7*szbcen;
      f_t w8=tsyg4*sxbcen+tsyg5*sybcen+tsyg8*szbcen;
      f_t w9=tszg3*sxbcen+tszg8*sybcen+tszg9*szbcen;
      mys.clear();
      mys.reserve(y2box-y1box+1);
      for (int ky = y1box; ky <= y2box; ky++) {
        int my = ky % ny;
        if (my < 0) my += ny;
        mys.push_back(my);
      }
      mzs.clear();
      mzs.reserve(z2box-z1box+1);
      for (int kz = z1box; kz <= z2box; kz++) {
        int mz = kz % nz;
        if (mz < 0) mz += nz;
        mzs.push_back(mz);
      }
      f_t distsx = distsm;
      f_t s1xx = sxsq - w7;
      f_t s1xy = sysq - w8;
      f_t s1xz = szsq - w9;
      for (int kx = x1box; kx <= x2box; kx++) {
        int mx = kx % nx;
        if (mx < 0) mx += nx;
        int mxny = mx * ny;
        f_t s2yz = s1xz;
        f_t s2_incr = s1xy;
        f_t s2 = distsx;
        std::vector<int>::const_iterator mye = mys.end();
        for (std::vector<int>::const_iterator
               myi=mys.begin();
               myi!=mye;
               myi++) {
          f_t s3_incr = s2yz;
          f_t dist = s2;
          data_type* data_mxnymynz = &data_ref[(mxny + (*myi)) * nz];
          std::vector<int>::const_iterator mze = mzs.end();
          for (std::vector<int>::const_iterator
                 mzi=mzs.begin();
                 mzi!=mze;
                 mzi++) {
            f_t dist_c  = dist;
            if( explicit_distance ) {
              // neglect previous stepwise distance calculation
              // use formula instead
              const f_t dx = xfi-kx*sx;
              const f_t dy = yfi-(myi-mys.begin()+y1box)*sy;
              const f_t dz = zfi-(mzi-mzs.begin()+z1box)*sz;
              dist_c = mr1*dx*dx+mr5*dy*dy+mr9*dz*dz
                    +tmr2*dx*dy+tmr3*dx*dz+tmr6*dy*dz;
              if( debug ) {
                MMTBX_ASSERT( dist_c>=0.0 );
                MMTBX_ASSERT( std::fabs(dist-dist_c)<0.001 );
              }
            }

            if (dist_c < cutoffsq) {
              data_type& dr = data_mxnymynz[*mzi];
              if (dr > 0 ) n_solvent--;
              if (dist_c < radsq) dr =  0;
              else if(dr>0)
              {
                dr = -dr;
                MMTBX_ASSERT( dr < 0 );
              }
            }
            dist += s3_incr;
            s3_incr += w6;
          }
          s2 += s2_incr;
          s2_incr += w4;
          s2yz += w5;
        }
        distsx += s1xx;
        s1xx += w1;
        s1xy += w2;
        s1xz += w3;
      }
    }
    accessible_surface_fraction = approx_surface_fraction_under_symmetry(
      data_ref.size(), n_solvent, space_group_order_z);
    return n_solvent;
  }


  struct shrink_neighbors
  {
    typedef double f_t;
    shrink_neighbors() {}

    shrink_neighbors(
      cctbx::uctbx::unit_cell const& unit_cell,
      af::c_grid<3>::index_type const& gridding_n_real,
      f_t const& shrink_truncation_radius)
    {
      int low[3];
      int high[3];
      for(unsigned i=0;i<3;i++) {
        double x = shrink_truncation_radius
                 * unit_cell.reciprocal_parameters()[i]
                 * gridding_n_real[i];
        low[i] = ifloor(-x);
        high[i] = iceil(x);
      }
      int n0 = static_cast<int>(gridding_n_real[0]);
      int n1 = static_cast<int>(gridding_n_real[1]);
      int n2 = static_cast<int>(gridding_n_real[2]);
      f_t shrink_truncation_radius_sq = shrink_truncation_radius
                                      * shrink_truncation_radius;
      cctbx::fractional<f_t> frac;
      for(int p0=low[0];p0<=high[0];p0++) {
        int m0 = scitbx::math::mod_positive(p0, n0);
        frac[0] = static_cast<f_t>(p0) / n0;
      for(int p1=low[1];p1<=high[1];p1++) {
        int m1 = scitbx::math::mod_positive(p1, n1);
        frac[1] = static_cast<f_t>(p1) / n1;
      for(int p2=low[2];p2<=high[2];p2++) {
        frac[2] = static_cast<f_t>(p2) / n2;
        f_t dist_sq = unit_cell.length_sq(frac);
        if (dist_sq < shrink_truncation_radius_sq) {
          int m2 = scitbx::math::mod_positive(p2, n2);
          table[m0][m1].push_back(m2);
        }
      }}}
    }

    typedef std::vector<int> dim2;
    typedef std::map<int, dim2> dim1;
    typedef std::map<int, dim1> dim0;
    dim0 table;
  };

  scitbx::int3 closest_grid_point(const rvector3_t &pos, const scitbx::int3 &grid_size)
  {
    rvector3_t tmp(pos);
    scitbx::mul(tmp, grid_size);
    scitbx::int3 result = iround( tmp );
    return result;
  }

  inline scitbx::int3 closest_grid_point(const scitbx::double3 &pos, const scitbx::int3 &grid_size)
  {
    scitbx::double3 tmp(pos);
    scitbx::mul(tmp, grid_size);
    return scitbx::vec3_cast<int,double>( scitbx::round( tmp ) );
  }


  typedef af::const_ref< cctbx::sgtbx::rt_mx > symop_array;

  void atom_mask::compute_contact_surface( std::size_t n_solvent)
  {
    typedef double f_t;
    typedef data_type data_type;
    cctbx::uctbx::unit_cell const& unit_cell = this->cell;
    const std::size_t space_group_order_z = this->group_order_z;
    // const symop_array sym_ops = group.all_ops().const_ref();
    const af::shared< cctbx::sgtbx::rt_mx > sym_ops_ = group.all_ops();
    const symop_array sym_ops = sym_ops_.const_ref();
    MMTBX_ASSERT( sym_ops.size() == group.order_z() );

    af::ref<data_type, af::c_grid<3> > data_ref = data.ref();
    std::size_t data_size = data_ref.size();
    if (shrink_truncation_radius == 0) {
      for(std::size_t ilxyz=0;ilxyz<data_size;ilxyz++) {
        if (data_ref[ilxyz] < 0 ) data_ref[ilxyz] = 0;
      }
      contact_surface_fraction = accessible_surface_fraction;
      return;
    }
    int nx = static_cast<int>(data_ref.accessor()[0]);
    int ny = static_cast<int>(data_ref.accessor()[1]);
    int nz = static_cast<int>(data_ref.accessor()[2]);
    const scitbx::int3 grid_sz(nx,ny,nz);
    af::versa<data_type, af::c_grid<3> > datacopy = data.deep_copy();
    af::const_ref<data_type, af::c_grid<3> >
      datacopy_ref = datacopy.const_ref();
    shrink_neighbors neighbors(
      unit_cell,
      data_ref.accessor(),
      shrink_truncation_radius);
    const data_type* datacopy_ptr = datacopy.begin();
    for(std::size_t ilxyz=0;ilxyz<data_size;ilxyz++,datacopy_ptr++) {
      if (*datacopy_ptr < 0 ) {
        int ly = static_cast<int>(ilxyz / nz);
        int lz = static_cast<int>(ilxyz - ly * nz);
        int lx = ly / ny;
        ly -= lx * ny;
        shrink_neighbors::dim0::const_iterator
          tab_i_end = neighbors.table.end();
        for(shrink_neighbors::dim0::const_iterator
          tab_i = neighbors.table.begin();
          tab_i != tab_i_end;
          tab_i++) {
          int mx = lx + tab_i->first;
          const double dmx = static_cast<double>(mx)/nx;
          shrink_neighbors::dim1::const_iterator
            tab_j_end = tab_i->second.end();
          for(shrink_neighbors::dim1::const_iterator
            tab_j = tab_i->second.begin();
            tab_j != tab_j_end;
            tab_j++) {
            int my = ly + tab_j->first;
            const double dmy = static_cast<double>(my)/ny;
          shrink_neighbors::dim2::const_iterator
            tab_k_end = tab_j->second.end();
          for(shrink_neighbors::dim2::const_iterator
            tab_k = tab_j->second.begin();
            tab_k != tab_k_end;
            tab_k++) {
              int mz = (*tab_k) + lz;
              const scitbx::double3 d_tst( dmx, dmy, static_cast<double>(mz)/nz );
              for(symop_array::const_iterator op=sym_ops.begin(); op!=sym_ops.end(); ++op)
              {
                // const cctbx::sgtbx::rt_mx &op = sym_ops[i_sym];
                scitbx::double3 sv = (*op) * d_tst;
                sv -= scitbx::floor(sv);
                scitbx::int3 tst = closest_grid_point(sv, grid_sz);
                if( tst[0] >= nx ) tst[0] -= nx;
                if( tst[1] >= ny ) tst[1] -= ny;
                if( tst[2] >= nz ) tst[2] -= nz;
                if( datacopy_ref( tst[0], tst[1], tst[2]) > 0 )
                {
                  data_ref[ilxyz] = -data_ref[ilxyz]; // 1;
                  MMTBX_ASSERT( data_ref[ilxyz] > 0 );
                  n_solvent++;
                  goto end_of_neighbors_loop;
                }
              }
            }
          }
        }
        data_ref[ilxyz] = 0;
        end_of_neighbors_loop:;
      }
    }
    contact_surface_fraction = approx_surface_fraction_under_symmetry(
      data_ref.size(), n_solvent, space_group_order_z);
  }

}} // namespace mmtbx::masks
