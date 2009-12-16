// Original implementation by
// Pavel Afonine, with speed optimizations by
// Ralf Grosse-Kunstleve.

#ifndef CCTBX_MASKS_AROUND_ATOMS_H
#define CCTBX_MASKS_AROUND_ATOMS_H

#include <cctbx/error.h>
#include <cctbx/uctbx.h>
#include <scitbx/math/modulo.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <map>
#include <vector>

namespace cctbx { namespace masks {

  namespace af = scitbx::af;

  template <typename DataType=int, typename FloatType=float>
  class around_atoms
  {
    public:
      typedef FloatType f_t;

      around_atoms(
        cctbx::uctbx::unit_cell const& unit_cell,
        std::size_t space_group_order_z,
        af::shared<scitbx::vec3<double> > const& sites_frac,
        af::shared<double> const& atom_radii,
        af::c_grid<3>::index_type const& gridding_n_real,
        f_t const& solvent_radius_,
        f_t const& shrink_truncation_radius_,
        bool explicit_distance_=false,
        bool debug_=false)
      :
        solvent_radius(solvent_radius_),
        shrink_truncation_radius(shrink_truncation_radius_),
        accessible_surface_fraction(-1),
        contact_surface_fraction(-1),
        debug(debug_),
        explicit_distance(explicit_distance_)
      {
        CCTBX_ASSERT(sites_frac.size() == atom_radii.size());
        CCTBX_ASSERT(solvent_radius >= 0);
        CCTBX_ASSERT(shrink_truncation_radius >= 0);
        CCTBX_ASSERT(gridding_n_real.const_ref().all_gt(0));
        data.resize(af::c_grid<3>(gridding_n_real), static_cast<DataType>(1));
        std::size_t n_solvent = compute_accessible_surface(
          unit_cell,
          space_group_order_z,
          sites_frac.const_ref(),
          atom_radii.const_ref());
        if( debug ) {
          n_atom_points = std::count(data.begin(), data.end(), 0);
          const size_t n1bar = std::count(data.begin(), data.end(), -1);
          const size_t n1 = std::count(data.begin(), data.end(), 1);
          CCTBX_ASSERT( n1 == n_solvent );
          CCTBX_ASSERT( n1 + n_atom_points + n1bar == data.size() );
        }
        else {
          n_atom_points = 0;
        }
        compute_contact_surface(
          unit_cell,
          space_group_order_z,
          n_solvent);
      }

      f_t solvent_radius;

      f_t shrink_truncation_radius;

      af::versa<DataType, af::c_grid<3> > data;

      //! Volume accessible for centers of solvent molecules.
      /*! Fraction of volume bounded by accessible surface.
          If the sphere with solvent_radius is rolled over the protein
          surface, the center of this sphere will be on the accessible
          surface.
       */
      double accessible_surface_fraction;

      /*! "contact_surf_fract" corresponds to the protein surface over
          which the sphere of radius solvent_radius is rolling.
          This is the contact
          surface. The relative volume bounded by the solvent is a ratio
          Vcontact/Vcell, where Vcontact is a volume encompassed by contact
          surface.
       */
      double contact_surface_fraction;

      size_t n_atom_points;  // for debugging purpose only

    protected:
      const bool debug;
      const bool explicit_distance; // for debugging purpose only

      static
      int
      ifloor(f_t const& x)
      {
        return scitbx::math::float_int_conversions<f_t, int>::ifloor(x);
      }

      static
      int
      iceil(f_t const& x)
      {
        return scitbx::math::float_int_conversions<f_t, int>::iceil(x);
      }

      static
      double
      approx_surface_fraction_under_symmetry(
        std::size_t n,
        std::size_t n_solvent,
        std::size_t space_group_order_z)
      {
        std::size_t n_non_solvent = (n - n_solvent)
                                  * space_group_order_z;
        if (n_non_solvent >= n) return 0;
        return static_cast<double>(n - n_non_solvent) / n;
      }

      std::size_t
      compute_accessible_surface(
        cctbx::uctbx::unit_cell const& unit_cell,
        std::size_t space_group_order_z,
        af::const_ref<scitbx::vec3<double> > const& sites_frac,
        af::const_ref<double> const& atom_radii)
      {
        // Severe code duplication: cctbx/maptbx/average_densities.h
        af::ref<DataType, af::c_grid<3> > data_ref = data.ref();
        std::size_t n_solvent = data_ref.size();
        int nx = static_cast<int>(data_ref.accessor()[0]);
        int ny = static_cast<int>(data_ref.accessor()[1]);
        int nz = static_cast<int>(data_ref.accessor()[2]);
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
          CCTBX_ASSERT( atmrad >= 0.0 );
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
              DataType* data_mxnymynz = &data_ref[(mxny + (*myi)) * nz];
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
                    CCTBX_ASSERT( dist_c>=0.0 );
                    CCTBX_ASSERT( std::fabs(dist-dist_c)<0.001 );
                  }
                }

                if (dist_c < cutoffsq) {
                  DataType& dr = data_mxnymynz[*mzi];
                  if (dr == 1) n_solvent--;
                  if (dist_c < radsq) dr =  0;
                  else if(dr!=0)    dr = -1;
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

      struct shrink_neighbors;
      friend struct shrink_neighbors;

      struct shrink_neighbors
      {
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

      void
      compute_contact_surface(
        cctbx::uctbx::unit_cell const& unit_cell,
        std::size_t space_group_order_z,
        std::size_t n_solvent)
      {
        af::ref<DataType, af::c_grid<3> > data_ref = data.ref();
        std::size_t data_size = data_ref.size();
        if (shrink_truncation_radius == 0) {
          for(std::size_t ilxyz=0;ilxyz<data_size;ilxyz++) {
            if (data_ref[ilxyz] == -1) data_ref[ilxyz] = 0;
          }
          contact_surface_fraction = accessible_surface_fraction;
          return;
        }
        int nx = static_cast<int>(data_ref.accessor()[0]);
        int ny = static_cast<int>(data_ref.accessor()[1]);
        int nz = static_cast<int>(data_ref.accessor()[2]);
        af::versa<DataType, af::c_grid<3> > datacopy = data.deep_copy();
        af::const_ref<DataType, af::c_grid<3> >
          datacopy_ref = datacopy.const_ref();
        shrink_neighbors neighbors(
          unit_cell,
          data_ref.accessor(),
          shrink_truncation_radius);
        const DataType* datacopy_ptr = datacopy.begin();
        for(std::size_t ilxyz=0;ilxyz<data_size;ilxyz++,datacopy_ptr++) {
          if (*datacopy_ptr == -1) {
            int ly = static_cast<int>(ilxyz / nz);
            int lz = static_cast<int>(ilxyz - ly * nz);
            int lx = ly / ny;
            ly -= lx * ny;
            typename shrink_neighbors::dim0::const_iterator
              tab_i_end = neighbors.table.end();
            for(typename shrink_neighbors::dim0::const_iterator
              tab_i = neighbors.table.begin();
              tab_i != tab_i_end;
              tab_i++) {
              int mx = lx + tab_i->first;
              while (mx >= nx) mx -= nx;
              int mxny = mx * ny;
              typename shrink_neighbors::dim1::const_iterator
                tab_j_end = tab_i->second.end();
              for(typename shrink_neighbors::dim1::const_iterator
                tab_j = tab_i->second.begin();
                tab_j != tab_j_end;
                tab_j++) {
                int my = ly + tab_j->first;
                while (my >= ny) my -= ny;
                const DataType*
                  datacopy_mxnymynz = &datacopy_ref[(mxny + my) * nz];
              typename shrink_neighbors::dim2::const_iterator
                tab_k_end = tab_j->second.end();
              for(typename shrink_neighbors::dim2::const_iterator
                tab_k = tab_j->second.begin();
                tab_k != tab_k_end;
                tab_k++) {
                  int mz = (*tab_k) + lz;
                  while (mz >= nz) mz -= nz;
                  if (datacopy_mxnymynz[mz] == 1) {
                    data_ref[ilxyz] = 1;
                    n_solvent++;
                    goto end_of_neighbors_loop;
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
  };

}} // namespace cctbx::masks

#endif // CCTBX_MASKS_AROUND_ATOMS_H
