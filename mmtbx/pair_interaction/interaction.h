#ifndef MMTBX_PAIR_INTERACTION_H
#define MMTBX_PAIR_INTERACTION_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <vector>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/versa_algebra.h>
#include <cctbx/import_scitbx_af.h>
#include <cmath>
#include <cctbx/adptbx.h>
#include <cctbx/xray/scattering_type_registry.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <mmtbx/error.h>
#include <cctbx/xray/targets.h>
#include <scitbx/matrix/outer_product.h>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>

#include <numeric>
#include <vector>
#include <algorithm>

#include <scitbx/array_family/accessors/c_grid.h>


using namespace std;
namespace mmtbx { namespace pair_interaction {
namespace af=scitbx::af;
using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;
using scitbx::vec2;


mat3<double> hessian(vec3<double> const& distanceVector,
                     double distanceReciprocal,
                     vec3<double> const& distanceUnitVector,
                     double fac1,
                     double fac2)
  {
    mat3<double> hessian_m = mat3<double>(0,0,0, 0,0,0, 0,0,0);
    vec3<double> distanceUnitVector2 = vec3<double> (
      distanceUnitVector[0] * distanceUnitVector[0],
      distanceUnitVector[1] * distanceUnitVector[1],
      distanceUnitVector[2] * distanceUnitVector[2]);
    double distanceReciprocal2 = distanceReciprocal * distanceReciprocal;
    double fac11 = fac1 * distanceReciprocal;
    double fac3 = fac2 + fac11;
    std::size_t increment = 4;
    int cntr = 0;
    for(std::size_t j=0; j < 9; j+=increment) {
      hessian_m[j] = fac3 * distanceUnitVector2[cntr] - fac11;
      cntr+=1;
    }
    hessian_m[1]= distanceReciprocal2 * distanceVector[0] * distanceVector[1] * fac3;
    hessian_m[2]= distanceReciprocal2 * distanceVector[0] * distanceVector[2] * fac3;
    hessian_m[5]= distanceReciprocal2 * distanceVector[1] * distanceVector[2] * fac3;
    hessian_m[3]=hessian_m[1];
    hessian_m[6]=hessian_m[2];
    hessian_m[7]=hessian_m[5];
    return hessian_m;
  }

template <typename FloatType=double>
class density_props
{
  public:
    FloatType density;
    FloatType gradient; // XXX it is actually dot product!
    vec3<FloatType> gradient_vector;
    mat3<FloatType> hessian;

    density_props()
    :
      density(0),
      gradient_vector(vec3<FloatType>(0,0,0)),
      hessian( mat3<FloatType>(0,0,0, 0,0,0, 0,0,0) ),
      gradient(0)
    {}

    density_props(
      FloatType const& density_,
      vec3<FloatType> const& gradient_vector_,
      mat3<FloatType> const& hessian_
      )
    :
      density(density_),
      gradient_vector(gradient_vector_),
      hessian(hessian_),
      gradient(dot(gradient_vector_))
    {}

    void add(density_props<> density_props_obj)
    {
      density += density_props_obj.density;
      gradient_vector += density_props_obj.gradient_vector;
      hessian += density_props_obj.hessian;
      gradient = dot(gradient_vector);
    }

    FloatType dot(vec3<FloatType> v) {
      return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    }

    FloatType cal_silva()
    {
      FloatType silva = 0;
      std::size_t increment = 3;
      FloatType cntr = 0;
      for(std::size_t j=0; j < 9; j+=increment) {
        FloatType temp = 0;
        for(std::size_t i=0; i < 3; i++) {
          temp += gradient_vector[i] * hessian[i+j];
        }
        FloatType diff = (density * temp - gradient_vector[cntr] * gradient);
        silva += diff*diff;
        cntr += 1;
      }
      return silva;
    }

  double get_dori_value() {
    FloatType silva = cal_silva();
    silva *=(4.0 / (gradient*gradient*gradient));
    silva /=(1.0 + silva);
    return silva;
  }

  double get_sedd_value() {
    FloatType silva = cal_silva();
    silva *= (4.0 / (std::pow(density, 8)));
    silva = std::log((1.0 + silva));
    return silva;
  }

  bool has_silva_interaction(std::string const & silva_type)
    {
    if(silva_type.compare("dori")==0) {
      if(density<0.001) return false;
      if(std::abs(gradient) < 1.e-9) return false;
      FloatType dori = get_dori_value();
      if(dori>=0.8 && dori <= 1) { // This breaks mmtbx/pair_interaction/tst_02.py
      //if(silva>=0.9) {
        return true; // it is 0.8 in original Java version.
      }
      else return false;
    }
    else if(silva_type.compare("sedd")==0) {
      if(density<0.1) return false;
      if(get_sedd_value()<=5) return true;
      else                    return false;
    }
    else {
      return false;
    }
    }

};

class wfc
{
  public:
    af::shared<vec3<int> > node_offsets;
    af::shared<vec3<int> > coefficients_of_first_derivative;
    af::shared<vec3<int> > coefficients_of_second_derivative;
    double                 prefactor_of_first_derivative;
    double                 prefactor_of_second_derivative;
    double                 a;
    double                 b;
    double                 position_max;
    double                 square_position_max;
    int                    ngrid;
    af::shared<double>     grid_positions;
    af::shared<double>     grid_values;
    af::shared<double>     first_derivative_of_grid_values;
    af::shared<double>     second_derivative_of_grid_values;
    double                 core_cutdens;
    double dx;
    int zz;
    af::shared<double> r_array;
    boost::python::list wfcin_array;
    af::shared<double> occ_electrons;

    wfc() {}

    wfc(int const& ngrid_,
        double const& zz_,
        af::shared<double>        const& r_array_,
        boost::python::list const& wfcin_array_,
        af::shared<double> occ_electrons_
        )
    :
      prefactor_of_first_derivative(1.0 / 120.0),
      prefactor_of_second_derivative(2.0 / 120.0),
      core_cutdens(1E-12),
      ngrid(ngrid_), zz(zz_),
      dx(0.002), // From current data set (.wfc files)
      r_array(r_array_), wfcin_array(wfcin_array_),
      occ_electrons(occ_electrons_)
    {

    af::shared<af::shared<double> > wfcin_array_flex;
    for(std::size_t i=0;i<boost::python::len(wfcin_array);i++) {
      wfcin_array_flex.push_back(
        boost::python::extract<af::shared<double> >(wfcin_array[i])() );
    }

    node_offsets = af::shared<vec3<int> >(6);
    node_offsets[0] = vec3<int>( 0, -2, -5 );
    node_offsets[1] = vec3<int>( 1, -1, -4 );
    node_offsets[2] = vec3<int>( 2, 0, -3  );
    node_offsets[3] = vec3<int>( 3, 1, -2  );
    node_offsets[4] = vec3<int>( 4, 2, -1  );
    node_offsets[5] = vec3<int>(5,3,0      );

    coefficients_of_first_derivative = af::shared<vec3<int> >(6);
    coefficients_of_first_derivative[0] = vec3<int>( -274, 6, -24 );
    coefficients_of_first_derivative[1] = vec3<int>( 600, -60, 150 );
    coefficients_of_first_derivative[2] = vec3<int>( -600, -40, -400 );
    coefficients_of_first_derivative[3] = vec3<int>( 400, 120, 600 );
    coefficients_of_first_derivative[4] = vec3<int>( -150, -30, -600 );
    coefficients_of_first_derivative[5] = vec3<int>(24, 4, 274);

    coefficients_of_second_derivative = af::shared<vec3<int> >(6);
    coefficients_of_second_derivative[0] = vec3<int>( 225, -5, -50 );
    coefficients_of_second_derivative[1] = vec3<int>( -770, 80, 305 );
    coefficients_of_second_derivative[2] = vec3<int>( 1070, -150, -780);
    coefficients_of_second_derivative[3] = vec3<int>( -780, 80, 1070 );
    coefficients_of_second_derivative[4] = vec3<int>( 305, -5, -770 );
    coefficients_of_second_derivative[5] = vec3<int>( -50, 0, 225 );

    first_derivative_of_grid_values = af::shared<double>(ngrid);
    second_derivative_of_grid_values = af::shared<double>(ngrid);
    grid_values = af::shared<double>(ngrid);

    af::shared<vec3<int> > noef  = node_offsets;
    af::shared<vec3<int> > coef1 = coefficients_of_first_derivative;
    af::shared<vec3<int> > coef2 = coefficients_of_second_derivative;
    double fac1  = prefactor_of_first_derivative;
    double fac2  = prefactor_of_second_derivative;
    double PI = scitbx::constants::pi;
    af::shared<vec3<double> > rr_array = af::shared<vec3<double> >(
      ngrid, vec3<double>(0,0,0));
    for(std::size_t i=0; i < ngrid; i++) {
       af::shared<double> wfcin_array_flex_i = wfcin_array_flex[i];
       af::shared<double> wfcin_array_flex_i_sq =
         wfcin_array_flex_i * wfcin_array_flex_i;
       for(std::size_t k=0; k < wfcin_array_flex_i.size(); k++) {
         rr_array[i][0] += wfcin_array_flex_i_sq[k]*occ_electrons[k];
       }
       if(rr_array[i][0]/(4.0*PI*(r_array[i]*r_array[i])) < core_cutdens) {
        ngrid=i+1;
        break;
       }
     }
    for(std::size_t i=0; i < ngrid; i++) {
      int ic=1;
      if     (i<=1)       ic=0;
      else if(i>=ngrid-3) ic=2;
      for(std::size_t j=0; j < 6; j++) {
        rr_array[i][1] = rr_array[i][1] + coef1[j][ic]* rr_array[i + noef[j][ic]][0];
        rr_array[i][2] = rr_array[i][2] + coef2[j][ic]* rr_array[i + noef[j][ic]][0];
      }
      rr_array[i][1] = rr_array[i][1] * fac1;
      rr_array[i][2] = rr_array[i][2] * fac2;
      double r = r_array[i];
      MMTBX_ASSERT(r != 0);
      double r1 = 1.0 / r;
      double r2 = r1 * r1;
      double r3 = r2 * r1;
      double r4 = r3 * r1;
      double delta = 1.0 / dx;
      double delta2 = delta * delta;
      grid_values[i] = rr_array[i][0] * r2 / (4.0 * PI);
      first_derivative_of_grid_values[i] =
        (rr_array[i][1] * delta - 2.0 * rr_array[i][0])* r3 / (4.0 * PI);
      second_derivative_of_grid_values[i] = (rr_array[i][2] * delta2- 5.0 *
        rr_array[i][1] * delta + 6.0 * rr_array[i][0])* r4 / (4.0 * PI);
    }
    grid_positions      = r_array;
    a                   = std::exp(-6.)/zz; // -6 is xmin in .wfc files
    b                   = dx;
    position_max        = r_array[ngrid-1];
    square_position_max = position_max*position_max;
  }
};

density_props<double> atom_density_props(
  vec3<double> const& p,
  vec3<double> const& a_xyz,
  wfc          const& wfc_obj)
{
  vec3<double> d_vector = p - a_xyz;
  double dx = d_vector[0];
  double dy = d_vector[1];
  double dz = d_vector[2];
  double d = std::sqrt(dx*dx + dy*dy + dz*dz); // norm
  if(d<1.e-10) d = 1.e-10;
  double d_reciprocal = 1.0/d;
  vec3<double> d_unit_vector=d_vector*d_reciprocal;
  double f  = 0;
  double fp = 0;
  double fpp= 0;
  if(d<wfc_obj.position_max) {
    int ir=0;
    double r=0;
    if(d<=wfc_obj.grid_positions[0]) {
      ir=1;
      r=wfc_obj.grid_positions[0];
    }
    else {
      ir=int(1+std::floor(std::log(d/wfc_obj.a)/wfc_obj.b));
      r=d;
    }
    af::tiny<double, 4>  rr(0,0,0,0);
    af::tiny<double, 4> dr1(0,0,0,0);
    af::versa<double, af::c_grid<2> > x1dr12(af::c_grid<2>(4,4), 0);
    for(std::size_t i=0; i < 4; i++) {
      int ii = std::min(std::max(ir, 2), wfc_obj.ngrid) - 3 + i;
      rr[i] = wfc_obj.grid_positions[ii];
      dr1[i] = r - rr[i];
      for(std::size_t j=0; j < i; j++) {
        x1dr12(i,j) = 1.0 / (rr[i] - rr[j]);
        x1dr12(j,i) = -x1dr12(i,j);
      }
    }
    for(std::size_t i=0; i < 4; i++) {
      int ii = std::min(std::max(ir, 2), wfc_obj.ngrid) - 3 + i;
      double prod = 1.0;
      for(std::size_t j=0; j < 4; j++) {
        if(i == j) continue;
        prod = prod * dr1[j] * x1dr12(i,j);
      }
      f   = f + wfc_obj.grid_values[ii] * prod;
      fp  = fp+wfc_obj.first_derivative_of_grid_values[ii]*prod;
      fpp = fpp+wfc_obj.second_derivative_of_grid_values[ii]*prod;
    }
  }
  return density_props<double>(
    f,
    d_unit_vector*fp,
    hessian(d_vector, d_reciprocal,d_unit_vector,-fp,fpp)
    );
}

density_props<double> build_density_props_obj(
  vec3<double> const&                 p,
  af::const_ref<vec3<double> > const& a_xyz,
  af::const_ref<int> const&           element_flags,
  std::vector<wfc> const&             wfcs,
  std::vector<int> const&             candidates)
{
  // candidates: atoms in ascending order, a superset of those within 10 Bohr
  density_props<double> density_props_obj = density_props<double>();
  for(std::size_t c=0; c < candidates.size(); c++) {
    int i = candidates[c];
    vec3<double> diff = a_xyz[i] - p;
    double dist_sq = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
    if(dist_sq<100) {
      density_props_obj.add(
        atom_density_props(p, a_xyz[i], wfcs[element_flags[i]]));
    }
  }
  density_props_obj.density  = std::max(density_props_obj.density,1.0E-30);
  return density_props_obj;
}

/*
  Atoms binned in cubic cells of the cutoff length: for a point, the atoms
  of its cell and the 26 neighbouring cells are all the atoms within the
  cutoff of it. The candidates of a cell are those atoms, in ascending order,
  so that a loop over them visits the atoms in the same order as a loop over
  all atoms.
*/
class cell_list
{
  public:
    vec3<int> n_cells;
    double cell_size;
    vec3<double> origin;
    std::vector<std::vector<int> > atoms_in_cell;
    std::vector<std::vector<int> > candidates_of_cell;

    cell_list(
      af::const_ref<vec3<double> > const& xyz,
      vec3<double> const& origin_,
      double const& cell_size_)
    :
      cell_size(cell_size_), origin(origin_)
    {
      vec3<double> xyz_max = origin;
      for(std::size_t j=0; j < xyz.size(); j++) {
        for(std::size_t k=0; k < 3; k++) xyz_max[k] = std::max(xyz_max[k], xyz[j][k]);
      }
      for(std::size_t k=0; k < 3; k++) {
        n_cells[k] = int(std::floor((xyz_max[k]-origin[k])/cell_size))+1;
      }
      atoms_in_cell.resize(n_cells[0]*n_cells[1]*n_cells[2]);
      for(std::size_t j=0; j < xyz.size(); j++) {
        atoms_in_cell[cell_index(cell_of(xyz[j]))].push_back(j);
      }
      candidates_of_cell.resize(atoms_in_cell.size());
      for(int ix=0; ix < n_cells[0]; ix++) {
      for(int iy=0; iy < n_cells[1]; iy++) {
      for(int iz=0; iz < n_cells[2]; iz++) {
        std::vector<int>& c = candidates_of_cell[cell_index(vec3<int>(ix,iy,iz))];
        for(int jx=std::max(ix-1,0); jx <= std::min(ix+1,n_cells[0]-1); jx++) {
        for(int jy=std::max(iy-1,0); jy <= std::min(iy+1,n_cells[1]-1); jy++) {
        for(int jz=std::max(iz-1,0); jz <= std::min(iz+1,n_cells[2]-1); jz++) {
          std::vector<int> const& a = atoms_in_cell[cell_index(vec3<int>(jx,jy,jz))];
          c.insert(c.end(), a.begin(), a.end());
        }}}
        std::sort(c.begin(), c.end());
      }}}
    }

    vec3<int> cell_of(vec3<double> const& p) const
    {
      vec3<int> c;
      for(std::size_t k=0; k < 3; k++) {
        c[k] = int(std::floor((p[k]-origin[k])/cell_size));
        c[k] = std::min(std::max(c[k], 0), n_cells[k]-1);
      }
      return c;
    }

    int cell_index(vec3<int> const& c) const
    {
      return (c[0]*n_cells[1] + c[1])*n_cells[2] + c[2];
    }

    std::vector<int> const& candidates(vec3<double> const& p) const
    {
      return candidates_of_cell[cell_index(cell_of(p))];
    }
};

/*
  Residue pairs (1-based residue indices, as atom_in_residue) having a
  density interaction of type silva_type. Grid points step_size apart from
  xyz_min: a point is assigned to the pair of residues of its two nearest
  atoms (within 10 Bohr), and the pair interacts if the density at any of its
  points does.
  residue_in_core: indexed by residue index; when not empty, only pairs with
  at least one residue in the core are examined. The other pairs are never
  used by the callers, and a point farther than 10 Bohr from every core atom
  cannot have a core atom among its nearest two, so the grid is restricted
  to the core region (the grid points themselves do not move).
*/
af::shared<vec3<int> > points_and_pairs(
  vec3<int> const& ngrid,
  double const& step_size,
  af::shared<vec3<double> > const& xyz,
  vec3<double> const& xyz_min,
  af::shared<int> const& atom_in_residue,
  af::shared<int> const& element_flags,
  boost::python::list const& wfc_obj,
  std::string const & silva_type,
  af::shared<bool> const& residue_in_core
  )
{
  const double cutoff = 10.0;
  // Create 2x2 bool map with lengths being max residue number, init with false
  int max_air = af::max(atom_in_residue.const_ref());
  af::versa<bool, af::c_grid<2> > pairs;
  pairs.resize(af::c_grid<2>(max_air+1, max_air+1), false);
  bool use_core = residue_in_core.size() > 0;
  if(use_core) MMTBX_ASSERT(residue_in_core.size() > max_air);
  std::vector<wfc> wfcs;
  for(std::size_t i=0; i < boost::python::len(wfc_obj); i++) {
    wfcs.push_back(boost::python::extract<wfc>(wfc_obj[i])());
  }
  // Grid index ranges: the whole grid, or the core region
  vec3<int> lo(0,0,0);
  vec3<int> hi(ngrid[0]-1, ngrid[1]-1, ngrid[2]-1);
  if(use_core) {
    bool any = false;
    vec3<double> cmin(0,0,0), cmax(0,0,0);
    for(std::size_t j=0; j < xyz.size(); j++) {
      if(!residue_in_core[atom_in_residue[j]]) continue;
      if(!any) { cmin = xyz[j]; cmax = xyz[j]; any = true; }
      for(std::size_t k=0; k < 3; k++) {
        cmin[k] = std::min(cmin[k], xyz[j][k]);
        cmax[k] = std::max(cmax[k], xyz[j][k]);
      }
    }
    if(!any) return af::shared<vec3<int> >();
    for(std::size_t k=0; k < 3; k++) {
      lo[k] = std::max(lo[k],
        int(std::floor((cmin[k]-cutoff-xyz_min[k])/step_size)));
      hi[k] = std::min(hi[k],
        int(std::ceil((cmax[k]+cutoff-xyz_min[k])/step_size)));
    }
  }

  cell_list cells(xyz.const_ref(), xyz_min, cutoff);
  af::shared<vec3<int> > interacting_pairs;
  for(int ix=lo[0]; ix <= hi[0]; ix++) {
    double stx = xyz_min[0]+ix*step_size;
    for(int iy=lo[1]; iy <= hi[1]; iy++) {
      double sty = xyz_min[1]+iy*step_size;
      for(int iz=lo[2]; iz <= hi[2]; iz++) {
        vec3<double> point = vec3<double>(
          stx,
          sty,
          xyz_min[2]+iz*step_size);
        // int: the distances are compared truncated to whole Bohr^2 (atoms
        // closer than that count as equally near; kept as it always was)
        int first  = 999999;
        int second = 999999;
        int atom_id_1 = -1;
        int atom_id_2 = -1;
        std::vector<int> const& candidates = cells.candidates(point);
        for(std::size_t c=0; c < candidates.size(); c++) {
          int j = candidates[c];
          vec3<double> diff = xyz[j] - point;
          double dist_sq = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
          if(dist_sq >= cutoff*cutoff) continue;
          if(dist_sq < first) {
            second = first;
            first = dist_sq;
            atom_id_2 = atom_id_1;
            atom_id_1 = j;
          }
          else if(dist_sq < second && dist_sq != first) {
            if(atom_id_1 != atom_id_2) {
              second = dist_sq;
              atom_id_2 = j;
            }
          }
        }
        if(atom_id_1 == -1) continue;
        if(atom_id_2 == -1) continue;
        if(atom_in_residue[atom_id_1] == atom_in_residue[atom_id_2]) continue;
        int ia1 = atom_in_residue[atom_id_1];
        int ia2 = atom_in_residue[atom_id_2];
        if(use_core && !residue_in_core[ia1] && !residue_in_core[ia2]) continue;

        vec3<int> pair;
        if(ia1<ia2) pair = vec3<int>(ia1, ia2, 0);
        else        pair = vec3<int>(ia2, ia1, 0);
        if(pairs(pair[0], pair[1]) == true) continue;

        density_props<double> density_props_obj = build_density_props_obj(
          point, xyz.const_ref(), element_flags.const_ref(), wfcs, candidates);
        bool has_interaction = density_props_obj.has_silva_interaction(silva_type);

        if(has_interaction) {
          pairs(pair[0], pair[1])=true;
          interacting_pairs.push_back(pair);
        }
  }}}
  return interacting_pairs;
}

}} // namespace mmtbx::pair_interaction

#endif
