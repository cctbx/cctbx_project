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

#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>

#include <numeric>
#include <vector>


using namespace std;
namespace mmtbx { namespace pair_interaction {
namespace af=scitbx::af;
using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;
using scitbx::vec2;


//af::versa<double, af::c_grid<2> > hessian(vec3<double> const& distanceVector,
mat3<double> hessian(vec3<double> const& distanceVector,
                     double distanceReciprocal,
                     vec3<double> const& distanceUnitVector,
                     double fac1,
                     double fac2)
  {
    //af::versa<double, af::c_grid<2> > hessian_m(af::c_grid<2>(3,3), 0);
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
      hessian( mat3<FloatType>(0,0,0, 0,0,0, 0,0,0) )
    {}

    density_props(
      FloatType const& density_,
      vec3<FloatType> const& gradient_vector_,
      mat3<FloatType> const& hessian_
      )
    :
      density(density_),
      gradient_vector(gradient_vector_),
      hessian(hessian_)
    {}

    void add(density_props<> density_props_obj)
    {
      density += density_props_obj.density;
      gradient_vector += density_props_obj.gradient_vector;
      hessian += density_props_obj.hessian;
    }

    FloatType cal_silva()
    {
      FloatType silva = 0;
      std::size_t increment = 3;
      FloatType cntr = 0;
      for(std::size_t j=0; j < 9; j+=increment) {
        FloatType temp = 0;
        for(std::size_t i=0; i < 3; i++) {
          temp += gradient_vector[j] * hessian[i+j];
        }
        FloatType diff = (density * temp - gradient_vector[cntr] * gradient);
        silva += diff*diff;
        cntr += 1;
      }
      return silva;
    }

  bool has_silva_interaction()
    {
    if(density<0.0001) return false;
    FloatType silva = cal_silva();
    MMTBX_ASSERT(gradient != 0);
    silva *=(4.0 / (gradient*gradient*gradient));
    silva /=(1.0 + silva);
    if(silva>=0.9) return true;
    else           return false;
    }

};

template <typename FloatType=double>
class wfc
{
  public:
    af::shared<vec3<int> > node_offsets;
    af::shared<vec3<int> > coefficients_of_first_derivative;
    af::shared<vec3<int> > coefficients_of_second_derivative;
    FloatType              prefactor_of_first_derivative;
    FloatType              prefactor_of_second_derivative;
    FloatType              a;
    FloatType              b;
    FloatType              position_max;
    FloatType              square_position_max;
    int                    ngrid;
    af::shared<FloatType>  grid_positions;
    af::shared<FloatType>  grid_values;
    af::shared<FloatType>  first_derivative_of_grid_values;
    af::shared<FloatType>  second_derivative_of_grid_values;
    double                 core_cutdens;
    af::shared<vec3<double> > rr_array;
    af::shared<double>        r_array;
    double dx;
    int zz;


    wfc() {}

    wfc(af::shared<vec3<int> > const& node_offsets_,
        af::shared<vec3<int> > const& coefficients_of_first_derivative_,
        af::shared<vec3<int> > const& coefficients_of_second_derivative_,
        FloatType              const& prefactor_of_first_derivative_,
        FloatType              const& prefactor_of_second_derivative_,
        double                 const& core_cutdens_,
        af::shared<vec3<double> > const& rr_array_,
        af::shared<double>        const& r_array_,
        int const& ngrid_,
        double const& zz_
        )
    :
      node_offsets(node_offsets_),
      coefficients_of_first_derivative(coefficients_of_first_derivative_),
      coefficients_of_second_derivative(coefficients_of_second_derivative_),
      prefactor_of_first_derivative(prefactor_of_first_derivative_),
      prefactor_of_second_derivative(prefactor_of_second_derivative_),
      core_cutdens(core_cutdens_), rr_array(rr_array_), r_array(r_array_),
      ngrid(ngrid_), zz(zz_),
      dx(0.002) // From current data set (.wfc files)
    {
    first_derivative_of_grid_values = af::shared<FloatType>(ngrid);
    second_derivative_of_grid_values = af::shared<FloatType>(ngrid);
    grid_values = af::shared<FloatType>(ngrid);

    af::shared<vec3<int> > noef  = node_offsets;
    af::shared<vec3<int> > coef1 = coefficients_of_first_derivative;
    af::shared<vec3<int> > coef2 = coefficients_of_second_derivative;
    FloatType fac1  = prefactor_of_first_derivative;
    FloatType fac2  = prefactor_of_second_derivative;

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
      double PI = scitbx::constants::pi;
      grid_values[i] = rr_array[i][0] * r2 / (4.0 * PI);
      first_derivative_of_grid_values[i] = (rr_array[i][1] * delta - 2.0 * rr_array[i][0])* r3 / (4.0 * PI);
      second_derivative_of_grid_values[i] = (rr_array[i][2] * delta2- 5.0 * rr_array[i][1] * delta + 6.0 * rr_array[i][0])* r4 / (4.0 * PI);
    }

    grid_positions      = r_array;
    a                   = std::exp(-6.)/zz; // -6 is xmin in .wfc files
    b                   = dx;
    position_max        = r_array[ngrid-1];
    square_position_max = position_max*position_max;

  }

};

template <typename FloatType=double>
density_props<FloatType> atom_density_props(
  vec3<FloatType> const& p,
  vec3<FloatType> const& a_xyz,
  wfc<FloatType>  const& wfc_obj)
{
  vec3<FloatType> d_vector = a_xyz - p;
  FloatType dx = d_vector[0];
  FloatType dy = d_vector[1];
  FloatType dz = d_vector[2];
  FloatType d = std::sqrt(dx*dx + dy*dy + dz*dz); // norm
  if(d<1.e-10) d = 1.e-10;
  FloatType d_reciprocal = 1.0/d;
  vec3<FloatType> d_unit_vector=d_vector*d_reciprocal;
  FloatType f  = 0;
  FloatType fp = 0;
  FloatType fpp= 0;
  if(d<wfc_obj.position_max) {
    int ir=0;
    FloatType r=0;
    if(d<=wfc_obj.grid_positions[0]) {
      ir=1;
      r=wfc_obj.grid_positions[0];
    }
    else {
      ir=int(1+std::floor(std::log(d/wfc_obj.a)/wfc_obj.b));
      r=d;
    }
    af::tiny<FloatType, 4>  rr(0,0,0,0);
    af::tiny<FloatType, 4> dr1(0,0,0,0);
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
      FloatType prod = 1.0;
      for(std::size_t j=0; j < 4; j++) {
        if(i == j) continue;
        prod = prod * dr1[j] * x1dr12(i,j);
      }
      f   = f + wfc_obj.grid_values[ii] * prod;
      fp  = fp+wfc_obj.first_derivative_of_grid_values[ii]*prod;
      fpp = fpp+wfc_obj.second_derivative_of_grid_values[ii]*prod;
    }
  }
  return density_props<FloatType>(
    f,
    d_unit_vector*fp,
    hessian(d_vector, d_reciprocal,d_unit_vector,-fp,fpp)
    );
}

bool has_interaction_at_point(
  vec3<double> const&              p,
  af::shared<vec3<double> > const& a_xyz,
  af::shared<int> const&           element_flags,
  boost::python::list const& wfc_obj
  )
{
  density_props<double> density_props_obj = density_props<double>();
  for(std::size_t i=0; i < a_xyz.size(); i++) {
    vec3<double> d_vector = a_xyz[i] - p;
    double dx = d_vector[0];
    double dy = d_vector[1];
    double dz = d_vector[2];
    double d = std::sqrt(dx*dx + dy*dy + dz*dz); // norm
    if(d<15) {
      wfc<double> tmp = boost::python::extract<wfc<double> >(wfc_obj[element_flags[i]])();
      density_props_obj.add(atom_density_props(p, a_xyz[i], tmp));
    }
  }
  density_props_obj.density  = std::max(density_props_obj.density,1.0E-30);
  vec3<double> gv = density_props_obj.gradient_vector;
  double dot = gv[0]*gv[0]+gv[1]*gv[1]+gv[2]*gv[2];
  density_props_obj.gradient = dot;
  return density_props_obj.has_silva_interaction();
}

template <typename FloatType=double>
class point_and_pair
{
  public:
    vec3<double> point;
    int i;
    int j;

    point_and_pair() {}

    point_and_pair(vec3<double> const& point_, int const& i_, int const& j_)
    :
      point(point_), i(i_), j(j_)
    {}
};

//template <typename datatype=int>
af::shared<vec3<int> > points_and_pairs(
  vec3<int> const& ngrid,
  double const& step_size,
  af::shared<vec3<double> > const& xyz,
  vec3<double> const& xyz_min,
  af::shared<int> const& atom_in_residue,
  af::shared<int> const& element_flags,
  boost::python::list const& wfc_obj
  )
{
  af::shared<vec3<int> > interacting_pairs;
  for(std::size_t ix=0; ix < ngrid[0]; ix++) {
    for(std::size_t iy=0; iy < ngrid[1]; iy++) {
      for(std::size_t iz=0; iz < ngrid[2]; iz++) {
        vec3<double> point = vec3<double>(
          xyz_min[0]+ix*step_size,
          xyz_min[1]+iy*step_size,
          xyz_min[2]+iz*step_size);
        int first  = 999999;
        int second = 999999;
        int atom_id_1 = -1;
        int atom_id_2 = -1;
        for(std::size_t j=0; j < xyz.size(); j++) {
          vec3<double> diff = xyz[j] - point;
          double dist_sq = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
          if(dist_sq >= 200) continue;
          if(dist_sq < first) {
            second = first;
            first = dist_sq;
            atom_id_2 = atom_id_1;
            atom_id_1 = j;
          }
          else if(dist_sq < second and dist_sq != first) {
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
        bool has_interaction = has_interaction_at_point(
          point, xyz, element_flags, wfc_obj);
        if(has_interaction) {
          vec3<int> pair;
          if(ia1<ia2) pair = vec3<int>(ia1, ia2, 0);
          else        pair = vec3<int>(ia2, ia1, 0);
          interacting_pairs.push_back(pair);
        }
  }}}
  return interacting_pairs;
}



}} // namespace mmtbx::pair_interaction

#endif
