#ifndef CCTBX_EXAMOKE_EXAMPLES_PRB_H
#define CCTBX_EXAMOKE_EXAMPLES_PRB_H
#include <cctbx/examples/merging/xscale_prototype_core.h>
#include <scitbx/constants.h>

using std::size_t;

namespace cctbx{
namespace merging{

class postrefine_base: public xscale6e {
  typedef scitbx::vec3<double> vec3;
  typedef scitbx::mat3<double> mat3;
  typedef scitbx::af::shared<vec3> vec3vec;
  typedef scitbx::af::shared<mat3> mat3vec;

  public:
    postrefine_base(int n_parameters):
      Iptr(0),Gptr(0),Bptr(0),Dptr(0),Axptr(0),Ayptr(0),
      xscale6e(n_parameters)
      {
      }

    const double* Iptr;
    const double* Gptr;
    const double* Bptr;
    const double* Dptr;
    const double* Axptr;
    const double* Ayptr;
    mat3vec Rx;
    mat3vec Ry;
    mat3vec dRx_dT;
    mat3vec dRy_dT;
    vec3vec s0array;
    scitbx::af::shared<double> Rh_Rs_ratio;

    void
    get_index_into_array(const scitbx::af::shared<double>& array){
      Iptr = array.begin();
      Gptr = Iptr + N_I;
      const double* last = Gptr;
      if ((bitflags & Bfactor) == Bfactor) {
        Bptr = last + N_G;
        last = Bptr;
      }
      if ((bitflags & Deff) == Deff) {
        Dptr = last + N_G;
        last = Dptr;
      }
      if ((bitflags & Rxy) == Rxy) {
        Axptr = last + N_G;
        Ayptr = last + 2*N_G;
        last = Ayptr;
      }
    }

    void
    set_rotational_derivatives() {
      //rotational operators for rotation around x, then y
      Rx = mat3vec(N_G);
      Ry = mat3vec(N_G);
      dRx_dT = mat3vec(N_G);
      dRy_dT = mat3vec(N_G);
      s0array = vec3vec(N_G);
      mat3 I3 = mat3(1,0,0,0,1,0,0,0,1); //identity

      vec3 unit_axis = vec3(1.,0.,0.);
      mat3 OP = scitbx::mat3<double>::outer_product(unit_axis,unit_axis);
      mat3 CP = scitbx::mat3<double>::cross_product_matrix(unit_axis);
      for (int idxf = 0; idxf < N_G; ++idxf){
        double c = std::cos(scitbx::constants::pi_180*Axptr[idxf]);
        double s = std::sin(scitbx::constants::pi_180*Axptr[idxf]);
        Rx[idxf]     = (  I3*c + OP*(1.-c) + CP*s );
        dRx_dT[idxf] = ( -I3*s + OP*s + CP*c )/scitbx::constants::pi_180;
      }
      unit_axis = vec3(0.,1.,0.);
      OP = scitbx::mat3<double>::outer_product(unit_axis,unit_axis);
      CP = scitbx::mat3<double>::cross_product_matrix(unit_axis);
      for (int idxf = 0; idxf < N_G; ++idxf){
        double c = std::cos(scitbx::constants::pi_180*Ayptr[idxf]);
        double s = std::sin(scitbx::constants::pi_180*Ayptr[idxf]);
        Ry[idxf]     = (  I3*c + OP*(1.-c) + CP*s );
        dRy_dT[idxf] = ( -I3*s + OP*s + CP*c )/scitbx::constants::pi_180;
        double wave = wavelength[ idxf ];
        s0array[idxf] = scitbx::vec3<double>(0,0,-1/wave);
      }
    }

    scitbx::af::shared<double>
    get_rh_rs_ratio(const scitbx::af::shared<double>& array) {
      get_index_into_array(array);//indices into the array of parameters
      scitbx::af::shared<double>result;
      for (int ix = 0; ix < fsim.raw_obs.size(); ++ix) {
        size_t frameno = fsim.frame[ix];
        double Rs  = (Dptr==0) ? 1./domain_size[ frameno ] : 1./Dptr[frameno];
        scitbx::vec3<double> q = Astar_matrix[frameno] * fsim.origHKL[ix];
        if (Axptr!=0) { q = Rx[frameno] * Ry[frameno] * q;} // q is rotated by small rotations about x & y
        double wave = wavelength[ frameno ];
        scitbx::vec3<double> s0 = scitbx::vec3<double>(0,0,-1/wave);
        double Rh  = (q+s0).length() - (1./wave);
        result.push_back(Rh/Rs);
      }
      Rh_Rs_ratio = result;
      return result;
    }

    inline
    void
    fvec_callable(scitbx::af::shared<double> &current_values) {
      // might be slower to get indices this way, gain is that it handles arbitrary parameter combinations
      get_index_into_array(current_values);//indices into the array of parameters
      scitbx::af::shared<double> rh_rs;
      if ((bitflags & PartialityDeff) == PartialityDeff){
        rh_rs = get_rh_rs_ratio(current_values);
      }
      for (int ix = 0; ix < fsim.raw_obs.size(); ++ix) {
        double Gitem  = Gptr[ fsim.frame[ix] ];
        double Bitem  = 1.;
        if ((bitflags & Bfactor) == Bfactor){ //note the necessary parentheses
          Bitem  = std::exp(-2.* Bptr[ fsim.frame[ix] ] * fsim.stol_sq[ix]); // :=exp(beta)
        }
        double Iitem  = Iptr[ fsim.miller[ix] ];

        double Pitem  = 1.;
        if ((bitflags & PartialityDeff) == PartialityDeff){
          Pitem -= rh_rs[ix]*rh_rs[ix];
        }
        residuals[ix] = fsim.raw_obs[ix] - Pitem * Gitem * Bitem * Iitem;
      }
    }

    void access_cpp_build_up_directly_eigen_eqn(bool objective_only, scitbx::af::shared<double> current_values) {
        get_index_into_array(current_values);//indices into the array of parameters
        if ((bitflags & Rxy) == Rxy){ set_rotational_derivatives(); }
        fvec_callable(current_values); // also calls Rh_rs ratio
        if (objective_only){
          add_residuals(residuals.const_ref(), weights.const_ref());
          return;
        }
        // add one of the normal equations per each observation
        for (int ix = 0; ix < fsim.raw_obs.size(); ++ix) {

          double Gitem  = Gptr[ fsim.frame[ix] ];
          double Bitem  = 1.;
          if ((bitflags & Bfactor) == Bfactor){ //note the necessary parentheses
          Bitem  = std::exp(-2.* Bptr[ fsim.frame[ix] ] * fsim.stol_sq[ix]); // :=exp(beta)
          }
          double Iitem  = Iptr[ fsim.miller[ix] ];
          double Pitem  = 1.;
          if ((bitflags & PartialityDeff) == PartialityDeff){
            Pitem -= Rh_Rs_ratio[ix]*Rh_Rs_ratio[ix];
          }

          scitbx::af::shared<std::size_t> jacobian_one_row_indices;
          scitbx::af::shared<double> jacobian_one_row_data;

          //derivative with respect to I ... must be indexed first
          jacobian_one_row_indices.push_back( fsim.miller[ix] );
          jacobian_one_row_data.push_back(-Bitem * Gitem * Pitem);

          //derivative with respect to G ... must be indexed second
          std::size_t jidx = N_I + fsim.frame[ix];
          jacobian_one_row_indices.push_back( jidx );
          jacobian_one_row_data.push_back(-Bitem * Iitem * Pitem);

          //derivative with respect to B ... must be indexed third
          if ((bitflags & Bfactor) == Bfactor){ //note the necessary parentheses
          jidx += N_G;
          jacobian_one_row_indices.push_back( jidx );
          jacobian_one_row_data.push_back(Bitem * Gitem * Iitem * Pitem * 2. * fsim.stol_sq[ix]);
          }

          //derivative with respect to Deff ... must be indexed fourth
          if ((bitflags & Deff) == Deff){ //note the necessary parentheses
          jidx += N_G;
          jacobian_one_row_indices.push_back( jidx );
          if ((bitflags & PartialityDeff) == PartialityDeff) {
            size_t frameno = fsim.frame[ix];
            scitbx::vec3<double> q = Astar_matrix[frameno] * fsim.origHKL[ix];
            double wave = wavelength[ frameno ];
            scitbx::vec3<double> s0 = scitbx::vec3<double>(0,0,-1/wave);
            double Rh  = (q+s0).length() - (1./wave);
            jacobian_one_row_data.push_back(2. * Bitem * Gitem * Iitem * Dptr[frameno] * Rh * Rh);
          } else {
            jacobian_one_row_data.push_back(0.);
          }
          }

          //derivatives with respect to Ax, Ay ... must be indexed 5th, 6th
          if ((bitflags & Rxy) == Rxy){ //note the necessary parentheses
            jidx += N_G;
            jacobian_one_row_indices.push_back( jidx );
            jidx += N_G;
            jacobian_one_row_indices.push_back( jidx );
            if ((bitflags & PartialityDeff) == PartialityDeff) {
              size_t frameno = fsim.frame[ix];
              scitbx::vec3<double> q = Astar_matrix[frameno] * fsim.origHKL[ix];
              vec3 qprime = Rx[frameno] * Ry[frameno] * q; // q is rotated by small rotations about x & y
              vec3 s1 = qprime+s0array[frameno];
              double dRh_dAx = (s1/s1.length()) * (dRx_dT[frameno] * Ry[frameno] * q);
              double dRh_dAy = (s1/s1.length()) * (Rx[frameno] * dRy_dT[frameno] * q);
              jacobian_one_row_data.push_back(2. * Bitem * Gitem * Iitem * Rh_Rs_ratio[ix] * dRh_dAx);
              jacobian_one_row_data.push_back(2. * Bitem * Gitem * Iitem * Rh_Rs_ratio[ix] * dRh_dAy);
            } else {
              jacobian_one_row_data.push_back(0.);
              jacobian_one_row_data.push_back(0.);
            }
          }

          //add_equation(residuals[ix], jacobian_one_row.const_ref(), weights[ix]);
          add_residual(residuals[ix], weights[ix]);
          add_equation_eigen(-residuals[ix],
              jacobian_one_row_indices.const_ref(), jacobian_one_row_data.const_ref(), weights[ix]);

        }
    }
};

}}

#endif // CCTBX_EXAMOKE_EXAMPLES_PRB_H
