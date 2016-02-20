#ifndef CCTBX_EXAMOKE_EXAMPLES_H
#define CCTBX_EXAMOKE_EXAMPLES_H
#include <scitbx/examples/bevington/prototype_core.h>
#include <scitbx/mat3.h>
#include <cctbx/miller.h>
#include <scitbx/vec3.h>

using std::size_t;

namespace cctbx{
namespace merging{
struct intensity_data {
  typedef scitbx::af::shared<std::size_t> veci;
  typedef scitbx::af::shared<double> vecd;
  typedef scitbx::af::shared<cctbx::miller::index<int> > milleri;

  veci frame;   // integer pointer to the frame or experiment (frame number)
  veci miller;  // integer pointer to the Miller index
  vecd raw_obs; // raw observation Bragg spot intensity
  vecd exp_var; // estimate of the measurement variance
  vecd stol_sq; // invariant array of (sin theta over lambda)**2
  milleri origHKL; // original index

  inline
  void
  reset_mem(){
     frame = veci(); miller = veci();
     raw_obs = vecd(); exp_var = vecd(); stol_sq = vecd();
  }

  inline
  boost::python::tuple
  estimate_G(const size_t& Nframe, const double& inv_d_sq_max, const double& inv_d_sq_min) const {
    /* simple weighted sum of observations to give the scale factor G */
    vecd result(Nframe);
    vecd weights(Nframe);
    veci visited(Nframe);
    for (size_t i = 0; i<raw_obs.size(); ++i){
      if ((inv_d_sq_max > 0.) && (inv_d_sq_max > 4.*stol_sq[i])) {continue;}
      if ((inv_d_sq_min > 0.) && (inv_d_sq_min < 4.*stol_sq[i])) {continue;}
      double weight = 1./exp_var[i];
      result[ frame[i] ] += weight * raw_obs[i];
      weights[ frame[i] ] += weight;
      visited[ frame[i] ] = 1;
    }
    for (size_t r = 0; r<result.size(); ++r){
      if (weights [r] > 0.) {
        result[ r ] /= weights [r];
      }
    }
    return boost::python::make_tuple(result,visited);
  }

  inline
  boost::python::tuple
  estimate_I(const size_t& Nmiller, const double& inv_d_sq_max, const double& inv_d_sq_min) const {
    /* simple weighted sum of observations to give the scale factor G */
    vecd result(Nmiller);
    vecd weights(Nmiller);
    veci visited(Nmiller);
    for (size_t i = 0; i<raw_obs.size(); ++i){
      if ((inv_d_sq_max > 0.) && (inv_d_sq_max > 4.*stol_sq[i])) {continue;}
      if ((inv_d_sq_min > 0.) && (inv_d_sq_min < 4.*stol_sq[i])) {continue;}
      double weight = 1./exp_var[i];
      result[ miller[i] ] += weight * raw_obs[i];
      weights[ miller[i] ] += weight;
      visited[ miller[i] ] = 1;
    }
    for (size_t r = 0; r<result.size(); ++r){
      if (weights [r] > 0.) {
        result[ r ] /= weights [r];
      }
    }
    return boost::python::make_tuple(result,visited);
  }

};

enum ParameterFlags {
  PartialityDeff = (1 << 0),     // Partiality is expressed in terms of mosaic block size only
  PartialityEtaDeff = (1 << 1), // Partiality is expressed in terms of both eta and Deff
  Bfactor = (1 << 2),             // Refine Bfactors for each lattice
  Deff = (1 << 3),                // Refine mosaic block size for each lattice
  Eta  = (1 << 4),                // Refine mosaic rotation for each lattice
  Rxy  = (1 << 5),                // Refine rotation angles on two axes perpendicular to the beam
};

struct scaling_common_functions {

    inline
    void set_cpp_data(intensity_data& f, int ni, int ng) // pass reference, avoid big copy
    {
        fsim = f;
        residuals = scitbx::af::shared<double>(fsim.raw_obs.size());
        weights = scitbx::af::shared<double>(fsim.raw_obs.size());
        for (int ix = 0; ix < fsim.raw_obs.size(); ++ix) {
          weights[ix] = 1./fsim.exp_var[ix];
        }
        N_I = ni; N_G = ng;
    }

    inline void set_parameter_flags(const int& b){
      bitflags = b;
    }
    inline void set_wavelength(scitbx::af::shared<double> b){
      wavelength = b;
    }
    inline void set_domain_size(scitbx::af::shared<double> b){
      domain_size = b;
    }
    inline void set_Astar_matrix(scitbx::af::shared<scitbx::mat3<double> > b){
      Astar_matrix = b;
    }

    scitbx::af::shared<double>
    get_rh_rs_ratio() const {
      scitbx::af::shared<double>result;
      for (int ix = 0; ix < fsim.raw_obs.size(); ++ix) {
        size_t frameno = fsim.frame[ix];
        double Rs  = 1./domain_size[ frameno ];
        scitbx::vec3<double> q = Astar_matrix[frameno] * fsim.origHKL[ix];
        double wave = wavelength[ frameno ];
        scitbx::vec3<double> s0 = scitbx::vec3<double>(0,0,-1/wave);
        double Rh  = (q+s0).length() - (1./wave);
        result.push_back(Rh/Rs);
      }
      return result;
    }

    void
    get_index_into_array(const scitbx::af::shared<double>& array){
      Iptr = array.begin();
      Gptr = Iptr + N_I;
      const double* last = Gptr;
      if ((bitflags & Bfactor) == Bfactor) {
        Bptr = last + N_G;
        last = Bptr;
      }
    }

    inline
    void
    fvec_callable(scitbx::af::shared<double> &current_values) {
      get_index_into_array(current_values);//indices into the array of parameters
      scitbx::af::shared<double> rh_rs;
      if ((bitflags & PartialityDeff) == PartialityDeff){
        rh_rs = get_rh_rs_ratio();
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

  protected:
    scitbx::af::shared<double> residuals, weights;
    intensity_data fsim;
    int N_I, N_G;
    int bitflags;
    scitbx::af::shared<double> wavelength;
    scitbx::af::shared<double> domain_size;
    scitbx::af::shared<scitbx::mat3<double> > Astar_matrix;
  protected:
    const double* Iptr;
    const double* Gptr;
    const double* Bptr;

};

class xscale6e: public scitbx::example::non_linear_ls_eigen_wrapper, public scaling_common_functions {
  public:
    xscale6e(int n_parameters):
      non_linear_ls_eigen_wrapper(n_parameters)
      {
      }

    void reset_mem(){
      fsim.reset_mem();
      eigen_wrapper.eigen_normal_matrix.resize(0,0);
      eigen_wrapper.eigen_normal_matrix = scitbx::example::linear_ls_eigen_wrapper::sparse_matrix_t();
      eigen_wrapper.solution_ = scitbx::af::shared<double>();
      eigen_wrapper.right_hand_side_ = scitbx::af::shared<double>();
    }

    void access_cpp_build_up_directly_eigen_eqn(bool objective_only, scitbx::af::shared<double> current_values) {
        fvec_callable(current_values);
        if (objective_only){
          add_residuals(residuals.const_ref(), weights.const_ref());
          return;
        }
        get_index_into_array(current_values);//indices into the array of parameters
        scitbx::af::shared<double> rh_rs;
        if ((bitflags & PartialityDeff) == PartialityDeff){
          rh_rs = get_rh_rs_ratio();
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
            Pitem -= rh_rs[ix]*rh_rs[ix];
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

          //add_equation(residuals[ix], jacobian_one_row.const_ref(), weights[ix]);
          add_residual(residuals[ix], weights[ix]);
          add_equation_eigen(-residuals[ix], jacobian_one_row_indices.const_ref(), jacobian_one_row_data.const_ref(), weights[ix]);
        }
    }
};

}}

#endif // CCTBX_EXAMOKE_EXAMPLES_H
