#ifndef CCTBX_EXAMOKE_EXAMPLES_H
#define CCTBX_EXAMOKE_EXAMPLES_H
#include <scitbx/examples/bevington/prototype_core.h>

using std::size_t;

namespace cctbx{
namespace merging{
struct intensity_data {
  typedef scitbx::af::shared<std::size_t> veci;
  typedef scitbx::af::shared<double> vecd;

  veci frame;   // integer pointer to the frame or experiment (frame number)
  veci miller;  // integer pointer to the Miller index
  vecd raw_obs; // raw observation Bragg spot intensity
  vecd exp_var; // estimate of the measurement variance
  vecd stol_sq; // invariant array of (sin theta over lambda)**2

  inline
  void
  reset_mem(){
     frame = veci(); miller = veci();
     raw_obs = vecd(); exp_var = vecd(); stol_sq = vecd();
  }

  inline
  boost::python::tuple
  estimate_G(const size_t& Nframe) const {
    /* simple weighted sum of observations to give the scale factor G */
    vecd result(Nframe);
    vecd weights(Nframe);
    veci visited(Nframe);
    for (size_t i = 0; i<raw_obs.size(); ++i){
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
  estimate_I(const size_t& Nmiller) const {
    /* simple weighted sum of observations to give the scale factor G */
    vecd result(Nmiller);
    vecd weights(Nmiller);
    veci visited(Nmiller);
    for (size_t i = 0; i<raw_obs.size(); ++i){
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

    inline
    void
    fvec_callable(scitbx::af::shared<double> &current_values) {
      const double* Iptr = current_values.begin(); //indices into the array of parameters
      const double* Gptr = Iptr + N_I;
      const double* Bptr = Gptr + N_G;
      for (int ix = 0; ix < fsim.raw_obs.size(); ++ix) {
        double Gitem  = Gptr[ fsim.frame[ix] ];
        double Bitem  = std::exp(-2.* Bptr[ fsim.frame[ix] ] * fsim.stol_sq[ix]); // :=exp(beta)
        double Iitem  = Iptr[ fsim.miller[ix] ];

        residuals[ix] = fsim.raw_obs[ix] - Gitem * Bitem * Iitem;
      }
    }

  protected:
    scitbx::af::shared<double> residuals, weights;
    intensity_data fsim;
    int N_I, N_G;
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
        const double* Iptr = current_values.begin(); //indices into the array of parameters
        const double* Gptr = Iptr + N_I;
        const double* Bptr = Gptr + N_G;
        // add one of the normal equations per each observation
        for (int ix = 0; ix < fsim.raw_obs.size(); ++ix) {

          double Gitem  = Gptr[ fsim.frame[ix] ];
          double Bitem  = std::exp(-2.* Bptr[ fsim.frame[ix] ] * fsim.stol_sq[ix]); // :=exp(beta)
          double Iitem  = Iptr[ fsim.miller[ix] ];

          scitbx::af::shared<std::size_t> jacobian_one_row_indices;
          scitbx::af::shared<double> jacobian_one_row_data;

          //derivative with respect to I ... must be indexed first
          jacobian_one_row_indices.push_back( fsim.miller[ix] );
          jacobian_one_row_data.push_back(-Bitem * Gitem);

          //derivative with respect to G ... must be indexed second
          jacobian_one_row_indices.push_back(N_I + fsim.frame[ix]);
          jacobian_one_row_data.push_back(-Bitem * Iitem);

          //derivative with respect to B ... must be indexed third
          jacobian_one_row_indices.push_back(N_I + N_G + fsim.frame[ix]);
          jacobian_one_row_data.push_back(Bitem * Gitem * Iitem * 2. * fsim.stol_sq[ix]);

          //add_equation(residuals[ix], jacobian_one_row.const_ref(), weights[ix]);
          add_residual(residuals[ix], weights[ix]);
          add_equation_eigen(-residuals[ix], jacobian_one_row_indices.const_ref(), jacobian_one_row_data.const_ref(), weights[ix]);
        }
        set_from_triplets();
        wipe_triplets();
    }
};

}}

#endif // CCTBX_EXAMOKE_EXAMPLES_H
