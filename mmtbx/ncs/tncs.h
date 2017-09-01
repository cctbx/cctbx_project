#ifndef MMTBX_TNCS_H
#define MMTBX_TNCS_H

#include <cctbx/import_scitbx_af.h>
#include <boost/python/list.hpp>
#include <scitbx/sym_mat3.h>
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/miller.h>
#include <scitbx/math/g_function.h>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
// #include <cctbx/maptbx/real_space_gradients_simple.h> // import dependency
#include <cctbx/sgtbx/space_group.h>
#include <mmtbx/error.h>


using namespace std;
namespace mmtbx { namespace ncs {
namespace af=scitbx::af;
namespace bp=boost::python;
using scitbx::vec3;
using scitbx::mat3;

/*
Class-container of parameters for NCS pair.
For usage from Python see /cctbx_project/mmtbx/ncs/tst_ncs_pair.py
*/

template <typename FloatType=double>
class pair
{
  public:
    scitbx::mat3<FloatType> r;
    scitbx::vec3<FloatType> t;
    FloatType radius;
    FloatType radius_estimate;
    FloatType fracscat;
    af::shared<FloatType> rho_mn;
    int id;
    pair() {
      r.fill(0);
      t.fill(0);
      radius=0;
      radius_estimate=0;
      fracscat=0;
      rho_mn.fill(0);
      id=-1;
    }
    pair(scitbx::mat3<FloatType> const& r_,
         scitbx::vec3<FloatType> const& t_,
         FloatType radius_,
         FloatType radius_estimate_,
         FloatType fracscat_,
         af::shared<FloatType> rho_mn_,
         int id_)
    :
      r(r_),t(t_),radius(radius_),radius_estimate(radius_estimate_),
      fracscat(fracscat_), rho_mn(rho_mn_), id(id_)
    {}

    void set_radius(FloatType _) { radius = _; }
    void set_rhoMN(af::shared<FloatType> _) { rho_mn = _; }
    void set_id(int _) { id = _; }
};

// --------------------------------------------------------------------------

template <typename FloatType>
class tncs_eps_factor_refinery { // Gfunction in Phaser
public:
  af::shared<af::double6> GfunTensorArray;
  af::shared<scitbx::vec3<FloatType> > ncsDeltaT;
  af::shared<FloatType> dEps_by_drho,dEps_by_dradius;
  FloatType refl_epsfac;
  af::shared<pair<FloatType> > pairs;
  int dim;
  int n_pairs;
  scitbx::mat3<FloatType> fractionalization_matrix;
  af::shared<scitbx::mat3<FloatType> > sym_mat;
  af::shared<FloatType> tncs_epsfac;
  cctbx::sgtbx::space_group space_group;
  af::shared<bool> centric_flags;
  af::shared<int>  epsilon;
  af::const_ref<cctbx::miller::index<int> > miller_indices;
  af::const_ref<FloatType> f_obs;
  af::const_ref<FloatType> sig_f_obs;
  af::const_ref<FloatType> SigmaN;
  af::const_ref<int> rbin;
  // resulting gradient
  af::shared<FloatType> Gradient_rhoMN;
  af::shared<FloatType> Gradient_radius;
  bool do_radius;
  bool do_rho_mn;
  bool target_called;

  tncs_eps_factor_refinery(
    bp::list const& pairs_,
    af::const_ref<FloatType> const& f_obs_,
    af::const_ref<FloatType> const& sig_f_obs_,
    af::const_ref<int> const& rbin_,
    af::const_ref<FloatType> const& SigmaN_,
    cctbx::sgtbx::space_group space_group_,
    af::const_ref<cctbx::miller::index<int> > miller_indices_,
    scitbx::mat3<FloatType> fractionalization_matrix_,
    af::shared<scitbx::mat3<FloatType> > sym_mat_)
  :
  refl_epsfac(0),
  sym_mat(sym_mat_),
  miller_indices(miller_indices_),
  space_group(space_group_),
  fractionalization_matrix(fractionalization_matrix_),
  f_obs(f_obs_),
  sig_f_obs(sig_f_obs_),
  rbin(rbin_),
  SigmaN(SigmaN_),
  do_radius(false),
  do_rho_mn(false),
  target_called(false)
  {
    MMTBX_ASSERT(f_obs.size()==sig_f_obs.size());
    MMTBX_ASSERT(f_obs.size()==miller_indices.size());
    MMTBX_ASSERT(SigmaN.size()==f_obs.size());
    int n_bins = af::max(rbin)+1;
    for(std::size_t i=0;i<bp::len(pairs_);i++) {
      pairs.push_back(bp::extract<pair<FloatType> >(pairs_[i])());
      MMTBX_ASSERT(pairs[i].rho_mn.size()==n_bins);
    }
    n_pairs = pairs.size();
    dim = n_pairs * sym_mat.size();
    GfunTensorArray.resize(dim,af::double6(0,0,0,0,0,0));
    ncsDeltaT.resize(dim,scitbx::vec3<FloatType>(0,0,0));
    centric_flags = space_group.is_centric(miller_indices);
    epsilon = space_group.epsilon(miller_indices);
    calcArrays();
  }

  void set_compute_gradients_radius() {
    do_radius = true;
    do_rho_mn = false;
  }

  void set_compute_gradients_rho_mn() {
    do_radius = false;
    do_rho_mn = true;
  }

  void update_pairs(bp::list const& new_pairs) {
    for(std::size_t i=0;i<bp::len(new_pairs);i++) {
      pairs[i] = bp::extract<pair<FloatType> >(new_pairs[i])();
    }
  }

  void calcArrays()
  {
  /*
   Prepare array of metric tensors needed to quickly calculate the G-functions
   for tNCS. This code assumes that the G-function corresponding to a spherical
   envelope will be used. Handling non-spherical envelopes would require
   interpolating a finely-sampled transform instead.
  */
    scitbx::mat3<FloatType> identity(1,0,0,0,1,0,0,0,1);
    // Pre-compute some results that don't depend on reflection number
    for(int ipair = 0; ipair < n_pairs; ipair++) {
      FloatType radSqr = scitbx::fn::pow2(pairs[ipair].radius);
      scitbx::mat3<FloatType> ncsR_t = pairs[ipair].r.transpose();
      scitbx::vec3<FloatType> ncsT = pairs[ipair].t;
      for (int isym = 0; isym < sym_mat.size(); isym++) {
  /*
    Metric matrix (represented as tensor) allows quick computation of
    |s_klmm|^2 and therefore rs^2 for G-function
    If F=fractionalization_matrix, R=Rotsym and (T) indicates transpose, then
    s_klmm = (identity-ncsR(T))*F(T)*R(T)*h
    |s_klmm|^2 = s_klmm(T).s_klmm =
                         h(T)*R*F*(identity-ncsR)*(identity-ncsR(T))*F(T)*R(T)*h
    The metric matrix is the part between h(T) and h.  In the
    GfunTensorArray, factors of 2 are included so
    the unique terms need to be computed only once.
  */
        int i = ipair*sym_mat.size()+isym;
        ncsDeltaT[i] = sym_mat[isym].transpose()*ncsT; // symmetry translation cancels
        scitbx::mat3<FloatType> metric_matrix =
          sym_mat[isym].transpose()*fractionalization_matrix*(identity-ncsR_t);
        metric_matrix = radSqr*metric_matrix*metric_matrix.transpose();
        GfunTensorArray[i] = af::double6(
            metric_matrix(0,0),  metric_matrix(1,1),  metric_matrix(2,2),
          2*metric_matrix(0,1),2*metric_matrix(0,2),2*metric_matrix(1,2));
      }
    }
  }

  void calcRefineTerms(cctbx::miller::index<int> const& miller, int sbin)
  {
    /* Calculate tNCS-related epsilon factor (to multiply by normal
       symmetry-related epsilon factor), plus derivatives that can be used to
       refine parameters characterising the tNCS
    */
    if(do_rho_mn) dEps_by_drho.resize(n_pairs);
    if(do_radius) dEps_by_dradius.resize(n_pairs);
    refl_epsfac = 1.; // Epsilon before modulation by tNCS
    /*Use half-triple summation rather than full quadruple summation, i.e. only
      paying attention to terms from molecules that are in similar orientations.
    */
    for (int ipair = 0; ipair < n_pairs; ipair++) {
      double wtfac(pairs[ipair].fracscat/sym_mat.size());
      if(do_rho_mn) dEps_by_drho[ipair] = 0.;
      if(do_radius) dEps_by_dradius[ipair] = 0.;
      for (int isym = 0; isym < sym_mat.size(); isym++) {
        // NCS-related contribution to epsilon, weighted by interference,
        // G-function and rhoMN terms
        int hh(miller[0]), kk(miller[1]), ll(miller[2]);
        int i = ipair*sym_mat.size()+isym;
        af::double6 GfunTensor = GfunTensorArray[i];
        FloatType rsSqr =
          GfunTensor[0]*hh*hh + GfunTensor[1]*kk*kk + GfunTensor[2]*ll*ll +
          GfunTensor[3]*hh*kk + GfunTensor[4]*hh*ll + GfunTensor[5]*kk*ll;
        if(rsSqr<0) {
          if(std::abs(rsSqr)<1.e-9) rsSqr=0;
          SCITBX_ASSERT(rsSqr>=0);
        }
        FloatType Gf = scitbx::math::g_function::GfuncOfRSsqr_approx(rsSqr);
        FloatType h_dot_T = scitbx::constants::two_pi*(miller*ncsDeltaT[i]);
        FloatType cosTerm = std::cos(h_dot_T);
        FloatType GcosTerm = 2.*Gf*cosTerm;
        refl_epsfac += GcosTerm*wtfac*pairs[ipair].rho_mn[sbin]; // stored parameter
        if(do_rho_mn) dEps_by_drho[ipair] += GcosTerm*wtfac;
        if(do_radius) {
          double sklmn = std::sqrt(rsSqr)/pairs[ipair].radius;
          dEps_by_dradius[ipair] += 2.*cosTerm*wtfac*pairs[ipair].rho_mn[sbin] *
            scitbx::math::g_function::dGfunc_by_dR(pairs[ipair].radius, sklmn);
        }
      }
    }
  }

  double target_gradient()
  {
    target_called=true;
    int nbins = pairs[0].rho_mn.size();
    int npars_ref = n_pairs*nbins;

    // Tabulate lists of equivalent pairs for constrained gradients
    int maxid(0);
    for (int ipair = 0; ipair < n_pairs; ipair++)
      maxid = std::max(maxid,pairs[ipair].id);
    std::vector<af::shared<FloatType> > pairs_with_id(maxid+1);
    if (do_rho_mn || do_radius)
      for (int ipair = 0; ipair < n_pairs; ipair++)
        pairs_with_id[pairs[ipair].id].push_back(ipair);

    // Reset gradient
    if(do_rho_mn) {
      Gradient_rhoMN.resize(npars_ref);
      Gradient_rhoMN.fill(0);
    }
    if(do_radius) {
      Gradient_radius.resize(n_pairs);
      Gradient_radius.fill(0);
    }
    //function, and analytic gradient for rhoMN parameters
    double minusLL(0);
    tncs_epsfac.resize(f_obs.size());
    for (unsigned r = 0; r < f_obs.size(); r++) {
      FloatType dLL_by_dEps(0);
      /*
      Look-up array of integer numbers for each reflection telling which
      resolution bin it belongs to (eg.: this Fobs belongs to bin number 12).
      */
      int s = rbin(r);
      //recalculate tncs_epsn from changed parameters
      calcRefineTerms(miller_indices[r], s);
      tncs_epsfac[r] = refl_epsfac; // result - tNCS epsilon factor
      /* Compute effective gfun with half-triple summation rather than full
         quadruple summation:
           assume NCS chosen closest to pure translation,
             ignore other NCS+symm combinations
           pair off-diagonal contributions from jncs>incs with jncs<incs
      */
      // epsn = eps = epsilon - symmetry factors, integers.
      double epsnSigmaN = epsilon[r]*tncs_epsfac[r]*SigmaN[r];
      // Code based on Wilson distribution with inflated variance for
      // measurement errors ==>
      double cent_fac = centric_flags[r] ? 0.5 : 1.0;
      double SigmaFactor = 2.*cent_fac;
      double ExpSig(SigmaFactor*scitbx::fn::pow2(sig_f_obs[r]));
      double V = epsnSigmaN + ExpSig;
      double f_obs_sq = scitbx::fn::pow2(f_obs[r]);
      // For simplicity, leave constants out of log-likelihood, which will not
      // affect refinement
      MMTBX_ASSERT(V>0);
      minusLL += cent_fac*(std::log(V) + f_obs_sq/V);
      dLL_by_dEps = cent_fac*SigmaN[r]*(V-f_obs_sq)/scitbx::fn::pow2(V);
      // <==
      for (int ipair = 0; ipair < n_pairs; ipair++) {
        int pair_id(pairs[ipair].id);
        if(do_rho_mn) {
          double gterm = dLL_by_dEps*dEps_by_drho[ipair];
          for (int jpair = 0; jpair < pairs_with_id[pair_id].size(); jpair++) {
            int this_pair = pairs_with_id[pair_id][jpair];
            Gradient_rhoMN[this_pair*nbins+s] += gterm;
          }
        }
        if(do_radius) {
          double gterm = dLL_by_dEps*dEps_by_dradius[ipair];
          for (int jpair = 0; jpair < pairs_with_id[pair_id].size(); jpair++) {
            int this_pair = pairs_with_id[pair_id][jpair];
            Gradient_radius[this_pair] += gterm;
          }
        }
      }
    }
    // sigma of 0.05 for rhoMN bin smoothness
    const double rhoMNbinwt(1./(2*scitbx::fn::pow2(0.05)));
    if(do_rho_mn || do_radius) {
      for (int ipair = 0; ipair < n_pairs; ipair++) {
        // Weakly restrain radii to initial estimates, with sigma of estimate/4
        int pair_id(pairs[ipair].id);
        if(do_radius) {
          double radiuswt(1./(2*scitbx::fn::pow2(pairs[ipair].radius_estimate/4.)));
          double delta_radius = pairs[ipair].radius - pairs[ipair].radius_estimate;
          minusLL += radiuswt*scitbx::fn::pow2(delta_radius);
          double gterm(2.*radiuswt*delta_radius);
          for (int jpair = 0; jpair < pairs_with_id[pair_id].size(); jpair++) {
            int this_pair = pairs_with_id[pair_id][jpair];
            Gradient_radius[this_pair] += gterm;
          }
        }
        if(do_rho_mn) {
          for(int s = 1; s < nbins-1; s++) { // restraints over inner bins
            double dmean = (pairs[ipair].rho_mn[s-1] + pairs[ipair].rho_mn[s+1])/2.;
            double delta_rhoMN = pairs[ipair].rho_mn[s] - dmean;
            minusLL += rhoMNbinwt*scitbx::fn::pow2(delta_rhoMN);
            double gterm(rhoMNbinwt*delta_rhoMN);
            for (int jpair = 0; jpair < pairs_with_id[pair_id].size(); jpair++) {
              int this_pair = pairs_with_id[pair_id][jpair];
              Gradient_rhoMN[this_pair*nbins+s-1] -= gterm;
              Gradient_rhoMN[this_pair*nbins+s]   += 2.*gterm;
              Gradient_rhoMN[this_pair*nbins+s+1] -= gterm;
            }
          }
        }
      }
    }
    return minusLL;
  }

  af::shared<FloatType> tncs_epsfac_result() {
    return tncs_epsfac;
  }

  af::shared<FloatType> gradient_rhoMN() {
    MMTBX_ASSERT(do_rho_mn && !do_radius && target_called);
    return Gradient_rhoMN;
  }

  af::shared<FloatType> gradient_radius() {
    MMTBX_ASSERT(!do_rho_mn && do_radius && target_called);
    return Gradient_radius;
  }

};

}} // namespace mmtbx::ncs

#endif // MMTBX_TNCS_H
