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

#include <cctbx/maptbx/real_space_gradients_simple.h> // import dependency


using namespace std;
namespace mmtbx { namespace ncs {
namespace af=scitbx::af;
namespace bp=boost::python;
using scitbx::vec3;
using scitbx::mat3;

template <typename FloatType=double>
class pair
{
  public:
    scitbx::mat3<FloatType> r;
    scitbx::vec3<FloatType> t;
    FloatType radius;
    FloatType weight;
    FloatType fracscat;
    af::shared<FloatType> rho_mn;
    pair() {
      r.fill(0);
      t.fill(0);
      radius=0;
      weight=0;
      fracscat=0;
      rho_mn.fill(0);
    }
    pair(scitbx::mat3<FloatType> const& r_,
         scitbx::vec3<FloatType> const& t_,
         FloatType radius_,
         FloatType weight_,
         FloatType fracscat_,
         af::shared<FloatType> rho_mn_)
    :
      r(r_),t(t_),radius(radius_),weight(weight_), fracscat(fracscat_),
      rho_mn(rho_mn_)
    {}
};

// --------------------------------------------------------------------------

template <typename FloatType>
class Gfunction {
public:
  std::vector<af::double6> GfunTensorArray;
  std::vector<scitbx::vec3<FloatType> > ncsDeltaT;
  std::vector<FloatType> dEps_by_drho;
  FloatType refl_epsfac;
  af::shared<pair<FloatType> > pairs;
  int n_sym_p;
  int dim;
  int n_pairs;
  scitbx::mat3<FloatType> fractionalization_matrix;
  af::shared<mat3<FloatType> > sym_mat;
  af::shared<FloatType> tncs_epsfac;

  Gfunction(
    bp::list const& pairs_,
    scitbx::mat3<FloatType> fractionalization_matrix_,
    af::shared<mat3<double> > sym_mat_,
    int n_sym_p_)
  :
  refl_epsfac(0), sym_mat(sym_mat_),
  fractionalization_matrix(fractionalization_matrix_)
  {
    for(std::size_t i=0;i<bp::len(pairs_);i++) {
      pairs.push_back(bp::extract<pair<FloatType> >(pairs_[i])());
    }
    n_pairs = pairs.size();
    n_sym_p = sym_mat.size();
    dim = n_pairs * n_sym_p;
    GfunTensorArray.resize(dim,af::double6(0,0,0,0,0,0));
    ncsDeltaT.resize(dim,scitbx::vec3<FloatType>(0,0,0));
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
      FloatType radSqr = std::pow(pairs[ipair].radius, 2);
      scitbx::mat3<FloatType> ncsR_t = pairs[ipair].r.transpose();
      scitbx::vec3<FloatType> ncsT = pairs[ipair].t;
      for (int isym = 0; isym < n_sym_p; isym++) {
        /*
          Metric matrix (represented as tensor) allows quick computation of
          |s_klmm|^2 and therefore rs^2 for G-function
          If F=fractionalization_matrix, R=Rotsym and (T) indicates transpose, then
          s_klmm = (identity-ncsR(T))*F(T)*R(T)*h
          |s_klmm|^2 = s_klmm(T).s_klmm = h(T)*R*F*(identity-ncsR)*(identity-ncsR(T))*F(T)*R(T)*h
          The metric matrix is the part between h(T) and h.  In the
          GfunTensorArray, factors of 2 are included so
          the unique terms need to be computed only once.
        */
        int i = ipair*n_sym_p+isym;
        ncsDeltaT[i] = sym_mat(isym).transpose()*ncsT; // symmetry translation cancels
        scitbx::mat3<FloatType> metric_matrix =
          sym_mat(isym).transpose()*fractionalization_matrix*(identity-ncsR_t);
        metric_matrix = radSqr*metric_matrix*metric_matrix.transpose();
        GfunTensorArray[i] = af::double6(
            metric_matrix(0,0),  metric_matrix(1,1),  metric_matrix(2,2),
          2*metric_matrix(0,1),2*metric_matrix(0,2),2*metric_matrix(1,2));
      }
    }
  }

  void calcRefineTerms(bool& do_gradient, bool& do_hessian,
                       cctbx::miller::index<int>& miller, int sbin)
  {
    /* Calculate tNCS-related epsilon factor (to multiply by normal
       symmetry-related epsilon factor), plus derivatives that can be used to
       refine parameters characterising the tNCS
    */
    dEps_by_drho.resize(n_pairs);
    refl_epsfac = 1.; // Epsilon before modulation by tNCS
    /*Use half-triple summation rather than full quadruple summation, i.e. only
      paying attention to terms from molecules that are in similar orientations.
    */
    for (int ipair = 0; ipair < n_pairs; ipair++) {
      double wtfac(pairs[ipair].fracscat);
      dEps_by_drho[ipair] = 0.;
      for (int isym = 0; isym < n_sym_p; isym++) {
        double Gf,GcosTerm;
        // NCS-related contribution to epsilon, weighted by interference,
        // G-function and rhoMN terms
        int hh(miller[0]), kk(miller[1]), ll(miller[2]);
        int i = ipair*n_sym_p+isym;
        af::double6 GfunTensor = GfunTensorArray[i];
        double rsSqr =
          GfunTensor[0]*hh*hh + GfunTensor[1]*kk*kk + GfunTensor[2]*ll*ll +
          GfunTensor[3]*hh*kk + GfunTensor[4]*hh*ll + GfunTensor[5]*kk*ll;
        Gf = scitbx::math::g_function::GfuncOfRSsqr_approx(rsSqr);
        double h_dot_T = scitbx::constants::two_pi*(miller*ncsDeltaT[i]);
        double cosTerm = std::cos(h_dot_T);
        GcosTerm = 2.*Gf*cosTerm;
        refl_epsfac += GcosTerm*wtfac*pairs[ipair].rhoMN[sbin]; // stored parameter
        if (do_gradient || do_hessian)
          dEps_by_drho[ipair] += GcosTerm*wtfac;
      }
    }
  }

  double target_gradient_hessian(
           af::const_ref<double> const& f_obs,
           af::const_ref<double> const& sig_f_obs,
           af::const_ref<double> const& SigmaN,
           af::const_ref<cctbx::miller::index<int> > const& miller_indices,
           af::const_ref<int> const& rbin, // PVA: see below for definition
           af::const_ref<int> const& epsn,
           af::const_ref<bool> const& centric_flags,
           af::shared<double>& Gradient,
           bool do_target, bool do_gradient,
           bool do_hessian)
  {
    static double EPS(std::numeric_limits<double>::min());
    int nbins = pairs[0].rho_mn.size();
    int npars_ref = n_pairs*nbins;
    if(do_gradient) {
      Gradient.resize(npars_ref);
      Gradient.fill(0);
      for (int g = 0; g < Gradient.size(); g++) Gradient[g] = 0;
    }
    //function, and analytic gradient and hessian for rhoMN parameters
    double minusLL(0);
    // Pre-compute some results that don't depend on reflection number
    /* XXX PVA: done in constructor?
    SpaceGroup sg = this->getSpaceGroup();
    UnitCell uc = this->getUnitCell();
    Gfunction gfun(sg,uc,pairs.size());
              gfun.calcArrays(pairs);
    */
    tncs_epsfac.resize(f_obs.size());
    for (unsigned r = 0; r < f_obs.size(); r++) {
      double dLL_by_dEps(0),d2LL_by_dEps2(0);
      /*
        PVA: look-up array of integer numbers for each refletion telling which
             resolution bin it belongs to. Say this Fobs belongs to bin number 12.
      */
      unsigned s = rbin(r);

      //recalculate tncs_epsn from changed parameters
      // XXX PVA: this requires update method for rho_mn done elsewhere
      calcRefineTerms(do_gradient,do_hessian,miller_indices[r],s);
      tncs_epsfac[r] = refl_epsfac; // result - tNCS epsilon factor

      if(do_target || do_gradient || do_hessian) {
        /* Compute effective gfun with half-triple summation rather than full
           quadruple summation: assume NCS chosen closest to pure translation,
           ignore other NCS+symm combinations pair off-diagonal contributions
           from jncs>incs with jncs<incs
        */
        // XXX PVA: epsn = eps = epsilon - symmetry factors, integers.
        double epsnSigmaN = epsn(r)*tncs_epsfac[r]*SigmaN[r];
        bool do_grad = (do_gradient || do_hessian); // Need pIobs gradient for both
        double grad,hess;
        // Code based on Wilson distribution with inflated variance for measurement errors ==>
        double cent_fac = centric_flags(r) ? 0.5 : 1.0; // XXX PVA: why not True/False?
        double SigmaFactor = 2.*cent_fac;
        double ExpSig(SigmaFactor*scitbx::fn::pow2(sig_f_obs[r]));
        double V = epsnSigmaN + ExpSig;
        double f_obs_sq = scitbx::fn::pow2(f_obs[r]);
        // For simplicity, leave constants out of log-likelihood, which will not affect refinement
        minusLL += cent_fac*(std::log(V) + f_obs_sq/V);
        if (do_gradient || do_hessian)
        {
          dLL_by_dEps = cent_fac*SigmaN[r]*(V-f_obs_sq)/scitbx::fn::pow2(V);
        }
        // <==
      }
      if(do_gradient) {
        for(int ipair = 0; ipair < n_pairs; ipair++) {
          Gradient[ipair*nbins+s] += dLL_by_dEps*dEps_by_drho[ipair];
        }
      }
    } // loop over reflections
    if (do_target || do_gradient || do_hessian)
    {
    const double rhoMNbinwt(1./(2*scitbx::fn::pow2(0.05))); // sigma of 0.05 for rhoMN bin smoothness
    if (do_gradient || do_hessian)
    for (int ipair = 0; ipair < n_pairs; ipair++)
    for (int s = 1; s < nbins-1; s++) // restraints over inner bins
    {
      double dmean = (pairs[ipair].rhoMN[s-1] + pairs[ipair].rhoMN[s+1])/2.;
      double delta = pairs[ipair].rhoMN[s] - dmean;
      minusLL += rhoMNbinwt*scitbx::fn::pow2(delta);
      if (do_gradient)
      {
        Gradient[ipair*nbins+s-1] -= rhoMNbinwt*delta;
        Gradient[ipair*nbins+s]   += 2.*rhoMNbinwt*delta;
        Gradient[ipair*nbins+s+1] -= rhoMNbinwt*delta;
      }
    }
    }
    return minusLL;
  }


};
// --------------------------------------------------------------------------

}} // namespace mmtbx::ncs

#endif // MMTBX_TNCS_H
