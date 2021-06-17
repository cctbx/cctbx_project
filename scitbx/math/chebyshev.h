#ifndef SCITBX_MATH_CHEBYSHEV_H
#define SCITBX_MATH_CHEBYSHEV_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>

namespace scitbx { namespace math {
namespace chebyshev{

  //! Chebyshev Base class
  template <typename FloatType = double>
  class chebyshev_base
  {
    public:
    /*! Default constructor */
    chebyshev_base() {}
    /*! Constructor that inits all members to nought */
    chebyshev_base(std::size_t const& n_terms,
                   FloatType const& low_limit,
                   FloatType const& high_limit)
    :
    n_terms_(n_terms),
    high_limit_(high_limit),
    low_limit_(low_limit),
    cheb_coefs_(n_terms,0)
    {
      SCITBX_ASSERT (n_terms>=2);
    }
    /*! Same as before, but load coefficients as well */
    chebyshev_base(std::size_t const& n_terms,
                   FloatType const& low_limit,
                   FloatType const& high_limit,
                   scitbx::af::const_ref<FloatType> const& cheb_coefs)
    :
    n_terms_(n_terms),
    high_limit_(high_limit),
    low_limit_(low_limit),
    cheb_coefs_(n_terms,0)
    {
      SCITBX_ASSERT (n_terms>=2);
      std::size_t terms=n_terms;

      if (n_terms>=cheb_coefs.size()){ terms = cheb_coefs.size();}

      for (std::size_t ii=0;ii<terms;ii++){
        cheb_coefs_[ii]=cheb_coefs[ii];
      }
      for (std::size_t ii=terms;ii<n_terms;ii++){
        cheb_coefs_[ii]=0;
      }
    }

    /*! Return chebyshev approximation for the given coefficients */
    FloatType
    f(FloatType const& x_in)
    {
      return( cheb_base_f(x_in) );
    }

    /*! Array version of previos function
     */
    scitbx::af::shared<FloatType>
    f(scitbx::af::const_ref<FloatType> const& x_in)
    {
      return( cheb_base_f(x_in) );
    }

    /*! Return the coefficients
     */
    scitbx::af::shared<FloatType>
    coefs()
    {
      return( cheb_coefs_ );
    }

    /*! Replace the coefficients with a new set.
     *  The number of coefficients can be made smaller,
     *  but not expanded
     */
    void
    replace(scitbx::af::const_ref<FloatType> const& new_coefs);

    protected:
    std::size_t n_terms_;
    FloatType high_limit_;
    FloatType low_limit_;
    scitbx::af::shared<FloatType> cheb_coefs_;


    FloatType
    transform(FloatType const& x_in);

    FloatType
    cheb_base_f(FloatType const& x_in);

    scitbx::af::shared<FloatType>
    cheb_base_f(scitbx::af::const_ref<FloatType> const& x_in);
  };

  //------------------

  template <typename FloatType>
  void
  chebyshev_base<FloatType>
  ::replace( scitbx::af::const_ref<FloatType> const& new_coefs)
  {
    std::size_t terms = new_coefs.size();
    if (terms>=n_terms_){
      terms = n_terms_;
    }
    for (unsigned ii=0;ii<terms;ii++){
      cheb_coefs_[ii]=new_coefs[ii];
    }
    for (unsigned ii=terms;ii<n_terms_;ii++){
      cheb_coefs_[ii]=0.0;
    }
  }

  //------------------

  template <typename FloatType>
  FloatType
  chebyshev_base<FloatType>
  ::transform(FloatType const& x_in)
  {
    typedef FloatType f_t;
    f_t epsilon;
    epsilon = 1.0E-12;
    f_t result=0;
    if(high_limit_ - low_limit_ != 0) {
      result = (x_in - (low_limit_ + high_limit_)*0.5)
        / (0.5*(high_limit_ - low_limit_));
    }
    SCITBX_ASSERT (result<=1+epsilon);
    SCITBX_ASSERT (result>=-1.0-epsilon);
    return(result);
  }

  //------------------

  template <typename FloatType>
  FloatType
  chebyshev_base<FloatType>
  ::cheb_base_f(FloatType const& x_in)
  {
    typedef FloatType f_t;
    f_t x;
    x = transform(x_in);
    f_t x2 = 2*x;
    f_t result = 0.0;
    f_t d=0, dd=0, sv=0;

    for (int ii=cheb_coefs_.size()-1;ii>=1;ii--){
      sv = d;
      d = x2*d-dd+cheb_coefs_[ii];
      dd = sv;
    }
    result = x*d -dd+0.5*cheb_coefs_[0];
    return (result);
  }

  //-------------------

  template <typename FloatType>
  scitbx::af::shared<FloatType>
  chebyshev_base<FloatType>
  ::cheb_base_f(scitbx::af::const_ref<FloatType> const& x_in)
  {
    scitbx::af::shared<FloatType> result(x_in.size(),0);
    for (unsigned ii=0;ii<x_in.size();ii++){
      result[ii]=cheb_base_f(x_in[ii]);
    }

    return (result);
  }

  //-----------------------------------------------------------------
  //! Chebyshev polynome, with derivative with respect to x

  /*! For this class, it only makes sense to
   *  have a constructor where coefficients are loaded.
   *  The prime use of this class is to have an approximation
   *  to a function and its derivative
   *
   *  If desired, this class could be expanded
   *  to include integrals as well.
   *
   *  Note that the 'replace' option is not available from python
   *  (see the corresponding _ext file)
   *  as this would mean that one would need to recompute the derivatives
   *  in principle, it could be done, but maybe at a later stage,
   *  when/if the need arises.
   */
  template <typename FloatType = double>
  class chebyshev_polynome: public chebyshev_base<FloatType>
  {
    public:
    chebyshev_polynome(){}

    chebyshev_polynome(std::size_t const& n_terms,
                       FloatType const& low_limit,
                       FloatType const& high_limit,
                       scitbx::af::const_ref<FloatType> const& cheb_coefs)
    :
    chebyshev_base<FloatType>(n_terms,
                              low_limit,
                              high_limit,
                              cheb_coefs),
    cheb_coefs_derivative_(n_terms,0),
    derivative_(n_terms,
                low_limit,
                high_limit)

    {
      // put the coefficients in to place
      std::size_t terms=n_terms;
      if (n_terms>=cheb_coefs.size()){ terms = cheb_coefs.size();}

      for (std::size_t ii=0;ii<terms;ii++){
        this->cheb_coefs_[ii]=cheb_coefs[ii];
      }
      for (std::size_t ii=terms;ii<n_terms;ii++){
        this->cheb_coefs_[ii]=0.0;
      }
      // compute the coefficients for the derivative
      compute_derivative_coefs();
    }

    // Give the derivative of a function at x
    FloatType
    dfdx(FloatType const& x)
    {
      return derivative_.f(x);
    }

    scitbx::af::shared<FloatType>
    dfdx(scitbx::af::const_ref<FloatType> const& x)
    {
      scitbx::af::shared<FloatType> result(x.size(),0);
      for (unsigned ii=0;ii<x.size();ii++){
        result[ii] = derivative_.f(x[ii]);
      }
      return (result);
    }


    //! Return the coefficients of the derivative
    scitbx::af::shared<FloatType>
    dfdx_coefs()
    {
      return (cheb_coefs_derivative_);
    }


    protected:
    /*! coefficient for the derivative */
    scitbx::af::shared<FloatType> cheb_coefs_derivative_;
    /*! a chebyshev polynome base class that computes
     *  the function value of the derivatives
     */
    chebyshev_base<FloatType> derivative_;


    /*! A function that computes the coefficients for the derivative
     */
    void
    compute_derivative_coefs();

  };

  //! Computes the coefficients for the derivative of a chebyshev polynome
  template <typename FloatType>
  void
  chebyshev_polynome<FloatType>
  ::compute_derivative_coefs()
  {

    cheb_coefs_derivative_[this->n_terms_-1]=0.0;
    cheb_coefs_derivative_[this->n_terms_-2]=2.0*(this->n_terms_-1)*
      this->cheb_coefs_[this->n_terms_-1];
    for (int ii= this->n_terms_ -3;ii>=0;ii--){
      cheb_coefs_derivative_[ii] = cheb_coefs_derivative_[ii+2] +
        2.0*(ii+1)*
        this->cheb_coefs_[ii+1];
    }
    FloatType domain_scale_part;
    domain_scale_part = 2.0/(this->high_limit_-this->low_limit_);
    for (int ii=0;ii< this->n_terms_; ii++){
      cheb_coefs_derivative_[ii]*=domain_scale_part;
    }
    derivative_.replace(cheb_coefs_derivative_.const_ref());

  }


  //------------------------------------------------------------
  //! A chebyshev polynome for coefficient fitting purposes
  /*! Same functionality as the base class
   *  but also returns derivatives of the function with respect to
   *  the coefficients
   */

  template <typename FloatType = double>
  class chebyshev_fitter: public chebyshev_base<FloatType>
  {
    public:
    //! Default constructor
    chebyshev_fitter(){}

    //! constructor with coefficient set to zero
    chebyshev_fitter(std::size_t const& n_terms,
                       FloatType const& low_limit,
                       FloatType const& high_limit)
    :
    chebyshev_base<FloatType>(n_terms,
                              low_limit,
                              high_limit)
    {}

    //! constructor with parameter array loading
    chebyshev_fitter(std::size_t const& n_terms,
                       FloatType const& low_limit,
                       FloatType const& high_limit,
                       scitbx::af::const_ref<FloatType> const& cheb_coefs)
    :
    chebyshev_base<FloatType>(n_terms,
                              low_limit,
                              high_limit,
                              cheb_coefs)
    {}

    scitbx::af::shared<FloatType>
    dfdcoefs(FloatType const& x_in)
    {
      scitbx::af::shared<FloatType> derivatives(this->cheb_coefs_.size(),0);
      typedef FloatType f_t;
      f_t x = this->transform(x_in);
      derivatives[0]=1.0;
      derivatives[1]=x;
      f_t x2 = 2.0*x;
      for (unsigned ii=2;ii<this->cheb_coefs_.size();ii++){
        derivatives[ii]=derivatives[ii-1]*x2 - derivatives[ii-2];
      }
      derivatives[0]=0.5;
      return(derivatives);
    }


  };



  //-------------------------------------------------------------

  //! A chebyshev polynome with gibbs damping
  /*! The Gibbs damping in this chebyshev polynome takes the
   *  form of g = 0.5*(1-tanh( (z-1)/(z*(1-z))))
   *  where z = ii/(n_terms_+1),
   *  and
   *  f(x) \approx 0.5*g_0c_0 + sum_{ii=0}^{M} g_ii c_ii T_ii(x)
   *  It ensures that the effective last coefficients are small
   *  This will large oscilations of a small period
   *  in the final function in the hope that it will prevent overfitting
   *  of data.
   *
   */

  template <typename FloatType = double>
  class chebyshev_smooth: public chebyshev_base<FloatType>
  {
    public:
    chebyshev_smooth(){}

    chebyshev_smooth(std::size_t const& n_terms,
                     FloatType const& low_limit,
                     FloatType const& high_limit)
    :
    chebyshev_base<FloatType>(n_terms,
                              low_limit,
                              high_limit),
    cheb_coefs_basic_(n_terms,0),
    cheb_coefs_smooth_(n_terms,0),
    dampening_(n_terms,0)
    {
      FloatType z;
      for (unsigned ii=0;ii<n_terms;ii++){
        // compute the damping factors
        z = ii/(n_terms+1.0);
        dampening_[ii] = 0.5*(1.0-std::tanh( (z-0.5)/(z*(1-z)) ) );
        // compute the 'new' chebyshev coefficients
      }
    }


    chebyshev_smooth(std::size_t const& n_terms,
                     FloatType const& low_limit,
                     FloatType const& high_limit,
                     scitbx::af::const_ref<FloatType> const& cheb_coefs)
    :
    chebyshev_base<FloatType>(n_terms,
                              low_limit,
                              high_limit),
    cheb_coefs_basic_(n_terms,0),
    cheb_coefs_smooth_(n_terms,0),
    dampening_(n_terms,0)
    {
      FloatType z;
      for (unsigned ii=0;ii<n_terms;ii++){
        // Basic coefficients
        cheb_coefs_basic_[ii]=cheb_coefs[ii];
        // compute the damping factors
        z = ii/(n_terms+1.0);
        if ( (z!=0) || (z!=1) ){
          dampening_[ii] = 0.5*(1.0-std::tanh( (z-0.5)/(z*(1-z)) ) );
        }
        if (z==0){
          dampening_[ii]=1.0;
        }
        if (z==1){
           dampening_[ii]=0.0;
        }
        // compute the 'new' chebyshev coefficients
        cheb_coefs_smooth_[ii] = cheb_coefs_basic_[ii]*dampening_[ii];
      }
      this->replace(cheb_coefs_smooth_.const_ref());
    }

    void
    replace_and_smooth(scitbx::af::const_ref<FloatType> const& coefs)
    {
      for (unsigned ii=0;ii<this->n_terms_;ii++){
        cheb_coefs_basic_[ii]=coefs[ii];
        cheb_coefs_smooth_[ii]=coefs[ii]*dampening_[ii];
        this->cheb_coefs_[ii] = cheb_coefs_smooth_[ii];
      }

    }

    //! In python this is accesible as smooth_coefs
    scitbx::af::shared<FloatType>
    smooth_coefs()
    {
      return (cheb_coefs_smooth_);
    }



    protected:
    scitbx::af::shared<FloatType> cheb_coefs_basic_;
    scitbx::af::shared<FloatType> cheb_coefs_smooth_;
    scitbx::af::shared<FloatType> dampening_;

  };





  //! The derivative of the chebyshev polynomial wrt to the coeficients.
  /*! Chebyshev polynome with Gibbs damping*/
  template <typename FloatType = double>
  class chebyshev_smooth_fitter: public chebyshev_base<FloatType>
  {
    public:
    //! Default constructor
    chebyshev_smooth_fitter(){}

    //! constructor with coefficient set to zero
    chebyshev_smooth_fitter(std::size_t const& n_terms,
                            FloatType const& low_limit,
                            FloatType const& high_limit)
    :
    chebyshev_base<FloatType>(n_terms,
                              low_limit,
                              high_limit),
    cheb_coefs_basic_(n_terms,0),
    cheb_coefs_smooth_(n_terms,0),
    dampening_(n_terms,0)
    {
      FloatType z;
      for (unsigned ii=0;ii<n_terms;ii++){
        // Basic coefficients
        cheb_coefs_basic_[ii]=0.0;
        // compute the damping factors
        z = ii/(n_terms+1.0);
        if ( (z!=0) || (z!=1) ){
          dampening_[ii] = 0.5*(1.0-std::tanh( (z-0.5)/(z*(1-z)) ) );
        }
        if (z==0){
          dampening_[ii]=1.0;
        }
        if (z==1){
          dampening_[ii]=0.0;
        }
        // compute the 'new' chebyshev coefficients
        cheb_coefs_smooth_[ii] = cheb_coefs_basic_[ii]*dampening_[ii];
      }
      this->replace(cheb_coefs_smooth_.const_ref());
    }
    //! constructor with parameter array loading
    chebyshev_smooth_fitter(std::size_t const& n_terms,
                       FloatType const& low_limit,
                       FloatType const& high_limit,
                       scitbx::af::const_ref<FloatType> const& cheb_coefs)
    :
    chebyshev_base<FloatType>(n_terms,
                              low_limit,
                              high_limit),
    cheb_coefs_basic_(n_terms,0),
    cheb_coefs_smooth_(n_terms,0),
    dampening_(n_terms,0)
    {
      FloatType z;
      for (unsigned ii=0;ii<n_terms;ii++){
        // Basic coefficients
        cheb_coefs_basic_[ii]=cheb_coefs[ii];
        // compute the damping factors
        z = ii/(n_terms+1.0);
        if ( (z!=0) || (z!=1) ){
          dampening_[ii] = 0.5*(1.0-std::tanh( (z-0.5)/(z*(1-z)) ) );
        }
        if (z==0){
          dampening_[ii]=1.0;
        }
        if (z==1){
          dampening_[ii]=0.0;
        }
        // compute the 'new' chebyshev coefficients
        cheb_coefs_smooth_[ii] = cheb_coefs_basic_[ii]*dampening_[ii];
      }
      this->replace(cheb_coefs_smooth_.const_ref());
    }

    scitbx::af::shared<FloatType>
    dfdcoefs(FloatType const& x_in)
    {
      scitbx::af::shared<FloatType> derivatives(this->cheb_coefs_.size(),0);
      typedef FloatType f_t;
      f_t x = this->transform(x_in);
      derivatives[0]=1.0;
      derivatives[1]=x;
      f_t x2 = 2.0*x;
      for (unsigned ii=2;ii<cheb_coefs_basic_.size();ii++){
        derivatives[ii]=derivatives[ii-1]*x2 - derivatives[ii-2];
      }
      derivatives[0]=0.5;
      for (unsigned ii=0;ii<cheb_coefs_basic_.size();ii++){
        derivatives[ii]*=dampening_[ii];
      }

      return(derivatives);
    }

    void
    replace_and_smooth(scitbx::af::const_ref<FloatType> const& new_coefs )
    {
      for (unsigned ii=0; ii< this->n_terms_; ii++){
        cheb_coefs_basic_[ii]=new_coefs[ii];
        cheb_coefs_smooth_[ii] = new_coefs[ii]*dampening_[ii];
      }
      this->replace(cheb_coefs_smooth_.const_ref());
    }
    scitbx::af::shared<FloatType>
    smooth_coefs()
    {
      return (cheb_coefs_smooth_);
    }

    protected:
    scitbx::af::shared<FloatType> cheb_coefs_basic_;
    scitbx::af::shared<FloatType> cheb_coefs_smooth_;
    scitbx::af::shared<FloatType> dampening_;

  };













  //! Some functionality for least squares fitting

  template <typename FloatType = double>
  class chebyshev_lsq
  {
    public:
    chebyshev_lsq(std::size_t const& n_terms,
                  FloatType const& low_limit,
                  FloatType const& high_limit,
                  scitbx::af::const_ref<FloatType> const& x_obs,
                  scitbx::af::const_ref<FloatType> const& y_obs,
                  scitbx::af::const_ref<FloatType> const& w_obs,
                  scitbx::af::const_ref<bool> const& free_flags)
    :
    shrd_x_obs_(x_obs.size(),0),
    shrd_y_obs_(x_obs.size(),0),
    shrd_w_obs_(x_obs.size(),0),
    shrd_free_flags_(x_obs.size(),0),
    n_terms_(n_terms),
    cheby_(n_terms, low_limit, high_limit)
    {
      SCITBX_ASSERT ( x_obs.size()==y_obs.size() );
      SCITBX_ASSERT ( x_obs.size()==w_obs.size() );
      SCITBX_ASSERT ( x_obs.size()==free_flags.size() );

      for (unsigned ii=0; ii<x_obs.size();ii++){
        shrd_x_obs_[ii]=x_obs[ii];
        shrd_y_obs_[ii]=y_obs[ii];
        shrd_w_obs_[ii]=w_obs[ii];
        shrd_free_flags_[ii] = free_flags[ii];
      }
      x_obs_ = shrd_x_obs_.const_ref();
      y_obs_ = shrd_y_obs_.const_ref();
      w_obs_ = shrd_w_obs_.const_ref();
      free_flags_ = shrd_free_flags_.const_ref();
    }

    FloatType residual();
    FloatType free_residual();

    scitbx::af::shared<FloatType> gradient();

    void replace(scitbx::af::const_ref<FloatType> const& new_coefs)
    {
      cheby_.replace(new_coefs);
    }

    scitbx::af::shared<FloatType>
    coefs()
    {
      return( cheby_.coefs() );
    }


    private:
    scitbx::af::shared<FloatType> shrd_x_obs_;
    scitbx::af::shared<FloatType> shrd_y_obs_;
    scitbx::af::shared<FloatType> shrd_w_obs_;
    scitbx::af::shared<bool> shrd_free_flags_;

    scitbx::af::const_ref<FloatType> x_obs_;
    scitbx::af::const_ref<FloatType> y_obs_;
    scitbx::af::const_ref<FloatType> w_obs_;
    scitbx::af::const_ref<bool> free_flags_;


    unsigned n_terms_;
    chebyshev_fitter<FloatType> cheby_;



  };

  template <typename FloatType>
  FloatType
  chebyshev_lsq<FloatType>
  ::residual()
  {
    FloatType result=0.0,tmp;
    for (unsigned ii=0;ii<x_obs_.size();ii++){
      if (free_flags_[ii]){
        tmp = (y_obs_[ii]- cheby_.f(x_obs_[ii]))/w_obs_[ii];
        result += tmp*tmp;
      }
    }
    return (result);
  }

  template <typename FloatType>
  FloatType
  chebyshev_lsq<FloatType>
  ::free_residual()
  {
    FloatType result=0.0,tmp;
    for (unsigned ii=0;ii<x_obs_.size();ii++){
      if (!free_flags_[ii]){
        if(w_obs_[ii]!=0) {
          tmp = (y_obs_[ii]- cheby_.f(x_obs_[ii]))/w_obs_[ii];
          result += tmp*tmp;
        }
      }
    }
    return (result);
  }

  template <typename FloatType>
  scitbx::af::shared<FloatType>
  chebyshev_lsq<FloatType>
  ::gradient()
  {

    scitbx::af::shared<FloatType> result(n_terms_,0);
    scitbx::af::shared<FloatType> tmp_grad;
    FloatType tmp=0.;
    for (unsigned ii=0; ii<x_obs_.size(); ii++){
      if (free_flags_[ii]){
        if(w_obs_[ii]!=0) {
          tmp = (y_obs_[ii] - cheby_.f(x_obs_[ii]))/(w_obs_[ii]*w_obs_[ii]);
          tmp_grad = cheby_.dfdcoefs(x_obs_[ii]);
          for (unsigned jj=0;jj<n_terms_;jj++){
            result[jj] += tmp_grad[jj]*tmp*(-2.0);
          }
        }
      }
    }
    // done
    return (result) ;

  }









  //! Some functionality for least squares fitting
  /*! This is to fit a chebyshev polynome with gibbs damping.
   *  Although elegant, it does not really seem to make
   *  a huge difference over the normal chebyshev polynome
   *  in the cases I have tested. There is !no! python interface
   *  for this, and this function might be removed at one stage or another
   */
  template <typename FloatType = double>
  class chebyshev_smooth_lsq
  {
    public:
    chebyshev_smooth_lsq(std::size_t const& n_terms,
                         FloatType const& low_limit,
                         FloatType const& high_limit,
                         scitbx::af::const_ref<FloatType> const& x_obs,
                         scitbx::af::const_ref<FloatType> const& y_obs,
                         scitbx::af::const_ref<FloatType> const& w_obs,
                         scitbx::af::const_ref<bool> const& free_flags)
    :
    n_terms_(n_terms),
    x_obs_(x_obs.size(),0),
    y_obs_(x_obs.size(),0),
    w_obs_(x_obs.size(),0),
    free_flags_(x_obs.size(),0),
    cheby_(n_terms, low_limit, high_limit)
    {
      SCITBX_ASSERT ( x_obs.size()==y_obs.size() );
      SCITBX_ASSERT ( x_obs.size()==w_obs.size() );
      SCITBX_ASSERT ( x_obs.size()==free_flags.size() );

      for (unsigned ii=0; ii<x_obs.size();ii++){
        x_obs_[ii]=x_obs[ii];
        y_obs_[ii]=y_obs[ii];
        w_obs_[ii]=w_obs[ii];
        free_flags_[ii] = free_flags[ii];
      }

    }

    FloatType residual();
    FloatType free_residual();

    scitbx::af::shared<FloatType> gradient();

    void replace(scitbx::af::const_ref<FloatType> const& new_coefs)
    {
      cheby_.replace_and_smooth(new_coefs);
    }

    scitbx::af::shared<FloatType>
    coefs()
    {
      return( cheby_.smooth_coefs() );
    }


    private:
    scitbx::af::shared<FloatType> x_obs_;
    scitbx::af::shared<FloatType> y_obs_;
    scitbx::af::shared<FloatType> w_obs_;
    scitbx::af::shared<bool> free_flags_;
    unsigned n_terms_;
    chebyshev_smooth_fitter<FloatType> cheby_;
  };

  template <typename FloatType>
  FloatType
  chebyshev_smooth_lsq<FloatType>
  ::residual()
  {
    FloatType result=0.0,tmp;
    for (unsigned ii=0;ii<x_obs_.size();ii++){
      if (free_flags_[ii]){
        tmp = (y_obs_[ii]- cheby_.f(x_obs_[ii]))/w_obs_[ii];
        result += tmp*tmp;
      }
    }
    return (result);
  }

  template <typename FloatType>
  FloatType
  chebyshev_smooth_lsq<FloatType>
  ::free_residual()
  {
    FloatType result=0.0,tmp;
    for (unsigned ii=0;ii<x_obs_.size();ii++){
      if (!free_flags_[ii]){
        tmp = (y_obs_[ii]- cheby_.f(x_obs_[ii]))/w_obs_[ii];
        result += tmp*tmp;
      }
    }
    return (result);
  }

  template <typename FloatType>
  scitbx::af::shared<FloatType>
  chebyshev_smooth_lsq<FloatType>
  ::gradient()
  {

    scitbx::af::shared<FloatType> result(n_terms_,0);
    scitbx::af::shared<FloatType> tmp_grad;
    FloatType tmp;
    for (unsigned ii=0; ii<x_obs_.size(); ii++){
      if (free_flags_[ii]){
        tmp = (y_obs_[ii] - cheby_.f(x_obs_[ii]))/(w_obs_[ii]*w_obs_[ii]);
        tmp_grad = cheby_.dfdcoefs(x_obs_[ii]);
        for (unsigned jj=0;jj<n_terms_;jj++){
          result[jj] += tmp_grad[jj]*tmp*(-2.0);
        }
      }
    }
    // done
    return (result) ;

  }

}}} // namespace scitbx::math::chebyshev

#endif // SCITBX_MATH_CHEBYSHEV_H
