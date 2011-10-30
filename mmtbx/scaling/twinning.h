//! Peter Zwart April 21, 2005
#ifndef MMTBX_SCALING_TWINNING_H
#define MMTBX_SCALING_TWINNING_H

#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/match_indices.h>
#include <cctbx/miller/asu.h>
#include <cctbx/miller/lookup_utils.h>
#include <cctbx/sgtbx/reciprocal_space_asu.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <scitbx/array_family/sort.h>
#include <scitbx/random.h>
#include <scitbx/math/erf.h>
#include <scitbx/math/bessel.h>
#include <scitbx/math/quadrature.h>
#include <cstdio>


namespace mmtbx { namespace scaling {
namespace twinning {

  template<typename FloatType>
  inline FloatType fabs_value(FloatType const& x)
  {
    FloatType eps=1e-13;
    if ( x>eps ){
      return(x);
    }
    if(x<-eps){
      return(-x);
    }
    return(+0);
  }

  template<typename FloatType>
  inline FloatType max_value(FloatType const& x, FloatType const& y )
  {
    if (x>y){
      return(x);
    }
    if (y>x){
      return(y);
    }
  }

  template<typename FloatType>
  inline FloatType substitute_when_zero(FloatType const& x, FloatType const& y )
  {
    FloatType eps=1e-13;
    if ( x<0 ){
      if (x>-eps){
        return(-y);
      }
    }
    if ( x>0 ){
      if (x<eps){
        return(y);
      }
    }
  }




  // These two lookup tables might be moved to the scitbx actually ...
  template<typename FloatType>
  class very_quick_erf
  {  // This lookup table might be less accurate, but is 40 times faster than the one in scitbx::math
  public:
    very_quick_erf(FloatType const& step_size)
    {
      SCITBX_ASSERT(step_size > 0);
      high_lim_=5.0; // hard wired high limit
      one_over_x_step_ = 1.0/step_size;
      unsigned n = static_cast<unsigned>(high_lim_*one_over_x_step_+0.5)+1;
      erf_table_.reserve(n);
      for(unsigned i=0;i<n;i++) {
        erf_table_.push_back( scitbx::math::erf(step_size * i) );
      }
    }

    FloatType erf( FloatType const& x )
    {
      FloatType ax,sign;
      if (x<0){
        sign=-1.0;
        ax = -x;
      } else{
        sign=1.0;
        ax = x;
      }
      if (ax < high_lim_){
        unsigned x_bin = static_cast<unsigned>(ax*one_over_x_step_+0.5);
        return( sign*erf_table_[x_bin] );
      } else {
        return( sign );
      }
    }

    FloatType
    loop_for_timings(std::size_t number_of_iterations, bool optimized)
    {
      FloatType result = 0;
      FloatType denom = static_cast<FloatType>(number_of_iterations / 10);
      if (optimized) {
        for (std::size_t i=0;i<number_of_iterations;i++) {
          result += erf(static_cast<FloatType>(i) / denom);
          result -= erf(static_cast<FloatType>(i) / denom);
        }
      }
      else {
        for (std::size_t i=0;i<number_of_iterations;i++) {
          result += scitbx::math::erf(static_cast<FloatType>(i) / denom);
          result -= scitbx::math::erf(static_cast<FloatType>(i) / denom);
        }
      }
      return result;
    }

  protected:
    scitbx::af::shared<FloatType> erf_table_;
    FloatType one_over_x_step_;
    FloatType high_lim_;

  };


  template<typename FloatType>
  class quick_log_ei0{
  public:
    quick_log_ei0(int const& n_points)
    {
      SCITBX_ASSERT( n_points> 50 );   // we need at least 50 points i think, although 5000 is more realistic.
      SCITBX_ASSERT( n_points< 50000); // no problems below 50000, did not check for larger values. most likely not needed
      n_ = n_points;                   // a factor 5 in timings is gained over the full computation
      FloatType t;
      t_step_ = 1.0/static_cast<FloatType>(n_);
      t_table_.reserve(n_);
      log_ei0_table_.reserve(n_);
      for (int ii=0;ii<n_-1;ii++){
        t = ii*t_step_;
        t_table_.push_back( t );
        t = -t/(1-t) + scitbx::math::bessel::ln_of_i0( t/(1-t) );
        log_ei0_table_.push_back( t );
      }
      t_table_.push_back(1.0);
      log_ei0_table_.push_back( log_ei0_table_[ n_-2 ]/2.0 ); // i don't want a zero here
    }

    FloatType log_ei0( FloatType const& x )
    {

      FloatType t,xx;
      xx = std::fabs(x);
      t = xx/(1.0+xx);
      if (t<0){
        t=1e-13;
      }
      int t_bin_low, t_bin_high;
      t_bin_low = int( std::floor( t*n_ ) );
      t_bin_high = t_bin_low+1;

      SCITBX_ASSERT( t >= 0);
      SCITBX_ASSERT( t_bin_low>= 0);
      FloatType f0 = log_ei0_table_[ t_bin_low ];
      FloatType f1 = log_ei0_table_[ t_bin_high ];
      FloatType t0 = t_table_[t_bin_low];
      //FloatType t1 = t_table_[t_bin_high];
      // linear interpolation
      FloatType alpha = (t-t0)*n_;///(t1-t0);
      FloatType result = f0*(1-alpha)+alpha*f1;
      return( result );
    }

    FloatType
    loop_for_timings(std::size_t number_of_iterations, bool optimized)
    {
      FloatType result = 0;
      FloatType denom = static_cast<FloatType>(number_of_iterations / 10);
      if (optimized) {
        for (std::size_t i=0;i<number_of_iterations;i++) {
          FloatType x = static_cast<FloatType>(i) / denom;
          result += log_ei0(x);
          result -= log_ei0(x);
        }
      }
      else {
        for (std::size_t i=0;i<number_of_iterations;i++) {
          FloatType x = static_cast<FloatType>(i) / denom;
          result += std::exp(-x)*scitbx::math::bessel::ln_of_i0(x);
          result -= std::exp(-x)*scitbx::math::bessel::ln_of_i0(x);
        }
      }
      return result;
    }

  protected:
    scitbx::af::shared<FloatType> t_table_;
    scitbx::af::shared<FloatType> log_ei0_table_;
    int n_;
    FloatType t_step_;
    FloatType one_over_t_step;

  };





  template<typename FloatType>
  class quick_ei0{
  public:
    quick_ei0(int const& n_points)
    {
      SCITBX_ASSERT( n_points> 50 );   // we need at least 50 points i think, although 5000 is more realistic.
      SCITBX_ASSERT( n_points< 50000); // no problems belwo 50000, did not check for larger values. most likely not needed
      n_ = n_points;                   // a factor 5 in timings is gained over the full computation
      FloatType t;
      t_step_ = 1.0/static_cast<FloatType>(n_);
      t_table_.reserve(n_);
      ei0_table_.reserve(n_);
      for (int ii=0;ii<n_-1;ii++){
        t = ii*t_step_;
        t_table_.push_back( t );
        t = -t/(1-t) + scitbx::math::bessel::ln_of_i0( t/(1-t) );
        ei0_table_.push_back( std::exp(t) );
      }
      t_table_.push_back(1.0);
      ei0_table_.push_back( ei0_table_[ n_-2 ]/2.0 ); // i don't want a zero here
    }

    FloatType ei0( FloatType const& x )
    {

      FloatType t,xx;
      xx = std::fabs(x);
      t = xx/(1.0+xx);
      if (t<0){
        t=1e-13;
      }
      int t_bin_low, t_bin_high;
      t_bin_low = int( std::floor( t*n_ ) );
      t_bin_high = t_bin_low+1;

      SCITBX_ASSERT( t >= 0);
      SCITBX_ASSERT( t_bin_low>= 0);
      FloatType f0 = ei0_table_[ t_bin_low ];
      FloatType f1 = ei0_table_[ t_bin_high ];
      FloatType t0 = t_table_[t_bin_low];
      //FloatType t1 = t_table_[t_bin_high];
      // linear interpolation
      FloatType alpha = (t-t0)*n_;///(t1-t0);
      FloatType result = f0*(1-alpha)+alpha*f1;
      return( result );
    }

    FloatType
    loop_for_timings(std::size_t number_of_iterations, bool optimized)
    {
      FloatType result = 0;
      FloatType denom = static_cast<FloatType>(number_of_iterations / 10);
      if (optimized) {
        for (std::size_t i=0;i<number_of_iterations;i++) {
          FloatType x = static_cast<FloatType>(i) / denom;
          result += ei0(x);
          result -= ei0(x);
        }
      }
      else {
        for (std::size_t i=0;i<number_of_iterations;i++) {
          FloatType x = static_cast<FloatType>(i) / denom;
          result += std::exp(-x)*scitbx::math::bessel::i0(x);
          result -= std::exp(-x)*scitbx::math::bessel::i0(x);
        }
      }
      return result;
    }

  protected:
    scitbx::af::shared<FloatType> t_table_;
    scitbx::af::shared<FloatType> ei0_table_;
    int n_;
    FloatType t_step_;
    FloatType one_over_t_step;

  };



  template<typename FloatType>
  inline cctbx::miller::index<> twin_mate( cctbx::miller::index<> const& hkl,
                                           scitbx::mat3<FloatType> const& twin_law )
  {
          int ht,kt,lt;
          ht = scitbx::math::iround(twin_law[0]*hkl[0] +
                                    twin_law[3]*hkl[1] +
                                    twin_law[6]*hkl[2]);

          kt = scitbx::math::iround(twin_law[1]*hkl[0] +
                                    twin_law[4]*hkl[1] +
                                    twin_law[7]*hkl[2]);

          lt = scitbx::math::iround(twin_law[2]*hkl[0] +
                                    twin_law[5]*hkl[1] +
                                    twin_law[8]*hkl[2]);

          cctbx::miller::index<> hkl_twin(ht,kt,lt);
          return( hkl_twin );
  }




  template<typename FloatType>
  FloatType twin_correlation( cctbx::miller::index<> const& hkl,
                              scitbx::af::const_ref<FloatType> const& data,
                              cctbx::sgtbx::space_group const& space_group,
                              bool const& anomalous_flag,
                              scitbx::mat3<FloatType> const& twin_law )
  {
    SCITBX_ASSERT( hkl.size() == data.size() );
    int twin_index;
    // make a lookup table please
    cctbx::miller::lookup_utils::lookup_tensor<FloatType> hkl_lookup(hkl,space_group,anomalous_flag);
    FloatType mean1=0,mean2=0,sig1=0,sig2=0,cov=0,x,y,n=0;
    for (int ii=0;ii<hkl.size();ii++){
      twin_index = hkl_lookup.find_hkl( twin_mate(hkl[ii],twin_law) );
      if (twin_index>=0){
        x =  data[twin_index];
        y = data[ii];
        mean1+=x;
        mean2+=y;
        cov+=x*y;
        sig1+=x*x;
        sig2+=y*y;
        n++;
      }
    }
    mean1/=n;
    mean2/=n;
    sig1/=n;
    sig2/=n;
    sig1=sig1-mean1*mean1;
    sig2=sig2-mean2*mean2;
    cov/=n;
    cov=cov-mean1*mean2;
    sig1=std::sqrt(sig1*sig2);
    if(sig1<=0){
      return(1.0);
    }
    return( cov/sig1 );

  }






  template <typename FloatType=double>
  class twin_r
  {
    public:
    /*! default constructor */
    twin_r() {}

    twin_r(scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
           scitbx::af::const_ref< FloatType > const& intensity,
           cctbx::sgtbx::space_group const& space_group,
           bool const& anomalous_flag,
           scitbx::mat3<FloatType> twin_law)
    :
    hkl_(),
    hkl_twin_(),
    location_(),
    intensity_(),
    r_sq_top_(0),
    r_sq_bottom_(0),
    r_abs_top_(0),
    r_abs_bottom_(0),
    cc_(0),
    twin_law_(twin_law),
    space_group_(space_group),
    hkl_lookup_(hkl, space_group, anomalous_flag)
    {
      SCITBX_ASSERT( hkl.size() == intensity.size() );
      int ht,kt,lt;
      for (unsigned ii=0;ii<hkl.size();ii++){
        hkl_.push_back( hkl[ii] );
        intensity_.push_back( intensity[ii] );
        // using iround to prevent mess ups
        ht = scitbx::math::iround(
               twin_law[0]*hkl[ii][0] +
               twin_law[3]*hkl[ii][1] +
               twin_law[6]*hkl[ii][2]);

        kt = scitbx::math::iround(
               twin_law[1]*hkl[ii][0] +
               twin_law[4]*hkl[ii][1] +
               twin_law[7]*hkl[ii][2]);

        lt = scitbx::math::iround(
               twin_law[2]*hkl[ii][0] +
               twin_law[5]*hkl[ii][1] +
               twin_law[8]*hkl[ii][2]);

        cctbx::miller::index<> tmp_ht(ht,kt,lt);
        // map the thing to the ASU please


        hkl_twin_.push_back( tmp_ht );
        location_.push_back( hkl_lookup_.find_hkl( tmp_ht ) );

      }

        compute_r_abs_value();
        compute_r_sq_value();
      compute_correlation();
    }


    void compute_correlation()
    {
      FloatType mean1=0,mean2=0,sig1=0,sig2=0,cov=0,x,y,n=0;

      int tmp_loc;
      for (unsigned ii=0;ii<hkl_.size();ii++){
        tmp_loc =  location_[ ii ];
        if (tmp_loc >= 0){
          if (tmp_loc != ii ){
            x = intensity_[ ii ];
            y = intensity_[ tmp_loc ];
            mean1+=x;
            mean2+=y;
            cov+=x*y;
            sig1+=x*x;
            sig2+=y*y;
            n+=1.0;
          }
        }
      }
      if (n>0){
        CCTBX_ASSERT(n != 0);
        mean1/=n;
        mean2/=n;
        sig1/=n;
        sig2/=n;
        sig1=sig1-mean1*mean1;
        sig2=sig2-mean2*mean2;
        cov/=n;
        cov=cov-mean1*mean2;
        sig1=std::sqrt(sig1*sig2);

        if (sig1<=0){
          cc_=1.0;
        } else{
          cc_=cov/sig1;
          CCTBX_ASSERT(sig1 != 0);
        }
      }
      else{
        cc_ = 0.0;
      }

    }



    void compute_r_abs_value()
    {
      FloatType top=0, bottom=0;
      FloatType x1,x2;
      int tmp_loc;
      for (unsigned ii=0;ii<hkl_.size();ii++){
        tmp_loc =  location_[ ii ];
        if (tmp_loc >= 0){
          if (tmp_loc != ii ){
            x1 = intensity_[ ii ];
            x2 = intensity_[ tmp_loc ];
            top += std::fabs(x1-x2);
            bottom += std::fabs( x1+x2 );
          }
        }
      }

      if (top>0){
        if (bottom>0){
          r_abs_bottom_ = bottom;
          r_abs_top_ = top;
        }
      }
    }

    void compute_r_sq_value()
    {
      FloatType top=0, bottom=0;
      FloatType x1,x2;
      int tmp_loc;
      for (unsigned ii=0;ii<hkl_.size();ii++){
        tmp_loc =  location_[ ii ];
        if (tmp_loc >= 0){
          if (tmp_loc != ii ){
            x1 = intensity_[ ii ];
            x2 = intensity_[ tmp_loc ];
            top += (x1-x2)*(x1-x2);
            bottom += ( x1+x2 )*( x1+x2 );
          }
        }
      }

      if (top>0){
        if (bottom>0){
          r_sq_bottom_ = bottom;
          r_sq_top_ = top;
        }
      }
    }



    FloatType r_abs_value()
    {
      FloatType result=0;
      if (r_abs_bottom_>0){
          result=r_abs_top_/r_abs_bottom_;
      }
      return (result);
    }

    FloatType r_sq_value()
    {
      FloatType result=0;
      if (r_sq_bottom_>0){
          result=r_sq_top_/r_sq_bottom_;
      }

      return (result);
    }

    FloatType correlation(){ return(cc_); }

    scitbx::af::tiny<FloatType,2> r_abs_pair()
    {
      scitbx::af::tiny<FloatType, 2> result;
      result[0] = r_abs_top_;
      result[1] = r_abs_bottom_;
      return (result);
    }

    scitbx::af::tiny<FloatType,2> r_sq_pair()
    {
      scitbx::af::tiny<FloatType, 2> result;
      result[0] = r_sq_top_;
      result[1] = r_sq_bottom_;
      return(result);
    }





    protected:

    FloatType r_abs_top_;
    FloatType r_abs_bottom_;
    FloatType r_sq_top_;
    FloatType r_sq_bottom_;
    FloatType cc_;

    scitbx::af::shared< cctbx::miller::index<> > hkl_;
    scitbx::af::shared< cctbx::miller::index<> > hkl_twin_;
    scitbx::af::shared< FloatType > intensity_;

    scitbx::af::shared< int > location_;
    scitbx::mat3<FloatType> twin_law_;
    cctbx::sgtbx::space_group space_group_;
    cctbx::miller::lookup_utils::lookup_tensor<
      FloatType> hkl_lookup_;

  };



  template <typename FloatType=double>
  class l_test
  {
    public:
    /*! Default constructor */
    l_test() {}

    l_test(scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
           scitbx::af::const_ref< FloatType > const& intensity,
           cctbx::sgtbx::space_group const& space_group,
           bool const& anomalous_flag,
           long parity_h,
           long parity_k,
           long parity_l,
           std::size_t max_delta_h)
    :
    parity_h_(parity_h),
    parity_k_(parity_k),
    parity_l_(parity_l),
    anomalous_flag_(anomalous_flag),
    max_delta_h_(max_delta_h),
    l_values_(),
    mean_l_(0),
    mean_l2_(0),
    hkl_lookup_( hkl, space_group, anomalous_flag ),
    hkl_(),
    intensity_(),
    diff_vectors_(),
    generator_(0), // always use the seed 0
    cumul_(50,0)
    {
      SCITBX_ASSERT ( hkl.size() == intensity.size() );
      SCITBX_ASSERT ( int(max_delta_h_)>=2 );


      for (unsigned ii=0;ii<hkl.size();ii++){
        intensity_.push_back(intensity[ii]);
        hkl_.push_back(hkl[ii]);
      }
      setup_diff_vectors();
      generate_pairs_and_compute_l_values();
      make_cumul();
      ml_estimate_alpha();
    }

    FloatType mean_l(){ return(mean_l_);}
    FloatType mean_l2(){ return(mean_l2_);}
    FloatType ml_alpha(){ return(ml_alpha_); }
    scitbx::af::shared<FloatType> cumul(){ return(cumul_);}


    protected:
    void setup_diff_vectors()
    {
      int tmp = max_delta_h_;

      for (int ii =-tmp; ii<=tmp; ii++){
        for (int jj =-tmp; jj<=tmp; jj++){
          for (int kk=-tmp; kk<=tmp;kk++){
            bool accept=true;

            if (ii%parity_h_!=0){  accept=false; }
            if (jj%parity_k_!=0){  accept=false; }
            if (kk%parity_l_!=0){  accept=false; }

            if ( abs(ii)+abs(jj)+abs(kk) >= tmp ){
              accept = false;
            }
            if ( abs(ii)+abs(jj)+abs(kk) < 2){ // not direct neighbours or self
              accept=false;
            }

            if (accept){
              cctbx::miller::index<> d_hkl(ii,jj,kk);
              diff_vectors_.push_back( d_hkl );
            }

          } // kk
        } // jj
      } // ii

    }

    void generate_pairs_and_compute_l_values()
    {
      scitbx::af::shared< size_t > diff_pick(hkl_.size(),0);
      diff_pick = generator_.random_size_t(
         hkl_.size(), diff_vectors_.size()  );

      long ht,kt,lt,index;
      FloatType p,q, ll,count=0;
      for (unsigned ii=0;ii<hkl_.size();ii++){
        ht = diff_vectors_[diff_pick[ii]][0] + hkl_[ii][0];
        kt = diff_vectors_[diff_pick[ii]][1] + hkl_[ii][1];
        lt = diff_vectors_[diff_pick[ii]][2] + hkl_[ii][2];

        cctbx::miller::index<> picked_pair(ht,kt,lt);
        index = hkl_lookup_.find_hkl(picked_pair);
        if (index >=0){
          p = intensity_[ii];
          q = intensity_[index];
          ll = (p-q)/(p+q);
          l_values_.push_back(ll);
          mean_l_+= std::fabs(ll);
          mean_l2_+=ll*ll;
          count++;
        }
      }

      mean_l_/=count;
      mean_l2_/=count;
    }

    void make_cumul()
    {
      FloatType l_current;
      for (unsigned ii=0;ii<50;ii++){
        l_current=ii/50.0;
        for (unsigned jj=0; jj<l_values_.size(); jj++){
          if (std::fabs(l_values_[jj])<=l_current){
            cumul_[ii]+=1.0;
          }
        }
      }
      for (unsigned ii=0;ii<50;ii++){
        cumul_[ii]/=static_cast<FloatType>(l_values_.size());
      }

    }


    FloatType log_p_alpha(FloatType l_value, FloatType alpha)
    {
      SCITBX_ASSERT( fabs(l_value) <=1);
      SCITBX_ASSERT( alpha < 0.5);
      SCITBX_ASSERT( alpha >=0 );
      FloatType P=0.;
      FloatType den = (l_value*l_value-1)*
      (-1+l_value*l_value+2.0*alpha*
       (-1+alpha+(-1+alpha)*(3+4*(-1+alpha)*alpha)*l_value*l_value));
      if(den == 0) return P;
      den = ((-1+(1-2*alpha)*(1-2*alpha)*l_value*l_value)*
      (-1+(1-2*alpha)*(1-2*alpha)*l_value*l_value));
      if(den == 0) return den;
      P/=den;
      if (P < 1e-10){
        P=1e-10;
      }
      return ( log(P) );
    }

    void ml_estimate_alpha()
    {
      scitbx::af::shared<FloatType> nll_alpha(500,0);
      FloatType trial_alpha;
      FloatType max_nll=-1E20;
      FloatType max_nll_alpha = 0.0;
      for (unsigned ii=0;ii<500;ii++){
        trial_alpha = ii/1001.0;
        for (unsigned jj=0;jj<l_values_.size();jj++){
          nll_alpha[ii]+=log_p_alpha( l_values_[jj], trial_alpha );
        }
        if (nll_alpha[ii] > max_nll){
          max_nll = nll_alpha[ii];
          max_nll_alpha =trial_alpha;
        }
      }
      ml_alpha_ = max_nll_alpha;
    }

    long parity_h_;
    long parity_k_;
    long parity_l_;
    bool anomalous_flag_;
    std::size_t max_delta_h_;
    scitbx::af::shared<FloatType> l_values_;
    FloatType mean_l_;
    FloatType mean_l2_;
    FloatType ml_alpha_;
    cctbx::miller::lookup_utils::lookup_tensor<
     FloatType> hkl_lookup_;
    scitbx::af::shared< cctbx::miller::index<> > hkl_;
    scitbx::af::shared<FloatType> intensity_;
    scitbx::af::shared< cctbx::miller::index<> > diff_vectors_;
    scitbx::random::mersenne_twister generator_;
    scitbx::af::shared<FloatType> cumul_;

  };



  template <typename FloatType=double>
  class detwin
  {
    public:
    detwin() {}
    detwin(scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
           scitbx::af::const_ref< FloatType > const& intensity,
           scitbx::af::const_ref< FloatType > const& sigma_i,
           cctbx::sgtbx::space_group const& space_group,
           bool const& anomalous_flag,
           scitbx::mat3<FloatType> twin_law)
    :
    hkl_(),
    hkl_twin_(),
    location_(),
    intensity_(),
    sigma_(),
    detwinned_hkl_(),
    detwinned_i_(),
    detwinned_sigi_(),
    twin_law_(twin_law),
    space_group_(space_group),
    hkl_lookup_(hkl, space_group, anomalous_flag),
    completeness_(0.0)
    {
      SCITBX_ASSERT( hkl.size() == intensity.size() );
      SCITBX_ASSERT( (sigma_i.size()==0) ||  (hkl.size() == sigma_i.size()) );
      int ht,kt,lt;
      for (unsigned ii=0;ii<hkl.size();ii++){
        hkl_.push_back( hkl[ii] );
        intensity_.push_back( intensity[ii] );
        if ( sigma_i.size() > 0 ){
          sigma_.push_back( sigma_i[ii] );
        }
        if ( sigma_i.size() == 0 ){
          sigma_.push_back( 0.0 );
        }
        //
        // using iround to prevent mess ups
        ht = scitbx::math::iround(
               twin_law[0]*hkl[ii][0] +
               twin_law[3]*hkl[ii][1] +
               twin_law[6]*hkl[ii][2]);

        kt = scitbx::math::iround(
               twin_law[1]*hkl[ii][0] +
               twin_law[4]*hkl[ii][1] +
               twin_law[7]*hkl[ii][2]);

        lt = scitbx::math::iround(
               twin_law[2]*hkl[ii][0] +
               twin_law[5]*hkl[ii][1] +
               twin_law[8]*hkl[ii][2]);

        cctbx::miller::index<> tmp_ht(ht,kt,lt);
        // map the thing to the ASU please


        hkl_twin_.push_back( tmp_ht );
        location_.push_back( hkl_lookup_.find_hkl( tmp_ht ) );
        if (location_[ii]<0){
          completeness_++;
        }
      }

      completeness_/=double(hkl_.size());
      completeness_=1-completeness_;
    }

    FloatType
    detwin_with_alpha(FloatType alpha)
    {
      return( detwin_with_alpha_(alpha) );
    }

    scitbx::af::shared< FloatType >
    detwinned_i(){
      SCITBX_ASSERT ( detwinned_i_.size() >0 );
      return( detwinned_i_ );
    }

    scitbx::af::shared< FloatType >
    detwinned_sigi(){
      SCITBX_ASSERT ( detwinned_sigi_.size() >0 );
      return( detwinned_sigi_);
    }

    scitbx::af::shared< cctbx::miller::index<> >
    detwinned_hkl(){
      SCITBX_ASSERT ( detwinned_hkl_.size() >0 );
      return( detwinned_hkl_);
    }


    scitbx::af::shared< long >
    location(){
      return( location_);
    }



    FloatType
    completeness(){ return(completeness_); }

    protected:
    FloatType
    detwin_with_alpha_(FloatType alpha)
    {
      SCITBX_ASSERT(alpha>=0.0);
      SCITBX_ASSERT(alpha<0.5);
      SCITBX_ASSERT(completeness_>0);
      detwinned_i_.clear();
      detwinned_sigi_.clear();
      detwinned_hkl_.clear();

      FloatType total=0;
      FloatType negatives=0;
      FloatType Iold1,Iold2,Sold1,Sold2,Inew,Snew;
      FloatType tmp_mult = std::sqrt( 1-2*alpha +2*alpha*alpha)/(1-2.0*alpha);

      for (unsigned ii=0;ii<hkl_.size();ii++){
        int tmp1 = ii;
        int tmp2 = location_[ii];
        if (tmp2 >=0){
          Iold1 = intensity_[tmp1];
          Iold2 = intensity_[tmp2];
          Sold1 = sigma_[tmp1];
          Sold2 = sigma_[tmp2];

          Inew = ((1.0-alpha)*Iold1 - alpha*Iold2)/(1-2.0*alpha);

          Snew =  tmp_mult*std::sqrt((Sold1*Sold1 + Sold1*Sold2)/2.0);


          detwinned_i_.push_back( Inew );
          detwinned_sigi_.push_back(Snew);
          detwinned_hkl_.push_back( hkl_[ii] );

          total++;
          if (Inew<0){
            negatives++;
          }
        }
      }
      return (negatives/total);
    }


    scitbx::af::shared< cctbx::miller::index<> > hkl_;
    scitbx::af::shared< cctbx::miller::index<> > hkl_twin_;
    scitbx::af::shared< long > location_;
    scitbx::af::shared< FloatType > intensity_;
    scitbx::af::shared< FloatType > sigma_;
    scitbx::af::shared< cctbx::miller::index<> > detwinned_hkl_;
    scitbx::af::shared< FloatType > detwinned_i_;
    scitbx::af::shared< FloatType > detwinned_sigi_;
    scitbx::mat3<FloatType> twin_law_;
    cctbx::sgtbx::space_group space_group_;
    cctbx::miller::lookup_utils::lookup_tensor<
      FloatType> hkl_lookup_;
    FloatType completeness_;

  };



  template <typename FloatType=double >
  class h_test
  {
    public:
    h_test(){}

    h_test(scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
           scitbx::af::const_ref< FloatType > const& intensity,
           scitbx::af::const_ref< FloatType > const& sigma_i,
           cctbx::sgtbx::space_group const& space_group,
           bool const& anomalous_flag,
           scitbx::mat3<FloatType> twin_law,
           FloatType const& fraction)
    :
    detwin_object_(hkl,
                   intensity,
                   sigma_i,
                   space_group,
                   anomalous_flag,
                   twin_law),
    location_(),
    intensity_(),
    h_values_(),
    h_array_(),
    h_cumul_obs_(),
    h_cumul_fit_(),
    fraction_(fraction),
    distance_(),
    fitted_alpha_(),
    mean_h_(),
    mean_h2_()
    {
      location_ = detwin_object_.location();
      for (unsigned ii=0;ii<hkl.size();ii++){
        intensity_.push_back( intensity[ii] );
      }
      FloatType p,q;
      scitbx::af::shared<FloatType> temp_p_values;
      scitbx::af::shared<FloatType> temp_q_values;
      scitbx::af::shared<FloatType> temp_sum_values;
      long index;
      for (unsigned ii=0;ii<intensity_.size();ii++){
        index = location_[ii];
        if (index >=0){
          p = intensity_[ii];
          q = intensity_[index];
          temp_p_values.push_back(p);
          temp_q_values.push_back(q);
          temp_sum_values.push_back(p+q);
        }
      }
      // Sort the temp_sum_values array
      scitbx::af::shared<std::size_t> sort_permutation;
      sort_permutation = scitbx::af::sort_permutation(
        temp_sum_values.const_ref(),
        true);
      unsigned limit_slot=unsigned( temp_sum_values.size()*fraction_ );
      if (limit_slot == 0) {
        throw std::invalid_argument("limit_slot == 0");
      }
      for (unsigned ii=0;ii<limit_slot;ii++){
        p = temp_p_values[ sort_permutation[ii] ];
        q = temp_q_values[ sort_permutation[ii] ];
        h_values_.push_back( std::abs(p-q) / (p+q) );
        mean_h_+= std::abs(p-q) / (p+q);
        mean_h2_+= (std::abs(p-q) / (p+q))*(std::abs(p-q) / (p+q));
      }
      mean_h_/=static_cast<FloatType>(limit_slot);
      mean_h2_/=static_cast<FloatType>(limit_slot);
      for (unsigned ii=0;ii<50;ii++){
        h_array_.push_back(ii/50.0);
        h_cumul_obs_.push_back(0.0);
        h_cumul_fit_.push_back(0.0);
      }
      make_cumul_();
      fit_cumul_();
    }

    scitbx::af::shared<FloatType>
    h_array()
    {
      return( h_array_ );
    }

    scitbx::af::shared<FloatType>
    h_values()
    {
      return( h_values_ );
    }
    scitbx::af::shared<FloatType>
    h_cumul_obs()
    {
      return( h_cumul_obs_ );
    }
    scitbx::af::shared<FloatType>
    h_cumul_fit()
    {
      return( h_cumul_fit_ );
    }
    FloatType
    distance()
    {
      return (distance_);
    }
    FloatType
    alpha()
    {
      return (fitted_alpha_);
    }
    FloatType
    mean_h()
    {
      return (mean_h_);
    }
    FloatType
    mean_h2()
    {
      return (mean_h2_);
    }

    protected:
    void
    make_cumul_()
    {
      FloatType h;
      for (unsigned ii=0;ii<h_array_.size();ii++){
        h = h_array_[ii];
        for (unsigned jj=0;jj<h_values_.size();jj++){
          if (h>=h_values_[jj]){
            h_cumul_obs_[ii]++;
          }
        }
      }
      for (unsigned ii=0;ii<h_array_.size();ii++){
        h_cumul_obs_[ii]/=static_cast<FloatType>(h_values_.size());
      }
    }

    void
    fit_cumul_()
    {
      // test 500 values of alpha please
      FloatType alpha, best_alpha=0.0;
      FloatType local_max;
      FloatType global_min=20.0;
      FloatType sh;
      for (unsigned ii=0;ii<500;ii++){
        alpha = ii/1001.0;
        local_max = 0.0;
        for (unsigned jj=0;jj<h_cumul_obs_.size();jj++){
          sh = h_array_[jj]/(1-2.0*alpha);
          if (sh > 1.0 ){
            sh = 1.0;
          }
          if ( std::abs(sh - h_cumul_obs_[jj]) > local_max ){
            local_max = std::abs(sh - h_cumul_obs_[jj]);
          }
        }
        if (local_max <= global_min){
          global_min = local_max;
          best_alpha = alpha;
        }
      }
      distance_= global_min;
      fitted_alpha_ = best_alpha;

      for (unsigned jj=0;jj<h_array_.size();jj++){
        sh = h_array_[jj]/(1-2.0*best_alpha);
        if (sh > 1.0 ){
            sh = 1.0;
        }
        h_cumul_fit_[jj]=(sh);
      }
    }

    detwin<FloatType> detwin_object_;
    scitbx::af::shared< long > location_;
    scitbx::af::shared< FloatType > intensity_;
    scitbx::af::shared< FloatType > h_values_;
    scitbx::af::shared< FloatType > h_array_;
    scitbx::af::shared< FloatType > h_cumul_obs_;
    scitbx::af::shared< FloatType > h_cumul_fit_;
    FloatType fraction_;
    FloatType distance_;
    FloatType fitted_alpha_;
    FloatType mean_h_;
    FloatType mean_h2_;
  };




  /* This is an implementation using Results from Randy Rewads Notes.
   * This provides basic routines needed to compute a likelihood for a twin fraction
   * A numerical integration is accried out. In limited tests, aboput 16 points gave resonable results
   * Note that th jhumber of points used is always twice what one specifies.
   */
  template < typename FloatType >
  class ml_murray_rust
  {
  public:
    // input normalized structure factors.
    ml_murray_rust( scitbx::af::const_ref<FloatType> const& z,
                    scitbx::af::const_ref<FloatType> const& sig_z,
                    scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
                    cctbx::sgtbx::space_group const& space_group,
                    bool const& anomalous_flag,
                    scitbx::mat3<FloatType> const& twin_law,
                    long const& n_hermite):
    qerf_(0.001), // an erf LUT
    n_hermite_(n_hermite) // gauss hermite related stuff
    {
      SCITBX_ASSERT( z.size() == sig_z.size() );
      SCITBX_ASSERT( z.size() == hkl.size() );
      cctbx::miller::index<> tmp_index;
      cctbx::miller::lookup_utils::lookup_tensor<FloatType>
        tmp_lut(hkl, space_group, anomalous_flag);

      for (std::size_t ii=0;ii<z.size();ii++){
        tmp_index = twin_mate( hkl[ii], twin_law );
        twin_mate_index_.push_back( tmp_lut.find_hkl( tmp_index ) );
        z_.push_back( z[ii] );
        sig_z_.push_back( sig_z[ii] );
      }

      // settting up the gauss hermite quadrature
      scitbx::math::quadrature::gauss_hermite_engine<FloatType> tmp_ghe(n_hermite_);
      x_ = tmp_ghe.x();
      wexs_ = tmp_ghe.w_exp_x_squared();

    }

    inline
    FloatType part_1( FloatType const& ito1, FloatType const& sito1,
                      FloatType const& ito2, FloatType const& sito2,
                      FloatType const& it2,  FloatType const& t)
    {
      FloatType tmp_a,tmp_b,tmp_c,tot=0;

      if (it2>=0){
        FloatType tmp1 = sito2*sito2;
        if(tmp1 != 0) {
          tmp_a = -(it2-ito2)*(it2-ito2)/tmp1
                  +(sito1*sito1)/(1.0*1.0)
                  -(2.0*(it2+ito1) )/(1.0);
          if(tmp_a > 300.0) tmp_a=300.0;
          tmp_a = std::exp(0.5*tmp_a);

          tmp_b =  it2*1.0*(t-1.0)
                  -sito1*sito1*t
                  +ito1*1.0*t;
          FloatType tmp2 = std::sqrt(2.0)*sito1*1.0*t;
          if(tmp1 != 0) {
            tmp_b = tmp_b / tmp2;
            tmp_b = qerf_.erf( tmp_b );
            if(t != 1.) {
              tmp_c = -ito1 + sito1*sito1/1.0 +it2*t/(1-t) ;
              if(sito1 != 0) {
                tmp_c = tmp_c/( std::sqrt(2.0) * sito1 );
                tmp_c = qerf_.erf( tmp_c );

                tot = tmp_a*(tmp_b+tmp_c);
                FloatType tmp3 = (2.0*std::sqrt(2.0*scitbx::constants::pi)*sito2*1.0*1.0*(2.0*t-1.0));
                if(tmp3 != 0) tot/=tmp3;
              }
            }
          }
        }
      }
      return(tot);
    }


    FloatType num_int( FloatType const& ito1,    FloatType const& sito1,
                       FloatType const& ito2,    FloatType const& sito2,
                       FloatType const& low_sig, FloatType const& high_sig,
                       FloatType const& t,
                       int const& n)
    {
      // carry out numerical integration using extended Simpson rule
      FloatType tmp_a=0, tmp_b=0, h, start, tmp_it2;
      h = (high_sig - low_sig)*sito2/( 2.0*n + 1 );
      start = low_sig*sito2 + ito2;
      for (int ii=1 ; ii<n; ii++){
        tmp_it2 = start + ii*2.0*h ;
        tmp_a += part_1(ito1,sito1,ito2,sito2,tmp_it2,t);
        tmp_it2 = start + (ii*2.0+1.0)*h ;
        tmp_b += part_1(ito1,sito1,ito2,sito2,tmp_it2,t);
      }
      FloatType result=0;
      result = 4.0*tmp_b + 2.0*tmp_a;
      tmp_a = part_1(ito1,sito1,ito2,sito2,start,t);
      tmp_b = part_1(ito1,sito1,ito2,sito2,start+h*(2.0*n+2.0),t);
      result = result + tmp_a+tmp_b;
      result = h*(result)/3.0;
      return ( result );
    }


    FloatType num_int_hermite( FloatType const& ito1,    FloatType const& sito1,
                               FloatType const& ito2,    FloatType const& sito2,
                               FloatType const& t)
    {
      FloatType result=0, tmp_it2, p, s2=std::sqrt(2.0);

      for (int ii=0 ; ii<4; ii++){
        if(ii<x_.size()) {
          tmp_it2 = x_[ii]*sito2*s2+ito2;
          p=0.0;
          if (tmp_it2>=0){
            p = part_1(ito1,sito1,ito2,sito2,tmp_it2,t);
          }
          if(ii<wexs_.size()) p*=wexs_[ii];
          result+=p;
        }
      }
      return ( result*sito2*s2 );
    }



    FloatType p_raw(FloatType const& it1, FloatType const& it2, FloatType const& t)
    {
       FloatType result=0;
       if (it2 <= (1-t)/t*it1){
         if (it2 >= t/(1-t)*it1){
           result = std::exp( -it2-it1 )/(1-2*t);
         }
       }
       return (result);
    }


    FloatType log_p_given_t(FloatType const& t, int const& n)
    {
       FloatType result=0, tmp_result, ito1,sito1,ito2,sito2;
       long tmp_index;
       for (long ii=0;ii<z_.size();ii++){
         tmp_index = twin_mate_index_[ii];
         if (tmp_index>=0){ // this means twin mate is available
           ito1 = z_[ii];
           sito1 = sig_z_[ii];
           ito2 = z_[tmp_index];
           sito2 = sig_z_[tmp_index];
           tmp_result = num_int(ito1,sito1,ito2,sito2,-5,5,t,n);
           if (tmp_result > 0){
             result += std::log( tmp_result );
           }
           else{
             result += std::log( 1e-36 );
           }

         }

       }
       return(result);
    }

    FloatType fast_log_p_given_t(FloatType const& t)
    {
       FloatType result=0, tmp_result, ito1,sito1,ito2,sito2;
       long tmp_index;
       for (long ii=0;ii<z_.size();ii++){
         if(ii<twin_mate_index_.size()&&ii<sig_z_.size()) {
           tmp_index = twin_mate_index_[ii];
           if (tmp_index>=0){ // this means twin mate is available
             ito1 = z_[ii];
             sito1 = sig_z_[ii];
             if(tmp_index<z_.size()&&tmp_index<sig_z_.size()) {
               ito2 = z_[tmp_index];
               sito2 = sig_z_[tmp_index];
               tmp_result = num_int_hermite(ito1,sito1,ito2,sito2,t);
               if (tmp_result > 0){
                 result += std::log( tmp_result );
               }
               else{
                 result += std::log( 1e-36 );
               }
             }
           }
         }
       }
       return(result);
    }


  protected:
    scitbx::af::shared< FloatType > z_;
    scitbx::af::shared< FloatType > sig_z_;
    scitbx::af::shared< long > twin_mate_index_;
    very_quick_erf<FloatType> qerf_;
    long n_hermite_;
    scitbx::af::shared< FloatType > x_;
    scitbx::af::shared< FloatType > wexs_;

  };




  template <typename FloatType>
  class ml_twin_with_ncs
  {
  public:
    // input normalized structure factors.
    ml_twin_with_ncs( scitbx::af::const_ref<FloatType> const& z,
                      scitbx::af::const_ref<FloatType> const& sig_z,
                      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
                      scitbx::af::const_ref< std::size_t > const& bins,
                      cctbx::sgtbx::space_group const& space_group,
                      bool const& anomalous_flag,
                      scitbx::mat3<FloatType> const& twin_law,
                      cctbx::uctbx::unit_cell const& unit_cell,
                      long const& n_hermite):
      qei0_(5000),
      n_hermite_(n_hermite),
      sig_lim_( 3.5 )
    {
      SCITBX_ASSERT( z.size() == sig_z.size() );
      SCITBX_ASSERT( z.size() == hkl.size() );
      SCITBX_ASSERT( z.size() == bins.size() );

      cctbx::miller::index<> tmp_index;
      cctbx::miller::lookup_utils::lookup_tensor<FloatType>
        tmp_lut(hkl, space_group, anomalous_flag);

      long tmp_bin;
      for (std::size_t ii=0;ii<z.size();ii++){
        tmp_index = twin_mate( hkl[ii], twin_law );
        twin_mate_index_.push_back( tmp_lut.find_hkl( tmp_index ) );
        z_.push_back( z[ii] );
        tmp_bin = bins[ii];
        if (tmp_bin>0){
          tmp_bin-=1;
        }
        bins_.push_back( tmp_bin );
        sig_z_.push_back( sig_z[ii] );
        d_star_sq_.push_back( unit_cell.d_star_sq( hkl[ii] ) );
      }

      // settting up the gauss hermite quadrature
      scitbx::math::quadrature::gauss_hermite_engine<FloatType> tmp_ghe(n_hermite_);
      xh_ = tmp_ghe.x();
      wh_ = tmp_ghe.w();

      // settting up the gauss hermite quadrature
      scitbx::math::quadrature::gauss_legendre_engine<FloatType> tmp_gle(n_hermite_);
      xl_ = tmp_gle.x();
      wl_ = tmp_gle.w();

      sqrt2_ = std::sqrt(2.0);
    }


    FloatType p_raw(FloatType const& z1, FloatType const& z2,
                    FloatType const& D_ncs, FloatType const& t)
    {
      //std::cout << "INHERE " << std::endl;
      FloatType part_1, part_2, x, divide_by, result=0;
      if (z2 <= (1-t)/t*z1){
         if (z2 >= t/(1-t)*z1){
           part_1 = std::exp( -(z1+z2)/(1-D_ncs*D_ncs) )/( (1-D_ncs*D_ncs)*(1-2.0*t) );

           divide_by = (1-2*t)*(1.0 - D_ncs*D_ncs );

           x = std::sqrt( ( (z1*(1-t) -t*z2) )*
                          ( (z2*(1-t) -t*z1) )
                          )*2.0*D_ncs/( divide_by );
           part_2 = qei0_.ei0(x)*std::exp(x);
           result = part_1*part_2;
         }
      }
      return( result );
    }


    FloatType num_int_z2( FloatType const& zo1, FloatType const& so1,
                          FloatType const& zo2, FloatType const& so2,
                          FloatType const& t,   FloatType const& D_ncs,
                          FloatType const& z1)
    {
      // please set up the limits of z2 for a given value of z1
      FloatType z2_low, z2_high, o_low, o_high;
      z2_low = t/(1-t)*z1;
      z2_high = (1-t)/t*z1;

      o_low = zo2 - sig_lim_*so2;
      o_high = zo2 + sig_lim_*so2;
      // now we have to modify the integration domain
      bool hermite=true; // swith maybe to gauss hermite quadrature?

      // modify upper and lower limits if it goes beyond the limits
      if ( o_low < z2_low ){
        hermite=false;
        o_low = z2_low;
        if (o_high < z2_low){// the intersting limit lies completely ouitside the domain of interest
          return( 0 );  // return 0
        }
      }
      if ( o_high > z2_high ){
        hermite=false;
        o_high = z2_high;
        if ( o_low > z2_high ){ // the intersting limit lies completely ouitside the domain of interest
          return(0); // return 0
        }
      }

      // the [o_low o_high] range lies completely in the domain [z2_low z2_high]
      // use guass hermite quadrature rule please

      FloatType result=0, z2;
      FloatType mid_point=(o_high+o_low)/2.0;
      FloatType span = (o_high-o_low)/2.0;

      if (hermite){
        for (int ii=0;ii<xh_.size();ii++){
          z2 = xh_[ii]*so2*sqrt2_+zo2;
          result+=p_raw(z1,z2,D_ncs,t)*wh_[ii];
        }
        result *= sqrt2_*so2;
        return(result);
      }

      if (!hermite){
        for (int ii=0;ii<xl_.size();ii++){
          z2=xl_[ii]*span + mid_point;
          result+= std::exp( -(z2-zo2)*(z2-zo2) /(2.0*so2*so2) )*p_raw(z1,z2,D_ncs,t)*wl_[ii];
        }
        result*=span;
        return(result);
      }

      return(0); // this should never happen!
    }

    FloatType num_int(FloatType const& zo1, FloatType const& so1,
                      FloatType const& zo2, FloatType const& so2,
                      FloatType const& t,   FloatType const& D_ncs
                     )
    {
      FloatType result = 0;
       // we have to figure out what we want to do
      FloatType z1_low, o_low, o_high;
      z1_low = 0;
      o_low = zo1 - sig_lim_*so1;
      o_high = zo1 + sig_lim_*so1;

      bool hermite=true;


      if ( o_high < z1_low ){
        return( 1e-32 );
      }
      if (o_low < z1_low ){
        o_low = z1_low;
        hermite=false;
      }


      FloatType mid_point=(o_high+o_low)/2.0;
      FloatType span = (o_high-o_low)/2.0;
      FloatType z1;

      if (hermite){ // use gauss-hermite integration
        for (int ii=0;ii<xh_.size();ii++){
          z1 = xh_[ii]*so1*sqrt2_+zo1;
          result+=num_int_z2(zo1,so1,zo2,so2,t,D_ncs,z1)*wh_[ii];
        }
        result*=so1*sqrt2_;
      }

      if(!hermite){
        for (int ii=0;ii<xl_.size();ii++){
          z1=xl_[ii]*span + mid_point;
          result+= std::exp( -(z1-zo1)*(z1-zo1)/(2.0*so1*so1) )*num_int_z2(zo1,so1,zo2,so2,t,D_ncs,z1)*wl_[ii];
        }
        result*=span;
      }
      if (result < 1e-32 ){
        result = 1e-32;
      }
      return (result);
    }

    scitbx::af::shared<FloatType> p_tot_given_t_and_coeff( FloatType const& t_coeff,
                                                           scitbx::af::const_ref<FloatType> d_coeffs )
    {
      // we use D_ncs = exp( -(coeff*coeff + 0.01) )
      FloatType result=0,tmp_result, z1,z2,s1,s2,D_ncs;
      FloatType twin_fraction,twin_fraction_tmp;
      FloatType d_t, td_t=0, tmp1;

      // setting up arrays for fd
      scitbx::af::shared<FloatType> score_coeffs;
      scitbx::af::shared<FloatType> score_coeffs_fd;
      for (int ii=0; ii<d_coeffs.size();ii++){
        score_coeffs.push_back( 0.0 );
        score_coeffs_fd.push_back( 0.0 );
      }


      long tmp_index;
      FloatType h=0.0001;

      FloatType twin_fraction_cap=0.48;
      twin_fraction = twin_fraction_cap/(1.0+std::exp(-t_coeff));
      twin_fraction_tmp = twin_fraction_cap/(1.0+std::exp( -(t_coeff+h) ));

      FloatType d_ncs_cap=0.95;

      for (long ii=0;ii<z_.size();ii++){
        tmp_index = twin_mate_index_[ii];
        if ( tmp_index >=0 ){
          z1=z_[ii];
          s1=sig_z_[ii];
          z2=z_[tmp_index];
          s2=sig_z_[tmp_index];

          D_ncs = d_ncs_cap/(1.0+std::exp(-d_coeffs[ bins_[ii] ]) );
          tmp_result = std::log( num_int(z1,s1,
                                         z2,s2,
                                         twin_fraction, D_ncs)
                                 );

          // fin diffs for twin fraction parameter
          d_t = std::log( num_int(z1,s1,
                                  z2,s2,
                                  twin_fraction_tmp, D_ncs)
                          );

          // fin difss for d_coeffs
          score_coeffs[ bins_[ii]  ]+=tmp_result;

          D_ncs = d_ncs_cap/(1.0+std::exp(-(d_coeffs[ bins_[ii] ] + h ) ) );
          tmp1 = std::log( num_int(z1,s1,
                                   z2,s2,
                                   twin_fraction, D_ncs)
                           );
          score_coeffs_fd[ bins_[ii]  ]+=tmp1;
          td_t += d_t ;
          result+= tmp_result;
        }
      }


      // get the gradients
      scitbx::af::shared<FloatType> result_vector;

      // compute fd's
      result_vector.push_back( -result );
      result_vector.push_back(( result-td_t  )/h );

      for (int ii=0;ii<d_coeffs.size();ii++){
        result_vector.push_back( -(score_coeffs_fd[ii]-score_coeffs[ii])/h  );
      }

      return (result_vector);
    }





  protected:
    scitbx::af::shared< FloatType > z_;
    scitbx::af::shared< FloatType > sig_z_;
    scitbx::af::shared< FloatType > d_star_sq_;
    scitbx::af::shared< long > bins_;
    scitbx::af::shared< long > twin_mate_index_;

    quick_ei0<FloatType> qei0_;

    long n_hermite_;
    scitbx::af::shared< FloatType > xh_;
    scitbx::af::shared< FloatType > wh_;

    scitbx::af::shared< FloatType > xl_;
    scitbx::af::shared< FloatType > wl_;

    FloatType sqrt2_;
    FloatType sig_lim_;

  };
















}}} // namespace mmtbx::scaling::detwinning
#endif // MMTBX_SCALING_DETWINING_H
