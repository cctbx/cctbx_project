//! Peter Zwart April 21, 2005
#ifndef MMTBX_SCALING_TWINNING_H
#define MMTBX_SCALING_TWINNING_H

#include <cstdio>
#include <iostream>
#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/match_indices.h>
#include <cctbx/miller/asu.h>
#include <cctbx/miller/lookup_utils.h>
#include <cctbx/sgtbx/reciprocal_space_asu.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <scitbx/array_family/sort.h>
#include <scitbx/random.h>
#include <map>
#include <vector>


namespace mmtbx { namespace scaling {
namespace twinning {


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
    }

    FloatType r_value()
    {
      FloatType top=0, bottom=0;
      FloatType x1,x2;
      unsigned tmp_loc;
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
      FloatType result=0;
      if (top>0){
        if (bottom>0){
          result=top/bottom;
        }
      }
      return (result);
    }

    protected:
    scitbx::af::shared< cctbx::miller::index<> > hkl_;
    scitbx::af::shared< cctbx::miller::index<> > hkl_twin_;
    scitbx::af::shared< FloatType > intensity_;

    scitbx::af::shared< long > location_;
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
      FloatType P;

      P = (l_value*l_value-1)*
      (-1+l_value*l_value+2.0*alpha*
       (-1+alpha+(-1+alpha)*(3+4*(-1+alpha)*alpha)*l_value*l_value));
      P/=((-1+(1-2*alpha)*(1-2*alpha)*l_value*l_value)*
      (-1+(1-2*alpha)*(1-2*alpha)*l_value*l_value));

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

      for (unsigned ii=0;ii<hkl_.size();ii++){
        int tmp1 = ii;
        int tmp2 = location_[ii];
        if (tmp2 >=0){
          Iold1 = intensity_[tmp1];
          Iold2 = intensity_[tmp2];
          Sold1 = sigma_[tmp1];
          Sold2 = sigma_[tmp2];

          Inew = ((1.0-alpha)*Iold1 - alpha*Iold2)/(1-2.0*alpha);
          Snew = pow(((1-alpha)/(1-2.0*alpha)),2.0)*Sold1*Sold1
            + pow(((alpha)/(1-2.0*alpha)),2.0)*Sold2*Sold2;

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
      for (unsigned ii=0;ii<limit_slot;ii++){
        p = temp_p_values[ sort_permutation[ii] ];
        q = temp_q_values[ sort_permutation[ii] ];
        h_values_.append( std::abs(p-q) / (p+q) );
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






















}}} // namespace mmtbx::scaling::detwinning
#endif // MMTBX_SCALING_DETWINING_H
