//! Peter Zwart April 19, 2005
#ifndef MMTBX_SCALING_RELATIVE_SCALING_H
#define MMTBX_SCALING_RELATIVE_SCALING_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/lookup_utils.h>
#include <mmtbx/scaling/scaling.h>
#include <mmtbx/scaling/absolute_scaling.h>
#include <mmtbx/scaling/twinning.h>
#include <scitbx/math/bessel.h>



namespace mmtbx { namespace scaling{
namespace relative_scaling{


  //! Computes LS relative scaling target function for a given pair of I's
  /*! Relative scaling is always performed anisotropically
   *  for obvious reasons. The isotropic case is NOT
   *  available.
   *
   *  Note that (in contrast with absolute scaling)
   *  the relative scaling is performed on intensities,as
   *  a Gaussian approximation is better justified in this case.
   *
   *  The scale will be applied on iobs2.
   */
  template <typename FloatType>
  FloatType
  rel_scale_single_ls_target(cctbx::miller::index<> const& hkl,
                             FloatType const& iobs1,
                             FloatType const& sigma_iobs1,
                             FloatType const& iobs2,
                             FloatType const& sigma_iobs2,
                             FloatType const& p_scale,
                             cctbx::uctbx::unit_cell const& uc,
                             scitbx::sym_mat3<FloatType> const& u)
  {
    typedef FloatType f_t;
    f_t result = 0.0;
    f_t V_star = pow( 1.0/(uc.volume()), 2.0/3.0);//RWGK's magic scalar
    // The scale factor is the sqrts from the one used in absolute scaling
    f_t k = mmtbx::scaling::absolute_scaling::wilson_get_aniso_scale(
                                               hkl,p_scale,V_star,u);
    result = (iobs1-k*k*iobs2)*(iobs1-k*k*iobs2)/
      (sigma_iobs1*sigma_iobs1+k*k*k*k*sigma_iobs2*sigma_iobs2);
    return(result);
  }

  //! Array version of rel_scale_single_ls_target.
  /*! See rel_scale_single_ls_target for details.
   */
  template <typename FloatType>
  FloatType
  rel_scale_total_ls_target(
            scitbx::af::const_ref<cctbx::miller::index<> > const& hkl,
            scitbx::af::const_ref<FloatType> const& iobs1,
            scitbx::af::const_ref<FloatType> const& sigma_iobs1,
            scitbx::af::const_ref<FloatType> const& iobs2,
            scitbx::af::const_ref<FloatType> const& sigma_iobs2,
            FloatType const& p_scale,
            cctbx::uctbx::unit_cell const& uc,
            scitbx::sym_mat3<FloatType> const& u)
  {

    SCITBX_ASSERT(hkl.size() == iobs1.size());
    SCITBX_ASSERT(hkl.size() == iobs2.size());
    SCITBX_ASSERT(hkl.size() == sigma_iobs1.size());
    SCITBX_ASSERT(hkl.size() == sigma_iobs2.size());

    typedef FloatType f_t;
    f_t result = 0.0;
    unsigned n_points=hkl.size();
    for (unsigned ii=0; ii<n_points; ii++){
      result +=rel_scale_single_ls_target(hkl[ii],
                                          iobs1[ii],
                                          sigma_iobs1[ii],
                                          iobs2[ii],
                                          sigma_iobs2[ii],
                                          p_scale,
                                          uc,
                                          u);
    }
    return(result);

  }

  //! Gradient of the LS target function
  template <typename FloatType>
  scitbx::af::shared<FloatType>
  rel_scale_single_ls_gradient(cctbx::miller::index<> const& hkl,
                               FloatType const& iobs1,
                               FloatType const& sigma_iobs1,
                               FloatType const& iobs2,
                               FloatType const& sigma_iobs2,
                               FloatType const& p_scale,
                               cctbx::uctbx::unit_cell const& uc,
                               scitbx::sym_mat3<FloatType> const& u)
  {
    typedef FloatType f_t;
    scitbx::af::shared<FloatType> result(7,0); // scale factor + u tensor
    f_t V_star = pow( 1.0/(uc.volume()), 2.0/3.0);//RWGK's magic scalar
    // The scale factor is the square from the one used in absolute scaling
    f_t k = mmtbx::scaling::absolute_scaling::wilson_get_aniso_scale(
                                               hkl,p_scale,V_star,u);
    f_t tmp = -4.0*k*iobs2*(iobs1-k*k*iobs2)/
      (sigma_iobs1*sigma_iobs1+k*k*k*k*sigma_iobs2*sigma_iobs2);

    f_t tmp_b = -4.0*k*k*(iobs1-k*k*iobs2)*(iobs1-k*k*iobs2)*
      sigma_iobs2*sigma_iobs2/
      (
       (sigma_iobs1*sigma_iobs1+k*k*k*k*sigma_iobs2*sigma_iobs2)*
       (sigma_iobs1*sigma_iobs1+k*k*k*k*sigma_iobs2*sigma_iobs2)
      );


    tmp = tmp+tmp_b;
    f_t tmp2 = scitbx::constants::pi*scitbx::constants::pi*V_star;

    result[0]=-tmp*k;
    result[1]=2.0*tmp*k*tmp2*hkl[0]*hkl[0];
    result[2]=2.0*tmp*k*tmp2*hkl[1]*hkl[1];
    result[3]=2.0*tmp*k*tmp2*hkl[2]*hkl[2];
    result[4]=4.0*tmp*k*tmp2*hkl[0]*hkl[1];
    result[5]=4.0*tmp*k*tmp2*hkl[0]*hkl[2];
    result[6]=4.0*tmp*k*tmp2*hkl[1]*hkl[2];

    return(result);
  }

  //! Array version of the gradient calculations
  template <typename FloatType>
  scitbx::af::shared<FloatType>
  rel_scale_total_ls_gradient(
            scitbx::af::const_ref<cctbx::miller::index<> > const& hkl,
            scitbx::af::const_ref<FloatType> const& iobs1,
            scitbx::af::const_ref<FloatType> const& sigma_iobs1,
            scitbx::af::const_ref<FloatType> const& iobs2,
            scitbx::af::const_ref<FloatType> const& sigma_iobs2,
            FloatType const& p_scale,
            cctbx::uctbx::unit_cell const& uc,
            scitbx::sym_mat3<FloatType> const& u)


  {
    SCITBX_ASSERT(hkl.size() == iobs1.size());
    SCITBX_ASSERT(hkl.size() == iobs2.size());
    SCITBX_ASSERT(hkl.size() == sigma_iobs1.size());
    SCITBX_ASSERT(hkl.size() == sigma_iobs2.size());

    scitbx::af::shared<FloatType> result(7,0); // scale factor + u tensor
    scitbx::af::shared<FloatType> tmp(7,0);

    unsigned n_points=hkl.size();
    for (unsigned ii=0; ii<n_points; ii++){
      tmp =rel_scale_single_ls_gradient(hkl[ii],
                                        iobs1[ii],
                                        sigma_iobs1[ii],
                                        iobs2[ii],
                                        sigma_iobs2[ii],
                                        p_scale,
                                        uc,
                                        u);
      result[0]+=tmp[0];
      result[1]+=tmp[1];
      result[2]+=tmp[2];
      result[3]+=tmp[3];
      result[4]+=tmp[4];
      result[5]+=tmp[5];
      result[6]+=tmp[6];
    }
    return(result);
  }



  template <typename FloatType=double>
  class least_squares_on_i_proper
  {
    public:
    /*! Default constructor */
    least_squares_on_i_proper() {};
    /*! Pass in all that is needed */
    least_squares_on_i_proper(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,     
      scitbx::af::const_ref< FloatType > const& i_nat,
      scitbx::af::const_ref< FloatType > const& sig_nat,
      scitbx::af::const_ref< FloatType > const& i_der,
      scitbx::af::const_ref< FloatType > const& sig_der,
      FloatType const& p_scale,
      cctbx::uctbx::unit_cell const& unit_cell,
      scitbx::sym_mat3<FloatType> const& u_rwgk)
    :    
    p_scale_(p_scale),
    unit_cell_( unit_cell ),
    u_rwgk_( u_rwgk ),
    scale_rwgk_( pow( unit_cell_.volume(), -2.0/3.0) )
    {
      // please check if all arrays have the same number of members      
      SCITBX_ASSERT ( hkl.size() ==  i_nat.size() );
      SCITBX_ASSERT ( hkl.size() ==  sig_nat.size() );
      SCITBX_ASSERT ( hkl.size() ==  i_der.size() );
      SCITBX_ASSERT ( hkl.size() ==  sig_der.size() );
      
      // fill up the arrays please
      for (unsigned ii=0;ii<hkl.size();ii++){
	hkl_.push_back( hkl[ii] );
	i_nat_.push_back( i_nat[ii] );
	sig_nat_.push_back( sig_nat[ii] );
	i_der_.push_back( i_der[ii] );
	sig_der_.push_back( sig_der[ii] );
      }
      // that should be it.   
    }


    
    //------------------------
    scitbx::af::shared<FloatType> 
    get_gradient()
    {
      return( gradient() );
    }
    
    FloatType
    get_function()
    {
      return( function() );
    }    
    //------------------------

    void
    set_p_scale(FloatType p_scale)
    {
      p_scale_=p_scale;
    }

    void 
    set_u_rwgk( scitbx::sym_mat3<FloatType> const& u_rwgk )
    {
      u_rwgk_[0] = u_rwgk[0];
      u_rwgk_[1] = u_rwgk[1];
      u_rwgk_[2] = u_rwgk[2];
      u_rwgk_[3] = u_rwgk[3];
      u_rwgk_[4] = u_rwgk[4];    
      u_rwgk_[5] = u_rwgk[5];      
    }

    void 
    set_params(FloatType p_scale,
	       scitbx::sym_mat3<FloatType> const& u_rwgk )
    {
      set_p_scale(p_scale);
      set_u_rwgk(u_rwgk);
    }

    
    protected:

    FloatType 
    get_k( unsigned index )
    {
      FloatType result,h,k,l;
      h = hkl_[index][0];
      k = hkl_[index][1];      
      l = hkl_[index][2];
      
      result = h*(h*u_rwgk_[0] + k*u_rwgk_[3] + l*u_rwgk_[4] )+
               k*(h*u_rwgk_[3] + k*u_rwgk_[1] + l*u_rwgk_[5] )+
               l*(h*u_rwgk_[4] + k*u_rwgk_[5] + l*u_rwgk_[2] );
      result = result*scitbx::constants::pi*scitbx::constants::pi
               *2.0*scale_rwgk_-p_scale_;   
      result = std::exp(result);
      return( result );     
    }
    
    FloatType 
    function( unsigned index )
    {
      FloatType result=0, k, weight;
      k=get_k( index );
      result = i_nat_[index] - i_der_[index]*k*k;
      result=result*result;
      weight = sig_nat_[index]*sig_nat_[index];
      weight = weight+sig_der_[index]*sig_der_[index]*k*k*k*k;
      result = result/weight;
      return(result);
    }

    FloatType
    function()
    {
      FloatType result=0;
      for (unsigned ii=0;ii<hkl_.size();ii++){
	result += function(ii);
      }
      return (result);
    }

    scitbx::af::shared<FloatType> 
    gradient( unsigned index )
    {
      scitbx::af::shared<FloatType> single_grad(7,0);
      FloatType tmp1,tmp2,tmp3,tmp4, result, k;
      k=get_k( index );

      tmp3 = -k*k*i_der_[index] + i_nat_[index];
      tmp4 = k*k*k*k*sig_der_[index]*sig_der_[index]
                   + sig_nat_[index]*sig_nat_[index];
      
      tmp1 = -4.0*k*k*tmp3*tmp3*sig_der_[index]*sig_der_[index]/(tmp4*tmp4);
      tmp2 = -4.0*i_der_[index]*k*tmp3/tmp4;

      result = tmp1+tmp2;

      FloatType tmp_const=2.0*scitbx::constants::pi*scitbx::constants::pi*
        scale_rwgk_;

      single_grad[0] = result* 
        -k;      
      single_grad[1] = result* 
        tmp_const*hkl_[index][0]*hkl_[index][0];
      single_grad[2] = result* 
        tmp_const*hkl_[index][1]*hkl_[index][1];
      single_grad[3] = result* 
        tmp_const*hkl_[index][1]*hkl_[index][2];
      single_grad[4] = result*
        tmp_const*hkl_[index][0]*hkl_[index][1]*2.0;      
      single_grad[5] = result*
        tmp_const*hkl_[index][0]*hkl_[index][2]*2.0;
      single_grad[6] = result* 
        tmp_const*hkl_[index][1]*hkl_[index][2]*2.0;      

      return( single_grad );
    }
    
    scitbx::af::shared<FloatType> 
    gradient()
    {
      scitbx::af::shared<FloatType> result(7,0);
      scitbx::af::shared<FloatType> tmp_gradient(7,0);
      for (unsigned ii=0;ii<hkl_.size();ii++){
	tmp_gradient = gradient( ii );
	for (unsigned jj=0; jj<7; jj++){
	  result[jj]+=tmp_gradient[jj];
	}
      }
      return(result);    
    }

    
    
    scitbx::af::shared< cctbx::miller::index<> > hkl_;
    scitbx::af::shared< FloatType > i_nat_;
    scitbx::af::shared< FloatType > sig_nat_;
    scitbx::af::shared< FloatType > i_der_;
    scitbx::af::shared< FloatType > sig_der_;
    FloatType p_scale_;
    cctbx::uctbx::unit_cell unit_cell_;
    scitbx::sym_mat3<FloatType> u_rwgk_;
    FloatType scale_rwgk_;

  };
    

















  // This might be a more appropriate target function
    template <typename FloatType=double>
  class least_squares_on_i
  {
    public:
    /*! Default constructor */
    least_squares_on_i() {};
    /*! Pass in all that is needed */
    least_squares_on_i(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,     
      scitbx::af::const_ref< FloatType > const& i_nat,
      scitbx::af::const_ref< FloatType > const& sig_nat,
      scitbx::af::const_ref< FloatType > const& i_der,
      scitbx::af::const_ref< FloatType > const& sig_der,
      FloatType const& p_scale,
      cctbx::uctbx::unit_cell const& unit_cell,
      scitbx::sym_mat3<FloatType> const& u_rwgk)
    :    
    p_scale_(p_scale),
    unit_cell_( unit_cell ),
    // u_rwgk_(u_rwgk[0],u_rwgk[1],u_rwgk[2],u_rwgk[3],u_rwgk[4],u_rwgk[5])
    u_rwgk_( u_rwgk ),
    scale_rwgk_( pow( unit_cell_.volume(), -2.0/3.0) )
    {
      // please check if all arrays have the same number of members      
      SCITBX_ASSERT ( hkl.size() ==  i_nat.size() );
      SCITBX_ASSERT ( hkl.size() ==  sig_nat.size() );
      SCITBX_ASSERT ( hkl.size() ==  i_der.size() );
      SCITBX_ASSERT ( hkl.size() ==  sig_der.size() );
      
      // fill up the arrays please
      for (unsigned ii=0;ii<hkl.size();ii++){
	hkl_.push_back( hkl[ii] );
	i_nat_.push_back( i_nat[ii] );
	sig_nat_.push_back( sig_nat[ii] );
	i_der_.push_back( i_der[ii] );
	sig_der_.push_back( sig_der[ii] );
      }
      // that should be it.   
    }


    
    //------------------------
    scitbx::af::shared<FloatType> 
    get_gradient()
    {
      return( gradient() );
    }
    
    FloatType
    get_function()
    {
      return( function() );
    }    
    //------------------------

    void
    set_p_scale(FloatType p_scale)
    {
      p_scale_=p_scale;
    }

    void 
    set_u_rwgk( scitbx::sym_mat3<FloatType> const& u_rwgk )
    {
      u_rwgk_[0] = u_rwgk[0];
      u_rwgk_[1] = u_rwgk[1];
      u_rwgk_[2] = u_rwgk[2];
      u_rwgk_[3] = u_rwgk[3];
      u_rwgk_[4] = u_rwgk[4];    
      u_rwgk_[5] = u_rwgk[5];      
    }

    void 
    set_params(FloatType p_scale,
	       scitbx::sym_mat3<FloatType> const& u_rwgk )
    {
      set_p_scale(p_scale);
      set_u_rwgk(u_rwgk);
    }

    
    protected:

    FloatType 
    get_k( unsigned index )
    {
      FloatType result,h,k,l;
      h = hkl_[index][0];
      k = hkl_[index][1];      
      l = hkl_[index][2];
      
      result = h*(h*u_rwgk_[0] + k*u_rwgk_[3] + l*u_rwgk_[4] )+
               k*(h*u_rwgk_[3] + k*u_rwgk_[1] + l*u_rwgk_[5] )+
               l*(h*u_rwgk_[4] + k*u_rwgk_[5] + l*u_rwgk_[2] );
      result = result*scitbx::constants::pi*scitbx::constants::pi
               *2.0*scale_rwgk_-p_scale_;   
      result = std::exp(result);
      return( result );     
    }
    
    FloatType 
    function( unsigned index )
    {
      FloatType result=0, k, weight;
      k=get_k( index );
      weight = i_nat_[index]*sig_der_[index]/i_der_[index];
      weight = weight*weight + sig_nat_[index]*sig_nat_[index];
      result = i_nat_[index] - i_der_[index]*k*k;
      result = result*result;
      result = result/weight;
      return(result);
    }

    FloatType
    function()
    {
      FloatType result=0;
      for (unsigned ii=0;ii<hkl_.size();ii++){
	result += function(ii);
      }
      return (result);
    }

    scitbx::af::shared<FloatType> 
    gradient( unsigned index )
    {
      scitbx::af::shared<FloatType> single_grad(7,0);
      FloatType result=0, k, weight;
      k=get_k( index );
      weight = i_nat_[index]*sig_der_[index]/i_der_[index];
      weight = weight*weight + sig_nat_[index]*sig_nat_[index];
      weight = 1.0/weight;
      result = i_nat_[index] - i_der_[index]*k*k;
      
      FloatType tmp_const=2.0*scitbx::constants::pi*scitbx::constants::pi*
        scale_rwgk_;

      single_grad[0] = 2.0*result*weight* -2.0*k*i_der_[index]* 
        -k;      
      single_grad[1] = 2.0*result*weight* -2.0*k*i_der_[index]* 
        tmp_const*hkl_[index][0]*hkl_[index][0];
      single_grad[2] = 2.0*result*weight* -2.0*k*i_der_[index]* 
        tmp_const*hkl_[index][1]*hkl_[index][1];
      single_grad[3] = 2.0*result*weight* -2.0*k*i_der_[index]* 
        tmp_const*hkl_[index][1]*hkl_[index][2];
      single_grad[4] = 2.0*result*weight* -2.0*k*i_der_[index]*
        tmp_const*hkl_[index][0]*hkl_[index][1]*2.0;      
      single_grad[5] = 2.0*result*weight* -2.0*k*i_der_[index]*
        tmp_const*hkl_[index][0]*hkl_[index][2]*2.0;
      single_grad[6] = 2.0*result*weight* -2.0*k*i_der_[index]* 
        tmp_const*hkl_[index][1]*hkl_[index][2]*2.0;      

      return( single_grad );
    }
    
    scitbx::af::shared<FloatType> 
    gradient()
    {
      scitbx::af::shared<FloatType> result(7,0);
      scitbx::af::shared<FloatType> tmp_gradient(7,0);
      for (unsigned ii=0;ii<hkl_.size();ii++){
	tmp_gradient = gradient( ii );
	for (unsigned jj=0; jj<7; jj++){
	  result[jj]+=tmp_gradient[jj];
	}
      }
      return(result);    
    }

    scitbx::af::shared< cctbx::miller::index<> > hkl_;
    scitbx::af::shared< FloatType > i_nat_;
    scitbx::af::shared< FloatType > sig_nat_;
    scitbx::af::shared< FloatType > i_der_;
    scitbx::af::shared< FloatType > sig_der_;
    FloatType p_scale_;
    cctbx::uctbx::unit_cell unit_cell_;
    scitbx::sym_mat3<FloatType> u_rwgk_;
    FloatType scale_rwgk_;

  };
    





































  //------------------------------------------------------------------
  // This routine is need for generating a 'property'
  // needed for the local_area determination
  //
  // The property that is generated in this instance
  // makes sure that only neighbours are generated
  // that are both present in the natiove and derivative
  

  template <typename FloatType=double>
  class property_matching_indices
  {
    public:
    /*! Default constructor */
    property_matching_indices(){}
    /*! This constyructuc a property (matching indices) */
    property_matching_indices(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_master,
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_set,
      cctbx::sgtbx::space_group const& space_group,
      bool const& anomalous_flag)
    :
    master_lookup_( hkl_master, space_group , anomalous_flag ),
    combined_property_( hkl_master.size(), false )
    {

      set_lut_ = master_lookup_.find_hkl( hkl_set );
      // we now have to lookup tables. Please make sure
      // that all reflections are present in the master
      // array
      for (unsigned ii=0;ii<hkl_set.size();ii++){
        SCITBX_ASSERT( set_lut_[ii]>=0 ); // set_a has to be a subset of
        combined_property_[set_lut_[ii] ]=true;
      }

    }

    cctbx::miller::lookup_utils::lookup_tensor<FloatType> master_lookup_;
    scitbx::af::shared< long > set_lut_; // lookup table
    scitbx::af::shared< bool > combined_property_;


  };



  //-----------------------------------------------------------
  // local scaling routines
  template <typename FloatType=double>
  class local_scaling_moment_based
  {
    public:
    /*! default constructor */
    local_scaling_moment_based(){}
    /*! This does (allmost) everything */
    local_scaling_moment_based(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_master,
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_sets,

      scitbx::af::const_ref< FloatType > const& data_set_a,
      scitbx::af::const_ref< FloatType > const& sigma_set_a,
      scitbx::af::const_ref< FloatType > const& data_set_b,
      scitbx::af::const_ref< FloatType > const& sigma_set_b,

      cctbx::sgtbx::space_group const& space_group,
      bool const& anomalous_flag,
      std::size_t const& radius,
      std::size_t const& depth,
      std::size_t const& at_least_this_number_of_neighbours,

      bool const& weights)
    :
    use_weights_( weights ),
    local_scales_( hkl_sets.size(), 1),
    local_scales_sigmas_( hkl_sets.size(), 0.01),
    property_generator_( hkl_master,
                         hkl_sets,
                         space_group,
                         anomalous_flag ),
    area_generator_( hkl_master,
                     property_generator_.combined_property_.const_ref(),
                     space_group,
                     anomalous_flag,
                     radius,
                     depth,
                     at_least_this_number_of_neighbours )
    {

      for (unsigned ii=0;ii<hkl_master.size();ii++){
        hkl_master_.push_back( hkl_master[ii] );
      }

      for (unsigned ii=0;ii<hkl_sets.size();ii++){
        hkl_sets_.push_back( hkl_sets[ii] );
        data_set_a_.push_back( data_set_a[ii] );
        sigma_set_a_.push_back( sigma_set_a[ii] );
        data_set_b_.push_back( data_set_b[ii] );
        sigma_set_b_.push_back( sigma_set_b[ii] );
      }

      // now eveything is set, please do the local scaling
      scale_it();
    }

    void
    scale_it()
    {
      // loop over all reflections of the native
      unsigned nb_index ;
      FloatType multiplier=0.0, constant=1.0, weight;
      if (use_weights_){
        multiplier=1.0;
        constant=0.0;
      }

      for (unsigned ii=0;ii<hkl_sets_.size();ii++){


        FloatType nat_mean=0, der_mean=0;
        FloatType nat_var=0, der_var=0;
        // compute the mean for the native and derivativesummats

        for (unsigned jj=1;
             jj<area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
             jj++){

          nb_index = area_generator_.area_[
            property_generator_.set_lut_[ ii ] ][jj];

          weight = sigma_set_b_[ nb_index ]* data_set_a_[ nb_index  ]/
            data_set_b_[ nb_index  ];
          weight = weight*weight
            + sigma_set_a_[ nb_index ]*sigma_set_a_[ nb_index ];
          weight=weight*multiplier+constant;
	  weight=1.0/weight;

          nat_mean += data_set_a_[ nb_index  ]*weight;
          nat_var += (data_set_a_[ nb_index  ]*weight*
                      data_set_a_[ nb_index  ]*weight);

          der_mean += data_set_b_[ nb_index  ]*weight;
          der_var += (data_set_b_[ nb_index  ]*weight*
                      data_set_b_[ nb_index  ]*weight);
        }
        if (der_mean>0){
          local_scales_[ii]=nat_mean/der_mean;
        }

      }

    }

    scitbx::af::shared<FloatType>
    get_scales()
    {
      return( local_scales_ );
    }



    protected:
    property_matching_indices<FloatType> property_generator_;
    cctbx::miller::lookup_utils::local_area<FloatType> area_generator_;

    scitbx::af::shared< cctbx::miller::index<> > hkl_master_;
    scitbx::af::shared< cctbx::miller::index<> > hkl_sets_;

    scitbx::af::shared< FloatType > data_set_a_;
    scitbx::af::shared< FloatType > sigma_set_a_;
    scitbx::af::shared< FloatType > data_set_b_;
    scitbx::af::shared< FloatType > sigma_set_b_;

    scitbx::af::shared< FloatType > local_scales_;
    scitbx::af::shared< FloatType > local_scales_sigmas_;
    bool use_weights_;


  };







    //-----------------------------------------------------------
  // local scaling routines
  template <typename FloatType=double>
  class local_scaling_ls_based
  {
    public:
    /*! default constructor */
    local_scaling_ls_based(){}
    /*! This does (allmost) everything */
    local_scaling_ls_based(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_master,
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_sets,

      scitbx::af::const_ref< FloatType > const& data_set_a,
      scitbx::af::const_ref< FloatType > const& sigma_set_a,
      scitbx::af::const_ref< FloatType > const& data_set_b,
      scitbx::af::const_ref< FloatType > const& sigma_set_b,

      cctbx::sgtbx::space_group const& space_group,
      bool const& anomalous_flag,
      std::size_t const& radius,
      std::size_t const& depth,
      std::size_t const& at_least_this_number_of_neighbours,

      bool const& weights)
    :
    use_weights_( weights ),
    local_scales_( hkl_sets.size(), 1),
    local_scales_sigmas_( hkl_sets.size(), 0.01),
    property_generator_( hkl_master,
                         hkl_sets,
                         space_group,
                         anomalous_flag ),
    area_generator_( hkl_master,
                     property_generator_.combined_property_.const_ref(),
                     space_group,
                     anomalous_flag,
                     radius,
                     depth,
                     at_least_this_number_of_neighbours )
    {

      for (unsigned ii=0;ii<hkl_master.size();ii++){
        hkl_master_.push_back( hkl_master[ii] );
      }

      for (unsigned ii=0;ii<hkl_sets.size();ii++){
        hkl_sets_.push_back( hkl_sets[ii] );
        data_set_a_.push_back( data_set_a[ii] );
        sigma_set_a_.push_back( sigma_set_a[ii] );
        data_set_b_.push_back( data_set_b[ii] );
        sigma_set_b_.push_back( sigma_set_b[ii] );
      }

      // now eveything is set, please do the local scaling
      scale_it();
    }

    void
    scale_it()
    {
      // loop over all reflections of the native
      unsigned nb_index ;
      FloatType multiplier=0.0, constant=1.0, weight;
      if (use_weights_){
        multiplier=1.0;
        constant=0.0;
      }

      for (unsigned ii=0;ii<hkl_sets_.size();ii++){


        FloatType top=0, bottom=0;
        for (unsigned jj=1;
             jj<area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
             jj++){

          nb_index = area_generator_.area_[
            property_generator_.set_lut_[ ii ] ][jj];

          weight = sigma_set_b_[ nb_index ]* data_set_a_[ nb_index  ]/
            data_set_b_[ nb_index  ];
          weight = weight*weight
            + sigma_set_a_[ nb_index ]*sigma_set_a_[ nb_index ];
          weight=weight*multiplier+constant;
	  weight=1.0/weight;

	  top+=weight*data_set_a_[nb_index]*data_set_b_[nb_index];
	  bottom+=weight*data_set_b_[nb_index]*data_set_b_[nb_index];

         }
        if (bottom>0){
          local_scales_[ii]=top/bottom;
        }

      }

    }

    scitbx::af::shared<FloatType>
    get_scales()
    {
      return( local_scales_ );
    }

    protected:
    property_matching_indices<FloatType> property_generator_;
    cctbx::miller::lookup_utils::local_area<FloatType> area_generator_;

    scitbx::af::shared< cctbx::miller::index<> > hkl_master_;
    scitbx::af::shared< cctbx::miller::index<> > hkl_sets_;

    scitbx::af::shared< FloatType > data_set_a_;
    scitbx::af::shared< FloatType > sigma_set_a_;
    scitbx::af::shared< FloatType > data_set_b_;
    scitbx::af::shared< FloatType > sigma_set_b_;

    scitbx::af::shared< FloatType > local_scales_;
    scitbx::af::shared< FloatType > local_scales_sigmas_;
    bool use_weights_;


  };













}}} // namespace mmtbx::scaling::relative_scaling
#endif // MMTBX_SCALING_RELATIVE_SCALING_H
