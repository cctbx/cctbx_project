//! Peter Zwart April 19, 2005
#ifndef MMTBX_SCALING_RELATIVE_SCALING_H
#define MMTBX_SCALING_RELATIVE_SCALING_H

#include <scitbx/constants.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/versa.h>
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



  /* relative scaling via minimising
     Residual=sum (Inat-k*k*Ider)**2/(wnat+k**4*wder)
  */
  template <typename FloatType=double>
  class least_squares_on_i_wt
  {
    public:
    /*! Default constructor */
    least_squares_on_i_wt() {};
    /*! Pass in all that is needed */
    least_squares_on_i_wt(
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

    scitbx::af::shared<FloatType>
    get_gradient(unsigned index)
    {
      return( gradient(index) );
    }

    FloatType
    get_function(unsigned index)
    {
      return( function(index) );
    }

    scitbx::af::shared<FloatType>
    hessian_as_packed_u(unsigned index)
    {
      return( hessian(index) );
    }

    scitbx::af::shared<FloatType>
    hessian_as_packed_u()
    {
      return( hessian() );
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
      if(result > 40.) result = std::exp(40.);
      else result = std::exp(result);
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
      tmp1 = -4.0*k*k*k*tmp3*tmp3*sig_der_[index]*sig_der_[index]/(tmp4*tmp4);
      tmp2 = -4.0*i_der_[index]*k*tmp3/tmp4;

      result = tmp1+tmp2;
      FloatType tmp_const=2.0*scitbx::constants::pi*scitbx::constants::pi*
        scale_rwgk_;

      single_grad[0] = result*
        -k;
      single_grad[1] = result*
        tmp_const*hkl_[index][0]*hkl_[index][0]*k;
      single_grad[2] = result*
        tmp_const*hkl_[index][1]*hkl_[index][1]*k;
      single_grad[3] = result*
        tmp_const*hkl_[index][2]*hkl_[index][2]*k;
      single_grad[4] = result*
        tmp_const*hkl_[index][0]*hkl_[index][1]*2.0*k;
      single_grad[5] = result*
        tmp_const*hkl_[index][0]*hkl_[index][2]*2.0*k;
      single_grad[6] = result*
        tmp_const*hkl_[index][1]*hkl_[index][2]*2.0*k;
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

    scitbx::af::shared< FloatType >
    hessian(unsigned index)
    {
      // first we need to know the pre-multiplier
      FloatType tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8, k;
      FloatType sndder,frstder;
      k = get_k( index );
      tmp1=i_nat_[index]-k*k*i_der_[index];
      tmp2=sig_der_[index]*sig_der_[index];
      tmp3=k*k*k*k*tmp2+sig_nat_[index]*sig_nat_[index];
      //-------------------------------------------
      tmp4 = 32*k*k*k*k*k*k*tmp1*tmp1*tmp2*tmp2/(tmp3*tmp3*tmp3);
      tmp5 = 32*k*k*k*k*i_der_[index]*tmp1*tmp2/(tmp3*tmp3);
      tmp6 =-12*k*k*tmp1*tmp1*tmp2/(tmp3*tmp3);
      tmp7 = 8*i_der_[index]*i_der_[index]*k*k/tmp3;
      tmp8 =-4*i_der_[index]*tmp1/tmp3;
      sndder=tmp4+tmp5+tmp6+tmp7+tmp8;
      //-------------------------------------------
      tmp4 = -4.0*k*k*k*tmp1*tmp1*tmp2/(tmp3*tmp3);
      tmp5 = -4.0*i_der_[index]*k*tmp1/tmp3;
      frstder=tmp4+tmp5;

      scitbx::af::shared<FloatType> tmp_hes(
        7*(7+1)/2, scitbx::af::init_functor_null<FloatType>() );


      FloatType dp,d0,d1,d2,d3,d4,d5,tps;
      tps=scitbx::constants::pi*scitbx::constants::pi*scale_rwgk_;
      dp=-1;
      d0=-2*tps*hkl_[index][0]*hkl_[index][0];
      d1=-2*tps*hkl_[index][1]*hkl_[index][1];
      d2=-2*tps*hkl_[index][2]*hkl_[index][2];
      d3=-4*tps*hkl_[index][0]*hkl_[index][1];
      d4=-4*tps*hkl_[index][0]*hkl_[index][2];
      d5=-4*tps*hkl_[index][1]*hkl_[index][2];
      //-------------------------------------------

      FloatType* h = tmp_hes.begin();
      // scale factor related stuff
      *h++=frstder*k + sndder*dp*dp*k*k;//0,0
      *h++=frstder*k*d0 + sndder*d0*k*k;//0,1
      *h++=frstder*k*d1 + sndder*d1*k*k;//0,2
      *h++=frstder*k*d2 + sndder*d2*k*k;//0,3
      *h++=frstder*k*d3 + sndder*d3*k*k;//0,4
      *h++=frstder*k*d4 + sndder*d4*k*k;//0,5
      *h++=frstder*k*d5 + sndder*d5*k*k;//0,6

      //------------------------------------------
      // u0 --> u0, u1, u2, u3, u4, u5
      *h++=frstder*k*d0*d0 + sndder*k*k*d0*d0;//
      *h++=frstder*k*d0*d1 + sndder*k*k*d0*d1;//
      *h++=frstder*k*d0*d2 + sndder*k*k*d0*d2;//
      *h++=frstder*k*d0*d3 + sndder*k*k*d0*d3;//
      *h++=frstder*k*d0*d4 + sndder*k*k*d0*d4;//
      *h++=frstder*k*d0*d5 + sndder*k*k*d0*d5;//
      //------------------------------------------
      // u1 --> u1, u2, u3, u4, u5
      *h++=frstder*k*d1*d1 + sndder*k*k*d1*d1;//
      *h++=frstder*k*d1*d2 + sndder*k*k*d1*d2;//
      *h++=frstder*k*d1*d3 + sndder*k*k*d1*d3;//
      *h++=frstder*k*d1*d4 + sndder*k*k*d1*d4;//
      *h++=frstder*k*d1*d5 + sndder*k*k*d1*d5;//
      //------------------------------------------
      // u2 --> u2, u3, u4, u5
      *h++=frstder*k*d2*d2 + sndder*k*k*d2*d2;//
      *h++=frstder*k*d2*d3 + sndder*k*k*d2*d3;//
      *h++=frstder*k*d2*d4 + sndder*k*k*d2*d4;//
      *h++=frstder*k*d2*d5 + sndder*k*k*d2*d5;//
      //------------------------------------------
      // u3 --> u3, u4, u5
      *h++=frstder*k*d3*d3 + sndder*k*k*d3*d3;//
      *h++=frstder*k*d3*d4 + sndder*k*k*d3*d4;//
      *h++=frstder*k*d3*d5 + sndder*k*k*d3*d5;//
      //------------------------------------------
      // u4 --> u4, u5
      *h++=frstder*k*d4*d4 + sndder*k*k*d4*d4;//
      *h++=frstder*k*d4*d5 + sndder*k*k*d4*d5;//

      //------------------------------------------
      // u5 --> u5
      *h++=frstder*k*d5*d5 + sndder*k*k*d5*d5;//

      return (tmp_hes);

    }

    scitbx::af::shared<FloatType>
    hessian()
    {
      scitbx::af::shared<FloatType> hes_mat(
         7*(7+1)/2, 0 );

      scitbx::af::shared<FloatType> tmp_hes_mat(
         7*(7+1)/2, scitbx::af::init_functor_null<FloatType>() );

      // loop over all reflections and get second derivatives
      for (unsigned index=0;index<hkl_.size();index++){
        tmp_hes_mat = hessian( index );
        hes_mat+=tmp_hes_mat;
      }
      return(hes_mat);
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


  /* relative scaling via minimising
     Residual=sum (Fnat-k*Fder)**2/(wnat+k**2*wder)
  */

  template <typename FloatType=double>
  class least_squares_on_f_wt
  {
    public:
    /*! Default constructor */
    least_squares_on_f_wt() {};
    /*! Pass in all that is needed */
    least_squares_on_f_wt(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
      scitbx::af::const_ref< FloatType > const& f_nat,
      scitbx::af::const_ref< FloatType > const& sig_nat,
      scitbx::af::const_ref< FloatType > const& f_der,
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
      SCITBX_ASSERT ( hkl.size() ==  f_nat.size() );
      SCITBX_ASSERT ( hkl.size() ==  sig_nat.size() );
      SCITBX_ASSERT ( hkl.size() ==  f_der.size() );
      SCITBX_ASSERT ( hkl.size() ==  sig_der.size() );

      // fill up the arrays please
      for (unsigned ii=0;ii<hkl.size();ii++){
        hkl_.push_back( hkl[ii] );
        f_nat_.push_back( f_nat[ii] );
        sig_nat_.push_back( sig_nat[ii] );
        f_der_.push_back( f_der[ii] );
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

    scitbx::af::shared<FloatType>
    get_gradient(unsigned index)
    {
      return( gradient(index) );
    }

    FloatType
    get_function(unsigned index)
    {
      return( function(index) );
    }

    scitbx::af::shared<FloatType>
    hessian_as_packed_u(unsigned index)
    {
      return( hessian(index) );
    }

    scitbx::af::shared<FloatType>
    hessian_as_packed_u()
    {
      return( hessian() );
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
      if(result > 40.) result = std::exp(40.);
      else result = std::exp(result);
      return( result );
    }

    FloatType
    function( unsigned index )
    {
      FloatType result=0, k, weight;
      k=get_k( index );
      result = f_nat_[index] - f_der_[index]*k;
      result=result*result;
      weight = sig_nat_[index]*sig_nat_[index];
      weight = weight+sig_der_[index]*sig_der_[index]*k*k;
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

      tmp3 = -k*f_der_[index] + f_nat_[index];
      tmp4 = k*k*sig_der_[index]*sig_der_[index]
                   + sig_nat_[index]*sig_nat_[index];

      tmp1 = -2.0*k*tmp3*tmp3*sig_der_[index]*sig_der_[index]/(tmp4*tmp4);
      tmp2 = -2.0*f_der_[index]*tmp3/tmp4;

      result = tmp1+tmp2;

      FloatType tmp_const=2.0*scitbx::constants::pi*scitbx::constants::pi*
        scale_rwgk_;

      single_grad[0] = result*
        -k;
      single_grad[1] = result*
        tmp_const*hkl_[index][0]*hkl_[index][0]*k;
      single_grad[2] = result*
        tmp_const*hkl_[index][1]*hkl_[index][1]*k;
      single_grad[3] = result*
        tmp_const*hkl_[index][2]*hkl_[index][2]*k;
      single_grad[4] = result*
        tmp_const*hkl_[index][0]*hkl_[index][1]*2.0*k;
      single_grad[5] = result*
        tmp_const*hkl_[index][0]*hkl_[index][2]*2.0*k;
      single_grad[6] = result*
        tmp_const*hkl_[index][1]*hkl_[index][2]*2.0*k;

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




    scitbx::af::shared< FloatType >
    hessian(unsigned index)
    {
      // first we need to know the pre-multiplier
      FloatType tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7, k;
      FloatType sndder,frstder;
      k = get_k( index );
      tmp1=f_nat_[index]-k*f_der_[index];
      tmp2=sig_der_[index]*sig_der_[index];
      tmp3=k*k*tmp2+sig_nat_[index]*sig_nat_[index];
      //------sndder-------------------------------------
      tmp4 = 8*k*k*tmp1*tmp1*tmp2*tmp2/(tmp3*tmp3*tmp3);
      tmp5 = 8*k*f_der_[index]*tmp1*tmp2/(tmp3*tmp3);
      tmp6 =-2*tmp1*tmp1*tmp2/(tmp3*tmp3);
      tmp7 = 2*f_der_[index]*f_der_[index]/tmp3;
      sndder=tmp4+tmp5+tmp6+tmp7;
      //------frsder-------------------------------------
      tmp4 = -2.0*k*tmp1*tmp1*tmp2/(tmp3*tmp3);
      tmp5 = -2.0*f_der_[index]*tmp1/tmp3;
      frstder=tmp4+tmp5;

      scitbx::af::shared<FloatType> tmp_hes(
        7*(7+1)/2, scitbx::af::init_functor_null<FloatType>() );


      FloatType dp,d0,d1,d2,d3,d4,d5,tps;
      tps=scitbx::constants::pi*scitbx::constants::pi*scale_rwgk_;
      dp=-1;
      d0=-2*tps*hkl_[index][0]*hkl_[index][0];
      d1=-2*tps*hkl_[index][1]*hkl_[index][1];
      d2=-2*tps*hkl_[index][2]*hkl_[index][2];
      d3=-4*tps*hkl_[index][0]*hkl_[index][1];
      d4=-4*tps*hkl_[index][0]*hkl_[index][2];
      d5=-4*tps*hkl_[index][1]*hkl_[index][2];
      //-------------------------------------------

      FloatType* h = tmp_hes.begin();
      // scale factor related stuff
      *h++=frstder*k + sndder*dp*dp*k*k;//0,0
      *h++=frstder*k*d0 + sndder*d0*k*k;//0,1
      *h++=frstder*k*d1 + sndder*d1*k*k;//0,2
      *h++=frstder*k*d2 + sndder*d2*k*k;//0,3
      *h++=frstder*k*d3 + sndder*d3*k*k;//0,4
      *h++=frstder*k*d4 + sndder*d4*k*k;//0,5
      *h++=frstder*k*d5 + sndder*d5*k*k;//0,6

      //------------------------------------------
      // u0 --> u0, u1, u2, u3, u4, u5
      *h++=frstder*k*d0*d0 + sndder*k*k*d0*d0;//
      *h++=frstder*k*d0*d1 + sndder*k*k*d0*d1;//
      *h++=frstder*k*d0*d2 + sndder*k*k*d0*d2;//
      *h++=frstder*k*d0*d3 + sndder*k*k*d0*d3;//
      *h++=frstder*k*d0*d4 + sndder*k*k*d0*d4;//
      *h++=frstder*k*d0*d5 + sndder*k*k*d0*d5;//
      //------------------------------------------
      // u1 --> u1, u2, u3, u4, u5
      *h++=frstder*k*d1*d1 + sndder*k*k*d1*d1;//
      *h++=frstder*k*d1*d2 + sndder*k*k*d1*d2;//
      *h++=frstder*k*d1*d3 + sndder*k*k*d1*d3;//
      *h++=frstder*k*d1*d4 + sndder*k*k*d1*d4;//
      *h++=frstder*k*d1*d5 + sndder*k*k*d1*d5;//
      //------------------------------------------
      // u2 --> u2, u3, u4, u5
      *h++=frstder*k*d2*d2 + sndder*k*k*d2*d2;//
      *h++=frstder*k*d2*d3 + sndder*k*k*d2*d3;//
      *h++=frstder*k*d2*d4 + sndder*k*k*d2*d4;//
      *h++=frstder*k*d2*d5 + sndder*k*k*d2*d5;//
      //------------------------------------------
      // u3 --> u3, u4, u5
      *h++=frstder*k*d3*d3 + sndder*k*k*d3*d3;//
      *h++=frstder*k*d3*d4 + sndder*k*k*d3*d4;//
      *h++=frstder*k*d3*d5 + sndder*k*k*d3*d5;//
      //------------------------------------------
      // u4 --> u4, u5
      *h++=frstder*k*d4*d4 + sndder*k*k*d4*d4;//
      *h++=frstder*k*d4*d5 + sndder*k*k*d4*d5;//

      //------------------------------------------
      // u5 --> u5
      *h++=frstder*k*d5*d5 + sndder*k*k*d5*d5;//

      return (tmp_hes);

    }

    scitbx::af::shared<FloatType>
    hessian()
    {
      scitbx::af::shared<FloatType> hes_mat(
         7*(7+1)/2, 0 );

      scitbx::af::shared<FloatType> tmp_hes_mat(
         7*(7+1)/2, scitbx::af::init_functor_null<FloatType>() );

      // loop over all reflections and get second derivatives
      for (unsigned index=0;index<hkl_.size();index++){
        tmp_hes_mat = hessian( index );
        hes_mat+=tmp_hes_mat;
      }
      return(hes_mat);
    }








    scitbx::af::shared< cctbx::miller::index<> > hkl_;
    scitbx::af::shared< FloatType > f_nat_;
    scitbx::af::shared< FloatType > sig_nat_;
    scitbx::af::shared< FloatType > f_der_;
    scitbx::af::shared< FloatType > sig_der_;
    FloatType p_scale_;
    cctbx::uctbx::unit_cell unit_cell_;
    scitbx::sym_mat3<FloatType> u_rwgk_;
    FloatType scale_rwgk_;

  };














  /* relative scaling via minimising
     Residual=sum (Inat-k*k*Ider)**2/(wcombined)
  */
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


    scitbx::af::shared<FloatType>
    get_gradient(unsigned index)
    {
      return( gradient(index) );
    }

    FloatType
    get_function(unsigned index)
    {
      return( function(index) );
    }

    scitbx::af::shared<FloatType>
    hessian_as_packed_u(unsigned index)
    {
      return( hessian(index) );
    }

    scitbx::af::shared<FloatType>
    hessian_as_packed_u()
    {
      return( hessian() );
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
      if(result > 40.) result = std::exp(40.);
      else result = std::exp(result);
      return( result );
    }

    FloatType
    function( unsigned index )
    {
      FloatType result=0, k, weight;
      k=get_k( index );
      weight = (i_nat_[index]/i_der_[index])*sig_der_[index];
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
        tmp_const*hkl_[index][0]*hkl_[index][0]*k;
      single_grad[2] = 2.0*result*weight* -2.0*k*i_der_[index]*
        tmp_const*hkl_[index][1]*hkl_[index][1]*k;
      single_grad[3] = 2.0*result*weight* -2.0*k*i_der_[index]*
        tmp_const*hkl_[index][2]*hkl_[index][2]*k;
      single_grad[4] = 2.0*result*weight* -2.0*k*i_der_[index]*
        tmp_const*hkl_[index][0]*hkl_[index][1]*2.0*k;
      single_grad[5] = 2.0*result*weight* -2.0*k*i_der_[index]*
        tmp_const*hkl_[index][0]*hkl_[index][2]*2.0*k;
      single_grad[6] = 2.0*result*weight* -2.0*k*i_der_[index]*
        tmp_const*hkl_[index][1]*hkl_[index][2]*2.0*k;

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


    scitbx::af::shared< FloatType >
    hessian(unsigned index)
    {
      // first we need to know the pre-multipliers
      FloatType tmp1,tmp2,tmp3,tmp4, k;
      FloatType sndder,frstder;
      k = get_k( index );
      tmp1=i_nat_[index]-k*k*i_der_[index];

      tmp2 = i_nat_[index]*sig_der_[index]/i_der_[index];
      tmp2 = tmp2*tmp2 + sig_nat_[index]*sig_nat_[index];

      //------sndder-------------------------------------
      tmp3 = 8*i_der_[index]*i_der_[index]*k*k/tmp2;
      tmp4 = -4*i_der_[index]*tmp1/tmp2;
      sndder=tmp3+tmp4;
      //------frsder-------------------------------------
      frstder=-4*i_der_[index]*k*tmp1/tmp2;

      scitbx::af::shared<FloatType> tmp_hes(
        7*(7+1)/2, scitbx::af::init_functor_null<FloatType>() );


      FloatType dp,d0,d1,d2,d3,d4,d5,tps;
      tps=scitbx::constants::pi*scitbx::constants::pi*scale_rwgk_;
      dp=-1;
      d0=-2*tps*hkl_[index][0]*hkl_[index][0];
      d1=-2*tps*hkl_[index][1]*hkl_[index][1];
      d2=-2*tps*hkl_[index][2]*hkl_[index][2];
      d3=-4*tps*hkl_[index][0]*hkl_[index][1];
      d4=-4*tps*hkl_[index][0]*hkl_[index][2];
      d5=-4*tps*hkl_[index][1]*hkl_[index][2];
      //-------------------------------------------

      FloatType* h = tmp_hes.begin();
      // scale factor related stuff
      *h++=frstder*k + sndder*dp*dp*k*k;//0,0
      *h++=frstder*k*d0 + sndder*d0*k*k;//0,1
      *h++=frstder*k*d1 + sndder*d1*k*k;//0,2
      *h++=frstder*k*d2 + sndder*d2*k*k;//0,3
      *h++=frstder*k*d3 + sndder*d3*k*k;//0,4
      *h++=frstder*k*d4 + sndder*d4*k*k;//0,5
      *h++=frstder*k*d5 + sndder*d5*k*k;//0,6

      //------------------------------------------
      // u0 --> u0, u1, u2, u3, u4, u5
      *h++=frstder*k*d0*d0 + sndder*k*k*d0*d0;//
      *h++=frstder*k*d0*d1 + sndder*k*k*d0*d1;//
      *h++=frstder*k*d0*d2 + sndder*k*k*d0*d2;//
      *h++=frstder*k*d0*d3 + sndder*k*k*d0*d3;//
      *h++=frstder*k*d0*d4 + sndder*k*k*d0*d4;//
      *h++=frstder*k*d0*d5 + sndder*k*k*d0*d5;//
      //------------------------------------------
      // u1 --> u1, u2, u3, u4, u5
      *h++=frstder*k*d1*d1 + sndder*k*k*d1*d1;//
      *h++=frstder*k*d1*d2 + sndder*k*k*d1*d2;//
      *h++=frstder*k*d1*d3 + sndder*k*k*d1*d3;//
      *h++=frstder*k*d1*d4 + sndder*k*k*d1*d4;//
      *h++=frstder*k*d1*d5 + sndder*k*k*d1*d5;//
      //------------------------------------------
      // u2 --> u2, u3, u4, u5
      *h++=frstder*k*d2*d2 + sndder*k*k*d2*d2;//
      *h++=frstder*k*d2*d3 + sndder*k*k*d2*d3;//
      *h++=frstder*k*d2*d4 + sndder*k*k*d2*d4;//
      *h++=frstder*k*d2*d5 + sndder*k*k*d2*d5;//
      //------------------------------------------
      // u3 --> u3, u4, u5
      *h++=frstder*k*d3*d3 + sndder*k*k*d3*d3;//
      *h++=frstder*k*d3*d4 + sndder*k*k*d3*d4;//
      *h++=frstder*k*d3*d5 + sndder*k*k*d3*d5;//
      //------------------------------------------
      // u4 --> u4, u5
      *h++=frstder*k*d4*d4 + sndder*k*k*d4*d4;//
      *h++=frstder*k*d4*d5 + sndder*k*k*d4*d5;//

      //------------------------------------------
      // u5 --> u5
      *h++=frstder*k*d5*d5 + sndder*k*k*d5*d5;//

      return (tmp_hes);

    }

    scitbx::af::shared<FloatType>
    hessian()
    {
      scitbx::af::shared<FloatType> hes_mat(
         7*(7+1)/2, 0 );

      scitbx::af::shared<FloatType> tmp_hes_mat(
         7*(7+1)/2, scitbx::af::init_functor_null<FloatType>() );

      // loop over all reflections and get second derivatives
      for (unsigned index=0;index<hkl_.size();index++){
        tmp_hes_mat = hessian( index );
        hes_mat+=tmp_hes_mat;
      }
      return(hes_mat);
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





  /* relative scaling via minimising
     Residual=sum (Fnat-k*Fder)**2/(wcombined)
  */
  template <typename FloatType=double>
  class least_squares_on_f
  {
    public:
    /*! Default constructor */
    least_squares_on_f() {};
    /*! Pass in all that is needed */
    least_squares_on_f(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
      scitbx::af::const_ref< FloatType > const& f_nat,
      scitbx::af::const_ref< FloatType > const& sig_nat,
      scitbx::af::const_ref< FloatType > const& f_der,
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
      SCITBX_ASSERT ( hkl.size() ==  f_nat.size() );
      SCITBX_ASSERT ( hkl.size() ==  sig_nat.size() );
      SCITBX_ASSERT ( hkl.size() ==  f_der.size() );
      SCITBX_ASSERT ( hkl.size() ==  sig_der.size() );

      // fill up the arrays please
      for (unsigned ii=0;ii<hkl.size();ii++){
        hkl_.push_back( hkl[ii] );
        f_nat_.push_back( f_nat[ii] );
        sig_nat_.push_back( sig_nat[ii] );
        f_der_.push_back( f_der[ii] );
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


    scitbx::af::shared<FloatType>
    get_gradient(unsigned index)
    {
      return( gradient(index) );
    }

    FloatType
    get_function(unsigned index)
    {
      return( function(index) );
    }

    scitbx::af::shared<FloatType>
    hessian_as_packed_u(unsigned index)
    {
      return( hessian(index) );
    }

    scitbx::af::shared<FloatType>
    hessian_as_packed_u()
    {
      return( hessian() );
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
      if(result > 40.) result = std::exp(40.);
      else result = std::exp(result);
      return( result );
    }

    FloatType
    function( unsigned index )
    {
      FloatType result=0, k, weight, denom, nom;
      k=get_k( index );
      denom = f_der_[index];
      nom = f_nat_[index]*sig_der_[index];
      if (denom <= 0 ){
        denom = 1.0;
        nom = sig_der_[index];
      }
      weight = nom/denom;
      weight = weight*weight + sig_nat_[index]*sig_nat_[index];
      result = f_nat_[index] - f_der_[index]*k;
      if(std::abs(result)<1.e+50) {
        result = result*result;
      }
      if (weight > 0 ){
        result = result/weight;
      }
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
      FloatType result=0, k, weight=0., denom, nom;
      k=get_k( index );
      nom = f_nat_[index]*sig_der_[index];
      denom = f_der_[index];
      if ( denom <= 0){
        denom = 1.0;
        nom = sig_der_[index];
      }
      if(denom != 0) weight = nom/denom;
      weight = weight*weight + sig_nat_[index]*sig_nat_[index];
      if (weight >=0){
      weight = 1.0/weight;
      } else {
        weight = 0.0;
      }
      result = f_nat_[index] - f_der_[index]*k;
      FloatType tmp_const=2.0*scitbx::constants::pi*scitbx::constants::pi*
        scale_rwgk_;
      cctbx::miller::index<int> ind = hkl_[index];
      FloatType rw2 = -2.*result*weight*k*f_der_[index];
      FloatType rw2t = rw2*tmp_const;
      single_grad[0] = -1.*rw2;
      single_grad[1] = rw2t*ind[0]*ind[0];
      single_grad[2] = rw2t*ind[1]*ind[1];
      single_grad[3] = rw2t*ind[2]*ind[2];
      single_grad[4] = rw2t*ind[0]*ind[1]*2.0;
      single_grad[5] = rw2t*ind[0]*ind[2]*2.0;
      single_grad[6] = rw2t*ind[1]*ind[2]*2.0;
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



        scitbx::af::shared< FloatType >
    hessian(unsigned index)
    {
      // first we need to know the pre-multipliers
      FloatType tmp1,tmp2,k, nom, denom;
      FloatType sndder,frstder;
      k = get_k( index );
      tmp1=f_nat_[index]-k*f_der_[index];

      nom = f_nat_[index]*sig_der_[index];
      denom = f_der_[index];
      if (denom <=0){
        nom = sig_der_[index];
        denom = 1.0;
      }
      tmp2 = nom/denom;
      tmp2 = tmp2*tmp2 + sig_nat_[index]*sig_nat_[index];
      if (tmp2 <=0){
        tmp2 = 1.0;
      }

      //------sndder-------------------------------------
      sndder=2*f_der_[index]*f_der_[index]/tmp2;
      //------frsder-------------------------------------
      frstder=-2*f_der_[index]*tmp1/tmp2;

      scitbx::af::shared<FloatType> tmp_hes(
        7*(7+1)/2, scitbx::af::init_functor_null<FloatType>() );


      FloatType dp,d0,d1,d2,d3,d4,d5,tps;
      tps=scitbx::constants::pi*scitbx::constants::pi*scale_rwgk_;
      dp=-1;
      d0=-2*tps*hkl_[index][0]*hkl_[index][0];
      d1=-2*tps*hkl_[index][1]*hkl_[index][1];
      d2=-2*tps*hkl_[index][2]*hkl_[index][2];
      d3=-4*tps*hkl_[index][0]*hkl_[index][1];
      d4=-4*tps*hkl_[index][0]*hkl_[index][2];
      d5=-4*tps*hkl_[index][1]*hkl_[index][2];
      //-------------------------------------------

      FloatType* h = tmp_hes.begin();
      // scale factor related stuff
      *h++=frstder*k + sndder*dp*dp*k*k;//0,0
      *h++=frstder*k*d0 + sndder*d0*k*k;//0,1
      *h++=frstder*k*d1 + sndder*d1*k*k;//0,2
      *h++=frstder*k*d2 + sndder*d2*k*k;//0,3
      *h++=frstder*k*d3 + sndder*d3*k*k;//0,4
      *h++=frstder*k*d4 + sndder*d4*k*k;//0,5
      *h++=frstder*k*d5 + sndder*d5*k*k;//0,6

      //------------------------------------------
      // u0 --> u0, u1, u2, u3, u4, u5
      *h++=frstder*k*d0*d0 + sndder*k*k*d0*d0;//
      *h++=frstder*k*d0*d1 + sndder*k*k*d0*d1;//
      *h++=frstder*k*d0*d2 + sndder*k*k*d0*d2;//
      *h++=frstder*k*d0*d3 + sndder*k*k*d0*d3;//
      *h++=frstder*k*d0*d4 + sndder*k*k*d0*d4;//
      *h++=frstder*k*d0*d5 + sndder*k*k*d0*d5;//
      //------------------------------------------
      // u1 --> u1, u2, u3, u4, u5
      *h++=frstder*k*d1*d1 + sndder*k*k*d1*d1;//
      *h++=frstder*k*d1*d2 + sndder*k*k*d1*d2;//
      *h++=frstder*k*d1*d3 + sndder*k*k*d1*d3;//
      *h++=frstder*k*d1*d4 + sndder*k*k*d1*d4;//
      *h++=frstder*k*d1*d5 + sndder*k*k*d1*d5;//
      //------------------------------------------
      // u2 --> u2, u3, u4, u5
      *h++=frstder*k*d2*d2 + sndder*k*k*d2*d2;//
      *h++=frstder*k*d2*d3 + sndder*k*k*d2*d3;//
      *h++=frstder*k*d2*d4 + sndder*k*k*d2*d4;//
      *h++=frstder*k*d2*d5 + sndder*k*k*d2*d5;//
      //------------------------------------------
      // u3 --> u3, u4, u5
      *h++=frstder*k*d3*d3 + sndder*k*k*d3*d3;//
      *h++=frstder*k*d3*d4 + sndder*k*k*d3*d4;//
      *h++=frstder*k*d3*d5 + sndder*k*k*d3*d5;//
      //------------------------------------------
      // u4 --> u4, u5
      *h++=frstder*k*d4*d4 + sndder*k*k*d4*d4;//
      *h++=frstder*k*d4*d5 + sndder*k*k*d4*d5;//

      //------------------------------------------
      // u5 --> u5
      *h++=frstder*k*d5*d5 + sndder*k*k*d5*d5;//

      return (tmp_hes);

    }

    scitbx::af::shared<FloatType>
    hessian()
    {
      scitbx::af::shared<FloatType> hes_mat(
         7*(7+1)/2, 0 );

      scitbx::af::shared<FloatType> tmp_hes_mat(
         7*(7+1)/2, scitbx::af::init_functor_null<FloatType>() );

      // loop over all reflections and get second derivatives
      for (unsigned index=0;index<hkl_.size();index++){
        tmp_hes_mat = hessian( index );
        hes_mat+=tmp_hes_mat;
      }
      return(hes_mat);
    }


    scitbx::af::shared< cctbx::miller::index<> > hkl_;
    scitbx::af::shared< FloatType > f_nat_;
    scitbx::af::shared< FloatType > sig_nat_;
    scitbx::af::shared< FloatType > f_der_;
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
  //
  // Note that both routines can be used for
  // intensity or amplitude based local scaling

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
                     at_least_this_number_of_neighbours ),
    location_of_master_in_set_( hkl_sets, space_group, anomalous_flag ),
    min_neighbours_( hkl_sets.size() ),
    max_neighbours_( 0 ),
    mean_neighbours_( 0 ),
    local_scales_( hkl_sets.size(), 1),
    local_scales_sigmas_( hkl_sets.size(), 0.01),
    use_weights_( weights )
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

      scitbx::af::shared<long> locations;
      locations = location_of_master_in_set_.find_hkl( hkl_master_.const_ref() );

      for (unsigned ii=0;ii<hkl_sets_.size();ii++){

        if (area_generator_.area_[property_generator_.set_lut_[ ii ] ].size() < min_neighbours_ ){
          min_neighbours_ = area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
        }
        if (area_generator_.area_[property_generator_.set_lut_[ ii ] ].size() > max_neighbours_ ){
          max_neighbours_ = area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
        }
        mean_neighbours_ += area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();

        FloatType nat_mean=0, der_mean=0;
        FloatType nat_var=0, der_var=0;
        // compute the mean for the native and derivativesummats

        for (unsigned jj=1;
             jj<area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
             jj++){

          nb_index = area_generator_.area_[
            property_generator_.set_lut_[ ii ] ][jj];
          nb_index = locations[ nb_index ];

          SCITBX_ASSERT(nb_index<hkl_sets_.size());

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
      mean_neighbours_/=hkl_sets_.size();
    }

    scitbx::af::shared<FloatType>
    get_scales()
    {
      return( local_scales_ );
    }

    scitbx::af::tiny<FloatType,3>
    stats()
    {
      scitbx::af::tiny<FloatType,3> result;
      result[0]=min_neighbours_;
      result[1]=max_neighbours_;
      result[2]=mean_neighbours_;
      return(result);
    }


    protected:
    property_matching_indices<FloatType> property_generator_;
    cctbx::miller::lookup_utils::local_area<FloatType> area_generator_;
    cctbx::miller::lookup_utils::lookup_tensor<FloatType> location_of_master_in_set_;


    FloatType min_neighbours_, max_neighbours_, mean_neighbours_;

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
                     at_least_this_number_of_neighbours ),
    location_of_master_in_set_( hkl_sets, space_group, anomalous_flag ),
    min_neighbours_( hkl_sets.size() ),
    max_neighbours_( 0 ),
    mean_neighbours_( 0 ),
    local_scales_( hkl_sets.size(), 1),
    local_scales_sigmas_( hkl_sets.size(), 0.01),
    use_weights_( weights )
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

      scitbx::af::shared<long> locations;
      locations = location_of_master_in_set_.find_hkl( hkl_master_.const_ref() );

      for (unsigned ii=0;ii<hkl_sets_.size();ii++){


        if (area_generator_.area_[property_generator_.set_lut_[ ii ] ].size() < min_neighbours_ ){
          min_neighbours_ = area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
        }
        if (area_generator_.area_[property_generator_.set_lut_[ ii ] ].size() > max_neighbours_ ){
          max_neighbours_ = area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
        }
        mean_neighbours_ += area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();


        FloatType top=0, bottom=0;
        for (unsigned jj=1;
             jj<area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
             jj++){

          nb_index = area_generator_.area_[
            property_generator_.set_lut_[ ii ] ][jj];
          nb_index = locations[ nb_index ];

          SCITBX_ASSERT(nb_index<hkl_sets_.size());

          weight = sigma_set_b_[ nb_index ]* data_set_a_[ nb_index  ]/
            data_set_b_[ nb_index  ];
          weight = weight*weight
           +  sigma_set_a_[ nb_index ]*sigma_set_a_[ nb_index ];

          weight=weight*multiplier+constant;
          weight=1.0/weight;

          top+=weight*data_set_a_[nb_index]*data_set_b_[nb_index];
          bottom+=weight*data_set_b_[nb_index]*data_set_b_[nb_index];

        }

        if (std::fabs(bottom)>1e-9){
          local_scales_[ii]=top/bottom;
        }

      }
      mean_neighbours_ /= hkl_sets_.size();
    }

    scitbx::af::shared<FloatType>
    get_scales()
    {
      return( local_scales_ );
    }

    scitbx::af::tiny<FloatType,3>
    stats()
    {
      scitbx::af::tiny<FloatType,3> result;
      result[0]=min_neighbours_;
      result[1]=max_neighbours_;
      result[2]=mean_neighbours_;
      return(result);
    }



    protected:
    property_matching_indices<FloatType> property_generator_;
    cctbx::miller::lookup_utils::local_area<FloatType> area_generator_;
    cctbx::miller::lookup_utils::lookup_tensor<FloatType> location_of_master_in_set_;

    FloatType min_neighbours_, max_neighbours_, mean_neighbours_;

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
  /*! Local scaling where scale factors are determined
   *  using an algorithm devised by Nikonov,
   *  Acta Cryst A39, 693-697.
   */

  template <typename FloatType=double>
  class local_scaling_nikonov
  {
    public:
    /*! default constructor */
    local_scaling_nikonov(){}
    /*! This does (allmost) everything */
    local_scaling_nikonov(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_master,
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_sets,

      scitbx::af::const_ref< FloatType > const& data_set_a,
      scitbx::af::const_ref< FloatType > const& data_set_b,
      scitbx::af::const_ref< FloatType > const& epsilons,
      scitbx::af::const_ref< bool > const& centric,
      FloatType const& threshold,

      cctbx::sgtbx::space_group const& space_group,
      bool const& anomalous_flag,
      std::size_t const& radius,
      std::size_t const& depth,
      std::size_t const& at_least_this_number_of_neighbours)
    :
    threshold_(threshold),
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
                     at_least_this_number_of_neighbours ),
    location_of_master_in_set_( hkl_sets, space_group, anomalous_flag ),
    min_neighbours_( hkl_sets.size() ),
    max_neighbours_( 0 ),
    mean_neighbours_( 0 ),
    local_scales_( hkl_sets.size(), 1),
    local_scales_sigmas_( hkl_sets.size(), 0.01)
    {

      for (unsigned ii=0;ii<hkl_master.size();ii++){
        hkl_master_.push_back( hkl_master[ii] );
      }

      for (unsigned ii=0;ii<hkl_sets.size();ii++){
        hkl_sets_.push_back( hkl_sets[ii] );
        data_set_a_.push_back( data_set_a[ii] );
        data_set_b_.push_back( data_set_b[ii] );
        epsilons_.push_back( epsilons[ii] );
        centric_.push_back( centric[ii] );
      }

      // now eveything is set, please do the local scaling
      scale_it();
    }

    void
    scale_it()
    {
      // loop over all reflections of the native
      unsigned nb_index ;

      scitbx::af::shared<long> locations;
      locations = location_of_master_in_set_.find_hkl( hkl_master_.const_ref() );

      for (unsigned ii=0;ii<hkl_sets_.size();ii++){

        // get some statistics regarding number of neighbours
        if (area_generator_.area_[property_generator_.set_lut_[ ii ] ].size() < min_neighbours_ ){
          min_neighbours_ = area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
        }
        if (area_generator_.area_[property_generator_.set_lut_[ ii ] ].size() > max_neighbours_ ){
          max_neighbours_ = area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
        }
        mean_neighbours_ += area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();


        FloatType mean_intensity_native=0;
        FloatType count=0;

        FloatType part1_ac=0, part1_c=0, part2_ac=0, n_ac=0, n_c=0, a,b, k_ac,k_c;

        // loop over all neighbours
        for (unsigned jj=1;
             jj<area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
             jj++){

          nb_index = area_generator_.area_[
            property_generator_.set_lut_[ ii ] ][jj];
          nb_index = locations[ nb_index ];

          mean_intensity_native+=data_set_a_[ nb_index  ]/epsilons_[nb_index];
          count++;
        }
        mean_intensity_native/=count;

        // the exact threshold is now known
        n_ac=0;
        n_c=0;
        for (unsigned jj=1;
             jj<area_generator_.area_[property_generator_.set_lut_[ ii ] ].size();
             jj++){

          nb_index = area_generator_.area_[
            property_generator_.set_lut_[ ii ] ][jj];
          nb_index = locations[ nb_index ];

          if (std::sqrt(mean_intensity_native)*threshold_*epsilons_[nb_index] <= data_set_a_[nb_index]  ){
            a=data_set_a_[nb_index]/epsilons_[nb_index];
            b=data_set_b_[nb_index]/epsilons_[nb_index];
            if (centric_[nb_index]){ // centric
              part1_c += b/a;
              n_c++;
            }
            else{
              part1_ac += b/a;
              part2_ac += (b*b)/(a*a);
              n_ac++;
            }
          }
        }
        k_c = 0.0;
        k_ac = 0.0;
        if (n_c > 1e-9){
          k_c = n_c/part1_c;
        }
        if (n_ac > 1e-9){
          k_ac = 2.0*part1_ac/n_ac - std::sqrt( part2_ac/n_ac );
          k_ac = 1.0/k_ac;
        }
        if (n_c+n_ac>1e-9){
          local_scales_[ii] = ( n_c*k_c + n_ac*k_ac)/(n_c+n_ac);
        }
        else{
          local_scales_[ii] = 1.0;
        }
      }
      mean_neighbours_ /= hkl_sets_.size();
    }

    scitbx::af::shared<FloatType>
    get_scales()
    {
      return( local_scales_ );
    }

    scitbx::af::tiny<FloatType,3>
    stats()
    {
      scitbx::af::tiny<FloatType,3> result;
      result[0]=min_neighbours_;
      result[1]=max_neighbours_;
      result[2]=mean_neighbours_;
      return(result);
    }



    protected:
    FloatType threshold_;
    property_matching_indices<FloatType> property_generator_;
    cctbx::miller::lookup_utils::local_area<FloatType> area_generator_;
    cctbx::miller::lookup_utils::lookup_tensor<FloatType> location_of_master_in_set_;

    FloatType min_neighbours_, max_neighbours_, mean_neighbours_;

    scitbx::af::shared< cctbx::miller::index<> > hkl_master_;
    scitbx::af::shared< cctbx::miller::index<> > hkl_sets_;

    scitbx::af::shared< FloatType > data_set_a_;
    scitbx::af::shared< FloatType > data_set_b_;
    scitbx::af::shared< FloatType > epsilons_;
    scitbx::af::shared< bool > centric_;

    scitbx::af::shared< FloatType > local_scales_;
    scitbx::af::shared< FloatType > local_scales_sigmas_;

  };

}}} // namespace mmtbx::scaling::relative_scaling
#endif // MMTBX_SCALING_RELATIVE_SCALING_H
