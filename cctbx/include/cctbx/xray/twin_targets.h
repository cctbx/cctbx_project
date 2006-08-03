#ifndef CCTBX_XRAY_TWIN_TARGETS_H
#define CCTBX_XRAY_TWIN_TARGETS_H

#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <complex>
#include <cmath>
#include <scitbx/math/bessel.h>
#include <cctbx/hendrickson_lattman.h>
#include <cctbx/miller/lookup_utils.h>
#include <cctbx/uctbx.h>
#include <cctbx/xray/f_model.h>
#include <scitbx/math/halton.h>

namespace cctbx { namespace xray { namespace twin_targets {

  template<typename FloatType>
  inline cctbx::miller::index<> twin_mate( cctbx::miller::index<> hkl,
                                           scitbx::mat3<FloatType> twin_law )
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


  template<typename FloatType> class least_squares_hemihedral_twinning_on_i{
  public:
  // You want to use this constructor
    least_squares_hemihedral_twinning_on_i(
      scitbx::af::const_ref< cctbx::miller::index<> >  const& hkl_obs,       //1 indices for calculated data
      scitbx::af::const_ref< FloatType >               const& i_obs,         //2 f calc
      scitbx::af::const_ref< FloatType >               const& w_obs,         //3 f bulk solvent
      scitbx::af::const_ref< cctbx::miller::index<> >  const& hkl_calc,      //4 f_model; not const to avoid CV issues
      sgtbx::space_group                               const& space_group,   //5 space group
      bool                                             const& anomalous_flag,//6 anomalous_flag
      FloatType                                        const& alpha,         //7 twin fraction
      scitbx::mat3<FloatType>                          const& twin_law       //8 twin law
      ):
      space_group_( space_group ),
      twin_law_(twin_law),
      alpha_(alpha)
      {
        CCTBX_ASSERT( (alpha >=0) && (alpha<=0.50) );
        CCTBX_ASSERT( hkl_obs.size() > 0);
        CCTBX_ASSERT( hkl_obs.size() == i_obs.size() );
        CCTBX_ASSERT( (hkl_obs.size() == w_obs.size()) || (w_obs.size()==0) );

        cctbx::miller::lookup_utils::lookup_tensor<FloatType>
          tmp_lookup_object( hkl_calc, space_group, anomalous_flag  );

        int tmp_loc;
        for (std::size_t ii=0;ii<hkl_obs.size();ii++){
          i_obs_.push_back( i_obs[ii] );
          if (w_obs.size() > 0){
            w_obs_.push_back( w_obs[ii] );
          }
          else {
            w_obs_.push_back( 1.0 );
          }
          tmp_loc = tmp_lookup_object.find_hkl( hkl_obs[ii] );
          CCTBX_ASSERT( tmp_loc >= 0 );
          calc_ori_lookup_table_.push_back( tmp_loc );
          tmp_loc = tmp_lookup_object.find_hkl( twin_mate( hkl_obs[ii],twin_law ) );
          CCTBX_ASSERT( tmp_loc >= 0 );
          calc_twin_lookup_table_ .push_back( tmp_loc );
        }
      }


      FloatType target(scitbx::af::const_ref<std::complex<FloatType> >
                       const& f_model) const
      {
        FloatType result=0,aa,ba,ab,bb,obs,calc;
        long calc_index_a, calc_index_b;
        for (std::size_t ii=0;ii<i_obs_.size();ii++){
          calc_index_a = calc_ori_lookup_table_[ ii ];
          calc_index_b = calc_twin_lookup_table_[ ii ];
          aa = f_model[calc_index_a].real();
          ba = f_model[calc_index_a].imag();
          ab = f_model[calc_index_b].real();
          bb = f_model[calc_index_b].imag();
          calc = (1-alpha_)*(aa*aa + ba*ba) + alpha_*(ab*ab + bb*bb);
          obs = i_obs_[ii];
          //std::cout << ii << " " << calc << " " << obs <<  " " << std::endl;
          result += w_obs_[ii]*(obs-calc)*(obs-calc);
        }
        return( result );
      }


      scitbx::af::tiny<scitbx::af::shared<FloatType>, 2 > d_target_d_ab
      (scitbx::af::const_ref<std::complex<FloatType> > const& f_model) const
        {
        scitbx::af::shared<FloatType> dtda(f_model.size(), 0 );
        scitbx::af::shared<FloatType> dtdb(f_model.size(), 0 );
        CCTBX_ASSERT ( f_model.size() == calc_ori_lookup_table_.size() );

        FloatType aa,ba,ab,bb,obs,calc;
        FloatType t1,dqdaa,dqdba,dqdab,dqdbb;
        FloatType dqdt1, dt1daa,dt1dba,dt1dab,dt1dbb;

        long calc_index_a, calc_index_b;
        for (std::size_t ii=0;ii<i_obs_.size();ii++){
          calc_index_a = calc_ori_lookup_table_[ ii ];
          calc_index_b = calc_twin_lookup_table_[ ii ];
          aa = f_model[calc_index_a].real();
          ba = f_model[calc_index_a].imag();
          ab = f_model[calc_index_b].real();
          bb = f_model[calc_index_b].imag();
          calc = (1-alpha_)*(aa*aa + ba*ba) + alpha_*(ab*ab + bb*bb);
          obs = i_obs_[ii];
          t1 = (obs-calc);
          dt1daa = 2.0*aa*(1-alpha_);
          dt1dba = 2.0*ba*(1-alpha_);
          dt1dab = 2.0*ab*(alpha_);
          dt1dbb = 2.0*bb*(alpha_);
          dqdaa = -2.0*t1*dt1daa;
          dqdba = -2.0*t1*dt1dba;
          dqdab = -2.0*t1*dt1dab;
          dqdbb = -2.0*t1*dt1dbb;
          // place them in the correct positions please
          dtda[ calc_index_a ] += dqdaa;
          dtdb[ calc_index_a ] += dqdba;
          dtda[ calc_index_b ] += dqdab;
          dtdb[ calc_index_b ] += dqdbb;
        }
        scitbx::af::tiny<scitbx::af::shared<FloatType>,2> result(dtda,dtdb);
        return( result  );
      }

      scitbx::af::shared< std::complex<FloatType> > d_target_d_fmodel
      (scitbx::af::const_ref<std::complex<FloatType> > const& f_model){
        scitbx::af::shared<std::complex<FloatType> > result;

        CCTBX_ASSERT ( f_model.size() == calc_ori_lookup_table_.size() );
        scitbx::af::tiny<scitbx::af::shared<FloatType>, 2 > derivs;
        derivs =  d_target_d_ab( f_model );

        for (std::size_t ii=0;ii<f_model.size();ii++){
          std::complex<FloatType> tmp(derivs[0][ii],-derivs[1][ii] );
          result.push_back( tmp );
        }
        return result;
      }

      FloatType d_target_d_alpha
      (scitbx::af::const_ref<std::complex<FloatType> > const& f_model) const
      {
        FloatType result=0,aa,ba,ab,bb,obs,ia,ib;
        long calc_index_a, calc_index_b;
        for (std::size_t ii=0;ii<i_obs_.size();ii++){
          calc_index_a = calc_ori_lookup_table_[ ii ];
          calc_index_b = calc_twin_lookup_table_[ ii ];
          aa = f_model[calc_index_a].real();
          ba = f_model[calc_index_a].imag();
          ab = f_model[calc_index_b].real();
          bb = f_model[calc_index_b].imag();
          ia=aa*aa+ba*ba;
          ib=ab*ab+bb*bb;
          obs = i_obs_[ii];
          //std::cout << obs << " " << ia << " " << ib << " " << ( -(1.0-alpha_)*ia - alpha_*ib + obs ) << std::endl;
          result += 2.0*(ia-ib)*( -(1.0-alpha_)*ia - alpha_*ib + obs )*w_obs_[ii];
        }
        return result;
      }

      void alpha( FloatType tmp_alpha )
      {
         alpha_ = tmp_alpha;
      }

      FloatType alpha()
      {
         return(alpha_);
      }


 protected:
      scitbx::af::shared<FloatType> i_obs_;
      scitbx::af::shared<FloatType> w_obs_;

      scitbx::mat3<FloatType> twin_law_;
      cctbx::sgtbx::space_group space_group_;

      //scitbx::af::shared<cctbx::miller::index<> > hkl_calc_;
      FloatType alpha_;

      scitbx::af::shared<long> calc_ori_lookup_table_;
      scitbx::af::shared<long> calc_twin_lookup_table_;

 };










template<typename FloatType> class least_squares_hemihedral_twinning_on_f{
 public:
  // You want to use this constructor
    least_squares_hemihedral_twinning_on_f(
      scitbx::af::const_ref< cctbx::miller::index<> >  const& hkl_obs,       //1 indices for calculated data
      scitbx::af::const_ref< FloatType >               const& f_obs,         //2 f calc
      scitbx::af::const_ref< FloatType >               const& w_obs,         //3 f bulk solvent
      scitbx::af::const_ref< cctbx::miller::index<> >  const& hkl_calc,      //4 f_model; not const to avoid CV issues
      sgtbx::space_group                               const& space_group,   //5 space group
      bool                                             const& anomalous_flag,//6 anomalous_flag
      FloatType                                        const& alpha,         //7 twin fraction
      scitbx::mat3<FloatType>                          const& twin_law       //8 twin law
      ):
      space_group_( space_group ),
      twin_law_(twin_law),
      alpha_(alpha)
      {
        CCTBX_ASSERT( (alpha >=0) && (alpha<=0.50) );
        CCTBX_ASSERT( hkl_obs.size() > 0);
        CCTBX_ASSERT( hkl_obs.size() == f_obs.size() );
        CCTBX_ASSERT( (hkl_obs.size() == w_obs.size()) || (w_obs.size()==0) );

        cctbx::miller::lookup_utils::lookup_tensor<FloatType>
          tmp_lookup_object( hkl_calc, space_group, anomalous_flag  );
        int tmp_loc;
        for (std::size_t ii=0;ii<hkl_obs.size();ii++){
          f_obs_.push_back( f_obs[ii] );
          if (w_obs.size() > 0){
            w_obs_.push_back( w_obs[ii] );
          }
          else {
            w_obs_.push_back( 1.0 );
          }
          tmp_loc = tmp_lookup_object.find_hkl( hkl_obs[ii] );
          CCTBX_ASSERT( tmp_loc >= 0 );
          calc_ori_lookup_table_.push_back( tmp_loc );
          tmp_loc = tmp_lookup_object.find_hkl( twin_mate( hkl_obs[ii],twin_law ) );
          CCTBX_ASSERT( tmp_loc >= 0 );
          calc_twin_lookup_table_ .push_back( tmp_loc );



        }
        CCTBX_ASSERT( hkl_obs.size() <= hkl_calc.size() );
      }


      FloatType target(scitbx::af::const_ref<std::complex<FloatType> >
                       const& f_model) const
      {
        FloatType result=0,aa,ba,ab,bb,obs,calc;
        long calc_index_a, calc_index_b;
        for (std::size_t ii=0;ii<f_obs_.size();ii++){
          calc_index_a = calc_ori_lookup_table_[ ii ];
          calc_index_b = calc_twin_lookup_table_[ ii ];
          aa = f_model[calc_index_a].real();
          ba = f_model[calc_index_a].imag();
          ab = f_model[calc_index_b].real();
          bb = f_model[calc_index_b].imag();
          calc = std::sqrt((1-alpha_)*(aa*aa + ba*ba) + alpha_*(ab*ab + bb*bb));
          obs = f_obs_[ii];
          //std::cout << ii << " " << calc << " " << obs <<  " " << std::endl;
          result += w_obs_[ii]*(obs-calc)*(obs-calc);
        }
        return( result );
      }


      scitbx::af::tiny<scitbx::af::shared<FloatType>, 2 > d_target_d_ab
      (scitbx::af::const_ref<std::complex<FloatType> > const& f_model) const
        {
        scitbx::af::shared<FloatType> dtda(f_model.size(), 0 );
        scitbx::af::shared<FloatType> dtdb(f_model.size(), 0 );
        CCTBX_ASSERT ( f_model.size() == calc_ori_lookup_table_.size() );

        FloatType aa,ba,ab,bb,obs,calc;
        FloatType t1,dqdaa,dqdba,dqdab,dqdbb;
        FloatType dqdt1, dt1daa,dt1dba,dt1dab,dt1dbb;

        long calc_index_a, calc_index_b;
        for (std::size_t ii=0;ii<f_obs_.size();ii++){
          calc_index_a = calc_ori_lookup_table_[ ii ];
          calc_index_b = calc_twin_lookup_table_[ ii ];
          aa = f_model[calc_index_a].real();
          ba = f_model[calc_index_a].imag();
          ab = f_model[calc_index_b].real();
          bb = f_model[calc_index_b].imag();
          calc = std::sqrt( (1-alpha_)*(aa*aa + ba*ba) + alpha_*(ab*ab + bb*bb) );
          obs = f_obs_[ii];
          t1 = (obs-calc);
          if (calc>0){
            dt1daa = -aa*(1-alpha_)/calc;
            dt1dba = -ba*(1-alpha_)/calc;
            dt1dab = -ab*(alpha_)/calc;
            dt1dbb = -bb*(alpha_)/calc;
            dqdaa = 2.0*t1*dt1daa;
            dqdba = 2.0*t1*dt1dba;
            dqdab = 2.0*t1*dt1dab;
            dqdbb = 2.0*t1*dt1dbb;
          }else{
            dqdaa = 0;
            dqdba = 0;
            dqdab = 0;
            dqdbb = 0;
          }
          // place them in the correct positions please
          dtda[ calc_index_a ] += dqdaa;
          dtdb[ calc_index_a ] += dqdba;
          dtda[ calc_index_b ] += dqdab;
          dtdb[ calc_index_b ] += dqdbb;
        }
        scitbx::af::tiny<scitbx::af::shared<FloatType>,2> result(dtda,dtdb);
        return( result  );
      }

      scitbx::af::shared< std::complex<FloatType> > d_target_d_fmodel
      (scitbx::af::const_ref<std::complex<FloatType> > const& f_model){
        scitbx::af::shared<std::complex<FloatType> > result;

        CCTBX_ASSERT ( f_model.size() == calc_ori_lookup_table_.size() );
        scitbx::af::tiny<scitbx::af::shared<FloatType>, 2 > derivs;
        derivs =  d_target_d_ab( f_model );

        for (std::size_t ii=0;ii<f_model.size();ii++){
          std::complex<FloatType> tmp(derivs[0][ii],-derivs[1][ii] );
          result.push_back( tmp );
        }
        return result;
      }

      FloatType d_target_d_alpha
      (scitbx::af::const_ref<std::complex<FloatType> > const& f_model) const
      {
        FloatType result=0,aa,ba,ab,bb,obs,ia,ib,t1,dtda,calc;
        long calc_index_a, calc_index_b;
        for (std::size_t ii=0;ii<f_obs_.size();ii++){
          calc_index_a = calc_ori_lookup_table_[ ii ];
          calc_index_b = calc_twin_lookup_table_[ ii ];
          aa = f_model[calc_index_a].real();
          ba = f_model[calc_index_a].imag();
          ab = f_model[calc_index_b].real();
          bb = f_model[calc_index_b].imag();
          ia=aa*aa+ba*ba;
          ib=ab*ab+bb*bb;
          obs = f_obs_[ii];
          calc= std::sqrt( (1-alpha_)*ia + alpha_*ib );
          t1 = obs-calc;
          dtda=0;
          if (calc>0){
            dtda = -0.5*(ia-ib)/calc;
          }
          result += -2.0*t1*dtda;
        }
        return result;
      }

      void alpha( FloatType tmp_alpha )
      {
         alpha_ = tmp_alpha;
      }

      FloatType alpha()
      {
         return(alpha_);
      }


 protected:
      scitbx::af::shared<FloatType> f_obs_;
      scitbx::af::shared<FloatType> w_obs_;

      scitbx::mat3<FloatType> twin_law_;
      cctbx::sgtbx::space_group space_group_;

      //scitbx::af::shared<cctbx::miller::index<> > hkl_calc_;
      FloatType alpha_;

      scitbx::af::shared<long> calc_ori_lookup_table_;
      scitbx::af::shared<long> calc_twin_lookup_table_;

 };

 template<typename FloatType>
 class hemihedral_r_values
 {
   public:
   hemihedral_r_values(scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_obs,       //1 obs indices
                       scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_calc,      //2 calc indices
                       cctbx::sgtbx::space_group                       const& space_group,   //3 space_group
                       bool                                            const& anomalous_flag,//4 anomalous flag
                       scitbx::mat3<FloatType>                         const& twin_law       //5 twin law
                      )
   :
   obs_size_( hkl_obs.size() ),
   calc_size_( hkl_calc.size() )
   {
      cctbx::miller::lookup_utils::lookup_tensor<FloatType> tmp_lookup(hkl_calc, space_group, anomalous_flag);
      obs_in_calc_lookup_ = tmp_lookup.find_hkl( hkl_obs );
      int ht,kt,lt,tmp_location;
      for (long ii=0;ii<hkl_obs.size();ii++){
        CCTBX_ASSERT( obs_in_calc_lookup_[ii] >= 0 );
        cctbx::miller::index<> tmp_miller_index=twin_mate(hkl_obs[ii], twin_law);
        tmp_location = tmp_lookup.find_hkl( tmp_miller_index );
        CCTBX_ASSERT( tmp_location>=0 );
        twin_related_obs_in_calc_lookup_.push_back( tmp_location );
      }

   }

   FloatType r_intensity_abs( scitbx::af::const_ref<FloatType> const& i_obs,
                              scitbx::af::const_ref< std::complex<FloatType> > const& f_model,
                              FloatType const& twin_fraction )
   {
     CCTBX_ASSERT( obs_size_ == i_obs.size() );
     CCTBX_ASSERT( calc_size_ == f_model.size() );
     FloatType top=0,bottom=0,tmp_a,tmp_b, i_calc, tmp_location,tmp_twin_location;
     for (long ii=0;ii<obs_size_;ii++){
       tmp_location = obs_in_calc_lookup_[ii];
       tmp_a  = f_model[ tmp_location ].real();
       tmp_b  = f_model[ tmp_location ].imag();
       i_calc = (tmp_a*tmp_a + tmp_b*tmp_b)*(1-twin_fraction);

       tmp_location = twin_related_obs_in_calc_lookup_[ii];
       tmp_a  = f_model[ tmp_location ].real();
       tmp_b  = f_model[ tmp_location ].imag();
       i_calc+= (tmp_a*tmp_a + tmp_b*tmp_b)*twin_fraction;

       top+= std::fabs( i_calc - i_obs[ii]  );
       bottom+= std::fabs(i_obs[ii]);
     }

     FloatType result=0.0;

     if (bottom>0){
       result = top/bottom;
     }

     return (result);
   }



   FloatType r_intensity_sq( scitbx::af::const_ref<FloatType> const& i_obs,
                             scitbx::af::const_ref< std::complex<FloatType> > const& f_model,
                             FloatType const& twin_fraction )
   {
     CCTBX_ASSERT( obs_size_ == i_obs.size() );
     CCTBX_ASSERT( calc_size_ == f_model.size() );
     FloatType top=0,bottom=0,tmp_a,tmp_b, i_calc, tmp_location,tmp_twin_location;
     for (long ii=0;ii<obs_size_;ii++){
       tmp_location = obs_in_calc_lookup_[ii];
       tmp_a  = f_model[ tmp_location ].real();
       tmp_b  = f_model[ tmp_location ].imag();
       i_calc = (tmp_a*tmp_a + tmp_b*tmp_b)*(1-twin_fraction);

       tmp_location = twin_related_obs_in_calc_lookup_[ii];
       tmp_a  = f_model[ tmp_location ].real();
       tmp_b  = f_model[ tmp_location ].imag();
       i_calc+= (tmp_a*tmp_a + tmp_b*tmp_b)*twin_fraction;

       top+= (i_calc-i_obs[ii])*(i_calc-i_obs[ii]);
       bottom+= i_obs[ii]*i_obs[ii];
     }

     FloatType result=0.0;

     if (bottom>0){
       result = top/bottom;
     }
     return (result);
   }




   FloatType r_amplitude_abs( scitbx::af::const_ref<FloatType> const& f_obs,
                             scitbx::af::const_ref< std::complex<FloatType> > const& f_model,
                             FloatType const& twin_fraction )
   {
     CCTBX_ASSERT( obs_size_ == f_obs.size() );
     CCTBX_ASSERT( calc_size_ == f_model.size() );
     FloatType top=0,bottom=0,tmp_a,tmp_b, f_calc, tmp_location,tmp_twin_location;
     for (long ii=0;ii<obs_size_;ii++){
       tmp_location = obs_in_calc_lookup_[ii];
       tmp_a  = f_model[ tmp_location ].real();
       tmp_b  = f_model[ tmp_location ].imag();
       f_calc = (tmp_a*tmp_a + tmp_b*tmp_b)*(1-twin_fraction);

       tmp_location = twin_related_obs_in_calc_lookup_[ii];
       tmp_a  = f_model[ tmp_location ].real();
       tmp_b  = f_model[ tmp_location ].imag();
       f_calc+= (tmp_a*tmp_a + tmp_b*tmp_b)*twin_fraction;

       top+= std::fabs( std::sqrt(f_calc)-f_obs[ii] );
       bottom+= f_obs[ii]; // allways positive anyway
     }

     FloatType result=0.0;

     if (bottom>0){
       result = top/bottom;
     }

     return (result);
   }



   FloatType r_amplitude_sq( scitbx::af::const_ref<FloatType> const& f_obs,
                             scitbx::af::const_ref< std::complex<FloatType> > const& f_model,
                             FloatType const& twin_fraction )
   {
     CCTBX_ASSERT( obs_size_ == f_obs.size() );
     CCTBX_ASSERT( calc_size_ == f_model.size() );
     FloatType top=0,bottom=0,tmp_a,tmp_b, f_calc, tmp_location,tmp_twin_location;
     for (long ii=0;ii<obs_size_;ii++){
       tmp_location = obs_in_calc_lookup_[ii];
       tmp_a  = f_model[ tmp_location ].real();
       tmp_b  = f_model[ tmp_location ].imag();
       f_calc = (tmp_a*tmp_a + tmp_b*tmp_b)*(1-twin_fraction);

       tmp_location = twin_related_obs_in_calc_lookup_[ii];
       tmp_a  = f_model[ tmp_location ].real();
       tmp_b  = f_model[ tmp_location ].imag();
       f_calc+= (tmp_a*tmp_a + tmp_b*tmp_b)*twin_fraction;

       top+= (std::sqrt(f_calc)-f_obs[ii])*(std::sqrt(f_calc)-f_obs[ii]);
       bottom+= f_obs[ii]*f_obs[ii];
     }

     FloatType result=0.0;

     if (bottom>0){
       result = top/bottom;
     }

     return (result);
   }

   protected:

   scitbx::af::shared<long> obs_in_calc_lookup_;
   scitbx::af::shared<long> twin_related_obs_in_calc_lookup_;
   long obs_size_;
   long calc_size_;
 };





  template<typename FloatType>
  class hemihedral_detwinner
  {
    public:
    hemihedral_detwinner( scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_obs,
                          scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_calc,
                          cctbx::sgtbx::space_group                       const& space_group,
                          bool                                            const& anomalous_flag,
                          scitbx::mat3<FloatType>                         const& twin_law
                        )
    :
    twin_completeness_(0)
    {
       CCTBX_ASSERT( (hkl_obs.size() <= hkl_calc.size()) || (hkl_calc.size()==0) );
       cctbx::miller::lookup_utils::lookup_tensor<FloatType> tmp_obs(hkl_obs, space_group, anomalous_flag);
       cctbx::miller::lookup_utils::lookup_tensor<FloatType> tmp_calc(hkl_calc, space_group, anomalous_flag);
       long tmp_loc;
       // map out where the twin related amplitudes are
       for (std::size_t ii=0;ii<hkl_obs.size();ii++){
          tmp_loc = tmp_obs.find_hkl( twin_mate(hkl_obs[ii], twin_law) );
          // tmp_loc is allowed to be negative, a twin-complete data set is not gauranteed
          if( tmp_loc < 0 ){
            twin_completeness_ += 1.0;
          }
          obs_to_twin_obs_.push_back( tmp_loc );

          tmp_loc = tmp_calc.find_hkl( hkl_obs[ii] );
          CCTBX_ASSERT( tmp_loc >= 0 );
          obs_to_calc_.push_back( tmp_loc );

          tmp_loc = tmp_calc.find_hkl( twin_mate(hkl_obs[ii], twin_law) );
          CCTBX_ASSERT( tmp_loc >= 0 );
          obs_to_twin_calc_.push_back( tmp_loc );
       }
       twin_completeness_/=FloatType(hkl_obs.size());

       // do similar stuff for calculated data
       for (std::size_t ii=0;ii<hkl_calc.size();ii++){
         tmp_loc = tmp_calc.find_hkl( twin_mate(hkl_calc[ii],twin_law) );
         CCTBX_ASSERT( tmp_loc >=0 );
         calc_to_twin_calc_.push_back( tmp_loc );
       }

    }

    scitbx::af::tiny< scitbx::af::shared<FloatType>, 2 >
    detwin_with_twin_fraction(scitbx::af::const_ref<FloatType> const& i_obs,
                              scitbx::af::const_ref<FloatType> const& sig_obs,
                              FloatType const& twin_fraction) const
    {
      scitbx::af::shared<FloatType> i_detwin;
      scitbx::af::shared<FloatType> s_detwin;

      CCTBX_ASSERT( i_obs.size() == sig_obs.size() );
      CCTBX_ASSERT( i_obs.size() == obs_to_twin_obs_.size() );

      FloatType i_a,s_a,i_b,s_b, n_i, n_s;
      int tmp_loc;

      FloatType tmp_mult = std::sqrt( 1-2*twin_fraction +2*twin_fraction*twin_fraction)/(1.0-2.0*twin_fraction);

      for (std::size_t ii=0;ii<i_obs.size();ii++){
        tmp_loc = obs_to_twin_obs_[ii];
        n_i = 0.0; // new intensity
        n_s = 0.0; // new sigma
        if (tmp_loc>=0){
           i_a = i_obs[ ii ];
           s_a = sig_obs[ ii ];
           i_b = i_obs[ tmp_loc ];
           s_b = sig_obs[ tmp_loc ];

           n_i = ((1.0-twin_fraction)*i_a - twin_fraction*i_b)/(1-2.0*twin_fraction);
           n_s =  tmp_mult*std::sqrt((s_a*s_a*(1-twin_fraction) + s_b*s_b*twin_fraction));
        } else { // twin related reflection is not there. do 'equipartitioning'
           i_a = i_obs[ii];
           s_a = sig_obs[ii];
           n_i = i_a*(1-twin_fraction);
           n_s = s_a*(1-twin_fraction);
        }
        i_detwin.push_back( n_i );
        s_detwin.push_back( n_s );

      }
      CCTBX_ASSERT( i_detwin.size() == i_obs.size() );
      CCTBX_ASSERT( s_detwin.size() == sig_obs.size() );
      scitbx::af::tiny< scitbx::af::shared<FloatType>, 2> result( i_detwin, s_detwin );

      return( result );
    }


    scitbx::af::tiny< scitbx::af::shared<FloatType>, 2 >
    detwin_with_model_data(scitbx::af::const_ref<FloatType> const& i_obs,
                           scitbx::af::const_ref<FloatType> const& sig_obs,
                           scitbx::af::const_ref< std::complex<FloatType> > const& f_model,
                           FloatType const& twin_fraction) const
    {
       CCTBX_ASSERT( i_obs.size() == obs_to_twin_obs_.size() );
       CCTBX_ASSERT( i_obs.size() == sig_obs.size() );
       CCTBX_ASSERT( f_model.size() == calc_to_twin_calc_.size() );

       scitbx::af::shared<FloatType> detwinned_i;
       scitbx::af::shared<FloatType> detwinned_s;

       FloatType a, b, o_a, s_a, o_b, s_b, c_a, c_b, frac1, frac2, n_i, n_s;
       int loc_twin_obs, loc_calc, loc_twin_calc;
       for (std::size_t ii=0;ii<i_obs.size();ii++){
         loc_twin_obs = obs_to_twin_obs_[ ii ];
         loc_calc = obs_to_calc_[ ii ];
         loc_twin_calc = obs_to_twin_calc_[ ii ];
         o_a = i_obs[ii];
         o_b = i_obs[ loc_twin_obs ];
         s_a = sig_obs[ii];
         s_b = sig_obs[ loc_twin_obs ];

         a = f_model[ loc_calc ].real();
         b = f_model[ loc_calc ].imag();
         c_a = (a*a+b*b);

         a = f_model[ loc_twin_calc ].real();
         b = f_model[ loc_twin_calc ].imag();
         c_b = (a*a+b*b);

         frac1 = c_a * (1-twin_fraction)/ ( c_a*(1.0-twin_fraction) + c_b*twin_fraction );
         frac2 = c_a *twin_fraction/ ( c_b*(1.0-twin_fraction) + c_a*twin_fraction );

         n_i = o_a*frac1 + o_b*frac2;
         n_s = std::sqrt( s_a*s_a*frac1*frac1 + s_b*s_b*frac2*frac2 );

         detwinned_i.push_back( n_i );
         detwinned_s.push_back( n_s );
       }
       scitbx::af::tiny< scitbx::af::shared<FloatType>, 2 > result( detwinned_i, detwinned_s );
       return( result );
    }

    protected:
    scitbx::af::shared<long> obs_to_twin_obs_;
    scitbx::af::shared<long> obs_to_calc_;
    scitbx::af::shared<long> obs_to_twin_calc_;
    scitbx::af::shared<long> calc_to_twin_calc_;
    FloatType twin_completeness_;

  };





}}}

#endif // CCTBX_XRAY_TWIN_TARGETS
