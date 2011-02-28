#ifndef CCTBX_XRAY_TWIN_TARGETS_H
#define CCTBX_XRAY_TWIN_TARGETS_H

#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <complex>
#include <cmath>
#include <scitbx/math/bessel.h>
#include <scitbx/math/quadrature.h>
#include <scitbx/math/erf.h>
#include <scitbx/line_search/more_thuente_1994.h>
#include <cctbx/hendrickson_lattman.h>
#include <cctbx/miller/lookup_utils.h>
#include <cctbx/uctbx.h>
#include <cctbx/miller/asu.h>
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


  template<typename FloatType>
  class twin_completion{
  protected:
    scitbx::mat3<FloatType> twin_law_;
    bool anomalous_flag_;
    cctbx::sgtbx::space_group space_group_;
    scitbx::af::shared<cctbx::miller::index<> > hkl_;
    scitbx::af::shared<cctbx::miller::index<> > twin_hkl_;
    cctbx::miller::lookup_utils::lookup_tensor<FloatType> ori_lookup_table_;
  public:
    twin_completion(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
      sgtbx::space_group const& space_group,
      bool const& anomalous_flag,
      scitbx::mat3<FloatType> const& twin_law)
    :
      twin_law_(twin_law),
      anomalous_flag_(anomalous_flag),
      space_group_( space_group ),
      ori_lookup_table_(hkl,space_group,anomalous_flag)
    {
      CCTBX_ASSERT( hkl.size() > 0 );
      for (std::size_t ii=0; ii<hkl.size(); ii++){
        hkl_.push_back( hkl[ii] );
        twin_hkl_.push_back( twin_mate( hkl[ii], twin_law ) );
      }

    }

    scitbx::af::shared< cctbx::miller::index<> > twin_complete()
    {
      scitbx::af::shared< cctbx::miller::index<> > tmp;
      for (std::size_t ii=0;ii<hkl_.size();ii++){
        tmp.push_back( hkl_[ii] );
        tmp.push_back( twin_hkl_[ii] );
      }
      scitbx::af::shared< std::size_t > unique;
      unique = cctbx::miller::unique_under_symmetry_selection(
        sgtbx::space_group_type( space_group_ ),
        anomalous_flag_,
        tmp.const_ref() );

      scitbx::af::shared< cctbx::miller::index<> > unique_index;
      for (std::size_t ii=0;ii<unique.size();ii++){
        unique_index.push_back( tmp[ unique[ii] ] );
      }
      return(unique_index);
    }


    bool check_free_flags(scitbx::af::const_ref< bool > const& flags )
    {
      CCTBX_ASSERT( flags.size() == hkl_.size() );
      // loop over all flags
      bool ori,twin;
      for (std::size_t ii=0; ii<hkl_.size();ii++){
        ori = flags[ii];
        long tmp_loc = ori_lookup_table_.find_hkl( twin_hkl_[ii] );
        if (tmp_loc >= 0){
          twin = flags[ tmp_loc ];
          if (ori != twin ){ // they are not equal. This is a problem
            return (false);
          }
        }
      }
      return( true );
    }

    scitbx::af::shared<bool>  get_free_model_selection(scitbx::af::const_ref< cctbx::miller::index<> > hkl_calc,
                                                       scitbx::af::const_ref< bool > const& flags )
    {
      // Declare an array with results
      scitbx::af::shared<bool> result(hkl_calc.size(),0);
      for( std::size_t ii=0;ii<hkl_calc.size();ii++){
        long index = ori_lookup_table_.find_hkl( hkl_calc[ii] );
        if (index < 0 ){
          index = ori_lookup_table_.find_hkl( twin_mate(hkl_calc[ii], twin_law_) );
        }
        if (index < 0){
          // this means that neither hkl_calc or its twin mate is in the observed data
          // this means that in all reasonability, it is a 'free' reflection. Free reflections are marked as 'True'
          result[ii]=true;
        }
        else{
          CCTBX_ASSERT( index < flags.size() );
          result[ii] = flags[index];
        }
      }
      return(result);
    }

    scitbx::af::shared<FloatType> twin_sum(scitbx::af::const_ref< FloatType > data, FloatType const& alpha)
    {
      scitbx::af::shared<FloatType> result(hkl_.size(),0);
      FloatType a,b;
      for (std::size_t ii=0;ii<hkl_.size();ii++){
        a = data[ii];
        long indx = ori_lookup_table_.find_hkl( twin_hkl_[ii] );
        if (indx>=0){
          b = data[indx];
        } else {
          b = a;
        }
        result[ii] = (1-alpha)*a + alpha*b;
      }
        return(result);
    }
  };

  template<typename FloatType>
  class least_squares_hemihedral_twinning_on_i {
  protected:
    scitbx::af::shared<FloatType> i_obs_;
    scitbx::af::shared<FloatType> w_obs_;

    scitbx::mat3<FloatType> twin_law_;
    cctbx::sgtbx::space_group space_group_;

    FloatType alpha_;

    scitbx::af::shared<long> calc_ori_lookup_table_;
    scitbx::af::shared<long> calc_twin_lookup_table_;

  public:
    // You want to use this constructor
    least_squares_hemihedral_twinning_on_i(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_obs,
        //1 indices for calculated data
      scitbx::af::const_ref< FloatType > const& i_obs, //2 f calc
      scitbx::af::const_ref< FloatType > const& w_obs, //3 f bulk solvent
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl_calc,
        //4 f_model; not const to avoid CV issues
      sgtbx::space_group const& space_group,   //5 space group
      bool const& anomalous_flag, //6 anomalous_flag
      FloatType const& alpha, //7 twin fraction
      scitbx::mat3<FloatType> const& twin_law) //8 twin law
    :
      twin_law_(twin_law),
      space_group_( space_group ),
      alpha_(alpha)
    {
      CCTBX_ASSERT( (alpha >=0) && (alpha<=1.00) );
      CCTBX_ASSERT( hkl_obs.size() > 0);
      CCTBX_ASSERT( hkl_obs.size() == i_obs.size() );
      CCTBX_ASSERT( (hkl_obs.size() == w_obs.size()) || (w_obs.size()==0) );

      cctbx::miller::lookup_utils::lookup_tensor<FloatType>
        tmp_lookup_object( hkl_calc, space_group, anomalous_flag  );

      for (std::size_t ii=0;ii<hkl_obs.size();ii++){
        i_obs_.push_back( i_obs[ii] );
        if (w_obs.size() > 0){
          w_obs_.push_back( w_obs[ii] );
        }
        else {
          w_obs_.push_back( 1.0 );
        }
        long tmp_loc = tmp_lookup_object.find_hkl( hkl_obs[ii] );
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
      for (std::size_t ii=0;ii<i_obs_.size();ii++){
        long calc_index_a = calc_ori_lookup_table_[ ii ];
        long calc_index_b = calc_twin_lookup_table_[ ii ];
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

      FloatType aa,ba,ab,bb,obs,calc;
      FloatType t1,dqdaa,dqdba,dqdab,dqdbb;
      FloatType dt1daa,dt1dba,dt1dab,dt1dbb;

      for (std::size_t ii=0;ii<i_obs_.size();ii++){
        long calc_index_a = calc_ori_lookup_table_[ ii ];
        long calc_index_b = calc_twin_lookup_table_[ ii ];
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

      scitbx::af::tiny<scitbx::af::shared<FloatType>, 2 > derivs;
      derivs =  d_target_d_ab( f_model );

      for (std::size_t ii=0;ii<f_model.size();ii++){
        std::complex<FloatType> tmp(derivs[0][ii],derivs[1][ii] );
        result.push_back( tmp );
      }
      return result;
    }

    FloatType d_target_d_alpha
    (scitbx::af::const_ref<std::complex<FloatType> > const& f_model) const
    {
      FloatType result=0,aa,ba,ab,bb,obs,ia,ib;
      for (std::size_t ii=0;ii<i_obs_.size();ii++){
        long calc_index_a = calc_ori_lookup_table_[ ii ];
        long calc_index_b = calc_twin_lookup_table_[ ii ];
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

    void set_weights( scitbx::af::const_ref<FloatType> const& weights  ){
      for (std::size_t ii=0;ii<w_obs_.size();ii++){
        w_obs_[ii] = weights[ii];
      }
    }
  };

  template<typename FloatType> class least_squares_hemihedral_twinning_on_f{
  protected:
    scitbx::af::shared<FloatType> f_obs_;
    scitbx::af::shared<FloatType> w_obs_;

    scitbx::mat3<FloatType> twin_law_;
    cctbx::sgtbx::space_group space_group_;

    FloatType eps_;
    FloatType alpha_;

    scitbx::af::shared<long> calc_ori_lookup_table_;
    scitbx::af::shared<long> calc_twin_lookup_table_;
  public:
    // You want to use this constructor
    least_squares_hemihedral_twinning_on_f(
      scitbx::af::const_ref< cctbx::miller::index<> >  const& hkl_obs,       //1 indices for calculated data
      scitbx::af::const_ref< FloatType >               const& f_obs,         //2 f calc
      scitbx::af::const_ref< FloatType >               const& w_obs,         //3 weights
      scitbx::af::const_ref< cctbx::miller::index<> >  const& hkl_calc,      //4 f_model
      sgtbx::space_group                               const& space_group,   //5 space group
      bool                                             const& anomalous_flag,//6 anomalous_flag
      FloatType                                        const& alpha,         //7 twin fraction
      scitbx::mat3<FloatType>                          const& twin_law       //8 twin law
      ):
      twin_law_(twin_law),
      space_group_( space_group ),
      eps_(1e-5),
      alpha_(alpha)
      {
        CCTBX_ASSERT( (alpha >=0) && (alpha<=1.00) );
        CCTBX_ASSERT( hkl_obs.size() > 0);
        CCTBX_ASSERT( hkl_obs.size() == f_obs.size() );
        CCTBX_ASSERT( (hkl_obs.size() == w_obs.size()) || (w_obs.size()==0) );
        cctbx::miller::lookup_utils::lookup_tensor<FloatType>
          tmp_lookup_object( hkl_calc, space_group, anomalous_flag  );
        for (std::size_t ii=0;ii<hkl_obs.size();ii++){
          f_obs_.push_back( f_obs[ii] );
          if (w_obs.size() > 0){
            w_obs_.push_back( w_obs[ii] );
          }
          else {
            w_obs_.push_back( 1.0 );
          }
          long tmp_loc = tmp_lookup_object.find_hkl( hkl_obs[ii] );
          CCTBX_ASSERT( tmp_loc >= 0 );
          calc_ori_lookup_table_.push_back( tmp_loc );
          //--------------------------------------------------------------------------------------//
          //-- calc_ori_lookup[ ii ] -> jj :: ii: obs index ii has index equal to calc index jj --//
          //--------------------------------------------------------------------------------------//
          tmp_loc = tmp_lookup_object.find_hkl( twin_mate( hkl_obs[ii],twin_law ) );
          CCTBX_ASSERT( tmp_loc >= 0 ); // If this assertion fails, it means that a twin related calculated miler index is not
                                        // in the list of calculated / model indices. This definently should not happen!
          calc_twin_lookup_table_ .push_back( tmp_loc );
          //----------------------------------------------------------------------------------------------------//
          //-- calc_twin_lookup[ ii ] -> jj :: ii: obs index ii has twin related index equal to calc index jj --//
          //----------------------------------------------------------------------------------------------------//
        }
        CCTBX_ASSERT( hkl_obs.size() <= hkl_calc.size() );
      }


      FloatType target(scitbx::af::const_ref<std::complex<FloatType> >
                       const& f_model) const
      {
        FloatType result=0,aa=0,ba=0,ab=0,bb=0,obs=0,calc=0;
        for (std::size_t ii=0;ii<f_obs_.size();ii++){
          long calc_index_a = calc_ori_lookup_table_[ ii ];
          long calc_index_b = calc_twin_lookup_table_[ ii ];
          aa = f_model[calc_index_a].real();
          ba = f_model[calc_index_a].imag();
          ab = f_model[calc_index_b].real();
          bb = f_model[calc_index_b].imag();
          FloatType sqrt_arg = 0.;
          if(std::abs(aa)<1.e+10 &&
             std::abs(bb)<1.e+10 &&
             std::abs(ab)<1.e+10 &&
             std::abs(ba)<1.e+10) {
            sqrt_arg = (1-alpha_)*(aa*aa + ba*ba) + alpha_*(ab*ab + bb*bb);
          }
          obs = f_obs_[ii];
          if(sqrt_arg>0) {
            calc = std::sqrt(sqrt_arg);
            result += w_obs_[ii]*(obs-calc)*(obs-calc);
          }
        }
        return( result );
      }

      scitbx::af::tiny<scitbx::af::shared<FloatType>, 2 > d_target_d_ab
      (scitbx::af::const_ref<std::complex<FloatType> > const& f_model) const
        {
          scitbx::af::shared<FloatType> dtda(f_model.size(), 0 );
          scitbx::af::shared<FloatType> dtdb(f_model.size(), 0 );
          FloatType aa=0,ba=0,ab=0,bb=0,obs=0,calc=0;
          FloatType t1=0,dqdaa=0,dqdba=0,dqdab=0,dqdbb=0;
          FloatType dt1daa=0,dt1dba=0,dt1dab=0,dt1dbb=0;

          for (std::size_t ii=0;ii<f_obs_.size();ii++){
            long calc_index_a = calc_ori_lookup_table_[ ii ]; // we try to find calculated indices. They are complete.
            long calc_index_b = calc_twin_lookup_table_[ ii ];//  we try to find calculated indices. They are complete.
            CCTBX_ASSERT( calc_index_a >-1 );
            CCTBX_ASSERT( calc_index_b >-1 );

            aa = f_model[calc_index_a].real();
            ba = f_model[calc_index_a].imag();
            ab = f_model[calc_index_b].real();
            bb = f_model[calc_index_b].imag();
            FloatType sqrt_arg = 0.;
            if(std::abs(aa)<1.e+50 &&
               std::abs(ba)<1.e+50 &&
               std::abs(ab)<1.e+50 &&
               std::abs(bb)<1.e+50) {
              sqrt_arg = (1-alpha_)*(aa*aa + ba*ba) + alpha_*(ab*ab + bb*bb);
            }
            calc = 0;
            if(sqrt_arg>0) calc = std::sqrt(sqrt_arg);
            obs = f_obs_[ii];
            t1 = (obs-calc);

            if (calc>eps_){
              calc = std::sqrt(sqrt_arg);
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
            dtda[ calc_index_a ] += dqdaa*w_obs_[ii];
            dtdb[ calc_index_a ] += dqdba*w_obs_[ii];
            dtda[ calc_index_b ] += dqdab*w_obs_[ii];
            dtdb[ calc_index_b ] += dqdbb*w_obs_[ii];
          }
          scitbx::af::tiny<scitbx::af::shared<FloatType>,2> result(dtda,dtdb);
          return( result  );
        }

      scitbx::af::shared< std::complex<FloatType> > d_target_d_fmodel
      (scitbx::af::const_ref<std::complex<FloatType> > const& f_model){
        scitbx::af::shared<std::complex<FloatType> > result;

        scitbx::af::tiny<scitbx::af::shared<FloatType>, 2 > derivs;
        derivs =  d_target_d_ab( f_model );

        for (std::size_t ii=0;ii<f_model.size();ii++){
          std::complex<FloatType> tmp(derivs[0][ii],derivs[1][ii] );
          result.push_back( tmp );
        }
        return result;
      }

      FloatType d_target_d_alpha
      (scitbx::af::const_ref<std::complex<FloatType> > const& f_model) const
      {
        FloatType result=0,aa=0,ba=0,ab=0,bb=0,obs=0,ia=0,ib=0,t1=0,dtda=0,calc=0;
        for (std::size_t ii=0;ii<f_obs_.size();ii++){
          long calc_index_a = calc_ori_lookup_table_[ ii ];
          long calc_index_b = calc_twin_lookup_table_[ ii ];
          aa = f_model[calc_index_a].real();
          ba = f_model[calc_index_a].imag();
          ab = f_model[calc_index_b].real();
          bb = f_model[calc_index_b].imag();
          if(std::abs(aa)<1.e+50 &&
             std::abs(ba)<1.e+50 &&
             std::abs(ab)<1.e+50 &&
             std::abs(bb)<1.e+50) {
            ia=aa*aa+ba*ba;
            ib=ab*ab+bb*bb;
          }
          obs = f_obs_[ii];
          FloatType sqrt_arg = (1-alpha_)*ia + alpha_*ib;
          dtda=0;
          if (sqrt_arg>0){
            calc= std::sqrt(sqrt_arg);
            t1 = obs-calc;
            dtda = -0.5*(ia-ib)/calc;
            result += -2.0*t1*dtda*w_obs_[ii];
          }
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

      void set_weights(scitbx::af::const_ref<FloatType> const& weights  ){
        for (std::size_t ii=0;ii<w_obs_.size();ii++){
          w_obs_[ii] = weights[ii];
        }
      }
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
      for (long ii=0;ii<hkl_obs.size();ii++){
        CCTBX_ASSERT( obs_in_calc_lookup_[ii] >= 0 );
        cctbx::miller::index<> tmp_miller_index=twin_mate(hkl_obs[ii], twin_law);
        long tmp_location = tmp_lookup.find_hkl( tmp_miller_index );
        CCTBX_ASSERT( tmp_location>=0 );
        twin_related_obs_in_calc_lookup_.push_back( tmp_location );
      }

   }

   FloatType r_intensity_abs( scitbx::af::const_ref<FloatType> const& i_obs,
                              scitbx::af::const_ref< std::complex<FloatType> > const& f_model,
                              scitbx::af::const_ref<bool> const& selection,
                              FloatType const& twin_fraction )
   {
     CCTBX_ASSERT( obs_size_  == i_obs.size() );
     CCTBX_ASSERT( calc_size_ == f_model.size() );
     CCTBX_ASSERT( (obs_size_ == selection.size()) || (selection.size()==0)  );

     FloatType top=0,bottom=0,tmp_a,tmp_b, i_calc;
     bool use;
     for (long ii=0;ii<obs_size_;ii++){
       use = true;
       if (selection.size()>0){
           use = selection[ii];
       }
       if (use){
         long tmp_location = obs_in_calc_lookup_[ii];
         CCTBX_ASSERT( tmp_location>=0 );
         tmp_a  = f_model[ tmp_location ].real();
         tmp_b  = f_model[ tmp_location ].imag();
         i_calc = (tmp_a*tmp_a + tmp_b*tmp_b)*(1-twin_fraction);

         tmp_location = twin_related_obs_in_calc_lookup_[ii];
         CCTBX_ASSERT( tmp_location>=0 );
         tmp_a  = f_model[ tmp_location ].real();
         tmp_b  = f_model[ tmp_location ].imag();
         i_calc+= (tmp_a*tmp_a + tmp_b*tmp_b)*twin_fraction;

         top+= std::fabs( i_calc - i_obs[ii]  );
         bottom+= std::fabs(i_obs[ii]);
       }
     }

     FloatType result=0.0;

     if (bottom>0){
       result = top/bottom;
     }

     return (result);
   }



   FloatType r_intensity_sq( scitbx::af::const_ref<FloatType> const& i_obs,
                             scitbx::af::const_ref< std::complex<FloatType> > const& f_model,
                             scitbx::af::const_ref<bool> const& selection,
                             FloatType const& twin_fraction )
   {
     CCTBX_ASSERT( obs_size_ == i_obs.size() );
     CCTBX_ASSERT( calc_size_ == f_model.size() );
     CCTBX_ASSERT( (obs_size_ == selection.size()) || (selection.size()==0)  );

     FloatType top=0,bottom=0,tmp_a,tmp_b, i_calc;
     bool use;
     for (long ii=0;ii<obs_size_;ii++){
       use = true;
       if (selection.size()>0){
           use = selection[ii];
       }
       if (use){
         long tmp_location = obs_in_calc_lookup_[ii];
         CCTBX_ASSERT( tmp_location>=0 );
         tmp_a  = f_model[ tmp_location ].real();
         tmp_b  = f_model[ tmp_location ].imag();
         i_calc = (tmp_a*tmp_a + tmp_b*tmp_b)*(1-twin_fraction);

         tmp_location = twin_related_obs_in_calc_lookup_[ii];
         CCTBX_ASSERT( tmp_location>=0 );
         tmp_a  = f_model[ tmp_location ].real();
         tmp_b  = f_model[ tmp_location ].imag();
         i_calc+= (tmp_a*tmp_a + tmp_b*tmp_b)*twin_fraction;

         top+= (i_calc-i_obs[ii])*(i_calc-i_obs[ii]);
         bottom+= i_obs[ii]*i_obs[ii];
       }
     }

     FloatType result=0.0;

     if (bottom>0){
       result = top/bottom;
     }
     return (result);
   }




   FloatType r_amplitude_abs( scitbx::af::const_ref<FloatType> const& f_obs,
                              scitbx::af::const_ref< std::complex<FloatType> > const& f_model,
                              scitbx::af::const_ref<bool> const& selection,
                              FloatType const& twin_fraction )
   {
     CCTBX_ASSERT( obs_size_ == f_obs.size() );
     CCTBX_ASSERT( calc_size_ == f_model.size() );
     CCTBX_ASSERT( (obs_size_ == selection.size()) || (selection.size()==0)  );

     FloatType top=0,bottom=0,tmp_a,tmp_b=0, f_calc=0;
     bool use;
     for (long ii=0;ii<obs_size_;ii++){
       use = true;
       if (selection.size()>0){
           use = selection[ii];
       }
       if (use){
         long tmp_location = obs_in_calc_lookup_[ii];
         CCTBX_ASSERT( tmp_location>=0 );
         tmp_a  = f_model[ tmp_location ].real();
         tmp_b  = f_model[ tmp_location ].imag();
         if(std::abs(tmp_a)<1.e+50 && std::abs(tmp_b)<1.e+50) {
           f_calc = (tmp_a*tmp_a + tmp_b*tmp_b)*(1.0-twin_fraction);
         }
         tmp_location = twin_related_obs_in_calc_lookup_[ii];
         CCTBX_ASSERT( tmp_location>=0 );
         tmp_a  = f_model[ tmp_location ].real();
         tmp_b  = f_model[ tmp_location ].imag();
         if(std::abs(tmp_a)<1.e+50 && std::abs(tmp_b)<1.e+50) {
           f_calc+= (tmp_a*tmp_a + tmp_b*tmp_b)*twin_fraction;
         }
         if(f_calc >= 0) {
           top+= std::fabs( std::sqrt(f_calc)-f_obs[ii] );
         }
         bottom+= f_obs[ii]; // allways positive anyway
       }
     }

     FloatType result=0.0;

     if (bottom>0){
       result = top/bottom;
     }
     return (result);
   }



   FloatType r_amplitude_sq( scitbx::af::const_ref<FloatType> const& f_obs,
                             scitbx::af::const_ref< std::complex<FloatType> > const& f_model,
                             scitbx::af::const_ref<bool> const& selection,
                             FloatType const& twin_fraction )
   {
     CCTBX_ASSERT( obs_size_ == f_obs.size() );
     CCTBX_ASSERT( calc_size_ == f_model.size() );
     CCTBX_ASSERT( (obs_size_ == selection.size()) || (selection.size()==0)  );
     FloatType top=0,bottom=0,tmp_a,tmp_b, f_calc;
     bool use;
     for (long ii=0;ii<obs_size_;ii++){
       use = true;
       if (selection.size()>0){
           use = selection[ii];
       }
       if (use){
         long tmp_location = obs_in_calc_lookup_[ii];
         CCTBX_ASSERT( tmp_location>=0 );
         tmp_a  = f_model[ tmp_location ].real();
         tmp_b  = f_model[ tmp_location ].imag();
         f_calc = (tmp_a*tmp_a + tmp_b*tmp_b)*(1-twin_fraction);

         tmp_location = twin_related_obs_in_calc_lookup_[ii];
         CCTBX_ASSERT( tmp_location>=0 );
         tmp_a  = f_model[ tmp_location ].real();
         tmp_b  = f_model[ tmp_location ].imag();
         f_calc+= (tmp_a*tmp_a + tmp_b*tmp_b)*twin_fraction;

         top+= (std::sqrt(f_calc)-f_obs[ii])*(std::sqrt(f_calc)-f_obs[ii]);
         bottom+= f_obs[ii]*f_obs[ii];
       }
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
       obs_size_  = hkl_obs.size();
       calc_size_ = hkl_calc.size();

       cctbx::miller::lookup_utils::lookup_tensor<FloatType> tmp_obs(hkl_obs, space_group, anomalous_flag);
       cctbx::miller::lookup_utils::lookup_tensor<FloatType> tmp_calc(hkl_calc, space_group, anomalous_flag);
       // map out where the twin related amplitudes are
       for (std::size_t ii=0;ii<hkl_obs.size();ii++){
          long tmp_loc = tmp_obs.find_hkl( twin_mate(hkl_obs[ii], twin_law) );
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
       CCTBX_ASSERT(hkl_obs.size() != 0);
       twin_completeness_/=FloatType(hkl_obs.size());
       // do similar stuff for calculated data
       for (std::size_t ii=0;ii<hkl_calc.size();ii++){
         long tmp_loc = tmp_calc.find_hkl( twin_mate(hkl_calc[ii],twin_law) );
         //if this hkl_calc is observed
         //if (  (tmp_obs.find_hkl( hkl_calc[ii] )>0) ||
         //      (tmp_obs.find_hkl( twin_mate(hkl_calc[ii],twin_law))>=0) ){
         //  CCTBX_ASSERT( tmp_loc >=0 );
         // }
         calc_to_twin_calc_.push_back( tmp_loc );
       }

    }

    scitbx::af::shared<long> obs_to_twin_obs()
    {
      return ( obs_to_twin_obs_ );
    }

    scitbx::af::shared<long> obs_to_calc()
    {
      return ( obs_to_calc_ );
    }

    scitbx::af::shared<long> obs_to_twin_calc()
    {
      return ( obs_to_twin_calc_ );
    }

    scitbx::af::shared<long> calc_to_twin_calc()
    {
      return (calc_to_twin_calc_ );
    }



    scitbx::af::tiny< scitbx::af::shared<FloatType>, 2 >
    detwin_with_twin_fraction(scitbx::af::const_ref<FloatType> const& i_obs,
                              scitbx::af::const_ref<FloatType> const& sig_obs,
                              FloatType const& twin_fraction) const
    {
      scitbx::af::shared<FloatType> i_detwin;
      scitbx::af::shared<FloatType> s_detwin;

      CCTBX_ASSERT( i_obs.size() == obs_size_ );
      CCTBX_ASSERT( (sig_obs.size() == obs_size_) || (sig_obs.size()==0) );
      FloatType eps=1e-3;
      CCTBX_ASSERT( (twin_fraction<0.5-eps) || (twin_fraction>0.5+eps) );
      CCTBX_ASSERT( twin_fraction >= 0.0 );
      CCTBX_ASSERT( twin_fraction <= 1.0 );


      FloatType i_a,s_a,i_b,s_b, n_i, n_s;

      FloatType tmp_mult = std::sqrt( 1-2*twin_fraction +2*twin_fraction*twin_fraction)/(1.0-2.0*twin_fraction);

      for (std::size_t ii=0;ii<i_obs.size();ii++){
        long tmp_loc = obs_to_twin_obs_[ii];
        n_i = 0.0; // new intensity
        n_s = 0.0; // new sigma
        if (tmp_loc>=0){
           i_a = i_obs[ ii ];
           i_b = i_obs[ tmp_loc ];
           s_a = 0.0;
           s_b = 0.0;
           if (sig_obs.size()!=0){
             s_a = sig_obs[ ii ];
             s_b = sig_obs[ tmp_loc ];
           }
           n_i = ((1.0-twin_fraction)*i_a - twin_fraction*i_b)/(1.0-2.0*twin_fraction);
           n_s =  tmp_mult*std::sqrt((s_a*s_a*(1.0-twin_fraction) + s_b*s_b*twin_fraction));
        } else { // twin related reflection is not there. set things to a negative value
           i_a = i_obs[ii];
           s_a = 0;
           if ( sig_obs.size() > 0 ){
             s_a = sig_obs[ii];
           }
           n_i = -100000.0; // i_a*(1.0-twin_fraction)/(1-2.0*twin_fraction);
           n_s = s_a*10;
        }
        i_detwin.push_back( n_i );
        s_detwin.push_back( n_s );

      }
      CCTBX_ASSERT( i_detwin.size() == i_obs.size() );
      CCTBX_ASSERT( s_detwin.size() == i_obs.size() );
      scitbx::af::tiny< scitbx::af::shared<FloatType>, 2> result( i_detwin, s_detwin );

      return( result );
    }

    scitbx::af::tiny< scitbx::af::shared<FloatType>, 2 >
    twin_with_twin_fraction( scitbx::af::const_ref<FloatType> const& i_obs,
                             scitbx::af::const_ref<FloatType> const& sig_obs,
                             FloatType const& twin_fraction ) const
    {
      scitbx::af::shared<FloatType> i_twinned;
      scitbx::af::shared<FloatType> s_twinned;

      CCTBX_ASSERT( i_obs.size() == obs_size_ );
      CCTBX_ASSERT( (sig_obs.size() == 0) || (sig_obs.size() == obs_size_) );
      CCTBX_ASSERT( twin_fraction >= 0 );
      CCTBX_ASSERT( twin_fraction <= 1 );

      FloatType i_out, s_out, s_a, s_b;
      for (std::size_t ii=0; ii<obs_size_; ii++){
        long tmp_loc = obs_to_twin_obs_[ii];
        //CCTBX_ASSERT(tmp_loc<=obs_size_);
        i_out = i_obs[ii];
        s_out = 100;
        if (tmp_loc>=0){
          s_a = 0.0;
          s_b = 0.0;
          if ( sig_obs.size() > 0 ){
            s_a = sig_obs[ii];
            s_b = sig_obs[ tmp_loc ];
          }
          i_out = (1.0-twin_fraction)*i_obs[ii] + twin_fraction*i_obs[ tmp_loc ];
          s_out = std::sqrt( (1.0-twin_fraction)*s_a*s_a + twin_fraction*s_b*s_b );
        }
        i_twinned.push_back( i_out );
        s_twinned.push_back( s_out );
      }
      scitbx::af::tiny< scitbx::af::shared<FloatType>, 2 > result( i_twinned, s_twinned );
      return ( result );
    }


     scitbx::af::tiny< scitbx::af::shared<FloatType>, 2 >
     detwin_with_model_data(scitbx::af::const_ref<FloatType> const& i_obs,
                            scitbx::af::const_ref<FloatType> const& sig_obs,
                            scitbx::af::const_ref<FloatType> const& f_model,
                            FloatType const& twin_fraction) const
     {
        CCTBX_ASSERT( ( i_obs.size() == sig_obs.size() ) || ( sig_obs.size()==0  ) );
        CCTBX_ASSERT( f_model.size() == calc_size_ );
        CCTBX_ASSERT( i_obs.size() == obs_size_ );

        scitbx::af::shared<FloatType> detwinned_i;
        scitbx::af::shared<FloatType> detwinned_s;

        FloatType o_a, s_a, o_b, s_b, c_a, c_b, frac1, frac2, n_i, n_s;
        for (std::size_t ii=0;ii<i_obs.size();ii++){
          long loc_twin_obs = obs_to_twin_obs_[ ii ];
          long loc_calc = obs_to_calc_[ ii ];
          long loc_twin_calc = obs_to_twin_calc_[ ii ];
          n_i =  i_obs[ii]; //-100000;
          n_s = 1000000.0;
          if (loc_twin_obs >=0){
            if (loc_calc>=0){
              if (loc_twin_calc>=0){
                CCTBX_ASSERT(loc_twin_obs >=0 && loc_calc >=0 && loc_twin_calc >= 0); // XXX trap for long standing seg. fault bug
                o_a = i_obs[ii];
                CCTBX_ASSERT(i_obs.size() > loc_twin_obs); // XXX trap for long standing seg. fault bug
                o_b = i_obs[ loc_twin_obs ];
                s_a = 0;
                s_b = 0;
                if (sig_obs.size()>0){
                  s_a = sig_obs[ii];
                  CCTBX_ASSERT(sig_obs.size() > loc_twin_obs); // XXX trap for long standing seg. fault bug
                  s_b = sig_obs[ loc_twin_obs ];
                }
                CCTBX_ASSERT(f_model.size() > loc_calc); // XXX trap for long standing seg. fault bug
                c_a = f_model[ loc_calc ];
                c_a = c_a*c_a;
                CCTBX_ASSERT(f_model.size() > loc_twin_calc); // XXX trap for long standing seg. fault bug
                c_b = f_model[ loc_twin_calc ];
                c_b = c_b*c_b;

                frac1 = c_a * (1-twin_fraction) / ( c_a*(1.0-twin_fraction) + c_b*twin_fraction );
                frac2 = c_a * twin_fraction / ( c_b*(1.0-twin_fraction) + c_a*twin_fraction );

                n_i = o_a*frac1 + o_b*frac2;
                n_s = std::sqrt( s_a*s_a*frac1*frac1 + s_b*s_b*frac2*frac2 );
              }
            }
          }
          detwinned_i.push_back( n_i );
          detwinned_s.push_back( n_s );

        }
        scitbx::af::tiny< scitbx::af::shared<FloatType>, 2 > result( detwinned_i, detwinned_s );
        return( result );
     }


    scitbx::af::tiny< scitbx::af::shared<FloatType>, 2 >
    detwin_with_model_data(scitbx::af::const_ref<FloatType> const& i_obs,
                           scitbx::af::const_ref<FloatType> const& sig_obs,
                           scitbx::af::const_ref< std::complex<FloatType> > const& f_model,
                           FloatType const& twin_fraction) const
    {
       CCTBX_ASSERT( ( i_obs.size() == sig_obs.size() ) || ( sig_obs.size()==0 )  );
       CCTBX_ASSERT( f_model.size() == calc_size_ );
       CCTBX_ASSERT( i_obs.size() == obs_size_ );

       CCTBX_ASSERT( twin_fraction >= 0 );
       CCTBX_ASSERT( twin_fraction <= 1 );

       scitbx::af::shared<FloatType> detwinned_i;
       scitbx::af::shared<FloatType> detwinned_s;

       FloatType a, b, o_a, s_a, o_b, s_b, c_a, c_b, frac1, frac2, n_i, n_s;
       for (std::size_t ii=0;ii<i_obs.size();ii++){
         long loc_twin_obs = obs_to_twin_obs_[ ii ];
         long loc_calc = obs_to_calc_[ ii ];
         long loc_twin_calc = obs_to_twin_calc_[ ii ];
         n_i = i_obs[ii];
         n_s = 10000.0;
         if (loc_twin_obs >=0){
           if (loc_calc >=0){
             if (loc_twin_calc >= 0){
                CCTBX_ASSERT(loc_twin_obs >=0 && loc_calc >=0 && loc_twin_calc >= 0); // XXX trap for long standing seg. fault bug
               o_a = i_obs[ii];
                  CCTBX_ASSERT(i_obs.size() > loc_twin_obs); // XXX trap for long standing seg. fault bug
               o_b = i_obs[ loc_twin_obs ];
               s_a = 0;
               s_b = 0;
               if ( sig_obs.size() > 0 ){
                 s_a = sig_obs[ii];
                  CCTBX_ASSERT(sig_obs.size() > loc_twin_obs); // XXX trap for long standing seg. fault bug
                 s_b = sig_obs[ loc_twin_obs ];
               }
            CCTBX_ASSERT(f_model.size() > loc_calc); // XXX trap for long standing seg. fault bug
               a = f_model[ loc_calc ].real();
               b = f_model[ loc_calc ].imag();
               c_a = (a*a+b*b);
                  CCTBX_ASSERT(f_model.size() > loc_twin_calc); // XXX trap for long standing seg. fault bug
               a = f_model[ loc_twin_calc ].real();
               b = f_model[ loc_twin_calc ].imag();
               c_b = (a*a+b*b);

               frac1 = c_a * (1-twin_fraction) / ( c_a*(1.0-twin_fraction) + c_b*twin_fraction );
               frac2 = c_a * twin_fraction / ( c_b*(1.0-twin_fraction) + c_a*twin_fraction );
               n_i = o_a*frac1 + o_b*frac2;
               n_s = std::sqrt( s_a*s_a*frac1*frac1 + s_b*s_b*frac2*frac2 );
             }
           }
         }
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
    std::size_t calc_size_;
    std::size_t obs_size_;

  };


  template<typename FloatType>
  class single_twin_likelihood
  {
  public:
    single_twin_likelihood( FloatType const& i_obs1,  FloatType const& s_obs1,
                            FloatType const& i_obs2,  FloatType const& s_obs2,
                            FloatType const& f_calc1, FloatType const& f_calc2,
                            FloatType const& eps1,    FloatType const& eps2,
                            bool      const& centric1,bool      const& centric2,
                            FloatType const& a,       FloatType const& b,
                            FloatType const& tf,      int       const& n_quad
                           )
      {
        n_quad_ = n_quad;
        io1_=i_obs1;
        so1_=s_obs1;
        io2_=i_obs2;
        so2_=s_obs2;
        fc1_=f_calc1;
        fc2_=f_calc2;
        a_=a;
        b_=b;
        e1_=eps1;
        e2_=eps2;
        centric1_=centric1;
        centric2_=centric2;
        tf_=tf;

        // quadrature needed for numerical integration
        scitbx::math::quadrature::gauss_legendre_engine<FloatType> tmp_gle(n_quad_);
        x_quad_ = tmp_gle.x();
        w_quad_ = tmp_gle.w();
      }

    FloatType log_p(FloatType f1,FloatType f2){
      FloatType result;
      result = log_likelihood_single(io1_,io2_,so1_,so2_,
                                     f1,f2,fc1_,fc2_,tf_,
                                     centric1_,centric2_,
                                     a_,b_,e1_,e2_);
      return(result);
    }


    scitbx::af::tiny<FloatType,2>
    d_log_p_d_f(FloatType f1, FloatType f2)
    {
      scitbx::af::tiny<FloatType,2> result(0,0);
      result = first_der_single(io1_, io2_,
                                so1_, so2_,
                                f1,   f2,
                                fc1_, fc2_,
                                tf_,
                                centric1_,centric2_,
                                a_,   b_,
                                e1_,  e2_);
      return(result);
    }

    scitbx::af::tiny<FloatType,3>
    dd_log_p_dd_f(FloatType f1, FloatType f2)
    {
      scitbx::af::tiny<FloatType,3> result(0,0,0);
      result = snd_der_single(io1_, io2_,
                              so1_, so2_,
                              f1,   f2,
                              fc1_, fc2_,
                              tf_,
                              centric1_,centric2_,
                              a_,   b_,
                              e1_,  e2_);
      return(result);
    }

    FloatType num_integrate(FloatType fm1, FloatType s1, FloatType fm2, FloatType s2, FloatType sigma_level)
    {



      FloatType mid1, span1, a1, b1;
      FloatType mid2, span2, a2, b2;
      mid1=fm1;
      mid2=fm2;
      span1=s1*sigma_level;
      span2=s2*sigma_level;
      FloatType tmp_max;
      tmp_max = log_p( fm1, fm2 );

      // make sure we do not pass zero
      a1 = fm1-span1;
      if (a1<0){
        a1 = 0;
        b1 = fm1+span1;
        fm1 = (a1+b1)/2.0;
        span1 = fm1;
      }
      a2 = fm2-span2;
      if (a2<0){
        a2 = 0;
        b2 = fm1+span1;
        fm2 = (a2+b2)/2.0;
        span2 = fm2;
      }

      FloatType y1,y2, result1,result2,tmp;
      result1=0;
      for (int ii=0;ii<n_quad_;ii++){
        y1 = x_quad_[ii]*span1 + mid1;
        result2=0;
        for (int jj=0;jj<n_quad_;jj++){
          y2       = x_quad_[jj]*span2 + mid2;
          tmp      = std::exp( log_p( y1, y2 )-tmp_max );
          result2 += tmp*w_quad_[jj];
        }
        result2 = result2*span2;
        result1+= result2*w_quad_[ii];
      }

      result1 = result1*span1;
      result1 = std::exp(tmp_max)*result1;
      return( result1 );
    }

    FloatType laplace_integrate(FloatType fm1, FloatType fm2)
    {
      // first we need the second derivatives
      scitbx::af::tiny<FloatType,3> snd_der;
      snd_der = dd_log_p_dd_f(fm1,fm2);
      FloatType det, result;
      det = snd_der[0]*snd_der[1]-snd_der[2]*snd_der[2];
      det = std::fabs(det);
      result=std::exp( log_p(fm1,fm2) )*scitbx::constants::pi*2.0/std::sqrt(det);
      return result;
    }





  protected:
    FloatType log_likelihood_single(FloatType io1, FloatType io2,
                                    FloatType so1, FloatType so2,
                                    FloatType f1,  FloatType f2,
                                    FloatType fc1, FloatType fc2,
                                    FloatType tf,
                                    bool centric1, bool centric2,
                                    FloatType a,   FloatType b,
                                    FloatType eps1, FloatType eps2)
    {
      float pm1=0, pm2=0, o12=0;
      if (centric1){
        pm1 = centric_log_likelihood_model(f1,fc1,a,b,eps1);
      } else{
        pm1 = acentric_log_likelihood_model(f1,fc1,a,b,eps1);
      }

      if (centric2){
        pm2 = centric_log_likelihood_model(f2,fc2,a,b,eps2);
      } else {
        pm2 = acentric_log_likelihood_model(f2,fc2,a,b,eps2);
      }

      o12 = q(io1,io2,so1,so2,f1,f2,tf);
      return( pm1+pm2+o12 );
    }

    scitbx::af::tiny<FloatType,2>
    first_der_single(FloatType io1, FloatType io2,
                     FloatType so1, FloatType so2,
                     FloatType f1,  FloatType f2,
                     FloatType fc1, FloatType fc2,
                     FloatType tf,
                     bool centric1, bool centric2,
                     FloatType a,   FloatType b,
                     FloatType eps1,FloatType eps2)
    {
      scitbx::af::tiny<FloatType,2> result(0,0);
      FloatType result1=0,result2=0;

      if (centric1){
        result1 = calc_fst_der_centric_model(f1,fc1,a,b,eps1);
      } else {
        result1 = calc_fst_der_acentric_model(f1,fc1,a,b,eps1);
      }

      if (centric2){
        result2 = calc_fst_der_centric_model(f2,fc2,a,b,eps2);
      } else {
        result2 = calc_fst_der_acentric_model(f2,fc2,a,b,eps2);
      }

      result1 += dq1(io1,io2,so1,so2,f1,f2,tf);
      result2 += dq2(io1,io2,so1,so2,f1,f2,tf);

      result[0]=result1;
      result[1]=result2;
      return( result );
    }


    scitbx::af::tiny<FloatType,3>
    snd_der_single(FloatType io1, FloatType io2,
                   FloatType so1, FloatType so2,
                   FloatType f1,  FloatType f2,
                   FloatType fc1, FloatType fc2,
                   FloatType tf,
                   bool centric1, bool centric2,
                   FloatType a,   FloatType b,
                   FloatType eps1,FloatType eps2)
    {
      scitbx::af::tiny<FloatType,3> result(0,0,0);
      FloatType result11=0,result22=0, result12=0;

      if (centric1){
        result11 = calc_snd_der_centric_model(f1,fc1,a,b,eps1);
      } else {
        result11 = calc_snd_der_acentric_model(f1,fc1,a,b,eps1);
      }

      if (centric2){
        result22 = calc_snd_der_centric_model(f2,fc2,a,b,eps2);
      } else {
        result22 = calc_snd_der_acentric_model(f2,fc2,a,b,eps2);
      }

      result11 += dq11(io1,io2,so1,so2,f1,f2,tf);
      result22 += dq22(io1,io2,so1,so2,f1,f2,tf);
      result12  = dq12(so1,so2,f1,f2,tf);

      result[0]=result11;
      result[1]=result22;
      result[2]=result12;
      return( result );
    }

   inline FloatType q(FloatType io1, FloatType io2,
                      FloatType so1, FloatType so2,
                      FloatType f1,  FloatType f2,
                      FloatType tf)
   {
     FloatType result=0,norma,tmp1,tmp2;
     norma=2.0*scitbx::constants::pi*so1*so2;
     FloatType ic1, ic2;
     ic1=(1-tf)*f1*f1+tf*f2*f2;
     ic2=(1-tf)*f2*f2+tf*f1*f1;
     tmp1=-(io1-ic1)*(io1-ic1)/(2.0*so1*so1);
     tmp2=-(io2-ic2)*(io2-ic2)/(2.0*so2*so2);
     result = -std::log(norma)+tmp1+tmp2;
     return( result );
   }

   inline FloatType dq1(FloatType io1, FloatType io2,
                        FloatType so1, FloatType so2,
                        FloatType f1,  FloatType f2,
                        FloatType tf)
   {
     FloatType result=0,tmp1,tmp2;
     FloatType ic1, ic2;
     ic1=(1-tf)*f1*f1+tf*f2*f2;
     ic2=(1-tf)*f2*f2+tf*f1*f1;
     tmp1=2*(io1-ic1)*(1-tf)*f1/(so1*so1);
     tmp2=2*(io2-ic2)*tf*f1/(so2*so2);
     result = tmp1+tmp2;
     return( result );
   }

   inline FloatType dq2(FloatType io1, FloatType io2,
                        FloatType so1, FloatType so2,
                        FloatType f1,  FloatType f2,
                        FloatType tf)
   {
     FloatType result=0,tmp1,tmp2;
     FloatType ic1, ic2;
     ic1=(1-tf)*f1*f1+tf*f2*f2;
     ic2=(1-tf)*f2*f2+tf*f1*f1;
     tmp1=2*(io1-ic1)*tf*f2/(so1*so1);
     tmp2=2*(io2-ic2)*(1-tf)*f2/(so2*so2);
     result = tmp1+tmp2;
     return( result );
   }


   inline FloatType dq11(FloatType io1, FloatType io2,
                         FloatType so1, FloatType so2,
                         FloatType f1,  FloatType f2,
                         FloatType tf)
   {
     FloatType result=0,tmp1,tmp2;
     tmp1= -2.0*(-1+tf)*(io1+3*f1*f1*(-1+tf)-f2*f2*tf)/(so1*so1);
     tmp2=  2.0*tf*(io2+f2*f2*(-1+tf)-3*f1*f1*tf)/(so2*so2);
     result = tmp1+tmp2;
     return( result );
   }


   inline FloatType dq22(FloatType io1, FloatType io2,
                         FloatType so1, FloatType so2,
                         FloatType f1,  FloatType f2,
                         FloatType tf)
   {
     FloatType result=0,tmp1,tmp2;
     tmp1= -2.0*(-1+tf)*(io2+3*f2*f2*(-1+tf)-f1*f1*tf)/(so2*so2);
     tmp2=  2.0*tf*(io1+f1*f1*(-1+tf)-3*f2*f2*tf)/(so1*so1);
     result = tmp1+tmp2;
     return( result );
   }


   inline FloatType dq12(FloatType so1, FloatType so2,
                         FloatType f1,  FloatType f2,
                         FloatType tf)
   {
     FloatType result;
     result = 4*f1*f2*(so1*so1+so2*so2)*(-1+tf)*tf/(so1*so1*so2*so2 );
     return( result );
   }

   inline FloatType acentric_log_likelihood_model(FloatType fo, FloatType fc, FloatType a, FloatType b, FloatType e){
      FloatType result;
      FloatType eb=e*b;
      if (fo<=1e-13){
        fo=1e-13;
      }
      FloatType x=2.0*a*fo*fc/eb;
      FloatType exparg; // = (fo*fo + alpha_[ii]*alpha_[ii]*f_calc_[ii]*f_calc_[ii]);
      //exparg = exparg/eb;
      //result = std::log(2.0) + std::log(fo) - std::log(eb)  -exparg  + std::log( scitbx::math::bessel::i0(x) );
      //FloatType result2;
      exparg = fo -  a*fc;
      exparg = exparg*exparg/eb;
      result = std::log(2.0) + std::log(fo) - std::log(eb)  -exparg + std::log( scitbx::math::bessel::ei0(x));
      return (result);
    }


    inline FloatType centric_log_likelihood_model(FloatType fo,  FloatType fc, FloatType a, FloatType b, FloatType e){
      FloatType result;
      FloatType eb=e*b;
      if (fo<=1e-13){
        fo=1e-13;
      }
      FloatType x=a*fo*fc/eb;
      FloatType exparg = (fo*fo + a*a*fc*fc)/(2.0*eb);
      FloatType tmp;
      if (x>40){
         tmp = x*0.999921 - 0.65543;
      } else {
         tmp = std::log( std::cosh(x) );
      }

      result = 0.5*std::log(2.0)-0.5*std::log(scitbx::constants::pi) - 0.5*std::log(eb)
        -exparg + tmp;
      return(result);
    }

    //--------------------------------------------------------------
    // compute the first derivative of the loglikelihood function
    inline FloatType
    calc_fst_der_acentric_model( FloatType fo, FloatType fc, FloatType a, FloatType b, FloatType e)
    {
      FloatType result;
      FloatType eb=e*b;
      if (fo<=1e-13){
        fo=1e-13;
      }
      FloatType x = 2.0*a*fo*fc/(eb);

      FloatType m = scitbx::math::bessel::i1_over_i0(x);

      result = (1.0/fo) - (2.0*fo/eb) +(2.0*a*fc/eb)*m;
      return (result);
    }

    // derivative of centrics
    inline FloatType
    calc_fst_der_centric_model( FloatType fo, FloatType fc, FloatType a, FloatType b, FloatType e)
    {
      FloatType result;
      FloatType eb=e*b;
      if (fo<=1e-13){
        fo=1e-13;
      }
      FloatType x = a*fo*fc/(eb);
      if (x<1e-13){
        x=1e-13;
      }
      FloatType m = std::tanh(x)*a*fc/(eb);;
      result = -fo/eb + m;
      return (result);
    }

    // second derivatives
    inline FloatType
    calc_snd_der_acentric_model( FloatType fo, FloatType fc, FloatType a, FloatType b, FloatType e)
    {
      FloatType result;
      FloatType eb=e*b;
      if (fo<=1e-13){
        fo=1e-13;
      }
      FloatType x = 2.0*a*fo*fc/(eb);
      FloatType m = scitbx::math::bessel::i1_over_i0(x);
      if (x<1e-13){
        x = 1e-13;
      }
      result = - (1.0/(fo*fo))
               - (2.0/eb)
               + (fc*4.0*a*a/(eb*eb))*
        (1.0-m/(x)-m*m);
      return (result);
    }

    inline FloatType calc_snd_der_centric_model( FloatType fo, FloatType fc, FloatType a, FloatType b, FloatType e )
    {
      FloatType result=0;
      FloatType eb=e*b;
      if (fo<=1e-13){
        fo=1e-13;
      }
      FloatType x = a*fo*fc/(eb);
      FloatType m = std::tanh( x );
      result = -1.0/eb + a*a*fc*fc*(1.0-m*m)/(eb*eb);
      return (result);
    }


    // the computation of the mean can be used as a good starting point for maximisation of the (partial) likelihood
    inline FloatType compute_mean_model( FloatType fo, FloatType fc, FloatType a, FloatType b, FloatType e, bool centric  ){
      FloatType result, tmp;
      if (centric){
        tmp = a*fc/std::sqrt(2.0*e*b);
        result = std::exp(-tmp*tmp)*std::sqrt(2.0*e*b/scitbx::constants::pi);
        result = result + a*fc*scitbx::math::erf(tmp);
      }
      else {
        FloatType x;
        x = a*fc*a*fc/(2.0*e*b);
        result = (e*b + a*a*fc*fc)*scitbx::math::bessel::ei0(x);
        result+= a*a*fc*fc*scitbx::math::bessel::ei1(x);
        result = result*0.5*std::sqrt(scitbx::constants::pi/(e*b));
      }
      return result;

    }
    //--------------------------------------------------------------
    FloatType io1_, so1_;
    FloatType io2_, so2_;
    FloatType fc1_, fc2_;
    FloatType a_, b_, e1_, e2_,tf_;
    bool centric1_, centric2_;
    int n_quad_;
    scitbx::af::shared<FloatType> x_quad_;
    scitbx::af::shared<FloatType> w_quad_;
  };

}}}

#endif // CCTBX_XRAY_TWIN_TARGETS
