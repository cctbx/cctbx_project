#ifndef CCTBX_XRAY_GROUPED_DATA_H
#define CCTBX_XRAY_GROUPED_DATA_H

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

namespace cctbx { namespace xray { namespace grouped_data {

  template<typename FloatType>
  class unmerged_data{
    protected:
      scitbx::af::shared< cctbx::miller::index<> > hkl_obs_;
      scitbx::af::shared< cctbx::miller::index<> > asu_hkl_;
        // these are ONLY indices in the ASU, in the same order as
        // Fmodel/Fcalc for instance
      scitbx::af::shared< long int > map_hkl_obs_to_asu_hkl_;
        // maps hkl_obs_ position to asu_hkl_ position
        // (includes symmetry ops to asu)
      scitbx::af::shared< std::vector<long int> > map_asu_hkl_to_hkl_obs_;
        // an array of observations that is associated with a particular
      cctbx::miller::lookup_utils::lookup_tensor<FloatType> asu_lookup_table_;
        // a lookup function that finds for any index the corresponding
        // asu index in a supplied array
      sgtbx::space_group space_group_;
      bool anomalous_flag_;

    public:
      unmerged_data(
        scitbx::af::const_ref< cctbx::miller::index<> > hkl_obs,
        scitbx::af::const_ref< cctbx::miller::index<> > asu_hkl,
        sgtbx::space_group const& space_group,
        bool const& anomalous_flag)
      :
        asu_lookup_table_(asu_hkl,space_group,anomalous_flag),
        space_group_(space_group),
        anomalous_flag_(anomalous_flag)
      {
        long int tmp_index;

        // push back asu hkl's and init map_asu_hkl_to_hkl_obs
        for (int ii=0;ii<asu_hkl.size();ii++ ){
          asu_hkl_.push_back( asu_hkl[ii] );
          std::vector<long> tmp_tab;
          map_asu_hkl_to_hkl_obs_.push_back( tmp_tab );
        }

        // make hkl_obs and associated lookup tables
        for (int ii=0; ii<hkl_obs.size();ii++){
          // push back the original index
          hkl_obs_.push_back( hkl_obs[ii] );
          // find this index in the asu please
          tmp_index = asu_lookup_table_.find_hkl( hkl_obs[ii] );
          CCTBX_ASSERT( tmp_index >= 0);
          CCTBX_ASSERT( tmp_index < asu_hkl.size() );
          map_hkl_obs_to_asu_hkl_.push_back( tmp_index );
          map_asu_hkl_to_hkl_obs_[ tmp_index ].push_back( ii );
        }

        // all lookup tables have been initialzed, that's it
        //for (int ii=0;ii<asu_hkl_.size();ii++){
          //std::cout << ii << "--->  ";
          //for (int jj=0;jj<map_asu_hkl_to_hkl_obs_[ii].size();jj++){
            //std::cout << map_asu_hkl_to_hkl_obs_[ii][jj] << " ";
          //}
          //std::cout << std::endl;
        //}


      }

      scitbx::af::shared< std::vector<long int> > map_asu_hkl_to_hkl_obs()
      {
        return (map_asu_hkl_to_hkl_obs_);
      }
  };

  /*
   * Below some applications for the unmerged data objects are found
   */

  template <typename FloatType>
  scitbx::af::shared<FloatType>
  kernel_mean(scitbx::af::const_ref<FloatType> const& d_star_sq_hkl,
              scitbx::af::const_ref<FloatType> const& I_hkl,
              scitbx::af::const_ref<FloatType> const& epsilon_hkl,
              scitbx::af::const_ref<FloatType> const& d_star_sq_array,
              FloatType const& kernel_width)
  {
    SCITBX_ASSERT(d_star_sq_hkl.size() == I_hkl.size() );
    SCITBX_ASSERT(d_star_sq_hkl.size() == epsilon_hkl.size() );

    scitbx::af::shared<FloatType> norma_I_array(d_star_sq_array.size(),0);
    scitbx::af::shared<FloatType> weights_array(d_star_sq_array.size(),0);
    FloatType x,dx,tmp_norm,result;
    FloatType eps=1E-8;

    // Use a simple kernel for binning purposes
    for (unsigned jj=0;jj<d_star_sq_hkl.size();jj++){
      x = d_star_sq_hkl[jj];
      for (unsigned ii=0 ;ii<d_star_sq_array.size();ii++){
        if (I_hkl[jj]>0){
          dx = x-d_star_sq_array[ii];
          dx = (dx*dx)/(2.0*kernel_width*kernel_width);
          result = std::exp(-dx);
          weights_array[ii] +=  result;
          norma_I_array[ii] +=  I_hkl[jj]*result/epsilon_hkl[jj];
        }
      }
    }
    // Now we just have to 'normalise' via the weights we have obtained
    for (unsigned ii=0 ;ii<d_star_sq_array.size();ii++){
      tmp_norm = weights_array[ii];
      if (tmp_norm <= eps){
        tmp_norm=eps;
      }
      norma_I_array[ii]/=tmp_norm;
    }
    return(norma_I_array);
  }

  template<typename FloatType>
  class merger{
    public:
      merger(scitbx::af::const_ref< cctbx::miller::index<> > hkl_obs,
               scitbx::af::const_ref< FloatType >              i_obs,
               scitbx::af::const_ref< FloatType >              s_obs,
               sgtbx::space_group                              const& space_group,
               bool                                            const& anomalous_flag,
               cctbx::uctbx::unit_cell                         const& uc):
      space_group_( space_group ),
      anomalous_flag_( anomalous_flag ),
      unit_cell_( uc )
      {
        CCTBX_ASSERT( hkl_obs.size() ==  i_obs.size() );
        CCTBX_ASSERT( hkl_obs.size() ==  s_obs.size() );

        // first copy over the data please
        for (int ii=0;ii<hkl_obs.size();ii++){
          hkl_obs_.push_back( hkl_obs[ii] );
          i_obs_.push_back( i_obs[ii] );
          s_obs_.push_back( s_obs[ii] );
          epsi_.push_back( space_group_.epsilon( hkl_obs_[ii] ) );
          centric_.push_back( space_group_.is_centric(hkl_obs_[ii]) );
          d_star_sq_.push_back( unit_cell_.d_star_sq( hkl_obs_[ii])  );
        }


        // now we need unique data
        scitbx::af::shared< std::size_t > unique;
        unique = cctbx::miller::unique_under_symmetry_selection(
          sgtbx::space_group_type( space_group_ ),
          anomalous_flag_,
          hkl_obs_.const_ref() );
        for (int ii=0; ii<unique.size();ii++){
          asu_hkl_.push_back( hkl_obs_[ unique[ii] ] );
        }

        // and please make the lookup tables
        unmerged_data<FloatType> umd_lut( hkl_obs_.const_ref(), asu_hkl_.const_ref(), space_group_, anomalous_flag_ );
        map_asu_to_obs_ = umd_lut.map_asu_hkl_to_hkl_obs();
      }

      // Likelihood bit
      FloatType log_lik_merge_single( FloatType F_true, long int this_index )
      {
        FloatType result=0, I_true, this_i, this_s;
        I_true = F_true*F_true;
        for (int ii=0;ii<map_asu_to_obs_[this_index].size();ii++){
          this_i  = i_obs_[ map_asu_to_obs_[this_index][ii] ];
          this_s  = s_obs_[ map_asu_to_obs_[this_index][ii] ];
          result += -(this_i-I_true)*(this_i-I_true)/(this_s*this_s*2.0) - std::log( scitbx::constants::pi*2.0 )*0.5 - std::log(this_s) ;
        }
        return( result );
      }

      std::vector<FloatType> combine_quadratics(FloatType ma, FloatType sa, FloatType mb, FloatType sb)
      {
        FloatType new_mean, new_var, weight;
        FloatType denominator = sa*sa+sb*sb;
        SCITBX_ASSERT(denominator != 0);
        new_mean = (sa*sa*mb + sb*sb*ma)/denominator;
        new_var  = (sa*sa*sb*sb)/denominator;
        weight   = -(ma-mb)*(ma-mb)/(2.0*denominator);
        std::vector<FloatType> result;
        result.push_back( new_mean );
        result.push_back( std::sqrt(new_var) );
        result.push_back( weight );
        return(result);
      }

      std::vector<FloatType> combine_obs(long int this_index)
      {
        FloatType new_mean, new_sig, weight=0.0, norma=0,this_i, this_s, dev_abs=0, dev_sq=0;

        new_mean = i_obs_[ map_asu_to_obs_[this_index][0] ];
        new_sig  = s_obs_[ map_asu_to_obs_[this_index][0] ];
        SCITBX_ASSERT(new_sig > 0); // FIXME rjgildea 2012-08-28
        norma = -std::log(2.0*scitbx::constants::pi)/2.0-std::log(new_sig);
        std::vector<FloatType> tmp_result;
        for (int ii=1;ii<map_asu_to_obs_[this_index].size();ii++){
          this_i  = i_obs_[ map_asu_to_obs_[this_index][ii] ];
          this_s  = s_obs_[ map_asu_to_obs_[this_index][ii] ];
          tmp_result = combine_quadratics(new_mean, new_sig, this_i, this_s);
          weight += tmp_result[2];
          new_mean = tmp_result[0];
          new_sig = tmp_result[1];
          SCITBX_ASSERT(this_s > 0);
          norma+=-std::log(2.0*scitbx::constants::pi)/2.0-std::log(this_s);
        }
        std::vector<FloatType> result;
        result.push_back( new_mean );
        result.push_back( new_sig );
        result.push_back( weight );
        result.push_back( norma );
        result.push_back( map_asu_to_obs_[this_index].size() );

        FloatType tmp;
        for (int ii=1;ii<map_asu_to_obs_[this_index].size();ii++){
          tmp = i_obs_[ map_asu_to_obs_[this_index][ii] ]-new_mean;
          dev_abs += std::abs( tmp );
          dev_sq  += tmp*tmp;
        }
        result.push_back( dev_abs );
        result.push_back( dev_sq );

        return( result );
        // output map
        // 0: combined mean intensity
        // 1: combined sigma
        // 2: log weight (chi squares)
        // 3: log normalisation constant from gaussian
        // 4: number of contributors
        // 5: sum |I-<I>|
        // 6: sum (I-<I>)^2

      }

      FloatType bic()
      {
        // Here we compute the 'Schwartz's Bayesian Information Criterion'
        // to get an approximate estimate of the marginal log-likelihood
        FloatType result=0;
        std::vector<FloatType> tmp_result;
        for (long int ii=0;ii<asu_hkl_.size();ii++){
          tmp_result = combine_obs( ii );
          result += tmp_result[2];// - 0.5*std::log( tmp_result[4] );
        }
        SCITBX_ASSERT(i_obs_.size() > 0);
        return (result-asu_hkl_.size()
          * std::log(static_cast<FloatType>(i_obs_.size()))*0.5 );
      }

      FloatType r_abs()
      {
        FloatType top=0, bottom=0;

        std::vector<FloatType> tmp_result;
        for (long int ii=0;ii<asu_hkl_.size();ii++){
          tmp_result = combine_obs( ii );
          SCITBX_ASSERT(tmp_result[4] != 0);
          top += tmp_result[5]/tmp_result[4];
          bottom += tmp_result[0];
        }
        return top/(std::max(bottom,1e-12));
      }

    protected:
      scitbx::af::shared< cctbx::miller::index<> >           hkl_obs_;
      scitbx::af::shared< FloatType >                          i_obs_;
      scitbx::af::shared< FloatType >                          s_obs_;
      scitbx::af::shared< FloatType >                           epsi_;
      //scitbx::af::shared< FloatType >                         mean_i_;
      scitbx::af::shared< bool >                             centric_;
      scitbx::af::shared< FloatType >                      d_star_sq_;

      scitbx::af::shared< cctbx::miller::index<> >           asu_hkl_;
      scitbx::af::shared< std::vector<long int> >     map_asu_to_obs_;

      sgtbx::space_group                                 space_group_;
      bool                                            anomalous_flag_;
      cctbx::uctbx::unit_cell                              unit_cell_;
  };

}}}

#endif // CCTBX_XRAY_GROUPED_DATA
