#ifndef SCITBX_MATH_RESAMPLE_H
#define SCITBX_MATH_RESAMPLE_H

#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/constants.h>
#include <scitbx/random.h>

#include <cmath>
#include <cstdio>

namespace scitbx { namespace math {
namespace resample{

  //! resampling class
  template <typename FloatType>
  class non_parametric_bootstrap
  {
    public:
    /*! Default constructor */
    non_parametric_bootstrap() {}
    /*! Constructor that inits all members to nought */
    non_parametric_bootstrap(
       scitbx::af::const_ref<FloatType> const& observations,
       long const& seed)
    :
    generator_(seed)
    {
      for (unsigned ii=0;ii<observations.size();ii++){
        observations_.push_back( observations[ii] );
      }
    }

    scitbx::af::shared<FloatType>
    draw(std::size_t n_new_obs)
    {
      // We have to draw n_new_obs numbers between
      // 0 and observations_.size() (0 inclusive, the latter one exclusive)
      scitbx::af::shared< std::size_t > random_integers(n_new_obs,0);
      random_integers = generator_.random_size_t(
        n_new_obs,
        observations_.size() );

      // The indices have been drawn, it is now straightforward
      // to make a new sample
      scitbx::af::shared<FloatType> new_data;
      for (unsigned ii=0;ii<n_new_obs;ii++){
        new_data.push_back( observations_[ random_integers[ ii] ] );
      }
      return( new_data );
    }

    scitbx::af::shared<FloatType>
    draw_from_random_jack_knifed_sample(std::size_t n_new_obs,
                                       std::size_t jack)
    {
      SCITBX_ASSERT( jack < observations_.size() );
      // I assume one knows what one is doing, this is just an safety check
      // first draw  observations_.size()-jack items
      scitbx::af::shared<std::size_t> jacked, sample;

      jacked = generator_.random_size_t(
        observations_.size()-jack,
        observations_.size() );
      sample = generator_.random_size_t(
        n_new_obs,
        observations_.size()-jack );

      scitbx::af::shared<FloatType> new_data;
      for (unsigned ii=0;ii<n_new_obs;ii++){
        new_data.push_back( observations_[ jacked[ sample[ ii] ] ] );
      }
      return( new_data );
    }

    protected:
    scitbx::af::shared<FloatType> observations_;
    scitbx::random::mersenne_twister generator_;
  };

  //------------------------------------------------------------

  //! resampling class
  template <typename FloatType = std::size_t>
  class non_parametric_bootstrap_as_int
  {
    public:
    /*! Default constructor */
    non_parametric_bootstrap_as_int() {}
    /*! Constructor that inits all members to nought */
    non_parametric_bootstrap_as_int(
       scitbx::af::const_ref<FloatType> const& observations,
       long const& seed)
    :
    generator_(seed)
    {
      for (unsigned ii=0;ii<observations.size();ii++){
        observations_.push_back( observations[ii] );
      }
    }

    scitbx::af::shared<FloatType>
    draw(std::size_t n_new_obs)
    {
      // We have to draw n_new_obs numbers between
      // 0 and observations_.size() (0 inclusive, the latter one exclusive)
      scitbx::af::shared< std::size_t > random_integers(n_new_obs,0);
      random_integers = generator_.random_size_t(
        n_new_obs,
        observations_.size() );

      // The indices have been drawn, it is now straightforward
      // to make a new sample
      scitbx::af::shared<FloatType> new_data;
      for (unsigned ii=0;ii<n_new_obs;ii++){
        new_data.push_back( observations_[ random_integers[ ii] ] );
      }
      return( new_data );
    }

    scitbx::af::shared<FloatType>
    draw_from_random_jack_knifed_sample(std::size_t n_new_obs,
                                       std::size_t jack)
    {
      SCITBX_ASSERT( jack < observations_.size() );
      // I assume one knows what one is doing, this is just an safety check
      // first draw  observations_.size()-jack items
      scitbx::af::shared<std::size_t> jacked, sample;

      jacked = generator_.random_size_t(
        observations_.size()-jack,
        observations_.size() );
      sample = generator_.random_size_t(
        n_new_obs,
        observations_.size()-jack );

      scitbx::af::shared<FloatType> new_data;
      for (unsigned ii=0;ii<n_new_obs;ii++){
        new_data.push_back( observations_[ jacked[ sample[ ii] ] ] );
      }
      return( new_data );
    }

    protected:
    scitbx::af::shared<FloatType> observations_;
    scitbx::random::mersenne_twister generator_;
  };





  //------------------------------------------------------------

  //! resampling class
  template <typename FloatType = double>
  class smooth_bootstrap
  {
    public:
    /*! Default constructor */
    smooth_bootstrap() {}
    /*! Constructor that inits all members to nought */
    smooth_bootstrap(
       scitbx::af::const_ref<FloatType> const& observations,
       long const& seed)
    :
    generator_(seed),
    observed_mean_(0),
    observed_std_(0)
    {
      for (unsigned ii=0;ii<observations.size();ii++){
        observations_.push_back( observations[ii] );
        observed_mean_+=observations[ii];
        observed_std_+=observations[ii]*observations[ii];
      }
      observed_mean_ /= FloatType( observations.size() );
      observed_std_ /= FloatType( observations.size()-1 );
      observed_std_ -= (observed_mean_ * observed_mean_);
      observed_std_ = std::sqrt( observed_std_ );
    }

    scitbx::af::shared<FloatType>
    draw(std::size_t n_new_obs)
    {
      // We have to draw n_new_obs numbers between
      // 0 and observations_.size() (0 inclusive, the latter one exclusive)
      scitbx::af::shared< std::size_t > random_integers(n_new_obs,0);
      random_integers = generator_.random_size_t(
        n_new_obs,
        observations_.size() );

      // The indices have been drawn, it is now straightforward
      // to make a new sample
      scitbx::af::shared<FloatType> new_data;
      FloatType add;
      for (unsigned ii=0;ii<n_new_obs;ii++){
        add=gauss( observed_std_/std::sqrt( FloatType(n_new_obs) ) );
        new_data.push_back( observations_[ random_integers[ ii] ] + add );
      }
      return( new_data );
    }

    scitbx::af::shared<FloatType>
    draw_from_random_jack_knifed_sample(std::size_t n_new_obs,
                                       std::size_t jack)
    {
      SCITBX_ASSERT( jack < observations_.size() );
      // I assume one knows what one is doing, this is just an safety check
      // first draw  observations_.size()-jack items
      scitbx::af::shared<std::size_t> jacked, sample;

      jacked = generator_.random_size_t(
        observations_.size()-jack,
        observations_.size() );
      sample = generator_.random_size_t(
        n_new_obs,
        observations_.size()-jack );

      scitbx::af::shared<FloatType> new_data;
      FloatType add;
      for (unsigned ii=0;ii<n_new_obs;ii++){
        add = gauss( observed_std_/std::sqrt( FloatType(n_new_obs) ) );
        new_data.push_back( observations_[ jacked[ sample[ ii] ] ] + add);
      }
      return( new_data );
    }

    protected:

    FloatType
    gauss(FloatType h)
    {
      // This is a basic
      FloatType u1,u2,result;
      u1 = generator_.random_double();
      u2 = generator_.random_double();
      result = std::sqrt( -2.0*std::log(u1) )*
        std::cos( 2.0* scitbx::constants::pi*u2 ) * h;
      return (result);

    }



    scitbx::af::shared<FloatType> observations_;
    scitbx::random::mersenne_twister generator_;
    FloatType observed_mean_;
    FloatType observed_std_;

  };

}}} // namespace scitbx::math::resample

#endif // SCITBX_MATH_RESAMPLE_H
