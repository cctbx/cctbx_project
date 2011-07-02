#ifndef CCTBX_HENDRICKSON_LATTMAN_H
#define CCTBX_HENDRICKSON_LATTMAN_H

#include <scitbx/array_family/tiny_plain.h>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/math/bessel.h>
#include <scitbx/math/atanh.h>
#include <scitbx/constants.h>
#include <boost/shared_array.hpp>
#include <complex>

namespace cctbx {

  //! Grouping of Hendrickson-Lattman coefficients.
  /*! Reference: Jan Drenth,
                 Principles of Protein X-Ray Crystallography,
                 Second edition, 1999, Chapter 14.
   */
  template<typename FloatType = double>
  class hendrickson_lattman : public af::tiny_plain<FloatType, 4>
  {
    public:
      typedef af::tiny_plain<FloatType, 4> base_type;

      //! Default constructor. The coefficients are not initialized!
      hendrickson_lattman() {}

      //! Initializtion from tiny array.
      hendrickson_lattman(base_type const& coeff)
      : base_type(coeff)
      {}

      //! Initializtion from plain pointer to coefficients.
      explicit
      hendrickson_lattman(const FloatType* coeff)
      {
        std::copy(coeff, coeff + 4, this->begin());
      }

      //! Initializtion given individual coefficients.
      hendrickson_lattman(
        FloatType const& a,
        FloatType const& b,
        FloatType const& c,
        FloatType const& d)
      {
        (*this)[0] = a;
        (*this)[1] = b;
        (*this)[2] = c;
        (*this)[3] = d;
      }

      /*! \brief Initialization given a phase integral (complex
          representation of centroid phase and figure of merit).
       */
      /*! The absolute value of the phase integral is truncated
          at max_figure_of_merit to avoid singularities.
          <p>
          See also: cctbx::miller::phase_integrator
       */
      hendrickson_lattman(
        bool centric_flag,
        std::complex<FloatType> const& phase_integral,
        FloatType const& max_figure_of_merit)
      {
        FloatType fom = std::min(std::abs(phase_integral),max_figure_of_merit);
        FloatType weight;
        if (centric_flag) {
          weight = scitbx::math::atanh(fom);
        }
        else {
          weight = scitbx::math::bessel::inverse_i1_over_i0(fom);
        }
        FloatType angle = std::arg(phase_integral);
        this->elems[0] = weight * std::cos(angle);
        this->elems[1] = weight * std::sin(angle);
        this->elems[2] = 0;
        this->elems[3] = 0;
      }

      //! Coefficients a,b,c,d as array.
      base_type const&
      coeff() const { return *this; }

      //! Coefficients a,b,c,d as array.
      base_type&
      coeff()       { return *this; }

      //! Individual coefficient a.
      FloatType const& a() const { return (*this)[0]; }

      //! Individual coefficient b.
      FloatType const& b() const { return (*this)[1]; }

      //! Individual coefficient c.
      FloatType const& c() const { return (*this)[2]; }

      //! Individual coefficient d.
      FloatType const& d() const { return (*this)[3]; }

      /*! \brief Coefficients for Friedel opposite (similar to conjugate
          complex of structure factor).
       */
      /*! Formula used: a, -b, c, -d
          <p>
          See also: cctbx::miller::sym_equiv_index
       */
      hendrickson_lattman
      conj() const
      {
        return hendrickson_lattman(a(), -b(), c(), -d());
      }

      //! Coefficients for symmetrically equivalent reflections.
      /*! The phase shift delta_phi must be given in radians.
          <p>
          See also: cctbx::miller::sym_equiv_index
       */
      hendrickson_lattman
      shift_phase(FloatType const& delta_phi) const
      {
        FloatType c1 = std::cos(delta_phi);
        FloatType s1 = std::sin(delta_phi);
        FloatType c2 = std::cos(2. * delta_phi);
        FloatType s2 = std::sin(2. * delta_phi);
        return hendrickson_lattman(
          a()*c1 - b()*s1,
          a()*s1 + b()*c1,
          c()*c2 - d()*s2,
          c()*s2 + d()*c2);
      }

      //! Phase combination.
      hendrickson_lattman
      operator+(hendrickson_lattman const& rhs) const
      {
        hendrickson_lattman result;
        for(unsigned i=0;i<4;i++) {
          result[i] = this->elems[i] + rhs[i];
        }
        return result;
      }

      //! In-place phase combination.
      hendrickson_lattman&
      operator+=(hendrickson_lattman const& rhs)
      {
        for(unsigned i=0;i<4;i++) {
          this->elems[i] += rhs[i];
        }
        return *this;
      }

      //! Division by integer (for averaging).
      hendrickson_lattman
      operator/(FloatType const& rhs) const
      {
        hendrickson_lattman result;
        for(unsigned i=0;i<4;i++) {
          result[i] = this->elems[i] / rhs;
        }
        return result;
      }

      //! Test for exact equality.
      bool
      operator==(hendrickson_lattman const& rhs) const
      {
        for(unsigned i=0;i<4;i++) {
          if (this->elems[i] != rhs[i]) return false;
        }
        return true;
      }

      //! Test for inequality.
      bool
      operator!=(hendrickson_lattman const& rhs) const
      {
        return !(*this == rhs);
      }

      struct phase_integration_cos_sin_table
      {
        unsigned n_steps;
        FloatType angular_step;
        boost::shared_array<af::tiny_plain<FloatType, 4> > data;

        phase_integration_cos_sin_table() {}

        phase_integration_cos_sin_table(
          unsigned n_steps_)
        :
          n_steps(n_steps_),
          angular_step(scitbx::constants::two_pi / n_steps),
          data(new af::tiny_plain<FloatType, 4>[n_steps])
        {
          af::tiny_plain<FloatType, 4>* d = data.get();
          for(unsigned i_step=0;i_step<n_steps_;i_step++) {
            FloatType angle = i_step * angular_step;
            *d++ = af::tiny_plain<FloatType, 4>(
              std::cos(angle),
              std::sin(angle),
              std::cos(angle+angle),
              std::sin(angle+angle));
          }
        }
      };
  };

} // namespace cctbx

#endif // CCTBX_HENDRICKSON_LATTMAN_H
