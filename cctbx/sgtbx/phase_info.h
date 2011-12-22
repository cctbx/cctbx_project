#ifndef CCTBX_SGTBX_PHASE_INFO_H
#define CCTBX_SGTBX_PHASE_INFO_H

#include <cctbx/miller.h>

#if defined(__llvm__) && defined(__APPLE_CC__) && __APPLE_CC__ <= 5658
#define CCTBX_SGTBX_PHASE_INFO_APPLE_LLVM2335_WORKAROUND
#endif

namespace cctbx { namespace sgtbx {

  class space_group; // forward declaration

  /*! \brief Handling of phase restrictions and optional evaluation
      of conditions for systematically absent reflections.
   */
  /*! A reflection with the Miller index h is "centric" if
      there is a symmetry operation with rotation part r such that
      h*r = -h. The phase of a centric reflection is restricted to
      two phase angels (modulo pi).
      <p>
      A reflection with the Miller index h is "systematically absent"
      if there is a symmetry operation with the rotation part r and
      the translation part t such that h*r = h and h*t != 0 mod 1.
   */
  class phase_info
  {
    public:
      //! Default constructor. Some data members are not initialized!
      phase_info() {}

      //! Determines the phase restriction for a given Miller index.
      /*! If no_test_sys_absent = false, it is also tested if a
          reflection with the given Miller index is systematically
          absent. If no_test_sys_absent = true, a faster algorithm
          is used that only determines the phase restriction. In the
          latter case the is_sys_absent() member function must not
          be used (an exception is thrown otherwise).
       */
      phase_info(
        sgtbx::space_group const& space_group,
        miller::index<> const& miller_index,
        bool no_test_sys_absent=false);

      /*! \brief Initialization with known product h*t and given
          translation denominator.
       */
      /*! sys_abs_was_tested indicates if it is known if a reflection
          with the given Miller index is systematically absent.
          The is_sys_absent() member function must not be used
          if sys_abs_was_tested = false (an exception is thrown
          otherwise).
          <p>
          Not available in Python.
       */
      phase_info(int ht, int t_den, bool sys_abs_was_tested)
      :
        ht_(ht), t_den_(t_den), sys_abs_was_tested_(sys_abs_was_tested)
      {}

      //! Tests if is_sys_absent() can be used.
      bool
      sys_abs_was_tested() const { return sys_abs_was_tested_; }

      //! Tests for systematically absent reflection.
      /*! Use sys_abs_was_tested() to determine if this test can be used.
          If the result of sys_abs_was_tested() = false, an exception
          is thrown if this test is used.
       */
      bool
      is_sys_absent() const
      {
        CCTBX_ASSERT(sys_abs_was_tested_);
        return ht_ == -2;
      }

      //! Tests if there actually is a phase restriction.
      /*! See class details.
       */
      bool
      is_centric() const { return ht_ >= 0; }

      //! Phase shift h*t (mod 1) corresponding to h*r = -h.
      /*! Low-level information for computing the restricted phases.
          ht() is multiplied by a factor t_den() in order to obtain
          an integer value.
          <p>
          See also: ht_angle()
       */
      int
      ht() const { return ht_; }

      //! Translation part denominator.
      /*! This is the factor by which ht() is multiplied.
       */
      int
      t_den() const { return t_den_; }

      //! Phase restriction in radians or degrees.
      /*! The return value is -1 if the phase is not restricted,
          and >= 0 and < pi or 180 otherwise.
       */
      double
      ht_angle(bool deg=false) const
      {
        if (!is_centric()) return -1.;
        return (pi_unit(deg) * ht_) / t_den_;
      }

      //! Tests if phase phi is compatible with restriction.
      /*! The tolerance compensates for rounding errors.
       */
      bool
      is_valid_phase(
        double phi,
        bool deg=false,
        double tolerance=1e-5) const;

      //! Nearest valid phase.
      /*! For acentric reflections equivalent to the input phase phi.
          For centric reflections, the restricted phase which is
          nearest to the input phase phi.
       */
      double
      nearest_valid_phase(double phi, bool deg=false) const;

      //! Projection of structure factor on vector with restricted phase.
      /*! For acentric reflections equivalent to input structure factor.
       */
      template <typename FloatType>
      std::complex<FloatType>
      valid_structure_factor(std::complex<FloatType> const& f) const
      {
        if (!is_centric()) return f;
        double a = ht_angle();
        double c = std::cos(a);
        double s = std::sin(a);
        double x = f.real() * c + f.imag() * s;
        return std::complex<FloatType>(x * c, x * s);
      }

    private:
      int ht_;
      int t_den_;
      bool sys_abs_was_tested_;

      double pi_unit(bool deg) const
      {
        if (deg) return 180.;
        return scitbx::constants::pi;
      }
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_PHASE_INFO_H
