#ifndef CCTBX_MILLER_SYM_EQUIV_H
#define CCTBX_MILLER_SYM_EQUIV_H

#include <cctbx/hendrickson_lattman.h>
#include <cctbx/sgtbx/space_group.h>

namespace cctbx { namespace miller {

  //! Class for symmetrically equivalent Miller indices.
  class sym_equiv_index
  {
    public:
      //! Default constructor. Some data members are not initialized!
      sym_equiv_index() {}

      //! Constructor.
      /*! hr is the product of the input index h and the
          rotation part r of a symmetry operation.
          <p>
          ht is the product of the input index h and the
          translation part of the symmetry operation,
          multiplied by the translation-part denominator t_den
          in order to obtain an integer value for ht.
          <p>
          friedel_flag indicates if Friedel's law was applied
          to arrive at h().
          <p>
          Not available in Python.
       */
      sym_equiv_index(index<> const& hr, int ht, int t_den, bool friedel_flag)
      : hr_(hr), ht_(ht) , t_den_(t_den), friedel_flag_(friedel_flag)
      {}

      //! The symmetrically equivalent index.
      index<>
      h() const
      {
        if (friedel_flag_) return -hr_;
        return hr_;
      }

      //! Product of Miller index and rotation part of symmetry operation.
      index<> const&
      hr() const { return hr_; }

      //! Product of Miller index and translation part of symmetry operation.
      int
      ht() const { return ht_; }

      //! Translation-part denominator.
      /*! This is the factor by which ht() is multiplied.
       */
      int
      t_den() const { return t_den_; }

      //! Phase shift h*t in radians or degrees.
      double
      ht_angle(bool deg = false) const
      {
        if (deg) return ht_angle_(360.);
        return ht_angle_(scitbx::constants::two_pi);
      }

      //! Flag for application of Friedel's law.
      /*! For centric reflections, friedel_flag() is always false.
       */
      bool
      friedel_flag() const { return friedel_flag_; }

      /*! \brief Returns new sym_equiv_index with flipped
          friedel_flag() if i_mate != 0.
       */
      sym_equiv_index
      mate(std::size_t i_mate=1) const
      {
        if (i_mate) return sym_equiv_index(hr_, ht_, t_den_, !friedel_flag_);
        return *this;
      }

      /*! \brief Phase for equivalent index, given phase for
          input Miller index.
       */
      /*! Formula used:<br>
          deg == 0: phi_eq = phi_in - (2 * pi * ht()) / t_den();<br>
          deg != 0: phi_eq = phi_in - (360 * ht()) / t_den();<br>
          if (friedel_flag()) phi_eq = -phi_eq;
       */
      template <class FloatType>
      FloatType
      phase_eq(FloatType const& phi_in, bool deg=false) const
      {
        FloatType phi_eq = phi_in - ht_angle(deg);
        if (friedel_flag_) return -phi_eq;
        return phi_eq;
      }

      /*! \brief Phase for input index, given phase for
          equivalent Miller index.
       */
      /*! Formula used:<br>
          if (friedel_flag()) phi_eq = -phi_eq;<br>
          deg == 0: phi_in = phi_eq + (2 * pi * ht()) / t_den();<br>
          deg != 0: phi_in = phi_eq + (360 * ht()) / t_den();
       */
      template <class FloatType>
      FloatType
      phase_in(FloatType phi_eq, bool deg=false) const
      {
        if (friedel_flag_) phi_eq = -phi_eq;
        return phi_eq + ht_angle(deg);
      }

      /*! \brief Complex value for equivalent index, given complex
          value for input index.
       */
      /*! Formula used:<br>
          f_eq = f_in * exp(-2 * pi * j * ht() / t_den());<br>
          where j is the imaginary number.<br>
          if (friedel_flag()) f_eq = conj(f_eq);
       */
      template <class FloatType>
      std::complex<FloatType>
      complex_eq(std::complex<FloatType> const& f_in) const
      {
        std::complex<FloatType>
        f_eq = f_in * std::polar(FloatType(1), -ht_angle());
        if (friedel_flag_) return std::conj(f_eq);
        return f_eq;
      }

      /*! \brief Complex value for input index, given complex
          value for equivalent index.
       */
      /*! Formula used:<br>
          if (friedel_flag()) f_eq = conj(f_eq);<br>
          f_in = f_eq * exp(2 * pi * j * ht() / t_den());<br>
          where j is the imaginary number.
       */
      template <class FloatType>
      std::complex<FloatType>
      complex_in(std::complex<FloatType> const& f_eq) const
      {
        std::complex<FloatType> shift = std::polar(FloatType(1), ht_angle());
        if (friedel_flag_) return std::conj(f_eq) * shift;
        return f_eq * shift;
      }

      /*! \brief Hendrickson-Lattman coefficients for equivalent index,
          given coefficients for input index.
       */
      template <class FloatType>
      hendrickson_lattman<FloatType>
      hendrickson_lattman_eq(hendrickson_lattman<FloatType> const& hl_in) const
      {
        hendrickson_lattman<FloatType> hl_eq = hl_in.shift_phase(-ht_angle());
        if (friedel_flag_) return hl_eq.conj();
        return hl_eq;
      }

      /*! \brief Hendrickson-Lattman coefficients for input index,
          given coefficients for equivalent index.
       */
      template <class FloatType>
      hendrickson_lattman<FloatType>
      hendrickson_lattman_in(hendrickson_lattman<FloatType> hl_eq) const
      {
        if (friedel_flag_) hl_eq = hl_eq.conj();
        return hl_eq.shift_phase(ht_angle());
      }

    protected:
      double ht_angle_(double period) const
      {
        return (period * ht_) / t_den_;
      }

      index<> hr_;
      int   ht_;
      int   t_den_;
      bool  friedel_flag_;
  };

  //! class for the handling of symmetrically equivalent Miller indices.
  class sym_equiv_indices
  {
    public:
      //! Default constructor. Some data members are not initialized!
      sym_equiv_indices() {}

      //! Computation of the symmetrically equivalent Miller indices.
      /*! The Miller index passed to the constructor is referred
          to as the "input Miller index."
          <p>
          The conditions for systematically absent reflections are
          NOT tested.
       */
      sym_equiv_indices(sgtbx::space_group const& sg, index<> const& h_in);

      //! Phase restriction (if any) for the input Miller index.
      sgtbx::phase_info
      phase_restriction() const
      {
        return sgtbx::phase_info(ht_restriction_, t_den_, false);
      }

      //! Tests if a reflection with the input Miller index is centric.
      /*! A reflection with the Miller index h is "centric" if
          there is a symmetry operation with rotation part r such
          that h*r = -h.<br>
          See also: phase_info
       */
      bool
      is_centric() const { return ht_restriction_ >= 0; }

      //! Raw array of symmetrically equivalent Miller indices.
      /*! Note that indices().size() is not in general equal to
          multiplicity().
       */
      af::shared<sym_equiv_index> const&
      indices() const { return indices_; }

      //! Multiplicity of the input Miller index.
      /*! For acentric reflections and in the presence of Friedel symmetry
          (anomalous_flag == false), multiplicity() is twice the number
          of symmetrically equivalent Miller indices().size().
          <p>
          For centric reflections or in the absence of Friedel symmetry
          (anomalous_flag == true), the multiplicity is equal to the
          number of symmetrically equivalent Miller indices().size().
       */
      int
      multiplicity(bool anomalous_flag) const
      {
        if (!anomalous_flag && !is_centric()) return 2 * indices_.size();
        return indices_.size();
      }

      //! Factor for loop over Friedel mates.
      /*! Useful for looping over all symmetrically equivalent reflections
          including Friedel mates (if anomalous_flag == false):
          <p>
          f_mates() == multiplicity(anomalous_flag) / indices().size()
          <p>
          See also: operator()
       */
      int
      f_mates(bool anomalous_flag) const
      {
        if (!anomalous_flag && !is_centric()) return 2;
        return 1;
      }

      //! Determines "epsilon" for the given Miller index.
      /*! The factor epsilon counts the number of times a Miller
          index h is mapped onto itself by symmetry. This factor
          is used for "statistical averaging" and in direct methods
          formulae.
          <p>
          Note that epsilon is directly related to the number
          of symmetrically equivalent indices().size():
          <p>
          epsilon == sgtbx::space_group::order_p() / indices().size()
       */
      int
      epsilon() const { return order_p_ / indices_.size(); }

      //! Medium-level access to the symmetrically equivalent Miller indices.
      /*! Intended use:<pre>
       sgtbx::space_group sg = ... // define space group
       miller::index<> h = ... // define input Miller index.
       bool anomalous_flag = ...
       miller::sym_equiv_indices h_eq(sg, h);
       for(std::size_t i_indices=0;i_indices<h_eq.indices().size();i_indices++)
         for(std::size_t i_mate=0;i_mate<h_eq.f_mates(anomalous_flag);i_mate++)
           miller::index<> k = h_eq(i_mate, i_indices).h();
          </pre>
          Note that it is possible and often more convenient to have a
          one-deep loop with multiplicity() iterations.
          <p>
          See also: operator()(std::size_t i)
          <p>
          Not available in Python.
       */
      sym_equiv_index
      operator()(std::size_t i_mate, std::size_t i_indices) const;

      //! High-level access to the symmetrically equivalent Miller indices.
      /*! Intended use:<pre>
          sgtbx::space_group sg = ... // define space group
          miller::index<> h = ... // define input Miller index.
          bool anomalous_flag = ...
          miller::sym_equiv_indices h_eq(sg, h);
          for(std::size_t i=0;i<h_eq.multiplicity(anomalous_flag);i++)
            miller::index<> k = h_eq(i).h();
          </pre>
       */
      sym_equiv_index
      operator()(std::size_t i) const;

      //! Tests if the phase phi is compatible with phase_restriction().
      /*! tolerance compensates for rounding errors.
          <p>
          See also: phase_info
       */
      bool
      is_valid_phase(double phi,
                     bool deg=false,
                     double tolerance=1e-5) const
      {
        return phase_restriction().is_valid_phase(phi, deg, tolerance);
      }

      //! Equivalent indices for a P1 listing.
      /*! If anomalous_flag == true the number of
          indices in P1 is always equal to indices().size().
          If anomalous_flag == false the number of indices in P1 is
          indices().size() for non-centric reflections and
          indices().size()/2 for centric reflections.
       */
      af::shared<sym_equiv_index>
      p1_listing(bool anomalous_flag) const;

    protected:
      void
      add(sym_equiv_index const& eq);

      struct index_mate_indices_decomposition
      {
        index_mate_indices_decomposition(
          std::size_t i_mate_,
          std::size_t i_indices_)
        :
          i_mate(i_mate_), i_indices(i_indices_)
        {}

        std::size_t i_mate, i_indices;
      };

      index_mate_indices_decomposition
      decompose_index_mate_indices(std::size_t i) const;

      int t_den_;
      std::size_t order_p_;
      int ht_restriction_;
      af::shared<sym_equiv_index> indices_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_SYM_EQUIV_H
