#ifndef CCTBX_ELTBX_XRAY_SCATTERING_H
#define CCTBX_ELTBX_XRAY_SCATTERING_H

#include <cctbx/eltbx/basic.h>
#include <scitbx/math/gaussian/sum.h>
#include <cctbx/import_scitbx_af.h>
#include <boost/config.hpp>
#include <ctype.h>

namespace cctbx { namespace eltbx { namespace xray_scattering {

  class gaussian : public scitbx::math::gaussian::sum<double>
  {
    public:
      typedef scitbx::math::gaussian::sum<double> base_t;

      //! Default constructor. Some data members are not initialized!
      gaussian() {}

      //! Initialization given an instance of the base type.
      gaussian(base_t const& gaussian_sum)
      :
        base_t(gaussian_sum)
      {}

      //! Initialization of the constant.
      explicit
      gaussian(double c, bool use_c=true)
      :
        base_t(c, use_c)
      {}

      //! Initialization of the terms and optionally the constant.
      /*! If c is different from zero use_c will automatically be
          set to true.
       */
      gaussian(
        af::small<double, base_t::max_n_terms> const& a,
        af::small<double, base_t::max_n_terms> const& b,
        double c=0,
        bool use_c=false)
      :
        base_t(a, b, c, use_c)
      {}

      /*! \brief Sum of Gaussian terms at the point stol
          (sin-theta-over-lambda), given stol^2.
       */
      /*! See also: at_stol(), at_d_star(), at_d_star_sq(),
                    uctbx::unit_cell::stol()
       */
      double
      at_stol_sq(double stol_sq) const
      {
        return at_x_sq(stol_sq);
      }

      /*! \brief Sum of Gaussian terms at the point stol
          (sin-theta-over-lambda).
       */
      /*! See also: at_stol_sq(), at_d_star(), at_d_star_sq(),
       */
      double
      at_stol(double stol) const
      {
        return at_x_sq(stol * stol);
      }

      /*! \brief Sum of Gaussian terms at the point d_star
          (1/d), given d_star^2.
       */
      /*! See also: at_stol_sq(), at_stol(), at_d_star(),
                    uctbx::unit_cell::d_star_sq()
       */
      double
      at_d_star_sq(double d_star_sq) const
      {
        return at_x_sq(d_star_sq / 4);
      }

      /*! \brief Sum of Gaussian terms at the point d_star
          (1/d).
       */
      /*! See also: at_stol_sq(), at_stol(), at_d_star_sq(),
                    uctbx::unit_cell::d_star_sq()
       */
      double
      at_d_star(double d_star) const
      {
        return at_x_sq(d_star * d_star / 4);
      }
  };

  namespace detail {

    static const float undefined = -99999.;

    template <std::size_t N>
    struct raw_table_entry
    {
      const char* label;
      float a[N], b[N], c;
    };

  }

  //! Coefficients for the Analytical Approximation to the Scattering Factor.
  /*! Used to work with coefficients from the International
      Tables Volume C (1992) (template parameter N = 4) and Waasmaier &
      Kirfel (1995), Acta Cryst. A51, 416-431 (template parameter N = 5).
   */
  template <std::size_t N>
  class base
  {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      base() : entry_(0), table_(0) {}

      //! Constructor. For internal use only.
      base(const detail::raw_table_entry<N>* table_raw,
           const char* table,
           std::string const& label,
           bool exact);

      //! Tests if the instance is constructed properly.
      /*! Shorthand for: label() != 0
          <p>
          Not available in Python.
       */
      bool
      is_valid() const { return entry_->label != 0; }

      //! Label of table. Either "IT1992" or "WK1995".
      const char* table() const { return table_; }

      //! Scattering factor label. E.g. "Si4+".
      const char* label() const { return entry_->label; }

      //! Fetches the Gaussian terms from the static table.
      gaussian
      fetch() const
      {
        typedef af::small<double, gaussian::max_n_terms> small_t;
        return gaussian(
          small_t(entry_->a, entry_->a+N),
          small_t(entry_->b, entry_->b+N),
          entry_->c, true);
      }

    protected:
      const char *table_;
      const detail::raw_table_entry<N>* entry_;
      void next_entry();
  };

  // non-inline constructor
  template <std::size_t N>
  base<N>::base(const detail::raw_table_entry<N>* table_raw,
                const char* table,
                std::string const& label,
                bool exact)
  :
    entry_(0), table_(table)
  {
    if (label == "const") {
      throw error("Reserved scattering type label: const");
    }
    std::string work_label = basic::strip_label(label, exact);
    int m = 0;
    for (const detail::raw_table_entry<N>* entry = table_raw;
         entry->label;
         entry++) {
      int i = basic::match_labels(work_label, entry->label);
      if (i < 0) {
        if (entry->c == detail::undefined) {
          throw error("Analytical approximation for given label not known.");
        }
        entry_ = entry;
        return;
      }
      if (i > m && !isdigit(entry->label[i - 1])) {
        m = i;
        entry_ = entry;
      }
    }
    if (exact || !entry_) {
      throw error("Unknown scattering type label: " + std::string(label));
    }
  }

  // non-inline member function
  template <std::size_t N>
  void
  base<N>::next_entry()
  {
    if (is_valid()) {
      for(;;) {
        entry_++;
        if (entry_->c != detail::undefined) break;
      }
    }
  }

  //! Coefficients for the Analytical Approximation to the Scattering Factor.
  /*! Coefficients from the International Tables for Crystallography,
      Volume C, Mathematical, Physical and Chemical Tables,
      Edited by A.J.C. Wilson, Kluwer Academic Publishers,
      Dordrecht/Boston/London, 1992.
      <p>
      Table 6.1.1.4 (pp. 500-502):
      Coefficients for analytical approximation to the scattering factors
      of Tables 6.1.1.1 and 6.1.1.3.
      <p>
      Table 6.1.1.4 is a reprint of Table 2.2B, pp. 99-101,
      International Tables for X-ray Crystallography, Volume IV,
      The Kynoch Press: Birmingham, England, 1974.
      There is just one difference for <tt>Tl3+</tt>:<pre>
                IT Vol IV 1974: a2 = 18.3841
                IT Vol C  1992: a2 = 18.3481</pre>
      <p>
      The value from IT Vol IV is presumed to be the correct one
      and used here.
      <p>
      If possible, the class xray_scattering::wk1995 should be used
      instead of this class.
      The coefficients of Waasmaier & Kirfel give more precise
      approximations than the coefficients of Volume C
      and some errors are corrected.
      <p>
      See also:
        xray_scattering::it1992_iterator,
        xray_scattering::wk1995,
        xray_scattering::base
   */
  class it1992: public base<4>
  {
    public:
      //! Facilitates meta-programming.
      typedef base<4> base_type;

      //! Default constructor. Calling certain methods may cause crashes!
      it1992() {}

      //! Looks up coefficients for the given scattering factor label.
      /*! If exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          E.g., "SI4" will be matched with "Si".<br>
          "Si4+" and "Si+4" will be matched with "Si4+".<br>
          See also: eltbx::basic::strip_label()
          <p>
          Note that the other methods of this class are inherited from
          class xray_scattering::base.
       */
      it1992(std::string const& label, bool exact=false);

    protected:
      friend class it1992_iterator;
  };

  /*! \brief Iterator over table of Coefficients for the Analytical
      Approximation to the Scattering Factor, International Tables 1992.
   */
  class it1992_iterator
  {
    public:
      //! Initialization of the iterator.
      it1992_iterator();

      //! Retrieves the next entry from the internal table.
      /*! Use it1992::is_valid() to detect end-of-iteration.
       */
      it1992
      next();

    private:
      it1992 current_;
  };

  //! Coefficients for the Analytical Approximation to the Scattering Factor.
  /*! Coefficients from Waasmaier & Kirfel (1995), Acta Cryst. A51, 416-431.
      <p>
      The original data can be found at:
      <pre>
      ftp://wrzx02.rz.uni-wuerzburg.de/pub/local/Crystallography/sfac.dat
      </pre>
      File picked up Jul  4, 1995.
      File verified  Sep 12, 2001.
      <p>
      Header of "sfac.dat":
      <pre>
        - Fitparameters of all atoms/ions (with the excepetion of O1-)
          from publication "New Analytical Scattering Factor Functions
          for Free Atoms and Ions", D. Waasmaier & A. Kirfel, Acta Cryst.
          A 1995, in press)

        - Fit for O1- based on the tabulated values of Table 2 (D.REZ,
          P.REZ & I.GRANT, Acta Cryst. (1994), A50, 481-497).

        - Fits for all other atoms/ions based on the tabulated values
          of Table 6.1.1.1 (atoms) and Table 6.1.1.3 (ions)
          (International Tables for Crystallography, Vol. C, 1992).
      </pre>
      <p>
      See also:
        xray_scattering::wk1995_iterator,
        xray_scattering::it1992,
        xray_scattering::base
   */
  class wk1995: public base<5>
  {
    public:
      //! Facilitates meta-programming.
      typedef base<5> base_type;

      //! Default constructor. Calling certain methods may cause crashes!
      wk1995() {}

      //! Looks up coefficients for the given scattering factor label.
      /*! If exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          E.g., "SI4" will be matched with "Si".<br>
          "Si4+" and "Si+4" will be matched with "Si4+".<br>
          See also: eltbx::basic::strip_label()
          <p>
          Note that the other methods of this class are inherited from
          class base.
       */
      wk1995(std::string const& label, bool exact=false);

    protected:
      friend class wk1995_iterator;
  };

  /*! \brief Iterator over table of Coefficients for the Analytical
      Approximation to the Scattering Factor, Waasmaier & Kirfel 1995.
   */
  class wk1995_iterator
  {
    public:
      //! Initialization of the iterator.
      wk1995_iterator();

      //! Retrieves the next entry from the internal table.
      /*! Use wk1995::is_valid() to detect end-of-iteration.
       */
      wk1995
      next();

    private:
      wk1995 current_;
  };

}}} // cctbx::eltbx::xray_scattering

#endif // CCTBX_ELTBX_XRAY_SCATTERING_H
