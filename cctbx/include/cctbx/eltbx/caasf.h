#ifndef CCTBX_ELTBX_CAASF_H
#define CCTBX_ELTBX_CAASF_H

#include <cctbx/eltbx/basic.h>
#include <scitbx/array_family/small.h>
#include <cctbx/import_scitbx_af.h>
#include <boost/config.hpp>
#include <cstddef>
#include <cmath>
#include <ctype.h>

namespace cctbx { namespace eltbx { namespace caasf {

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
  /*! Currently used to work with coefficients from the International
      Tables Volume C (1992) (template parameter N = 4) and Waasmaier &
      Kirfel (1995), Acta Cryst. A51, 416-431 (template parameter N = 5).
   */
  template <std::size_t N>
  class base
  {
    public:
      //! Enables compile-time array allocation in other procedures.
      BOOST_STATIC_CONSTANT(std::size_t, n_plus_1 = N + 1);

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

      //! Label of table. Currently either "IT1992" or "WK1995".
      const char* table() const { return table_; }

      //! Scattering factor label. E.g. "Si4+".
      const char* label() const { return entry_->label; }

      //! Number of a and b coefficients. Currently either 4 or 5.
      /*! Note that the total number of coefficients is 2*n_ab()+1.
          <p>
          Not available in Python.
       */
      static std::size_t n_ab() { return N; }

      //! Coefficient a(i), with 0 <= i < n_ab().
      /*! No range checking (for runtime efficiency).
          <p>
          Python: a() returns a tuple of N floats.
       */
      float a(int i) const { return entry_->a[i]; }

      //! Coefficient b(i), with 0 <= i < n_ab().
      /*! No range checking (for runtime efficiency).
          <p>
          Python: b() returns a tuple of N floats.
       */
      float b(int i) const { return entry_->b[i]; }

      //! Coefficient c.
      float c() const { return entry_->c; }

#     include <cctbx/eltbx/caasf_at_d_star_sq.h>

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
    std::string work_label = basic::strip_label(label, exact);
    int m = 0;
    for (const detail::raw_table_entry<N>* entry = table_raw;
         entry->label;
         entry++) {
      if (std::string(entry->label) == "const" && work_label != "CONST") {
        continue;
      }
      if (std::string(entry->label) == "custom" && work_label != "CUSTOM") {
        continue;
      }
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
      throw error("Unknown scattering factor label.");
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
      If possible, the class caasf::wk1995 should be used
      instead of this class.
      The coefficients of Waasmaier & Kirfel give more precise
      approximations than the coefficients of Volume C
      and some errors are corrected.
      <p>
      See also:
        caasf::it1992_iterator,
        caasf::wk1995,
        caasf::base
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
          class caasf::base.
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
        caasf::wk1995_iterator,
        caasf::it1992,
        caasf::base
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

#ifndef CCTBX_ELTBX_CAASF_CUSTOM_MAX_N_AB
#define CCTBX_ELTBX_CAASF_CUSTOM_MAX_N_AB 5
#endif

  class custom
  {
    public:
      //! Maximum number of a,b pairs.
      BOOST_STATIC_CONSTANT(std::size_t,
        max_n_ab=CCTBX_ELTBX_CAASF_CUSTOM_MAX_N_AB);

      //! Default constructor. Some data members are not initialized!
      custom() {}

      //! Initialization with label and coefficients.
      custom(
        std::string const& label,
        af::small<float, max_n_ab> const& a,
        af::small<float, max_n_ab> const& b,
        float c)
      :
        label_(label),
        a_(a),
        b_(b),
        c_(c)
      {
        CCTBX_ASSERT(a_.size() == b_.size());
      }

      //! Scattering factor label as passed to the constructor.
      std::string const&
      label() const { return label_; }

      //! Number of a and b coefficients as passed to the constructor.
      std::size_t
      n_ab() const { return a_.size(); }

      //! Array of coefficients a.
      af::small<float, max_n_ab> const&
      a() const { return a_; }

      //! Array of coefficients b.
      af::small<float, max_n_ab> const&
      b() const { return b_; }

      //! Coefficient a(i), with 0 <= i < n_ab().
      /*! No range checking (for runtime efficiency).
          <p>
          Not available in Python.
       */
      float
      a(int i) const { return a_[i]; }

      //! Coefficient b(i), with 0 <= i < n_ab().
      /*! No range checking (for runtime efficiency).
          <p>
          Not available in Python.
       */
      float
      b(int i) const { return b_[i]; }

      //! Coefficient c.
      float
      c() const { return c_; }

#     include <cctbx/eltbx/caasf_at_d_star_sq.h>

    private:
      std::string label_;
      af::small<float, max_n_ab> a_;
      af::small<float, max_n_ab> b_;
      float c_;
  };

}}} // cctbx::eltbx::caasf

#endif // CCTBX_ELTBX_CAASF_H
