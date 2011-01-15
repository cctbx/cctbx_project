#ifndef CCTBX_ELTBX_XRAY_SCATTERING_H
#define CCTBX_ELTBX_XRAY_SCATTERING_H

#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/xray_scattering/gaussian.h>
#include <boost/optional.hpp>
#include <boost/config.hpp>
#include <stdexcept>
#include <ctype.h>

namespace cctbx { namespace eltbx { namespace xray_scattering {

  static const char *standard_labels[] = {
    "H", "D", "T",
    "Hsds",
    "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
    "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
    "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
    "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "H1-", "D1-", "T1-",
    "Li1+", "Be2+", "Cval", "O1-", "O2-", "F1-", "Na1+", "Mg2+",
    "Al3+", "Sival", "Si4+", "Cl1-", "K1+", "Ca2+", "Sc3+", "Ti2+",
    "Ti3+", "Ti4+", "V2+", "V3+", "V5+", "Cr2+", "Cr3+", "Mn2+",
    "Mn3+", "Mn4+", "Fe2+", "Fe3+", "Co2+", "Co3+", "Ni2+", "Ni3+",
    "Cu1+", "Cu2+", "Zn2+", "Ga3+", "Ge4+", "Br1-", "Rb1+", "Sr2+",
    "Y3+", "Zr4+", "Nb3+", "Nb5+", "Mo3+", "Mo5+", "Mo6+", "Ru3+",
    "Ru4+", "Rh3+", "Rh4+", "Pd2+", "Pd4+", "Ag1+", "Ag2+", "Cd2+",
    "In3+", "Sn2+", "Sn4+", "Sb3+", "Sb5+", "I1-", "Cs1+", "Ba2+",
    "La3+", "Ce3+", "Ce4+", "Pr3+", "Pr4+", "Nd3+", "Pm3+", "Sm3+",
    "Eu2+", "Eu3+", "Gd3+", "Tb3+", "Dy3+", "Ho3+", "Er3+", "Tm3+",
    "Yb2+", "Yb3+", "Lu3+", "Hf4+", "Ta5+", "W6+", "Os4+", "Ir3+",
    "Ir4+", "Pt2+", "Pt4+", "Au1+", "Au3+", "Hg1+", "Hg2+", "Tl1+",
    "Tl3+", "Pb2+", "Pb4+", "Bi3+", "Bi5+", "Ra2+", "Ac3+", "Th4+",
    "U3+", "U4+", "U6+", "Np3+", "Np4+", "Np6+", "Pu3+", "Pu4+",
    "Pu6+", 0
  };

  inline
  bool
  is_reserved_scattering_type_label(
    std::string const& label)
  {
    if (label == "const") return true;
    if (label == "unknown") return true;
    return false;
  }

  inline
  void
  throw_if_reserved_scattering_type_label(
    std::string const& label)
  {
    if (is_reserved_scattering_type_label(label)) {
      throw std::runtime_error(
        "Reserved scattering type label: \""+label+"\"");
    }
  }

  inline
  boost::optional<std::string>
  get_standard_label(
    std::string const& label,
    bool exact=false,
    bool optional=false)
  {
    if (label == "const") return boost::optional<std::string>(label);
    std::string work_label = basic::strip_label(label, exact);
    const char* result = 0;
    int m = 0;
    for (const char **std_lbl = standard_labels; *std_lbl; std_lbl++) {
      int i = basic::match_labels(work_label, *std_lbl);
      if (i < 0 /* exact match */) {
        return boost::optional<std::string>(*std_lbl);
      }
      if (i > m && !isdigit((*std_lbl)[i-1])) {
        m = i;
        result = *std_lbl;
      }
    }
    if (exact || result == 0) {
      if (optional) return boost::optional<std::string>();
      throw std::runtime_error(
        "Unknown scattering type label: \"" + label + "\"");
    }
    return boost::optional<std::string>(result);
  }

  inline
  std::string
  replace_hydrogen_isotype_labels(std::string standard_label)
  {
    if      (standard_label == "D"   || standard_label == "T"  ) return "H";
    else if (standard_label == "D1-" || standard_label == "T1-") return "H1-";
    return standard_label;
  }

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
      base() : table_(0), entry_(0) {}

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
    table_(table), entry_(0)
  {
    throw_if_reserved_scattering_type_label(label);
    std::string std_lbl = replace_hydrogen_isotype_labels(
      *get_standard_label(label, exact));
    for (const detail::raw_table_entry<N>* entry = table_raw;
         entry->label;
         entry++) {
      if (entry->label == std_lbl) {
        entry_ = entry;
        break;
      }
    }
    if (entry_ == 0) {
      throw std::runtime_error(
        "Unknown scattering type label: \"" + label + "\"");
    }
  }

  // non-inline member function
  template <std::size_t N>
  void
  base<N>::next_entry()
  {
    if (is_valid()) entry_++;
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
