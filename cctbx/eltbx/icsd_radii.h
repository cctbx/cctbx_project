#ifndef CCTBX_ELTBX_ICSD_RADII_H
#define CCTBX_ELTBX_ICSD_RADII_H

#include <string>

namespace cctbx { namespace eltbx { namespace icsd_radii {

  namespace detail {

    struct raw_record
    {
      const char* label;
      float       radius;
    };

  }

  //! Access to table of ionic radii.
  /*! Reference:<pre>
                          U s e r ' s  M a n u a l
                         I C S D  -  C R Y S T I N
                         =========================
                    Inorganic Crystal Structure Database
                            in conjunction with
                   Crystal Structure Information System
                           and its application to
                CCDF - Cambridge Crystallographic Data File
                                    and
                           MDF - Metal Data File
                    G.Bergerhoff, B.Kilger, C.Witthauer,
                            R.Hundt, R.Sievers
                                    Bonn
                                    1986
             Institut fuer Anorganische Chemie der Universitaet
            -------------------------------------------------------------
            ICSD/CRYSTIN User's Manual.
            English Version. Translated by Ruth Schubert.
            Updated Dec. 1986.</pre>

      These radii are also the basis of the distance tests,  which are
      routinely  carried  out when collecting data in the ICSD system.
      In  this  connection,   negative increments are obtained for the
      specially small atoms C+4, D+1, H+1 and N+5.
   */
  class table
  {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      table() : record_(0)  {}

      //! Looks up ionic radius for the given ion label.
      /*! If exact == true, the ion label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.
          <p>
          E.g., "SI4" will be matched with "Si".<br>
          "Si4+" and "Si+4" will be matched with "Si4+".<br>
          <p>
          See also: eltbx::basic::strip_label()
       */
      explicit
      table(std::string const& label, bool exact=false);

      //! Tests if the instance is constructed properly.
      /*! Shorthand for: label() != 0
          <p>
          Not available in Python.
       */
      bool
      is_valid() const { return record_->label != 0; }

      //! Label from table.
      const char*
      label() const { return record_->label; }

      //! Radius [Angstrom] from table.
      float
      radius() const { return record_->radius; }

    private:
      const detail::raw_record* record_;
      friend class table_iterator;
  };

  /*! \brief Iterator over table of radii.
   */
  class table_iterator
  {
    public:
      //! Initialization of the iterator.
      table_iterator();

      //! Retrieves the next entry from the internal table.
      /*! Use table::is_valid() to detect end-of-iteration.
       */
      table
      next();

    private:
      table current_;
  };

}}} // cctbx::eltbx::icsd_radii

#endif // CCTBX_ELTBX_ICSD_RADII_H
