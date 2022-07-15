#ifndef CCTBX_ELTBX_COVALENT_RADII_H
#define CCTBX_ELTBX_COVALENT_RADII_H

#include <string>

namespace cctbx { namespace eltbx { namespace covalent_radii {

  namespace detail {

    struct raw_record
    {
      const char* label;
      float       radius;
      float       esd;
    };

  }

  //! Access to table of covalent radii.
  /*! Reference:<pre>
      Beatriz Cordero, Veronica Gomez, Ana E. Platero-Prats, Marc Reves,
      Jorge Echeverria, Eduard Cremades, Flavia Barragan and Santiago Alvarez.
      Covalent radii revisited. Dalton Trans., 2008, 2832-2838,
      doi:10.1039/b801115j
   */
  class table
  {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      table() : record_(0)  {}

      //! Looks up covalent radius for the given atom label.
      /*! If exact == true, the covalent label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.
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

      //! E.s.d. [Angstrom] from table.
      float
      esd() const { return record_->esd; }

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

}}} // cctbx::eltbx::covalent_radii

#endif // CCTBX_ELTBX_COVALENT_RADII_H
