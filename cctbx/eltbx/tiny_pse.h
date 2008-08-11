#ifndef CCTBX_ELTBX_TINY_PSE_H
#define CCTBX_ELTBX_TINY_PSE_H

#include <string>

namespace cctbx { namespace eltbx { namespace tiny_pse {

  namespace detail
  {
    struct raw_record {
      int z;
      const char* symbol;
      const char* name;
      float weight;
    };
  }

  //! Tiny Periodic System of Elements.
  /*! Reference:<br>
      CRC Handbook of Chemistry & Physics, 63rd edition, 1982-1983<br>
      CRC Handbook of Chemistry & Physics, 70th edition, 1989-1990
   */
  class table
  {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      table() : record_(0) {}

      //! Looks up table entry by element symbol.
      /*! If exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.
          <p>
          See also: eltbx::basic::strip_label()
       */
      explicit
      table(std::string const& label, bool exact=false);

      //! Looks up table entry by atomic number.
      explicit
      table(int atomic_number);

      //! Tests if the instance is constructed properly.
      /*! Shorthand for: atomic_number() != 0
          <p>
          Not available in Python.
       */
      bool
      is_valid() const { return record_->z != 0; }

      //! Atomic number.
      int
      atomic_number() const { return record_->z; }

      //! Element symbol.
      const char*
      symbol() const { return record_->symbol; }

      //! Element name.
      const char*
      name() const { return record_->name; }

      //! Atomic weight.
      float
      weight() const { return record_->weight; }

    private:
      const detail::raw_record* record_;
      friend class table_iterator;
  };

  /*! \brief Iterator over table entries.
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

}}} // cctbx::eltbx::tiny_pse

#endif // CCTBX_ELTBX_TINY_PSE_H
