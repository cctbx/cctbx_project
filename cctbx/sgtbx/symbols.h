#ifndef CCTBX_SGTBX_SYMBOLS_H
#define CCTBX_SGTBX_SYMBOLS_H

#include <cctbx/sgtbx/group_codes.h>
#include <string>

namespace cctbx { namespace sgtbx {

  namespace symbols {
    namespace tables {

      struct main_symbol_dict_entry
      {
        int         sg_number;
        const char* qualifier;
        const char* hermann_mauguin;
        const char* hall;
      };

    } // namespace tables
  } // namespace symbols

  //! class for the handling of space group symbols of various types.
  /*! The purpose of this class is to convert several conventional
      space group notations to hall() symbols by using lookup tables.
      The Hall symbols can then be used to initialize objects
      of class SpaceGroup.
      <p>
      Supported space group notations are:
      <ul>
      <li>
        <b>Space group numbers</b> as defined in the International Tables for
        Crystallography Vol. A (1983). Optionally, the space group numbers
        can be followed by a colon and one of the following characters:
        <p>
        <ul>
          <li><b><tt>1</tt></b>: Origin choice 1,
            for space groups with two origin choices
          <li><b><tt>2</tt></b>: Origin choice 2,
            for space groups with two origin choices
          <li><b><tt>H</tt></b>: Hexagonal basis,
            for rhombohedral space groups
          <li><b><tt>R</tt></b>: Rhombohedral basis,
            for rhombohedral space groups
        </ul>
        <p>
        By default, origin choice 2 or a hexagonal basis is used.
        <p>
        Examples:
        <ul>
          <li><tt>48</tt>
          <li><tt>48:1</tt>
          <li><tt>155</tt>
          <li><tt>48:R</tt>
        </ul>
        <p>
      <li>
        <b>Hermann-Mauguin symbols</b> as defined in the International Tables
        for Crystallography Vol. A (1983). Subscripts are entered without
        formatting. Optionally, subscripts can be surrounded by parentheses.
        Spaces are not required. Optionally, the Hermann-Mauguin symbols
        can be followed by a colon and one of the characters
        <tt>1</tt>,
        <tt>2</tt>,
        <tt>H</tt>, and
        <tt>R</tt>
        as explained above.
        <p>
        Examples:
        <ul>
          <li><tt>P 21 21 21</tt>
          <li><tt>R 3:R</tt>
        </ul>
        <p>
      <li>
        <b>Schoenflies symbols</b> as defined in the International Tables for
        Crystallography Vol. A (1983). Superscripts are entered by
        prepending the character '<b><tt>^</tt></b>'.
        Optionally, the Schoenflies symbols can be followed by a colon
        and one of the characters
        <tt>1</tt>,
        <tt>2</tt>,
        <tt>H</tt>, and
        <tt>R</tt>
        as explained above.
        <p>
        Examples:
        <ul>
          <li><tt>D2^4</tt>
          <li><tt>D3^7:R</tt>
        </ul>
        <p>
      <li>
        <b>Hall symbols</b> as defined in the International Tables for
        Crystallography Vol. B (2001).
        In contrast to Hermann-Mauguin symbols, Hall symbols can
        be used to specify any space group representation.
        Hall symbols are entered by prepending '<tt>Hall:</tt>'.
        <p>
        Examples:
        <ul>
          <li><tt>Hall: P 41</tt>
          <li><tt>Hall: C 2y (x,y,-x+z)</tt>
        </ul>
        <p>
        When Hall symbols are used, all the lookup algorithms provided
        by this class are bypassed. This feature is provided for
        generality and convenience. Note that number(),
        hermann_mauguin() etc. are not defined if a Hall symbol is
        used!
        <p>
      </ul>
      <p>
      See also: \ref page_multiple_cell
   */
  class space_group_symbols
  {
    public:
      //! Default constructor. Some data members are not initialized!
      space_group_symbols() {}

      //! Lookup space group Symbol.
      /*! For a general introduction see class details.<br>
          table_id is one of "" (the empty string), "I1952",
          "A1983", or "Hall".<br>
          <p>
          The default table lookup preferences are described
          in the class details section.
          <p>
          <b>I1952</b> selects the preferences of the
          International Tables for Crystallography, Volume I, 1952:<br>
          monoclinic unique axis c (example: P 2),<br>
          origin choice 1 (example: P n n n),<br>
          rhombohedral axes (example: R 3).
          <p>
          <b>A1983</b> selects the preferences of the
          International Tables for Crystallography, Volume A, 1983:<br>
          monoclinic unique axis b (example: P 2),<br>
          origin choice 1 (example: P n n n),<br>
          hexagonal axes (example: R 3).
          <p>
          <b>Hall</b> signals that Symbol is a Hall symbols
          without a leading "Hall: ". The table lookup is bypassed.
          Note that number(), hermann_mauguin() etc. are not defined
          if a Hall symbol is used!
       */
      explicit
      space_group_symbols(
        std::string const& symbol,
        std::string const& table_id="");

      //! Lookup space group number.
      /*! For a general introduction see class details.<br>
          table_id is one of "", "I1952", "A1983", or "Hall". See the
          other constructor for details.<br>
          See also: Extension()
       */
      explicit
      space_group_symbols(
        int space_group_number,
        std::string const& extension="",
        std::string const& table_id="");

      //! For internal use only.
      space_group_symbols(
        const symbols::tables::main_symbol_dict_entry* entry,
        char extension);

      //! Tests if the instance is constructed properly.
      /*! Shorthand for: number() != 0
          <p>
          Not available in Python.
       */
      bool
      is_valid() const { return number_ != 0; }

      //! Space group number according to the International Tables.
      /*! A number in the range 1 - 230. This number uniquely defines
          the space group type.<br>
          Note the distinction between "space group type" and space
          group representation" (i.e. setting, origin choice, cell
          choice, ...). For many of the 230 space group types there
          are multiple space group representations listed in the
          International Tables.
       */
      int
      number() const { return number_; }

      //! Schoenflies symbol.
      /*! One of the 230 unique Schoenflies symbols defined in the
          International Tables. A Schoenflies symbol uniquely defines
          the space group type.<br>
          Note the distinction between "space group type" and space
          group representation" (i.e. setting, origin choice, cell
          choice, ...). For many of the 230 space group types there
          are multiple space group representations listed in the
          International Tables.
       */
      std::string const&
      schoenflies() const { return schoenflies_; }

      //! A qualifier for the classification of alternative representations.
      /*! A qualifier for monoclinic and orthorhombic space groups.<br>
          For monoclinic space groups, the qualifier takes the
          form "x" or "xn", where x is one of {a, b, c, -a, -b, -c},
          and n is one of {1, 2, 3}.
          The letters define the "unique axis" according to Table 4.3.1
          in the International Tables Volume A (1983), and the numbers
          define the "cell choice."<br>
          For orthorhombic space groups, the qualifier is one of {abc,
          ba-c, cab, -cab, bca, a-cb}, according to Table 4.3.1 in the
          International Tables Volume A (1983).<br>
          Note that this qualifier is purely informational and not
          actively used in any of the symbol lookup algorithms.
       */
      std::string const&
      qualifier() const { return qualifier_; }

      //! Hermann-Mauguin symbol as defined in the International Tables.
      /*! Hermann-Mauguin (H-M) symbols were originally designed as a
          convenient description of given space-group representations.
          While it is natural to derive a H-M symbol for a given list
          of symmetry operations, it is problematic to derive the
          symmetry operations from a H-M symbol. In particular, for
          a number of space groups there is an ambiguity in the
          selection of the location of the origin with respect to the
          symmetry elements. For the conventional space group
          representations listed in the International Tables, the
          ambiguity in the origin selection is overcome by using an
          Extension().
       */
      std::string const&
      hermann_mauguin() const { return hermann_mauguin_; }

      //! Extension to the Hermann-Mauguin symbol.
      /*! For some space groups, the extension is used to distinguish
          between origin choices, or the choice of hexagonal or
          rhombohedral axes:<br>
          Extension "1": Origin choice 1.<br>
          Extension "2": Origin choice 2.<br>
          Extension "H": Hexagonal axes.<br>
          Extension "R": Rhombohedral axes.<br>
          The extension is '\0' (the null character) otherwise.<br>
          See also: hermann_mauguin()
       */
      char
      extension() const { return extension_; }

      //! Change-of-basis symbol part of universal Hermann-Mauguin symbol.
      /*! The parentheses are not included. For example, if the universal
          Hermann-Mauguin symbol is "P 3 (y,z,x)", the return value
          is "y,z,x".
       */
      std::string const&
      change_of_basis_symbol() { return change_of_basis_symbol_; }

      /*! Hermann-Mauguin symbol with extension and change-of-basis
          symbol appended (if any).
       */
      /*! If the extension is '\0' (the null character)
          and the change-of-basis operator is the identity matrix,
          universal_hermann_mauguin() is equivalent to
          hermann_mauguin().
          <p>
          The universal Hermann-Mauguin symbol uniquely identifies a
          tabulated space group representation.
       */
      std::string const&
      universal_hermann_mauguin() const
      {
        return universal_hermann_mauguin_;
      }

      //! Hall symbol.
      /*! The space group notation of Hall was designed to be "computer
          adapted". Hall symbols have some similarities with
          hermann_mauguin() symbols, but define the space group
          representation without ambiguities. Another advantage is that
          any 3-dimensional crystallographic space group representation
          can be described by a Hall symbol.<br>
          The most common use of Hall symbols in this implementation
          is to initialize objects of class SpaceGroup.
       */
      std::string const&
      hall() const { return hall_; }

      //! Determines the point group type.
      /*! The code returned is a matrix group code. There are
          exactly 32 possible return values, corresponding to
          the 32 crystallographic point group types.
          <p>
          Python: returns a string representing the point group type.
       */
      matrix_group::code
      point_group_type() const;

      //! Determines the Laue group type.
      /*! The code returned is a matrix group code. There are
          exactly 11 possible return values, corresponding to
          the 11 Laue group types.
          <p>
          Python: returns a string representing the Laue group type.
       */
      matrix_group::code
      laue_group_type() const
      {
        return point_group_type().laue_group_type();
      }

      //! Determines the crystal system.
      /*! There are exactly 7 possible return values.
          <p>
          Python: returns a string representing the crystal system.
       */
      crystal_system::code
      crystal_system() const
      {
        return point_group_type().crystal_system();
      }

    private:
      int         number_;
      std::string schoenflies_;
      std::string qualifier_;
      std::string hermann_mauguin_;
      char        extension_;
      std::string change_of_basis_symbol_;
      std::string universal_hermann_mauguin_;
      std::string hall_;
      bool
      set_all(
        const symbols::tables::main_symbol_dict_entry* entry,
        char work_extension,
        std::string const& std_table_id);
      int hall_pass_through(std::string const& symbol);
      void clear();
  };

  //! Iterator for the 530 tabulated space group representations.
  class space_group_symbol_iterator
  {
    public:
      //! Initialize the iterator.
      space_group_symbol_iterator();

      //! Get symbols for next space group representation.
      /*! result.number() == 0 signals that the end of the internal
          table was reached.
       */
      space_group_symbols next();

    private:
      const symbols::tables::main_symbol_dict_entry* entry_;
      int n_hall_;
      int i_hall_;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SYMBOLS_H
