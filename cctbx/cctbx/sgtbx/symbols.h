// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_SYMBOLS_H
#define CCTBX_SGTBX_SYMBOLS_H

#include <string>

namespace sgtbx {

  namespace symbols {
    namespace tables {

      struct Main_Symbol_Dict_Entry {
        int         SgNumber;
        const char* Qualifier;
        const char* Hermann_Mauguin;
        const char* Hall;
      };

    } // namespace tables
  } // namespace symbols

  //! class for the handling of space group symbols of various types.
  /*! The purpose of this class is to convert several conventional
      space group notations to Hall() symbols by using lookup tables.
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
        generality and convenience. Note that SgNumber(),
        Hermann_Mauguin() etc. are not defined if a Hall symbol is
        used!
        <p>
      </ul>
   */
  class SpaceGroupSymbols {
    public:
      //! Lookup space group Symbol.
      /*! For a general introduction see class details.<br>
          TableId is one of "" (the empty string), "I1952",
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
          Note that SgNumber(), Hermann_Mauguin() etc. are not defined
          if a Hall symbol is used!
       */
      SpaceGroupSymbols(const std::string& Symbol,
                        const std::string& TableId = "");
      //! Lookup space group number.
      /*! For a general introduction see class details.<br>
          TableId is one of "", "I1952", "A1983", or "Hall". See the
          other constructor for details.<br>
          See also: Extension()
       */
      SpaceGroupSymbols(int SgNumber, const std::string& Extension = "",
                        const std::string& TableId = "");
      //! For internal use only.
      SpaceGroupSymbols(const symbols::tables::Main_Symbol_Dict_Entry* Entry,
                        char Extension);
      //! Space group number according to the International Tables.
      /*! A number in the range 1 - 230. This number uniquely defines
          the space group type.<br>
          Note the distinction between "space group type" and space
          group representation" (i.e. setting, origin choice, cell
          choice, ...). For many of the 230 space group types there
          are multiple space group representations listed in the
          International Tables.
       */
      inline int SgNumber() const { return m_SgNumber; }
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
      inline const std::string& Schoenflies() const { return m_Schoenflies; }
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
      inline const std::string& Qualifier() const { return m_Qualifier; }
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
      inline const std::string& Hermann_Mauguin() const {
        return m_Hermann_Mauguin; }
      //! Extension to the Hermann-Mauguin symbol.
      /*! For some space groups, the extension is used to distinguish
          between origin choices, or the choice of hexagonal or
          rhombohedral axes:<br>
          Extension "1": Origin choice 1.<br>
          Extension "2": Origin choice 2.<br>
          Extension "H": Hexagonal axes.<br>
          Extension "R": Rhombohedral axes.<br>
          The extension is '\0' (the null character) otherwise.<br>
          See also: Hermann_Mauguin()
       */
      inline char Extension() const { return m_Extension; }
      //! Hermann-Mauguin symbol with extension appended (if any).
      /*! If the extension is '\0' (the null character),
          ExtendedHermann_Mauguin() is equivalent to
          Hermann_Mauguin(). Otherwise a colon and the one-character
          extension is appended to the Hermann-Mauguin symbol.
          <p>
          The extended Hermann-Mauguin symbol uniquely identifies a
          tabulated space group representation.
       */
      const std::string& ExtendedHermann_Mauguin() const {
        return m_ExtendedHermann_Mauguin; }
      //! Hall symbol.
      /*! The space group notation of Hall was designed to be "computer
          adapted". Hall symbols have some similarities with
          Hermann_Mauguin() symbols, but define the space group
          representation without ambiguities. Another advantage is that
          any 3-dimensional crystallographic space group representation
          can be described by a Hall symbol.<br>
          The most common use of Hall symbols in this implementation
          is to initialize objects of class SpaceGroup.
       */
      inline const std::string& Hall() const { return m_Hall; }
    private:
      int         m_SgNumber;
      std::string m_Schoenflies;
      std::string m_Qualifier;
      std::string m_Hermann_Mauguin;
      char        m_Extension;
      std::string m_ExtendedHermann_Mauguin;
      std::string m_Hall;
      void SetAll(const symbols::tables::Main_Symbol_Dict_Entry* Entry,
                  char WorkExtension,
                  const std::string& StdTableId);
      int HallPassThrough(const std::string& Symbol);
      void Clear();
  };

  //! Iterator for the 656 tabulated space group representations.
  class SpaceGroupSymbolIterator {
    public:
      //! Initialize the iterator.
      SpaceGroupSymbolIterator();
      //! Get symbols for next space group representation.
      /*! result.SgNumber() == 0 signals that the end of the internal
          table was reached.
       */
      SpaceGroupSymbols next();
    private:
      const symbols::tables::Main_Symbol_Dict_Entry* Entry;
      int nHall;
      int iHall;
  };

} // namespace sgtbx

#endif // CCTBX_SGTBX_SYMBOLS_H
