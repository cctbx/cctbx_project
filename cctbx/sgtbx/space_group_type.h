#ifndef CCTBX_SGTBX_SPACE_GROUP_TYPE_H
#define CCTBX_SGTBX_SPACE_GROUP_TYPE_H

#include <cctbx/sgtbx/space_group.h>

namespace cctbx { namespace sgtbx {

  //! Computation and use of change-of-basis matrix to reference setting.
  /*! A space group type is characterized by the space group number
      according to the International Tables for Crystallography,
      Volume A, and a change-of-basis matrix that transforms
      the coordinates of a given setting to a reference setting.
   */
  class space_group_type
  {
    public:
      //! Default constructor.
      /*! Equivalent to space_group_type(space_group())
       */
      space_group_type() : number_(1), cb_op_is_tidy_(true) {}

      //! Initializer using space group symbols.
      /*! Equivalent to
          space_group_type(
            space_group(space_group_symbols(symbol, table_id)),
            tidy_cb_op)
       */
      explicit
      space_group_type(
        std::string const& symbol,
        std::string const& table_id="",
        bool tidy_cb_op=true);

      //! Determines the space group type.
      /*! The algorithm for the determination of the space group
          type is published in
          <A HREF="http://journals.iucr.org/a/issues/1999/02/02/au0146/"
          ><I>Acta Cryst.</I> 1999, <B>A55</B>:383-395</A>.
          The result is the space group number according to the
          International Tables for Crystallography, and a
          change-of-basis matrix that transforms the given space group
          to its reference setting. r_den and t_den are the rotation part
          denominator and the translation part denominator, respectively,
          for the change-of-basis matrix.<p>
          In general, there are several change-of-basis matrices that
          could be used. If tidy_cb_op == true, the operations of the
          affine normalizer are used to determine a "standard"
          change-of-basis matrix. This is, the change-of-basis matrix
          is independent of the order of the given symmetry operations.
          In principle, the space group number and the "standard"
          change-of-basis matrix could be used to test for the
          equivalence of two groups of symmetry operations (but the
          algorithm used for operator== is more efficient).
          On average, the additional time required for the
          determination of the "standard" change-of-basis matrix
          is about 70% of the time for the main determination
          of the space group type. This is, the tidy_cb_op == true
          option is relatively expensive. However, the absolute
          runtime is only a small fraction of a second.
       */
      explicit
      space_group_type(
        space_group const& group,
        bool tidy_cb_op=true,
        int r_den=cb_r_den,
        int t_den=cb_t_den);

      //! Access to space group passed to the constructor.
      space_group const&
      group() const { return group_; }

      //! Space group number according to the International Tables.
      int
      number() const { return number_; }

      //! Change-of-basis operator.
      change_of_basis_op const&
      cb_op() const { return cb_op_; }

      //! Setting as passed to the constructor via tidy_cp_op.
      bool
      cb_op_is_tidy() const { return cb_op_is_tidy_; }

      //! Gets the additional generators of the Euclidean normalizer.
      /*! See International Tables for Crystallography Volume A,
          1983, Table 15.3.2. The generators are tabulated for
          reference settings and transformed to the given setting
          using cb_op().
          <p>
          flag_k2l = true requests the additional generator from
          column "Inversion through a centre at" of Table 15.3.2.
          <p>
          flag_l2n = true requests the additional generators from
          column "Further generators" of Table 15.3.2.
          <p>
          See also: class StructureSeminvariant
       */
      af::shared<rt_mx>
      addl_generators_of_euclidean_normalizer(
        bool flag_k2l,
        bool flag_l2n) const;

      /*! \brief Adds the additional generators of the Euclidean normalizer
          to the space group.
       */
      /*! See also: addl_generators_of_euclidean_normalizer(),
                    space_group::expand_smx()
       */
      space_group
      expand_addl_generators_of_euclidean_normalizer(
        bool flag_k2l,
        bool flag_l2n) const
      {
        space_group result = group_;
        result.expand_smx(
          addl_generators_of_euclidean_normalizer(flag_k2l, flag_l2n)
            .const_ref());
        return result;
      }

      //! Tests for the 22 (11 pairs) enantiomorphic space groups.
      /*! A space group G is enantiomorphic if G and -I.G.-I have
          two different space group types. I is the unit matrix.<br>
          There are 11 pairs of enantiomorphic space groups,
          for example P41 and P43.
          <p>
          The notion of enantiomorphic space groups is connected
          to the notion of <i>affine</i> space group types. There
          are 219 affine space group types, compared to the
          230 conventional space group types. A pair of
          enantiomorphic space groups corresponds to a single
          affine space group type.
          <p>
          See also: change_of_hand_op()
       */
      bool
      is_enantiomorphic() const;

      //! Simple table lookup.
      /*! A space group is symmorphic if, apart from the lattice translations,
          all symmetry operations leave one common point fixed. This is
          true for 73 space group types, e.g. C222.
          <p>
          See also:
            http://reference.iucr.org/dictionary/Symmorphic_space_groups
       */
      bool
      is_symmorphic() const;

      //! Determines a change-of-hand matrix.
      /*! If the space group is centro-symmetric, the change-of-hand
          matrix is the identity matrix.
          <p>
          If the space group belongs to one of the 22 enantiomorphic
          space group types, the change-of-hand matrix is determined as
          a centre of inversion that is located at the origin of the
          reference setting.
          <p>
          If the space group is not enantiomorphic and not
          centro-symmetric, the change-of-hand matrix is determined as
          a centre of inversion of the Euclidean normalizer.
          <p>
          The change-of-hand matrix can be used to transform the
          symmetry operations to obtain the enantiomorph symmetry
          operations (use space_group::change_basis()), and to transform
          fractional coordinates. For example:<pre>
          space_group sg = ...;
          space_group_type sg_type(sg);
          change_of_basis_op cb_op = sg_type.change_of_hand_op();
          space_group enantiomorph_sg = sg.change_basis(cb_op);
          fractional<> x = ...;
          fractional<> enatiomorph_x = cb_op(x);</pre>
          <p>
          See also: addl_generators_of_euclidean_normalizer(),
                    is_enantiomorphic(), space_group::change_basis()
       */
      change_of_basis_op
      change_of_hand_op() const;

      //! Builds a Hall symbol for the given symmetry operations.
      /*! For a given group of symmetry operations, there are in
          general several plausible choices for a Hall symbol.
          The Hall symbol returned by this algorithm is derived
          from the Hall symbol for the reference setting. A change-of-basis
          operator is attached if necessary (i.e. "P 2 (y,z,x)").
          <p>
          If tidy_cb_op == true, the returned Hall symbol
          is a reproducible representation of a
          given setting that is independent of the order of
          the symmetry operations.
       */
      std::string
      hall_symbol(bool tidy_cb_op = true) const;

      /*! \brief Builds a universal Hermann-Mauguin symbol for the given
          symmetry operations.
       */
      /*! If tidy_cb_op == true, the returned symbol
          is a reproducible representation of a
          given setting that is independent of the order of
          the symmetry operations.
       */
      std::string
      universal_hermann_mauguin_symbol(bool tidy_cb_op = true) const;

      //! Determines conventional Hermann-Mauguin symbol or Hall symbol.
      /*! First, group().match_tabulated_settings() is called. If the given
          symmetry operations correspond to one of the 530 tabulated
          settings, the extended Hermann-Mauguin symbol is returned
          (e.g. "P n n n :2").  Otherwise the string "Hall: " + the
          result of hall_symbol() is returned (e.g. "Hall:  P 2 2
          -1n").
          <p>
          The result of lookup_symbol() can be used as the input
          for the constructor of class space_group_symbols.
          <p>
          If ad_hoc_1992 is true, symbols for space groups
          39, 41, 64, 67, and 68 will contain the "e" character:
          <p>
          Acta Cryst. (1992). A48, 727-732.
          Symbols for symmetry elements and symmetry operations.
          Final report of the International Union of Crystallography
          Ad-Hoc Committee on the nomenclature of symmetry.
          <p>
          http://www.iucr.org/iucr-top/comm/cnom/symsym/node7.html
          <pre>
          Space group No.   39   41   64   67   68
          Symbol in ITA83: Abm2 Aba2 Cmca Cmma Ccca
               New symbol: Aem2 Aea2 Cmce Cmme Ccce
          </pre>
          <p>
          The "e" symbols are ambiguous for space groups 67 and 68:
          <pre>
          No. 67:
            Cmme: Cmma or Cmmb
            Aemm: Abmm or Acmm
            Bmem: Bmcm or Bmam
          No. 68
            Ccce: Ccca or Cccb
            Aeaa: Abaa or Acaa
            Bbeb: Bbcb or Bbab
          </pre>
          This function produces the first choice, e.g. Cmme = Cmma.
       */
      std::string
      lookup_symbol(
        bool ad_hoc_1992=false) const;

    protected:
      space_group group_;
      int number_;
      change_of_basis_op cb_op_;
      bool cb_op_is_tidy_;
      mutable std::string hall_symbol_tidy_true_;
      mutable std::string hall_symbol_tidy_false_;
      mutable std::string uhm_symbol_tidy_true_;
      mutable std::string uhm_symbol_tidy_false_;
      mutable std::string lookup_symbol_;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SPACE_GROUP_TYPE_H
