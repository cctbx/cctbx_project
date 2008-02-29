#ifndef IOTBX_PDB_INPUT_H
#define IOTBX_PDB_INPUT_H

#include <iotbx/pdb/hierarchy_v1.h>
#include <iotbx/pdb/hierarchy_v2.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/tiny.h>

namespace iotbx { namespace pdb {

  static const double anisou_factor = 1.e-4;

  static const char blank_altloc_char = ' ';

  //! Helper for looping over index ranges.
  template <typename ElementType>
  struct range_loop
  {
    range_loop() {}

    range_loop(
      af::const_ref<ElementType> const& indices,
      unsigned begin=0)
    :
      i_end(indices.end()),
      i(indices.begin()),
      end(begin)
    {}

    range_loop(
      std::vector<ElementType> const& indices,
      unsigned begin=0)
    :
      i_end(indices.size() == 0 ? 0 : (&*indices.begin()) + indices.size()),
      i(indices.size() == 0 ? 0 : &*indices.begin()),
      end(begin)
    {}

    bool
    next()
    {
      if (i == i_end) return false;
      begin = end;
      end = static_cast<unsigned>(*i++);
      size = end - begin;
      return true;
    }

    void
    skip_to_last()
    {
      if (i != i_end) i = i_end-1;
    }

    protected:
      const ElementType* i_end;
      const ElementType* i;
    public:
      unsigned begin;
      unsigned end;
      unsigned size;
  };

  //! Facilitates fast processing and comprehensive error messages.
  class line_info
  {
    public:
      const char* source_info;
      unsigned line_number;
      const char* data;
      unsigned size;

    protected:
      std::string error_source_info_;
      unsigned error_line_number_;
      std::string error_line_;
      unsigned error_column_;
      std::string error_message_;

    public:
      line_info() {}

      line_info(const char* source_info_)
      :
        source_info(source_info_),
        line_number(0),
        error_column_(0)
      {}

      template <typename StringType>
      void
      set_error(
        unsigned error_column,
        StringType error_message)
      {
        if (error_column_ != 0) return;
        error_source_info_ = (source_info ? source_info : "");
        error_line_number_ = line_number;
        error_line_ = std::string(data, size);
        error_column_ = error_column;
        error_message_ = error_message;
      }

      bool
      error_occured() const { return error_column_ != 0; }

      std::string
      format_exception_message() const;

      void
      check_and_throw_runtime_error() const
      {
        if (error_column_ == 0) return;
        throw std::runtime_error(format_exception_message());
      }

      bool
      is_blank_line() const
      {
        for(unsigned i=0;i<size;i++) {
          if (data[i] != ' ') return false;
        }
        return true;
      }

      std::string
      strip_data(unsigned start_at_column=0) const
      {
        unsigned sz = size;
        while (sz > start_at_column) {
          if (data[--sz] != ' ') {
            ++sz;
            break;
          }
        }
        if (sz <= start_at_column) return std::string();
        return std::string(data+start_at_column, data+sz);
      }
  };

  //! Detects old-style PDB files with the PDB access code in columns 73-76.
  struct columns_73_76_evaluator
  {
    const char* finding;
    bool is_old_style;
    unsigned number_of_atom_and_hetatm_lines;

    typedef af::tiny<char, 4> columns_73_76_t;

    struct columns_73_76_t_lexical_less_than
    {
      bool operator()(columns_73_76_t const& a, columns_73_76_t const& b) const
      {
        if (a[0] < b[0]) return true;
        if (a[0] > b[0]) return false;
        if (a[1] < b[1]) return true;
        if (a[1] > b[1]) return false;
        if (a[2] < b[2]) return true;
        if (a[2] > b[2]) return false;
        if (a[3] < b[3]) return true;
        return false;
      }
    };

    typedef std::map<
      columns_73_76_t,
      unsigned,
      columns_73_76_t_lexical_less_than> columns_73_76_dict_t;

    static
    bool
    is_record_type(const char* name, const char* line_data)
    {
      for(unsigned i=0;i<6;i++) {
        if (line_data[i] != name[i]) return false;
      }
      return true;
    }

    columns_73_76_evaluator() {}

    columns_73_76_evaluator(
      af::const_ref<std::string> const& lines,
      unsigned is_frequent_threshold_atom_records=1000,
      unsigned is_frequent_threshold_other_records=100);
  };

  //! Efficient processing of input atom labels.
  struct input_atom_labels
  {
    static const unsigned compacted_size = 4+3+2+4+1+4+1;
    char compacted[compacted_size];

    char*       name_begin()       { return compacted; }
    const char* name_begin() const { return compacted; }
    str4        name_small() const { return str4(name_begin(), true); }
    std::string name()       const { return std::string(name_begin(),4); }

    char*       resname_begin()       { return compacted+4; }
    const char* resname_begin() const { return compacted+4; }
    str3        resname_small() const { return str3(resname_begin(), true); }
    std::string resname()       const { return std::string(resname_begin(),3);}

    char*       chain_begin()       { return compacted+7; }
    const char* chain_begin() const { return compacted+7; }
    str2        chain_small() const
    {
      if (chain_begin()[0] == ' ') return str2(chain_begin()[1]);
      return str2(chain_begin(), true);
    }
    std::string chain()       const
    {
      if (chain_begin()[0] == ' ') return std::string(chain_begin()+1,1);
      return std::string(chain_begin(),2);
    }

    char*       resseq_begin()       { return compacted+9; }
    const char* resseq_begin() const { return compacted+9; }
    str4        resseq_small() const { return str4(resseq_begin(), true); }
    std::string resseq()       const { return std::string(resseq_begin(),4); }

    char*       icode_begin()       { return compacted+13; }
    const char* icode_begin() const { return compacted+13; }
    str1        icode_small() const { return str1(icode_begin(), true); }
    std::string icode()       const { return std::string(icode_begin(),1); }

    char*       resid_begin()       { return compacted+9; }
    const char* resid_begin() const { return compacted+9; }
    str5        resid_small() const { return str5(resid_begin(), true); }
    std::string resid()       const { return std::string(resid_begin(),5); }

    char*       segid_begin()       { return compacted+14; }
    const char* segid_begin() const { return compacted+14; }
    str4        segid_small() const { return str4(segid_begin(), true); }
    std::string segid()       const { return std::string(segid_begin(),4); }

    char*       altloc_begin()       { return compacted+18; }
    const char* altloc_begin() const { return compacted+18; }
    str1        altloc_small() const { return str1(altloc_begin(), true); }
    std::string altloc()       const { return std::string(altloc_begin(),1); }

    str4 altloc_resname_small() const
    {
      str4 result(compacted+3, true);
      result.elems[0] = *altloc_begin();
      return result;
    }

    input_atom_labels() {}

    input_atom_labels(pdb::line_info& line_info)
    {
      //  7 - 11  Integer       serial   Atom serial number.
      // 13 - 16  Atom          name     Atom name.
      // 17       Character     altLoc   Alternate location indicator.
      // 18 - 20  Residue name  resName  Residue name.
      // 21 - 22                chainID  Chain identifier.
      // 23 - 26  Integer       resSeq   Residue sequence number.
      // 27       AChar         iCode    Code for insertion of residues.
      // 73 - 76  LString(4)    segID    Segment identifier, left-justified.
      extract(line_info,12,4,name_begin());
      extract(line_info,16,1,altloc_begin());
      extract(line_info,17,3,resname_begin());
      extract(line_info,20,2,chain_begin());
      extract(line_info,22,4,resseq_begin());
      extract(line_info,26,1,icode_begin());
      extract(line_info,72,4,segid_begin());
    }

    static
    af::tiny<char, 5>
    extract_serial_number_string(
      pdb::line_info& line_info)
    {
      af::tiny<char, 5> result;
      extract(line_info,6,5,result.begin());
      return result;
    }

    static
    void
    extract(
      pdb::line_info& line_info,
      unsigned i_begin,
      unsigned n,
      char* target)
    {
      unsigned j = 0;
      while (i_begin < line_info.size) {
        if (j == n) return;
        target[j++] = line_info.data[i_begin++];
      }
      while (j < n) target[j++] = ' ';
    }

    static
    bool
    are_equal(
      pdb::line_info& line_info,
      unsigned i_begin,
      unsigned n,
      const char* target)
    {
      unsigned j = 0;
      while (i_begin < line_info.size) {
        if (j == n) return true;
        if (target[j++] != line_info.data[i_begin++]) return false;
      }
      while (j < n) if (target[j++] != ' ') return false;
      return true;
    }

    std::string
    pdb_format() const;

    void
    check_equivalence(pdb::line_info& line_info) const;

    int
    compare(input_atom_labels const& other) const
    {
      const char* s = compacted;
      const char* o = other.compacted;
      if (*  s < *  o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      if (*++s < *++o) return -1; if (*s > *o) return 1;
      return 0;
    }

    bool
    operator!=(input_atom_labels const& other) const
    {
      const char* s = compacted;
      const char* o = other.compacted;
      if (*  s != *  o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      if (*++s != *++o) return true;
      return false;
    }

    bool
    equal_ignoring_altloc(input_atom_labels const& other) const
    {
      const char* s = compacted;
      const char* o = other.compacted;
      if (*  s != *  o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      return true;
    }

    bool
    is_in_same_residue(input_atom_labels const& other) const
    {
      const char* s = resname_begin();
      const char* o = other.resname_begin();
      if (*  s != *  o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      s = resseq_begin();
      o = other.resseq_begin();
      if (*  s != *  o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      if (*++s != *++o) return false;
      return (*icode_begin() == *other.icode_begin());
    }
  };

  //! Processing of PDB strings.
  class input
  {
    public:
      typedef std::map<str6, unsigned> record_type_counts_t;
      typedef std::map<std::string, std::vector<unsigned> > str_sel_cache_t;

      input() {}

      input(std::string const& file_name);

      input(
        const char* source_info,
        af::const_ref<std::string> const& lines);

    protected:
      void
      process(af::const_ref<std::string> const& lines);

    public:
      std::string const&
      source_info() const { return source_info_; }

      record_type_counts_t const&
      record_type_counts() const { return record_type_counts_; }

      af::shared<std::string> const&
      unknown_section() const { return unknown_section_; }

      af::shared<std::string> const&
      title_section() const { return title_section_; }

      af::shared<std::string> const&
      remark_section() const { return remark_section_; }

      af::shared<std::string> const&
      primary_structure_section() const { return primary_structure_section_; }

      af::shared<std::string> const&
      heterogen_section() const { return heterogen_section_; }

      af::shared<std::string> const&
      secondary_structure_section() const
      {
        return secondary_structure_section_;
      }

      af::shared<std::string> const&
      connectivity_annotation_section() const
      {
        return connectivity_annotation_section_;
      }

      af::shared<std::string> const&
      miscellaneous_features_section() const
      {
        return miscellaneous_features_section_;
      }

      af::shared<std::string> const&
      crystallographic_section() const { return crystallographic_section_; }

      af::shared<input_atom_labels> const&
      input_atom_labels_list() const { return input_atom_labels_list_; }

      af::shared<std::string> const&
      atom_serial_number_strings() const { return atom_serial_number_strings_;}

      af::shared<hierarchy_v1::atom> const&
      atoms() const { return atoms_; }

      af::shared<hierarchy_v2::atom>
      atoms_v2();

      af::shared<int> const&
      model_numbers() const { return model_numbers_; }

      af::shared<std::size_t> const&
      model_indices() const { return model_indices_; }

      af::shared<std::size_t> const&
      ter_indices() const { return ter_indices_; }

      af::shared<std::vector<unsigned> > const&
      chain_indices() const { return chain_indices_; }

      af::shared<std::size_t> const&
      break_indices() const { return break_indices_; }

      af::shared<std::string> const&
      connectivity_section() const { return connectivity_section_; }

      af::shared<std::string> const&
      bookkeeping_section() const { return bookkeeping_section_; }

#define IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(attr, type) \
      str_sel_cache_t const& \
      attr##_selection_cache() const \
      { \
        if (!attr##_selection_cache_is_up_to_date_) { \
          /* intermediate map to avoid frequent construction of string */ \
          typedef std::map<type, std::vector<unsigned> > interm_map_t; \
          interm_map_t sel_cache; \
          unsigned i_seq = 0; \
          const input_atom_labels* iall_end = input_atom_labels_list_.end();\
          for(const input_atom_labels* ial=input_atom_labels_list_.begin(); \
                                       ial!=iall_end;ial++) { \
            sel_cache[ial->attr##_small()].push_back(i_seq++); \
          } \
          for(interm_map_t::iterator sc=sel_cache.begin(); \
                                     sc!=sel_cache.end();sc++) { \
            attr##_selection_cache_[     /* swap avoids deep copy */ \
              std::string(sc->first.elems)].swap(sc->second); \
          } \
          attr##_selection_cache_is_up_to_date_ = true; \
        } \
        return attr##_selection_cache_; \
      }

      IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(name, str4)
      IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(altloc, str1)
      IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(resname, str3)
      IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(chain, str2)
      IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(resseq, str4)
      IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(icode, str1)
      IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(segid, str4)

      bool
      model_numbers_are_unique() const;

      af::shared<std::size_t>
      model_atom_counts() const;

      af::shared<std::vector<unsigned> >
      find_duplicate_atom_labels() const;

      //! not const because atom parents are modified.
      hierarchy_v1::root
      construct_hierarchy_v1();

      //! not const because atom parents are modified.
      hierarchy_v2::root
      construct_hierarchy_v2();

      unsigned
      number_of_alternative_groups_with_blank_altloc() const
      {
        SCITBX_ASSERT(construct_hierarchy_v1_was_called_before);
        return number_of_alternative_groups_with_blank_altloc_;
      }

      unsigned
      number_of_alternative_groups_without_blank_altloc() const
      {
        SCITBX_ASSERT(construct_hierarchy_v1_was_called_before);
        return number_of_alternative_groups_without_blank_altloc_;
      }

      unsigned
      number_of_chains_with_altloc_mix() const
      {
        SCITBX_ASSERT(construct_hierarchy_v1_was_called_before);
        return number_of_chains_with_altloc_mix_;
      }

      af::shared<std::size_t> const&
      i_seqs_alternative_group_with_blank_altloc() const
      {
        SCITBX_ASSERT(construct_hierarchy_v1_was_called_before);
        return i_seqs_alternative_group_with_blank_altloc_;
      }

      af::shared<std::size_t> const&
      i_seqs_alternative_group_without_blank_altloc() const
      {
        SCITBX_ASSERT(construct_hierarchy_v1_was_called_before);
        return i_seqs_alternative_group_without_blank_altloc_;
      }

      std::map<std::string, unsigned>
      atom_element_counts() const;

#define IOTBX_PDB_INPUT_EXTRACT_ATOM_ATTR(attr, attr_type) \
      af::shared<attr_type > \
      extract_atom_##attr() const \
      { \
        af::shared<attr_type > result( \
          atoms_.size(), af::init_functor_null<attr_type >()); \
        attr_type* r = result.begin(); \
        const hierarchy_v1::atom* atoms_end = atoms_.end(); \
        for(const hierarchy_v1::atom* a=atoms_.begin();a!=atoms_end;a++) { \
          *r++ = a->data->attr; \
        } \
        return result; \
      }

      IOTBX_PDB_INPUT_EXTRACT_ATOM_ATTR(xyz, vec3)
      IOTBX_PDB_INPUT_EXTRACT_ATOM_ATTR(sigxyz, vec3)
      IOTBX_PDB_INPUT_EXTRACT_ATOM_ATTR(occ, double)
      IOTBX_PDB_INPUT_EXTRACT_ATOM_ATTR(sigocc, double)
      IOTBX_PDB_INPUT_EXTRACT_ATOM_ATTR(b, double)
      IOTBX_PDB_INPUT_EXTRACT_ATOM_ATTR(sigb, double)
      IOTBX_PDB_INPUT_EXTRACT_ATOM_ATTR(uij, sym_mat3)
      IOTBX_PDB_INPUT_EXTRACT_ATOM_ATTR(siguij, sym_mat3)

      af::shared<std::size_t>
      extract_atom_hetero() const;

      af::shared<std::size_t>
      extract_atom_flag_altloc() const;

      void
      reset_atom_tmp(int first_value=0, int increment=1) const;

    protected:
      std::string source_info_;
      record_type_counts_t record_type_counts_;
      af::shared<std::string> unknown_section_;
      af::shared<std::string> title_section_;
      af::shared<std::string> remark_section_;
      af::shared<std::string> primary_structure_section_;
      af::shared<std::string> heterogen_section_;
      af::shared<std::string> secondary_structure_section_;
      af::shared<std::string> connectivity_annotation_section_;
      af::shared<std::string> miscellaneous_features_section_;
      af::shared<std::string> crystallographic_section_;
      af::shared<input_atom_labels> input_atom_labels_list_;
      af::shared<std::string> atom_serial_number_strings_;
      af::shared<hierarchy_v1::atom> atoms_;
      af::shared<hierarchy_v2::atom> atoms_v2_;
      af::shared<int>         model_numbers_;
      af::shared<std::size_t> model_indices_;
      af::shared<std::size_t> ter_indices_;
      af::shared<std::vector<unsigned> > chain_indices_;
      af::shared<std::size_t> break_indices_;
      af::shared<unsigned>    break_record_line_numbers;
      af::shared<std::string> connectivity_section_;
      af::shared<std::string> bookkeeping_section_;
      bool construct_hierarchy_v1_was_called_before;
      unsigned number_of_alternative_groups_with_blank_altloc_;
      unsigned number_of_alternative_groups_without_blank_altloc_;
      unsigned number_of_chains_with_altloc_mix_;
      af::shared<std::size_t> i_seqs_alternative_group_with_blank_altloc_;
      af::shared<std::size_t> i_seqs_alternative_group_without_blank_altloc_;
      //
      mutable bool            name_selection_cache_is_up_to_date_;
      mutable str_sel_cache_t name_selection_cache_;
      mutable bool            altloc_selection_cache_is_up_to_date_;
      mutable str_sel_cache_t altloc_selection_cache_;
      mutable bool            resname_selection_cache_is_up_to_date_;
      mutable str_sel_cache_t resname_selection_cache_;
      mutable bool            chain_selection_cache_is_up_to_date_;
      mutable str_sel_cache_t chain_selection_cache_;
      mutable bool            resseq_selection_cache_is_up_to_date_;
      mutable str_sel_cache_t resseq_selection_cache_;
      mutable bool            icode_selection_cache_is_up_to_date_;
      mutable str_sel_cache_t icode_selection_cache_;
      mutable bool            segid_selection_cache_is_up_to_date_;
      mutable str_sel_cache_t segid_selection_cache_;
  };

}} // namespace iotbx::pdb

#endif // IOTBX_PDB_INPUT_H
