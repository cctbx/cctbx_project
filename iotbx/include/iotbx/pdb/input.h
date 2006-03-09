#ifndef IOTBX_PDB_INPUT_H
#define IOTBX_PDB_INPUT_H

#include <iotbx/pdb/hierarchy.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/auto_array.h>
#include <set>

namespace iotbx { namespace pdb {

  static const int min_resseq = -999;
  static const int max_resseq = 9999;
  static const double anisou_factor = 1.e-4;

  static const char        blank_altloc_char = ' ';
  static const std::string blank_altloc_stl(" ");

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
      format_exception_message() const
      {
        std::string result;
        if (error_source_info_.size() != 0) {
          result += error_source_info_;
          if (error_line_number_ != 0) {
            result += ", ";
          }
        }
        else if (error_line_number_ != 0) {
          result += "input ";
        }
        if (error_line_number_ != 0) {
          char buf[64];
          std::sprintf(buf, "line %u", error_line_number_);
          result += buf;
        }
        if (result.size() == 0) {
          result = "input line";
        }
        return result + ":\n  " + error_line_
          + "\n  " + std::string(std::max(1u,error_column_)-1, '-')
          + "^\n  " + error_message_;
      }

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

  inline
  int
  field_as_int(
    pdb::line_info& line_info,
    unsigned i_begin,
    unsigned i_end)
  {
    int result = 0;
    int sign = 0;
    bool have_non_blank = false;
    unsigned j_end = std::min(i_end, line_info.size);
    while (i_begin < j_end) {
      char c = line_info.data[i_begin++];
      switch (c) {
        case '+':
          if (have_non_blank || sign != 0) {
            line_info.set_error(i_begin, "unexpected plus sign.");
            return 0;
          }
          sign = 1;
          break;
        case '-':
          if (have_non_blank || sign != 0) {
            line_info.set_error(i_begin, "unexpected minus sign.");
            return 0;
          }
          sign = -1;
          break;
        case ' ':
        case '0': result *= 10;              break;
        case '1': result *= 10; result += 1; break;
        case '2': result *= 10; result += 2; break;
        case '3': result *= 10; result += 3; break;
        case '4': result *= 10; result += 4; break;
        case '5': result *= 10; result += 5; break;
        case '6': result *= 10; result += 6; break;
        case '7': result *= 10; result += 7; break;
        case '8': result *= 10; result += 8; break;
        case '9': result *= 10; result += 9; break;
        default:
          line_info.set_error(i_begin, "unexpected character.");
          return 0;
      }
      if (c != ' ') have_non_blank = true;
    }
    while (i_begin++ < i_end) result *= 10;
    if (sign < 0) return -result;
    return result;
  }

  // i_end-i_begin must be less than 32. This is not checked to
  // maximize performance.
  inline
  double
  field_as_double(
    pdb::line_info& line_info,
    unsigned i_begin,
    unsigned i_end)
  {
    char buf[32];
    char* b = buf;
    bool have_non_blank = false;
    unsigned i_col = i_begin;
    unsigned j_end = std::min(i_end, line_info.size);
    while (i_col < j_end) {
      char c = line_info.data[i_col++];
      if (c == ' ') {
        if (have_non_blank) *b++ = '0';
        else                i_begin++;
      }
      else {
        if (c == 'x' || c == 'X') {
          // suppress conversion of hexadecimal literals
          c = '?'; // any character that triggers an error in strtod()
        }
        *b++ = c;
        have_non_blank = true;
      }
    }
    if (!have_non_blank) return 0;
    while (i_col++ < i_end) *b++ = '0';
    *b = '\0';
    char *endptr;
    double result = std::strtod(buf, &endptr);
    if (endptr == buf) {
      line_info.set_error(i_begin+1, "not a floating-point number.");
    }
    if (endptr != b) {
      line_info.set_error(i_begin+1+(endptr-buf), "unexpected character.");
    }
    return result;
  }

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
      unsigned is_frequent_threshold_other_records=100)
    :
      finding("Undecided."),
      is_old_style(false),
      number_of_atom_and_hetatm_lines(0)
    {
      columns_73_76_dict_t atom_columns_73_76_dict;
      columns_73_76_dict_t other_columns_73_76_dict;
      for(const std::string* line=lines.begin();line!=lines.end();line++) {
        if (line->size() < 80) continue;
        const char* line_data = line->data();
        if (line_data[79] == ' ') continue;
        af::tiny<char, 4> columns_73_76;
        {
          bool have_non_blank = false;
          const char* l = line_data+72;
          unsigned i = 0;
          while (i < 4) {
            if (*l != ' ') have_non_blank = true;
            columns_73_76[i++] = *l++;
          }
          if (!have_non_blank) continue;
        }
        if (   is_record_type("ATOM  ", line_data)
            || is_record_type("HETATM", line_data)) {
          atom_columns_73_76_dict[columns_73_76]++;
          number_of_atom_and_hetatm_lines++;
        }
        else if (!(   is_record_type("SIGATM", line_data)
                   || is_record_type("ANISOU", line_data)
                   || is_record_type("SIGUIJ", line_data)
                   || is_record_type("TER   ", line_data))) {
          other_columns_73_76_dict[columns_73_76]++;
        }
      }
      if (atom_columns_73_76_dict.size() == 0) {
        finding = "Blank columns 73-76 on ATOM and HETATM records.";
        return;
      }
      if (   atom_columns_73_76_dict.size() == 1
          && other_columns_73_76_dict.size() == 1) {
        columns_73_76_t const& a = atom_columns_73_76_dict.begin()->first;
        columns_73_76_t const& o = other_columns_73_76_dict.begin()->first;
        if (a.all_eq(o) && a[0] != ' ' && a[3] != ' ') {
          finding = "Exactly one common label in columns 73-76.";
          is_old_style = true;
          return;
        }
      }
      unsigned sum_counts_common_four_character_field = 0;
      for(columns_73_76_dict_t::const_iterator
          a_item=atom_columns_73_76_dict.begin();
          a_item!=atom_columns_73_76_dict.end();
          a_item++) {
        columns_73_76_t const& a = a_item->first;
        if (a[0] == ' ' || a[3] == ' ') continue;
        columns_73_76_dict_t::const_iterator
          o_item = other_columns_73_76_dict.find(a);
        if (o_item != other_columns_73_76_dict.end()) {
          if (   a_item->second > is_frequent_threshold_atom_records
              && o_item->second > is_frequent_threshold_other_records) {
            finding = "Frequent common labels in columns 73-76.";
            is_old_style = true;
            return;
          }
          unsigned sum_counts = a_item->second + o_item->second;
          if (sum_counts_common_four_character_field < sum_counts) {
              sum_counts_common_four_character_field = sum_counts;
          }
        }
      }
      if (sum_counts_common_four_character_field == 0) {
        finding = "No common label in columns 73-76.";
        is_old_style = false;
        return;
      }
      // Compare first three characters only.
      columns_73_76_t a_begin = atom_columns_73_76_dict.begin()->first;
      if (a_begin[0] == ' ' || a_begin[2] == ' ') return;
      for(columns_73_76_dict_t::const_iterator
            a_item=atom_columns_73_76_dict.begin();
            a_item!=atom_columns_73_76_dict.end();
            a_item++) {
        columns_73_76_t const& a = a_item->first;
        if (   a[0] != a_begin[0]
            || a[1] != a_begin[1]
            || a[2] != a_begin[2]) return;
      }
      for(columns_73_76_dict_t::const_iterator
            o_item=other_columns_73_76_dict.begin();
            o_item!=other_columns_73_76_dict.end();
            o_item++) {
        columns_73_76_t const& o = o_item->first;
        if (   o[0] != a_begin[0]
            || o[1] != a_begin[1]
            || o[2] != a_begin[2]) return;
      }
      finding = "Exactly one common label in columns 73-76"
                " comparing only the first three characters.";
      is_old_style = true;
    }
  };

  namespace record_type {

    enum section {
      unknown,
      title,
      primary_structure,
      heterogen,
      secondary_structure,
      connectivity_annotation,
      miscellaneous_features,
      crystallographic,
      coordinate,
      connectivity,
      bookkeeping
    };

    enum id {
      xxxxxx,
      // Title Section
      header,
      onhold, // non-standard
      obslte,
      title_,
      caveat,
      compnd,
      source,
      keywds,
      expdta,
      author,
      revdat,
      sprsde,
      jrnl__,
      remark,
      ftnote, // non-standard
      // Primary Structure Section
      dbref_,
      seqadv,
      seqres,
      modres,
      // Heterogen Section
      het___,
      hetnam,
      hetsyn,
      formul,
      // Secondary Structure Section
      helix_,
      sheet_,
      turn__,
      // Connectivity Annotation Section
      ssbond,
      link__,
      hydbnd,
      sltbrg,
      cispep,
      // Miscellaneous Features Section
      site__,
      // Crystallographic and Coordinate Transformation Section
      cryst1,
      origx_,
      scale_,
      mtrix_,
      tvect_,
      // Coordinate Section
      model_,
      atom__,
      sigatm,
      anisou,
      siguij,
      ter___,
      break_, // non-standard
      hetatm,
      endmdl,
      // Connectivity Section
      conect,
      // Bookkeeping Section
      master,
      end___
    };

    /* requirements (not checked to maximize performance):
         line_size > 0
         name must be 6 characters exactly
         name[0] == line_data[0]
     */
    inline
    bool
    is_name(const char* name, const char* line_data, unsigned line_size)
    {
      if (line_size > 5) {
        if (*(++line_data) != *(++name)) return false;
        if (*(++line_data) != *(++name)) return false;
        if (*(++line_data) != *(++name)) return false;
        if (*(++line_data) != *(++name)) return false;
        if (*(++line_data) != *(++name)) return false;
        return true;
      }
      else {
        unsigned i = 0;
        while (++i < line_size) {
          if (*(++line_data) != *(++name)) return false;
        }
        while (++i < 6) {
          if (*(++name) != ' ') return false;
        }
      }
      return true;
    }

    inline
    bool
    is_name(
      const char* name,
      const char* line_data,
      unsigned line_size,
      const char* digits)
    {
      if (line_size < 6) return false;
      unsigned i = 0;
      while (++i < 5) {
        if (*(++line_data) != *(++name)) return false;
      }
      line_data++;
      while (*digits != '\0') {
        if (*line_data == *digits++) return true;
      }
      return false;
    }

    struct info
    {
      record_type::section section;
      record_type::id id;

      info() {}

      info(
        const char* line_data,
        unsigned line_size)
      {
        if (line_size == 0) {
          set(unknown, xxxxxx);
          return;
        }
        switch (line_data[0]) {
          case 'A':
            if (is_name("ATOM  ", line_data, line_size)) {
              set(coordinate, atom__); return;
            }
            if (is_name("ANISOU", line_data, line_size)) {
              set(coordinate, anisou); return;
            }
            if (is_name("AUTHOR", line_data, line_size)) {
              set(title, author); return;
            }
            break;
          case 'H':
            if (is_name("HETATM", line_data, line_size)) {
              set(coordinate, hetatm); return;
            }
            if (is_name("HELIX ", line_data, line_size)) {
              set(secondary_structure, helix_); return;
            }
            if (is_name("HEADER", line_data, line_size)) {
              set(title, header); return;
            }
            if (is_name("HYDBND", line_data, line_size)) {
              set(connectivity_annotation, hydbnd); return;
            }
            if (is_name("HET   ", line_data, line_size)) {
              set(heterogen, het___); return;
            }
            if (is_name("HETNAM", line_data, line_size)) {
              set(heterogen, hetnam); return;
            }
            if (is_name("HETSYN", line_data, line_size)) {
              set(heterogen, hetsyn); return;
            }
            break;
          case 'S':
            if (is_name("SIGATM", line_data, line_size)) {
              set(coordinate, sigatm); return;
            }
            if (is_name("SIGUIJ", line_data, line_size)) {
              set(coordinate, siguij); return;
            }
            if (is_name("SHEET ", line_data, line_size)) {
              set(secondary_structure, sheet_); return;
            }
            if (is_name("SSBOND", line_data, line_size)) {
              set(connectivity_annotation, ssbond); return;
            }
            if (is_name("SLTBRG", line_data, line_size)) {
              set(connectivity_annotation, sltbrg); return;
            }
            if (is_name("SOURCE", line_data, line_size)) {
              set(title, source); return;
            }
            if (is_name("SPRSDE", line_data, line_size)) {
              set(title, sprsde); return;
            }
            if (is_name("SEQADV", line_data, line_size)) {
              set(primary_structure, seqadv); return;
            }
            if (is_name("SEQRES", line_data, line_size)) {
              set(primary_structure, seqres); return;
            }
            if (is_name("SCALE ", line_data, line_size, "123")) {
              set(crystallographic, scale_); return;
            }
            if (is_name("SITE  ", line_data, line_size)) {
              set(miscellaneous_features, site__); return;
            }
          case 'C':
            if (is_name("CONECT", line_data, line_size)) {
              set(connectivity, conect); return;
            }
            if (is_name("CISPEP", line_data, line_size)) {
              set(connectivity_annotation, cispep); return;
            }
            if (is_name("COMPND", line_data, line_size)) {
              set(title, compnd); return;
            }
            if (is_name("CAVEAT", line_data, line_size)) {
              set(title, caveat); return;
            }
            if (is_name("CRYST1", line_data, line_size)) {
              set(crystallographic, cryst1); return;
            }
            break;
          case 'R':
            if (is_name("REMARK", line_data, line_size)) {
              set(title, remark); return;
            }
            if (is_name("REVDAT", line_data, line_size)) {
              set(title, revdat); return;
            }
            break;
          case 'T':
            if (is_name("TER   ", line_data, line_size)) {
              set(coordinate, ter___); return;
            }
            if (is_name("TURN  ", line_data, line_size)) {
              set(secondary_structure, turn__); return;
            }
            if (is_name("TITLE ", line_data, line_size)) {
              set(title, title_); return;
            }
            if (is_name("TVECT ", line_data, line_size)) {
              set(crystallographic, tvect_); return;
            }
            break;
          case 'B':
            if (is_name("BREAK ", line_data, line_size)) {
              set(coordinate, break_); return;
            }
            break;
          case 'L':
            if (is_name("LINK  ", line_data, line_size)) {
              set(connectivity_annotation, link__); return;
            }
            break;
          case 'M':
            if (is_name("MODEL ", line_data, line_size)) {
              set(coordinate, model_); return;
            }
            if (is_name("MODRES", line_data, line_size)) {
              set(primary_structure, modres); return;
            }
            if (is_name("MTRIX ", line_data, line_size, "123")) {
              set(crystallographic, mtrix_); return;
            }
            if (is_name("MASTER", line_data, line_size)) {
              set(bookkeeping, master); return;
            }
            break;
          case 'E':
            if (is_name("ENDMDL", line_data, line_size)) {
              set(coordinate, endmdl); return;
            }
            if (is_name("EXPDTA", line_data, line_size)) {
              set(title, expdta); return;
            }
            if (is_name("END   ", line_data, line_size)) {
              set(bookkeeping, end___); return;
            }
            break;
          case 'F':
            if (is_name("FORMUL", line_data, line_size)) {
              set(heterogen, formul); return;
            }
            if (is_name("FTNOTE", line_data, line_size)) {
              set(title, ftnote); return;
            }
          case 'O':
            if (is_name("ORIGX ", line_data, line_size, "123")) {
              set(crystallographic, origx_); return;
            }
            if (is_name("ONHOLD", line_data, line_size)) {
              set(title, onhold); return;
            }
            if (is_name("OBSLTE", line_data, line_size)) {
              set(title, obslte); return;
            }
          case 'K':
            if (is_name("KEYWDS", line_data, line_size)) {
              set(title, keywds); return;
            }
          case 'J':
            if (is_name("JRNL  ", line_data, line_size)) {
              set(title, jrnl__); return;
            }
          case 'D':
            if (is_name("DBREF ", line_data, line_size)) {
              set(primary_structure, dbref_); return;
            }
          default:
            break;
        }
        set(unknown, xxxxxx);
      }

      void
      set(record_type::section section_, record_type::id id_)
      {
        section = section_;
        id = id_;
      }
    };

  } // namespace record_type

  inline
  int
  read_model_number(pdb::line_info& line_info)
  {
    unsigned i_col = 6;
    for(;i_col<line_info.size;i_col++) {
      if (line_info.data[i_col] != ' ') break;
    }
    bool ok = true;
    char buf[16];
    unsigned i_buf = 0;
    for(;i_col<line_info.size;i_col++) {
      if (line_info.data[i_col] == ' ') break;
      buf[i_buf++] = line_info.data[i_col];
      if (i_buf == sizeof buf) {
        ok = false;
        break;
      }
    }
    if (ok && i_buf != 0) {
      buf[i_buf] = '\0';
      char *endptr;
      long result = std::strtol(buf, &endptr, 10);
      if (endptr == buf + i_buf) return result;
    }
    return field_as_int(line_info,10,14);
  }

  struct input_atom_labels
  {
    int32_t resseq;
    static const unsigned compacted_size = 4+4+1+1+4+1;
    char compacted[compacted_size];

    char*       name_begin()       { return compacted; }
    const char* name_begin() const { return compacted; }
    str4        name_small() const { return str4(name_begin()); }
    std::string name()       const { return std::string(name_begin(),4); }

    char*       resname_begin()       { return compacted+4; }
    const char* resname_begin() const { return compacted+4; }
    str4        resname_small() const { return str4(resname_begin()); }
    std::string resname()       const { return std::string(resname_begin(),4);}

    char*       chain_begin()       { return compacted+8; }
    const char* chain_begin() const { return compacted+8; }
    str1        chain_small() const { return str1(chain_begin()); }
    std::string chain()       const { return std::string(chain_begin(),1); }

    char*       icode_begin()       { return compacted+9; }
    const char* icode_begin() const { return compacted+9; }
    str1        icode_small() const { return str1(icode_begin()); }
    std::string icode()       const { return std::string(icode_begin(),1); }

    char*       segid_begin()       { return compacted+10; }
    const char* segid_begin() const { return compacted+10; }
    str4        segid_small() const { return str4(segid_begin()); }
    std::string segid()       const { return std::string(segid_begin(),4); }

    char*       altloc_begin()       { return compacted+14; }
    const char* altloc_begin() const { return compacted+14; }
    str1        altloc_small() const { return str1(altloc_begin()); }
    std::string altloc()       const { return std::string(altloc_begin(),1); }

    input_atom_labels() {}

    input_atom_labels(pdb::line_info& line_info)
    :
      //  7 - 11  Integer       serial   Atom serial number.
      // 13 - 16  Atom          name     Atom name.
      // 17       Character     altLoc   Alternate location indicator.
      // 18 - 20  Residue name  resName  Residue name.
      // 21       *** non-standard: appended to resName
      // 22       Character     chainID  Chain identifier.
      // 23 - 26  Integer       resSeq   Residue sequence number.
      // 27       AChar         iCode    Code for insertion of residues.
      // 73 - 76  LString(4)    segID    Segment identifier, left-justified.
      resseq(field_as_int(line_info,22,26))
    {
      extract(line_info,12,4,name_begin());
      extract(line_info,16,1,altloc_begin());
      extract(line_info,17,4,resname_begin());
      extract(line_info,21,1,chain_begin());
      extract(line_info,26,1,icode_begin());
      extract(line_info,72,4,segid_begin());
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
    pdb_format() const
    {
      SCITBX_ASSERT(resseq >= min_resseq);
      SCITBX_ASSERT(resseq <= max_resseq);
      char buf[128];
      std::sprintf(buf, "\"%4.4s%1.1s%4.4s%1.1s%4d%1.1s\"",
        name_begin(), altloc_begin(), resname_begin(), chain_begin(),
        resseq, icode_begin());
      std::string result;
      result += buf;
      if (str4(segid_begin()).stripped_size() != 0) {
        result += " segid=\"" + segid() + "\"";
      }
      return result;
    }

    void
    check_equivalence(pdb::line_info& line_info) const
    {
      if      (!are_equal(line_info,12,4,name_begin())) {
        line_info.set_error(13, "name mismatch.");
      }
      else if (!are_equal(line_info,16,1,altloc_begin())) {
        line_info.set_error(17, "altloc mismatch.");
      }
      else if (!are_equal(line_info,17,4,resname_begin())) {
        line_info.set_error(18, "resname mismatch.");
      }
      else if (!are_equal(line_info,21,1,chain_begin())) {
        line_info.set_error(22, "chain mismatch.");
      }
      else if (resseq != field_as_int(line_info,22,26)) {
        line_info.set_error(23, "resseq mismatch.");
      }
      else if (!are_equal(line_info,26,1,icode_begin())) {
        line_info.set_error(27, "icode mismatch.");
      }
      else if (*chain_begin() == ' '
            && !are_equal(line_info,72,4,segid_begin())) {
        line_info.set_error(74, "segid mismatch.");
      }
    }

    int
    compare(input_atom_labels const& other) const
    {
      if (resseq < other.resseq) return -1;
      if (resseq > other.resseq) return  1;
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
      return 0;
    }

    bool
    operator!=(input_atom_labels const& other) const
    {
      if (resseq != other.resseq) return true;
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
      return false;
    }

    bool
    not_equal_ignoring_altloc(input_atom_labels const& other) const
    {
      if (resseq != other.resseq) return true;
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
      return false;
    }

    bool
    is_alternate(input_atom_labels const& other) const
    {
      char self_altloc = *altloc_begin();
      if (self_altloc != blank_altloc_char) return true;
      if (self_altloc == *other.altloc_begin()) return false;
      return !not_equal_ignoring_altloc(other);
    }
  };

  struct i_seq_input_atom_labels
  {
    unsigned i_seq;
    input_atom_labels labels;

    bool
    operator<(i_seq_input_atom_labels const& other) const
    {
      switch (labels.compare(other.labels))
      {
        case -1: return true;
        case  1: return false;
        default: break;
      }
      return (i_seq < other.i_seq);
    }
  };

  scitbx::auto_array<i_seq_input_atom_labels>
  sort_input_atom_labels(
    const pdb::input_atom_labels* input_atom_labels_list,
    unsigned i_begin,
    unsigned i_end)
  {
    scitbx::auto_array<i_seq_input_atom_labels>
      result(new i_seq_input_atom_labels[i_end-i_begin]);
    i_seq_input_atom_labels* r = result.get();
    const pdb::input_atom_labels* ial = input_atom_labels_list + i_begin;
    for(unsigned i_seq=i_begin;i_seq<i_end;i_seq++) {
      r->i_seq = i_seq;
      (r++)->labels = *ial++;
    }
    r = result.get();
    std::sort(r, &r[i_end-i_begin]);
    return result;
  }

  atom*
  process_atom_record(pdb::line_info& line_info)
  {
    // 13 - 16  Atom          name          Atom name.
    // 31 - 38  Real(8.3)     x             Orthogonal coordinates for X in
    //                                      Angstroms.
    // 39 - 46  Real(8.3)     y             Orthogonal coordinates for Y in
    //                                      Angstroms.
    // 47 - 54  Real(8.3)     z             Orthogonal coordinates for Z in
    //                                      Angstroms.
    // 55 - 60  Real(6.2)     occupancy     Occupancy.
    // 61 - 66  Real(6.2)     tempFactor    Temperature factor.
    // 73 - 76  LString(4)    segID         Segment identifier, left-justified.
    // 77 - 78  LString(2)    element       Element symbol, right-justified.
    // 79 - 80  LString(2)    charge        Charge on the atom.
    return new atom(
      str4(line_info.data,line_info.size,12,' '), // name
      str4(line_info.data,line_info.size,72,' '), // segid
      str2(line_info.data,line_info.size,76,' '), // element
      str2(line_info.data,line_info.size,78,' '), // charge
      vec3(
        field_as_double(line_info,30,38),
        field_as_double(line_info,38,46),
        field_as_double(line_info,46,54)), // xyz
      vec3(0,0,0), // sigxyz
      field_as_double(line_info,54,60), // occ
      0, // sigocc
      field_as_double(line_info,60,66), // b
      0, // sigb
      sym_mat3(-1,-1,-1,-1,-1,-1), // uij
      sym_mat3(-1,-1,-1,-1,-1,-1)); // siguij
  }

  void
  process_sigatm_record(
    pdb::line_info& line_info,
    pdb::input_atom_labels const& input_atom_labels,
    pdb::atom& atom)
  {
    atom.sigxyz = vec3(
      field_as_double(line_info,30,38),
      field_as_double(line_info,38,46),
      field_as_double(line_info,46,54));
    atom.sigocc = field_as_double(line_info,54,60);
    atom.sigb = field_as_double(line_info,60,66);
    input_atom_labels.check_equivalence(line_info);
  }

  void
  process_anisou_record(
    pdb::line_info& line_info,
    pdb::input_atom_labels const& input_atom_labels,
    pdb::atom& atom)
  {
    // 29 - 35  Integer       U(1,1)
    // 36 - 42  Integer       U(2,2)
    // 43 - 49  Integer       U(3,3)
    // 50 - 56  Integer       U(1,2)
    // 57 - 63  Integer       U(1,3)
    // 64 - 70  Integer       U(2,3)
    atom.uij[0] = field_as_int(line_info,28,35);
    atom.uij[1] = field_as_int(line_info,35,42);
    atom.uij[2] = field_as_int(line_info,42,49);
    atom.uij[3] = field_as_int(line_info,49,56);
    atom.uij[4] = field_as_int(line_info,56,63);
    atom.uij[5] = field_as_int(line_info,63,70);
    atom.uij *= anisou_factor;
    input_atom_labels.check_equivalence(line_info);
  }

  void
  process_siguij_record(
    pdb::line_info& line_info,
    pdb::input_atom_labels const& input_atom_labels,
    pdb::atom& atom)
  {
    atom.siguij[0] = field_as_int(line_info,28,35);
    atom.siguij[1] = field_as_int(line_info,35,42);
    atom.siguij[2] = field_as_int(line_info,42,49);
    atom.siguij[3] = field_as_int(line_info,49,56);
    atom.siguij[4] = field_as_int(line_info,56,63);
    atom.siguij[5] = field_as_int(line_info,63,70);
    atom.siguij *= anisou_factor;
    input_atom_labels.check_equivalence(line_info);
  }

  class model_record_oversight
  {
    private:
      pdb::line_info& line_info;
      enum styles { unknown, bare, encapsulated };
      int style;
      bool active_block;

    public:
      model_record_oversight(pdb::line_info& line_info_)
      :
        line_info(line_info_),
        style(unknown),
        active_block(false)
      {}

      bool
      atom_is_allowed_here()
      {
        if (style == bare || active_block) return true;
        if (style == unknown) {
          style = bare;
          return true;
        }
        line_info.set_error(1,
          "ATOM or HETATM record is outside MODEL/ENDMDL block.");
        return false;
      }

      bool
      model_is_allowed_here()
      {
        if (style == encapsulated) {
          if (active_block) {
            line_info.set_error(1,
              "Missing ENDMDL for previous MODEL record.");
            active_block = false;
            return false;
          }
          active_block = true;
          return true;
        }
        if (style == unknown) {
          style = encapsulated;
          active_block = true;
          return true;
        }
        line_info.set_error(1,
          "MODEL record must appear before any ATOM or HETATM records.");
        return false;
      }

      void
      check_if_endmdl_is_allowed_here()
      {
        if (style != encapsulated || !active_block) {
          line_info.set_error(1,
            "No matching MODEL record.");
        }
        active_block = false;
      }

      bool
      expecting_endmdl_record() const { return active_block; }
  };

  struct model_range_loop
  {
    model_range_loop() {}

    model_range_loop(
      af::const_ref<std::size_t> const& model_indices)
    :
      mii_end(model_indices.end()),
      mii(model_indices.begin()),
      end(0)
    {}

    bool
    next()
    {
      if (mii == mii_end) return false;
      begin = end;
      end = static_cast<unsigned>(*mii++);
      size = end - begin;
      return true;
    }

    protected:
      const std::size_t* mii_end;
      const std::size_t* mii;
    public:
      unsigned begin;
      unsigned end;
      unsigned size;
  };

  class conformer_selections
  {
    public:
      std::vector<unsigned> common;
      std::map<std::string, std::vector<unsigned> > alternates;
  };

  class input
  {
    public:
      typedef std::map<str6, unsigned> record_type_counts_t;
      typedef std::map<std::string, std::vector<unsigned> > str_sel_cache_t;
      typedef af::shared<std::vector<unsigned> > int_sel_cache_t;

      input()
      :
        name_selection_cache_is_up_to_date_(false),
        altloc_selection_cache_is_up_to_date_(false),
        resname_selection_cache_is_up_to_date_(false),
        chain_selection_cache_is_up_to_date_(false),
        resseq_selection_cache_is_up_to_date_(false),
        icode_selection_cache_is_up_to_date_(false),
        segid_selection_cache_is_up_to_date_(false)
      {}

      input&
      process(
        const char* source_info,
        af::const_ref<std::string> const& lines)
      {
        columns_73_76_evaluator columns_73_76_eval(lines);
        input_atom_labels_list_.reserve(
            input_atom_labels_list_.size()
          + columns_73_76_eval.number_of_atom_and_hetatm_lines);
        atoms_.reserve(
            atoms_.size()
          + columns_73_76_eval.number_of_atom_and_hetatm_lines);
        input_atom_labels* current_input_atom_labels = 0;
        atom* current_atom = 0;
        bool expect_anisou;
        bool expect_sigatm;
        bool expect_siguij;
        const char* error_message_no_matching_atom
          = "no matching ATOM or HETATM record.";
        unsigned atom___counts = 0;
        unsigned hetatm_counts = 0;
        unsigned sigatm_counts = 0;
        unsigned anisou_counts = 0;
        unsigned siguij_counts = 0;
        pdb::line_info line_info(source_info);
        pdb::model_record_oversight model_record_oversight(line_info);
        for(unsigned i_line=0;i_line<lines.size();i_line++) {
          std::string const& line = lines[i_line];
          line_info.line_number++;
          line_info.size = static_cast<unsigned>(line.size());
          if (   columns_73_76_eval.is_old_style
              && line_info.size > 72) {
            line_info.size = 72;
          }
          line_info.data = line.data();
          record_type::info record_type_info(line_info.data, line_info.size);
          if      (record_type_info.id == record_type::atom__) {
            if (model_record_oversight.atom_is_allowed_here()) {
              input_atom_labels_list_.push_back(input_atom_labels(line_info));
              atoms_.push_back(atom_shptr(process_atom_record(line_info)));
              if (line_info.error_occured()) {
                input_atom_labels_list_.pop_back();
                atoms_.pop_back();
              }
              else {
                atom___counts++;
                current_input_atom_labels = &input_atom_labels_list_.back();
                current_atom = atoms_.back().get();
                expect_anisou = true;
                expect_sigatm = true;
                expect_siguij = true;
              }
            }
          }
          else if (record_type_info.id == record_type::hetatm) {
            if (model_record_oversight.atom_is_allowed_here()) {
              input_atom_labels_list_.push_back(input_atom_labels(line_info));
              atoms_.push_back(atom_shptr(process_atom_record(line_info)));
              if (line_info.error_occured()) {
                input_atom_labels_list_.pop_back();
                atoms_.pop_back();
              }
              else {
                hetatm_counts++;
                current_input_atom_labels = &input_atom_labels_list_.back();
                current_atom = atoms_.back().get();
                expect_anisou = true;
                expect_sigatm = true;
                expect_siguij = true;
              }
            }
          }
          else if (record_type_info.id == record_type::anisou) {
            if (!expect_anisou) {
              line_info.set_error(1, error_message_no_matching_atom);
            }
            expect_anisou = false;
            anisou_counts++;
            process_anisou_record(
              line_info, *current_input_atom_labels, *current_atom);
          }
          else if (record_type_info.id == record_type::sigatm) {
            if (!expect_sigatm) {
              line_info.set_error(1, error_message_no_matching_atom);
            }
            expect_sigatm = false;
            sigatm_counts++;
            process_sigatm_record(
              line_info, *current_input_atom_labels, *current_atom);
          }
          else if (record_type_info.id == record_type::siguij) {
            if (!expect_siguij) {
              line_info.set_error(1, error_message_no_matching_atom);
            }
            expect_siguij = false;
            siguij_counts++;
            process_siguij_record(
              line_info, *current_input_atom_labels, *current_atom);
          }
          else {
            record_type_counts_[
              str6(line_info.data, line_info.size, 0, ' ')]++;
            current_input_atom_labels = 0;
            current_atom = 0;
            expect_anisou = false;
            expect_sigatm = false;
            expect_siguij = false;
            if      (   record_type_info.id == record_type::remark
                     || record_type_info.id == record_type::ftnote) {
              remark_section_.push_back(line_info.strip_data());
            }
            else if (record_type_info.id == record_type::ter___) {
              ter_indices_.push_back(input_atom_labels_list_.size());
            }
            else if (record_type_info.id == record_type::break_) {
              break_indices_.push_back(input_atom_labels_list_.size());
            }
            else if (record_type_info.id == record_type::model_) {
              if (model_record_oversight.model_is_allowed_here()) {
                if (input_atom_labels_list_.size() != 0) {
                  model_indices_.push_back(input_atom_labels_list_.size());
                }
                model_numbers_.push_back(read_model_number(line_info));
              }
            }
            else if (record_type_info.id == record_type::endmdl) {
              model_record_oversight.check_if_endmdl_is_allowed_here();
            }
            else if (record_type_info.section == record_type::title) {
              title_section_.push_back(line_info.strip_data());
            }
            else if (record_type_info.section
                     == record_type::primary_structure) {
              primary_structure_section_.push_back(line_info.strip_data());
            }
            else if (record_type_info.section
                     == record_type::heterogen) {
              heterogen_section_.push_back(line_info.strip_data());
            }
            else if (record_type_info.section
                     == record_type::secondary_structure) {
              secondary_structure_section_.push_back(line_info.strip_data());
            }
            else if (record_type_info.section
                     == record_type::connectivity_annotation) {
              connectivity_annotation_section_.push_back(
                line_info.strip_data());
            }
            else if (record_type_info.section
                     == record_type::miscellaneous_features) {
              miscellaneous_features_section_.push_back(
                line_info.strip_data());
            }
            else if (record_type_info.section
                     == record_type::crystallographic) {
              crystallographic_section_.push_back(line_info.strip_data());
            }
            else if (record_type_info.section
                     == record_type::connectivity) {
              connectivity_section_.push_back(line_info.strip_data());
            }
            else if (record_type_info.section
                     == record_type::bookkeeping) {
              bookkeeping_section_.push_back(line_info.strip_data());
            }
            else if (!line_info.is_blank_line()) {
              unknown_section_.push_back(line_info.strip_data());
            }
          }
          line_info.check_and_throw_runtime_error();
        }
        if (atom___counts != 0) record_type_counts_["ATOM  "] += atom___counts;
        if (hetatm_counts != 0) record_type_counts_["HETATM"] += hetatm_counts;
        if (sigatm_counts != 0) record_type_counts_["SIGATM"] += sigatm_counts;
        if (anisou_counts != 0) record_type_counts_["ANISOU"] += anisou_counts;
        if (siguij_counts != 0) record_type_counts_["SIGUIJ"] += siguij_counts;
        if (model_indices_.size() != model_numbers_.size()) {
          model_indices_.push_back(input_atom_labels_list_.size());
        }
        else if (   model_indices_.size() == 0
                 || model_indices_.back() != input_atom_labels_list_.size()) {
          if (input_atom_labels_list_.size() != 0) {
            model_numbers_.push_back(0);
            model_indices_.push_back(input_atom_labels_list_.size());
          }
        }
        if (model_record_oversight.expecting_endmdl_record()) {
          throw std::runtime_error("ENDMDL record missing at end of input.");
        }
        SCITBX_ASSERT(model_indices_.size() == model_numbers_.size());
        return *this;
      }

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

      af::shared<atom_shptr> const&
      atoms() const { return atoms_; }

      af::shared<int> const&
      model_numbers() const { return model_numbers_; }

      af::shared<std::size_t> const&
      model_indices() const { return model_indices_; }

      af::shared<std::size_t> const&
      ter_indices() const { return ter_indices_; }

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
      IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(resname, str4)
      IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(chain, str1)

      int_sel_cache_t const&
      resseq_selection_cache() const
      {
        if (!resseq_selection_cache_is_up_to_date_) {
          resseq_selection_cache_.resize(max_resseq-min_resseq+1);
          unsigned i_seq = 0;
          const input_atom_labels* iall_end = input_atom_labels_list_.end();
          for(const input_atom_labels* ial=input_atom_labels_list_.begin();
                                       ial!=iall_end;ial++) {
            int resseq = ial->resseq;
            SCITBX_ASSERT(resseq >= min_resseq);
            SCITBX_ASSERT(resseq <= max_resseq);
            resseq_selection_cache_[resseq-min_resseq].push_back(i_seq++);
          }
          resseq_selection_cache_is_up_to_date_ = true;
        }
        return resseq_selection_cache_;
      }

      IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(icode, str1)
      IOTBX_PDB_INPUT_SELECTION_STR_SEL_CACHE_MEMBER_FUNCTION(segid, str4)

      bool
      model_numbers_are_unique() const
      {
        std::set<int> unique_numbers(
          model_numbers_.begin(), model_numbers_.end());
        return (unique_numbers.size() == model_numbers_.size());
      }

      af::shared<std::size_t>
      model_atom_counts() const
      {
        af::shared<std::size_t> result(af::reserve(model_indices_.size()));
        model_range_loop mr(model_indices_.const_ref());
        while (mr.next()) result.push_back(mr.size);
        return result;
      }

      af::shared<std::vector<unsigned> >
      find_duplicate_atom_labels() const
      {
        af::shared<std::vector<unsigned> > result;
        const input_atom_labels* iall = input_atom_labels_list_.begin();
        model_range_loop mr(model_indices_.const_ref());
        while (mr.next()) {
          if (mr.size < 2) continue;
          std::map<unsigned, std::vector<unsigned> > dup_map;
          scitbx::auto_array<i_seq_input_atom_labels>
            sorted_input_atom_labels = sort_input_atom_labels(
              iall, mr.begin, mr.end);
          const i_seq_input_atom_labels* sj = sorted_input_atom_labels.get();
          const i_seq_input_atom_labels* si = sj++;
          for(unsigned js=1;js<mr.size;js++) {
            if (sj->labels != si->labels) {
              si = sj++;
            }
            else {
              dup_map[si->i_seq].push_back((sj++)->i_seq);
            }
          }
          for(std::map<unsigned, std::vector<unsigned> >::iterator
                di=dup_map.begin(); di!=dup_map.end(); di++) {
            di->second.push_back(di->first);
            std::sort(di->second.begin(), di->second.end());
            result.push_back(std::vector<unsigned>());
            result.back().swap(di->second); // swap avoids deep copy
          }
        }
        return result;
      }

      af::small<unsigned, 2>
      find_false_blank_altloc() const
      {
        str_sel_cache_t const& alsc = altloc_selection_cache();
        if (alsc.size() < 2 || alsc.find(blank_altloc_stl) == alsc.end()) {
          return af::small<unsigned, 2>();
        }
        const input_atom_labels* iall = input_atom_labels_list_.begin();
        model_range_loop mr(model_indices_.const_ref());
        while (mr.next()) {
          if (mr.size < 2) continue;
          scitbx::auto_array<i_seq_input_atom_labels>
            sorted_input_atom_labels(
              sort_input_atom_labels(iall, mr.begin, mr.end));
          const i_seq_input_atom_labels* sj = sorted_input_atom_labels.get();
          const i_seq_input_atom_labels* si = sj++;
          for(unsigned js=1;js<mr.size;js++) {
            if (sj->labels.not_equal_ignoring_altloc(si->labels)) {
              si = sj++;
              continue;
            }
            if (   (*(sj->labels.altloc_begin()) == blank_altloc_char)
                != (*(si->labels.altloc_begin()) == blank_altloc_char)) {
              af::small<unsigned, 2> result;
              result.push_back(si->i_seq);
              result.push_back(sj->i_seq);
              return result;
            }
            sj++;
          }
        }
        return af::small<unsigned, 2>();
      }

      af::shared<pdb::conformer_selections>
      conformer_selections() const
      {
        af::shared<pdb::conformer_selections>
          result((af::reserve(model_indices_.size())));
        const input_atom_labels* iall = input_atom_labels_list_.begin();
        model_range_loop mr(model_indices_.const_ref());
        while (mr.next()) {
          std::vector<unsigned> common;
          std::map<str1, std::vector<unsigned> > alternates;
          scitbx::auto_array<i_seq_input_atom_labels>
            sorted_input_atom_labels(
              sort_input_atom_labels(iall, mr.begin, mr.end));
          const i_seq_input_atom_labels* sj = sorted_input_atom_labels.get();
          const i_seq_input_atom_labels* si = (mr.size == 1 ? sj : sj++);
          if (mr.size != 0) {
            if (sj->labels.is_alternate(si->labels)) {
              alternates[sj->labels.altloc_small()].push_back(sj->i_seq);
            }
            else {
              common.push_back(sj->i_seq);
            }
            for(unsigned js=1;js<mr.size;js++,si=sj++) {
              if (sj->labels.is_alternate(si->labels)) {
                alternates[sj->labels.altloc_small()].push_back(sj->i_seq);
              }
              else {
                common.push_back(sj->i_seq);
              }
            }
          }
          result.push_back(pdb::conformer_selections());
          pdb::conformer_selections& result_back = result.back();
          result_back.common.assign(common.begin(), common.end());
          for(std::map<str1, std::vector<unsigned> >::iterator
                ai=alternates.begin(); ai!=alternates.end(); ai++) {
            result_back.alternates[std::string(ai->first.elems)]
              .assign(ai->second.begin(), ai->second.end());
          }
        }
        return result;
      }

    protected:
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
      af::shared<atom_shptr>  atoms_;
      af::shared<int>         model_numbers_;
      af::shared<std::size_t> model_indices_;
      af::shared<std::size_t> ter_indices_;
      af::shared<std::size_t> break_indices_;
      af::shared<std::string> connectivity_section_;
      af::shared<std::string> bookkeeping_section_;
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
      mutable int_sel_cache_t resseq_selection_cache_;
      mutable bool            icode_selection_cache_is_up_to_date_;
      mutable str_sel_cache_t icode_selection_cache_;
      mutable bool            segid_selection_cache_is_up_to_date_;
      mutable str_sel_cache_t segid_selection_cache_;
  };

}} // namespace iotbx::pdb

#endif // IOTBX_PDB_INPUT_H
