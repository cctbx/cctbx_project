#include <iotbx/pdb/input.h>
#include <iotbx/error.h>
#include <scitbx/misc/file_utils.h>

namespace iotbx { namespace pdb {

  std::string
  line_info::format_exception_message() const
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

  //! Fortran-style conversion of an integer input field.
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

  //! Fortran-style conversion of a floating-point input field.
  /*! i_end-i_begin must be less than 32. This is not checked to
      maximize performance.
   */
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
        if (c == 'x' || c == 'X' || c == 'n' || c == 'N') {
          // suppress conversion of hexadecimal literals and "nan"
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

  columns_73_76_evaluator::columns_73_76_evaluator(
    af::const_ref<std::string> const& lines,
    unsigned is_frequent_threshold_atom_records,
    unsigned is_frequent_threshold_other_records)
  :
    finding("Undecided."),
    is_old_style(false),
    number_of_atom_and_hetatm_lines(0)
  {
    columns_73_76_dict_t atom_columns_73_76_dict;
    columns_73_76_dict_t other_columns_73_76_dict;
    af::tiny<char, 4> header_idcode(' ', ' ', ' ', ' ');
    for(const std::string* line=lines.begin();line!=lines.end();line++) {
      if (line->size() < 6) continue;
      const char* line_data = line->data();
      bool is_atom_or_hetatm_record = (
           is_record_type("ATOM  ", line_data)
        || is_record_type("HETATM", line_data));
      if (is_atom_or_hetatm_record) {
        number_of_atom_and_hetatm_lines++;
      }
      if (is_record_type("HEADER", line_data)) {
        if (line->size() < 66) continue;
        const char* l = line_data + 62;
        char* h = header_idcode.begin();
        *h++ = *l++;
        *h++ = *l++;
        *h++ = *l++;
        *h   = *l  ;
      }
      if (line->size() < 80) continue;
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
      if (is_atom_or_hetatm_record) {
        atom_columns_73_76_dict[columns_73_76]++;
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
    if (atom_columns_73_76_dict.size() == 1) {
      columns_73_76_t const& a = atom_columns_73_76_dict.begin()->first;
      if (a[0] != ' ' && a[3] != ' ') {
        if (other_columns_73_76_dict.size() == 0) {
          if (a.all_eq(header_idcode)) {
            is_old_style = true;
          }
        }
        if (other_columns_73_76_dict.size() == 1) {
          if (a.all_eq(other_columns_73_76_dict.begin()->first)) {
            is_old_style = true;
          }
        }
        if (is_old_style) {
          finding = "Exactly one common label in columns 73-76.";
          return;
        }
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

  //! Fast processing of record types.
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
      split_, // added in PDB V3.2
      caveat,
      compnd,
      source,
      keywds,
      expdta,
      nummdl, // added in PDB V3.2
      mdltyp, // added in PDB V3.2
      author,
      revdat,
      sprsde,
      jrnl__,
      remark,
      ftnote, // non-standard
      // Primary Structure Section
      dbref_,
      dbref1, // added in PDB V3.2
      dbref2, // added in PDB V3.2
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
            if (line_size >= 4) {
              // "ATOM  ": ignore columns 5 & 6 for compatibility with CNS 1.2
              if (   line_data[1] == 'T'
                  && line_data[2] == 'O'
                  && line_data[3] == 'M') {
                set(coordinate, atom__); return;
              }
              if (is_name("ANISOU", line_data, line_size)) {
                set(coordinate, anisou); return;
              }
              if (is_name("AUTHOR", line_data, line_size)) {
                set(title, author); return;
              }
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
            if (is_name("SPLIT ", line_data, line_size)) {
              set(title, split_); return;
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
            break;
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
            if (is_name("MDLTYP", line_data, line_size)) {
              set(title, mdltyp); return;
            }
            if (is_name("MASTER", line_data, line_size)) {
              set(bookkeeping, master); return;
            }
            break;
          case 'N':
            if (is_name("NUMMDL", line_data, line_size)) {
              set(title, nummdl); return;
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
            break;
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
            break;
          case 'K':
            if (is_name("KEYWDS", line_data, line_size)) {
              set(title, keywds); return;
            }
            break;
          case 'J':
            if (is_name("JRNL  ", line_data, line_size)) {
              set(title, jrnl__); return;
            }
            break;
          case 'D':
            if (is_name("DBREF ", line_data, line_size)) {
              set(primary_structure, dbref_); return;
            }
            if (is_name("DBREF1", line_data, line_size)) {
              set(primary_structure, dbref1); return;
            }
            if (is_name("DBREF2", line_data, line_size)) {
              set(primary_structure, dbref2); return;
            }
            break;
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

  //! Tolerant processing of MODEL records.
  str8
  read_model_id(pdb::line_info& line_info)
  {
    char blank = ' ';
    unsigned i_col = 6;
    for(;i_col<line_info.size&&i_col<10U;i_col++) {
      if (line_info.data[i_col] != blank) break;
    }
    str8 result;
    unsigned i = 0;
    for(;i_col<line_info.size&&i_col<14U;i_col++) {
      result.elems[i++] = line_info.data[i_col];
    }
    if (i < 4U) {
      unsigned n = 4U - i;
      std::memmove(result.elems+n, result.elems, i);
      std::fill_n(result.elems, n, blank);
      i = 4U;
    }
    result.elems[i] = '\0';
    return result;
  }

namespace detail {

  void
  input_atom_labels::check_equivalence(pdb::line_info& line_info) const
  {
    if      (!are_equal(line_info,12,4,name_begin())) {
      line_info.set_error(13, "name mismatch.");
    }
    else if (!are_equal(line_info,16,1,altloc_begin())) {
      line_info.set_error(17, "altloc mismatch.");
    }
    else if (!are_equal(line_info,17,3,resname_begin())) {
      line_info.set_error(18, "resname mismatch.");
    }
    else if (!are_equal(line_info,20,2,chain_begin())) {
      line_info.set_error(21, "chain mismatch.");
    }
    else if (!are_equal(line_info,22,4,resseq_begin())) {
      line_info.set_error(23, "resseq mismatch.");
    }
    else if (!are_equal(line_info,26,1,icode_begin())) {
      line_info.set_error(27, "icode mismatch.");
    }
    else if (   chain_begin()[1] == ' '
             && !are_equal(line_info,72,4,segid_begin())) {
      line_info.set_error(74, "segid mismatch.");
    }
  }

} // namespace detail

  hierarchy::atom
  process_atom_record(pdb::line_info& line_info, bool hetero)
  {
    // 13 - 16  Atom          name          Atom name.
    // 17       Character     altLoc        Alternate location indicator.
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
    return hierarchy::atom(
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
      sym_mat3(-1,-1,-1,-1,-1,-1), // siguij
      hetero,
      str5(line_info.data,line_info.size, 6,' '), // serial
      str4(line_info.data,line_info.size,12,' '), // name
      str4(line_info.data,line_info.size,72,' '), // segid
      str2(line_info.data,line_info.size,76,' '), // element
      str2(line_info.data,line_info.size,78,' ')); // charge
  }

  void
  process_sigatm_record(
    pdb::line_info& line_info,
    detail::input_atom_labels const& input_atom_labels,
    hierarchy::atom_data& atom_data)
  {
    atom_data.sigxyz = vec3(
      field_as_double(line_info,30,38),
      field_as_double(line_info,38,46),
      field_as_double(line_info,46,54));
    atom_data.sigocc = field_as_double(line_info,54,60);
    atom_data.sigb = field_as_double(line_info,60,66);
    input_atom_labels.check_equivalence(line_info);
  }

  void
  process_anisou_record(
    pdb::line_info& line_info,
    detail::input_atom_labels const& input_atom_labels,
    hierarchy::atom_data& atom_data)
  {
    // 29 - 35  Integer       U(1,1)
    // 36 - 42  Integer       U(2,2)
    // 43 - 49  Integer       U(3,3)
    // 50 - 56  Integer       U(1,2)
    // 57 - 63  Integer       U(1,3)
    // 64 - 70  Integer       U(2,3)
    atom_data.uij[0] = field_as_int(line_info,28,35);
    atom_data.uij[1] = field_as_int(line_info,35,42);
    atom_data.uij[2] = field_as_int(line_info,42,49);
    atom_data.uij[3] = field_as_int(line_info,49,56);
    atom_data.uij[4] = field_as_int(line_info,56,63);
    atom_data.uij[5] = field_as_int(line_info,63,70);
    atom_data.uij *= anisou_factor;
    input_atom_labels.check_equivalence(line_info);
  }

  void
  process_siguij_record(
    pdb::line_info& line_info,
    detail::input_atom_labels const& input_atom_labels,
    hierarchy::atom_data&
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
                             atom_data
#endif
                                      )
  {
#ifdef IOTBX_PDB_ENABLE_ATOM_DATA_SIGUIJ
    atom_data.siguij[0] = field_as_int(line_info,28,35);
    atom_data.siguij[1] = field_as_int(line_info,35,42);
    atom_data.siguij[2] = field_as_int(line_info,42,49);
    atom_data.siguij[3] = field_as_int(line_info,49,56);
    atom_data.siguij[4] = field_as_int(line_info,56,63);
    atom_data.siguij[5] = field_as_int(line_info,63,70);
    atom_data.siguij *= anisou_factor;
#endif
    input_atom_labels.check_equivalence(line_info);
  }

  //! Processing of MODEL & ENDMDL records.
  class model_record_oversight
  {
    protected:
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

      bool
      endmdl_is_allowed_here()
      {
        if (style != encapsulated || !active_block) {
          line_info.set_error(1,
            "No matching MODEL record.");
          return false;
        }
        active_block = false;
        return true;
      }

      bool
      expecting_endmdl_record() const { return active_block; }
  };

  //! Detection of chain boundaries.
  class chain_tracker
  {
    protected:
      af::shared<std::vector<unsigned> > chain_indices;
      std::vector<unsigned>* current_chain_indices;
      std::vector<unsigned> current_chain_indices_ignoring_segid;
      unsigned n_atoms;
      char previous_chain_and_segid[2+4];
      std::vector<str4> unique_segids;

      void
      evaluate_unique_segids()
      {
        // if unique_segids contains duplicates ignore segid
        std::set<str4> segid_set;
        std::vector<str4>::const_iterator us_end = unique_segids.end();
        for(std::vector<str4>::const_iterator
              us=unique_segids.begin();us!=us_end;us++) {
          if (!segid_set.insert(*us).second) {
            // we have a duplicate
            current_chain_indices->swap(current_chain_indices_ignoring_segid);
            break;
          }
        }
        current_chain_indices = 0;
        current_chain_indices_ignoring_segid.clear();
        unique_segids.clear();
      }

    public:
      chain_tracker()
      :
        current_chain_indices(0),
        n_atoms(0)
      {
        previous_chain_and_segid[0] = '\n';
      }

      void
      next_atom_labels(const detail::input_atom_labels& current_labels)
      {
        if (current_chain_indices == 0) {
          chain_indices.push_back(std::vector<unsigned>());
          current_chain_indices = &chain_indices.back();
        }
        char* p = previous_chain_and_segid;
        if (*p != '\n') {
          // test if previous labels and current labels belong to
          // different chains:
          //   - change in chainid
          //   - if common chainid is blank: change in segid
          if (   p[0] != current_labels.chain_begin()[0]
              || p[1] != current_labels.chain_begin()[1]) {
            current_chain_indices->push_back(n_atoms);
            current_chain_indices_ignoring_segid.push_back(n_atoms);
          }
          else if (p[1] == ' ') {
            ++p;
            const char* c = current_labels.segid_begin();
            if (   (*++p != *  c)
                || (*++p != *++c)
                || (*++p != *++c)
                || (*++p != *++c)) {
              current_chain_indices->push_back(n_atoms);
            }
            p = previous_chain_and_segid;
          }
        }
        // copy chain and segid
        *p++ = current_labels.chain_begin()[0];
        *p++ = current_labels.chain_begin()[1];
        const char* c = current_labels.segid_begin();
        *p++ = *c++;
        *p++ = *c++;
        *p++ = *c++;
        *p   = *c;
        // keep track of unique segids
        if (unique_segids.size() == 0) {
          unique_segids.push_back(current_labels.segid_small());
        }
        else {
          c = current_labels.segid_begin();
          p = unique_segids.back().elems;
          if (   *p++ != *c++
              || *p++ != *c++
              || *p++ != *c++
              || *p   != *c) {
            unique_segids.push_back(current_labels.segid_small());
          }
        }
        n_atoms++;
      }

      void
      transition()
      {
        if (*previous_chain_and_segid != '\n') {
          if (current_chain_indices != 0) {
            current_chain_indices->push_back(n_atoms);
            current_chain_indices_ignoring_segid.push_back(n_atoms);
          }
          *previous_chain_and_segid = '\n';
        }
      }

      void
      endmdl()
      {
        transition();
        if (current_chain_indices == 0) {
          chain_indices.push_back(std::vector<unsigned>());
        }
        else {
          evaluate_unique_segids();
        }
      }

      af::shared<std::vector<unsigned> >
      finish()
      {
        transition();
        if (current_chain_indices != 0) {
          evaluate_unique_segids();
        }
        return chain_indices;
      }
  };

  input::input(std::string const& file_name)
  :
    source_info_("file " + file_name)
  {
    process(scitbx::misc::file_to_lines(file_name).const_ref());
  }

  input::input(
    const char* source_info,
    af::const_ref<std::string> const& lines)
  :
    source_info_(source_info ? source_info : "")
  {
    process(lines);
  }

  void
  input::process(af::const_ref<std::string> const& lines)
  {
    columns_73_76_evaluator columns_73_76_eval(lines);
    input_atom_labels_list_.reserve(
        input_atom_labels_list_.size()
      + columns_73_76_eval.number_of_atom_and_hetatm_lines);
    atoms_.reserve(
        atoms_.size()
      + columns_73_76_eval.number_of_atom_and_hetatm_lines);
    detail::input_atom_labels* current_input_atom_labels = 0;
    hierarchy::atom_data* current_atom_data = 0;
    bool expect_anisou = false;
    bool expect_sigatm = false;
    bool expect_siguij = false;
    const char* error_message_no_matching_atom
      = "no matching ATOM or HETATM record.";
    unsigned atom___counts = 0;
    unsigned hetatm_counts = 0;
    unsigned sigatm_counts = 0;
    unsigned anisou_counts = 0;
    unsigned siguij_counts = 0;
    pdb::line_info line_info(source_info_.c_str());
    pdb::model_record_oversight model_record_oversight(line_info);
    pdb::chain_tracker chain_tracker;
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
          input_atom_labels_list_.push_back(
            detail::input_atom_labels(line_info));
          atoms_.push_back(process_atom_record(line_info,/*hetero*/false));
          if (line_info.error_occured()) {
            input_atom_labels_list_.pop_back();
            atoms_.pop_back();
          }
          else {
            atom___counts++;
            current_input_atom_labels = &input_atom_labels_list_.back();
            chain_tracker.next_atom_labels(*current_input_atom_labels);
            current_atom_data = atoms_.back().data.get();
            expect_anisou = true;
            expect_sigatm = true;
            expect_siguij = true;
          }
        }
      }
      else if (record_type_info.id == record_type::hetatm) {
        if (model_record_oversight.atom_is_allowed_here()) {
          input_atom_labels_list_.push_back(
            detail::input_atom_labels(line_info));
          atoms_.push_back(process_atom_record(line_info,/*hetero*/true));
          if (line_info.error_occured()) {
            input_atom_labels_list_.pop_back();
            atoms_.pop_back();
          }
          else {
            hetatm_counts++;
            current_input_atom_labels = &input_atom_labels_list_.back();
            chain_tracker.next_atom_labels(*current_input_atom_labels);
            current_atom_data = atoms_.back().data.get();
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
        else {
          expect_anisou = false;
          anisou_counts++;
          process_anisou_record(
            line_info, *current_input_atom_labels, *current_atom_data);
        }
      }
      else if (record_type_info.id == record_type::sigatm) {
        if (!expect_sigatm) {
          line_info.set_error(1, error_message_no_matching_atom);
        }
        else {
          expect_sigatm = false;
          sigatm_counts++;
          process_sigatm_record(
            line_info, *current_input_atom_labels, *current_atom_data);
        }
      }
      else if (record_type_info.id == record_type::siguij) {
        if (!expect_siguij) {
          line_info.set_error(1, error_message_no_matching_atom);
        }
        else {
          expect_siguij = false;
          siguij_counts++;
          process_siguij_record(
            line_info, *current_input_atom_labels, *current_atom_data);
        }
      }
      else {
        record_type_counts_[
          str6(line_info.data, line_info.size, 0, ' ')]++;
        current_input_atom_labels = 0;
        current_atom_data = 0;
        expect_anisou = false;
        expect_sigatm = false;
        expect_siguij = false;
        if      (   record_type_info.id == record_type::remark
                 || record_type_info.id == record_type::ftnote) {
          remark_section_.push_back(line_info.strip_data());
        }
        else if (record_type_info.id == record_type::ter___) {
          ter_indices_.push_back(input_atom_labels_list_.size());
          chain_tracker.transition();
        }
        else if (record_type_info.id == record_type::break_) {
          break_indices_.push_back(input_atom_labels_list_.size());
          break_record_line_numbers.push_back(line_info.line_number);
        }
        else if (record_type_info.id == record_type::model_) {
          if (model_record_oversight.model_is_allowed_here()) {
            model_ids_.push_back(read_model_id(line_info));
          }
        }
        else if (record_type_info.id == record_type::endmdl) {
          if (model_record_oversight.endmdl_is_allowed_here()) {
            model_indices_.push_back(input_atom_labels_list_.size());
            chain_tracker.endmdl();
          }
        }
        else if (record_type_info.section == record_type::title) {
          title_section_.push_back(line_info.strip_data());
          chain_tracker.transition();
        }
        else if (record_type_info.section
                 == record_type::primary_structure) {
          primary_structure_section_.push_back(line_info.strip_data());
          chain_tracker.transition();
        }
        else if (record_type_info.section
                 == record_type::heterogen) {
          heterogen_section_.push_back(line_info.strip_data());
          chain_tracker.transition();
        }
        else if (record_type_info.section
                 == record_type::secondary_structure) {
          secondary_structure_section_.push_back(line_info.strip_data());
          chain_tracker.transition();
        }
        else if (record_type_info.section
                 == record_type::connectivity_annotation) {
          connectivity_annotation_section_.push_back(
            line_info.strip_data());
          chain_tracker.transition();
        }
        else if (record_type_info.section
                 == record_type::miscellaneous_features) {
          miscellaneous_features_section_.push_back(
            line_info.strip_data());
          chain_tracker.transition();
        }
        else if (record_type_info.section
                 == record_type::crystallographic) {
          crystallographic_section_.push_back(line_info.strip_data());
          chain_tracker.transition();
        }
        else if (record_type_info.section
                 == record_type::connectivity) {
          connectivity_section_.push_back(line_info.strip_data());
          chain_tracker.transition();
        }
        else if (record_type_info.section
                 == record_type::bookkeeping) {
          bookkeeping_section_.push_back(line_info.strip_data());
          chain_tracker.transition();
        }
        else if (!line_info.is_blank_line()) {
          unknown_section_.push_back(line_info.strip_data());
        }
      }
      line_info.check_and_throw_runtime_error();
    }
    chain_indices_ = chain_tracker.finish();
    if (model_record_oversight.expecting_endmdl_record()) {
      throw std::invalid_argument("ENDMDL record missing at end of input.");
    }
    if (   model_indices_.size() == 0
        && input_atom_labels_list_.size() != 0) {
      model_ids_.push_back(str8());
      model_indices_.push_back(input_atom_labels_list_.size());
    }
    IOTBX_ASSERT(model_indices_.size() == model_ids_.size());
    IOTBX_ASSERT(model_indices_.size() == chain_indices_.size());
    if (atom___counts != 0) record_type_counts_["ATOM  "] += atom___counts;
    if (hetatm_counts != 0) record_type_counts_["HETATM"] += hetatm_counts;
    if (sigatm_counts != 0) record_type_counts_["SIGATM"] += sigatm_counts;
    if (anisou_counts != 0) record_type_counts_["ANISOU"] += anisou_counts;
    if (siguij_counts != 0) record_type_counts_["SIGUIJ"] += siguij_counts;
  }

  af::shared<std::string>
  input::model_ids() const
  {
    af::shared<std::string> result((af::reserve(model_ids_.size())));
    const str8* model_ids_end = model_ids_.end();
    for(const str8* i=model_ids_.begin();i!=model_ids_end;i++) {
      result.push_back(std::string(i->elems));
    }
    return result;
  }

  af::shared<std::size_t>
  input::model_atom_counts() const
  {
    af::shared<std::size_t> result(af::reserve(model_indices_.size()));
    range_loop<std::size_t> mr(model_indices_.const_ref());
    while (mr.next()) result.push_back(mr.size);
    return result;
  }

  af::shared<hierarchy::atom_with_labels>
  input::atoms_with_labels() const
  {
    input_atoms_with_labels_af_shared g;
    g.result.reserve(atoms_.size());
    g.run(*this);
    return g.result;
  }

}} // namespace iotbx::pdb
