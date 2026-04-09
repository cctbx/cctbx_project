// cctbx_project/xcif/xcif_ext.cpp
// Boost.Python bindings for the xcif CIF parser.
// IP4: aliasing shared_ptr, GIL release for parse, schema-aware interning.

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <xcif/data_model.h>
#include <xcif/numeric.h>
#include <xcif/mapped_file.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <memory>
#include <cmath>
#include <stdexcept>

namespace bp = boost::python;
using namespace xcif;

// ===== GIL RAII guard =============================================

struct GILRelease {
  PyThreadState* save;
  GILRelease() : save(PyEval_SaveThread()) {}
  ~GILRelease() { PyEval_RestoreThread(save); }
};

// ===== String helpers =============================================

static bp::object intern_sv(const string_view& sv) {
  PyObject* s = PyUnicode_FromStringAndSize(sv.data(), sv.size());
  if (!s) bp::throw_error_already_set();
  PyUnicode_InternInPlace(&s);
  return bp::object(bp::handle<>(s));
}

static bp::object plain_sv(const string_view& sv) {
  PyObject* s = PyUnicode_FromStringAndSize(sv.data(), sv.size());
  if (!s) bp::throw_error_already_set();
  return bp::object(bp::handle<>(s));
}

// Schema-aware categorical tags (IP4 decision 4C)
static const char* const CATEGORICAL_TAGS[] = {
  "_atom_site.label_comp_id",
  "_atom_site.type_symbol",
  "_atom_site.label_atom_id",
  "_atom_site.label_asym_id",
  "_atom_site.label_entity_id",
  "_atom_site.label_alt_id",
  "_atom_site.auth_comp_id",
  "_atom_site.auth_atom_id",
  "_atom_site.auth_asym_id",
};
static const std::size_t N_CAT =
  sizeof(CATEGORICAL_TAGS) / sizeof(CATEGORICAL_TAGS[0]);

static bool is_categorical(const string_view& tag) {
  CIEqual eq;
  for (std::size_t i = 0; i < N_CAT; ++i)
    if (eq(tag, string_view(CATEGORICAL_TAGS[i]))) return true;
  return false;
}

// ===== Wrapper types (shared_ptr keeps Document alive) ============

typedef std::shared_ptr<Document> DocPtr;

struct BlockWrap {
  DocPtr doc;
  const Block* p;
  BlockWrap() : p(0) {}
  BlockWrap(DocPtr d, const Block* b) : doc(d), p(b) {}
};

struct LoopWrap {
  DocPtr doc;
  const Loop* p;
  LoopWrap() : p(0) {}
  LoopWrap(DocPtr d, const Loop* l) : doc(d), p(l) {}
};

// ===== Parse (GIL released) =======================================

static DocPtr py_parse(const std::string& text) {
  DocPtr result;
  {
    GILRelease gil;
    result = std::make_shared<Document>(
      parse(text.c_str(), text.size(), "<string>"));
  }
  return result;
}

static DocPtr py_parse_file(const std::string& path) {
  DocPtr result;
  {
    GILRelease gil;
    result = std::make_shared<Document>(parse_file(path.c_str()));
  }
  return result;
}

// ===== Document wrappers ==========================================

static std::size_t doc_len(DocPtr self) { return self->size(); }

static BlockWrap doc_getitem(DocPtr self, int idx) {
  if (idx < 0 || static_cast<std::size_t>(idx) >= self->size()) {
    PyErr_SetString(PyExc_IndexError, "block index out of range");
    bp::throw_error_already_set();
  }
  return BlockWrap(self, &(*self)[static_cast<std::size_t>(idx)]);
}

static bp::object doc_find_block(DocPtr self, const std::string& name) {
  const Block* b = self->find_block(string_view(name));
  if (!b) return bp::object();
  return bp::object(BlockWrap(self, b));
}

// ===== Block wrappers =============================================

static bp::object block_name(const BlockWrap& self) {
  return intern_sv(self.p->name());
}

static bool block_has_tag(const BlockWrap& self, const std::string& t) {
  return self.p->has_tag(string_view(t));
}

static bp::object block_find_value(const BlockWrap& self,
                                   const std::string& tag) {
  string_view sv(tag);
  if (!self.p->has_tag(sv)) return bp::object();
  return intern_sv(self.p->find_value(sv));
}

static bp::object block_find_loop(const BlockWrap& self,
                                  const std::string& tag) {
  const Loop* lp = self.p->find_loop(string_view(tag));
  if (!lp) return bp::object();
  return bp::object(LoopWrap(self.doc, lp));
}

static bp::list block_loops(const BlockWrap& self) {
  bp::list result;
  const std::vector<Loop>& v = self.p->loops();
  for (std::size_t i = 0; i < v.size(); ++i)
    result.append(LoopWrap(self.doc, &v[i]));
  return result;
}

static bp::list block_pair_tags(const BlockWrap& self) {
  bp::list result;
  const std::vector<std::pair<string_view, string_view>>& p = self.p->pairs();
  for (std::size_t i = 0; i < p.size(); ++i)
    result.append(intern_sv(p[i].first));
  return result;
}

static bp::list block_pair_values(const BlockWrap& self) {
  bp::list result;
  const std::vector<std::pair<string_view, string_view>>& p = self.p->pairs();
  for (std::size_t i = 0; i < p.size(); ++i)
    result.append(plain_sv(p[i].second));
  return result;
}

static bp::list block_save_frames(const BlockWrap& self) {
  bp::list result;
  const std::vector<Block>& sfs = self.p->save_frames();
  for (std::size_t i = 0; i < sfs.size(); ++i)
    result.append(BlockWrap(self.doc, &sfs[i]));
  return result;
}

static bp::object block_find_save_frame(const BlockWrap& self,
                                        const std::string& name) {
  const Block* sf = self.p->find_save_frame(string_view(name));
  if (!sf) return bp::object();
  return bp::object(BlockWrap(self.doc, sf));
}

// ===== Loop wrappers ==============================================

static std::size_t loop_width(const LoopWrap& self) {
  return self.p->width();
}
static std::size_t loop_length(const LoopWrap& self) {
  return self.p->length();
}

static bp::list loop_tags(const LoopWrap& self) {
  bp::list result;
  const std::vector<string_view>& t = self.p->tags();
  for (std::size_t i = 0; i < t.size(); ++i)
    result.append(intern_sv(t[i]));
  return result;
}

static bool loop_has_tag(const LoopWrap& self, const std::string& t) {
  return self.p->has_tag(string_view(t));
}

static int loop_column_index(const LoopWrap& self,
                             const std::string& tag) {
  std::size_t ci = self.p->column_index(string_view(tag));
  return ci == std::size_t(-1) ? -1 : static_cast<int>(ci);
}

static bp::object loop_value(const LoopWrap& self,
                             std::size_t row, std::size_t col) {
  if (row >= self.p->length() || col >= self.p->width()) {
    PyErr_SetString(PyExc_IndexError, "loop cell out of range");
    bp::throw_error_already_set();
  }
  return plain_sv(self.p->value(row, col));
}

static bp::list loop_column(const LoopWrap& self,
                            const std::string& tag) {
  string_view sv_tag(tag);
  std::vector<string_view> col = self.p->column(sv_tag);
  bool do_intern = is_categorical(sv_tag);
  bp::list result;
  for (std::size_t i = 0; i < col.size(); ++i)
    result.append(do_intern ? intern_sv(col[i]) : plain_sv(col[i]));
  return result;
}

// ===== Flex column extraction =====================================

static scitbx::af::shared<double>
loop_column_flex_double(const LoopWrap& self, const std::string& tag) {
  string_view sv_tag(tag);
  std::size_t ci = self.p->column_index(sv_tag);
  if (ci == std::size_t(-1)) {
    PyErr_SetString(PyExc_KeyError, tag.c_str());
    bp::throw_error_already_set();
  }
  std::size_t n = self.p->length();
  scitbx::af::shared<double> result(n);
  for (std::size_t r = 0; r < n; ++r)
    result[r] = as_double(self.p->value(r, ci));
  return result;
}

static scitbx::af::shared<int>
loop_column_flex_int(const LoopWrap& self, const std::string& tag) {
  string_view sv_tag(tag);
  std::size_t ci = self.p->column_index(sv_tag);
  if (ci == std::size_t(-1)) {
    PyErr_SetString(PyExc_KeyError, tag.c_str());
    bp::throw_error_already_set();
  }
  std::size_t n = self.p->length();
  scitbx::af::shared<int> result(n);
  for (std::size_t r = 0; r < n; ++r)
    result[r] = as_int(self.p->value(r, ci));
  return result;
}

static scitbx::af::shared<std::string>
loop_column_flex_string(const LoopWrap& self, const std::string& tag) {
  string_view sv_tag(tag);
  std::size_t ci = self.p->column_index(sv_tag);
  if (ci == std::size_t(-1)) {
    PyErr_SetString(PyExc_KeyError, tag.c_str());
    bp::throw_error_already_set();
  }
  std::size_t n = self.p->length();
  scitbx::af::shared<std::string> result(n);
  for (std::size_t r = 0; r < n; ++r) {
    string_view v = self.p->value(r, ci);
    result[r] = std::string(v.data(), v.size());
  }
  return result;
}

// ===== Numeric wrappers ===========================================

static double py_as_double(const std::string& s) {
  return as_double(string_view(s));
}

static int py_as_int(const std::string& s) {
  return as_int(string_view(s));
}

static bp::tuple py_as_double_with_su(const std::string& s) {
  std::pair<double, double> r = as_double_with_su(string_view(s));
  return bp::make_tuple(r.first, r.second);
}

static bool py_is_null(const std::string& s) {
  return is_null(string_view(s));
}

static bool py_is_unknown(const std::string& s) {
  return is_unknown(string_view(s));
}

static bool py_is_inapplicable(const std::string& s) {
  return is_inapplicable(string_view(s));
}

// ===== Module definition ==========================================

static void translate_invalid_argument(const std::invalid_argument& e) {
  PyErr_SetString(PyExc_ValueError, e.what());
}

static void translate_overflow_error(const std::overflow_error& e) {
  PyErr_SetString(PyExc_OverflowError, e.what());
}

BOOST_PYTHON_MODULE(xcif_ext)
{
  bp::register_exception_translator<std::invalid_argument>(
    translate_invalid_argument);
  bp::register_exception_translator<std::overflow_error>(
    translate_overflow_error);
  bp::class_<BlockWrap>("Block", bp::no_init)
    .add_property("name", &block_name)
    .def("has_tag", &block_has_tag)
    .def("find_value", &block_find_value)
    .def("find_loop", &block_find_loop)
    .add_property("loops", &block_loops)
    .add_property("pair_tags", &block_pair_tags)
    .add_property("pair_values", &block_pair_values)
    .add_property("save_frames", &block_save_frames)
    .def("find_save_frame", &block_find_save_frame)
    ;

  bp::class_<LoopWrap>("Loop", bp::no_init)
    .add_property("width", &loop_width)
    .add_property("length", &loop_length)
    .add_property("tags", &loop_tags)
    .def("has_tag", &loop_has_tag)
    .def("column_index", &loop_column_index)
    .def("value", &loop_value)
    .def("column", &loop_column)
    .def("column_as_flex_double", &loop_column_flex_double)
    .def("column_as_flex_int", &loop_column_flex_int)
    .def("column_as_flex_string", &loop_column_flex_string)
    ;

  bp::class_<Document, DocPtr, boost::noncopyable>("Document", bp::no_init)
    .def("__len__", &doc_len)
    .def("__getitem__", &doc_getitem)
    .def("find_block", &doc_find_block)
    ;

  // Free functions
  bp::def("parse", &py_parse);
  bp::def("parse_file", &py_parse_file);
  bp::def("as_double", &py_as_double);
  bp::def("as_int", &py_as_int);
  bp::def("as_double_with_su", &py_as_double_with_su);
  bp::def("is_null", &py_is_null);
  bp::def("is_unknown", &py_is_unknown);
  bp::def("is_inapplicable", &py_is_inapplicable);
}
