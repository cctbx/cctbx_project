#include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <scitbx/array_family/accessors/flex_grid.h>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  struct flex_grid_wrappers : boost::python::pickle_suite
  {
    typedef flex_grid_default_index_type df_i_t;
    typedef flex_grid<> w_t;
    typedef flex_grid_default_index_type::value_type ivt;
    typedef boost::python::class_<w_t> c_w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      set_focus_convenience_overloads, set_focus, 1, 6)

    static boost::python::tuple
    getinitargs(w_t const& fg)
    {
      bool open_range = true;
      return boost::python::make_tuple(
        fg.origin(),
        fg.last(open_range),
        open_range);
    }

    static df_i_t
    getstate(w_t const& fg)
    {
      return fg.focus();
    }

    static void
    setstate(w_t& fg, df_i_t const& state)
    {
      fg.set_focus(state);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> copy_const_reference;
      c_w_t("grid")
        .def(init<df_i_t const&>((arg_("all"))))
        .def(init<ivt const&, optional<ivt const&, ivt const&,
                  ivt const&, ivt const&, ivt const&> >((
          arg_("all_0"), arg_("all_1"), arg_("all_2"),
          arg_("all_3"), arg_("all_4"), arg_("all_5"))))
        .def(init<df_i_t const&, df_i_t const&, bool>((
          arg_("origin"),
          arg_("last"),
          arg_("open_range")=true)))
        .def("set_focus",
          (w_t(w_t::*)(df_i_t const&, bool)) &w_t::set_focus, (
            arg_("focus"), arg_("open_range")=true))
        .def("set_focus",
          (w_t(w_t::*)(ivt const&, ivt const&, ivt const&,
                       ivt const&, ivt const&, ivt const&)) 0,
            set_focus_convenience_overloads())
        .def("nd", &w_t::nd)
        .def("size_1d", &w_t::size_1d)
        .def("is_0_based", &w_t::is_0_based)
        .def("origin", &w_t::origin)
        .def("all", &w_t::all, copy_const_reference())
        .def("last", (df_i_t(w_t::*)(bool)) &w_t::last, (
          arg_("open_range")=true))
        .def("is_padded", &w_t::is_padded)
        .def("focus", (df_i_t(w_t::*)(bool)) &w_t::focus, (
          arg_("open_range")=true))
        .def("focus_size_1d", &w_t::focus_size_1d)
        .def("is_trivial_1d", &w_t::is_trivial_1d)
        .def("shift_origin", &w_t::shift_origin)
        .def("is_valid_index", &w_t::is_valid_index, (arg_("index")))
        .def("__call__",
          (std::size_t(w_t::*)(df_i_t const&) const) &w_t::operator(), (
            arg_("index")))
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def_pickle(flex_grid_wrappers())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_flex_grid()
  {
    flex_grid_wrappers::wrap();
  }

}}} // namespace scitbx::af::boost_python
