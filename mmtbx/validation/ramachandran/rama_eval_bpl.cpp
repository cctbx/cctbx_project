#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <mmtbx/validation/ramachandran/rama_eval.h>


namespace mmtbx { namespace validation { namespace ramachandran {
namespace {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(eval_score_int, rama_eval::evaluate_score, 2, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(eval_score_str, rama_eval::evaluate_score, 2, 2)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_val_int, rama_eval::get_value, 3, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_val_str, rama_eval::get_value, 3, 3)

  void init_module()
  {
    using namespace boost::python;


    typedef return_value_policy<return_by_value> rbv;

    class_<rama_eval>("rama_eval", no_init)
      .def(init<>())

      .def("evaluate_score", static_cast< int(rama_eval::*)
          (std::string const&, double const&)>
          (&rama_eval::evaluate_score), eval_score_int())

      .def("evaluate_score", static_cast< int(rama_eval::*)
          (int const&, double const&)>
          (&rama_eval::evaluate_score), eval_score_str())

      .def("get_value", static_cast< double(rama_eval::*)
          (std::string const&, double const&, double const&)>
          (&rama_eval::get_value), get_val_str())

      .def("get_value", static_cast< double(rama_eval::*)
          (int const&, double const&, double const&)>
          (&rama_eval::get_value), get_val_int())

    ;
  }

} // namespace <anonymous>
}}} // namespace mmtbx::validation::ramachandran

BOOST_PYTHON_MODULE(mmtbx_validation_ramachandran_ext)
{
  mmtbx::validation::ramachandran::init_module();
}
