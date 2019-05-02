#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

#include <mmtbx/tls/utils.h>

// TLSMatrices Overloads
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mat_get_overloads, getValuesByString, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mat_set_overloads, setValuesByString, 1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mat_any_overloads, any, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mat_dec_overloads, decompose, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mat_val_overloads, isValid, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mat_nrm_overloads, normalise, 2, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mat_prm_overloads, paramCount, 0, 2)
// TLSAmplitudes Overloads
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(amp_any_overloads, any, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(amp_nrm_overloads, normalise, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(amp_prm_overloads, paramCount, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(amp_set_overloads, setValues, 1, 2)
// TLSMatricesAndAmplitudes Overloads
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(maa_nul_overloads, isNull, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(maa_iv1_overloads, isValid, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(maa_iv2_overloads, isValid, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(maa_nrm_overloads_a, normaliseByAmplitudes, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(maa_nrm_overloads_m, normaliseByMatrices, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(maa_prm_overloads, paramCount, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(maa_res_overloads, resetIfNull, 0, 2)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(maa_uij_overloads, uijs, 2, 3)
// TLSMatricesAndAmplitudesList Overloads
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mal_nul_overloads, isNull, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mal_nrm_overloads_a, normaliseByAmplitudes, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mal_nrm_overloads_m, normaliseByMatrices, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mal_prm_overloads, paramCount, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mal_res_overloads, resetNullModes, 0, 2)

namespace mmtbx { namespace tls { namespace utils {
  namespace bp = boost::python;
  namespace af = scitbx::af;

namespace {

  // Define Pickle Functions for classes
  struct tlsmal_pickle_suite : bp::pickle_suite
  {
    static
    boost::python::tuple
    getinitargs(const TLSMatricesAndAmplitudesList& l)
    {
      // Three arrays
      // 1) matrices (length * 21)
      // 2) amplitudes (length * n_amplitudes)

      // Extract the size of the list, and the number of amplitudes in each
      size_t n_set = l.size();
      size_t n_amp = l.getConst(0)->getAmplitudesConst()->size();
      // Create arrays of values
      dblArrNd mat_values(af::flex_grid<>(n_set, 21));
      dblArrNd amp_values(af::flex_grid<>(n_set, n_amp)); // do not need to initialise to zeros as all elements filled
      // Iterate through list
      for (int i = 0; i < n_set; i++) {
        TLSMatricesAndAmplitudes* ma = l.getConst(i);
        // Extract matrices and amplitudes
        dblArr1d ma_mats = ma->getMatricesConst()->getValuesByString("TLS", true);
        dblArr1d ma_amps = ma->getAmplitudesConst()->getValues();
        // Copy to output arrays
        memcpy(&mat_values(i,0), &ma_mats[0], sizeof(double) * 21);
        memcpy(&amp_values(i,0), &ma_amps[0], sizeof(double) * n_amp);
      }
      return bp::make_tuple( mat_values, amp_values );
    }
  };

  bp::list expand_wrapper(TLSMatricesAndAmplitudes& c)
  {
    bp::list result;
    // call original function
    af::shared<TLSMatrices> v = c.expand();
    // put all the strings inside the python list
    af::shared<TLSMatrices>::iterator it;
    for (it = v.begin(); it != v.end(); ++it){
      result.append(*it);
    }
    return result;
  }

  void init_module()
 {
    using namespace boost::python;
    using boost::python::arg;

    def("uij_eigenvalues", uij_eigenvalues);

    class_<TLSMatrices>(
        "TLSMatrices",
        init< scitbx::sym_mat3<double> const&,
              scitbx::sym_mat3<double> const&,
              scitbx::mat3<double> const& >
              ( ( arg("T"), arg("L"), arg("S") ) )
        )
      // Other constructors
      .def( init<> ()
      )
      .def( init< TLSMatrices const& >
          ( ( arg("other") ) )
      )
      .def( init< af::shared<double> >
          ( ( arg("values") ) )
      )

      // Properties
      .add_property("T", &TLSMatrices::getT, "Return the values of the T-matrix")
      .add_property("L", &TLSMatrices::getL, "Return the values of the L-matrix")
      .add_property("S", &TLSMatrices::getS, "Return the values of the S-matrix")

      // Class/static methods
      .def("get_tolerance", &TLSMatrices::getTolerance,
          "Return the current tolerance for which a number is considered to be non-zero."
          ).staticmethod("get_tolerance")
      .def("set_tolerance", &TLSMatrices::setTolerance,
          ( arg("tolerance") ),
          "Set the tolerance of the object for determining when values are non-zero."
          ).staticmethod("set_tolerance")
      .def("get_precision", &TLSMatrices::getPrecision,
          "Return the number of decimal places used for rounding."
          ).staticmethod("get_precision")
      .def("set_precision", &TLSMatrices::setPrecision,
          ( arg("precision") ),
          "Set the number of decimal places for rounding matrix values."
          ).staticmethod("set_precision")

      // Operator overloads
      .def(self + self)       // add operator
      .def(self * double())   // l-multiply operator
      .def(double() * self)   // r-multiply operator

      // Instance Methods
      .def("add", &TLSMatrices::add,
          ( arg("other") ),
          "Add the matrix values from another TLSMatrix object."
          )
      .def("any", &TLSMatrices::any,
          mat_any_overloads(
            ( arg("component_string"), arg("tolerance") ),
            "Returns True if any of the selected matrix values are non-zero, else False.\n"
            "Takes <component_string> as a string containing letters T-L-S or a combination, and a tolerance. \n"
            "If <tolerance> is not defined, the class tolerance is used."
            )
          )
      .def("copy", &TLSMatrices::copy,
          "Create a completely separate copy of the object.",
          return_value_policy<manage_new_object>()
          )
      .def("decompose", &TLSMatrices::decompose,
          mat_dec_overloads(
            ( arg("tolerance") ),
            "Perform TLS Decomposition into fundamental motions.\n"
            "Returns a TLSDecomposition object.\n"
            "If <tolerance> is not defined, the class tolerance is used."
            )
          )
      .def("get", &TLSMatrices::getValuesByString,
          mat_get_overloads(
            ( arg("component_string"), arg("include_szz") ),
            "Get values of selected matrices as a single array. <component_string> must be a string containing letters T, L, S or a combination. Letters must be in the order T-L-S."
            )
          ) // Get internal matrix values
      .def("is_valid", &TLSMatrices::isValid,
          mat_val_overloads(
            ( arg("tolerance") ),
            "Returns True if matrices can be decomposed into fundamental motions, and False otherwise.\n"
            "If <tolerance> is not defined, the class tolerance is used."
            )
          )
      .def("multiply", &TLSMatrices::multiply,
          ( arg("scalar") ),
          "Multiply all matrix values by a constant."
          )
      .def("n_params", &TLSMatrices::paramCount,
          mat_prm_overloads(
            ( args("free"), arg("non_zero") ),
            "Return the number of parameters in the object.\n"
            "free=True/False controls whether Szz is included in the calculation.\n"
            "Setting non_zero=True returns the number of parameters that are non-zero."
            )
          )
      .def("normalise", &TLSMatrices::normalise,
          mat_nrm_overloads(
            ( arg("sites_cart"), arg("origin"), arg("target"), arg("tolerance") ),
            "Scale the matrix values so that the average eigenvalue of the generated Uij (from <sites_cart> and <origin>) is <target>.\n"
            "Scaling will not be performed if the average eigenvalue before scaling is less than <tolerance>.\n"
            "Returns the multiplier value required to apply the inverse scaling to the TLSAmplitudes"
            "If <tolerance> is not defined, the class tolerance is used."
            )
          )
      .def("reset", &TLSMatrices::reset,
          "Set all matrix values to zero."
          )
      .def("set", &TLSMatrices::setValuesByString,
          mat_set_overloads(
            ( arg("values"), arg("component_string") = "TLS", arg("include_szz") = true ),
            "Set values of selected matrices from a single array. <component_string> must be a string containing letters T, L, S or a combination. Letters must be in the order T-L-S."
            )
          ) // Set internal matrix values
      .def("summary", &TLSMatrices::summary,
          "Return a summary string of the object."
          )
      .def("uijs", &TLSMatrices::uijs,
          ( arg("sites_cart"), arg("origin" ) ),
          "Return Uijs calculated from coordinates and an origin"
          )
    ;

    class_<TLSAmplitudes>(
        "TLSAmplitudes",
        init< size_t >
            ( ( arg("n") ) )
        )
      // Other constructors
      .def( init< TLSAmplitudes const& >
          ( ( arg("other") ) )
      )
      .def( init< af::shared<double> const& >
          ( ( arg("values") ) )
      )

      // Properties
      .add_property("description", &TLSAmplitudes::getDescription,
          "A short description of the intention behind the class"
          )
      .add_property("values", &TLSAmplitudes::getValues,
          "Get the current amplitude values"
          )

      // Class/static methods
      .def("get_tolerance", &TLSAmplitudes::getTolerance,
          "Return the current tolerance for which a number is considered to be non-zero."
          ).staticmethod("get_tolerance")
      .def("set_tolerance", &TLSAmplitudes::setTolerance,
          ( arg("tolerance") ),
          "Set the tolerance of the object for determining when values are non-zero."
          ).staticmethod("set_tolerance")
      .def("get_precision", &TLSAmplitudes::getPrecision,
          "Return the number of decimal places used for rounding."
          ).staticmethod("get_precision")
      .def("set_precision", &TLSAmplitudes::setPrecision,
          ( arg("precision") ),
          "Set the number of decimal places for rounding matrix values."
          ).staticmethod("set_precision")

      // Operator overloads
      .def(self + self)       // add operator
      .def(self * double())   // l-multiply operator
      .def(double() * self)   // r-multiply operator

      // Operators/special methods
      .def("__getitem__", &TLSAmplitudes::operator[],
          ( arg( "index" ) )
          )

      // Instance methods
      .def("add", &TLSAmplitudes::add,
          ( arg("other") ),
          "Add the matrix values from another TLSAmplitudes object of the same size."
          )
      .def("any", &TLSAmplitudes::any,
          amp_any_overloads(
            ( arg("tolerance") ),
            "Returns True if any of the amplitude values are non-zero, else False.\n"
            "If <tolerance> is not defined, the class tolerance is used."
            )
          )
      .def("copy", &TLSAmplitudes::copy,
          "Create a completely separate copy of the object.",
          return_value_policy<manage_new_object>()
          )
      .def("get", &TLSAmplitudes::getValues,
          "Get all amplitude values."
          )
      .def("get", &TLSAmplitudes::getValuesBySelection,
          ( arg("selection") ),
          "Get selected amplitude values. <selection> is a list of indices for which values to return."
          )
      .def("multiply", &TLSAmplitudes::multiply,
          ( arg("scalar") ),
          "Multiply all amplitude values by a constant."
          )
      .def("normalise", &TLSAmplitudes::normalise,
          amp_nrm_overloads(
            ( arg("target") ),
            "Scale the amplitude values so that the average amplitudes is <target>\n"
            "Returns the multiplier value required to apply the inverse scaling to the TLSMatrices"
            )
          )
      .def("n_params", &TLSAmplitudes::paramCount,
          amp_prm_overloads(
            ( arg("non_zero") ),
            "Return the number of parameters in the object.\n"
            "Setting non_zero=True returns the number of parameters that are non-zero."
            )
          )
      .def("reset", &TLSAmplitudes::reset,
          "Set all amplitude values to one."
          )
      .def("set", &TLSAmplitudes::setValues,
          ( arg("values") ),
          "Set all amplitude values from an array <values>."
          ) // Set internal matrix values
      .def("set", &TLSAmplitudes::setValuesBySelection,
          ( arg ("values"), arg("selection") ),
          "Set selected amplitude values from an array <values>. <selection> is a list of indices for which values to set."
          ) // Set internal matrix values
      .def("size", &TLSAmplitudes::size,
          "Get the number of amplitudes"
          )
      .def("summary", &TLSAmplitudes::summary,
          "Get summary string for the amplitudes in the class."
          )
      .def("zero_values", &TLSAmplitudes::zeroValues,
          "Set all values to zero."
          )
      .def("zero_negative_values", &TLSAmplitudes::zeroNegativeValues,
          "Set any negative values to zero."
          )
    ;

    bool (TLSMatricesAndAmplitudes::*MA_is_valid1)
      (double tolerance)
      = &TLSMatricesAndAmplitudes::isValid;
    bool (TLSMatricesAndAmplitudes::*MA_is_valid2)
      (const selArr1d &selection, double tolerance)
      = &TLSMatricesAndAmplitudes::isValid;

    symArrNd (TLSMatricesAndAmplitudes::*MA_uijs1)
      (const vecArrNd &sites_carts, const vecArr1d &origins)
      = &TLSMatricesAndAmplitudes::uijs;
    symArrNd (TLSMatricesAndAmplitudes::*MA_uijs2)
      (const vecArrNd &sites_carts, const vecArr1d &origins, const selArr1d &selection)
      = &TLSMatricesAndAmplitudes::uijs;

    class_<TLSMatricesAndAmplitudes>(
        "TLSMatricesAndAmplitudes",
        init< size_t >
            ( ( arg("n") ) )
        )
      // Other Constructors
      .def( init< TLSMatricesAndAmplitudes const& >
          ( ( arg("other") ) )
      )
      .def( init< TLSMatrices &, TLSAmplitudes & >
          ( ( arg("matrices"), arg("amplitudes") ) )
      )
      .def( init< dblArr1d const&, dblArr1d const& >
          ( ( arg("matrix_values"), arg("amplitude_values") ) )
      )

      // Properties
      .add_property("amplitudes", bp::make_function( &TLSMatricesAndAmplitudes::getAmplitudes, return_internal_reference<>() ) )
      .add_property("matrices",   bp::make_function( &TLSMatricesAndAmplitudes::getMatrices,   return_internal_reference<>() ) )

      // Instance methods
      .def("copy", &TLSMatricesAndAmplitudes::copy,
          "Create a completely separate copy of the object.",
          return_value_policy<manage_new_object>()
          )
      .def("expand", expand_wrapper,
          "Multiply the TLSMatrices object by each amplitude in the TLSAmplitudes object.\n"
          "Returns a list of the multiplied objects."
          )
      .def("is_null", &TLSMatricesAndAmplitudes::isNull,
          maa_nul_overloads(
            ( arg("matrices_tolerance"), arg("amplitudes_tolerance") ),
            "Returns True if all matrix values are 0.0 OR all amplitudes are 0.0 (i.e. whether this will return non-zero Uijs)"
            )
          )
      .def("is_valid", MA_is_valid1,
          maa_iv1_overloads(
            ( arg("tolerance") ),
            "Returns True if (a selection of) the amplitude-multiplied matrices can be decomposed into fundamental motions, and False otherwise.\n"
            "If <tolerance> is not defined, the class tolerance is used."
            )
          )
      .def("is_valid", MA_is_valid2,
          maa_iv2_overloads(
            ( arg("selection"), arg("tolerance") ),
            "Returns True if (a selection of) the amplitude-multiplied matrices can be decomposed into fundamental motions, and False otherwise.\n"
            "If <tolerance> is not defined, the class tolerance is used."
            )
          )
      .def("n_params", &TLSMatricesAndAmplitudes::paramCount,
          maa_prm_overloads(
            ( args("free"), arg("non_zero") ),
            "Return the number of parameters in the object.\n"
            "free=True/False controls whether Szz is included in the calculation.\n"
            "Setting non_zero=True returns the number of parameters that are non-zero."
            )
          )
      .def("normalise_by_amplitudes", &TLSMatricesAndAmplitudes::normaliseByAmplitudes,
          maa_nrm_overloads_a(
            ( arg("target") ),
            "Normalise the TLSAmplitudes to the values of 'target'. Apply the inverse scaling to the TLSMatrices."
            )
          )
      .def("normalise_by_matrices", &TLSMatricesAndAmplitudes::normaliseByMatrices,
          maa_nrm_overloads_m(
            ( arg("sites_carts"), arg("origins"), arg("target") ),
            "Normalise the TLSMatrices to the values of 'target'. Apply the inverse scaling to the Amplitudes."
            )
          )
      .def("reset", &TLSMatricesAndAmplitudes::reset,
          "Reset all matrix values to 0.0 and all amplitudes to 1.0"
          )
      //.def("set_label", &TLSMatricesAndAmplitudes::setLabel,
      //    ( arg("label") ),
      //    "Label the object with an integer."
      //    )
      .def("summary", &TLSMatricesAndAmplitudes::summary,
          "Return a summary string for the object"
          )
      .def("uijs", MA_uijs1,
          ( arg("sites_carts"), arg("origins") ),
          "Calculate Uijs for all of the contained amplitudes"
          )
      .def("uijs", MA_uijs2,
          ( arg("sites_carts"), arg("origins"), arg("selection") ),
          "Calculate Uijs for a selection of the contained amplitudes"
          )
      .def("reset_if_null", &TLSMatricesAndAmplitudes::resetIfNull,
          maa_res_overloads(
            ( arg("matrices_tolerance"), arg("amplitudes_tolerance") ),
            "Reset matrices and amplitudes if either is all zeros"
            )
          )
    ;

    symArrNd (TLSMatricesAndAmplitudesList::*MAL_uijs1)
      (const vecArrNd &sites_carts, const vecArr1d &origins)
      = &TLSMatricesAndAmplitudesList::uijs;
    symArrNd (TLSMatricesAndAmplitudesList::*MAL_uijs2)
      (const vecArrNd &sites_carts, const vecArr1d &origins, const selArr1d &selection)
      = &TLSMatricesAndAmplitudesList::uijs;

    class_<TLSMatricesAndAmplitudesList>(
        "TLSMatricesAndAmplitudesList",
        init< size_t, size_t >
            ( ( arg("length"), arg("n_amplitudes") ) )
            )

      // Constructor from arrays
      .def( init< dblArrNd const&, dblArrNd const& >
          ( ( arg("matrix_values"), arg("amplitude_values") ) )
      )

      // Define pickle
      .def_pickle(tlsmal_pickle_suite())

      // Operators/special methods
      .def("__getitem__", &TLSMatricesAndAmplitudesList::operator[],
          ( arg( "index" ) ),
          return_internal_reference<>()
          )

      // Instance methods
      .def("copy", &TLSMatricesAndAmplitudesList::copy,
          "Create a completely separate copy of the object.",
          return_value_policy<manage_new_object>()
          )
      .def("get", &TLSMatricesAndAmplitudesList::get,
          ( arg("index") ),
          "Get a TLSMatricesAndAmplitudes object by index",
          return_internal_reference<>()
          )
      .def("is_null", &TLSMatricesAndAmplitudesList::isNull,
          mal_nul_overloads(
            ( arg("matrices_tolerance"), arg("amplitudes_tolerance") ),
            "Returns True if all matrix values are 0.0 OR all amplitudes are 0.0 (i.e. whether this will return non-zero Uijs) for ALL modes."
            )
          )
      .def("n_params", &TLSMatricesAndAmplitudesList::paramCount,
          mal_prm_overloads(
            ( args("free"), arg("non_zero") ),
            "Return the number of parameters in the object.\n"
            "free=True/False controls whether Szz is included in the calculation.\n"
            "Setting non_zero=True returns the number of parameters that are non-zero."
            )
          )
      .def("normalise_by_amplitudes", &TLSMatricesAndAmplitudesList::normaliseByAmplitudes,
          mal_nrm_overloads_a(
            ( arg("target") ),
            "Normalise all TLSAmplitudes to the values of 'target'. Apply the inverse scaling to the TLSMatrices."
            )
          )
      .def("normalise_by_matrices", &TLSMatricesAndAmplitudesList::normaliseByMatrices,
          mal_nrm_overloads_m(
            ( arg("sites_carts"), arg("origins"), arg("target") ),
            "Normalise all TLSMatrices to the values of 'target'. Apply the inverse scaling to the Amplitudes."
            )
          )
      .def("reset", &TLSMatricesAndAmplitudesList::reset,
          "Reset all matrices to zeros and amplitudes to ones"
          )
      .def("reset_matrices", &TLSMatricesAndAmplitudesList::resetMatrices,
          "Reset all matrices to zeros"
          )
      .def("reset_null_modes", &TLSMatricesAndAmplitudesList::resetNullModes,
          mal_res_overloads(
            ( arg("matrices_tolerance"), arg("amplitudes_tolerance") ),
            "For each mode, reset amplitudes and matrices if either is approximately zero"
            )
          )
      .def("size", &TLSMatricesAndAmplitudesList::size,
          "Get the number of amplitudes"
          )
      .def("summary", &TLSMatricesAndAmplitudesList::summary,
          "Get summary string for the class."
          )
      .def("uijs", MAL_uijs1,
          ( args("sites_carts"), arg("origins") ),
          "Calculate total Uijs for all of the contained modes"
          )
      .def("uijs", MAL_uijs2,
          ( args("sites_carts"), arg("origins"), arg("selection") ),
          "Calculate total Uijs for all of the contained modes (for a selection of amplitudes)"
          )
      .def("zero_amplitudes", &TLSMatricesAndAmplitudesList::zeroAmplitudes,
          ( arg("selection") ),
          "Set amplitudes to zero for the selected modes"
          )
      .def("zero_negative_amplitudes", &TLSMatricesAndAmplitudesList::zeroNegativeAmplitudes,
          "Set negative amplitudes to zero for all containted classes"
          )
    ;
  }

} // namespace <anonymous>
}}} // namespace mmtbx::tls::utils

BOOST_PYTHON_MODULE(mmtbx_tls_utils_ext)
{
  mmtbx::tls::utils::init_module();
}
