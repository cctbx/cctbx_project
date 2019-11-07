#ifndef MMTBX_TLS_OPTIMISE_AMPLITUDES_H
#define MMTBX_TLS_OPTIMISE_AMPLITUDES_H

#include <string>

// Basic data types
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/sym_mat3.h>

// Allow arrays of the above
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>

#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>

namespace mmtbx { namespace tls { namespace optimise {

// These are used a lot...
namespace af = scitbx::af;
namespace bp = boost::python;

// Commonly used types
typedef scitbx::vec3<double> vec;
typedef scitbx::mat3<double> mat;
typedef scitbx::sym_mat3<double> sym;
// 1-dimensional arrays
typedef af::shared<double> dblArr1d;
typedef af::shared<size_t> selArr1d;
typedef af::shared<int>    intArr1d;
typedef af::shared<bool>   blnArr1d;
typedef af::shared<sym>    symArr1d;
// N-dimensional arrays
typedef af::versa<bool,   af::flex_grid<> > blnArrNd;
typedef af::versa<double, af::flex_grid<> > dblArrNd;
typedef af::versa<sym,    af::flex_grid<> > symArrNd;

class MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator {

  public:
    //! Main constructor
    MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
        const symArrNd &target_uijs,
        const dblArrNd &target_weights,
        const dblArr1d &base_amplitudes,
        const bp::list &base_uijs,
        const bp::list &base_atom_indices,
        const selArr1d &dataset_hash,
        const symArr1d &residual_uijs
        );

    dblArr1d getCurrentAmplitudes();
    void setCurrentAmplitudes(const dblArr1d &values);
    void printCurrentAmplitudes();

    void setResidualMask(const blnArr1d &mask);

    bp::tuple computeFunctionalAndGradients();

  private:
    // Input variables
    //
    const symArrNd target_uijs;
    const dblArrNd target_weights;
    af::shared<symArr1d*> base_u;  // Base Uijs
    af::shared<selArr1d*> base_i;  // Atom indices for base uijs
    const selArr1d dataset_hash; // maps base elements to datasets
    const symArr1d residual_uijs;
    blnArr1d residual_mask;

    // Quantities calculated from input variables
    //
    const size_t n_dst, n_atm, n_base, n_total;
    int residual_mask_total; // Number of datasets to use for residual optimisation

    // Internal intermediate variables
    //
    dblArr1d initial_amplitudes;
    dblArr1d current_amplitudes;
    symArrNd total_uijs;

    // Optimisation/output variables
    //
    double functional;
    dblArr1d gradients;

    size_t n_call;

    void zero();
    void sanitise_current_amplitudes();
    void calculate_total_uijs();
    void calculate_f_g_least_squares();

};

} } } // close namepsace mmtbx/tls/optimise

#endif
