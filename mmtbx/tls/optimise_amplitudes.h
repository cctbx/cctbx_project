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

//!  Object for calculating functional and gradients for fitting a sets of Uijs to a target (multi-dimensional) list of Uijs.
/*!
  Uses a least-squares target function with weights for each atom (in each dataset),
  and a combination of optional penalty terms on the amplitudes of each of the individual components.
*/
class MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator {

  public:
    /*!
      Main contructor
    */
    MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
      const symArrNd &target_uijs,                //!< Array of target symmetric Uij matrices for each dataset (shape: n_datasets, n_atoms).
      const dblArrNd &target_weights,             //!< Array of target weights during optimisation (shape: n_datasets, n_atoms).
      const dblArr1d &base_amplitudes,            //!< Array of starting amplitudes for each of the base Uij terms (length n_base).
      const bp::list &base_uijs,                  //!< Sets of symmetric Uij matrices that are the "components" for which the amplitudes are being optimised (length n_base).
      const bp::list &base_atom_indices,          //!< Array of indices which indicate which atom each of the components in base_uijs corresponds to (length n_base).
      const selArr1d &base_dataset_hash,          //!< Array of indices which indicate which dataset each of the components in base_uijs corresponds to (length n_base).
      const symArr1d &atomic_uijs,                //!< Array of symmetric Uij matrices that describe the individual Uij for each atom (in all datasets) (length n_atoms).
      double weight_sum_of_amplitudes,            //!< weight for lasso-like penalty term of the sum of all base_amplitudes and the sum of the magnitudes of atomic_uijs.
      double weight_sum_of_squared_amplitudes,    //!< weight for ridge regression-like penalty term of the sum of all squared base_amplitudes and the sum of the squared magnitudes of atomic_uijs.
      double weight_sum_of_amplitudes_squared     //!< weight for the lasso-like term squared.
      );

    /*!
      Returns the current set of amplitudes - concatenate the current base_amplitudes with the magnitudes of atomic_uijs.
      Returns 1D array of length (n_base + n_atoms)
    */
    dblArr1d getCurrentAmplitudes();

    /*!
      Set the current amplitudes.
      Requires 1D array of length (n_base + n_atoms)
    */
    void setCurrentAmplitudes(const dblArr1d &values);

    /*!
      Prints the current amplitudes.
    */
    void printCurrentAmplitudes();

    /*!
      Set which datasets are to be used for the optimisation of the magnitudes of the atomic_uijs.
      Requires 1D boolean array of indices of length n_datasets
    */
    void setAtomicOptimisationMask(const blnArr1d &mask);

    /*!
      Calculate function and gradients for the current set of amplitudes.
    */
    bp::tuple computeFunctionalAndGradients();

  private:

    // Input variables
    //
    const symArrNd target_uijs;
    const dblArrNd target_weights;      //!< Weights on input atoms
    af::shared<symArr1d*> base_u;       //!< list of list of base Uij elements
    af::shared<selArr1d*> base_i;       //!< list of list of indices mapping base_uij elements to atoms
    const selArr1d base_dataset_hash;   //!< maps base elements to datasets
    const symArr1d atomic_uijs;         //!< sym Uij for each of the atoms
    blnArr1d atomic_mask;               //!< mask of which datasets are used to optimise the atomic amplitudes

    // Weights
    //
    double weight_sum_amp;              //!< Weight term on sum(amplitudes)
    double weight_sum_sqr_amp;          //!< Weight term on sum(amplitudes^2)
    double weight_sum_amp_sqr;          //!< Weight term on sum(amplitudes)^2

    // Quantities calculated from input variables
    //
    const size_t n_dst, n_atm, n_base, n_total;
    int atomic_mask_total;              //!< Number of datasets to use for atomic optimisation
    dblArr1d dataset_weights;           //!< Average weight for each dataset
    double average_dataset_weight;      //!< Average weight over all datasets

    // Internal intermediate variables
    //
    dblArr1d initial_amplitudes;        //!< Copy of the input amplitudes
    dblArr1d current_amplitudes;        //!< Array of the current amplitudes
    symArrNd total_uijs;                //!< Array of the current total uijs (sum over all components)

    // Optimisation/output variables
    //
    double functional;                  //!< Current functional
    dblArr1d gradients;                 //!< Current gradients

    size_t n_call;

    void zero();                        //!< Reinitialises internal intermediate variables for next calculation
    void sanitiseCurrentAmplitudes();   //!< Zeros negative amplitudes
    void calculateDatasetWeights();     //!< Calculate the average Uij weight of each dataset
    void calculateTotalUijs();          //!< Apply the amplitudes to the base_uijs and sum over all elements
    void calculateFGLeastSquares();     //!< Calculate the least-squares contributions to the functional and the gradients
    void calculateFGSumAmp();           //!< Calculate the lasso-like sum-of-amplitudes contributions to the functional and the gradients
    void calculateFGSumSqrAmp();        //!< Calculate the ridge regression-like sum-of-squared-amplitudes contributions to the functional and the gradients
    void calculateFGSumAmpSqr();        //!< Calculate the sum-of-amplitudes-squared contributions to the functional and the gradients

};

} } } // close namepsace mmtbx/tls/optimise

#endif
