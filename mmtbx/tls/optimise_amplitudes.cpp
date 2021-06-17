#include <cmath>
#include <iostream>

#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>

#include <mmtbx/tls/utils.h>
#include <mmtbx/tls/optimise_amplitudes.h>

namespace mmtbx { namespace tls { namespace optimise {

template <typename T>
T find_max(const af::shared<T> &array)
{
  T max_val = array[0];
  for (size_t i=1; i<array.size(); i++)
  {
    if (max_val < array[i]) { max_val = array[i]; }
  }
  return max_val;
};

bool is_zero(const sym s, double tol=1e-12)
{
  if (
    (std::abs(s[0])<tol) &&
    (std::abs(s[1])<tol) &&
    (std::abs(s[2])<tol) &&
    (std::abs(s[3])<tol) &&
    (std::abs(s[4])<tol) &&
    (std::abs(s[5])<tol)
    )
  {
    return true;
  }
  return false;
}

//! Main constructor
MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
    const symArrNd &target_uijs,
    const dblArrNd &target_weights,
    const dblArr1d &base_amplitudes,
    const bp::list &base_uijs,
    const bp::list &base_atom_indices,
    const selArr1d &base_dataset_hash,
    const symArr1d &atomic_uijs,
    double weight_sum_of_amplitudes,
    double weight_sum_of_squared_amplitudes,
    double weight_sum_of_amplitudes_squared ) :
  target_uijs(target_uijs),
  target_weights(target_weights),
  base_dataset_hash(base_dataset_hash),
  atomic_uijs(atomic_uijs),
  weight_sum_amp(weight_sum_of_amplitudes),
  weight_sum_sqr_amp(weight_sum_of_squared_amplitudes),
  weight_sum_amp_sqr(weight_sum_of_amplitudes_squared),
  n_dst(target_uijs.accessor().all()[0]),
  n_atm(target_uijs.accessor().all()[1]),
  n_base(base_amplitudes.size()),
  n_total(base_amplitudes.size() + atomic_uijs.size()),
  atomic_mask_total(0),
  n_call(0)
{
  std::ostringstream errMsg;

  // ==========================
  // Check target uijs
  // ==========================
  if (target_uijs.accessor().nd() != 2) {
    errMsg << "invalid target_uijs: must be 2-dimensional flex array (currently " << target_uijs.accessor().nd() << ")";
    throw std::invalid_argument( errMsg.str() );
  }

  // ==========================
  // Check weights
  // ==========================
  if (target_uijs.accessor().nd() != target_weights.accessor().nd())
  {
    errMsg << "invalid dimension of target_weights (dimension " << target_weights.accessor().nd() << "): must be same dimension as target_uijs (dimension " << target_uijs.accessor().nd() << ")";
    throw std::invalid_argument( errMsg.str() );
  }
  for (size_t i=0; i<target_uijs.accessor().nd(); i++)
  {
    if (target_uijs.accessor().all()[i] != target_weights.accessor().all()[i])
    {
      errMsg << "incompatible dimension of target_weights (axis " << i << "): must be same size as target_uijs (" << target_weights.accessor().all()[i] << " != " << target_uijs.accessor().all()[i] <<")";
      throw std::invalid_argument( errMsg.str() );
    }
  }

  // ==========================
  // Check base_amplitudes, base_uijs and base_atom_indices (common error message)
  // ==========================
  errMsg << "invalid input base components. "
    << "base_amplitudes (length " << base_amplitudes.size()
    << "), base_uijs (length " << bp::len(base_uijs)
    << ") and base_atom_indices (length " << bp::len(base_atom_indices)
    << ") must all be the same length";
  if (base_amplitudes.size() != n_base) { throw std::invalid_argument( errMsg.str() ); }
  if (bp::len(base_uijs) != n_base) { throw std::invalid_argument( errMsg.str() ); }
  if (bp::len(base_atom_indices) != n_base) { throw std::invalid_argument( errMsg.str() ); }
  errMsg.str(""); // clear the previous error message

  // ==========================
  // Check weights
  // ==========================
  if (weight_sum_amp < 0.0)     { throw std::invalid_argument( "weight_sum_of_amplitudes must be positive" ); }
  if (weight_sum_sqr_amp < 0.0) { throw std::invalid_argument( "weight_sum_of_squared_amplitudes must be positive" ); }
  if (weight_sum_amp_sqr < 0.0) { throw std::invalid_argument( "weight_sum_of_amplitudes_squared must be positive" ); }

  // ==========================
  // Unpack base_uijs and base_atom_indices
  // ==========================
  base_u.reserve(bp::len(base_uijs));
  base_i.reserve(bp::len(base_atom_indices));
  for (std::size_t i = 0; i < n_base; ++i)
  {
    symArr1d* atm_u = new symArr1d(bp::extract<symArr1d>(base_uijs[i]));
    selArr1d* atm_i = new selArr1d(bp::extract<selArr1d>(base_atom_indices[i]));

    if (atm_u->size() != atm_i->size())
    {
      errMsg << "incompatible pair (element " << i << ") in base_uijs/base_atom_indices: pairwise elements must be the same length (" << atm_u->size() << " and " << atm_i->size() << ")";
      throw std::invalid_argument( errMsg.str() );
    }
    if (find_max(*atm_i) >= n_atm)
    {
      errMsg << "invalid selection in base_atom_indices (" << find_max(*atm_i) << "): attempting to select atom outside of array (size " << n_atm << ")";
      throw std::invalid_argument( errMsg.str() );
    }

    base_u.push_back(atm_u);
    base_i.push_back(atm_i);
  }

  // ==========================
  // Check dataset hash
  // ==========================
  if (base_dataset_hash.size() != n_base) {
    errMsg << "invalid base_dataset_hash (length " << base_dataset_hash.size() << "): must be same length as base_amplitudes, base_uijs & base_atom_indices (length " << base_amplitudes.size() << ")";
    throw std::invalid_argument( errMsg.str() );
  }
  if (find_max(base_dataset_hash) >= n_dst)
  {
    errMsg << "invalid value in base_dataset_hash (" << find_max(base_dataset_hash) << "): attempts to select element outside range of target_uijs (size " << n_dst << ")";
    throw std::invalid_argument( errMsg.str() );
  }
  for (size_t i_dst=0; i_dst < n_dst; i_dst++)
  {
    bool found = false;
    for (size_t i_base=0; i_base < n_base; i_base++)
    {
      if (base_dataset_hash[i_base] == i_dst)
      {
        found = true;
        break;
      }
    }
    if (found == false)
    {
      errMsg << "Dataset index " << i_dst << " is not present in base_dataset_hash -- this dataset has no base elements associated with it.";
      throw std::invalid_argument( errMsg.str() );
    }
  }

  // ==========================
  // Check atomic uijs
  // ==========================
  if (atomic_uijs.size() != n_atm)
  {
    errMsg << "invalid size of atomic_uijs (" << atomic_uijs.size() << "): must match 2nd dimension of target_uijs (" << n_atm << ")";
    throw std::invalid_argument( errMsg.str() );
  }

  // ==========================
  // Assign class members
  // ==========================

  // Initialise atomic amplitude array
  dblArr1d res_amplitudes(n_atm, 1.0);

  // Concatenate amplitude arrays
  initial_amplitudes.reserve(n_total);
  std::copy(base_amplitudes.begin(), base_amplitudes.end(), std::back_inserter(initial_amplitudes));
  std::copy(res_amplitudes.begin(), res_amplitudes.end(), std::back_inserter(initial_amplitudes));

  // Copy to current values
  current_amplitudes = dblArr1d(initial_amplitudes);

  // Total Uijs (summed over levels) (datasets * atoms)
  total_uijs = symArrNd(af::flex_grid<>(n_dst, n_atm), sym(0.,0.,0.,0.,0.,0.));

  // Calculate dataset weights (average target weight over each dataset)
  calculateDatasetWeights();

  // Initialise blank atomic mask
  setAtomicOptimisationMask(blnArr1d(n_dst, true));

  // Initialise output
  functional = 0.0;
  gradients = dblArr1d(n_total, 0.0);

}

dblArr1d MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::getCurrentAmplitudes()
{
  return current_amplitudes;
}

void MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::setCurrentAmplitudes(const dblArr1d &values)
{
  if (values.size() != current_amplitudes.size()) {
    std::ostringstream errMsg;
    errMsg << "Input array (size " << values.size() << ") must be the same length as current_amplitudes (size " << current_amplitudes.size() << ")";
    throw std::invalid_argument( errMsg.str() );
  }
  for (size_t i=0; i<current_amplitudes.size(); i++)
  {
    current_amplitudes[i] = values[i];
  }
}

void MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::printCurrentAmplitudes()
{
  for (size_t i=0; i<current_amplitudes.size(); i++)
  {
    std::cout << i << " - " << current_amplitudes[i] << std::endl;
  }
}

void MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::setAtomicOptimisationMask(const blnArr1d &mask)
{
  if (mask.size() != n_dst) {
    std::ostringstream errMsg;
    errMsg << "Input array (size " << mask.size() << ") must be the same length as number of datasets (" << n_dst << ")";
    throw std::invalid_argument( errMsg.str() );
  }
  atomic_mask = blnArr1d(mask);
  atomic_mask_total = 0;
  for (size_t i=0; i<atomic_mask.size(); i++)
  {
    atomic_mask_total += (int)atomic_mask[i];
  }
}

void MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::calculateDatasetWeights()
{
  // Overall weights for each dataset
  dataset_weights = dblArr1d(n_dst, 0.0);
  average_dataset_weight = 0.0;
  for (size_t i_dst=0; i_dst<n_dst; i_dst++)
  {
    for (size_t i_atm=0; i_atm<n_atm; i_atm++)
    {
      dataset_weights[i_dst] += target_weights(i_dst,i_atm) / (double)n_atm;
    }
    average_dataset_weight += dataset_weights[i_dst] / (double)n_dst;
  }
}

bp::tuple MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::computeFunctionalAndGradients()
{
  // Reset everything
  zero();

  n_call++;

  // Reset negative amplitudes
  //
  sanitiseCurrentAmplitudes();

  // Apply amplitudes to base uijs
  //
  calculateTotalUijs();

  // Calculate least-squares component of target function
  //
  calculateFGLeastSquares();

  // Penalty: Sum(amplitudes)
  if (weight_sum_amp > 0.0)
  {
    calculateFGSumAmp();
  }

  // Penalty: Sum(amplitudes^2)
  if (weight_sum_sqr_amp > 0.0)
  {
    calculateFGSumSqrAmp();
  }

  // Penalty: (Sum(amplitudes))^2
  if (weight_sum_amp_sqr > 0.0)
  {
    calculateFGSumAmpSqr();
  }

  // Return as python tuple for optimiser
  //
  return bp::make_tuple(functional, gradients);
}

void MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::zero()
{
  // Reset functional and gradients
  functional = 0.0;
  memset(&gradients[0], 0.0, sizeof(double) * gradients.size());

  // Zero-out the level uijs
  memset(&total_uijs[0], 0.0, sizeof(sym) * total_uijs.size());
}

void MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::sanitiseCurrentAmplitudes()
{
  // Set negative amplitudes to zero
  for (size_t i_opt=0; i_opt<current_amplitudes.size(); i_opt++)
  {
    if (current_amplitudes[i_opt] < 0.0)
    {
      // Constrain the value to be zero (and thereby negate any benefit to the functional that could have been gained)
      current_amplitudes[i_opt] = 0.0;
    }
  }
}

void MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::calculateTotalUijs()
{
  // ==========================================================
  // Sum over amplitudes to generate totals - !!! BASE TERMS !!!
  //
  for (size_t i_base=0; i_base<n_base; i_base++)
  {
    symArr1d &base_u_atom = *(base_u[i_base]);
    selArr1d &base_i_atom = *(base_i[i_base]);
    size_t i_dst = base_dataset_hash[i_base];
    size_t i_opt = i_base;

    // Iterate through atoms associated with this base element
    for (size_t i_atm_x=0; i_atm_x<base_u_atom.size(); i_atm_x++)
    {
      // apply multipliers
      sym m = current_amplitudes[i_opt] * base_u_atom[i_atm_x];
      // Skip if null
      if (is_zero(m))
      {
        continue;
      }
      // Add to total uijs
      total_uijs(i_dst, base_i_atom[i_atm_x]) += m;
    }
  }
  //
  // Sum over amplitudes to generate totals - !!! ATOMIC TERMS !!!
  //
  for (size_t i_atm=0; i_atm<n_atm; i_atm++)
  {
    size_t i_opt = n_base + i_atm;
    sym m = current_amplitudes[i_opt] * atomic_uijs[i_atm];
    for (size_t i_dst=0; i_dst<n_dst; i_dst++)
    {
      total_uijs(i_dst, i_atm) += m;
    }
  }
  // ==========================================================
}

void MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::calculateFGLeastSquares()
{
  // Normalisation term (none)
  double norm_all = 1.;

  // Normalisation term for atomic level
  double norm_res = 0.0; // if unchanged, atomic level will not be optimised
  if (atomic_mask_total > 0)
  {
    // Upweights by term n_all/n_calc to immitate being calculated over all datasets
    norm_res = norm_all * (double)(n_dst) / (double)(atomic_mask_total);
  }

  // ==========================================================
  // Calculate functional and gradients
  for (size_t i_dst=0; i_dst<n_dst; i_dst++)
  {
    // Extract weights for this dataset
    dblArr1d d_wgts(&target_weights(i_dst,0), &target_weights(i_dst, n_atm));
    // Extract total uijs for this dataset
    symArr1d d_total(&total_uijs(i_dst,0), &total_uijs(i_dst, n_atm));

    // Calculate difference to target
    symArr1d d_diffs(n_atm); // all elements will be populated
    for (size_t i_atm=0; i_atm<n_atm; i_atm++)
    {
      sym d_diff_i = target_uijs(i_dst,i_atm) - total_uijs(i_dst,i_atm);
      d_diffs[i_atm] = d_diff_i;

      // Calculate functional -- least-squares component
      for (size_t i_elem=0; i_elem<6; i_elem++)
      {
        functional += (
          d_diff_i[i_elem] *
          d_diff_i[i_elem] *
          d_wgts[i_atm] *
          norm_all
          );
      }
    }

    // Calculate gradients -- BASE TERMS
    for (size_t i_base=0; i_base<n_base; i_base++)
    {
      // Only calculate if this base term is related to this dataset
      if (i_dst != base_dataset_hash[i_base])
      {
        continue;
      }

      // i_ relative to full length of amplitudes
      size_t i_opt = i_base;

      // Extract the base elements
      symArr1d &base_u_atom = *(base_u[i_base]);
      selArr1d &base_i_atom = *(base_i[i_base]);

      // Iterate through the atoms in the base element
      for (size_t i_atm_x=0; i_atm_x<base_i_atom.size(); i_atm_x++)
      {
        // i_atm relative to n_atm
        size_t i_atm = base_i_atom[i_atm_x];

        // Extract the uij for this base element
        sym base_sym = base_u_atom[i_atm_x];

        // Skip this atom if the base_uij is zero at this position
        if (is_zero(base_sym))
        {
          continue;
        }

        // Corresponding u_diff
        sym diff_sym = d_diffs[i_atm];

        // Gradient from least-squares
        // -2*base*diffs*wgts
        for (size_t i_elem=0; i_elem<6; i_elem++)
        {
          gradients[i_opt] += (
            -2.0 *
            base_sym[i_elem] *
            diff_sym[i_elem] *
            d_wgts[i_atm] *
            norm_all
            );
        }
      }
    }

    // Calculate gradients -- ATOMIC TERMS
    if (atomic_mask[i_dst])
    {
      for (size_t i_atm=0; i_atm<n_atm; i_atm++)
      {
        size_t i_opt = n_base + i_atm;

        // Extract the uij for this base element
        sym base_sym = atomic_uijs[i_atm];

        // Skip this atom if the base_uij is zero at this position
        if (is_zero(base_sym))
        {
          continue;
        }

        // Corresponding u_diff
        sym diff_sym = d_diffs[i_atm];

        // Gradient from least-squares
        // -2*base*diffs*wgts
        for (size_t i_elem=0; i_elem<6; i_elem++)
        {
          gradients[i_opt] += (
            -2.0 *
            base_sym[i_elem] *
            diff_sym[i_elem] *
            d_wgts[i_atm] *
            norm_res
            );
        }
      }
    }
  }
}

void MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::calculateFGSumAmp()
{

  // common multiplier for all terms
  double multiplier_all = weight_sum_amp;
  double multiplier_res = multiplier_all * (double)n_dst;

  //
  // calculate functional and gradients
  //
  // base terms
  //
  for (size_t i_base=0; i_base<n_base; i_base++)
  {
    // Extract dataset
    size_t i_dst = base_dataset_hash[i_base];
    // weight by dataset weight
    functional += multiplier_all * dataset_weights[i_dst] * current_amplitudes[i_base];
    gradients[i_base] += multiplier_all * dataset_weights[i_dst];
  }
  //
  // atomic terms
  //
  double res_amp;
  for (size_t i_atm=0; i_atm<n_atm; i_atm++)
  {
    size_t i_opt = n_base + i_atm;
    sym m = current_amplitudes[i_opt] * atomic_uijs[i_atm];
    res_amp = (m[0] + m[1] + m[2]) / 3.0;
    // weight by average_dataset_weight
    functional += multiplier_res * average_dataset_weight * res_amp;
    gradients[i_opt] += multiplier_res * average_dataset_weight;
  }
}

void MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::calculateFGSumSqrAmp()
{

  // common multiplier for all terms
  double multiplier_all = weight_sum_sqr_amp;
  double multiplier_res = multiplier_all * (double)n_dst;

  //
  // calculate functional and gradients
  //
  // base terms
  //
  for (size_t i_base=0; i_base<n_base; i_base++)
  {
    // Extract dataset
    size_t i_dst = base_dataset_hash[i_base];
    // weight by dataset weight
    functional += multiplier_all * dataset_weights[i_dst] * current_amplitudes[i_base] * current_amplitudes[i_base];
    gradients[i_base] += 2.0 * multiplier_all * dataset_weights[i_dst] * current_amplitudes[i_base];
  }
  //
  // atomic terms
  //
  double res_amp;
  for (size_t i_atm=0; i_atm<n_atm; i_atm++)
  {
    size_t i_opt = n_base + i_atm;
    sym m = current_amplitudes[i_opt] * atomic_uijs[i_atm];
    res_amp = (m[0] + m[1] + m[2]) / 3.0;
    // weight by average_dataset_weight
    functional += multiplier_res * average_dataset_weight * res_amp * res_amp;
    gradients[i_opt] += 2.0 * multiplier_res * average_dataset_weight * res_amp;
  }
}

void MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator::calculateFGSumAmpSqr()
{

  // common multiplier for all terms
  double multiplier_all = weight_sum_amp_sqr;

  //
  // calculate amplitude sums
  //
  // base terms
  //
  dblArr1d dst_amp_sum(n_dst, 0.0); // Sum of amplitudes
  for (size_t i_base=0; i_base<n_base; i_base++)
  {
    size_t i_dst = base_dataset_hash[i_base];
    dst_amp_sum[i_dst] += current_amplitudes[i_base];
  }
  //
  // atomic terms
  //
  double res_amp_sum = 0.0;
  for (size_t i_atm=0; i_atm<n_atm; i_atm++)
  {
    size_t i_opt = n_base + i_atm;
    sym m = current_amplitudes[i_opt] * atomic_uijs[i_atm];
    res_amp_sum += (m[0] + m[1] + m[2]) / 3.0;
  }

  //
  // Calculate functionals
  //
  for (size_t i_dst=0; i_dst<n_dst; i_dst++)
  {
    functional += (multiplier_all * dataset_weights[i_dst]) * ( dst_amp_sum[i_dst] + res_amp_sum ) * ( dst_amp_sum[i_dst] + res_amp_sum );
  }

  // Calculate gradients
  //
  // base terms
  //
  for (size_t i_base=0; i_base<n_base; i_base++)
  {
    size_t i_dst = base_dataset_hash[i_base];
    gradients[i_base] += 2.0 * (multiplier_all * dataset_weights[i_dst]) * ( dst_amp_sum[i_dst] + res_amp_sum );
  }
  //
  // atomic terms
  //
  for (size_t i_atm=0; i_atm<n_atm; i_atm++)
  {
    size_t i_opt = n_base + i_atm;
    for (size_t i_dst=0; i_dst<n_dst; i_dst++)
    {
      gradients[i_opt] += 2.0 * (multiplier_all * dataset_weights[i_dst]) * ( dst_amp_sum[i_dst] + res_amp_sum );
    }
  }
}

} } } // close namepsace mmtbx/tls/optimise
