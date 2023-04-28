#ifndef CCTBX_OTHER_RESTRAINTS_SUMP_H
#define CCTBX_OTHER_RESTRAINTS_SUMP_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cctbx/restraints.h>

namespace cctbx { namespace other_restraints {

  struct sump_proxy {
    //! Default constructor. Some data members are not initialized!
    sump_proxy() {}

    /*! Constructor
    * arguments:
    * @param i_seqs - indices of components into site occupancies
    * @param coefficients - multiplication factor for particular component
    * @param weight - the restraint weight
    * @param target - the target value of sum_over_i(occu_i*coefficient_i)
    * @param labels - could be empty - then site labels or indices will be used
    * @param all_i_seqs - use this when restrain is over dependent occupancy
    *   constraint
    * @param group_sizes - defines structure of all_i_seq
    *   group_sizes.size() = i_seqs.size() &&
    *     sum_over_i(group_sizes_i) = all_i_seqs.size()
    */
    sump_proxy(
      af::shared<unsigned> const& i_seqs,
      af::shared<double> const& coefficients,
      double weight, double target,
      af::shared<std::string> const& labels,
      af::shared<unsigned> const& all_i_seqs,
      af::shared<unsigned> const& group_sizes)
    : i_seqs(i_seqs),
      coefficients(coefficients),
      weight(weight),
      target(target),
      labels(labels),
      all_i_seqs(all_i_seqs),
      group_sizes(group_sizes)
    {
      CCTBX_ASSERT(i_seqs.size() == coefficients.size());
      CCTBX_ASSERT(labels.size() == 0 || labels.size() == i_seqs.size());
      CCTBX_ASSERT(group_sizes.size() == i_seqs.size());
      CCTBX_ASSERT(all_i_seqs.size() >= i_seqs.size());
      CCTBX_ASSERT(af::sum(group_sizes.const_ref()) == all_i_seqs.size());
    }

    //! Indices into array of occupancies.
    af::shared<unsigned> i_seqs;
    af::shared<double> coefficients;
    double weight, target;

    af::shared<std::string> labels;
    af::shared<unsigned> all_i_seqs;
    af::shared<unsigned> group_sizes;
  };

  class sump {
  public:
    //! Constructor.
    sump(
      af::shared<double> const& occupancies,
      af::shared<double> const& coefficients,
      double weight, double target)
    : coefficients(coefficients),
      weight(weight),
      target(target)
    {
      double sum = 0;
      for (size_t i = 0; i < occupancies.size(); i++) {
        sum += occupancies[i] * coefficients[i];
      }
      delta_ = sum - target;
    }

    //! Constructor.
    sump(
      af::shared<double> const& occupancies,
      sump_proxy const& proxy)
      : coefficients(proxy.coefficients),
      weight(proxy.weight),
      target(proxy.target)
    {
      double sum = 0;
      for (size_t i = 0; i < proxy.i_seqs.size(); i++) {
        sum += occupancies[proxy.i_seqs[i]] * coefficients[i];
      }
      delta_ = sum - target;
    }

    //! weight * delta[i]**2.
    double
    residual() const {
      return weight * scitbx::fn::pow2(delta_);
    }


    void linearise(
      cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
      cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
      af::shared<unsigned> const& i_seqs) const
    {
      std::size_t row_i = linearised_eqns.next_row();
      linearised_eqns.weights[row_i] = weight;
      linearised_eqns.deltas[row_i] = delta_;
      for (std::size_t i = 0; i < i_seqs.size(); i++) {
        cctbx::xray::parameter_indices const &ids_i = parameter_map[i_seqs[i]];
        if (ids_i.occupancy == -1) {
          continue;
        }
        linearised_eqns.design_matrix(row_i, ids_i.occupancy) = coefficients[i];
      }
    }

    double delta() const { return delta_; }

  protected:
    double delta_;
  public:
    af::shared<double> coefficients;
    double weight, target;
  };

}} // namespace cctbx::adp_restraints

#endif // CCTBX_OTHER_RESTRAINTS_SUMP_H
