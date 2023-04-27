#ifndef CCTBX_OTHER_RESTRAINTS_SUMP_H
#define CCTBX_OTHER_RESTRAINTS_SUMP_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cctbx/restraints.h>

namespace cctbx { namespace other_restraints {

  struct sump_proxy {
    //! Default constructor. Some data members are not initialized!
    sump_proxy() {}

    //! Constructor.
    sump_proxy(
      af::shared<unsigned> const& i_seqs,
      af::shared<double> const& coefficients,
      double weight, double target)
    : i_seqs(i_seqs),
      coefficients(coefficients),
      weight(weight),
      target(target)
    {}

    //! Indices into array of sites.
    af::shared<unsigned> i_seqs;
    af::shared<double> coefficients;
    double weight, target;
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
