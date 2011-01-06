#ifndef SMTBX_REFINEMENT_CONSTRAINTS_SCATTERER_PARAMETERS_H
#define SMTBX_REFINEMENT_CONSTRAINTS_SCATTERER_PARAMETERS_H

#include <boost/iterator/filter_iterator.hpp>
#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

/// Helper to present crystallographic parameters in a user-defined order
/** The class constraints::scatterer_parameters below features named parameters,
    whereas this class presents them in a list. The former may therefore use
    the latter to define an ordering.
 */
struct ordered_scatterer_parameters
  : public af::tiny<asu_parameter *, 3>
{
  typedef af::tiny<asu_parameter *, 3> base_t;

  struct is_variable {
    bool operator()(asu_parameter const *p) const {
      return p->is_variable();
    }
  };

  typedef boost::filter_iterator<is_variable, base_t::const_iterator>
          const_iterator;

  ordered_scatterer_parameters(asu_parameter *p0,
                               asu_parameter *p1,
                               asu_parameter *p2)
  : base_t(p0, p1, p2)
  {}

  /// Those iterators skip invariable parameters
  //@{
  const_iterator begin() const {
    return const_iterator(is_variable(), base_t::begin(), base_t::end());
  }

  const_iterator end() const {
    return const_iterator(is_variable(), base_t::end(), base_t::end());
  }
  //@}
};


/// The constraints::parameters corresponding to a single scatterer
struct scatterer_parameters
{
  typedef asu_parameter::scatterer_type scatterer_type;

  scatterer_type const *scatterer;
  asu_parameter *site, *occupancy, *u;

  scatterer_parameters() {}

  scatterer_parameters(scatterer_type const *scatterer)
    : scatterer(scatterer),
      site(0), occupancy(0), u(0)
  {}

  scatterer_parameters(scatterer_type const *scatterer,
                       asu_parameter *site,
                       asu_parameter *occupancy,
                       asu_parameter *u)
    : scatterer(scatterer),
      site(site), occupancy(occupancy), u(u)
  {}

  /// Parameters in the same order as the derivatives in grad Fc
  /** This shall abide to the convention of smtbx::structure_factors.
   */
  ordered_scatterer_parameters ordered() const {
    return ordered_scatterer_parameters(site, u, occupancy);
  }
};


/// Map grad Fc derivative indices to constraint::parameter component indices
/** The mapping is returned as an array p such that a component index j
    of a parameter used in the constraint system which correspond to the
    i-th element of grad Fc comes as p[i] = j
 */
af::shared<std::size_t>
mapping_to_grad_fc(af::const_ref<scatterer_parameters> const &params) {
  af::shared<std::size_t> result((af::reserve(5*params.size()))); // heuristic
  for (std::size_t i=0; i<params.size(); ++i) {
    BOOST_FOREACH (asu_parameter const *p, params[i].ordered()) {
      index_range r = p->component_indices_for(params[i].scatterer);
      SMTBX_ASSERT(r.is_valid())(params[i].scatterer->label);
      for (std::size_t j=r.first(); j<r.last(); ++j) result.push_back(j);
    }
  }
  return result;
}


/// Write annotations for each components of each scatterers in grad Fc order
void
write_component_annotations(af::const_ref<scatterer_parameters> const &params,
                            std::ostream &output)
{
  for (std::size_t i=0; i<params.size(); ++i) {
    BOOST_FOREACH (asu_parameter const *p, params[i].ordered()) {
      p->write_component_annotations_for(params[i].scatterer, output);
    }
  }
}

}}}
#endif // GUARD
