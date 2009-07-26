#ifndef SCITBX_RIGID_BODY_BODY_T_H
#define SCITBX_RIGID_BODY_BODY_T_H

#include <scitbx/rigid_body/joint_t.h>
#include <scitbx/array_family/versa.h>

namespace scitbx { namespace rigid_body {

  template <typename FloatType>
  struct body_t
  {
    typedef FloatType ft;

    unsigned number_of_sites;
    ft sum_of_masses;
    shared_ptr<alignment_t<ft> > alignment;
    af::versa<ft, af::mat_grid> i_spatial;
    shared_ptr<joint_t<ft> > joint;
    rotr3<ft> cb_tree;
    int parent;

    virtual
    ~body_t() {}

    virtual
    af::const_ref<ft>
    qd() const = 0;
  };

}} // namespace scitbx::rigid_body

#endif // GUARD
