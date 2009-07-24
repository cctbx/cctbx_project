#ifndef SCITBX_RIGID_BODY_BODY_LIB_H
#define SCITBX_RIGID_BODY_BODY_LIB_H

#include <scitbx/rigid_body/joint_lib.h>

namespace scitbx { namespace rigid_body { namespace body_lib {

  template <typename FloatType>
  struct body_t
  {
    typedef FloatType ft;

    unsigned number_of_sites;
    ft sum_of_masses;
    boost::shared_ptr<joint_lib::alignment_t<ft> > alignment;
    af::versa<ft, af::mat_grid> i_spatial;
    boost::shared_ptr<joint_lib::joint_t<ft> > joint;
    rotr3<ft> cb_tree;
    int parent;

    virtual
    ~body_t() {}

    virtual
    af::const_ref<ft>
    qd() const = 0;
  };

}}} // namespace scitbx::rigid_body::tardy

#endif // GUARD
