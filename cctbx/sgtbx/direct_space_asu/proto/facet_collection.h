#ifndef CCTBX_SGTBX_FACET_COLLECTION_H
#define CCTBX_SGTBX_FACET_COLLECTION_H

//! \cond

#include <memory>
#include <iosfwd>

#include "cut.h"

namespace cctbx { namespace sgtbx { namespace asu {

  class facet_collection
  {
  public:
    typedef std::auto_ptr<facet_collection> pointer;

    virtual bool is_inside(const rvector3_t &p) const = 0;
    virtual bool is_inside(const scitbx::af::int3 &num,
        const scitbx::af::int3 &den) const = 0;
    virtual bool is_inside(const scitbx::af::int3 &num) const = 0;
    virtual bool is_inside_shape_only(const scitbx::af::double3 &point,
        double tol) const = 0;
    virtual short where_is(const scitbx::af::int3 &num,
        const scitbx::af::int3 &den) const = 0;
    virtual short where_is(const scitbx::af::int3 &num) const = 0;
    virtual void optimize_for_grid(const scitbx::af::int3 &grid_size) = 0;
    virtual void get_optimized_grid_limits(scitbx::af::long3 &max_p) const = 0;
    virtual pointer new_copy() const = 0;
    virtual pointer new_shape_only() const = 0;
    virtual pointer new_shape_only_keep_inclusive_flag() const = 0;
    virtual size_type size() const = 0;
    virtual void change_basis(const change_of_basis_op &) =0;
    virtual void get_nth_plane(size_type i, cut &plane) const = 0;
    virtual double get_tolerance(const scitbx::af::double3 &tol3d) const = 0;
    virtual void print(std::ostream &os) const = 0;

    // DO NOT USE!!!
    virtual pointer add_face(const cut &face) const = 0;

    virtual ~facet_collection() {};
  }; // class facet_collection

  typedef facet_collection::pointer (*asu_func)();

  extern asu_func asu_table[230];


}}}
//! \endcond

#endif

