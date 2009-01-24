#ifndef CCTBX_SGTBX_ABSTRACT_H 
#define CCTBX_SGTBX_ABSTRACT_H

#include <memory>

#include "cut.h"

namespace cctbx { namespace sgtbx { namespace asu {

  class abstract
  {
  public:
    typedef std::auto_ptr<abstract> ptr;

    virtual bool is_inside(const rvector3_t &p) const = 0;
    virtual ptr new_copy() const = 0;
    virtual ptr new_volume_only() const = 0;
    virtual size_type size() const = 0;
    virtual void change_basis(const change_of_basis_op &) =0;
    virtual void get_nth_plane(size_type i, cut &plane) const = 0;

    virtual ~abstract() {};
  }; // class abstract
  
  typedef abstract::ptr (*asu_func)();

  extern asu_func asu_table[230];


}}}
#endif

