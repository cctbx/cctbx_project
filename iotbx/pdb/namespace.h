#ifndef IOTBX_PDB_NAMESPACE_H
#define IOTBX_PDB_NAMESPACE_H

#include <scitbx/array_family/shared.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/misc/positive_getitem_index.h>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

namespace iotbx {

//! Handling of files in PDB format.
namespace pdb {

  namespace af = scitbx::af;

  typedef scitbx::vec3<double> vec3;
  typedef scitbx::sym_mat3<double> sym_mat3;

  using scitbx::positive_getitem_index;

  using boost::shared_ptr;
  using boost::weak_ptr;

  struct stream_write
  {
    virtual void
    operator()(const char*, unsigned) = 0;

    virtual
    ~stream_write() {}
  };

}} // namespace iotbx::pdb

#endif // IOTBX_PDB_NAMESPACE_H
