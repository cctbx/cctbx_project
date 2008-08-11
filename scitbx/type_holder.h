#ifndef SCITBX_TYPE_HOLDER_H
#define SCITBX_TYPE_HOLDER_H

namespace scitbx {

    template <class T> struct type_holder { typedef T type; };

} // namespace scitbx

#endif // SCITBX_TYPE_HOLDER_H
