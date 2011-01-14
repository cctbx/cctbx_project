#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/accessors/mat_grid.h>

namespace scitbx { namespace af {

  template <typename ElementType, int Rows, int Cols>
  class tiny_mat_ref : public ref<ElementType, mat_grid>
  {
  public:
    tiny_mat_ref(ElementType *storage)
    : ref<ElementType, mat_grid>(storage, mat_grid(Rows, Cols))
    {}
  };

  template <typename ElementType, int Rows, int Cols>
  class tiny_mat_const_ref : public const_ref<ElementType, mat_grid>
  {
  public:
    tiny_mat_const_ref(ElementType const *storage)
    : const_ref<ElementType, mat_grid>(storage, mat_grid(Rows, Cols))
    {}
  };

}}
