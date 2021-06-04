#ifndef SIMTBX_KOKKOS_INSTANCE_H
#define SIMTBX_KOKKOS_INSTANCE_H

#include <Kokkos_Core.hpp>

namespace simtbx {
namespace Kokkos {

class kokkos_instance {
  public:
    kokkos_instance();
    kokkos_instance(int const&);
    ~kokkos_instance();

    void finalize_kokkos();

  private:
    bool bFinalized;
   
};
} // Kokkos
} // simtbx
#endif // SIMTBX_KOKKOS_INSTANCE_H
