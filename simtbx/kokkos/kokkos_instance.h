#ifndef SIMTBX_KOKKOS_INSTANCE_H
#define SIMTBX_KOKKOS_INSTANCE_H

#include "simtbx/kokkos/kokkos_types.h"

namespace simtbx {
namespace Kokkos {

class kokkos_instance {
  public:
    kokkos_instance();
    kokkos_instance(int const&);
    ~kokkos_instance();

    void finalize_kokkos();

    int get_deviceID() const;

  private:
    bool bFinalized = false;
    int deviceID = -1;

};
} // Kokkos
} // simtbx
#endif // SIMTBX_KOKKOS_INSTANCE_H
