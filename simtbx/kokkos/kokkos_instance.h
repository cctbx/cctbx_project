#ifndef SIMTBX_KOKKOS_INSTANCE_H
#define SIMTBX_KOKKOS_INSTANCE_H

#include "kokkostbx/kokkos_types.h"

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
    static bool m_isInitialized;
    static int m_instances;
    bool m_isFinalized = false;
    int m_deviceID = -1;


};
} // Kokkos
} // simtbx
#endif // SIMTBX_KOKKOS_INSTANCE_H
