#include "simtbx/kokkos/kokkos_instance.h"

using Kokkos::InitializationSettings;
using Kokkos::initialize;
using Kokkos::finalize;

namespace simtbx {
namespace Kokkos {

  bool kokkos_instance::m_isInitialized = false;
  int kokkos_instance::m_instances = 0;

  kokkos_instance::kokkos_instance() {
    printf("NO OPERATION, NO DEVICE NUMBER");
  }

  kokkos_instance::kokkos_instance(int const& t_deviceID) {
    if (!m_isInitialized) {
      initialize(InitializationSettings()
                      .set_device_id(t_deviceID));

      m_isInitialized = true;
      m_isFinalized = false;
      m_deviceID = t_deviceID;
    }
    ++m_instances;
  }

  int
  kokkos_instance::get_deviceID() const {
    return m_deviceID;
  }

  void
  kokkos_instance::finalize_kokkos() {
    finalize();
    m_isFinalized = true;
  }

  kokkos_instance::~kokkos_instance() {
    --m_instances;
    if (!m_isFinalized && m_instances<1) { finalize(); }
  }


} // Kokkos
} // simtbx
