#include "simtbx/kokkos/kokkos_instance.h"

using Kokkos::InitArguments;
using Kokkos::initialize;
using Kokkos::finalize;

namespace simtbx {
namespace Kokkos {

  kokkos_instance::kokkos_instance() {
    printf("NO OPERATION, NO DEVICE NUMBER");
  }

  kokkos_instance::kokkos_instance(int const& t_deviceID) {
    InitArguments kokkos_init;
    kokkos_init.device_id = t_deviceID;

    initialize(kokkos_init);
    bFinalized = false;
    deviceID = t_deviceID;
  }

  int
  kokkos_instance::get_deviceID() const {
    return deviceID;
  }

  void
  kokkos_instance::finalize_kokkos() {
    finalize();
    bFinalized = true;
  }

  kokkos_instance::~kokkos_instance() {
    if (!bFinalized) { finalize(); }
  }


} // Kokkos
} // simtbx
