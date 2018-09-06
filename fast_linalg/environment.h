#pragma once

#include <fast_linalg/cblas.h>
#include <string>

namespace fast_linalg {

  /// Information about OpenBLAS
  /** All instances of this class obviously refer to the same information. */
  struct environment {
    /// Number of threads to be used in parallelisation
    int threads() { return openblas_get_num_threads(); }

    /// Set the number of threads to be used in parallelisation
    void set_threads(int n) { openblas_set_num_threads(n); }

    /// Number of physical cores on the machine
    /** E.g. a machine with 2 hexacore processors will reports 12 cores
        whether hyperthreading is enabled or not.
     */
    int physical_cores() { return openblas_get_num_procs(); }

    /// The family of the CPU (e.g. "Nehalem" for an Intel processor)
    std::string cpu_family() { return std::string(openblas_get_corename()); }

    /// A summary of the options used during the build of this library
    std::string build_config() { return std::string(openblas_get_config()); }
  };
}
