#ifndef OMPTBX_LOCK_H
#define OMPTBX_LOCK_H

#include <omptbx/omp_or_stubs.h>

namespace omptbx {

  class lock
  {
    private:
      omp_lock_t omp_lock_;

    public:
      lock() { omp_init_lock(&omp_lock_); }

      ~lock() { omp_destroy_lock(&omp_lock_); }

      void
      set() { omp_set_lock(&omp_lock_); }

      void
      unset() { omp_unset_lock(&omp_lock_); }

    private:
      lock(const lock&);
      const lock& operator=(const lock&);
  };

} // namespace omptbx

#endif /* OMPTBX_LOCK_H */
