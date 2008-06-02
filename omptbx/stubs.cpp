#if !defined(_OPENMP)

#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <omptbx/stubs.h>

extern "C" {

void
omp_set_num_threads(int num_threads)
{
}

int
omp_get_num_threads(void)
{
  return 1;
}

int
omp_get_max_threads(void)
{
  return 1;
}

int
omp_get_thread_num(void)
{
  return 0;
}

int
omp_get_num_procs(void)
{
  return 1;
}

void
omp_set_dynamic(int dynamic_threads)
{
}

int
omp_get_dynamic(void)
{
  return 0;
}

int
omp_in_parallel(void)
{
  return 0;
}

void
omp_set_nested(int nested)
{
}

int
omp_get_nested(void)
{
  return 0;
}

enum { UNLOCKED = -1, INIT, LOCKED };

void
omp_init_lock(omp_lock_t* lock)
{
  *lock = UNLOCKED;
}

void
omp_destroy_lock(omp_lock_t* lock)
{
  *lock = INIT;
}

void
omp_set_lock(omp_lock_t* lock)
{
  if (*lock == UNLOCKED) {
    *lock = LOCKED;
  }
  else if (*lock == LOCKED) {
    fprintf(stderr, "omptbx error: deadlock in using lock variable\n");
    exit(1);
  }
  else {
    fprintf(stderr, "omptbx error: lock not initialized\n");
    exit(1);
  }
}

void
omp_unset_lock(omp_lock_t* lock)
{
  if (*lock == LOCKED) {
    *lock = UNLOCKED;
  }
  else if (*lock == UNLOCKED) {
    fprintf(stderr, "omptbx error: lock not set\n");
    exit(1);
  }
  else {
    fprintf(stderr, "omptbx error: lock not initialized\n");
    exit(1);
  }
}

int
omp_test_lock(omp_lock_t* lock)
{
  if (*lock == UNLOCKED)
  {
    *lock = LOCKED;
    return 1;
  }
  else if (*lock == LOCKED)
  {
    return 0;
  }
  else
  {
    fprintf(stderr, "omptbx error: lock not initialized\n");
    exit(1);
  }
}

enum { NOOWNER = -1, MASTER = 0 };

void
omp_init_nest_lock(omp_nest_lock_t* nlock)
{
  nlock->owner = NOOWNER;
  nlock->count = 0;
}

void
omp_destroy_nest_lock(omp_nest_lock_t* nlock)
{
  nlock->owner = NOOWNER;
  nlock->count = UNLOCKED;
}

void
omp_set_nest_lock(omp_nest_lock_t* nlock)
{
  if (nlock->owner == MASTER && nlock->count >= 1) {
    nlock->count++;
  }
  else if (nlock->owner == NOOWNER && nlock->count == 0) {
    nlock->owner = MASTER;
    nlock->count = 1;
  }
  else {
    fprintf(stderr, "omptbx error: lock corrupted or not initialized\n");
    exit(1);
  }
}

void
omp_unset_nest_lock(omp_nest_lock_t* nlock)
{
  if (nlock->owner == NOOWNER && nlock->count >= 1) {
    nlock->count--;
    if (nlock->count == 0) {
      nlock->owner = NOOWNER;
    }
  }
  else if (nlock->owner == NOOWNER && nlock->count == 0) {
    fprintf(stderr, "omptbx error: lock not set\n");
    exit(1);
  }
  else {
    fprintf(stderr, "omptbx error: lock corrupted or not initialized\n");
    exit(1);
  }
}

int
omp_test_nest_lock(omp_nest_lock_t* nlock)
{
  omp_set_nest_lock(nlock);
  return nlock->count;
}

double
omp_get_wtime(void)
{
  /* This function does not provide a working
     wallclock timer. Replace it with a version
     customized for the target machine.
  */
  //return 0.0;
  throw std::runtime_error("omptbx: omp_get_wtime() not implemented.");
}

double
omp_get_wtick(void)
{
  /* This function does not provide a working
     clock tick function. Replace it with
     a version customized for the target machine.
   */
  // return 365. * 86400.;
  throw std::runtime_error("omptbx: omp_get_wtick() not implemented.");
}

} // extern "C"

#endif /* !defined(_OPENMP) */
