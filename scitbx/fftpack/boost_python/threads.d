#pragma D option quiet

#define HASH(k) k*(k+3) % 63

#define SECTION 1

#if SECTION == 1

pid$target::*real_to_complex*forward_compressed*:entry,
pid$target::*complex_to_complex*transform*:entry
{
  // arg0 = "this" pointer of the object these member functions are called upon
  @counts_per_thread[HASH(tid), probefunc, arg0] = count();
  @counts_per_cpu[cpu, probefunc, arg0] = count();
}

END {
  printf("\n ** Number of calls **\n");
  printf("\nPer thread\n");
  printa(@counts_per_thread);
  printf("\nPer cpu\n");
  printa(@counts_per_cpu);
  printf("\n** How threads used the cpus **\n");
}

pid$target::*real_to_complex_3d*forward*:entry
{
  profiling = 1;
}

profile:::profile-500
/pid == $target && profiling/
{
  @sample[HASH(tid)] = count();
}

#endif