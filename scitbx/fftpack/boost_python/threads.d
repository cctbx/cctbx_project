#pragma D option quiet

#define HASH(k) k*(k+3) % 63

#define SECTION 1

#if SECTION == 1

pid$target::*real_to_complex_3d*forward?scitbx??af??ref?double*:entry
{
  profiling = 1;
  start = timestamp;
}

pid$target::*real_to_complex_3d*forward?scitbx??af??ref?double*:return
{
  end = timestamp;
  profiling = 0;
}

pid$target::GOMP_parallel_start:entry
{
  start_parallel = timestamp;
}

pid$target::GOMP_parallel_end:return
{
  end_parallel = timestamp;
}

pid$target::*real_to_complex*forward_compressed*:entry,
pid$target::*complex_to_complex*transform*:entry
{
  // arg0 = "this" pointer of the object these member functions are called upon
  @counts_per_thread[HASH(tid), probefunc, arg0] = count();
  @counts_per_cpu[cpu, probefunc, arg0] = count();
}

END {
  printf("\nTotal runtime: %i ns, parallelizable runtime: %i ns\n", 
         end - start, end_parallel - start_parallel);
  printf("\n ** Number of calls **\n");
  printf("\nPer thread\n");
  printa(@counts_per_thread);
  printf("\nPer cpu\n");
  printa(@counts_per_cpu);
  printf("\n** How threads used the cpus **\n");
}

profile:::profile-500
/pid == $target && profiling/
{
  @sample[cpu, HASH(tid)] = count();
}

#endif