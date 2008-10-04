/*
This is best used as:
  libtbx.dtrace libtbx/dtrace/profile.d path/to/script.py

It is useful as it is to profile Python function calls but it is also a good
starting point for hacking to achieve more advanced needs (C++ functions profiling, statistics per thread, etc).
*/

#pragma D option quiet

BEGIN {
  cctbx = "cctbx_project/";
}

python$target:::function-entry
/   index(copyinstr(arg0), "/System") == -1
 && index(copyinstr(arg0), "/usr") == -1
/
{
  self->ts=timestamp;
  f = copyinstr(arg0);
  i = index(f, cctbx);
  self->file = i >= 0 ? substr(f, i + strlen(cctbx)) : f;
  self->function = copyinstr(arg1);
  self->line = arg2;
  self->traced=1;
}

python$target:::function-return
/ self->traced /
{
  @time[self->function, self->file, self->line]
    = sum((timestamp - self->ts)/1000);
  self->traced = 0;
}

#ifndef SHOW_TOP
#define SHOW_TOP 10
#endif

END {
  printf("Top function by time in micro-seconds\n");
  trunc(@time, SHOW_TOP);
  printa("%'@9i   %-35s%s (%i)\n", @time);
}
