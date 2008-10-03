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
/   stringof(copyin(arg0, 7)) != "/System" 
 && stringof(copyin(arg0, 4)) != "/usr"
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

END {
  printf("Top function by time in micro-seconds\n");
  trunc(@time, 10);
  printa("%'@9i   %-35s%s (%i)\n", @time);
}