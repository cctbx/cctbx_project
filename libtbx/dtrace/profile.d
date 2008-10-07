/*
This is best used as:
  libtbx.dtrace profile path/to/script.py

It may be useful as it is but DTrace is to tracing what the shell is to file
manipulation: best to write simple, short, focused scripts. So this script is 
best used as a starting demo!

Beware that recursive calls are not properly handled by the by.
*/

#pragma D option quiet

BEGIN {
  cctbx = "cctbx_project/";
}

python$target:::function-entry
{
  self->entry_time = vtimestamp;
  self->line = arg2;
  self->traced = 1;
}

python$target:::function-return
/ self->traced /
{
  this->file_path = copyinstr(arg0);
  this->file = basename(this->file_path);
  this->directory = dirname(this->file_path);
  this->function = copyinstr(arg1);
  @time[this->function, this->file, self->line, this->directory]
    = sum(vtimestamp - self->entry_time);
  self->traced = 0;
}

pid$target::scitbx*:entry,
pid$target::cctbx*:entry,
pid$target::mmtbx*:entry,
pid$target::smtbx*:entry
{
  self->cpp_entry_time = vtimestamp;
  self->cpp_traced = 1;
}

pid$target::scitbx*:return,
pid$target::cctbx*:return,
pid$target::mmtbx*:return,
pid$target::smtbx*:return
/ self->cpp_traced /
{
  @time[probefunc, probemod, 0, "-"] = sum(vtimestamp - self->cpp_entry_time);
  self->cpp_traced = 0;
}

#ifndef SHOW_TOP
#define SHOW_TOP 10
#endif

END {
  printf("Top function by time in micro-seconds\n");
  trunc(@time, SHOW_TOP);
  normalize(@time, 1000);
  printa("%'@9i   %-35s%-35s (%i: %s)\n", @time);
}
