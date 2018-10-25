from __future__ import absolute_import, division, print_function
from libtbx import easy_run

def marked_for_commit():
  lines = easy_run.fully_buffered(
    command="svn status").raise_if_errors().stdout_lines
  for li in lines:
    status = li[0] # c.f. [1]
    if status in ('A','M','R'):
      path = li[5:].strip()
      yield path

### References:

### [1] http://svnbook.red-bean.com/en/1.1/ch03s05.html#svn-ch-3-sect-5.3.1
### "In this output format svn status prints five columns of characters,
### followed by several whitespace characters, followed by a file or directory
### name. The first column tells the status of a file or directory
### and/or its contents. ...."

if (__name__ == "__main__"):
  print("\n".join(marked_for_commit()))
