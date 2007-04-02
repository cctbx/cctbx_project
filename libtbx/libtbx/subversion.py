import os

def marked_for_commit():
  pipe = os.popen('svn status')
  for li in pipe:
    status = li[0] # c.f. [1]
    path = li[5:].strip()
    if status in ('A','M','R'): yield path


### References:

### [1] http://svnbook.red-bean.com/en/1.1/ch03s05.html#svn-ch-3-sect-5.3.1
### "In this output format svn status prints five columns of characters,
### followed by several whitespace characters, followed by a file or directory
### name. The first column tells the status of a file or directory
### and/or its contents. ...."
