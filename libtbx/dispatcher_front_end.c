#include <process.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

int
is_exe(const char* path)
{
  int n, i;
  n = strlen(path);
  if (n < 4) return 0;
  n -= 4;
  for(i=0;i<4;i++,n++) {
    if (tolower(path[n]) != ".exe"[i]) return 0;
  }
  return 1;
}

int
is_python_exe(const char* path)
{
  int n, i, c;
  n = strlen(path);
  if (is_exe(path)) n -= 4;
  if (n < 6) return 0;
  n -= 6;
  for(i=0;i<6;i++,n++) {
    if (tolower(path[n]) != "python"[i]) return 0;
  }
  if (n > 6) {
    c = path[n-7];
    if (c != '.' && c != '\\') return 0;
  }
  return 1;
}

int
main(int argc, char *const argv[])
{
  char* libtbx_python_exe;
  char* libtbx_build;
  char* dispatcher;
  char** extended_argv;
  int n, i;
  const char* dispatcher_name = "\\dispatcher";
  _putenv("PYTHONHOME=");
  libtbx_python_exe = getenv("LIBTBX_PYTHON_EXE");
  assert(libtbx_python_exe != NULL);
  extended_argv = malloc((argc + 3) * sizeof(char*));
  assert(extended_argv != NULL);
  extended_argv[0] = libtbx_python_exe;
  n = 1;
  i = 1;
  if (!is_python_exe(argv[0])) {
    libtbx_build = getenv("LIBTBX_BUILD");
    assert(libtbx_build != NULL);
    dispatcher = malloc(strlen(libtbx_build) + strlen(dispatcher_name) + 1);
    assert(dispatcher != NULL);
    strcpy(dispatcher, libtbx_build);
    strcat(dispatcher, dispatcher_name);
    extended_argv[n++] = dispatcher;
    i = 0;
  }
  for(;i<argc;i++,n++) {
    extended_argv[n] = malloc(strlen(argv[i]) + 3);
    assert(extended_argv[n] != NULL);
    strcpy(extended_argv[n], "\"");
    strcat(extended_argv[n], argv[i]);
    strcat(extended_argv[n], "\"");
  }
  extended_argv[n] = NULL;
  spawnv(P_WAIT, extended_argv[0], extended_argv);
  if (errno) {
    perror(argv[0]);
    return 1;
  }
  return 0;
}
