#include <process.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

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

char*
getenv_certain(const char* argv0, const char* var_name)
{
  char* var_value = getenv(var_name);
  if (var_value == NULL) {
    fprintf(stderr, "%s: environment variable is not defined: %s\n",
      argv0, var_name);
    exit(1);
  }
  return var_value;
}

void*
malloc_certain(const char* argv0, unsigned long n)
{
  void* ptr;
  errno = 0;
  ptr = malloc(n);
  if (ptr == NULL) {
    fprintf(stderr, "%s: error allocating %lu bytes", argv0, n);
    if (errno) {
      fprintf(stderr, ": %s", strerror(errno));
    }
    fprintf(stderr, "\n");
    exit(2);
  }
  return ptr;
}

static const char* /* APPLICATION PREFIX interleaved with random digits */
application_prefix = "3A7P0P0L4I8C3A5T5I6O3N490P0R7E6F8I6X4";

int
main(int argc, char *const argv[])
{
  char env_var_name[256];
  char* env_python_exe;
  char* env_build;
  char* dispatcher;
  char** extended_argv;
  int n, i;
  const char* dispatcher_name = "\\dispatcher";
  _putenv("PYTHONHOME=");
  strcpy(env_var_name, application_prefix);
  strcat(env_var_name, "_PYTHON_EXE");
  env_python_exe = getenv_certain(argv[0], env_var_name);
  extended_argv = malloc_certain(argv[0], (argc + 3) * sizeof(char*));
  extended_argv[0] = env_python_exe;
  n = 1;
  i = 1;
  if (!is_python_exe(argv[0])) {
    strcpy(env_var_name, application_prefix);
    strcat(env_var_name, "_BUILD");
    env_build = getenv_certain(argv[0], env_var_name);
    dispatcher = malloc_certain(
      argv[0], strlen(env_build) + strlen(dispatcher_name) + 1);
    strcpy(dispatcher, env_build);
    strcat(dispatcher, dispatcher_name);
    extended_argv[n++] = dispatcher;
    i = 0;
  }
  for(;i<argc;i++,n++) {
    extended_argv[n] = malloc_certain(argv[0], strlen(argv[i]) + 3);
    strcpy(extended_argv[n], "\"");
    strcat(extended_argv[n], argv[i]);
    strcat(extended_argv[n], "\"");
  }
  extended_argv[n] = NULL;
  spawnv(P_WAIT, extended_argv[0], extended_argv);
  if (errno) {
    fprintf(stderr, "%s: error starting %s: %s\n",
      argv[0], extended_argv[0], strerror(errno));
    exit(3);
  }
  return 0;
}
