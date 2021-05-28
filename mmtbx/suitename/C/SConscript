import libtbx.load_env

Import("env_base", "env_etc")

env = env_base.Clone(LIBS=env_etc.libm)
if (libtbx.manual_date_stamp < 20090819):
  # XXX backward compatibility 2009-08-19
  env.Replace(CCFLAGS=env_etc.ccflags_base)
if (env_etc.compiler != "win32_cl"):
  env.Replace(LINK=env_base["CC"])

exe = env.Program(
  target=["#suitename/exe/suitename"],
  source=[
    "suitename.c",
    "suiteninit.c",
    "suiteninpt.c",
    "suitenout.c",
    "suitenscrt.c",
    "suitenutil.c"])
