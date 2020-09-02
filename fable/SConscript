import libtbx.load_env
Import("env_etc")

env_etc.fable_include = libtbx.env.dist_path("fable")

if not libtbx.env.module_is_installed("fable"):
  if (not env_etc.no_boost_python and hasattr(env_etc, "boost_adaptbx_include")):
    Import("env_no_includes_boost_python_ext")
    env = env_no_includes_boost_python_ext.Clone()
    env_etc.enable_more_warnings(env=env)
    env_etc.include_registry.append(
      env=env,
      paths=[
        env_etc.libtbx_include,
        env_etc.fable_include,
        env_etc.boost_adaptbx_include,
        env_etc.boost_include,
        env_etc.python_include])
    env.SharedLibrary(
      target="#lib/fable_ext",
      source=["ext.cpp"])
