import libtbx.command_line.configure
import libtbx.env
build_options = libtbx.command_line.configure.build_options_t()
build_options.get_from_libtbx_env(env=libtbx.env.cache)
build_options.report()
