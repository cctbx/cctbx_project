import libtbx
import libtbx.config
libtbx.env = libtbx.config.unpickle()
libtbx.env.set_os_environ_all_dist()
