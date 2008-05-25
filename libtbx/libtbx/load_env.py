import libtbx
import libtbx.config
import os
libtbx.env = libtbx.config.unpickle()
libtbx.env.set_os_environ_all_dist()
libtbx.env.dispatcher_name = os.environ.get("LIBTBX_DISPATCHER_NAME")
libtbx.env.full_testing = os.environ.get("LIBTBX_FULL_TESTING") is not None
