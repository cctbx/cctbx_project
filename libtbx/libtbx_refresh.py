# copy script to enable checkout even if libtbx sources are removed already
source = self.env.under_dist("libtbx", "libtbx/sourceforge_checkout.py")
copy = self.env.under_build("sourceforge_checkout.py")
open(copy, "w").write(open(source).read())
#
self.env.write_dispatcher_in_bin(
  source_file=copy,
  target_file="libtbx.sourceforge_checkout")
