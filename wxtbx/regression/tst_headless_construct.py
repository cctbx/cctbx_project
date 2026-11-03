from __future__ import absolute_import, division, print_function
import sys

from wxtbx.regression._warnings import install_wx_deprecation_filters

# Install at module load so subsequent imports of wx (and any wxtbx widget
# module reached during construction) operate under the strict filter.
install_wx_deprecation_filters()

def exercise():
  import wx
  try:
    app = wx.App(False)
  except (Exception, SystemExit) as e:
    # wx.App aborts on a headless box; on GTK this surfaces as a SystemExit
    # (a BaseException), which a plain "except Exception" would not catch.
    print("Skipped: cannot create wx.App ({})".format(e))
    print("OK")
    return
  frame = wx.Frame(None)
  from wxtbx.regression._headless_construct_cases import CASES
  ok = 0
  skipped = 0
  failed = []
  for label, factory, requires_data, skip_reason in CASES:
    if requires_data:
      print("SKIP {}: {}".format(label, skip_reason))
      skipped += 1
      continue
    install_wx_deprecation_filters()
    try:
      widget = factory(frame)
      try:
        widget.Destroy()
      except Exception:
        pass
      ok += 1
      print("OK   {}".format(label))
    except Exception as e:
      failed.append((label, e))
      print("FAIL {}: {}: {}".format(label, type(e).__name__, e))
  frame.Destroy()
  print("ok={} skipped={} failed={}".format(ok, skipped, len(failed)))
  if failed:
    sys.exit(1)
  print("OK")

if __name__ == "__main__":
  exercise()
