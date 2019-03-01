from __future__ import absolute_import, division, print_function

import os


def test_load_path():
    from dxtbx.serialize.filename import load_path

    os.environ["HELLO_WORLD"] = "EXPANDED"
    new_path = os.path.join("~", "$HELLO_WORLD", "path")
    path = load_path(new_path)
    assert path == os.path.join(os.path.expanduser("~"), "EXPANDED", "path")
    new_path = os.path.join("$HELLO_WORLD", "path")
    path = load_path(new_path)
    assert path == os.path.abspath(os.path.join("EXPANDED", "path"))
