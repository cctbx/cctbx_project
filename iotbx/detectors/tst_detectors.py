from __future__ import absolute_import, division, print_function
def tst_detectors():
    import libtbx.load_env
    import os

    image_directory = libtbx.env.find_in_repositories(
        relative_path = 'bpcx_regression/use_case_xds_method')
    if not image_directory:
        print('bpcx_regression not configured, skipping test')
        return

    import iotbx.detectors
    f = iotbx.detectors.ImageFactory(
        os.path.join(image_directory, 'thau2_O0_K0_P0_1_0001.cbf'))
    f.read()
    length = len(f.get_raw_data())
    fi = f.get_flex_image()
    assert(length == fi.size1() * fi.size2())
    print('OK')

if __name__ == '__main__':
    tst_detectors()
