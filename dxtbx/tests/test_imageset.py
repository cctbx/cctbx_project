from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
import glob
import os

import dxtbx.tests.imagelist
import pytest


@pytest.mark.parametrize(
    "image",
    dxtbx.tests.imagelist.smv_images
    + dxtbx.tests.imagelist.tiff_images
    + dxtbx.tests.imagelist.cbf_multitile_images
    + dxtbx.tests.imagelist.cbf_images,
    ids=(
        dxtbx.tests.imagelist.smv_image_ids
        + dxtbx.tests.imagelist.tiff_image_ids
        + dxtbx.tests.imagelist.cbf_multitile_image_ids
        + dxtbx.tests.imagelist.cbf_image_ids
    ),
)
def test_format(dials_regression, image):
    from dxtbx.format.Registry import Registry

    db_fail_count = 0
    print(image)
    image = os.path.join(dials_regression, *(image.split("/")))
    format_class = Registry.find(image)
    reader = format_class.get_reader()([image])
    masker = format_class.get_masker()([image])

    N = len(reader)

    for i in range(N):
        data = reader.read(i)
        mask = masker.get(i)

    iset = format_class.get_imageset([image])


def test_image_tile():
    from dxtbx.format.image import ImageTileInt
    from scitbx.array_family import flex

    data = flex.int(flex.grid(10, 10))
    name = "TileName"

    tile = ImageTileInt(data, name)

    assert tile.data().all_eq(data)
    assert tile.name() == name
    assert tile.empty() is False


def test_image():
    import dxtbx.format.image
    from scitbx.array_family import flex

    data = flex.int(flex.grid(10, 10))
    name = "TileName0"
    tile0 = dxtbx.format.image.ImageTileInt(data, name)
    image = dxtbx.format.image.ImageInt(tile0)
    for i in range(1, 4):
        data = flex.int(flex.grid(10, 10))
        name = "TileName%d" % i
        tile = dxtbx.format.image.ImageTileInt(data, name)
        image.append(tile)

    assert image.n_tiles() == 4
    for i in range(image.n_tiles()):
        tile = image.tile(i)
        assert tile.name() == "TileName%d" % i


def test_image_buffer():
    import dxtbx.format.image
    from scitbx.array_family import flex

    data = flex.int(flex.grid(10, 10))
    name = "TileName0"
    tile0 = dxtbx.format.image.ImageTileInt(data, name)
    image = dxtbx.format.image.ImageInt(tile0)

    b = dxtbx.format.image.ImageBuffer(image)
    assert b.is_int() is True
    assert b.is_double() is False
    assert b.is_empty() is False


def test_external_lookup():
    import dxtbx.format.image
    from dxtbx.imageset import ExternalLookup
    from scitbx.array_family import flex

    mask = flex.bool(flex.grid(10, 10), True)
    gain = flex.double(flex.grid(10, 10), 1)
    pedestal = flex.double(flex.grid(10, 10), 2)

    lookup = ExternalLookup()
    lookup.mask.data = dxtbx.format.image.ImageBool(
        dxtbx.format.image.ImageTileBool(mask)
    )
    lookup.gain.data = dxtbx.format.image.ImageDouble(
        dxtbx.format.image.ImageTileDouble(gain)
    )
    lookup.pedestal.data = dxtbx.format.image.ImageDouble(
        dxtbx.format.image.ImageTileDouble(pedestal)
    )

    mask2 = lookup.mask.data.tile(0).data()
    gain2 = lookup.gain.data.tile(0).data()
    pedestal2 = lookup.pedestal.data.tile(0).data()

    assert mask2.all_eq(mask)
    assert gain2.all_eq(gain)
    assert pedestal2.all_eq(pedestal)


def test_imagesetdata(dials_regression):
    from dxtbx.imageset import ImageSetData
    from dxtbx.format.image import (
        ImageBool,
        ImageTileBool,
        ImageDouble,
        ImageTileDouble,
    )
    from scitbx.array_family import flex
    from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus as FormatClass

    filenames = sorted(
        glob.glob(os.path.join(dials_regression, "centroid_test_data", "*.cbf"))
    )

    ReaderClass = FormatClass.get_reader()
    MaskerClass = FormatClass.get_masker()

    reader = ReaderClass(filenames)
    masker = MaskerClass(filenames)

    handle = ImageSetData(reader, masker)

    data = handle.get_data(0).as_int().tile(0).data()
    mask = handle.get_mask(0).tile(0).data()

    assert handle.has_single_file_reader() is False

    path = handle.get_path(0)
    assert path == filenames[0]

    master_path = handle.get_master_path()
    assert master_path == ""

    identifier = handle.get_image_identifier(0)
    assert identifier == filenames[0]

    beam = FormatClass(filenames[0]).get_beam()
    detector = FormatClass(filenames[0]).get_detector()
    goniometer = FormatClass(filenames[0]).get_goniometer()
    scan = FormatClass(filenames[0]).get_scan()

    handle.set_beam(beam, 0)
    handle.set_detector(detector, 0)
    handle.set_goniometer(goniometer, 0)
    handle.set_scan(scan, 0)

    beam2 = handle.get_beam(0)
    detector2 = handle.get_detector(0)
    goniometer2 = handle.get_goniometer(0)
    scan2 = handle.get_scan(0)

    assert beam2 == beam
    assert detector2 == detector
    assert goniometer2 == goniometer
    assert scan2 == scan

    mask = flex.bool(flex.grid(10, 10), True)
    gain = flex.double(flex.grid(10, 10), 1)
    pedestal = flex.double(flex.grid(10, 10), 2)

    handle.external_lookup.mask.data = ImageBool(ImageTileBool(mask))
    handle.external_lookup.gain.data = ImageDouble(ImageTileDouble(gain))
    handle.external_lookup.pedestal.data = ImageDouble(ImageTileDouble(pedestal))

    mask2 = handle.external_lookup.mask.data.tile(0).data()
    gain2 = handle.external_lookup.gain.data.tile(0).data()
    pedestal2 = handle.external_lookup.pedestal.data.tile(0).data()

    assert mask2.all_eq(mask)
    assert gain2.all_eq(gain)
    assert pedestal2.all_eq(pedestal)


@pytest.fixture(scope="session")
def centroid_files(dials_regression):
    return [
        os.path.join(dials_regression, "centroid_test_data", "centroid_%04d.cbf" % i)
        for i in range(1, 10)
    ]


@pytest.fixture
def centroid_files_and_imageset(centroid_files):
    from dxtbx.format.Registry import Registry

    # Create the format class
    format_class = Registry.find(centroid_files[0])

    # Create the reader
    imageset = format_class.get_imageset(centroid_files, as_imageset=True)

    return centroid_files, imageset


class TestImageSet(object):
    def test_imageset(self, centroid_files_and_imageset):
        filenames, imageset = centroid_files_and_imageset

        # Run a load of tests
        self.tst_get_item(imageset)
        assert len(imageset) == len(filenames)
        self.tst_iter(imageset)
        self.tst_paths(imageset, filenames)
        self.tst_get_detectorbase(imageset, range(len(filenames)), 9)
        self.tst_get_models(imageset, range(len(filenames)), 9)

    def tst_get_item(self, imageset):
        image = imageset[0]
        with pytest.raises(Exception):
            image = imageset[9]

        imageset2 = imageset[3:7]
        image = imageset2[0]
        with pytest.raises(Exception):
            image = imageset2[5]

        assert len(imageset2) == 4
        self.tst_get_detectorbase(imageset2, range(0, 4), 5)
        self.tst_get_models(imageset2, range(0, 4), 5)
        self.tst_paths(imageset2, imageset.paths()[3:7])
        self.tst_iter(imageset2)

        imageset2 = imageset[3:5]
        image = imageset2[0]
        with pytest.raises(Exception):
            image = imageset2[2]

        assert len(imageset2) == 2
        self.tst_get_detectorbase(imageset2, range(0, 2), 2)
        self.tst_get_models(imageset2, range(0, 2), 2)
        self.tst_paths(imageset2, imageset.paths()[3:5])
        self.tst_iter(imageset2)

    @staticmethod
    def tst_iter(imageset):
        for image in imageset:
            pass

    @staticmethod
    def tst_paths(imageset, filenames1):
        filenames2 = imageset.paths()
        for f1, f2 in zip(filenames1, filenames2):
            assert f1 == f2

    @staticmethod
    def tst_get_detectorbase(imageset, indices, outside_index):
        for i in indices:
            imageset.get_detectorbase(i)

        with pytest.raises(Exception):
            imageset.get_detectorbase(outside_index)

    def tst_get_models(self, imageset, indices, outside_index):
        for i in indices:
            self.tst_get_models_index(imageset, i)

        with pytest.raises(Exception):
            self.tst_get_models_index(imageset, outside_index)

    @staticmethod
    def tst_get_models_index(imageset, index=None):
        imageset.get_detector(index)
        imageset.get_beam(index)

    @staticmethod
    def tst_set_models(imageset):
        from dxtbx.model import Beam, Detector, Panel

        # Create some other models
        beam = Beam((1, 0, 0), 0.5)
        detector = Detector(
            Panel(
                "UNKNOWN",
                "Panel",
                (1, 0, 0),
                (0, 1, 0),
                (0, 0, 1),
                (0.1, 0.1),
                (1000, 1000),
                (0, 1),
            )
        )

        # Override sweep models
        imageset.set_beam(beam)
        imageset.set_detector(detector)

        # Ensure this doens't interfere with reading
        for i in imageset:
            pass

        # Get the models back and check they're ok
        beam2 = imageset.get_beam()
        detector2 = imageset.get_detector()
        assert beam2 == beam
        assert detector2 == detector

        # Get the models from an index back and check they're not the same
        beam2 = imageset.get_beam(0)
        detector2 = imageset.get_detector(0)
        assert beam2 != beam
        assert detector2 != detector


class TestImageSweep(object):
    def test(self, centroid_files):
        from dxtbx.imageset import ImageSweep
        from dxtbx.format.Registry import Registry

        # Create the format class
        format_class = Registry.find(centroid_files[0])

        # Create the sweep
        sweep = format_class.get_imageset(centroid_files)

        # Run a load of tests
        self.tst_get_item(sweep)
        assert len(sweep) == len(centroid_files)
        self.tst_iter(sweep)
        self.tst_paths(sweep, centroid_files)
        self.tst_get_detectorbase(sweep, range(len(centroid_files)), 9)
        self.tst_get_models(sweep, range(len(centroid_files)), 9)
        self.tst_get_array_range(sweep, (0, 9))
        self.tst_set_models(sweep)

    def tst_get_item(self, sweep):
        image = sweep[0]
        with pytest.raises(Exception):
            image = sweep[9]

        sweep2 = sweep[3:7]
        image = sweep2[0]
        with pytest.raises(Exception):
            image = sweep2[5]

        assert len(sweep2) == 4
        self.tst_get_detectorbase(sweep2, range(0, 4), 5)
        self.tst_get_models(sweep2, range(0, 4), 5)
        self.tst_paths(sweep2, sweep.paths()[3:7])
        self.tst_iter(sweep2)
        self.tst_get_array_range(sweep2, (3, 7))

        with pytest.raises(IndexError):
            sweep2 = sweep[3:7:2]

    @staticmethod
    def tst_iter(sweep):
        for image in sweep:
            pass

    @staticmethod
    def tst_paths(sweep, filenames1):
        filenames2 = sweep.paths()
        for f1, f2 in zip(filenames1, filenames2):
            assert f1 == f2

    @staticmethod
    def tst_get_detectorbase(sweep, indices, outside_index):
        for i in indices:
            sweep.get_detectorbase(i)

        with pytest.raises(Exception):
            sweep.get_detectorbase(outside_index)

    def tst_get_models(self, sweep, indices, outside_index):
        self.tst_get_models_index(sweep)
        for i in indices:
            self.tst_get_models_index(sweep, i)

    @staticmethod
    def tst_get_models_index(sweep, index=None):
        if index is not None:
            sweep.get_detector(index)
            sweep.get_beam(index)
            sweep.get_goniometer(index)
            sweep.get_scan(index)
        else:
            sweep.get_detector()
            sweep.get_beam()
            sweep.get_goniometer()
            sweep.get_scan()

        # Ensure state at zero
        sweep[0]
        scan1 = sweep.get_scan()
        # Put sweep to end
        sweep[len(sweep) - 1]
        scan2 = sweep.get_scan()
        assert scan1 == scan2

    @staticmethod
    def tst_get_array_range(sweep, array_range):
        assert sweep.get_array_range() == array_range

    @staticmethod
    def tst_set_models(sweep):
        from dxtbx.model import Beam, Detector, Panel

        # Get some models
        beam = sweep.get_beam()
        gonio = sweep.get_goniometer()
        detector = sweep.get_detector()

        # Modify the geometry
        assert len(detector) == 1
        beam.set_direction((1, 0, 0))
        gonio.set_rotation_axis((0, 1, 0))
        detector[0].set_local_frame((1, 0, 0), (0, 1, 0), (0, 0, 1))

        # Override sweep models
        sweep.set_beam(beam)
        sweep.set_goniometer(gonio)
        sweep.set_detector(detector)

        # Ensure this doens't interfere with reading
        for i in sweep:
            pass

        # Get the models back and check they're ok
        beam2 = sweep.get_beam()
        gonio2 = sweep.get_goniometer()
        detector2 = sweep.get_detector()
        assert beam2 == beam
        assert gonio2 == gonio
        assert detector2 == detector

        # Get the models from an index back and check they're the same
        beam2 = sweep.get_beam(0)
        gonio2 = sweep.get_goniometer(0)
        detector2 = sweep.get_detector(0)
        assert beam2 == beam
        assert gonio2 == gonio
        assert detector2 == detector

        # Get a sub sweep
        sub_sweep = sweep[3:7]

        # Get the models back and check they're ok
        beam2 = sub_sweep.get_beam()
        gonio2 = sub_sweep.get_goniometer()
        detector2 = sub_sweep.get_detector()
        assert beam2 == beam
        assert gonio2 == gonio
        assert detector2 == detector

        # Get the models from an index back and check they're not the same
        beam2 = sub_sweep.get_beam(0)
        gonio2 = sub_sweep.get_goniometer(0)
        detector2 = sub_sweep.get_detector(0)
        assert beam2 == beam
        assert gonio2 == gonio
        assert detector2 == detector


@pytest.mark.xfail(
    reason="This test is broken since the master h5 file is inconsistent with NeXus format after commit dbb0bf7"
)
def test_nexus_file(dials_regression):
    filename = os.path.join(
        dials_regression,
        "image_examples",
        "LCLS_cspad_nexus",
        "cxi78513_bslz4_r0014_subset4_master.h5",
    )

    from dxtbx.format.Registry import Registry

    format_class = Registry.find(filename)

    iset = format_class.get_imageset([filename])

    assert len(iset) == 2
    for i in range(len(iset)):
        data = iset.get_raw_data(i)
        mask = iset.get_mask(i)
        b = iset.get_beam(i)
        d = iset.get_detector(i)
        g = iset.get_goniometer(i)
        s = iset.get_scan(i)

    iset = format_class.get_imageset([filename], single_file_indices=[1])
    assert len(iset) == 1

    for i in range(len(iset)):
        data = iset.get_raw_data(i)
        mask = iset.get_mask(i)
        b = iset.get_beam(i)
        d = iset.get_detector(i)
        g = iset.get_goniometer(i)
        s = iset.get_scan(i)


@pytest.mark.parametrize("lazy", (True, False))
def test_SACLA_MPCCD_Cheetah_File(dials_regression, lazy):
    filename = os.path.join(
        dials_regression,
        "image_examples",
        "SACLA_MPCCD_Cheetah",
        "run266702-0-subset.h5",
    )

    from dxtbx.format.Registry import Registry

    format_class = Registry.find(filename)

    iset = format_class.get_imageset([filename], lazy=lazy)

    assert len(iset) == 4
    for i in range(len(iset)):
        assert iset.get_raw_data(i)
        #   assert iset.get_mask(i)
        assert iset.get_beam(i)
        assert iset.get_detector(i)
        assert iset.get_goniometer(i) is None
        assert iset.get_scan(i) is None

    iset = format_class.get_imageset([filename], single_file_indices=[1], lazy=lazy)
    assert len(iset) == 1

    for i in range(len(iset)):
        assert iset.get_raw_data(i)
        #   assert iset.get_mask(i)
        assert iset.get_beam(i)
        assert iset.get_detector(i)
        assert iset.get_goniometer(i) is None
        assert iset.get_scan(i) is None


def test_imagesetfactory(centroid_files, dials_regression):
    from dxtbx.imageset import ImageSetFactory, ImageSweep

    filenames = centroid_files

    sweep = ImageSetFactory.new(filenames)

    assert isinstance(sweep[0], ImageSweep)

    template = os.path.join(dials_regression, "centroid_test_data", "centroid_####.cbf")
    image_range = (3, 6)

    sweep = ImageSetFactory.from_template(template, image_range)

    assert isinstance(sweep[0], ImageSweep)
    assert len(sweep[0]) == 4
    assert sweep[0].paths()[0].endswith("3.cbf")
    assert sweep[0].paths()[-1].endswith("6.cbf")

    imageset = ImageSetFactory.make_imageset(filenames)
    assert len(imageset) == 9

    imageset = ImageSetFactory.make_imageset(filenames, check_format=False)
    assert len(imageset) == 9

    sweep = ImageSetFactory.make_sweep(template, list(range(1, 9 + 1)))
    assert len(sweep) == 9

    sweep = ImageSetFactory.make_sweep(template, list(range(3, 6 + 1)))
    assert len(sweep) == 4


def test_pickle_imageset(centroid_files):
    from dxtbx.imageset import ImageSetFactory, ImageSweep

    sweep = ImageSetFactory.new(centroid_files)[0]

    # Read the 5th image
    image = sweep[4]

    # Pickle, then unpickle
    sweep2 = pickle.loads(pickle.dumps(sweep))

    assert sweep.get_template() == sweep2.get_template()
    assert sweep.get_array_range() == sweep2.get_array_range()
    assert sweep.get_beam() == sweep2.get_beam()
    assert sweep.get_goniometer() == sweep2.get_goniometer()
    assert sweep.get_scan() == sweep2.get_scan()
    assert sweep.paths() == sweep2.paths()
    assert sweep == sweep2

    # Check auxiliary methods after pickling
    sweep3 = sweep2[0:4]
    sweep4 = sweep3[0:2]
    sweep4.get_detectorbase(0)
    sweep4[0]
