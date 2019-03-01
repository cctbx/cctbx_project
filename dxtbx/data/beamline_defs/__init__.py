from __future__ import absolute_import, division


def get_beamline_definition(detector_id, **kwargs):
    import unicodedata
    import string
    import types

    if not isinstance(detector_id, types.UnicodeType):
        detector_id = unicode(detector_id, "utf-8", "ignore")

    valid_chars = frozenset("_.%s%s" % (string.ascii_letters, string.digits))
    filename = unicodedata.normalize("NFKD", detector_id).encode("ASCII", "ignore")
    filename = "".join(c if c in valid_chars else "_" for c in filename)
    while "__" in filename:
        filename = filename.replace("__", "_")
        # http://stackoverflow.com/questions/1546226/a-simple-way-to-remove-multiple-spaces-in-a-string-in-python/15913564#15913564

    try:
        import importlib

        beamline = importlib.import_module("dxtbx.data.beamline_defs.%s" % filename)
        generator_object = beamline.get_definition(**kwargs)
    except ImportError:
        generator_object = Dummy()
    generator_object.set_detector_name(detector_id)
    generator_object.set_block_name(filename)
    return generator_object


class template(object):
    def CIF_block(self):
        """Interface function to generate a CIF block for this detector."""
        raise RuntimeError("This needs to be overridden.")

    def mmCIF_block(self):
        """Interface function to generate an mmCIF block for this detector."""
        raise RuntimeError("This needs to be overridden.")

    def _lookup(self, mmCIFsemantics):
        keys = {
            # Entries can either be
            #
            # a tuple
            #    ( CIF-string, mmCIF-string )
            #
            # or a string such as
            #    "_something?something_else"
            # where
            #   _something_something_else is CIF, and
            #   _something.something_else is mmCIF
            "df.detector": ("_diffrn_detector", "_diffrn.detector"),
            "df.rad.mono": "_diffrn_radiation?monochromator",
            "df.rad.source": "_diffrn_radiation?source",
            "df.rad.type": "_diffrn_radiation?type",
            "df.m.dev": "_diffrn?measurement_device",
            "df.m.dev_type": "_diffrn?measurement_device_type",
            "df.m.method": "_diffrn?measurement_method",
            "df.m.spec_supp": "_diffrn?measurement_specimen_support",
        }
        column = 1 if mmCIFsemantics else 0

        def __lookup(key):
            entry = keys[key]
            if isinstance(entry, basestring):
                return entry.replace("?", "." if mmCIFsemantics else "_")
            else:
                return entry[column]

        return __lookup

    def get_block_name(self):
        return self._block_name

    def set_block_name(self, block_name):
        self._block_name = block_name

    def get_detector_name(self):
        return self._detector_name

    def set_detector_name(self, detector_name):
        self._detector_name = detector_name

    def write_block_to_file(self, block, filename):
        import iotbx.cif.model

        c = iotbx.cif.model.cif()
        c[self.get_block_name()] = block
        with open(filename, "w") as fh:
            c.show(out=fh)

    def _date_to_epoch(self, year, month, day):
        import datetime

        return (
            datetime.datetime(year, month, day, 0, 0) - datetime.datetime(1970, 1, 1)
        ).total_seconds()

    def _generate_block(self):
        import iotbx.cif.model

        return iotbx.cif.model.block()

    def __str__(self):
        return "CIF block generator for %s (%s.py)" % (
            self._detector_name,
            self._block_name,
        )


class Dummy(template):
    def __init__(self):
        self.CIF_block = self._generate_block
        self.mmCIF_block = self._generate_block

    def __str__(self):
        return "Dummy CIF generator. No information for detector %s (%s.py)" % (
            self._detector_name,
            self._block_name,
        )
