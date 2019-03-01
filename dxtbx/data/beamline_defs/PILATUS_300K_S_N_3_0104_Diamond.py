from __future__ import absolute_import, division
import dxtbx.data.beamline_defs


class get_definition(dxtbx.data.beamline_defs.template):
    def __init__(self, timestamp=None, **kwargs):
        self._timestamp = timestamp

    def CIF_block(self):
        """Interface function to generate a CIF block for this detector."""
        return self._identify_time()(mmcif=False)

    def mmCIF_block(self):
        """Interface function to generate an mmCIF block for this detector."""
        return self._identify_time()(mmcif=True)

    def _identify_time(self):
        """Detector has always been on I19."""
        return self._at_I19

    def _base(self, mmcif=False):
        """Generates
        1. a CIF/mmCIF block that contains information that
        is always true about the detector.
        2. a lookup function for CIF/mmCIF strings."""
        # prepare string lookup table
        l = self._lookup(mmcif)

        import iotbx.cif.model

        b = iotbx.cif.model.block()
        b[l("df.detector")] = "Photon counting pixel array"
        b[l("df.rad.type")] = "Synchrotron"

        return b, l

    def _at_I19(self, mmcif=False):
        b, l = self._base(mmcif)

        #   b[l('df.m.dev')]       = 'Fixed \\c 3-circle diffractometer'
        #   b[l('df.m.dev_type')]  = 'Fluid Film Devices'
        b[l("df.m.method")] = "shutterless scans"
        #   b[l('df.m.spec_supp')] = 'MiTeGen MicroMount'
        b[l("df.rad.source")] = "Diamond Light Source Beamline I19-2"
        b[l("df.rad.mono")] = "Silicon 111"

        return b
