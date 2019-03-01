from __future__ import absolute_import, division, print_function

from dxtbx.format.Format import Format
from dxtbx.format.FormatStill import FormatStill
from dxtbx.format.FormatMultiImage import Reader
from dxtbx.format.FormatMultiImageLazy import FormatMultiImageLazy
from libtbx.phil import parse

locator_str = """
  experiment = None
    .type = str
    .help = Experiment identifier, e.g. mfxo1916
  run = None
    .type = ints
    .help = Run number or a list of runs to process
  mode = idx
    .type = str
    .help = Mode for reading the xtc data (see LCLS documentation)
  data_source = None
    .type = str
    .help = Complete LCLS data source.  Overrides experiment and run.  Example: \
            exp=mfxo1916:run=20:smd \
            More info at https://confluence.slac.stanford.edu/display/PSDM/Manual#Manual-Datasetspecification
  detector_address = None
    .type = str
    .multiple = True
    .help = detector used for collecting the data at LCLS
  calib_dir = None
    .type = str
    .help = Specifiy path to custom calib directory if needed
"""
locator_scope = parse(locator_str)


class XtcReader(Reader):
    def nullify_format_instance(self):
        """ No-op for XTC streams. No issue with multiprocessing. """
        pass


class FormatXTC(FormatMultiImageLazy, FormatStill, Format):
    def __init__(self, image_file, **kwargs):
        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        FormatMultiImageLazy.__init__(self, **kwargs)
        FormatStill.__init__(self, image_file, **kwargs)
        Format.__init__(self, image_file, **kwargs)
        self.current_index = None
        self.current_event = None
        self._psana_runs = {}  ## empty container, to prevent breaking other formats
        if "locator_scope" in kwargs:
            self.params = FormatXTC.params_from_phil(
                master_phil=kwargs["locator_scope"], user_phil=image_file, strict=True
            )
        else:
            self.params = FormatXTC.params_from_phil(
                master_phil=locator_scope, user_phil=image_file, strict=True
            )
        assert self.params.mode == "idx", "idx mode should be used for analysis"
        self._initialized = True

    @staticmethod
    def understand(image_file):
        """ Extracts the datasource and detector_address from the image_file and then feeds it to PSANA
        If PSANA fails to read it, then input may not be an xtc/smd file. If success, then OK.
        If detector_address is not provided, a command line promp will try to get the address
        from the user """
        try:
            import psana
        except ImportError:
            return False
        try:
            params = FormatXTC.params_from_phil(locator_scope, image_file)
        except Exception:
            return False
        if params is None:
            return False

        ds = FormatXTC._get_datasource(image_file, params)
        return True

    @staticmethod
    def params_from_phil(master_phil, user_phil, strict=False):
        """ Read the locator file """
        try:
            user_input = parse(file_name=user_phil)
            working_phil, unused = master_phil.fetch(
                sources=[user_input], track_unused_definitions=True
            )
            unused_args = ["%s=%s" % (u.path, u.object.words[0].value) for u in unused]
            if len(unused_args) > 0 and strict:
                for unused_arg in unused_args:
                    print(unused_arg)
                print(
                    "Incorrect of unused parameter in locator file. Please check and retry"
                )
                return None
            params = working_phil.extract()
            return params
        except Exception:
            return None

    @classmethod
    def get_reader(Class):
        """
    Return a reader class

    """
        obj = XtcReader
        # Note, need to set this on the parent class since it's a scoped global variable
        Reader._format_class_ = Class
        return obj

    def populate_events(self):
        """ Read the timestamps from the XTC stream.  Assumes the psana idx mode of reading data.
        Handles multiple LCLS runs by concatenating the timestamps from multiple runs together
        in a single list and creating a mapping. """
        if hasattr(self, "times") and len(self.times) > 0:
            return

        if not self._psana_runs:
            self._psana_runs = self._get_psana_runs(self._ds)

        self.times = []
        self.run_mapping = {}
        for run in self._psana_runs.values():
            times = run.times()
            self.run_mapping[run.run()] = (
                len(self.times),
                len(self.times) + len(times),
                run,
            )
            self.times.extend(times)

    def get_run_from_index(self, index=None):
        """ Look up the run number given an index """
        if index is None:
            index = 0
        for run_number in self.run_mapping:
            start, stop, run = self.run_mapping[run_number]
            if index >= start and index < stop:
                return run
        raise IndexError("Index is not within bounds")

    def _get_event(self, index):
        """ Retrieve a psana event given and index. This is the slow step for reading XTC streams,
        so implement a cache for the last read event. """
        if index == self.current_index:
            return self.current_event
        else:
            self.current_index = index
            self.current_event = self.get_run_from_index(index).event(self.times[index])
            return self.current_event

    @staticmethod
    def _get_datasource(image_file, params):
        """ Construct a psana data source object given the locator parameters """
        import psana

        if params.calib_dir is not None:
            psana.setOption("psana.calib-dir", params.calib_dir)
        if params.data_source is None:
            if (
                params.experiment is None
                or params.run is None
                or params.mode is None
                or len(params.run) == 0
            ):
                return False
            img = "exp=%s:run=%s:%s" % (
                params.experiment,
                ",".join(["%d" % r for r in params.run]),
                params.mode,
            )
        else:
            img = params.data_source
        return psana.DataSource(img)

    @staticmethod
    def _get_psana_runs(datasource):
        """
    Extracts the runs,
    These can only be extracted once,
    only call this method after datasource is set
    """
        # this is key,value = run_integer, psana.Run, e.g. {62: <psana.Run(@0x7fbd0e23c990)>}
        psana_runs = {r.run(): r for r in datasource.runs()}
        return psana_runs

    def get_psana_timestamp(self, index):
        """ Get the cctbx.xfel style event timestamp given an index """
        from xfel.cxi.cspad_ana import cspad_tbx
        import psana

        evt = self._get_event(index)
        time = evt.get(psana.EventId).time()
        fid = evt.get(psana.EventId).fiducials()

        sec = time[0]
        nsec = time[1]

        return cspad_tbx.evt_timestamp((sec, nsec / 1e6))


if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        print(FormatXTC.understand(arg))
