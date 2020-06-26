from read_geom import read_geom

from libtbx.phil import parse
from libtbx.utils import Sorry

import h5py
import os
import subprocess
import sys
from pathlib import Path
from lxml import etree, objectify
import logging

LOGGING_FORMAT = "%(levelname)s: %(message)s"
logging.basicConfig(level=logging.INFO, format=LOGGING_FORMAT)
logger = logging.getLogger()

p = Path(__file__)
path_to_defs = p.parent / 'definitions/applications/NXmx.nxdl.xml'


root = objectify.parse(str(path_to_defs))

for elem in root.getiterator():
    try:
        elem.tag = etree.QName(elem).localname
    except ValueError:
        print(f"Error with {elem.tag}")

phil_scope = parse("""
  cxi_file = None
    .type = str
    .help = cheetah file used to read in image data(.cxi).
  geom_file = None
    .type = str
    .help = geometry file to be read in for AGIPD detector (.geom).
  detector_distance = None
    .type = float
    .help = AGIPD Detector distance
  wavelength = None
    .type = float
    .help = AGIPD wavelength override
  mode = vds cxi
    .type = choice
    .optional = False
    .help = Input data file format. VDS: virtual data set. CXI: \
            Cheetah file format.
  
  nexus_details {
    instrument_name = SwissFEL ARAMIS BEAMLINE ESB
      .type = str
      .help = Name of instrument
    instrument_short_name = ESB
      .type = str
      .help = short name for instrument, perhaps the acronym
    source_name = SwissFEL ARAMIS
      .type = str
      .help = Name of the neutron or x-ray storage ring/facility
    source_short_name = SwissFEL ARAMIS
      .type = str
      .help = short name for source, perhaps the acronym
    start_time = None
      .type = str
      .help = ISO 8601 time/date of the first data point collected in UTC, \
              using the Z suffix to avoid confusion with local time
    end_time = None
      .type = str
      .help = ISO 8601 time/date of the last data point collected in UTC, \
              using the Z suffix to avoid confusion with local time. \
              This field should only be filled when the value is accurately \
              observed. If the data collection aborts or otherwise prevents \
              accurate recording of the end_time, this field should be omitted
    end_time_estimated = None
      .type = str
      .help = ISO 8601 time/date of the last data point collected in UTC, \
              using the Z suffix to avoid confusion with local time. \
              This field may be filled with a value estimated before an \
              observed value is avilable.
    sample_name = None
      .type = str
      .help = Descriptive name of sample
    total_flux = None
      .type = float
      .help = flux incident on beam plane in photons per second
  }

""")


def get_root_path(tree_elem: etree) -> str:
    l = list(x.attrib['type'].replace('NX', '') for x in tree_elem.iterancestors())
    return '/'.join(l[:-1][::-1])


def get_git_revision_hash() -> str:
    definitions_path = os.path.join(os.path.dirname(sys.argv[0]), 'definitions')
    current_dir = os.getcwd()
    os.chdir(definitions_path)
    definitions_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip()
    os.chdir(current_dir)
    return definitions_hash.decode()


class Agipd2nexus:

    def __init__(self, args):
        self.group_rules = {}
        self.field_rules = {}

    def params_from_phil(self, args):
        user_phil = []
        for arg in args:
            if os.path.isfile(arg):
                user_phil.append(parse(file_name=arg))
            else:
                try:
                    user_phil.append(parse(arg))
                except Exception as e:
                    raise Sorry("Unrecognized argument: %s" % arg)
        self.params = phil_scope.fetch(sources=user_phil).extract()

    def create_group(self, h5file: h5py.File, lx_elem: objectify.ObjectifiedElement, lx_parent: str) -> None:
        """
        Create a new h5.Group attached to the `parent_elem`
        """
        NXname = lx_elem.attrib['type']
        if NXname.startswith('NX'):
            name = NXname.replace('NX', '')
        else:
            logger.warning(f'name {NXname} is not a correct NX-name')
            return
        logger.debug(f"GROUP: {lx_parent} --> {NXname}")
        root_path = get_root_path(lx_elem)
        full_path = '/'.join([root_path, name])
        parent = h5file[root_path]
        if full_path in self.group_rules:
            new_group = parent.create_group(self.group_rules[full_path])
        else:
            new_group = parent.create_group(name)
        new_group.attrs['NX_class'] = NXname

    def create_field(self, h5file: h5py.File, lx_elem: objectify.ObjectifiedElement, lx_parent: str) -> None:
        NXname = lx_elem.attrib['name']
        name = NXname.replace('NX', '')
        root_path = get_root_path(lx_elem)
        full_path = '/'.join([root_path, name])
        # parent = h5file[path]
        logger.debug(f"FIELD: {lx_parent} --> {NXname}")
        if full_path in self.field_rules:
            h5file[full_path] = self.field_rules[full_path]
        else:
            h5file[full_path] = 'XXX'

    def create_nexus_master_file(self):
        self.output_file_name = os.path.splitext(self.params.cxi_file)[0] + '_master.h5'
        f = h5py.File(self.output_file_name, 'w')
        logger.info(f"file {self.output_file_name} was created")

        for elem in root.getiterator(('group', 'field')):
            parent = elem.getparent().attrib['type']
            if parent == 'group':
                # create root element of the file
                entry = f.create_group('entry')
                entry.attrs['NX_class'] = 'NXentry'
                continue

            if elem.tag == 'group':
                if 'minOccurs' in elem.keys() and elem.attrib['minOccurs'] == '0':
                    # print(f">>> GROUP {elem.attrib['type']} is optional")
                    continue
                self.create_group(f, elem, parent)
            elif elem.tag == 'field':
                if 'minOccurs' in elem.keys() and elem.attrib['minOccurs'] == '0':
                    # optional field, skip
                    # print(f">>> FIELD {elem.attrib['name']} is optional")
                    continue
                # field[parent] += [elem.attrib['name']]
                self.create_field(f, elem, parent)


class Ruleset(Agipd2nexus):

    def __init__(self, args):
        # constructor of the parent class first
        Agipd2nexus.__init__(self, args)

        self.params_from_phil(args)
        if self.params.detector_distance == None:
            # FIXME: should be taken from the `geom` file
            self.params.detector_distance = 177.0  # Set detector distance arbitrarily if nothing is provided
        self.hierarchy = read_geom(self.params.geom_file)
        self.n_quads = 4
        self.n_modules = 4
        self.n_asics = 8

        self.field_rules = {
            'entry/definition': f"NXmx:{get_git_revision_hash()}",      # TODO: _create_scalar
            'entry/start_time': str(self.params.nexus_details.start_time),
            'entry/end_time': str(self.params.nexus_details.start_time),
            'entry/end_time_estimated': str(self.params.nexus_details.start_time),
        }
        self.group_rules = {
            'NXdetector': ['ELE_D0'],
            'NXdetector_group': ['AGIPD'],

        }


if __name__ == '__main__':
  nexus_helper = Ruleset(sys.argv[1:])
  nexus_helper.create_nexus_master_file()
  # os.system(f'h5glance {nexus_helper.output_file_name} --attrs')