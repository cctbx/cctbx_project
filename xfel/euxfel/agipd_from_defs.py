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
import numpy as np
from datetime import datetime as dt
from typing import Union, Any
from enum import Enum, auto
from collections import Counter

LOGGING_FORMAT = "%(levelname)s: %(message)s"
logging.basicConfig(level=logging.INFO, format=LOGGING_FORMAT)
logger = logging.getLogger()

p = Path(__file__)
path_to_defs = p.parent / 'definitions/applications/NXmx.nxdl.xml'
path_to_cnxvalidate = "~/github/cnxvalidate/build/nxvalidate"
path_to_phil = p.parent / "AGIPD.phil"

root = objectify.parse(str(path_to_defs))

for elem in root.getiterator():
    try:
        elem.tag = etree.QName(elem).localname
    except ValueError:
        print(f"Error with {elem.tag}")

with open(path_to_phil) as phil:
    phil_scope = parse(phil.read())

class NxType(Enum):
    """ Enumeration for elements withing NeXus schema """
    group = auto()
    field = auto()
    attribute = auto()


class NexusElement:
    def __init__(self, full_path: str, value: Any, nxtype: NxType, dtype: str, parent: str = '',
                 attrs: dict = None) -> None:
        self.value = value
        self.dtype = dtype
        self.type = nxtype
        self.parent = parent
        self.full_path = full_path
        self.attrs = attrs

    def push(self, h5file):
        """ Write an element to the file """
        if self.full_path:
            h5file[self.full_path] = self.value
            if self.attrs:
                for k, v in self.attrs.items():
                    h5file[self.full_path].attrs[k] = v
        else:
            logger.error(f"Cannot push {self.name}")


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
        self.additional_elements = {}
        self.params_from_phil(args)
        self.output_file_name = os.path.splitext(self.params.cxi_file)[0] + '_master.h5'
        self.stat = Counter()

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

    @staticmethod
    def get_root_path(tree_elem: etree) -> str:
        """ return a NeXus path (e.g. `entry/instrument/detector_group`) from an lxml tree element """
        l = list(x.attrib['type'].replace('NX', '') for x in tree_elem.iterancestors())
        return '/'.join(l[:-1][::-1])

    def create_group(self, h5file: h5py.File, lx_elem: objectify.ObjectifiedElement, lx_parent: str) -> None:
        """
        Create a new h5.Group attached to the `parent_elem`

        The group name might be provided in the `group_rules`, otherwise it is created on the fly
        either by reducing `NXblah` to `blah` or from the `name` attribute

        Group rules might contain additional actions like creation specific content
        """
        NXname = lx_elem.attrib['type']     # TODO: the name might be contained in the `name` attribute
        if NXname.startswith('NX'):
            name = NXname.replace('NX', '')
        else:
            logger.warning(f'name {NXname} is not a correct NX-name')
            return
        logger.debug(f"GROUP: {lx_parent} --> {NXname}")
        root_path = self.get_root_path(lx_elem)
        full_path = '/'.join([root_path, name])
        parent = h5file[root_path]
        if full_path in self.group_rules:
            for n in self.group_rules[full_path]['names']:
                new_group = parent.create_group(n)
                new_group.attrs['NX_class'] = NXname
                logger.debug(f"group {n} was created")
                self.stat['groups from rules'] += 1
        else:
            new_group = parent.create_group(name)
            new_group.attrs['NX_class'] = NXname
            logger.debug(f"group {name} was created")
            self.stat['groups from defs'] += 1

    def create_field(self, h5file: h5py.File,
                     lx_elem: Union[objectify.ObjectifiedElement, dict]) -> None:
        NXname = lx_elem.attrib['name']
        name = NXname.replace('NX', '')
        if isinstance(lx_elem, objectify.ObjectifiedElement):
            root_path = self.get_root_path(lx_elem)
            full_path = '/'.join([root_path, name])

        if full_path in self.field_rules:
            logger.debug(f"Add {full_path} from a rule")
            if isinstance(self.field_rules[full_path], dict):
                vector = self.field_rules[full_path]
                h5file[full_path] = np.array(vector.pop('value'), dtype='f')
                for key, attribute in vector.items():
                    h5file[full_path].attrs[key] = attribute
            elif isinstance(self.field_rules[full_path], NexusElement):
                self.field_rules[full_path].push(h5file)
            else:
                h5file[full_path] = self.field_rules[full_path]
            logger.debug(f"field {full_path} was added")
            self.stat['fields from rules'] += 1

        else:
            logger.debug(f"Add {full_path} from definition")
            field = NexusElement(full_path=full_path, value='XXX', nxtype=NxType.field, dtype='f')
            field.push(h5file)
            self.stat['fields from defs'] += 1

    def create_attribute(self, h5file: h5py.File, elem):
        NXname = elem.attrib['name']
        name = NXname.replace('NX', '')
        if isinstance(elem, objectify.ObjectifiedElement):
            root_path = self.get_root_path(elem)
            full_path = '/'.join([root_path, name])

        h5file[full_path] = elem.attrib['value']
        self.stat['attr'] += 1

    def create_nexus_master_file(self):
        f = h5py.File(self.output_file_name, 'w')
        logger.info(f"file {self.output_file_name} was created")

        for k, v in self.global_attrs.items():
            f.attrs[k] = v

        for elem in root.getiterator(('group', 'field', 'attribute')):
            try:
                parent = elem.getparent().attrib['type']
            except KeyError:
                # that's an attribute
                parent = elem.getparent().attrib['name']

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
                    logger.debug(f">>> FIELD {elem.attrib['name']} is optional")
                    # optional field, skip
                    continue
                self.create_field(f, elem)
            elif elem.tag == 'attribute':
                if 'optional' in elem.keys() and elem.attrib['optional'] == 'true':
                    continue
                logger.warning(f"Adding attr {elem.attrib['name']}")
                # self.create_attribute(f, elem)
        for path, elem in self.additional_elements.items():
            if elem.type == NxType.field:
                field = elem
                setattr(field, 'full_path', path)
                field.push(f)
                self.stat['fields from add'] += 1

        #     elif elem.type == NxType.attribute:
        #         self.create_attribute(f, elem)


class Ruleset(Agipd2nexus):

    def __init__(self, args):
        # constructor of the parent class first
        Agipd2nexus.__init__(self, args)

        if self.params.detector_distance == None:
            # FIXME: should be taken from the `geom` file
            self.params.detector_distance = 177.0  # Set detector distance arbitrarily if nothing is provided
        self.hierarchy = read_geom(self.params.geom_file)
        self.n_quads = 4
        self.n_modules = 4
        self.n_asics = 8

        self.field_rules = {
            # 'entry/definition': np.str(f"NXmx:{get_git_revision_hash()}"),      # TODO: _create_scalar?
            'entry/definition': np.str(f"NXmx"),      # TODO: _create_scalar?
            'entry/file_name': np.str(self.output_file_name),
            # 'entry/start_time': np.str(self.params.nexus_details.start_time),
            'entry/start_time': np.str('2000-10-10T00:00:00.000Z'),     # FIXME: what is the real data?
            # 'entry/end_time': np.str(self.params.nexus_details.end_time),
            'entry/end_time': np.str('2000-10-10T01:00:00.000Z'),
            # 'entry/end_time_estimated': np.str(self.params.nexus_details.end_time_estimated),
            'entry/end_time_estimated': np.str('2000-10-10T01:00:00.000Z'),
            'entry/data/data': h5py.ExternalLink(self.params.cxi_file, "entry_1/data_1/data"),
            'entry/instrument/detector_group/group_index': np.array(list(range(1, 3)), dtype='i'), # XXX: why 16, not 2?
            'entry/instrument/detector_group/group_names': np.array([np.string_('AGIPD'), np.string_('ELE_D0')],
                                                                    dtype='S12'),
            'entry/instrument/detector_group/group_parent': np.array([-1, 1], dtype='i'),
            'entry/instrument/detector_group/trtansformations/AXIS_D0':
                {'value': 0.0, 'depends_on': 'AXIS_RAIL', 'equipment': 'detector',
                 'equipment_component': 'detector_arm',
                 'transformation_type': 'rotation', 'units': 'degrees', 'vector': (0., 0., -1.),
                 'offset': self.hierarchy.local_origin, 'offset_units': 'mm'},
            'entry/source/name': NexusElement(full_path='entry/source/name', value=self.params.nexus_details.source_name,
                                              nxtype=NxType.field, dtype='s',
                                              attrs={'short_name': self.params.nexus_details.source_short_name}),
        }
        self.group_rules = {
            'NXdetector': {'names': ['ELE_D0']},
            'NXdetector_group': {'names': ['AGIPD']},

        }
        self.additional_elements = {
            'entry/instrument/detector_group/group_type':
                NexusElement(full_path='entry/instrument/detector_group/group_type',
                             nxtype=NxType.field, value=[1, 2], dtype='i')
        }
        self.global_attrs = {
            'NXclass': 'NXroot',
            'file_name': self.output_file_name,
            'file_time': str(dt.now()),
            'HDF5_Version': '1.2.3'
        }

if __name__ == '__main__':
  nexus_helper = Ruleset(sys.argv[1:])
  nexus_helper.create_nexus_master_file()
  logger.info("Stats:\n\t" + "\n\t".join(f"{k}: {v}" for k, v in nexus_helper.stat.items()))
  os.system(f'h5glance {nexus_helper.output_file_name} --attrs')
  # os.system(f"{path_to_cnxvalidate} -l definitions ~/xfel/examples/swissFEL_example/spb/{nexus_helper.output_file_name}"
  #           f" | grep source")
  os.system(f"{path_to_cnxvalidate} -l definitions ~/xfel/examples/swissFEL_example/spb/{nexus_helper.output_file_name}")