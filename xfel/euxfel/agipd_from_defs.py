from __future__ import division
from read_geom import read_geom
from extra_data import RunDirectory

from libtbx.phil import parse
from libtbx.utils import Sorry

import h5py
from h5py import string_dtype as h5py_str
import os
import subprocess
import sys
from pathlib import Path
from lxml import etree, objectify
import logging
import numpy as np
from datetime import datetime as dt
from typing import Union, Any, List
from enum import Enum, auto
from collections import Counter
from itertools import product

LOGGING_FORMAT = "%(levelname)s: %(message)s"
logging.basicConfig(level=logging.INFO, format=LOGGING_FORMAT)
logger = logging.getLogger()

p = Path(__file__)
path_to_defs = p.parent / "definitions/applications/NXmx.nxdl.xml"
path_to_phil = p.parent / "AGIPD.phil"

try:
    root = objectify.parse(str(path_to_defs))
except OSError:
    # clone the submodule with definitions:
    os.system(
        f"cd {p.parent.parent.parent} && git submodule update --init --remote && cd -"
    )
    root = objectify.parse(str(path_to_defs))

for elem in root.getiterator():
    try:
        elem.tag = etree.QName(elem).localname
    except ValueError:
        logger.warning(f"Error with {elem.tag}")

with open(path_to_phil) as phil:
    phil_scope = parse(phil.read())


class NxType(Enum):
    """ Enumeration for elements withing NeXus schema """

    group = auto()
    field = auto()
    attribute = auto()


class LazyFunc:
    """ Very primitive version of lazy executor """

    def __init__(self, *args):
        self.func, self.loc_src, self.loc_dest = args

    def push(self, h5file: h5py.File):
        self.func(self.loc_src, h5file[self.loc_dest])
        # close the source file after usage:
        self.func.__self__.close()


class NexusElement:
    def __init__(
        self,
        full_path: str,
        value: Any,
        nxtype: NxType,
        dtype: str,
        parent: str = "",
        attrs: dict = None,
    ) -> None:
        self.value = value
        self.dtype = dtype
        self.type = nxtype
        self.parent = parent
        self.full_path = full_path
        self.attrs = attrs

    def push(self, h5file):
        """ Write an element to the file """
        if self.full_path:
            h5file.create_dataset(self.full_path, data=self.value, dtype=self.dtype)
            if self.attrs:
                for k, v in self.attrs.items():
                    h5file[self.full_path].attrs[k] = v
        else:
            logger.error(f"Cannot push {self.name}")

    def __str__(self):
        return f"{self.value} [type:{type(self.value)}, dtype:{self.dtype}]"


def get_git_revision_hash() -> str:
    definitions_path = os.path.join(os.path.dirname(sys.argv[0]), "definitions")
    current_dir = os.getcwd()
    os.chdir(definitions_path)
    definitions_hash = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip()
    os.chdir(current_dir)
    return definitions_hash.decode()


class Agipd2nexus:
    def __init__(self, args):
        self.group_rules = {}
        self.field_rules = {}
        self.additional_elements = {}
        self.params_from_phil(args)
        self.output_file_name = (
            os.path.splitext(self.params.cxi_file)[0] + "_master_from_defs.h5"
        )
        self.stat = Counter()

    def params_from_phil(self, args):
        """`args` are command line arguments"""
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

    def translate_groups(self, h5_path: str):
        local_path = h5_path[:]
        for NXname, local_name in self.group_rules.items():
            local_path = local_path.replace(NXname, local_name["names"])
        return local_path

    def get_root_path(self, tree_elem: etree) -> List[str]:
        """ return a NeXus path (e.g. `entry/instrument/detector_group`) from an lxml tree element

        The trickery is that we compare the path's elements with the group names which we want to rename
        and substitute those path's elements accordingly.
        And if for a certain element we need to create several groups, then paths get duplicated with those names:

        Example:
            we have `entry/instrument/detector_group` in the schema, and
            `detector_group` are [`DET_1`, `DET_2`]. Then we've got two paths as the result:
            [`entry/instrument/DET_1`, `entry/instrument/DET_2`]

        """
        ancestor_list = list(x.attrib["type"] for x in tree_elem.iterancestors())
        for j, elem in enumerate(ancestor_list):
            # an `elem` might occur in the `group_rules` no more than once!
            if elem in self.group_rules:
                # for NXname, local_names in self.group_rules.items():
                # check if we use our own names
                if "names" not in self.group_rules[elem]:
                    continue
                local_names = self.group_rules[elem]["names"]
                if len(local_names) == 1:
                    ancestor_list[j] = [local_names[0]]
                else:
                    ancestor_list[j] = []
                    for local_name in local_names:
                        # iterate through many identical groups
                        ancestor_list[j] += [local_name]
            else:
                ancestor_list[j] = [elem]

        ancestor_list = ancestor_list[:-1][
            ::-1
        ]  # drop the last element and reverse the order
        root_path = [
            "/" + "/".join(path).replace("NX", "") for path in product(*ancestor_list)
        ]  # concatenate the values
        return root_path

    @staticmethod
    def get_nxdlPath(h5_path: h5py.Group) -> str:
        par = h5_path.parent.attrs.get("NX_class")
        if par == b"NXentry":
            return par.decode()
        return f"{Agipd2nexus.get_nxdlPath(h5_path.parent)}/{par.decode()}"

    def create_group(
        self, h5file: h5py.File, lx_elem: objectify.ObjectifiedElement, lx_parent: str
    ) -> None:
        """
        Create a new h5.Group attached to the `parent_elem`

        The group name might be provided in the `group_rules`, otherwise it is created on the fly
        either by reducing `NXblah` to `blah` or from the `name` attribute

        Group rules might contain additional actions like creation specific content
        """
        NXname = lx_elem.attrib[
            "type"
        ]  # TODO: the name might be contained in the `name` attribute
        if NXname.startswith("NX"):
            name = NXname.replace("NX", "")
        else:
            logger.warning(f"name {NXname} is not a correct NX-name")
            return
        logger.debug(f"GROUP: {lx_parent} --> {NXname}")
        root_paths = self.get_root_path(lx_elem)
        for root_path in root_paths:
            parent = h5file[root_path]
            if NXname in self.group_rules and "names" in self.group_rules[NXname]:
                for n in self.group_rules[NXname]["names"]:
                    new_group = parent.create_group(n)
                    new_group.attrs["NX_class"] = NXname
                    logger.debug(f"group {n} was created from a rule")
                    self.stat["groups from rules"] += 1
            else:
                new_group = parent.create_group(name)
                new_group.attrs["NX_class"] = NXname
                logger.debug(f"group {name} was created")
                self.stat["groups from defs"] += 1

    def create_field(
        self,
        h5file: h5py.File,
        lx_elem: Union[objectify.ObjectifiedElement, dict],
        recommended: bool = False,
    ) -> None:
        """
        The function tries to add a field to the HDF5 hierarchy.

        If there is a `recommended` marker but the element is not in the `field_rules` then it will be skipped.
        If `recommended == False` then it means the element is mandatory, and even if it is not in the `field_rules`
        it will be filled with a dummy value `7.7777777`
        """
        NXname = lx_elem.attrib["name"]
        name = NXname.replace("NX", "")
        if isinstance(lx_elem, objectify.ObjectifiedElement):
            root_paths = self.get_root_path(lx_elem)
            full_paths = ["/".join([path, name]) for path in root_paths]

        for full_path in full_paths:
            if full_path in self.field_rules:
                logger.debug(
                    f"Add {full_path} from a rule: {self.field_rules[full_path]}"
                )
                if isinstance(self.field_rules[full_path], dict):
                    vector = self.field_rules[full_path]
                    h5file[full_path] = np.array(vector.pop("value"), dtype="f")
                    for key, attribute in vector.items():
                        h5file[full_path].attrs[key] = attribute
                elif isinstance(self.field_rules[full_path], (NexusElement, LazyFunc)):
                    self.field_rules[full_path].push(h5file)
                else:
                    h5file[full_path] = self.field_rules[full_path]
                logger.debug(f"field {full_path} was added")
                self.stat["fields from rules"] += 1

            else:
                if recommended:
                    logger.info(
                        f"Recommended element {full_path} is not in the field rules. Skipped"
                    )
                    continue
                logger.warning(f"Add {full_path} from definition as `7.7777777`")
                field = NexusElement(
                    full_path=full_path, value=7.7777777, nxtype=NxType.field, dtype="f"
                )
                field.push(h5file)
                self.stat["fields from defs"] += 1

    def create_attribute(self, h5file: h5py.File, elem):
        NXname = elem.attrib["name"]
        name = NXname.replace("NX", "")
        if isinstance(elem, objectify.ObjectifiedElement):
            root_path = self.get_root_path(elem)
            full_path = "/".join([root_path, name])

        h5file[full_path] = elem.attrib["value"]
        self.stat["attr"] += 1

    def create_nexus_master_file(self):
        self.out_file = h5py.File(self.output_file_name, "w")
        logger.info(f"file {self.output_file_name} was created")

        for k, v in self.global_attrs.items():
            self.out_file.attrs[k] = v

        for elem in root.getiterator(("group", "field", "attribute")):
            try:
                parent = elem.getparent().attrib["type"]
            except KeyError:
                # that's an attribute
                parent = elem.getparent().attrib["name"]

            if parent == "group":
                # create root element of the file
                entry = self.out_file.create_group("entry")
                entry.attrs["NX_class"] = "NXentry"
                continue

            if elem.tag == "group":
                if (
                    "minOccurs" in elem.keys()
                    and elem.attrib["minOccurs"] == "0"
                    and elem.attrib["type"] not in self.group_rules
                ):
                    continue
                self.create_group(self.out_file, elem, parent)
            elif elem.tag == "field":
                if ("minOccurs" in elem.keys() and elem.attrib["minOccurs"] == "0") or (
                    "optional" in elem.keys() and elem.attrib["optional"] == "true"
                ):
                    logger.debug(f">>> FIELD {elem.attrib['name']} is optional")
                    # an optional field, skip
                    continue
                elif (
                    "recommended" in elem.keys()
                    and elem.attrib["recommended"] == "true"
                ):
                    self.create_field(self.out_file, elem, recommended=True)
                else:
                    self.create_field(self.out_file, elem)
            elif elem.tag == "attribute":
                if "optional" in elem.keys() and elem.attrib["optional"] == "true":
                    continue
                logger.debug(f"Adding attr {elem.attrib['name']}")
        for path, elem in self.additional_elements.items():
            if elem.type == NxType.field:
                field = elem
                setattr(field, "full_path", path)
                field.push(self.out_file)
                logger.debug(f"Additional elem was added to {path}")
                self.stat["fields from add"] += 1
        self.out_file.close()


class Ruleset(Agipd2nexus):
    def __init__(self, args):
        # constructor of the parent class first
        Agipd2nexus.__init__(self, args)

        self.hierarchy = read_geom(self.params.geom_file)

        if self.params.cxi_file and Path(self.params.cxi_file).exists():
            cxi = h5py.File(self.params.cxi_file, "r")
        else:
            cxi = None

        if self.params.detector_distance is None:
            # try to take from the geometry file:
            if self.hierarchy.detector_distance:
                self.params.detector_distance = self.hierarchy.detector_distance
            else:
                raise Sorry(
                    "Detector distance is undefined! You should set it either in `.phil` or in `.geom` files, "
                    "or pass as a command line argument: `detector_distance=123.45` (in mm)"
                )
        if self.params.wavelength is None:
            # try to take from the geometry file:
            if self.hierarchy.incident_wavelength:
                self.params.wavelength = self.hierarchy.incident_wavelength
            else:
                raise Sorry(
                    "Incident wavelength is undefined! You should set it either in `.phil` or in `.geom` files, "
                    "or pass as a command line argument: `wavelength=1.2345` (in angstrom)"
                )
        self.n_quads = 4
        self.n_modules = 4
        self.n_asics = 8

        # ==== CREATE DETECTOR MODULES ====
        """
        Add 4 quadrants
        Nexus coordiate system, into the board            AGIPD detector
             o --------> (+x)                             Q3=(12,13,14,15) Q0=(0,1,2,3)
             |                                                        o
             |                                            Q2=(8,9,10,11)   Q1=(4,5,6,7)
             v
            (+y)
        """
        panels = []
        for q, quad in self.hierarchy.items():
            for m, module in quad.items():
                panels.extend([module[key] for key in module])
        fast = max([int(panel["max_fs"]) for panel in panels]) + 1
        slow = max([int(panel["max_ss"]) for panel in panels]) + 1
        pixel_size = panels[0]["pixel_size"]
        assert [
            pixel_size == panels[i + 1]["pixel_size"] for i in range(len(panels) - 1)
        ].count(False) == 0

        quad_fast = fast
        quad_slow = slow * self.n_quads
        module_fast = quad_fast
        module_slow = quad_slow // self.n_quads
        asic_fast = module_fast
        asic_slow = module_slow // self.n_asics

        self.group_rules = {
            "NXdetector": {"names": ["ELE_D0"]},
            "NXdetector_group": {"names": ["AGIPD"]},
            "NXtransformations": {},
            "NXdetector_module": {"names": []},  # 'names' will be populated below
        }
        array_name = "ARRAY_D0"
        det_path = "/entry/instrument/ELE_D0/"
        t_path = det_path + "transformations/"

        class Transform(NexusElement):
            def __init__(
                self, name: str, value: Any = [0.0], attrs: dict = None
            ) -> None:
                default_attrs = {
                    "equipment": "detector",
                    "transformation_type": "rotation",
                    "units": "degrees",
                    "offset_units": "mm",
                    "vector": (0.0, 0.0, -1.0),
                }
                NexusElement.__init__(
                    self,
                    full_path=t_path + name,
                    value=value,
                    nxtype=NxType.field,
                    dtype="f",
                    attrs={**default_attrs, **attrs},
                )

        det_field_rules = {}  # extends mandatory field rules
        det_additional_rules = (
            {}
        )  # additional transformation elements (fields) for the detector

        for quad in range(self.n_quads):  # iterate quadrants
            q_key = f"q{quad}"
            q_name = f"AXIS_D0Q{quad}"
            quad_vector = self.hierarchy[q_key].local_origin.elems

            q_elem = Transform(
                q_name,
                attrs={
                    "depends_on": t_path + "AXIS_D0",
                    "offset": quad_vector,
                    "equipment_component": "detector_quad",
                },
            )
            det_additional_rules[t_path + q_name] = q_elem
            for module_num in range(
                self.n_modules
            ):  # iterate modules within a quadrant
                m_key = f"p{(quad * self.n_modules) + module_num}"
                m_name = f"AXIS_D0Q{quad}M{module_num}"
                module_vector = self.hierarchy[q_key][m_key].local_origin.elems
                m_elem = Transform(
                    m_name,
                    attrs={
                        "depends_on": t_path + q_name,
                        "equipment_component": "detector_module",
                        "offset": module_vector,
                    },
                )
                det_additional_rules[t_path + m_name] = m_elem
                for asic_num in range(self.n_asics):  # iterate asics within a module
                    a_key = f"p{(quad * self.n_modules) + module_num}a{asic_num}"
                    a_name = f"AXIS_D0Q{quad}M{module_num}A{asic_num}"
                    asic_vector = self.hierarchy[q_key][m_key][a_key][
                        "local_origin"
                    ].elems

                    a_elem = Transform(
                        a_name,
                        attrs={
                            "depends_on": t_path + m_name,
                            "equipment_component": "detector_asic",
                            "offset": asic_vector,
                        },
                    )
                    det_additional_rules[t_path + a_name] = a_elem
                    det_module_name = array_name + f"Q{quad}M{module_num}A{asic_num}"
                    # populate ``group_rules`` with detector modules:
                    self.group_rules["NXdetector_module"]["names"] += [det_module_name]

                    def full_m_name(name: str) -> str:
                        return det_path + det_module_name + "/" + name

                    det_field_rules[full_m_name("data_origin")] = np.array(
                        [(quad * self.n_modules) + module_num, asic_slow * asic_num, 0],
                        dtype="i",
                    )
                    det_field_rules[full_m_name("data_size")] = np.array(
                        [1, asic_slow, asic_fast], dtype="i"
                    )

                    fast = self.hierarchy[q_key][m_key][a_key]["local_fast"].elems
                    slow = self.hierarchy[q_key][m_key][a_key]["local_slow"].elems

                    det_field_rules[full_m_name("fast_pixel_direction")] = NexusElement(
                        full_path=full_m_name("fast_pixel_direction"),
                        value=[pixel_size],
                        dtype="f",
                        nxtype=NxType.field,
                        attrs={
                            "depends_on": t_path
                            + f"AXIS_D0Q{quad}M{module_num}A{asic_num}",
                            "transformation_type": "translation",
                            "units": "mm",
                            "vector": fast,
                            "offset": (0.0, 0.0, 0.0),
                        },
                    )
                    det_field_rules[full_m_name("slow_pixel_direction")] = NexusElement(
                        full_path=full_m_name("slow_pixel_direction"),
                        value=[pixel_size],
                        dtype="f",
                        nxtype=NxType.field,
                        attrs={
                            "depends_on": t_path
                            + f"AXIS_D0Q{quad}M{module_num}A{asic_num}",
                            "transformation_type": "translation",
                            "units": "mm",
                            "vector": slow,
                            "offset": (0.0, 0.0, 0.0),
                        },
                    )
        self.field_rules = {
            # '/entry/definition': np.string_(f"NXmx:{get_git_revision_hash()}"),      # TODO: _create_scalar?
            "/entry/definition": np.string_(
                "NXmx"
            ),  # XXX: whoa! this is THE criteria of being a "nexus format"    !
            "/entry/file_name": np.string_(self.output_file_name),
            # '/entry/start_time': np.string_(self.params.nexus_details.start_time),
            "/entry/start_time": np.string_(
                "2000-10-10T00:00:00.000Z"
            ),  # FIXME: what is the real data?
            # '/entry/end_time': np.string_(self.params.nexus_details.end_time),
            "/entry/end_time": np.string_("2000-10-10T01:00:00.000Z"),
            # '/entry/end_time_estimated': np.string_(self.params.nexus_details.end_time_estimated),
            "/entry/end_time_estimated": np.string_("2000-10-10T01:00:00.000Z"),
            "/entry/data/data": LazyFunc(cxi.copy, "entry_1/data_1/data", "entry/data"),
            "/entry/instrument/name": NexusElement(
                full_path="/entry/instrument/name",
                value=self.params.nexus_details.instrument_name,
                nxtype=NxType.field,
                dtype=h5py_str(),
                attrs={"short_name": self.params.nexus_details.instrument_short_name},
            ),
            "/entry/instrument/AGIPD/group_index": np.array(
                list(range(1, 3)), dtype="i"
            ),  # XXX: why 16, not 2?
            "/entry/instrument/AGIPD/group_names": np.array(
                [np.string_("AGIPD"), np.string_("ELE_D0")], dtype="S12"
            ),
            "/entry/instrument/AGIPD/group_parent": np.array([-1, 1], dtype="i"),
            "/entry/instrument/beam/incident_wavelength": NexusElement(
                full_path="/entry/instrument/beam/incident_wavelength",
                value=self.params.wavelength,
                nxtype=NxType.field,
                dtype="f",
                attrs={"units": "angstrom"},
            ),
            "/entry/instrument/beam/total_flux": NexusElement(
                full_path="/entry/instrument/beam/total_flux",
                value=self.get_xgm_data(cxi),
                nxtype=NxType.field,
                dtype="f",
                attrs={"units": "Hz"},
            ),
            "/entry/instrument/ELE_D0/data": h5py.SoftLink("/entry/data/data"),
            "/entry/instrument/ELE_D0/sensor_material": "Si",  # FIXME: move to the `phil`-file
            "/entry/instrument/ELE_D0/sensor_thickness": NexusElement(
                full_path="/entry/instrument/ELE_D0/sensor_thickness",
                value=300.0,  # FIXME: move to the `phil`-file
                nxtype=NxType.field,
                dtype="f",
                attrs={"units": "microns"},
            ),
            "/entry/sample/depends_on": np.str("."),  # XXX: Why not `np.string_`??
            "/entry/sample/name": NexusElement(
                full_path="/entry/sample/name",
                value=self.params.sample_name,
                nxtype=NxType.field,
                dtype=h5py_str(),
            ),
            "/entry/source/name": NexusElement(
                full_path="/entry/source/name",
                value=self.params.nexus_details.source_name,
                nxtype=NxType.field,
                dtype=h5py_str(),
                attrs={"short_name": self.params.nexus_details.source_short_name},
            ),
        }

        self.field_rules = {**self.field_rules, **det_field_rules}
        self.additional_elements = {
            "/entry/instrument/AGIPD/group_type": NexusElement(
                full_path="/entry/instrument/AGIPD/group_type",
                value=[1, 2],
                nxtype=NxType.field,
                dtype="i",
            ),
            f"{t_path}/AXIS_D0": Transform(
                "AXIS_D0",
                attrs={
                    "depends_on": t_path + "AXIS_RAIL",
                    "equipment_component": "detector_arm",
                    "offset": self.hierarchy.local_origin,
                },
            ),
            f"{t_path}/AXIS_RAIL": NexusElement(
                full_path=t_path + "AXIS_RAIL",
                dtype="f",
                nxtype=NxType.field,
                value=[self.params.detector_distance],
                attrs={
                    "depends_on": np.string_("."),
                    "equipment": "detector",
                    "equipment_component": "detector_arm",
                    "transformation_type": "translation",
                    "units": "mm",
                    "vector": (0.0, 0.0, 1.0),
                },
            ),
        }
        self.additional_elements = {**self.additional_elements, **det_additional_rules}
        self.global_attrs = {
            "NX_class": "NXroot",
            "file_name": self.output_file_name,
            "file_time": str(dt.now()),
            "HDF5_Version": h5py.version.hdf5_version,
        }

    def get_xgm_data(self, cxi):
        run_raw = RunDirectory(self.params.xgm_dir)
        # Photon energy: E = hc/λ
        hc = 1.9864458571489e-25  # [joule * meter]
        # wavelength is in ångström, it gives the 1E10 factor
        pe = hc * 1e10 / self.params.wavelength  # photon energy in joules
        tids = cxi["entry_1/trainId"]
        intensity = run_raw.get_array(
            self.params.xgm_addr, "data.intensityTD", extra_dims=["pulseID"]
        )
        xgm_data = intensity[
            (intensity["trainId"] >= tids[0]) & (intensity["trainId"] <= tids[-1])
        ][:, : self.params.xgm_pulse_num]
        # XGM data is in μJ, it gives the 1E-6 factor
        # Finally, we normalize this per 220 nanoseconds (2.2E-7)
        return xgm_data * 1e-6 / pe / 2.2e-7


if __name__ == "__main__":
    nexus_helper = Ruleset(sys.argv[1:])
    nexus_helper.create_nexus_master_file()
    logger.info(
        "Stats:\n\t" + "\n\t".join(f"{k}: {v}" for k, v in nexus_helper.stat.items())
    )
