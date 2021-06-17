from __future__ import division
import sys, os
from scitbx.matrix import col
from libtbx.phil import parse
from libtbx.utils import Sorry

help_str = """Converts a CrystFEL file to DIALS json format."""

phil_scope = parse(
    """
  geom_file = None
    .type = str
    .help = CrystFEL geometry file to convert
  show_plot = False
    .type = bool
    .help = plot of detector geometry
"""
)


class PanelGroup(dict):
    def __init__(self):
        self.center = None
        self.detector_distance = None  # in mm
        self.incident_wavelength = None  # in angstrom
        self.local_origin = None
        self.local_fast = col((1, 0, 0))
        self.local_slow = col((0, 1, 0))

    def setup_centers(self) -> None:
        if self.center is not None:  # the center is already defined
            return
        center = col((0.0, 0.0, 0.0))
        for key, child in self.items():
            if isinstance(child, PanelGroup):
                if child.center is None:
                    child.setup_centers()
                center += child.center
            else:
                center += child["center"]
        center /= len(self)
        self.center = center

    def setup_local_frames(self) -> None:
        for key, child in self.items():
            if isinstance(child, PanelGroup):
                child.local_origin = child.center - self.center
                child.setup_local_frames()
            else:
                child["local_origin"] = child["origin"] - self.center


def read_geom(geom_file: str) -> PanelGroup:
    """Function to read in the CrystFEL .geom file."""
    panels = {}
    rigid_groups = {}
    collections = {}

    def known_panels() -> set:
        all_keys = []
        for value in rigid_groups.values():
            all_keys.extend(value)
        for value in collections.values():
            all_keys.extend(value)
        return set(all_keys)

    pixel_size = None

    with open(geom_file) as geom:
        lines = geom.readlines()

    lines = (line.split(";")[0] for line in lines)  # cut out comments
    lines = (line for line in lines if len(line.split("=")) == 2)
    geometry = dict(map(lambda x: x.strip(), line.split("=")) for line in lines)  # noqa

    if "res" in geometry:
        pixel_size = 1000 / float(geometry.pop("res"))  # mm
    else:
        raise KeyError("Pixel size is not defined!")

    for key, value in geometry.items():
        if "rigid_group" in key:
            if "collection" in key:
                collections[key.split("rigid_group_collection_")[1]] = value.split(",")
            else:
                rigid_groups[key.split("rigid_group_")[1]] = value.split(",")
        else:
            if "/" not in key:
                continue
            panel = key.split("/")[0].strip()
            key = key.split("/")[1].strip()
            if panel not in known_panels():
                continue
            if panel not in panels:
                panels[panel] = {}
            panels[panel][key] = value

    mapping = {}
    for panel in panels:
        mapping[panel] = {}
        for group in rigid_groups:
            if panel not in rigid_groups[group]:
                continue
            for collection in collections:
                if group in collections[collection]:
                    mapping[panel][collection] = group, len(rigid_groups[group])
    # example of mapping entry: mapping['p0a0'] = {'asics': ('p0', 8), 'quadrants': ('q0', 32)}
    parents = {}
    for panel in panels:
        parents[panel] = [
            mapping[panel][k][0]
            for k in sorted(mapping[panel], key=lambda x: x[1], reverse=True)
        ]
    # example of parents entry:  parents['p0a0'] = ['q0', 'p0']
    # IE parents are listed in reverse order of immediacy (p0 is the parent of p0a0 and q0 is the parent of p0)

    hierarchy = PanelGroup()

    if "clen" in geometry:
        hierarchy.detector_distance = float(geometry.pop("clen")) * 1000
    if (
        "photon_energy" in geometry
    ):  # h * c / e = 1.23984198E-6 [SI] -- eV to angstrom which itself is [1E-10 SI]
        hierarchy.incident_wavelength = 1.23984198e4 / float(
            geometry.pop("photon_energy")
        )

    def add_node(panel, parent, parents, depth):
        if depth == len(parents):
            parent[panel] = panels[panel]
        else:
            if parents[depth] not in parent:
                parent[parents[depth]] = PanelGroup()
            add_node(panel, parent[parents[depth]], parents, depth + 1)

    for panel in panels:
        add_node(panel, hierarchy, parents[panel], 0)

    # example of a hierarchy entry:
    # hierarchy['q0']['p0']['p0a0'] = full panel dictionary

    def parse_vector(vector):
        try:
            x, y, z = vector.split(" ")
        except ValueError:
            x, y = vector.split(" ")
            return col((float(x.rstrip("x")), float(y.rstrip("y")), 0.0))
        else:
            return col(
                (float(x.rstrip("x")), float(y.rstrip("y")), float(z.rstrip("z")))
            )

    # set up panel vectors in lab space
    for panel in panels:
        panels[panel]["origin"] = col(
            (
                float(panels[panel]["corner_x"]) * pixel_size,
                float(panels[panel]["corner_y"]) * pixel_size,
                0.0,
            )
        )
        if "coffset" in panels[panel]:
            panels[panel]["origin"] += col(
                (0, 0, 1000 * float(panels[panel]["coffset"]))
            )
        panels[panel]["fast"] = panels[panel]["local_fast"] = parse_vector(
            panels[panel]["fs"]
        ).normalize()
        panels[panel]["slow"] = panels[panel]["local_slow"] = parse_vector(
            panels[panel]["ss"]
        ).normalize()
        center_fast = (
            panels[panel]["fast"]
            * pixel_size
            * (int(panels[panel]["max_fs"]) - int(panels[panel]["min_fs"]) + 1)
            / 2.0
        )
        center_slow = (
            panels[panel]["slow"]
            * pixel_size
            * (int(panels[panel]["max_ss"]) - int(panels[panel]["min_ss"]) + 1)
            / 2.0
        )
        panels[panel]["center"] = panels[panel]["origin"] + center_fast + center_slow
        panels[panel]["pixel_size"] = pixel_size

    assert "pixel_size" not in panels

    hierarchy.setup_centers()
    hierarchy.local_origin = hierarchy.center
    hierarchy.setup_local_frames()
    return hierarchy


def run(args):
    if "-h" in args or "--help" in args or "-c" in args:
        print(help_str)
        phil_scope.show(attributes_level=2)
        return

    user_phil = []
    for arg in args:
        if os.path.isfile(arg):
            user_phil.append(parse(file_name=arg))
        else:
            try:
                user_phil.append(parse(arg))
            except Exception as e:
                raise Sorry("Unrecognized argument: %s" % arg)
    params = phil_scope.fetch(sources=user_phil).extract()

    hierarchy = read_geom(params.geom_file)

    # Plot the detector model highlighting the hierarchical structure of the detector
    def plot_node(cummulative, node, name):
        if isinstance(node, PanelGroup):
            plt.arrow(
                cummulative[0],
                cummulative[1],
                node.local_origin[0],
                node.local_origin[1],
            )
            for childname, child in node.items():
                plot_node(cummulative + node.local_origin, child, childname)
        else:
            plt.arrow(
                cummulative[0],
                cummulative[1],
                node["local_origin"][0],
                node["local_origin"][1],
            )

            ori = node["origin"]
            fast_at_zero = (
                node["fast"]
                * node["pixel_size"]
                * (int(node["max_fs"]) - int(node["min_fs"]) + 1)
            )
            slow_at_zero = (
                node["slow"]
                * node["pixel_size"]
                * (int(node["max_ss"]) - int(node["min_ss"]) + 1)
            )
            plt.arrow(ori[0], ori[1], fast_at_zero[0], fast_at_zero[1], color="blue")
            plt.arrow(ori[0], ori[1], slow_at_zero[0], slow_at_zero[1], color="red")

            plt.text(ori[0], ori[1], name)

    if params.show_plot:
        from matplotlib import pyplot as plt

        plot_node(col((0, 0, 0)), hierarchy, "root")
        plt.xlim(-200, 200)
        plt.ylim(200, -200)

        plt.show()


if __name__ == "__main__":
    run(sys.argv[1:])
