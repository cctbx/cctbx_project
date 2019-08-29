#!/usr/bin/env python

"""
Lists the module python dependencies and the state of install.
"""

from __future__ import absolute_import, division, print_function

import argparse
import copy
import re
import sys

import libtbx.load_env
import libtbx.pkg_utils
from libtbx.pkg_utils import pkg_resources
from libtbx.pkg_utils import packaging
from six.moves import range
from six.moves import zip


BOLD = "\033[1m"
RED = "\033[1;31m"
MAGENTA = "\033[1;35m"
GREEN = "\033[32m"
GRAY = "\033[37m"
NC = "\033[0m"

# A regex to check for requirement characters that are not shell-safe
reSafe = re.compile(r"[^a-zA-Z0-9,._+:@%/\[\]-]")


def _get_requirement_status(requirement):
    # type: (packaging.requirements.Requirement) -> Tuple[Any, str]
    """Get the requirement status of a package.

    Args:
        requirement (packaging.requirements.Requirement):
            The requirement object to check

    Returns:
        Tuple[str, str]:
            The status value, and - if it exists - the current version.
            The status value can currently be one of:
                - "SKIPPED": The package is excluded by env markers
                - "VALID":   Already installed with correct version
                - "UNKNOWNEXTRA": Installed, but missing extras
                - "MISSING": The package is not installed
                - "MISMATCH": An incompatible version is installed
    """
    requirement = copy.deepcopy(requirement)
    if requirement.marker and not requirement.marker.evaluate():
        # Try to get a version without the marker
        requirement.marker = None
        # Look to see if we have an installed version anyway
        try:
            version = pkg_resources.require(str(requirement))[0].version
            return ("SKIPPED", version)
        except pkg_resources.ResolutionError:
            return ("SKIPPED", None)

    # Try to find this package
    try:
        pkg = pkg_resources.require(str(requirement))[0]
        return ("VALID", pkg.version)
    except pkg_resources.UnknownExtra:
        # Package exists, but extra is missing. Get the base description.
        requirement.extras = set()
        version = pkg_resources.require(str(requirement))[0].version
        return ("UNKNOWNEXTRA", version)
    except pkg_resources.VersionConflict:
        requirement.specifier = packaging.specifiers.SpecifierSet()
        version = pkg_resources.require(str(requirement))[0].version
        return ("MISMATCH", version)
    except pkg_resources.DistributionNotFound:
        # Package not installed
        return ("MISSING", None)


def _shell_quote_requirement(req):
    # type: (Union[packaging.requirements.Requirement, str]) -> str
    """Returns a shell-quoted requirement string.

    Since we know the form of the requirement, it is safe to single-quote
    anything and convert internal single-quotes to double (we know no
    recursive quoting). But if it's not necessary to do this, it won't.

    Args:
        req (Requirement or str):
            The requirement, as an object or PEP508 string.

    Returns:
        str:
            A shell-safe, potentially escaped version of the string.
    """
    out = str(req)
    if reSafe.search(out):
        # Something needs to be quoted. Quote it, but we know the
        # general form and know we don't have recursively nested quotes
        out = "'" + out.replace("'", '"') + "'"
    return out


def _print_table(rows):
    # type: (List[List[str]]) -> None
    """Prints a pretty formatted table of dependencies.

    Args:
        rows (List[List[str]]):
            The row data, containing six columns in the order of
            Status, Name, Specifier, Available, Marker, From
    """

    titles = ["Status", "Name", "Specifier", "Available", "Marker", "From"]
    justification = ["l", "l", "l", "r", "l", "l"]

    # Remove marker row if no markers
    if not any(x[4] for x in rows):
        for row in rows + [titles] + [justification]:
            row.pop(4)

    column_widths = [
        max(len(str(row[i])) for row in rows + [titles]) for i in range(len(rows[0]))
    ]

    colors = {
        "VALID": GREEN,
        "MISSING": RED,
        "UNKNOWNEXTRA": MAGENTA,
        "MISMATCH": MAGENTA,
        "SKIPPED": GRAY,
        "Status": "",
    }
    # Print a formatted table
    for row in [titles] + rows:
        color = colors[row[0]]
        out_cols = []
        for col, just, width in zip(row[1:], justification[1:], column_widths[1:]):
            col_txt = col.ljust(width)
            if just == "r":
                col_txt = col.rjust(width)
            elif just == "c":
                col_txt = col.center(width)
            out_cols.append(col_txt)
        print(color + " ".join(out_cols) + (NC if color else ""))


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Show python package requirements")
    parser.add_argument(
        "--norefresh", action="store_true", help="Don't re-read module config"
    )
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument(
        "--all", action="store_true", help="Just show a list of all dependencies"
    )
    modes.add_argument(
        "--actions", action="store_true", help="Show packages to install/upgrade"
    )
    args = parser.parse_args()

    # Should we do the normal actions if nothing special has been asked for?
    normal_output = not args.all and not args.actions

    # Check for packaging
    if not libtbx.pkg_utils.pkg_resources:
        print(
            "\n".join(
                [
                    "  WARNING: Can not verify python package requirements - pip/setuptools out of date",
                    "  Please update pip and setuptools by running:",
                    "",
                    "    libtbx.python -m pip install pip setuptools --upgrade",
                    "",
                    "  or following the instructions at https://pip.pypa.io/en/stable/installing/",
                ]
            )
        )
        sys.exit(1)

    # Re-read all config files, without a refresh. Not entirely sure that
    # this is the proper thing to do as makes a slight race condition with
    # refresh, and cannot get a stage where module is listed but not refreshed?
    if not args.norefresh:
        for module in libtbx.env.module_list:
            module.process_libtbx_config()

    # Gather the requirements
    requirements = sorted(
        libtbx.pkg_utils.collate_python_requirements(libtbx.env.module_list),
        key=lambda x: x.name.lower(),
    )

    rows = []
    to_action = []

    for requirement in requirements:
        # Check the status of this requirement
        status, cur_ver = _get_requirement_status(requirement)
        name = requirement.name
        if requirement.extras:
            name = name + "[" + ",".join(requirement.extras) + "]"
        # Build a row object for table display
        rows.append(
            [
                status,
                name,
                str(requirement.specifier),
                str(cur_ver) if cur_ver is not None else "",
                str(requirement.marker) if requirement.marker else "",
                ", ".join(sorted(requirement.modules)),
            ]
        )
        if status in {"MISSING", "UNKNOWNEXTRA", "MISMATCH"}:
            to_action.append(requirement)

    # Decide how to do output
    if normal_output and rows:
        _print_table(rows)

        if to_action:
            print(
                "\nThere are packaging actions required. pip command to resolve\npackaging differences:\n"
            )
            print(
                " ".join(
                    ["    libtbx.pip", "install"]
                    + [_shell_quote_requirement(req) for req in to_action]
                )
            )
    elif normal_output and not rows:
        print("No dependencies found.")
    elif args.all:
        print(" ".join(_shell_quote_requirement(x) for x in requirements))
    elif args.actions:
        print(" ".join(_shell_quote_requirement(x) for x in to_action))
