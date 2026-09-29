#!/usr/bin/env python3
"""Make the adopted Outside Reviewer transport bundle from a frozen packet."""

import sys
if __name__ == "__main__" and not sys.flags.isolated:
    sys.path.pop(0)

import argparse
from contextlib import ExitStack
import hashlib
import importlib.util
import io
import os
from pathlib import Path
import tarfile
import tempfile

# Load the reviewed sibling explicitly. With -I, the tools directory is not
# on sys.path, so an unlisted local module cannot shadow the standard library.
_checker_path = Path(__file__).with_name("screen_check.py")
_checker_spec = importlib.util.spec_from_file_location("screen_check", _checker_path)
_checker = importlib.util.module_from_spec(_checker_spec)
_checker_spec.loader.exec_module(_checker)
regular_files, verify = _checker.regular_files, _checker.verify


def sha(data):
    return hashlib.sha256(data).hexdigest()


def add_bytes(archive, name, data):
    info = tarfile.TarInfo(name)
    info.size = len(data)
    info.mode = 0o644
    info.mtime = 0
    archive.addfile(info, io.BytesIO(data))


def directory_identity(stat_result):
    return stat_result.st_dev, stat_result.st_ino


def require_outside_companions(packet_root, companions):
    """Companions and every ancestor must be separate from the packet."""
    companion_root = companions.resolve(strict=True)
    for directory in (companion_root, *companion_root.parents):
        if directory.samefile(packet_root):
            raise ValueError("companions must be outside the frozen packet")


def make_bundle(packet, companions, output):
    # Only the packet's own parent can receive a bundle.  A separate output
    # directory could be moved under the packet after a successful preflight,
    # even if its descriptor stays open.  The packet's parent cannot move
    # under its own child while the packet remains in that parent.
    with ExitStack() as stack:
        packet_root = packet.resolve(strict=True)
        require_outside_companions(packet_root, companions)
        packet_parent_fd = os.open(packet_root.parent, os.O_RDONLY | os.O_DIRECTORY)
        stack.callback(os.close, packet_parent_fd)
        packet_identity = directory_identity(
            os.stat(packet_root.name, dir_fd=packet_parent_fd))
        parent_fd = os.open(output.parent, os.O_RDONLY | os.O_DIRECTORY)
        stack.callback(os.close, parent_fd)
        if directory_identity(os.fstat(parent_fd)) != directory_identity(
                os.fstat(packet_parent_fd)):
            raise ValueError("output must be beside the frozen packet")
        packet_id = verify(packet)
        files = regular_files(packet)
        names = {name for name, _ in files}
        if not {"README_FIRST.md", "PROOF_SUMMARY.md"}.issubset(names):
            raise ValueError("packet needs README_FIRST.md and PROOF_SUMMARY.md before freeze")
        if any(Path(name).name in {"APPROVALS.md", "CHECKER.md", "ERRATA.md"}
               for name in names):
            raise ValueError("companion file present inside frozen packet")
        companions_present = regular_files(companions)
        companion_names = {name for name, _ in companions_present}
        if not {"APPROVALS.md", "CHECKER.md"}.issubset(companion_names):
            raise ValueError("outside companions need APPROVALS.md and CHECKER.md")
        if companion_names - {"APPROVALS.md", "CHECKER.md", "ERRATA.md"}:
            raise ValueError("unrecognized outside companion")
        with tempfile.TemporaryDirectory(prefix="review-bundle-") as temporary:
            inner = Path(temporary) / "PACKET.tgz"
            with tarfile.open(inner, "w:gz") as archive:
                for name, path in files:
                    add_bytes(archive, name, path.read_bytes())
            payload = {"PACKET.tgz": inner.read_bytes()}
            for name, path in companions_present:
                payload[name] = path.read_bytes()
            checksums = "".join(f"{sha(payload[name])}  {name}\n" for name in sorted(payload))
            payload["SHA256SUMS"] = checksums.encode("utf-8")
            require_outside_companions(packet_root, companions)
            if directory_identity(os.stat(packet_root.name, dir_fd=parent_fd)) != packet_identity:
                raise ValueError("packet moved while bundling")
            try:
                output_fd = os.open(output.name, os.O_CREAT | os.O_EXCL | os.O_RDWR,
                                    0o644, dir_fd=parent_fd)
            except FileExistsError as error:
                raise ValueError("output already exists") from error
            with os.fdopen(output_fd, "w+b") as stream:
                with tarfile.open(fileobj=stream, mode="w:gz") as archive:
                    for name in sorted(payload):
                        add_bytes(archive, name, payload[name])
                stream.flush()
                stream.seek(0)
                bundle_hash = sha(stream.read())
        # Source files might change during assembly; reject that case.
        if verify(packet) != packet_id:
            raise ValueError("packet changed while bundling; do not send the archive")
    print(f"BUNDLE {bundle_hash}")
    print(f"Packet identity: {packet_id}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    parser.add_argument("companions", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        make_bundle(args.packet, args.companions, args.output)
    except (ValueError, OSError) as error:
        parser.exit(2, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
