#!/usr/bin/env python3
"""Normalized digests for a test installation's fingerprint (Revision 13). Report-only.

A rebuild that restores an installation to its recorded commits can still change a few
files' bytes only through known timestamp fields:
- static archives (for example lib/*.a): every member header carries a date, and macOS
  ranlib writes the current time into the index member's header (__.SYMDEF);
- a library that embeds its build time as text (for example phaser_ext.so, generated
  phaser/Version_shared).

This tool computes a digest of such a file that ignores only those fields, so that a project
fingerprint script can print, besides its exact FINGERPRINT, a NORMALIZED value. The
queue helper uses the normalized value only to say "timestamp-only differences: attended
validation needed"; it never makes an installation available from it, and nothing here
replaces a recorded fingerprint.

What is ignored, exactly:
- archive (`!<arch>` format, BSD or GNU): the 12-byte date field of each member header.
  Member names, owner, group, mode, size, contents and order, and the index member's
  contents, are all compared. Anything that does not parse as such an archive (a thin
  archive included) is compared byte for byte.
- build-time file: every text date of the form `Sun Oct 11 08:54:45 2026` (C asctime form)
  is masked; in a 64-bit Mach-O file also the LC_UUID value and the LC_CODE_SIGNATURE data,
  which the linker derives from the file's other bytes (so they change whenever the date
  does). The signature is masked only where the linker puts it (checkpoint 5): inside the
  __LINKEDIT segment's file range, ending at the end of that segment and of the file, and
  after every other segment's file range; a signature declared anywhere else makes the
  whole file compared byte for byte (nothing masked; stderr says so). Any other byte is
  compared. Other formats (ELF, universal binaries) get only the
  date masked, so a derived build id there still makes the value differ (conservative).
These rules guard against mistakes, not deliberate edits: the ranges are read from the
file's own load commands, so a file whose load commands were edited on purpose can be
given the same normalized value. That only yields "attended validation needed"; an
installation is available only when its exact recorded fingerprint matches, and a new
value is recorded only by an attended validation.

Usage:
  fingerprint_normalize.py file --kind archive|build-time PATH
      prints `<digest>  <PATH>`; exit 0 (a file that does not parse as that kind is
      digested byte for byte; stderr says so).
  fingerprint_normalize.py detail --root DIR [--archives GLOB]... [--build-time REL]... < DETAIL
      reads fingerprint detail lines (`<sha256>  <path relative to DIR>`, as printed by
      `shasum -a 256`; other lines are copied unchanged) and writes them with the digest of
      each matching file replaced by `archive:<digest>` or `build-time:<digest>`.
      A GLOB is matched part by part (`lib/*.a` matches lib/x.a, not lib/sub/x.a).
Exit 2 for a usage error or an unreadable listed file. Standard library only.
"""

import argparse
import fnmatch
import hashlib
import re
import struct
import sys
from pathlib import Path, PurePosixPath

AR_MAGIC = b"!<arch>\n"
HEADER = 60
DATE = re.compile(rb"(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun) "
                  rb"(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) "
                  rb"[ 0-3][0-9] [0-2][0-9]:[0-5][0-9]:[0-6][0-9] [0-9]{4}")
MACHO_64 = b"\xcf\xfa\xed\xfe"
LC_SEGMENT_64, LC_UUID, LC_CODE_SIGNATURE = 0x19, 0x1b, 0x1d


def archive_digest(data):
    """(digest, members) of an `!<arch>` archive ignoring only member header dates, or None."""
    if not data.startswith(AR_MAGIC):
        return None
    digest = hashlib.sha256(AR_MAGIC)
    pos, members = len(AR_MAGIC), 0
    while pos < len(data):
        header = data[pos:pos + HEADER]
        if len(header) < HEADER or header[58:60] != b"`\n":
            return None
        size_field = header[48:58].strip()
        if not size_field.isdigit():
            return None
        size = int(size_field)
        end = pos + HEADER + size
        if end > len(data):
            return None
        pad = data[end:end + (size % 2)]
        # name[0:16] date[16:28] uid[28:34] gid[34:40] mode[40:48] size[48:58] fmag[58:60]
        digest.update(header[:16] + b" " * 12 + header[28:])
        digest.update(data[pos + HEADER:end] + pad)
        pos = end + len(pad)
        members += 1
    return digest.hexdigest(), members


def build_time_digest(data):
    """(digest, dates masked, derived fields masked) of a file with an embedded build time."""
    masked = bytearray(data)
    dates = 0
    for match in DATE.finditer(data):
        masked[match.start():match.end()] = b"#" * (match.end() - match.start())
        dates += 1
    derived = []
    if data[:4] == MACHO_64 and len(data) >= 32:
        ncmds, sizeofcmds = struct.unpack_from("<II", data, 16)
        offset, end = 32, 32 + sizeofcmds
        segments, uuids, signatures = [], [], []
        for _ in range(ncmds):
            if offset + 8 > min(end, len(data)):
                break
            cmd, cmdsize = struct.unpack_from("<II", data, offset)
            if cmdsize < 8:
                break
            if cmd == LC_SEGMENT_64 and offset + 56 <= len(data):
                name = data[offset + 8:offset + 24].rstrip(b"\0")
                fileoff, filesize = struct.unpack_from("<QQ", data, offset + 40)
                segments.append((name, fileoff, filesize))
            elif cmd == LC_UUID and offset + 24 <= len(data):
                uuids.append(offset + 8)
            elif cmd == LC_CODE_SIGNATURE and offset + 16 <= len(data):
                signatures.append(struct.unpack_from("<II", data, offset + 8))
            offset += cmdsize
        for start in uuids:
            masked[start:start + 16] = b"\0" * 16
            derived.append("LC_UUID")
        # Checkpoint 5 (review CP3 finding 6): the signature is masked only where the
        # linker puts it: inside __LINKEDIT, ending at its end and at the end of the
        # file, after the load commands and after every other segment's file range.
        # Anything else (an edited or malformed file) is compared byte for byte.
        linkedit = [seg for seg in segments if seg[0] == b"__LINKEDIT"]
        for dataoff, datasize in signatures:
            ok = len(signatures) == 1 and len(linkedit) == 1 and datasize > 0
            if ok:
                _, fileoff, filesize = linkedit[0]
                others = [seg[1] + seg[2] for seg in segments
                          if seg[0] != b"__LINKEDIT" and seg[2]]
                ok = (fileoff <= dataoff and dataoff >= 32 + sizeofcmds
                      and dataoff + datasize == fileoff + filesize == len(data)
                      and all(dataoff >= other for other in others))
            if not ok:
                derived = ["LC_CODE_SIGNATURE not masked (range %d+%d is not at the end of "
                           "__LINKEDIT and of the file); nothing masked, compared byte for byte"
                           % (dataoff, datasize)]
                return hashlib.sha256(data).hexdigest(), dates, derived
            masked[dataoff:dataoff + datasize] = b"\0" * datasize
            derived.append("LC_CODE_SIGNATURE")
    return hashlib.sha256(bytes(masked)).hexdigest(), dates, derived


def file_digest(path, kind):
    data = Path(path).read_bytes()
    if kind == "archive":
        result = archive_digest(data)
        if result is None:
            print("%s: not a parsable archive; compared byte for byte" % path, file=sys.stderr)
            return hashlib.sha256(data).hexdigest()
        return result[0]
    digest, dates, derived = build_time_digest(data)
    if not dates:
        print("%s: no build-time date found" % path, file=sys.stderr)
    for note in derived:
        if "not masked" in note:
            print("%s: %s" % (path, note), file=sys.stderr)
    return digest


def matches(rel, pattern):
    parts, wanted = PurePosixPath(rel).parts, PurePosixPath(pattern).parts
    return len(parts) == len(wanted) and all(fnmatch.fnmatchcase(p, w) for p, w in zip(parts, wanted))


DETAIL_LINE = re.compile(r"([0-9a-f]{64})  (.+)")


def detail(args):
    root = Path(args.root)
    if not root.is_dir():
        print("no such directory: %s" % root, file=sys.stderr)
        return 2
    build_time = set(args.build_time or [])
    for line in sys.stdin.read().splitlines():
        match = DETAIL_LINE.fullmatch(line)
        if match:
            rel = match.group(2)
            kind = ("build-time" if rel in build_time else
                    "archive" if any(matches(rel, g) for g in args.archives or []) else None)
            if kind:
                try:
                    line = "%s:%s  %s" % (kind, file_digest(root / rel, kind), rel)
                except OSError as error:
                    print("cannot read %s: %s" % (rel, error), file=sys.stderr)
                    return 2
        print(line)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog="fingerprint_normalize.py", allow_abbrev=False)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("file", allow_abbrev=False)
    p.add_argument("--kind", choices=("archive", "build-time"), required=True)
    p.add_argument("path")
    p = sub.add_parser("detail", allow_abbrev=False)
    p.add_argument("--root", required=True)
    p.add_argument("--archives", action="append")
    p.add_argument("--build-time", action="append")
    args = parser.parse_args(argv)
    if args.cmd == "file":
        try:
            print("%s  %s" % (file_digest(args.path, args.kind), args.path))
        except OSError as error:
            print("cannot read %s: %s" % (args.path, error), file=sys.stderr)
            return 2
        return 0
    return detail(args)


if __name__ == "__main__":
    sys.exit(main())
