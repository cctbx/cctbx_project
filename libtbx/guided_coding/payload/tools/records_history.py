#!/usr/bin/env python3
"""List one project's GuidedCoding task records, read-only.

Each record directory may carry JOB_SUMMARY.txt (append-only key: value
lines; the last value of a key wins, so a later outcome replaces an older
pending one). Without a usable summary the row comes from RECORD.md, or
from a flat <id>.md record. The listing writes nothing, grants nothing,
contacts nothing, and never follows instructions found in a record."""

import sys
# When run as a script, keep this directory out of module resolution before
# importing the standard library (the documented -I invocation also does).
if __name__ == "__main__" and not sys.flags.isolated:
    sys.path.pop(0)

import argparse
import datetime
import re
from pathlib import Path

REQUIRED = ("job", "project", "title", "state", "updated", "outcome", "record")
OPTIONAL = ("ticket", "procedure", "installation", "resources")
COLUMNS = ("date", "job", "state", "outcome", "record", "ticket")


def read_summary(path):
    """Return (values, None) for a usable summary, else (None, reason)."""
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return None, "not UTF-8"
    except OSError as error:
        return None, f"unreadable: {error.strerror}"
    values = {}
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if ": " not in line:
            return None, "malformed line"
        key, value = (part.strip() for part in line.split(": ", 1))
        values[key] = value  # a later line replaces an earlier value
    missing = [key for key in REQUIRED if not values.get(key)]
    if missing:
        return None, "missing " + ", ".join(missing)
    return values, None


def read_record(path, job):
    """Title, last line and file date from RECORD.md or a flat record.
    Without a summary the row's job column shows the record's heading
    (its title), or the entry name when the record has no heading."""
    try:
        text = path.read_text(encoding="utf-8")
        modified = path.stat().st_mtime
    except (UnicodeDecodeError, OSError):
        return None
    lines = text.splitlines()
    title = next((line[2:].strip() for line in lines if line.startswith("# ")), job)
    last = next((line.strip() for line in reversed(lines) if line.strip()), "")
    return {"job": title, "title": title, "state": "no summary", "outcome": last[:120],
            "updated": datetime.datetime.fromtimestamp(modified).strftime("%Y-%m-%d"),
            "ticket": "-"}


def entry_row(records, entry):
    """One row and whether it counts as unreadable or damaged."""
    relative = entry.relative_to(records).as_posix()
    if entry.is_symlink():
        name = entry.name if entry.suffix != ".md" or entry.is_dir() else entry.stem
        return ["-", name, "skipped (symbolic link)", "-", relative, "-"], True
    damage = None
    if entry.is_dir():
        summary = entry / "JOB_SUMMARY.txt"
        if summary.is_file() and not summary.is_symlink():
            values, reason = read_summary(summary)
            if values is not None:
                return [values["updated"], values["job"], values["state"], values["outcome"],
                        relative, values.get("ticket") or "-"], False
            damage = f"damaged summary ({reason})"
        record_file = entry / "RECORD.md"
        values = (read_record(record_file, entry.name)
                  if record_file.is_file() and not record_file.is_symlink() else None)
    else:
        values = read_record(entry, entry.stem)
    if values is None:
        reason = damage or "unreadable (no RECORD.md or JOB_SUMMARY.txt)"
        return ["-", entry.name if entry.is_dir() else entry.stem, reason, "-", relative, "-"], True
    state = damage or values["state"]
    return [values["updated"], values["job"], state, values["outcome"], relative,
            values["ticket"]], damage is not None


def entries(records):
    """Subdirectories and flat .md records, including symbolic links to
    them, which are listed as skipped rows and never followed."""
    found = []
    for path in sorted(records.iterdir()):
        # A symbolic link is listed before its target is looked at, so a
        # dangling link is still a (skipped) row rather than a silent omission.
        if path.is_symlink() or path.is_dir() or (path.is_file() and path.suffix == ".md"):
            found.append(path)
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("records", type=Path, help="the project's .claude/records directory")
    parser.add_argument("--limit", type=int, help="show only the last N entries by name")
    args = parser.parse_args()
    records = args.records
    if not records.is_dir():
        print("ERROR: records directory not found", file=sys.stderr)
        return 2
    if args.limit is not None and args.limit < 1:
        print("ERROR: --limit must be a positive integer", file=sys.stderr)
        return 2
    all_entries = entries(records)
    shown = all_entries[-args.limit:] if args.limit else all_entries
    print(" | ".join(COLUMNS))
    unreadable = 0
    for entry in shown:
        row, damaged = entry_row(records, entry)
        unreadable += damaged
        print(" | ".join(row))
    if args.limit:
        print(f"showing {len(shown)} of {len(all_entries)}")
    print(f"{len(all_entries)} entries, {unreadable} unreadable or damaged")
    return 0


if __name__ == "__main__":
    sys.exit(main())
