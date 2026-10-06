"""Controls for the claims that records_history.py makes (SPEC_A section 3).

The tool is only ever run through subprocess (-I -B); nothing here imports
it. The fixture records directory is built under a temporary directory with
fixed modification times, so every expected row is deterministic.
"""

import datetime
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

TOOL = Path(__file__).resolve().parent.parent / "payload" / "tools" / "records_history.py"
HEADER = "date | job | state | outcome | record | ticket"
FOOTER = re.compile(r"^(\d+) entries, (\d+) unreadable or damaged$")
# Local noon: the date is the same whether the tool renders local or UTC time
# (for zones within twelve hours of UTC).
STAMP = datetime.datetime(2026, 9, 15, 12, 0, 0).timestamp()
LONG_LINE = "Gamma outcome " + "x" * 140
NAMES = ("2026-10-01-alpha", "2026-10-02-beta", "2026-10-03-gamma.md", "2026-10-04-delta",
         "2026-10-05-epsilon", "2026-10-06-zeta", "2026-10-07-empty", "2026-10-08-held",
         "2026-10-09-retired")
TITLES = {"2026-10-02-beta": "Beta record", "2026-10-03-gamma.md": "Gamma flat record",
          "2026-10-04-delta": "Delta record", "2026-10-05-epsilon": "Epsilon record",
          "2026-10-06-zeta": "Zeta record"}
LAST_LINES = {"2026-10-02-beta": "Beta final outcome line.",
              "2026-10-04-delta": "Delta fallback line.",
              "2026-10-05-epsilon": "Epsilon fallback line.",
              "2026-10-06-zeta": "Zeta fallback line."}
IMPORT = re.compile(r"^\s*(?:import\s+(?P<modules>[^#]+)|from\s+(?P<package>[\w.]+)\s+import\b)")


def summary(job, state, outcome, ticket="-"):
    return (f"job: {job}\nproject: phenix\ntitle: {job} title\nstate: {state}\n"
            f"updated: 2026-10-01T10:00Z\noutcome: {outcome}\nrecord: RECORD.md\n"
            f"ticket: {ticket}\nprocedure: r10 rev16\ninstallation: /opt/phenix\nresources: -\n")


def listing(root):
    """Names, modification times and sizes of everything under root."""
    entries = [("", root.stat().st_mtime_ns, 0)]
    for path in sorted(root.rglob("*")):
        info = path.lstat()
        entries.append((path.relative_to(root).as_posix(), info.st_mtime_ns, info.st_size))
    return entries


def imported_modules(source):
    found = set()
    for line in source.splitlines():
        match = IMPORT.match(line)
        if not match:
            continue
        if match["package"]:
            found.add(match["package"].split(".")[0])
        else:
            for item in match["modules"].split(","):
                found.add(item.strip().split()[0].split(".")[0])
    return found


class RecordsHistoryChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.records = Path(self.temp.name) / "records"
        self.records.mkdir()
        # (1) a valid summary whose state/updated/outcome were appended later
        alpha = self.directory("2026-10-01-alpha")
        (alpha / "JOB_SUMMARY.txt").write_text(
            summary("2026-10-01-alpha", "planned", "PLAN drafted",
                    ticket="/Users/dev/Downloads/alpha-ticket.md")
            + "# appended at integration\nstate: integrated\nupdated: 2026-10-02T11:00Z\n"
              "outcome: RESULT integrated by the Developer\n")
        self.record(alpha / "RECORD.md", "Alpha record", "Alpha last line.")
        # (2) RECORD.md only
        self.record(self.directory("2026-10-02-beta") / "RECORD.md", "Beta record",
                    "Beta final outcome line.")
        # (3) a flat record whose last line is longer than 120 characters
        self.record(self.records / "2026-10-03-gamma.md", "Gamma flat record", LONG_LINE)
        # (4) damaged summary: a required key is missing; RECORD.md fallback
        delta = self.directory("2026-10-04-delta")
        (delta / "JOB_SUMMARY.txt").write_text(
            summary("2026-10-04-delta", "building", "WORK in progress").replace(
                "record: RECORD.md\n", ""))
        self.record(delta / "RECORD.md", "Delta record", "Delta fallback line.")
        # (5) damaged summary: not UTF-8; RECORD.md fallback
        epsilon = self.directory("2026-10-05-epsilon")
        (epsilon / "JOB_SUMMARY.txt").write_bytes(
            summary("2026-10-05-epsilon", "result", "RESULT ready").encode("utf-8")
            + b"outcome: \xff\xfe broken\n")
        self.record(epsilon / "RECORD.md", "Epsilon record", "Epsilon fallback line.")
        # (6) damaged summary: a line without ': '; RECORD.md fallback
        zeta = self.directory("2026-10-06-zeta")
        (zeta / "JOB_SUMMARY.txt").write_text(
            summary("2026-10-06-zeta", "approved", "PLAN approved") + "state stopped\n")
        self.record(zeta / "RECORD.md", "Zeta record", "Zeta fallback line.")
        # (7) an empty directory
        self.directory("2026-10-07-empty")
        # (8) a held job and a retired job
        (self.directory("2026-10-08-held") / "JOB_SUMMARY.txt").write_text(
            summary("2026-10-08-held", "held", "HELD by the Developer pending review"))
        (self.directory("2026-10-09-retired") / "JOB_SUMMARY.txt").write_text(
            summary("2026-10-09-retired", "retired", "RETIRED after publication"))
        # files that are not records
        (self.records / "INTRO_SHOWN").write_text("2026-10-01\n")
        (self.records / "notes.txt").write_text("not a record\n")

    def directory(self, name):
        path = self.records / name
        path.mkdir()
        return path

    def record(self, path, title, last_line):
        path.write_text(f"# {title}\n\nSome body text.\n\n{last_line}\n\n")
        os.utime(path, (STAMP, STAMP))

    def run_tool(self, *arguments):
        return subprocess.run([sys.executable, "-I", "-B", str(TOOL), *arguments],
                              cwd=self.temp.name, text=True, capture_output=True)

    def rows(self, stdout):
        """The row lines: everything but the header, the footer and `showing`."""
        lines = stdout.splitlines()
        self.assertEqual(lines[0], HEADER)
        return [line for line in lines[1:]
                if not FOOTER.match(line) and not line.startswith("showing ")]

    def row(self, stdout, name):
        matches = [line for line in self.rows(stdout) if name in line]
        self.assertEqual(len(matches), 1, (name, matches))
        return matches[0]

    def fields(self, stdout, name):
        fields = self.row(stdout, name).split(" | ")
        self.assertEqual(len(fields), 6, fields)
        self.assertTrue(fields[4] == name or fields[4].startswith(name + "/"), fields)
        return fields

    def test_header_rows_in_name_order_and_counts(self):
        """SPEC_A 3: header line, one row per entry sorted by name, final
        counts line, exit 0 despite unreadable rows."""
        result = self.run_tool(str(self.records))
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.splitlines()
        self.assertEqual(lines[0], HEADER)
        self.assertEqual(lines[-1], "9 entries, 4 unreadable or damaged")
        rows = self.rows(result.stdout)
        self.assertEqual(len(rows), 9, rows)
        for index, name in enumerate(NAMES):
            self.assertIn(name, rows[index])

    def test_summary_rows_show_the_last_repeated_values(self):
        """SPEC_A 3: when state/updated/outcome repeat, the LAST values are shown."""
        stdout = self.run_tool(str(self.records)).stdout
        date, job, state, outcome, record, ticket = self.fields(stdout, "2026-10-01-alpha")
        self.assertTrue(date.startswith("2026-10-02"), date)
        self.assertEqual(job, "2026-10-01-alpha")
        self.assertEqual(state, "integrated")
        self.assertEqual(outcome, "RESULT integrated by the Developer")
        self.assertEqual(ticket, "/Users/dev/Downloads/alpha-ticket.md")
        self.assertNotIn("planned", self.row(stdout, "2026-10-01-alpha"))
        self.assertNotIn("PLAN drafted", self.row(stdout, "2026-10-01-alpha"))

    def test_record_only_directory_falls_back_to_record_md(self):
        """SPEC_A 3: no summary: title from the # heading, state `no summary`,
        outcome = last non-blank line, date = mtime, ticket `-`."""
        stdout = self.run_tool(str(self.records)).stdout
        date, job, state, outcome, record, ticket = self.fields(stdout, "2026-10-02-beta")
        self.assertEqual(date, "2026-09-15")
        self.assertEqual(state, "no summary")
        self.assertEqual(outcome, "Beta final outcome line.")
        self.assertEqual(ticket, "-")
        self.assertIn("Beta record", self.row(stdout, "2026-10-02-beta"))

    def test_flat_record_is_an_entry_with_a_truncated_outcome(self):
        """SPEC_A 3: a flat <id>.md is read as RECORD.md; the outcome is cut to
        120 characters; record is the file name."""
        stdout = self.run_tool(str(self.records)).stdout
        date, job, state, outcome, record, ticket = self.fields(stdout, "2026-10-03-gamma.md")
        self.assertEqual(record, "2026-10-03-gamma.md")
        self.assertEqual(date, "2026-09-15")
        self.assertEqual(state, "no summary")
        self.assertEqual(outcome, LONG_LINE[:120])
        self.assertEqual(ticket, "-")
        self.assertIn("Gamma flat record", self.row(stdout, "2026-10-03-gamma.md"))

    def test_damaged_summaries_say_so_and_fall_back(self):
        """SPEC_A 3: a missing required key, a non-UTF-8 file or a line without
        `: ` gives `damaged summary (<reason>)` with the RECORD.md values."""
        stdout = self.run_tool(str(self.records)).stdout
        for name in ("2026-10-04-delta", "2026-10-05-epsilon", "2026-10-06-zeta"):
            with self.subTest(entry=name):
                row = self.row(stdout, name)
                self.assertRegex(row, r"damaged summary \([^)]+\)")
                self.assertIn(TITLES[name], row)
                self.assertIn(LAST_LINES[name], row)
                self.assertIn("2026-09-15", row)
        # the damaged summaries' own values are not shown
        self.assertNotIn("WORK in progress", stdout)
        self.assertNotIn("PLAN approved", stdout)
        self.assertNotIn("RESULT ready", stdout)

    def test_empty_directory_is_unreadable_and_other_files_are_ignored(self):
        """SPEC_A 3: neither file -> `unreadable (no RECORD.md or
        JOB_SUMMARY.txt)`; INTRO_SHOWN and other non-.md files are no entries."""
        stdout = self.run_tool(str(self.records)).stdout
        self.assertIn("unreadable (no RECORD.md or JOB_SUMMARY.txt)",
                      self.row(stdout, "2026-10-07-empty"))
        self.assertNotIn("INTRO_SHOWN", stdout)
        self.assertNotIn("notes.txt", stdout)

    def test_held_and_retired_jobs_are_listed(self):
        """SPEC_A 3: held and retired jobs stay in the listing with their state."""
        stdout = self.run_tool(str(self.records)).stdout
        held = self.fields(stdout, "2026-10-08-held")
        self.assertEqual(held[2], "held")
        self.assertEqual(held[3], "HELD by the Developer pending review")
        self.assertEqual(held[5], "-")
        retired = self.fields(stdout, "2026-10-09-retired")
        self.assertEqual(retired[2], "retired")
        self.assertEqual(retired[3], "RETIRED after publication")

    def test_limit_keeps_the_last_entries_by_name(self):
        """SPEC_A 3: --limit N keeps the last N by name and prints showing N of M."""
        result = self.run_tool(str(self.records), "--limit", "2")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = self.rows(result.stdout)
        self.assertEqual(len(rows), 2, rows)
        self.assertIn("2026-10-08-held", rows[0])
        self.assertIn("2026-10-09-retired", rows[1])
        self.assertIn("showing 2 of 9", result.stdout)
        self.assertRegex(result.stdout, r"(?m)^9 entries, \d+ unreadable or damaged$")
        self.assertNotIn("2026-10-01-alpha", result.stdout)

    def test_missing_directory_exits_2(self):
        """SPEC_A 3: a missing path or a plain file is `records directory not found`."""
        for path in (self.records / "nowhere", self.records / "INTRO_SHOWN"):
            with self.subTest(path=path.name):
                result = self.run_tool(str(path))
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertIn("ERROR:", result.stderr)
                self.assertIn("records directory not found", result.stderr)
                self.assertEqual(result.stdout, "")

    def test_tool_writes_nothing(self):
        """SPEC_A 3: the records directory and the tool's own directory are
        identical (names, mtimes, sizes) before and after every invocation."""
        before = (listing(self.records), listing(TOOL.parent))
        self.assertEqual(self.run_tool(str(self.records)).returncode, 0)
        self.assertEqual(self.run_tool(str(self.records), "--limit", "2").returncode, 0)
        self.assertEqual(self.run_tool(str(self.records / "nowhere")).returncode, 2)
        self.assertEqual((listing(self.records), listing(TOOL.parent)), before)

    def test_imports_only_the_allowed_standard_modules(self):
        """SPEC_A 3: the source imports nothing outside argparse, datetime, os,
        re, sys and pathlib (a regex over import/from lines)."""
        found = imported_modules(TOOL.read_text(encoding="utf-8"))
        self.assertTrue(found, "no import lines found")
        self.assertLessEqual(found, {"argparse", "datetime", "os", "re", "sys", "pathlib"})

    # ---- SPEC_A 6: symbolic links and --limit -------------------------------

    def footer(self, stdout):
        """(entries, unreadable or damaged) from the last output line."""
        match = FOOTER.match(stdout.splitlines()[-1])
        self.assertIsNotNone(match, stdout)
        return int(match.group(1)), int(match.group(2))

    def test_symbolic_link_entries_are_skipped_and_never_followed(self):
        """SPEC_A 6 (3): a symbolic link in the records directory, whether to a
        directory holding a valid JOB_SUMMARY.txt or to a flat .md file, is its
        own row with state `skipped (symbolic link)`, counted as unreadable or
        damaged; nothing from the link targets appears in the output."""
        before = self.run_tool(str(self.records))
        self.assertEqual(before.returncode, 0, before.stderr)
        self.assertEqual(self.footer(before.stdout), (9, 4))
        outside = Path(self.temp.name) / "outside"
        target_dir = outside / "2026-10-10-linked-target"
        target_dir.mkdir(parents=True)
        (target_dir / "JOB_SUMMARY.txt").write_text(
            summary("2026-10-10-linked-target", "published", "SENTINEL-LINKED-DIR-OUTCOME",
                    ticket="/sentinel/linked-dir-ticket.md"))
        self.record(target_dir / "RECORD.md", "SENTINEL-LINKED-DIR-TITLE",
                    "SENTINEL-LINKED-DIR-LAST-LINE")
        self.record(outside / "linked-flat-target.md", "SENTINEL-LINKED-FLAT-TITLE",
                    "SENTINEL-LINKED-FLAT-LAST-LINE")
        link_dir = self.records / "2026-10-10-linked-dir"
        link_file = self.records / "2026-10-11-linked-file.md"
        link_dir.symlink_to(target_dir, target_is_directory=True)
        link_file.symlink_to(outside / "linked-flat-target.md")
        # the links resolve: a tool that followed them would see valid records
        self.assertTrue((link_dir / "JOB_SUMMARY.txt").is_file())
        self.assertTrue(link_file.is_file())
        result = self.run_tool(str(self.records))
        self.assertEqual(result.returncode, 0, result.stderr)
        output = result.stdout + result.stderr
        self.assertNotIn("SENTINEL", output)
        self.assertNotIn("/sentinel/", output)
        self.assertNotIn("published", output)
        for name in (link_dir.name, link_file.name):
            with self.subTest(entry=name):
                self.assertIn("skipped (symbolic link)", self.row(result.stdout, name))
        self.assertEqual(len(self.rows(result.stdout)), 11)
        self.assertEqual(result.stdout.splitlines()[-1], "11 entries, 6 unreadable or damaged")

    def test_symlinked_summary_or_record_inside_an_entry_is_absent(self):
        """SPEC_A 6 (3): a symlinked JOB_SUMMARY.txt is treated as absent (the
        row falls back to the entry's own RECORD.md with state `no summary`,
        and the target's values do not appear); a symlinked RECORD.md with no
        summary leaves the entry unreadable."""
        outside = Path(self.temp.name) / "outside"
        outside.mkdir()
        (outside / "eta_summary.txt").write_text(
            summary("2026-10-12-eta", "published", "SENTINEL-ETA-OUTCOME",
                    ticket="/sentinel/eta-ticket.md"))
        eta = self.directory("2026-10-12-eta")
        self.record(eta / "RECORD.md", "Eta record", "Eta fallback line.")
        (eta / "JOB_SUMMARY.txt").symlink_to(outside / "eta_summary.txt")
        self.record(outside / "theta_record.md", "SENTINEL-THETA-TITLE",
                    "SENTINEL-THETA-LAST-LINE")
        theta = self.directory("2026-10-13-theta")
        (theta / "RECORD.md").symlink_to(outside / "theta_record.md")
        self.assertTrue((eta / "JOB_SUMMARY.txt").is_file())
        self.assertTrue((theta / "RECORD.md").is_file())
        result = self.run_tool(str(self.records))
        self.assertEqual(result.returncode, 0, result.stderr)
        output = result.stdout + result.stderr
        self.assertNotIn("SENTINEL", output)
        self.assertNotIn("/sentinel/", output)
        self.assertNotIn("published", output)
        date, job, state, outcome, record, ticket = self.fields(result.stdout, "2026-10-12-eta")
        self.assertEqual(date, "2026-09-15")
        self.assertEqual(state, "no summary")
        self.assertEqual(outcome, "Eta fallback line.")
        self.assertEqual(ticket, "-")
        self.assertIn("Eta record", self.row(result.stdout, "2026-10-12-eta"))
        self.assertNotIn("damaged summary", self.row(result.stdout, "2026-10-12-eta"))
        self.assertIn("unreadable (no RECORD.md or JOB_SUMMARY.txt)",
                      self.row(result.stdout, "2026-10-13-theta"))
        self.assertEqual(len(self.rows(result.stdout)), 11)
        self.assertEqual(result.stdout.splitlines()[-1], "11 entries, 5 unreadable or damaged")

    def test_limit_must_be_a_positive_integer(self):
        """SPEC_A 6 (3): --limit 0 and --limit -1 exit 2 with `--limit must be
        a positive integer` on stderr and no rows on stdout; --limit 1 still
        works."""
        for value in ("0", "-1"):
            with self.subTest(limit=value):
                result = self.run_tool(str(self.records), "--limit", value)
                self.assertEqual(result.returncode, 2, result.stdout)
                self.assertIn("ERROR:", result.stderr)
                self.assertIn("--limit must be a positive integer", result.stderr)
                for name in NAMES:
                    self.assertNotIn(name, result.stdout)
                self.assertNotIn("unreadable or damaged", result.stdout)
                self.assertNotIn("showing ", result.stdout)
        result = self.run_tool(str(self.records), "--limit", "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(self.rows(result.stdout)), 1)
        self.assertIn("2026-10-09-retired", self.rows(result.stdout)[0])
        self.assertIn("showing 1 of 9", result.stdout)

    def test_dangling_symbolic_links_are_skipped_rows(self):
        """SPEC_A 8 (3): a symbolic link whose target does not exist, named like
        a record directory or like <id>.md, is its own row with `skipped
        (symbolic link)`, counted in the entries and in unreadable or damaged;
        the tool still exits 0."""
        missing = Path(self.temp.name) / "missing"
        link_dir = self.records / "2026-10-14-dangling"
        link_file = self.records / "2026-10-15-dangling.md"
        link_dir.symlink_to(missing / "2026-10-14-target", target_is_directory=True)
        link_file.symlink_to(missing / "target.md")
        for link in (link_dir, link_file):
            self.assertTrue(link.is_symlink())
            self.assertFalse(link.exists())
        result = self.run_tool(str(self.records))
        self.assertEqual(result.returncode, 0, result.stderr)
        for name in (link_dir.name, link_file.name):
            with self.subTest(entry=name):
                self.assertIn("skipped (symbolic link)", self.row(result.stdout, name))
        self.assertEqual(len(self.rows(result.stdout)), 11)
        self.assertEqual(result.stdout.splitlines()[-1], "11 entries, 6 unreadable or damaged")
        self.assertNotIn("Traceback", result.stderr)


if __name__ == "__main__":
    unittest.main()
