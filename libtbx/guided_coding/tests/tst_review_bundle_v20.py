"""Transport checks for the adopted Outside Reviewer bundle layout."""

import contextlib
import io
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).parent.parent / "payload" / "tools"))
import review_bundle
import screen_check


class ReviewBundleChecks(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        self.packet, self.companions = root / "packet", root / "companions"
        self.packet.mkdir()
        self.companions.mkdir()
        (self.packet / "README_FIRST.md").write_text("Problem and reading order\n")
        (self.packet / "PROOF_SUMMARY.md").write_text("Evidence, then provisional verdict\n")
        (self.packet / "RUN_LOG.txt").write_text("command exited 0\n")
        (self.companions / "APPROVALS.md").write_text("Synthetic approval; not a real gate\n")
        (self.companions / "CHECKER.md").write_text("Synthetic reading\n")
        self.out = root / "bundle.tgz"
        with contextlib.redirect_stdout(io.StringIO()):
            screen_check.freeze(self.packet)

    def build(self):
        with contextlib.redirect_stdout(io.StringIO()):
            review_bundle.make_bundle(self.packet, self.companions, self.out)

    def test_bundle_has_two_integrity_domains(self):
        self.build()
        with tarfile.open(self.out, "r:gz") as outer:
            files = {entry.name: outer.extractfile(entry).read() for entry in outer if entry.isfile()}
        self.assertEqual(set(files), {"PACKET.tgz", "APPROVALS.md", "CHECKER.md", "SHA256SUMS"})
        checksums = files["SHA256SUMS"].decode().splitlines()
        self.assertNotIn("SHA256SUMS", " ".join(checksums))
        for name in ("PACKET.tgz", "APPROVALS.md", "CHECKER.md"):
            self.assertIn(f"{screen_check.digest(files[name])}  {name}", checksums)
        with tarfile.open(fileobj=io.BytesIO(files["PACKET.tgz"]), mode="r:gz") as inner:
            manifest = inner.extractfile("MANIFEST.sha256").read()
        self.assertEqual(screen_check.digest(manifest), screen_check.verify(self.packet))

    def test_companions_nested_inside_packet_refused_before_write(self):
        nested = self.packet / "companions"
        nested.mkdir()
        for name in ("APPROVALS.md", "CHECKER.md"):
            (nested / name).write_bytes((self.companions / name).read_bytes())
        (self.packet / "MANIFEST.sha256").unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            screen_check.freeze(self.packet)
        with self.assertRaisesRegex(ValueError, "outside the frozen packet"):
            review_bundle.make_bundle(self.packet, nested, self.out)
        self.assertFalse(self.out.exists())
        screen_check.verify(self.packet)

        # Also reject duplicated companions when the passed directory is
        # outside but the already frozen packet contains copies.
        with self.assertRaisesRegex(ValueError, "companion file present inside"):
            self.build()
        self.assertFalse(self.out.exists())
        screen_check.verify(self.packet)

    def test_missing_summary_changed_packet_and_output_collision(self):
        self.build()
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.build()
        (self.packet / "RUN_LOG.txt").write_text("modified after freeze\n")
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            review_bundle.make_bundle(self.packet, self.companions, self.out.parent / "other.tgz")

    def test_output_inside_packet_refused_without_damaging_freeze(self):
        with self.assertRaisesRegex(ValueError, "beside the frozen packet"):
            review_bundle.make_bundle(self.packet, self.companions, self.packet / "stray.tgz")
        self.assertFalse((self.packet / "stray.tgz").exists())
        screen_check.verify(self.packet)

    def test_symlinked_parent_into_packet_refused_before_write(self):
        alias = self.packet.parent / "alias"
        alias.symlink_to(self.packet, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "beside the frozen packet"):
            review_bundle.make_bundle(self.packet, self.companions, alias / "stray.tgz")
        self.assertFalse((self.packet / "stray.tgz").exists())
        screen_check.verify(self.packet)

    def test_symlinked_packet_parent_can_receive_bundle(self):
        alias = self.packet.parent / "parent-alias"
        alias.symlink_to(self.packet.parent, target_is_directory=True)
        with contextlib.redirect_stdout(io.StringIO()):
            review_bundle.make_bundle(self.packet, self.companions, alias / "bundle.tgz")
        self.assertTrue((self.packet.parent / "bundle.tgz").is_file())
        screen_check.verify(self.packet)

    def test_case_variant_of_packet_refused_on_case_insensitive_volume(self):
        variant = self.packet.with_name(self.packet.name.upper())
        if not variant.exists() or not variant.samefile(self.packet):
            self.skipTest("requires a case-insensitive filesystem")
        with self.assertRaisesRegex(ValueError, "beside the frozen packet"):
            review_bundle.make_bundle(self.packet, self.companions,
                                      variant / "case-variant.tgz")
        self.assertFalse((self.packet / "case-variant.tgz").exists())
        screen_check.verify(self.packet)

    def test_replaced_output_alias_cannot_redirect_open_to_packet(self):
        alias = self.packet.parent / "alias"
        alias.symlink_to(self.packet.parent, target_is_directory=True)
        original_add = review_bundle.add_bytes
        replaced = False

        def replace_alias(archive, name, data):
            nonlocal replaced
            if not replaced:
                alias.unlink()
                alias.symlink_to(self.packet, target_is_directory=True)
                replaced = True
            return original_add(archive, name, data)

        with mock.patch.object(review_bundle, "add_bytes", side_effect=replace_alias):
            with contextlib.redirect_stdout(io.StringIO()):
                review_bundle.make_bundle(self.packet, self.companions,
                                          alias / "race.tgz")
        self.assertTrue(replaced)
        self.assertTrue((self.packet.parent / "race.tgz").is_file())
        self.assertFalse((self.packet / "race.tgz").exists())
        screen_check.verify(self.packet)

    def test_movable_sibling_directory_refused_before_move_or_write(self):
        separate = self.packet.parent / "separate"
        separate.mkdir()
        invoked = False
        original_add = review_bundle.add_bytes

        def move_directory(archive, name, data):
            nonlocal invoked
            invoked = True
            separate.rename(self.packet / "moved")
            return original_add(archive, name, data)

        with mock.patch.object(review_bundle, "add_bytes", side_effect=move_directory):
            with self.assertRaisesRegex(ValueError, "beside the frozen packet"):
                review_bundle.make_bundle(self.packet, self.companions,
                                          separate / "output.tgz")
        self.assertFalse(invoked)
        self.assertFalse((self.packet / "moved/output.tgz").exists())
        screen_check.verify(self.packet)

    def test_packet_parent_cannot_be_moved_under_its_packet(self):
        original_add = review_bundle.add_bytes
        attempted = False

        def attempt_move(archive, name, data):
            nonlocal attempted
            if not attempted:
                attempted = True
                with self.assertRaises(OSError):
                    self.packet.parent.rename(self.packet / "moved-parent")
            return original_add(archive, name, data)

        with mock.patch.object(review_bundle, "add_bytes", side_effect=attempt_move):
            self.build()
        self.assertTrue(attempted)
        self.assertTrue(self.out.is_file())
        self.assertFalse((self.packet / "moved-parent").exists())
        screen_check.verify(self.packet)

    def test_packet_moved_before_output_creation_refused(self):
        original_add = review_bundle.add_bytes
        moved_packet = self.packet.parent / "relocated-packet"
        moved = False

        def move_packet(archive, name, data):
            nonlocal moved
            if not moved:
                self.packet.rename(moved_packet)
                moved = True
            return original_add(archive, name, data)

        with mock.patch.object(review_bundle, "add_bytes", side_effect=move_packet):
            with self.assertRaises(OSError):
                self.build()
        self.assertTrue(moved)
        self.assertFalse(self.out.exists())
        screen_check.verify(moved_packet)

    def test_nested_output_directory_inside_packet_refused(self):
        (self.packet / "nested").mkdir()
        (self.packet / "MANIFEST.sha256").unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            screen_check.freeze(self.packet)
        with self.assertRaisesRegex(ValueError, "beside the frozen packet"):
            review_bundle.make_bundle(self.packet, self.companions,
                                      self.packet / "nested" / "stray.tgz")
        self.assertFalse((self.packet / "nested" / "stray.tgz").exists())
        screen_check.verify(self.packet)

    def test_requires_summary_and_refuses_linked_companion(self):
        (self.packet / "MANIFEST.sha256").unlink()
        (self.packet / "PROOF_SUMMARY.md").unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            screen_check.freeze(self.packet)
        with self.assertRaisesRegex(ValueError, "needs README_FIRST"):
            self.build()
        (self.packet / "MANIFEST.sha256").unlink()
        (self.packet / "PROOF_SUMMARY.md").write_text("Rebuilt synthetic summary\n")
        with contextlib.redirect_stdout(io.StringIO()):
            screen_check.freeze(self.packet)
        (self.companions / "link").symlink_to(self.companions / "CHECKER.md")
        with self.assertRaisesRegex(ValueError, "symbolic link"):
            self.build()


if __name__ == "__main__":
    unittest.main()
