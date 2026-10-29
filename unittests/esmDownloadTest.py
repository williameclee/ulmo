"""Offline tests for ESM archive selection and download provenance.

Run with: python3 -B -m unittest discover -s unittests -p esmDownloadTest.py
"""

import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

MODULE = Path(__file__).resolve().parents[1] / "data" / "download_esm_epoch.py"
SPEC = importlib.util.spec_from_file_location("download_esm_epoch", MODULE)
downloader = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(downloader)


class EsmDownloadTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.grace = self.root / "GSM-fixture"
        self.grace.write_text("GRCOF2 0 0 1 0 0 0 20030101.0000 20030101.1200 nnnn\n")

    def archive(self, component, hours):
        output = io.BytesIO()
        with tarfile.open(fileobj=output, mode="w:gz") as archive:
            for hour in hours:
                name = f"mtmshc_{component}_20030101_{hour:02d}.gfc"
                content = f"modelname {name[:-4]}\nend_of_head\ngfc 0 0 {hour} 0\n".encode()
                member = tarfile.TarInfo(f"./01/{name}")
                member.size = len(content)
                archive.addfile(member, io.BytesIO(content))
        return output.getvalue()

    def run_download(self, component, hours):
        arguments = [str(MODULE), "--grace-file", str(self.grace), "--output",
                     str(self.root / "data"), "--components", component]
        with patch("sys.argv", arguments), patch.object(downloader.urllib.request,
                "urlopen", return_value=io.BytesIO(self.archive(component, hours))):
            downloader.main()

    def test_unsorted_members_and_exclusive_end(self):
        self.run_download("H", [12, 6, 0])
        manifest = json.loads((self.root / "data/manifest.json").read_text())
        entries = manifest["components"]["H"]
        self.assertEqual([entry["timestamp"] for entry in entries],
                         ["2003-01-01T00:00:00", "2003-01-01T06:00:00"])
        self.assertEqual(len(list((self.root / "data/raw/H").glob("*.gfc"))), 2)
        for entry in entries:
            data = (self.root / "data" / entry["file"]).read_bytes()
            self.assertEqual(entry["sha256"], hashlib.sha256(data).hexdigest())

    def test_missing_and_duplicate_samples_fail(self):
        with self.assertRaisesRegex(ValueError, "Missing samples"):
            self.run_download("H", [0])
        with self.assertRaisesRegex(ValueError, "Duplicate timestamp"):
            self.run_download("H", [0, 0, 6])

    def test_add_component_preserves_manifest(self):
        self.run_download("H", [6, 0])
        self.run_download("O", [0, 6])
        manifest = json.loads((self.root / "data/manifest.json").read_text())
        self.assertEqual(set(manifest["components"]), {"H", "O"})

    def test_different_interval_cannot_reuse_destination(self):
        self.run_download("H", [0, 6])
        self.grace.write_text("GRCOF2 0 0 1 0 0 0 20030101.0000 20030101.0600 nnnn\n")
        with self.assertRaisesRegex(ValueError, "different GRACE"):
            self.run_download("O", [0])

    def test_modified_existing_sample_is_preserved(self):
        self.run_download("H", [0, 6])
        file = next((self.root / "data/raw/H").glob("*.gfc"))
        file.write_text("changed content")
        with self.assertRaisesRegex(ValueError, "Existing data differ"):
            self.run_download("H", [0, 6])
        self.assertEqual(file.read_text(), "changed content")


if __name__ == "__main__":
    unittest.main()
