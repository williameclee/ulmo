"""Offline ESM archive tests: python3 -B -m pytest unittests/test_esm_download.py."""

import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tarfile

import pytest

MODULE = Path(__file__).resolve().parents[1] / "data" / "download_esm_epoch.py"
SPEC = importlib.util.spec_from_file_location("download_esm_epoch", MODULE)
downloader = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(downloader)


def make_archive(component, hours):
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as archive:
        for hour in hours:
            name = f"mtmshc_{component}_20030101_{hour:02d}.gfc"
            content = f"modelname {name[:-4]}\nend_of_head\ngfc 0 0 {hour} 0\n".encode()
            member = tarfile.TarInfo(f"./01/{name}")
            member.size = len(content)
            archive.addfile(member, io.BytesIO(content))
    return output.getvalue()


@pytest.fixture
def grace_file(tmp_path):
    grace = tmp_path / "GSM-fixture"
    grace.write_text("GRCOF2 0 0 1 0 0 0 20030101.0000 20030101.1200 nnnn\n")
    return grace


@pytest.fixture
def run_download(tmp_path, grace_file, monkeypatch):
    monkeypatch.setenv("IFILES", str(tmp_path))
    def run(component, hours, source=None):
        arguments = [str(MODULE), "--grace-file", str(source or grace_file), "--components", component]
        monkeypatch.setattr(sys, "argv", arguments)
        monkeypatch.setattr(downloader.urllib.request, "urlopen",
                            lambda *args, **kwargs: io.BytesIO(make_archive(component, hours)))
        downloader.main()
    return run


def test_unsorted_members_and_exclusive_end(tmp_path, run_download):
    run_download("H", [12, 6, 0])
    manifest = json.loads((tmp_path / "ESM/ESAESM2/manifests/GSM-fixture.json").read_text())
    entries = manifest["components"]["H"]
    assert [entry["timestamp"] for entry in entries] == [
        "2003-01-01T00:00:00", "2003-01-01T06:00:00"]
    assert len(list((tmp_path / "ESM/ESAESM2/raw/H").glob("*.gfc"))) == 2
    for entry in entries:
        data = (tmp_path / "ESM/ESAESM2" / entry["file"]).read_bytes()
        assert entry["sha256"] == hashlib.sha256(data).hexdigest()


@pytest.mark.parametrize("hours, error", [([0], "Missing samples"),
                                        ([0, 0, 6], "Duplicate timestamp")])
def test_missing_and_duplicate_samples_fail(run_download, hours, error):
    with pytest.raises(ValueError, match=error):
        run_download("H", hours)


def test_add_component_preserves_manifest(tmp_path, run_download):
    run_download("H", [6, 0])
    run_download("O", [0, 6])
    manifest = json.loads((tmp_path / "ESM/ESAESM2/manifests/GSM-fixture.json").read_text())
    assert set(manifest["components"]) == {"H", "O"}


def test_changed_source_cannot_reuse_manifest(grace_file, run_download):
    run_download("H", [0, 6])
    grace_file.write_text("GRCOF2 0 0 1 0 0 0 20030101.0000 20030101.0600 nnnn\n")
    with pytest.raises(ValueError, match="different GRACE"):
        run_download("O", [0])


def test_modified_existing_sample_is_preserved(tmp_path, run_download):
    run_download("H", [0, 6])
    file = next((tmp_path / "ESM/ESAESM2/raw/H").glob("*.gfc"))
    file.write_text("changed content")
    with pytest.raises(ValueError, match="Existing data differ"):
        run_download("H", [0, 6])
    assert file.read_text() == "changed content"


def test_overlapping_epochs_share_raw_samples(tmp_path, run_download):
    run_download("H", [0, 6])
    sample = tmp_path / "ESM/ESAESM2/raw/H/mtmshc_H_20030101_06.gfc"
    modified = sample.stat().st_mtime_ns
    other = tmp_path / "GSM-other"
    other.write_text("GRCOF2 0 0 1 0 0 0 20030101.0600 20030101.1200 nnnn\n")
    run_download("H", [6], source=other)
    assert sample.stat().st_mtime_ns == modified
    assert len(list(sample.parent.glob("*.gfc"))) == 2
    manifests = tmp_path / "ESM/ESAESM2/manifests"
    assert len(list(manifests.glob("*.json"))) == 2
    entries = json.loads((manifests / "GSM-other.json").read_text())["components"]["H"]
    assert [entry["timestamp"] for entry in entries] == ["2003-01-01T06:00:00"]


def test_missing_data_root_fails(run_download, monkeypatch, capsys):
    monkeypatch.delenv("IFILES")
    with pytest.raises(SystemExit) as error:
        run_download("H", [0, 6])
    assert error.value.code == 2
    assert "Set IFILES" in capsys.readouterr().err
