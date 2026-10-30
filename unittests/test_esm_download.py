"""Offline ESM archive tests: python3 -B -m pytest unittests/test_esm_download.py."""

import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tarfile

import pytest

MODULE = Path(__file__).resolve().parents[1] / "data" / "download_esm.py"
SPEC = importlib.util.spec_from_file_location("download_esm", MODULE)
downloader = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(downloader)


def make_archive(component, hours, date="20030101"):
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as archive:
        for hour in hours:
            name = f"mtmshc_{component}_{date}_{hour:02d}.gfc"
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


def test_parsing_does_not_read_or_download(tmp_path, monkeypatch):
    monkeypatch.setenv("IFILES", str(tmp_path))
    missing_source = tmp_path / "not-downloaded-yet"
    args = downloader.parse_arguments(["--grace-file", str(missing_source)])
    assert args.grace_file == missing_source
    assert args.components == ["H", "I", "O"]
    assert args.data_root == tmp_path
    assert list(tmp_path.iterdir()) == []


def test_programmatic_download_uses_explicit_root(tmp_path, grace_file, monkeypatch):
    monkeypatch.delenv("IFILES", raising=False)
    monkeypatch.setattr(downloader.urllib.request, "urlopen",
                        lambda *args, **kwargs: io.BytesIO(make_archive("H", [6, 0])))
    manifest = downloader.download_epoch(grace_file, ["H"], tmp_path / "explicit")
    assert manifest == tmp_path / "explicit/ESM/ESAESM2/manifests/GSM-fixture.json"
    assert len(json.loads(manifest.read_text())["components"]["H"]) == 2


def test_completed_component_survives_later_download_failure(tmp_path, grace_file, monkeypatch):
    def open_archive(url, **kwargs):
        if url.endswith("_H.tar.gz"):
            return io.BytesIO(make_archive("H", [0, 6]))
        raise OSError("connection interrupted")
    monkeypatch.setattr(downloader.urllib.request, "urlopen", open_archive)
    with pytest.raises(OSError, match="connection interrupted"):
        downloader.download_epoch(grace_file, ["H", "O"], tmp_path)
    path = tmp_path / "ESM/ESAESM2/manifests/GSM-fixture.json"
    assert set(json.loads(path.read_text())["components"]) == {"H"}


def test_epoch_crossing_year_boundary_downloads_both_archives(tmp_path, grace_file, monkeypatch):
    grace_file.write_text("GRCOF2 0 0 1 0 0 0 20021231.1800 20030101.0600 nnnn\n")
    urls = []
    def open_archive(url, **kwargs):
        urls.append(url)
        if "mtm_2002_" in url:
            return io.BytesIO(make_archive("H", [18], date="20021231"))
        return io.BytesIO(make_archive("H", [0]))
    monkeypatch.setattr(downloader.urllib.request, "urlopen", open_archive)
    path = downloader.download_epoch(grace_file, ["H"], tmp_path)
    entries = json.loads(path.read_text())["components"]["H"]
    assert [entry["timestamp"] for entry in entries] == [
        "2002-12-31T18:00:00", "2003-01-01T00:00:00"]
    assert urls == [f"{downloader.BASE}/mtmshc/mtm_2002_H.tar.gz",
                    f"{downloader.BASE}/mtmshc/mtm_2003_H.tar.gz"]


@pytest.mark.parametrize("start,end", [("20030101.0100", "20030101.1200"),
                                      ("20030101.1200", "20030101.0000"),
                                      ("19941231.1800", "19950101.0600"),
                                      ("20061231.1800", "20070101.0600")])
def test_invalid_epoch_fails_before_download(tmp_path, grace_file, monkeypatch, start, end):
    grace_file.write_text(f"GRCOF2 0 0 1 0 0 0 {start} {end} nnnn\n")
    def unexpected_network(*args, **kwargs):
        pytest.fail("Invalid input must not open a network connection")
    monkeypatch.setattr(downloader.urllib.request, "urlopen", unexpected_network)
    with pytest.raises(ValueError):
        downloader.download_epoch(grace_file, ["H"], tmp_path)
    assert not (tmp_path / "ESM").exists()


def test_production_cli_defaults_to_full_record_without_grace(tmp_path, monkeypatch):
    monkeypatch.setenv("IFILES", str(tmp_path))
    calls = []
    monkeypatch.setattr(downloader, "download_dataset", lambda *args: calls.append(args))
    downloader.main([])
    assert calls == [(list(downloader.COMPONENTS), tmp_path, list(range(1995, 2007)))]
    downloader.main(["--years", "1996", "2000", "--components", "O"])
    assert calls[-1] == (["O"], tmp_path, [1996, 2000])
    assert list(tmp_path.iterdir()) == []


def test_full_record_schedule_covers_every_component_year(tmp_path, monkeypatch):
    calls = []
    def record_archive(output, component, year, expected_times):
        calls.append((component, year, sorted(expected_times)))
        return [{"timestamp": t.isoformat()} for t in sorted(expected_times)]
    monkeypatch.setattr(downloader, "download_archive", record_archive)
    paths = downloader.download_dataset(downloader.COMPONENTS, tmp_path)
    assert len(paths) == 12
    assert len(calls) == 72
    for component in downloader.COMPONENTS:
        times = [t for c, _, stamps in calls if c == component for t in stamps]
        assert len(times) == 17532
        assert times[0] == downloader.dt.datetime(1995, 1, 1)
        assert times[-1] == downloader.dt.datetime(2006, 12, 31, 18)
        assert all(b - a == downloader.dt.timedelta(hours=6) for a, b in zip(times, times[1:]))
    leap = json.loads(paths[1].read_text())
    assert "grace_file" not in leap
    assert len(leap["components"]["O"]) == 1464
    assert leap["end_exclusive"] == "1997-01-01T00:00:00"


def make_year_archive(component, year):
    """Small GFC payloads on the full annual time grid, in reverse order."""
    start = downloader.dt.datetime(year, 1, 1)
    end = downloader.dt.datetime(year + 1, 1, 1)
    step = downloader.dt.timedelta(hours=6)
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as archive:
        for k in reversed(range((end - start) // step)):
            date = start + k * step
            name = f"mtmshc_{component}_{date:%Y%m%d_%H}.gfc"
            content = f"modelname {name[:-4]}\nend_of_head\ngfc 0 0 {date.hour} 0\n".encode()
            member = tarfile.TarInfo(f"./{date:%m}/{name}")
            member.size = len(content)
            archive.addfile(member, io.BytesIO(content))
    return output.getvalue()


def test_year_resume_verifies_checkpoints_and_recovers_missing_files(tmp_path, monkeypatch):
    archive = make_year_archive("H", 1996)
    calls = []
    def open_archive(url, **kwargs):
        calls.append(url)
        return io.BytesIO(archive)
    monkeypatch.setattr(downloader.urllib.request, "urlopen", open_archive)
    path, = downloader.download_dataset(["H"], tmp_path, [1996])
    output = path.parent.parent
    entries = json.loads(path.read_text())["components"]["H"]
    assert len(entries) == 1464
    sample = output / entries[0]["file"]
    modified = sample.stat().st_mtime_ns
    downloader.download_dataset(["H"], tmp_path, [1996])
    assert len(calls) == 1  # All samples checked locally; no archive opened.
    (output / entries[-1]["file"]).unlink()
    downloader.download_dataset(["H"], tmp_path, [1996])
    assert len(calls) == 2
    assert sample.stat().st_mtime_ns == modified
    assert (output / entries[-1]["file"]).is_file()
    sample.write_text("modified locally")
    with pytest.raises(ValueError, match="differ from checkpoint"):
        downloader.download_dataset(["H"], tmp_path, [1996])
    assert len(calls) == 2
    assert sample.read_text() == "modified locally"


def test_year_checkpoints_survive_failure_and_preserve_demo_files(tmp_path, grace_file, monkeypatch):
    # The demo and production run share the same January sample paths.
    monkeypatch.setattr(downloader.urllib.request, "urlopen",
                        lambda *args, **kwargs: io.BytesIO(make_archive("H", [0, 6])))
    demo = downloader.download_epoch(grace_file, ["H"], tmp_path)
    sample = demo.parent.parent / "raw/H/mtmshc_H_20030101_00.gfc"
    modified = sample.stat().st_mtime_ns
    archive = make_year_archive("H", 2003)
    def interrupted(url, **kwargs):
        if "2003_H" in url:
            return io.BytesIO(archive)
        raise OSError("connection interrupted")
    monkeypatch.setattr(downloader.urllib.request, "urlopen", interrupted)
    with pytest.raises(OSError, match="connection interrupted"):
        downloader.download_dataset(["H"], tmp_path, [2003, 2004])
    annual = demo.parent / "ESAESM2_2003.json"
    assert len(json.loads(annual.read_text())["components"]["H"]) == 1460
    assert sample.stat().st_mtime_ns == modified
    assert len(json.loads(demo.read_text())["components"]["H"]) == 2
    calls = []
    next_archive = make_year_archive("H", 2004)
    def resume(url, **kwargs):
        calls.append(url)
        return io.BytesIO(next_archive)
    monkeypatch.setattr(downloader.urllib.request, "urlopen", resume)
    paths = downloader.download_dataset(["H"], tmp_path, [2003, 2004])
    assert len(paths) == 2
    assert calls == [f"{downloader.BASE}/mtmshc/mtm_2004_H.tar.gz"]


@pytest.mark.parametrize("years", [[], [1994], [2007], [2000.5], [True]])
def test_invalid_production_years_fail_before_io(tmp_path, monkeypatch, years):
    def unexpected_network(*args, **kwargs):
        pytest.fail("Invalid input must not open a network connection")
    monkeypatch.setattr(downloader.urllib.request, "urlopen", unexpected_network)
    with pytest.raises(ValueError, match="integer years"):
        downloader.download_dataset(["H"], tmp_path, years)
    assert list(tmp_path.iterdir()) == []


def test_grace_demo_and_year_selection_are_mutually_exclusive(tmp_path, monkeypatch):
    monkeypatch.setenv("IFILES", str(tmp_path))
    with pytest.raises(SystemExit) as error:
        downloader.parse_arguments(["--years", "2003", "--grace-file", "GSM"])
    assert error.value.code == 2
