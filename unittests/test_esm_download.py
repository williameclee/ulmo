"""Offline ESM archive tests: python3 -B -m pytest unittests/test_esm_download.py."""

import importlib.util
import io
import json
from pathlib import Path
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


def test_component_only_cli(tmp_path, monkeypatch):
    monkeypatch.setenv("IFILES", str(tmp_path))
    calls = []
    monkeypatch.setattr(downloader, "download_dataset", lambda *args: calls.append(args))
    downloader.main([])
    downloader.main(["--components", "H", "O"])
    assert calls == [(list(downloader.COMPONENTS), tmp_path), (["H", "O"], tmp_path)]
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("arguments", [["--years", "2003"], ["--grace-file", "GSM"],
                                       ["--components", "invalid"]])
def test_invalid_cli_arguments(tmp_path, monkeypatch, arguments):
    monkeypatch.setenv("IFILES", str(tmp_path))
    with pytest.raises(SystemExit) as error:
        downloader.parse_arguments(arguments)
    assert error.value.code == 2


def test_missing_data_root_fails(monkeypatch, capsys):
    monkeypatch.delenv("IFILES", raising=False)
    with pytest.raises(SystemExit) as error:
        downloader.parse_arguments([])
    assert error.value.code == 2
    assert "Set IFILES" in capsys.readouterr().err


@pytest.mark.parametrize("hours, error", [([0], "Missing samples"),
                                        ([0, 0, 6], "Duplicate timestamp")])
def test_missing_and_duplicate_samples_fail(tmp_path, monkeypatch, hours, error):
    monkeypatch.setattr(downloader.urllib.request, "urlopen",
                        lambda *args, **kwargs: io.BytesIO(make_archive("H", hours)))
    times = {downloader.dt.datetime(2003, 1, 1, h) for h in [0, 6]}
    with pytest.raises(ValueError, match=error):
        downloader.download_archive(tmp_path, "H", 2003, times)


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
    monkeypatch.setattr(downloader, "YEARS", (1996,))
    archive = make_year_archive("H", 1996)
    calls = []
    def open_archive(url, **kwargs):
        calls.append(url)
        return io.BytesIO(archive)
    monkeypatch.setattr(downloader.urllib.request, "urlopen", open_archive)
    path, = downloader.download_dataset(["H"], tmp_path)
    output = path.parent.parent
    entries = json.loads(path.read_text())["components"]["H"]
    assert len(entries) == 1464
    sample = output / entries[0]["file"]
    modified = sample.stat().st_mtime_ns
    downloader.download_dataset(["H"], tmp_path)
    assert len(calls) == 1  # All samples checked locally; no archive opened.
    (output / entries[-1]["file"]).unlink()
    downloader.download_dataset(["H"], tmp_path)
    assert len(calls) == 2
    assert sample.stat().st_mtime_ns == modified
    assert (output / entries[-1]["file"]).is_file()
    sample.write_text("modified locally")
    with pytest.raises(ValueError, match="differ from checkpoint"):
        downloader.download_dataset(["H"], tmp_path)
    assert len(calls) == 2
    assert sample.read_text() == "modified locally"


def test_year_checkpoints_survive_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(downloader, "YEARS", (2003, 2004))
    archive = make_year_archive("H", 2003)
    def interrupted(url, **kwargs):
        if "2003_H" in url:
            return io.BytesIO(archive)
        raise OSError("connection interrupted")
    monkeypatch.setattr(downloader.urllib.request, "urlopen", interrupted)
    with pytest.raises(OSError, match="connection interrupted"):
        downloader.download_dataset(["H"], tmp_path)
    annual = tmp_path / "ESM/ESAESM2/manifests/ESAESM2_2003.json"
    assert len(json.loads(annual.read_text())["components"]["H"]) == 1460
    calls = []
    next_archive = make_year_archive("H", 2004)
    def resume(url, **kwargs):
        calls.append(url)
        return io.BytesIO(next_archive)
    monkeypatch.setattr(downloader.urllib.request, "urlopen", resume)
    paths = downloader.download_dataset(["H"], tmp_path)
    assert len(paths) == 2
    assert calls == [f"{downloader.BASE}/mtmshc/mtm_2004_H.tar.gz"]


def test_add_component_preserves_manifest(tmp_path, monkeypatch):
    monkeypatch.setattr(downloader, "YEARS", (2003,))
    for component in ["H", "O"]:
        archive = make_year_archive(component, 2003)
        monkeypatch.setattr(downloader.urllib.request, "urlopen",
                            lambda *args, **kwargs: io.BytesIO(archive))
        path, = downloader.download_dataset([component], tmp_path)
    assert set(json.loads(path.read_text())["components"]) == {"H", "O"}
