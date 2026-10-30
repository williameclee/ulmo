"""Offline tests for ESM verification and its separate command-line interface."""

import hashlib
import importlib.util
import json
from pathlib import Path
import struct

import pytest

MODULE = Path(__file__).resolve().parents[1] / "data" / "verify_esm_epoch.py"
SPEC = importlib.util.spec_from_file_location("verify_esm_epoch", MODULE)
verifier = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(verifier)


@pytest.fixture
def manifest_file(tmp_path):
    folder = tmp_path / "ESM/ESAESM2"
    (folder / "raw/H").mkdir(parents=True)
    (folder / "manifests").mkdir()
    entries = []
    for hour, amplitude in [(0, 2), (6, 4)]:
        lines = [f"gfc {degree} {order} {amplitude * (degree + 1)} {-amplitude * order}\n"
                 for degree in range(181) for order in range(degree + 1)]
        data = ("end_of_head\n" + "".join(lines)).encode("ascii")
        relative = f"raw/H/mtmshc_H_20030101_{hour:02d}.gfc"
        (folder / relative).write_bytes(data)
        entries.append(dict(file=relative, timestamp=f"2003-01-01T{hour:02d}:00:00",
                            bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
    manifest = dict(start="2003-01-01T00:00:00", end_exclusive="2003-01-01T12:00:00",
                    components={"H": entries})
    path = folder / "manifests/GSM-fixture.json"
    path.write_text(json.dumps(manifest))
    return path


def test_parser_does_not_access_paths(tmp_path):
    manifest = tmp_path / "missing.json"
    output = tmp_path / "output"
    args = verifier.parse_arguments([str(manifest), "--output", str(output)])
    assert args.manifest == manifest
    assert args.output == output
    assert not output.exists()


def test_cli_saves_independent_means_and_report(tmp_path, manifest_file, capsys):
    output = tmp_path / "verified"
    verifier.main([str(manifest_file), "--output", str(output)])
    report = json.loads((output / "raw_validation.json").read_text())
    assert json.loads(capsys.readouterr().out) == report
    assert report["H"]["files"] == 2
    assert report["H"]["coefficient_rows"] == 16471
    assert report["H"]["checksums_verified"]
    expected = [value for degree in range(181) for order in range(degree + 1)
                for value in (3 * (degree + 1), -3 * order)]
    actual = struct.unpack(f"<{len(expected)}d", (output / "reference_H.bin").read_bytes())
    assert list(actual) == expected


@pytest.mark.parametrize("corruption", ["checksum", "duplicate_harmonic", "missing_harmonic", "timestamp"])
def test_verification_rejects_corrupted_input(tmp_path, manifest_file, corruption):
    manifest = json.loads(manifest_file.read_text())
    entry = manifest["components"]["H"][0]
    path = manifest_file.parent.parent / entry["file"]
    if corruption == "timestamp":
        entry["timestamp"] = "2003-01-01T06:00:00"
    else:
        data = path.read_text()
        if corruption == "checksum":
            data = data.replace("gfc 0 0 2 0", "gfc 0 0 9 0")
        elif corruption == "duplicate_harmonic":
            data += "gfc 0 0 2 0\n"
        else:
            data = data.replace("gfc 0 0 2 0\n", "")
        path.write_text(data)
        if corruption != "checksum":
            entry["bytes"] = path.stat().st_size
            entry["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_file.write_text(json.dumps(manifest))
    with pytest.raises(AssertionError):
        verifier.verify_epoch(manifest_file, tmp_path / "verified")
    assert not (tmp_path / "verified/reference_H.bin").exists()
