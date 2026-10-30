"""Download the full ESAESM2 six-hour, degree-180 record (1995-2006).

Stream annual component archives into $IFILES/ESM/ESAESM2/raw/COMPONENT.
Completed component-years are checkpointed and verified locally on restart.
Use --years/--components for subsets, or --grace-file for the single-epoch
proof-of-concept demo. Archives need not have chronologically sorted members.
Saved files have SHA-256 integrity records; archive streams are not retained
or verified as complete gzip files.
"""

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import tarfile
import urllib.request

BASE = "https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2"
COMPONENTS = ("A", "H", "I", "O", "S", "AOHIS")
YEARS = tuple(range(1995, 2007))


def parse_arguments(argv=None):
    """Parse CLI options and resolve the configured download root."""
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--grace-file", type=Path,
                           help="Optional single-epoch demo: read the GRCOF2 interval")
    selection.add_argument("--years", type=int, nargs="+", choices=YEARS,
                           help="Whole years to download; default: all 1995-2006")
    parser.add_argument("--components", nargs="+", choices=COMPONENTS,
                        help="Default: all components (H I O for the GRACE demo)")
    args = parser.parse_args(argv)
    if args.components is None:
        args.components = ["H", "I", "O"] if args.grace_file else list(COMPONENTS)
    if args.years is None:
        args.years = list(YEARS)
    data_root = os.environ.get("IFILES")
    if not data_root:
        parser.error("Set IFILES to the configured data root")
    args.data_root = Path(data_root)
    return args


def read_grace_epoch(grace_file):
    """Read and validate the source interval; return its manifest and samples."""
    grace_bytes = grace_file.read_bytes()
    intervals = {tuple(line.split()[7:9]) for line in grace_bytes.decode("utf-8").splitlines()
                 if line.startswith("GRCOF2 ")}
    if len(intervals) != 1 or len(next(iter(intervals))) != 2:
        raise ValueError("Expected one consistent GRCOF2 start/end interval")
    start, end = [dt.datetime.strptime(s, "%Y%m%d.%H%M") for s in intervals.pop()]
    if start >= end or any(t.minute or t.second or t.hour % 6 for t in [start, end]):
        raise ValueError("Interval must have increasing, six-hour-aligned boundaries")
    step = dt.timedelta(hours=6)
    if start.year < 1995 or (end - step).year > 2006:
        raise ValueError("ESAESM2 six-hour fields cover 1995 through 2006")
    times = [start + k * step for k in range((end - start) // step)]
    manifest = dict(dataset="ESAESM2", start=start.isoformat(), end_exclusive=end.isoformat(),
                    grace_file=grace_file.name,
                    grace_sha256=hashlib.sha256(grace_bytes).hexdigest(),
                    sampling_hours=6, components={})
    return manifest, times


def load_manifest(path, manifest):
    """Preserve completed components only when source identity matches."""
    if path.exists():
        previous = json.loads(path.read_text())
        if any(previous.get(key) != value for key, value in manifest.items() if key != "components"):
            raise ValueError("Existing manifest belongs to a different GRACE input or dataset interval")
        manifest["components"] = previous["components"]
    return manifest


def save_manifest(path, manifest):
    """Replace the manifest after a whole component has completed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".json.part")
    temporary.write_text(json.dumps(manifest, indent=2) + "\n")
    temporary.replace(path)


def download_archive(output, component, year, expected_times):
    """Stream one yearly archive and retain exactly the requested samples."""
    dest = output / "raw" / component
    dest.mkdir(parents=True, exist_ok=True)
    found = {}
    url = f"{BASE}/mtmshc/mtm_{year}_{component}.tar.gz"
    print(f"Reading {url}; need {len(expected_times)} samples", flush=True)
    with urllib.request.urlopen(url, timeout=120) as response:
        with tarfile.open(fileobj=response, mode="r|gz") as archive:
            for member in archive:
                name = Path(member.name).name
                match = re.fullmatch(r"mtmshc_" + re.escape(component)
                                     + r"_(\d{8}_\d{2})\.gfc", name)
                if not member.isfile() or not match:
                    continue
                timestamp = dt.datetime.strptime(match[1], "%Y%m%d_%H")
                if timestamp not in expected_times:
                    continue
                if timestamp in found:
                    raise ValueError(f"Duplicate timestamp: {member.name}")
                data = archive.extractfile(member).read()
                if len(data) != member.size or b"end_of_head" not in data:
                    raise ValueError(f"Incomplete ESM file: {member.name}")
                path = dest / name
                if path.exists():
                    if path.read_bytes() != data:
                        raise ValueError(f"Existing data differ: {path}")
                else:
                    temporary = path.with_suffix(".gfc.part")
                    temporary.write_bytes(data)
                    temporary.replace(path)
                found[timestamp] = dict(file=path.relative_to(output).as_posix(),
                    timestamp=timestamp.isoformat(), bytes=len(data),
                    sha256=hashlib.sha256(data).hexdigest(), archive=url,
                    member=member.name)
                if len(found) % 30 == 0:
                    print(f"  {component}: {len(found)}/{len(expected_times)}", flush=True)
                if found.keys() == expected_times:
                    break
    if found.keys() != expected_times:
        raise ValueError(f"Missing samples for {component}: {expected_times - found.keys()}")
    return [found[t] for t in sorted(found)]


def download_epoch(grace_file, components, data_root):
    """Download an epoch using explicit inputs, independent of CLI/environment.

    Return the manifest path. A completed component is recorded even if a
    later component fails; existing samples must match the archived bytes.
    """
    if not components or any(component not in COMPONENTS for component in components):
        raise ValueError(f"Choose one or more components from {COMPONENTS}")
    output = Path(data_root) / "ESM" / "ESAESM2"
    manifest, times = read_grace_epoch(Path(grace_file))
    manifest_file = output / "manifests" / (manifest["grace_file"] + ".json")
    manifest = load_manifest(manifest_file, manifest)
    for component in components:
        entries = []
        for year in sorted({t.year for t in times}):
            expected_times = {t for t in times if t.year == year}
            entries.extend(download_archive(output, component, year, expected_times))
        manifest["components"][component] = entries
        save_manifest(manifest_file, manifest)
        print(f"Completed {component}: {len(entries)} files", flush=True)
    return manifest_file


def verified_component(output, component, entries, times):
    """Check a completed checkpoint before skipping its archive on restart.

    Missing files trigger another archive read. Modified files or inconsistent
    manifests raise an error rather than silently replacing local data.
    """
    if [entry["timestamp"] for entry in entries] != [t.isoformat() for t in times]:
        raise ValueError(f"Incomplete or inconsistent checkpoint for {component}")
    complete = True
    for entry, timestamp in zip(entries, times):
        relative = f"raw/{component}/mtmshc_{component}_{timestamp:%Y%m%d_%H}.gfc"
        if Path(entry["file"]) != Path(relative):
            raise ValueError(f"Unexpected checkpoint path: {entry['file']}")
        path = output / relative
        if not path.is_file():
            complete = False
            continue
        data = path.read_bytes()
        if len(data) != entry["bytes"] or hashlib.sha256(data).hexdigest() != entry["sha256"]:
            raise ValueError(f"Existing data differ from checkpoint: {path}")
    return complete


def download_dataset(components, data_root, years=YEARS):
    """Download complete years without requiring a GRACE observation.

    Defaults to the full temporal record. Save one manifest per year, after
    each component completes, so failures do not discard earlier progress.
    Return the annual manifest paths. Reuse existing epoch-demo raw files.
    """
    if not components or any(c not in COMPONENTS for c in components):
        raise ValueError(f"Choose one or more components from {COMPONENTS}")
    if not years or any(isinstance(y, bool) or not isinstance(y, int) or y not in YEARS for y in years):
        raise ValueError("ESAESM2 six-hour fields cover integer years 1995 through 2006")
    output = Path(data_root) / "ESM" / "ESAESM2"
    manifest_files = []
    for year in sorted(set(years)):
        start, end = dt.datetime(year, 1, 1), dt.datetime(year + 1, 1, 1)
        step = dt.timedelta(hours=6)
        times = [start + k * step for k in range((end - start) // step)]
        manifest = dict(dataset="ESAESM2", start=start.isoformat(),
                        end_exclusive=end.isoformat(), sampling_hours=6, components={})
        path = output / "manifests" / f"ESAESM2_{year}.json"
        manifest = load_manifest(path, manifest)
        for component in dict.fromkeys(components):
            previous = manifest["components"].get(component)
            if previous is not None and verified_component(output, component, previous, times):
                print(f"Verified {component} {year}: {len(previous)} existing files", flush=True)
                continue
            entries = download_archive(output, component, year, set(times))
            manifest["components"][component] = entries
            save_manifest(path, manifest)
            print(f"Completed {component} {year}: {len(entries)} files", flush=True)
        manifest_files.append(path)
    return manifest_files


def main(argv=None):
    args = parse_arguments(argv)
    if args.grace_file:
        download_epoch(args.grace_file, args.components, args.data_root)
    else:
        download_dataset(args.components, args.data_root, args.years)


if __name__ == "__main__":
    main()
