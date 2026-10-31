"""Download the full ESAESM2 six-hour, degree-180 record (1995-2006).

Stream annual component archives into $IFILES/ESM/ESAESM2/raw/COMPONENT.
Completed component-years are checkpointed and verified locally on restart.
Use --components to select components; all years are always downloaded.
Archives need not have chronologically sorted members.
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
    parser.add_argument(
        "--components",
        nargs="+",
        choices=COMPONENTS,
        default=list(COMPONENTS),
        help="Components to download; default: all components",
    )
    args = parser.parse_args(argv)
    data_root = os.environ.get("IFILES")
    if not data_root:
        parser.error("Set IFILES to the configured data root")
    args.data_root = Path(data_root)
    return args


def load_manifest(path, manifest):
    """Preserve completed components only when source identity matches."""
    if path.exists():
        previous = json.loads(path.read_text())
        if any(
            previous.get(key) != value
            for key, value in manifest.items()
            if key != "components"
        ):
            raise ValueError(
                "Existing manifest belongs to a different dataset interval"
            )
        # Retain completed components while allowing new components to be added.
        manifest["components"] = previous["components"]
    return manifest


def save_manifest(path, manifest):
    """Replace the manifest after a whole component has completed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    # Publish a complete JSON document only after writing its temporary file.
    temporary = path.with_suffix(".json.part")
    temporary.write_text(json.dumps(manifest, indent=2) + "\n")
    temporary.replace(path)


def download_archive(output, comp, year, expected_times):
    """Stream one yearly archive and retain exactly the requested samples."""
    dest = output / "raw" / comp
    dest.mkdir(parents=True, exist_ok=True)
    found = {}
    url = f"{BASE}/mtmshc/mtm_{year}_{comp}.tar.gz"
    print(f"Reading {url}; need {len(expected_times)} samples", flush=True)
    with urllib.request.urlopen(url, timeout=120) as response:
        # Stream the compressed archive instead of saving a full annual download.
        with tarfile.open(fileobj=response, mode="r|gz") as archive:
            for member in archive:
                # Use only a validated basename; archive directory paths are not extracted.
                name = Path(member.name).name
                match = re.fullmatch(
                    r"mtmshc_" + re.escape(comp) + r"_(\d{8}_\d{2})\.gfc", name
                )
                if not member.isfile() or not match:
                    continue
                timestamp = dt.datetime.strptime(match[1], "%Y%m%d_%H")
                # Select by timestamp rather than archive order, which may be unsorted.
                if timestamp not in expected_times:
                    continue
                if timestamp in found:
                    raise ValueError(f"Duplicate timestamp: {member.name}")
                data = archive.extractfile(member).read()
                if len(data) != member.size or b"end_of_head" not in data:
                    raise ValueError(f"Incomplete ESM file: {member.name}")
                path = dest / name
                # Reuse identical samples from earlier runs; preserve and reject changed data.
                if path.exists():
                    if path.read_bytes() != data:
                        raise ValueError(f"Existing data differ: {path}")
                else:
                    tmp_path = path.with_suffix(".gfc.part")
                    tmp_path.write_bytes(data)
                    tmp_path.replace(path)
                # Record sample integrity and archive provenance using portable relative paths.
                found[timestamp] = dict(
                    file=path.relative_to(output).as_posix(),
                    timestamp=timestamp.isoformat(),
                    bytes=len(data),
                    sha256=hashlib.sha256(data).hexdigest(),
                    archive=url,
                    member=member.name,
                )
                if len(found) % 30 == 0:
                    print(f"  {comp}: {len(found)}/{len(expected_times)}", flush=True)
                # Stop once all requested timestamps are found; the archive tail is not checked.
                if found.keys() == expected_times:
                    break
    # Reject partial coverage before returning chronologically sorted entries.
    if found.keys() != expected_times:
        raise ValueError(f"Missing samples for {comp}: {expected_times - found.keys()}")
    return [found[t] for t in sorted(found)]


def verified_component(output, comp, entries, times):
    """Check a completed checkpoint before skipping its archive on restart.

    Missing files trigger another archive read. Modified files or inconsistent
    manifests raise an error rather than silently replacing local data.
    """
    # Require the exact ordered time grid before trusting a completed checkpoint.
    if [entry["timestamp"] for entry in entries] != [t.isoformat() for t in times]:
        raise ValueError(f"Incomplete or inconsistent checkpoint for {comp}")
    complete = True
    for entry, timestamp in zip(entries, times):
        # Derive the expected path independently of the checkpoint entry.
        relative = f"raw/{comp}/mtmshc_{comp}_{timestamp:%Y%m%d_%H}.gfc"
        if Path(entry["file"]) != Path(relative):
            raise ValueError(f"Unexpected checkpoint path: {entry['file']}")
        path = output / relative
        # Missing samples allow a reread, but continue checking other files for changes.
        if not path.is_file():
            complete = False
            continue
        data = path.read_bytes()
        if (
            len(data) != entry["bytes"]
            or hashlib.sha256(data).hexdigest() != entry["sha256"]
        ):
            raise ValueError(f"Existing data differ from checkpoint: {path}")
    return complete


def download_dataset(comps, data_root):
    """Download the full temporal record for the selected components.

    Save one manifest per year, after
    each component completes, so failures do not discard earlier progress.
    Return the annual manifest paths. Reuse existing epoch-demo raw files.
    """
    if not comps or any(c not in COMPONENTS for c in comps):
        raise ValueError(f"Choose one or more components from {COMPONENTS}")
    output = Path(data_root) / "ESM" / "ESAESM2"
    manifest_files = []
    for year in YEARS:
        start, end = dt.datetime(year, 1, 1), dt.datetime(year + 1, 1, 1)
        step = dt.timedelta(hours=6)
        # Derive the annual grid from calendar boundaries, including leap days.
        times = [start + k * step for k in range((end - start) // step)]
        manifest = dict(
            dataset="ESAESM2",
            start=start.isoformat(),
            end_exclusive=end.isoformat(),
            sampling_hours=6,
            components={},
        )
        path = output / "manifests" / f"ESAESM2_{year}.json"
        manifest = load_manifest(path, manifest)
        # Deduplicate component requests while preserving their requested order.
        for comp in dict.fromkeys(comps):
            prev = manifest["components"].get(comp)
            # Skip the archive only after every checkpointed sample passes local checks.
            if prev is not None and verified_component(output, comp, prev, times):
                print(
                    f"Verified {comp} {year}: {len(prev)} existing files",
                    flush=True,
                )
                continue
            entries = download_archive(output, comp, year, set(times))
            manifest["components"][comp] = entries
            # Save immediately so completed component-years survive interrupted runs.
            save_manifest(path, manifest)
            print(f"Completed {comp} {year}: {len(entries)} files", flush=True)
        manifest_files.append(path)
    return manifest_files


# Dispatch parsed options to processing functions with explicit inputs.
def main(argv=None):
    args = parse_arguments(argv)
    download_dataset(args.components, args.data_root)


if __name__ == "__main__":
    main()
