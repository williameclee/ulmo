"""Stream ESAESM2 yearly archives, retaining only one GRACE interval.

Read the interval from GRCOF2 records, not the inclusive day-of-year filename.
Archive members are not chronologically sorted. Stop only after every expected
six-hour sample has been found. Individual saved files have SHA-256 checksums;
the partial archive stream itself is not a verified full-year gzip download.
Inputs are saved under $IFILES/ESM/ESAESM2/raw/COMPONENT. Per-GRACE-source
manifests share that root, so overlapping epochs reuse the same raw files.
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grace-file", type=Path, required=True)
    parser.add_argument("--components", nargs="+", default=["H", "I", "O"],
                        choices=["A", "H", "I", "O", "S", "AOHIS"])
    args = parser.parse_args()
    data_root = os.environ.get("IFILES")
    if not data_root:
        parser.error("Set IFILES to the configured data root")
    output = Path(data_root) / "ESM" / "ESAESM2"
    grace_bytes = args.grace_file.read_bytes()
    grace_text = grace_bytes.decode('utf-8')
    intervals = {tuple(line.split()[7:9]) for line in grace_text.splitlines()
                 if line.startswith("GRCOF2 ")}
    if len(intervals) != 1 or len(next(iter(intervals))) != 2:
        raise ValueError("Expected one consistent GRCOF2 start/end interval")
    start, end = [dt.datetime.strptime(s, "%Y%m%d.%H%M") for s in intervals.pop()]
    if start >= end or any(t.minute or t.second or t.hour % 6 for t in [start, end]):
        raise ValueError("Interval must have increasing, six-hour-aligned boundaries")
    times = []
    t = start
    while t < end:
        times.append(t)
        t += dt.timedelta(hours=6)
    if start.year < 1995 or times[-1].year > 2006:
        raise ValueError("ESAESM2 six-hour fields cover 1995 through 2006")
    manifest_folder = output / "manifests"
    manifest_folder.mkdir(parents=True, exist_ok=True)
    manifest = dict(dataset="ESAESM2", start=start.isoformat(), end_exclusive=end.isoformat(),
                    grace_file=args.grace_file.name,
                    grace_sha256=hashlib.sha256(grace_bytes).hexdigest(),
                    sampling_hours=6, components={})
    manifest_file = manifest_folder / (args.grace_file.name + ".json")
    if manifest_file.exists():
        previous = json.loads(manifest_file.read_text())
        if any(previous.get(key) != value for key, value in manifest.items() if key != "components"):
            raise ValueError("Existing manifest belongs to a different GRACE input or interval")
        manifest["components"] = previous["components"]
    for component in args.components:
        dest = output / "raw" / component
        dest.mkdir(parents=True, exist_ok=True)
        entries = []
        for year in sorted({t.year for t in times}):
            expected_times = {t for t in times if t.year == year}
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
                        found[timestamp] = dict(file=str(path.relative_to(output)),
                            timestamp=timestamp.isoformat(), bytes=len(data),
                            sha256=hashlib.sha256(data).hexdigest(), archive=url,
                            member=member.name)
                        if len(found) % 30 == 0:
                            print(f"  {component}: {len(found)}/{len(expected_times)}", flush=True)
                        if found.keys() == expected_times:
                            break
            if found.keys() != expected_times:
                raise ValueError(f"Missing samples for {component}: {expected_times - found.keys()}")
            entries.extend(found[t] for t in sorted(found))
        manifest["components"][component] = entries
        # Persist each completed component even if a later connection fails.
        temporary_manifest = manifest_file.with_suffix(".json.part")
        temporary_manifest.write_text(json.dumps(manifest, indent=2) + "\n")
        temporary_manifest.replace(manifest_file)
        print(f"Completed {component}: {len(entries)} files", flush=True)


if __name__ == "__main__":
    main()
