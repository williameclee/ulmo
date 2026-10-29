"""Verify retained ESM files and independently average their Stokes records.

Writes little-endian doubles (C,S interleaved in degree/order sequence) for
comparison with MATLAB. Uses only the Python standard library.
"""

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import struct


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("--output", type=Path, required=True,
                        help="Destination for reference coefficients and verification report")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.data / "manifest.json").read_text())
    start = dt.datetime.fromisoformat(manifest["start"])
    end = dt.datetime.fromisoformat(manifest["end_exclusive"])
    expected_times = []
    t = start
    while t < end:
        expected_times.append(t.isoformat())
        t += dt.timedelta(hours=6)
    report = {}
    for component, entries in manifest["components"].items():
        assert [entry["timestamp"] for entry in entries] == expected_times
        nrows = 181 * 182 // 2
        mean = [0.0] * (2 * nrows)
        total_bytes = 0
        for entry in entries:
            data = (args.data / entry["file"]).read_bytes()
            assert len(data) == entry["bytes"]
            assert hashlib.sha256(data).hexdigest() == entry["sha256"]
            total_bytes += len(data)
            seen = set()
            for line in data.decode("ascii").splitlines():
                if not line.startswith("gfc "):
                    continue
                tag, degree, order, cosine, sine = line.split()
                degree, order = int(degree), int(order)
                assert 0 <= order <= degree <= 180
                index = degree * (degree + 1) // 2 + order
                assert index not in seen
                seen.add(index)
                # Divide each sample before summing, unlike the MATLAB processor.
                mean[2 * index] += float(cosine) / len(entries)
                mean[2 * index + 1] += float(sine) / len(entries)
            assert len(seen) == nrows
        reference = args.output / f"reference_{component}.bin"
        reference.write_bytes(struct.pack(f"<{len(mean)}d", *mean))
        report[component] = dict(files=len(entries), bytes=total_bytes,
                                 coefficient_rows=nrows, checksums_verified=True)
    (args.output / "raw_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
