"""Validate the complete release inventory before touching a draft release."""

import hashlib
import os
from pathlib import Path
import re

tag = os.environ["RELEASE_TAG"]
if not re.fullmatch(r"v\d+\.\d+\.\d+", tag):
    raise SystemExit("Expected a stable version tag such as v0.1.1")
version = tag[1:]
prefix = f"netOP_{version}"
expected = {f"{prefix}.tar.gz"}
for series in ("4.4", "4.5", "4.6"):
    for architecture, extension in (
        ("arm64", "tgz"), ("x86_64", "tgz"), ("x86_64", "zip")
    ):
        expected.add(f"{prefix}_R-{series}_{architecture}.{extension}")

directory = Path("release-assets")
actual = {path.name for path in directory.iterdir()}
if actual != expected:
    raise SystemExit(
        f"Asset mismatch: missing={sorted(expected - actual)}, "
        f"unexpected={sorted(actual - expected)}"
    )

checksums = []
for name in sorted(expected):
    path = directory / name
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"Missing or empty asset: {name}")
    checksums.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {name}\n")
(directory / "SHA256SUMS").write_text("".join(checksums))
print("Validated nine binaries and one source tarball; wrote SHA256SUMS.")
