#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


MANIFEST_NAME = "ANALYSIS_MANIFEST.csv"
CHECKSUM_NAME = "SHA256SUMS.csv"
TRANSIENT_PATTERNS = (
    "*.part",
    "*.range*",
    "*.pyc",
    ".DS_Store",
    ".Rhistory",
    ".RData",
    "*~",
)
TRANSIENT_DIRECTORIES = {"__pycache__", ".pytest_cache", ".ipynb_checkpoints"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_transient(path: Path) -> bool:
    return any(part in TRANSIENT_DIRECTORIES for part in path.parts) or any(
        path.match(pattern) for pattern in TRANSIENT_PATTERNS
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the analysis package manifest and SHA-256 list."
    )
    parser.add_argument(
        "--package",
        type=Path,
        required=True,
        help="Submission package root.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Directory for the manifest files. Defaults to "
            "supplementary/S13_reproducibility/manifests under the package."
        ),
    )
    parser.add_argument(
        "--path-prefix",
        default="",
        help="Optional prefix added to every manifest relative_path.",
    )
    args = parser.parse_args()

    package = args.package.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else package / "supplementary" / "S13_reproducibility" / "manifests"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / MANIFEST_NAME
    checksum_path = output_dir / CHECKSUM_NAME
    excluded = {manifest_path.resolve(), checksum_path.resolve()}

    rows = []
    zero_byte_files = []
    for path in sorted(package.rglob("*")):
        if (
            not path.is_file()
            or path.resolve() in excluded
            or is_transient(path)
        ):
            continue
        relative_path = path.relative_to(package)
        if any("Internal_Audit" in part for part in relative_path.parts):
            continue
        relative = relative_path.as_posix()
        if args.path_prefix:
            relative = f"{args.path_prefix.strip('/')}/{relative}"
        size = path.stat().st_size
        if size == 0:
            zero_byte_files.append(relative)
        rows.append((relative, size, sha256(path)))

    if zero_byte_files:
        raise SystemExit(
            "Refusing to build a manifest with zero-byte files: "
            + ", ".join(zero_byte_files)
        )

    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["relative_path", "size_bytes", "sha256"])
        writer.writerows(rows)

    with checksum_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["relative_path", "sha256"])
        writer.writerows((relative, digest) for relative, _, digest in rows)

    print(
        f"manifest_files={len(rows)} "
        f"manifest={manifest_path.relative_to(package)} "
        f"checksums={checksum_path.relative_to(package)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
