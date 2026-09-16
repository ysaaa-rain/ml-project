"""Download the official RegulonDB promoter and E. coli K-12 genome assets."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from datetime import date
from pathlib import Path
from typing import Any
from urllib.request import Request, urlopen


DEFAULT_FILES = {
    "promoters_gff3": "https://regulondb.ccg.unam.mx/media/raw/gff3/PromoterSet.gff3",
    "genome_fasta": "https://regulondb.ccg.unam.mx/media/raw/e_coli_k12.fna",
}

DEFAULT_OUTPUT_NAMES = {
    "promoters_gff3": "regulondb_promoters_20260916.gff3",
    "genome_fasta": "regulondb_e_coli_k12_20260916.fna",
}


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _fetch(url: str, *, retries: int = 3, timeout: int = 120) -> tuple[bytes, dict[str, str]]:
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            request = Request(url, headers={"User-Agent": "ml-project/PR01-02 source audit"})
            with urlopen(request, timeout=timeout) as response:
                headers = {
                    "etag": response.headers.get("ETag", ""),
                    "content_type": response.headers.get("Content-Type", ""),
                }
                return response.read(), headers
        except Exception as error:  # pragma: no cover - network-dependent branch
            last_error = error
            if attempt + 1 < retries:
                time.sleep(float(attempt + 1))
    raise RuntimeError(f"failed to download {url}: {last_error}")


def download_assets(
    *,
    raw_dir: str | Path,
    report_output: str | Path,
    files: dict[str, str] | None = None,
    output_names: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Snapshot the official RegulonDB raw assets and write a hash manifest."""

    raw_path = Path(raw_dir)
    raw_path.mkdir(parents=True, exist_ok=True)
    selected = dict(files or DEFAULT_FILES)
    names = dict(DEFAULT_OUTPUT_NAMES)
    names.update(output_names or {})
    assets: list[dict[str, Any]] = []
    for asset_id, url in selected.items():
        filename = names.get(asset_id, Path(url.rstrip("/")).name)
        if not filename or filename in {".", ".."} or "/" in filename or "\\" in filename:
            raise ValueError(f"unsafe RegulonDB asset filename: {filename!r}")
        payload, headers = _fetch(url)
        output = raw_path / filename
        output.write_bytes(payload)
        assets.append(
            {
                "asset_id": asset_id,
                "url": url,
                "local_path": str(output.resolve()),
                "bytes": len(payload),
                "sha256": sha256_bytes(payload),
                **headers,
            }
        )

    report = {
        "source_id": "regulondb",
        "source_version": "official portal raw-media snapshot",
        "access_date": date.today().isoformat(),
        "landing_page": "https://regulondb.ccg.unam.mx/",
        "documentation": "https://github.com/regulondbunam/regulondb-docs/blob/main/00_about_policies/about_us.md",
        "assets": assets,
        "note": "Raw files are ignored by Git; the report records provenance and hashes for local reproduction.",
    }
    report_path = Path(report_output)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument(
        "--report-output",
        default="data/interim/regulondb_download_report_20260916.json",
    )
    args = parser.parse_args()
    report = download_assets(raw_dir=args.raw_dir, report_output=args.report_output)
    for asset in report["assets"]:
        print(f"{asset['asset_id']}_bytes={asset['bytes']}")
        print(f"{asset['asset_id']}_sha256={asset['sha256']}")
    print(f"report={Path(args.report_output).resolve()}")


if __name__ == "__main__":
    main()
