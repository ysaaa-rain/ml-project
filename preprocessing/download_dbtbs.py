"""Download and parse the official DBTBS release 4.1 promoter pages."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import urljoin, urlparse
from urllib.request import Request, urlopen

import pandas as pd


DEFAULT_BASE_URL = "https://dbtbs.hgc.jp/ver4/"
DNA_RE = re.compile(r"^[ACGTN]+$")
POSITION_RE = re.compile(r"^(\d+)\.\.(\d+)$")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


class _HTMLTableParser(HTMLParser):
    """Collect table rows while retaining line breaks inside evidence cells."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.rows: list[list[str]] = []
        self._row: list[str] | None = None
        self._cell: list[str] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        tag = tag.lower()
        if tag == "tr":
            self._row = []
        elif tag in {"td", "th"} and self._row is not None:
            self._cell = []
        elif tag == "br" and self._cell is not None:
            self._cell.append("\n")

    def handle_endtag(self, tag: str) -> None:
        tag = tag.lower()
        if tag in {"td", "th"} and self._row is not None and self._cell is not None:
            self._row.append("".join(self._cell))
            self._cell = None
        elif tag == "tr" and self._row is not None:
            if self._row:
                self.rows.append(self._row)
            self._row = None

    def handle_data(self, data: str) -> None:
        if self._cell is not None:
            self._cell.append(data)


def _clean_text(value: str) -> str:
    lines = [re.sub(r"\s+", " ", part).strip() for part in value.split("\n")]
    return "; ".join(part for part in lines if part)


def _clean_sequence(value: str) -> str:
    return re.sub(r"\s+", "", value).upper()


def parse_index_links(html: str, *, base_url: str = DEFAULT_BASE_URL) -> list[str]:
    """Return unique promoter page URLs linked by the official index."""

    links = re.findall(r"href\s*=\s*[\"'](COG/prom/[^\"']+)[\"']", html, re.IGNORECASE)
    return sorted({urljoin(base_url, link) for link in links})


def extract_promoter_rows(html: str, *, page_name: str) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Extract experimentally listed sigma-promoter rows from one DBTBS page."""

    parser = _HTMLTableParser()
    parser.feed(html)
    rows = parser.rows
    header_index = None
    for index, row in enumerate(rows):
        normalized = [re.sub(r"[^a-z]", "", cell.lower()) for cell in row]
        if any("bindingfactor" in cell for cell in normalized) and any(
            "bindingseq" in cell for cell in normalized
        ):
            header_index = index
            break
    counts = {
        "promoter_rows_seen": 0,
        "non_sigma_promoters_skipped": 0,
        "invalid_or_missing_sequence_skipped": 0,
        "position_parse_failures": 0,
    }
    if header_index is None:
        return [], counts

    records: list[dict[str, Any]] = []
    for row_index, row in enumerate(rows[header_index + 1 :], start=1):
        if len(row) < 6:
            continue
        binding_factor = _clean_text(row[0])
        regulation = _clean_text(row[1])
        if regulation.lower() != "promoter":
            continue
        counts["promoter_rows_seen"] += 1
        if not re.match(r"^sig(?:ma)?[a-z0-9]+$", binding_factor, re.IGNORECASE):
            counts["non_sigma_promoters_skipped"] += 1
            continue
        location = _clean_text(row[2])
        absolute_position = _clean_text(row[3])
        sequence = _clean_sequence(row[4])
        evidence = _clean_text(row[5])
        if not sequence or not DNA_RE.fullmatch(sequence):
            counts["invalid_or_missing_sequence_skipped"] += 1
            continue
        position_match = POSITION_RE.fullmatch(absolute_position)
        if position_match is None:
            counts["position_parse_failures"] += 1
            source_start = pd.NA
            source_end = pd.NA
        else:
            source_start = int(position_match.group(1))
            source_end = int(position_match.group(2))
        source_record_id = f"dbtbs_v4.1:{page_name}:promoter:{row_index}"
        records.append(
            {
                "sequence_id": source_record_id,
                "species": "Bacillus subtilis",
                "taxon_id": pd.NA,
                "assembly_accession": pd.NA,
                "source_dataset": "dbtbs_v4.1",
                "sequence": sequence,
                "sequence_length": len(sequence),
                "tss_position": pd.NA,
                "strand": pd.NA,
                "sigma_factor_type": binding_factor,
                "known_element_annotations": pd.NA,
                "promoter_strength": pd.NA,
                "evidence_level": "experimental",
                "source_record_id": source_record_id,
                "source_page": page_name,
                "promoter_location": location,
                "source_tss_coordinate": source_start,
                "source_tss_end_coordinate": source_end,
                "source_evidence": evidence,
                "source_binding_factor": binding_factor,
            }
        )
    return records, counts


def _fetch(url: str, *, retries: int = 3, timeout: int = 60) -> bytes:
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            request = Request(url, headers={"User-Agent": "ml-project/PR01-02 source audit"})
            with urlopen(request, timeout=timeout) as response:
                return response.read()
        except Exception as error:  # pragma: no cover - network-dependent branch
            last_error = error
            if attempt + 1 < retries:
                time.sleep(1.0 * (attempt + 1))
    raise RuntimeError(f"failed to download {url}: {last_error}")


def _page_filename(url: str) -> str:
    name = Path(urlparse(url).path).name
    if not name or name in {".", ".."} or "/" in name or "\\" in name:
        raise ValueError(f"unsafe DBTBS page filename: {name!r}")
    return name


def _download_page(url: str, pages_dir: Path) -> tuple[str, bytes | None, str | None]:
    filename = _page_filename(url)
    try:
        payload = _fetch(url)
        (pages_dir / filename).write_bytes(payload)
        return url, payload, None
    except Exception as error:  # pragma: no cover - network-dependent branch
        return url, None, str(error)


def download_and_parse(
    *,
    raw_dir: str | Path,
    metadata_output: str | Path,
    report_output: str | Path,
    base_url: str = DEFAULT_BASE_URL,
    workers: int = 4,
    max_pages: int | None = None,
) -> dict[str, Any]:
    """Snapshot DBTBS index/pages and write a canonical promoter metadata TSV."""

    if workers < 1:
        raise ValueError("workers must be at least 1")
    raw_path = Path(raw_dir)
    pages_dir = raw_path / "pages"
    pages_dir.mkdir(parents=True, exist_ok=True)
    index_url = urljoin(base_url, "promtable.html")
    index_bytes = _fetch(index_url)
    index_path = raw_path / "promtable.html"
    index_path.write_bytes(index_bytes)
    page_urls = parse_index_links(index_bytes.decode("utf-8", errors="replace"), base_url=base_url)
    if max_pages is not None:
        if max_pages < 1:
            raise ValueError("max_pages must be at least 1 when supplied")
        page_urls = page_urls[:max_pages]

    all_records: list[dict[str, Any]] = []
    page_hashes: dict[str, str] = {}
    failed_pages: dict[str, str] = {}
    aggregate_counts = {
        "promoter_rows_seen": 0,
        "non_sigma_promoters_skipped": 0,
        "invalid_or_missing_sequence_skipped": 0,
        "position_parse_failures": 0,
    }
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_download_page, url, pages_dir) for url in page_urls]
        for future in as_completed(futures):
            url, payload, error = future.result()
            if error is not None or payload is None:
                failed_pages[url] = error or "unknown download error"
                continue
            page_name = _page_filename(url)
            page_hashes[page_name] = sha256_bytes(payload)
            records, counts = extract_promoter_rows(
                payload.decode("utf-8", errors="replace"),
                page_name=page_name,
            )
            all_records.extend(records)
            for key, value in counts.items():
                aggregate_counts[key] += value

    frame = pd.DataFrame(all_records)
    if frame.empty:
        raise RuntimeError("DBTBS download completed without any usable sigma-promoter rows")
    frame = frame.sort_values("sequence_id").reset_index(drop=True)
    metadata_path = Path(metadata_output)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(metadata_path, sep="\t", index=False)
    report = {
        "source_id": "dbtbs_v4.1",
        "source_version": "release 4.1",
        "access_date": date.today().isoformat(),
        "landing_page": urljoin(base_url, "home.html"),
        "index_url": index_url,
        "index_path": str(index_path.resolve()),
        "index_bytes": len(index_bytes),
        "index_sha256": sha256_bytes(index_bytes),
        "pages_requested": len(page_urls),
        "pages_downloaded": len(page_hashes),
        "pages_failed": failed_pages,
        "page_hashes": page_hashes,
        "raw_pages_dir": str(pages_dir.resolve()),
        "metadata_output": str(metadata_path.resolve()),
        "metadata_sha256": sha256_bytes(metadata_path.read_bytes()),
        "metadata_rows": len(frame),
        "sigma_counts": frame["sigma_factor_type"].value_counts().to_dict(),
        "parse_counts": aggregate_counts,
    }
    report_path = Path(report_output)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", default="data/raw/dbtbs_v4.1_20260916")
    parser.add_argument("--metadata-output", default="data/interim/dbtbs_metadata_20260916.tsv")
    parser.add_argument("--report-output", default="data/interim/dbtbs_download_report_20260916.json")
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--max-pages", type=int, default=None, help="Only for a bounded smoke test")
    args = parser.parse_args()
    report = download_and_parse(
        raw_dir=args.raw_dir,
        metadata_output=args.metadata_output,
        report_output=args.report_output,
        base_url=args.base_url,
        workers=args.workers,
        max_pages=args.max_pages,
    )
    print(f"pages_downloaded={report['pages_downloaded']}")
    print(f"pages_failed={len(report['pages_failed'])}")
    print(f"metadata_rows={report['metadata_rows']}")
    print(f"metadata_sha256={report['metadata_sha256']}")
    print(f"metadata_output={report['metadata_output']}")


if __name__ == "__main__":
    main()
