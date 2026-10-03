"""Check Python packages and MEME Suite tools for PR01-02 motif analysis."""

from __future__ import annotations

import argparse
from importlib import metadata
import json
import platform
from pathlib import Path
import shutil
import subprocess
import sys


PYTHON_PACKAGES = (
    "matplotlib",
    "numpy",
    "pandas",
    "pyyaml",
    "pytest",
)

MEME_TOOLS = ("meme", "streme", "dreme", "fimo", "tomtom")
LOCAL_MEME_BIN = Path(__file__).resolve().parents[1] / "tmp" / "meme-suite-5.5.9" / "bin"


def _package_version(package: str) -> str | None:
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def _tool_version(executable: str) -> str | None:
    for argument in ("--version", "-version"):
        try:
            completed = subprocess.run(
                [executable, argument],
                check=False,
                capture_output=True,
                text=True,
                timeout=10,
            )
        except (OSError, subprocess.TimeoutExpired):
            continue
        output = (completed.stdout or completed.stderr).strip()
        if completed.returncode == 0 and output:
            return output.splitlines()[0]
    return None


def _tool_path(tool: str) -> str | None:
    system_path = shutil.which(tool)
    if system_path:
        return system_path
    local_path = LOCAL_MEME_BIN / tool
    return str(local_path) if local_path.is_file() else None


def collect_environment() -> dict:
    """Collect a JSON-serializable experiment environment snapshot."""

    tools = {}
    for tool in MEME_TOOLS:
        executable = _tool_path(tool)
        tools[tool] = {
            "available": executable is not None,
            "path": executable,
            "version": _tool_version(executable) if executable else None,
        }
    package_versions = {package: _package_version(package) for package in PYTHON_PACKAGES}
    return {
        "python": {
            "executable": sys.executable,
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
        },
        "platform": platform.platform(),
        "packages": package_versions,
        "meme_suite": {
            "ready": all(tools[tool]["available"] for tool in ("meme", "streme", "fimo", "tomtom")),
            "tools": tools,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", help="Optional JSON output path")
    args = parser.parse_args()
    snapshot = collect_environment()
    payload = json.dumps(snapshot, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
