"""Notebook helpers that keep the workshop usable in Colab and locally.

Import :func:`install_colab_shims` in a notebook before any optional Colab
imports.  On Colab it does nothing.  In Jupyter it provides small equivalents
for the Drive mount and file-upload interfaces used by the workshop notebooks.
"""

from __future__ import annotations

import os
import sys
import types
from pathlib import Path


def is_colab() -> bool:
    """Return whether this kernel is running in Google Colab."""
    try:
        import google.colab  # noqa: F401
    except ImportError:
        return False
    return True


def project_root() -> Path:
    """Locate the repository root even when Jupyter starts in a subfolder."""
    for candidate in (Path.cwd(), *Path.cwd().parents):
        if (candidate / "utils").is_dir():
            return candidate
    return Path.cwd()


def shared_cache_root(workshop_name: str = "NOAA-workshop2") -> str:
    """Return the shared Drive location in Colab or a persistent local cache."""
    if is_colab():
        return f"/content/drive/Shareddrives/{workshop_name}"

    root = Path(os.environ.get("NOAA_WORKSHOP_DATA_DIR", project_root() / ".noaa-workshop-cache"))
    root.mkdir(parents=True, exist_ok=True)
    return str(root)


def content_path(name: str) -> str:
    """Map a Colab ``/content`` file to a project-local file outside Colab."""
    if is_colab():
        return f"/content/{name}"
    local_content = project_root() / ".noaa-workshop-cache" / "content"
    local_content.mkdir(parents=True, exist_ok=True)
    return str(local_content / name)


def install_colab_shims() -> None:
    """Install local fallbacks for the limited ``google.colab`` APIs we use."""
    if is_colab():
        return

    # Do not overwrite a real package if a user has one installed.
    if "google.colab" in sys.modules:
        return

    colab = types.ModuleType("google.colab")
    drive = types.ModuleType("google.colab.drive")
    files = types.ModuleType("google.colab.files")

    def mount(_mountpoint: str, force_remount: bool = False) -> None:
        del force_remount
        print(f"Local Jupyter detected; using cache at {shared_cache_root()} (Google Drive is not mounted).")

    def upload() -> dict[str, bytes]:
        """Prompt for a local path and mirror Colab's ``files.upload`` result."""
        source = input("Path of file to upload (leave blank to cancel): ").strip()
        if not source:
            return {}
        path = Path(source).expanduser()
        if not path.is_file():
            print(f"File not found: {path}")
            return {}
        destination = Path.cwd() / path.name
        if destination.resolve() != path.resolve():
            destination.write_bytes(path.read_bytes())
        return {path.name: destination.read_bytes()}

    drive.mount = mount
    files.upload = upload
    colab.drive = drive
    colab.files = files
    sys.modules["google.colab"] = colab
    sys.modules["google.colab.drive"] = drive
    sys.modules["google.colab.files"] = files
