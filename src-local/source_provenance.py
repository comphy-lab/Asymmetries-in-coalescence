"""Truthful source identity for contour-run metadata."""
from pathlib import Path
import subprocess


def source_commit(root: Path) -> str:
    """Return HEAD, adding -dirty or -unverified when appropriate."""
    try:
        head = subprocess.run(
            ["git", "rev-parse", "--verify", "HEAD"], cwd=root,
            capture_output=True, text=True, check=False,
        )
        if head.returncode != 0 or not head.stdout.strip():
            return "unavailable"
        commit = head.stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=normal"],
            cwd=root, capture_output=True, text=True, check=False,
        )
    except OSError:
        return "unavailable"
    if status.returncode != 0:
        return commit + "-unverified"
    return commit + ("-dirty" if status.stdout else "")
