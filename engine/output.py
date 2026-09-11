"""Where the scripts write figures and tables (next to the scripts, not the cwd)."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def figure_path(name):
    d = ROOT / "figures"
    d.mkdir(exist_ok=True)
    return d / name


def table_path(name):
    d = ROOT / "tables"
    d.mkdir(exist_ok=True)
    return d / name
