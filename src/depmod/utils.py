from __future__ import annotations

import shutil
from contextlib import contextmanager
from pathlib import Path


def rmtree(path: str | Path) -> None:
    """Remove a file/directory."""
    path = Path(path).absolute()
    if path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def mkdir(path: str | Path, overwrite: bool = False) -> None:
    """Create a directory and its parent if does not exist."""
    path = Path(path).absolute()
    if path.is_dir() and overwrite:
        rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def mkparent(filepath: str | Path) -> None:
    """Create the parent directory of a file."""
    parent = Path(filepath).absolute().parent
    if not parent.is_dir():
        mkdir(parent)


@contextmanager
def mkopen(filepath: str | Path, mode: str, **kwargs):
    """Same as 'open()' but create the parent directory if it does not exists"""
    filepath = Path(filepath).absolute()
    mkparent(filepath)
    with open(filepath, mode=mode, **kwargs) as fd:
        yield fd


__all__ = ["rmtree", "mkdir", "mkparent", "mkopen"]
