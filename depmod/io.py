"""IO Module.

Helper module that provide various I/O functions.

In particular functions to parse lammps file's header and extract the lattice vectors.
"""
from __future__ import annotations

import gzip
import re
import shutil
import tarfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import IO

import numpy as np

RE_COMPRESSED = re.compile(r"[.]gz")
RE_TARFILE = re.compile(r"[.]tar")
RE_TARGZ = re.compile(r"[.]tar|[.gz]")


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
def mkopen(filepath: str | Path, mode: str, **kwargs) -> IO:
    """Same as 'open()' but create the parent directory if it does not exists"""
    filepath = Path(filepath).absolute()
    mkparent(filepath)
    with open(filepath, mode=mode, **kwargs) as fd:
        yield fd


@dataclass
class SystemInfos:
    """
    Atomic system informations.

    Attributes
    ----------
    filepath : Path
        Path to the file.
    natom : int
        Number of atom in the system.
    cell : np.ndarray
        A 3x3 array that represent the lattice vector column wise.
    """

    filepath: Path
    natom: int
    cell: np.ndarray


def guess_atom_format(filepath: str | Path) -> str | None:
    """Guess the file format."""
    filepath = Path(filepath)

    try:
        fmt = Path(re.sub(r"[.]gz|[.]tar", "", filepath.name)).suffixes[0].replace(".", "")
    except IndexError:
        raise ValueError(f"Can't guess format from '{filepath.name}' (no suffix).") from None

    return fmt


def guess_filetype(filepath: str | Path) -> tuple[bool, bool]:
    """Guess the file type."""
    is_compressed = bool(RE_COMPRESSED.search(filepath.name))
    is_tarfile = bool(RE_TARFILE.search(filepath.name))
    return (is_compressed, is_tarfile)


def read_atom(
    filepath: str | Path,
    fmt: str | None = None,
    tar: bool | None = None,
    gz: bool | None = None,
    **kwargs,
) -> SystemInfos:
    """Read an atom file and extract informations on the system."""
    filepath = Path(filepath)

    if fmt is None:
        try:
            fmt = guess_atom_format(filepath)
        except ValueError as error:
            raise RuntimeError(
                "Could not guess the file format. " "Try to explicitely specify 'fmt=??'."
            ) from error

    try:
        reader = AtomFormat.get_reader(fmt)
    except ValueError as error:
        raise RuntimeError(
            f"'{fmt}' is not a supported file format."
            f"\nSupported file format are: {AtomFormat.list_readers()}"
        ) from error

    with open_atom_file(filepath, tar=tar, gz=gz) as fd:
        natom, cell = reader(fd, **kwargs)

    return SystemInfos(filepath, natom, cell)


@contextmanager
def open_atom_file(filepath: str | Path, tar: bool | None = None, gz: bool | None = None) -> IO:
    """Context manager for opening file."""
    filepath = Path(filepath)

    # check if the file exists
    if not filepath.is_file():
        raise FileNotFoundError(f"{filepath}")

    # check if the file is a tar and compressed
    is_compressed, is_tarfile = guess_filetype(filepath)

    # overrided by function argument
    if tar is not None:
        is_tarfile = tar
    if gz is not None:
        is_compressed = gz

    if is_compressed and is_tarfile:
        with tarfile.open(filepath, "r|gz") as fd:
            yield fd.extractfile(fd.next())
    elif is_compressed:
        with gzip.open(filepath, "rb") as fd:
            yield fd
    else:
        with open(filepath, "rb") as fd:
            yield fd


class AtomFormat:
    """Static class to manage reading functions."""

    __readers: dict = {}
    __labels: dict = {}

    @classmethod
    def reader(cls, extention: list[str], fmt_label) -> callable:
        """Decorate a function to be associate to a file format."""

        def wrapper(func):
            cls.__labels[fmt_label] = []
            for ext in extention:
                cls.__readers[ext] = func
                cls.__labels[fmt_label].append(ext)
            return func

        return wrapper

    @classmethod
    def get_reader(cls, fmt: str) -> callable:
        """Return the appropriate reading function from the given format."""
        try:
            reader_func = cls.__readers[fmt]
        except KeyError:
            raise ValueError(f"Uknown file format '{fmt}'") from None

        return reader_func

    @classmethod
    def list_readers(cls):
        """Return the file extension associated to registered formats."""
        return cls.__labels


@AtomFormat.reader(["lmp", "lmpdata", "lmp-data", "lammps-data"], "lammps data")
def read_lammps_data(fd) -> tuple[int, np.ndarray]:
    """Reader for lammps data file."""
    re_strnum = re.compile(rb"([+-]?\d*[.]?\d*[Eefd]?[+-]?\d+)")
    re_nat = re.compile(rb"\s*\d+\s+atoms")
    re_box = re.compile(rb"(.*xlo\s+xhi|.*ylo\s+yhi|.*zlo\s+zhi)")
    re_axis = re.compile(rb"\wlo")
    re_tilt = re.compile(rb"xy|xz|yz")
    re_section = re.compile(rb"(\s*Masses|\s*Atoms)")

    def __parse_natom(line):
        if re_nat.search(line):
            natom = int(re_strnum.search(line).group())
            return (True, "nat", natom)
        return (False, None, None)

    def __parse_axis(line):
        if re_box.search(line):
            axis = re_axis.search(line).group().decode()[0]
            xlo, xhi = map(float, re_strnum.findall(line))
            return (True, axis, xhi - xlo)
        return (False, None, None)

    def __parse_tilt(line):
        if re_tilt.search(line):
            order = {k: i for i, k in enumerate(["xy", "xz", "yz"])}
            axis = map(lambda x: x.decode(), re_tilt.findall(line))
            num = map(float, re_strnum.findall(line))
            tilt = [0.0] * 3
            for i, j in zip(axis, num, strict=True):
                tilt[order[i]] = j
            return (True, "tilt", tilt)
        return (False, None, None)

    def __parse_stop_flag(line):
        if re_section.search(line):
            return (True, "stop", True)
        return (False, None, None)

    _check = {"natom", "x", "y", "z", "tilt", "stop"}

    _parser = {
        "natom": __parse_natom,
        "x": __parse_axis,
        "y": __parse_axis,
        "z": __parse_axis,
        "stop": __parse_stop_flag,
        "tilt": __parse_tilt,
    }

    infos = {"stop": False}

    try:
        while True:
            line = next(fd)

            for key in _check:
                parsed, ikey, value = _parser[key](line)

                if parsed:
                    infos[ikey] = value
                    _check.remove(key)
                    break

                if infos["stop"]:
                    break

    except StopIteration:
        pass

    for k in ("nat", "x", "y", "z"):
        if not infos.get(k):
            raise KeyError(f"Fail to parse lammps data file. Field '{k}' is missing.")

    if "tilt" in infos:
        for i, k in enumerate(("xy", "xz", "yz")):
            infos[k] = infos["tilt"][i]

    for k in ("xy", "xz", "yz"):
        if not infos.get(k):
            infos[k] = 0.0

    cell = np.zeros((3, 3), dtype=float)

    cell[0, 0] = infos["x"]
    cell[1, 1] = infos["y"]
    cell[2, 2] = infos["z"]

    cell[0, 1] = infos["xy"]
    cell[0, 2] = infos["xz"]
    cell[1, 2] = infos["yz"]

    return (infos["nat"], cell)


@AtomFormat.reader(["dump", "lmpdump", "lmp-dump", "lammps-dump"], "lammps dump")
def read_lammps_dump(fd) -> tuple[int, np.ndarray]:
    pass
