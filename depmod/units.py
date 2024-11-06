"""Module to manage units conversion.

Mainly because LAMMPS might use several internal unit systems.
See: https://docs.lammps.org/units.html

Might be overkill... but whatever.
"""
#TODO: Might go for a simpler version ??
from __future__ import annotations

import warnings

from collections import namedtuple

import numpy as np


class UnitDimensionError(Exception):
    pass


class InvalidUnitSymbolError(Exception):
    pass


class UnitDimension:
    """Respresent the dimension of a Unit in the S.I. basis."""

    def __init__(
        self,
        T: int=0, # s: Time
        L: int=0, # m: Length
        M: int=0, # kg: Mass
        A: int=0, # A: Eletric Current
        K: int=0, # K: Temperature
        N: int=0, # mol: Matter Quantity
        J: int=0  # Cd: Light
    ) -> None:

        # represent the exponent associated to each s.i. unit
        self.dim = np.array([T, L, M, A, K, N, J], dtype=int)

    def copy(self):
        return self.__class__(*self.dim)
    
    def __hash__(self):
        return hash(tuple(self.dim))

    def __str__(self):
        _str = []
        _tmp = []
        for i, sym in zip(
            [2, 1, 0, 4, 5, 3, 6],
            ["kg", "m", "s", "K", "mol", "A", "Cd"]
        ):
            dim = self.dim[i]
            if dim == 0:
                continue
            elif dim == 1:
                _str.append("%s")
                _tmp.append(sym)
            elif (dim < 0 or dim > 1):
                _str.append("%s^%d")
                _tmp += [sym, dim]
        return (".".join(_str) % tuple(_tmp))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return np.all(self.dim == other.dim)

    def __add__(self, other):
        if self == other:
            return self.__class__(*self.dim)
        raise UnitDimensionError()

    def __sub__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        return self.__class__(*(self.dim + other.dim))

    def __rmul__(self, other):
        return self.__class__(*(other.dim + self.dim))
    
    def __truediv__(self, other):
        return self.__class__(*(self.dim - other.dim))

    def __rtruediv__(self, other):
        return self.__class__(*(other.dim - self.dim))

    def __pow__(self, other: int | float):
        if not isinstance(other, (int, float)):
            raise UnitDimensionError()
        return self.__class__(*(self.dim * other))


def cast_other_as_unit(func):
    def wrapper(self, value: int | float | Unit):
        if not issubclass(type(value), Unit):
            value = Unit(value)
        return func(self, value)
    return wrapper


class Unit:
    def __init__(self, value: float = 0., dim: UnitDimension | None = None):
        self.val = value
        self.dim = UnitDimension() if dim is None else dim

    @cast_other_as_unit
    def __eq__(self, other):
        return (self.val == other.val) and (self.dim == other.dim)

    def __hash__(self):
        return hash((self.val, hash(self.dim)))

    @cast_other_as_unit
    def __add__(self, other):
        return self.__class__(
            self.val + other.val,
            self.dim + other.dim
        )
    
    @cast_other_as_unit
    def __sub__(self, other):
        return self.__class__(
            self.val - other.val,
            self.dim - other.dim
        )

    @cast_other_as_unit
    def __mul__(self, other):
        return self.__class__(
            self.val * other.val,
            self.dim * other.dim
        )
    
    @cast_other_as_unit
    def __rmul__(self, other):
        return other.__mul__(self)
    
    @cast_other_as_unit
    def __truediv__(self, other):
        return self.__class__(
            self.val / other.val,
            self.dim / other.dim
        )
    
    @cast_other_as_unit
    def __rtruediv__(self, other):
        return other.__truediv__(self)
    
    def __pow__(self, other):
        return self.__class__(
            self.val ** other,
            self.dim ** other
        )

    def __str__(self):
        return f"Unit({self.val}, [{self.dim}])"

    def __repr__(self):
        return str(self)

    def __int__(self):
        return int(self.val)
    
    def __float__(self):
        return self.val


UnitDefinition = namedtuple("UnitDefinition", ["unit", "symbols", "name", "infos"])


class UnitContainer:
    def __init__(self):
        self.dct = {}
        self.sym = {}

    def register(self, udef: UnitDefinition):
        ukey = hash(udef)

        if ukey in self.dct:
            return False

        self.dct[ukey] = udef

        for symb in udef.symbols:
            self.sym[symb] = ukey

    def get(self, key):
        return self.dct.get(self.sym.get(key))

    def __getitem__(self, key: str):
        try:
            return self.dct[self.sym[key]]
        except KeyError:
            raise NotImplementedError(f"Symbol '{key}' is not associated to any units.")


class UnitTable:
    def __init__(self):
        self.__units = UnitContainer()
        self.__constants = UnitContainer()

    def __register(self, container, symbols, unit, name, infos) -> Unit:
        if not isinstance(symbols, (list, tuple)):
            symbols = [symbols]

        _symbols = []
        for symbol in symbols:
            if self.__units.get(symbol) is not None:
                warnings.warn(f"{self.__class__.__name__}: \"{symbol}\" already registered. Skipped.")
                continue
            _symbols.append(symbol)

        container.register(UnitDefinition(unit, tuple(_symbols), name, infos))

        return unit

    def create_unit(self, symbols, unit, name, infos = None) -> Unit:
        return self.__register(
            self.__units,
            symbols,
            unit,
            name,
            infos
        )

    def create_constant(self, symbols, unit, name, infos = None) -> Unit:
        return self.__register(
            self.__constants,
            symbols,
            unit,
            name,
            infos
        )

    
    @property
    def cte(self):
        return self.__constants

    def __getitem__(self, key):
        try:
            return self.__units[key].unit
        except NotImplementedError as err:
            raise InvalidUnitSymbolError from err


__codata_version__ = "2018"


CODATA = {
    "2018" : {

    }
}


MASS = UnitDimension(M=1)
LENGTH = DISTANCE = UnitDimension(L=1)
TIME = UnitDimension(T=1)
TEMPERATURE = UnitDimension(K=1)
MOLE = UnitDimension(N=1)
SURFACE = LENGTH * LENGTH
VOLUME = SURFACE * LENGTH
DENSITY = MASS / VOLUME
VELOCITY = DISTANCE / TIME
FORCE = MASS * LENGTH / TIME ** 2
ENERGY = FORCE * LENGTH
PRESSURE = FORCE / SURFACE

_QUANTITY = {
    "mass": MASS,
    "distance": LENGTH,
    "time": TIME,
    "energy": ENERGY,
    "velocity": VELOCITY,
    "force": ENERGY,
    "temperature": TEMPERATURE,
    "pressure": PRESSURE,
    "density": DENSITY
}

_RQUANTITY = {v: k for k, v in _QUANTITY.items()}


def create_unit_table(codata_version):
    ut = UnitTable()

    # Length units
    m = ut.create_unit("m", Unit(1., LENGTH), "metre")

    # Mass units
    kg = ut.create_unit("kg", Unit(1., MASS), "kilogram")
    g = ut.create_unit("g", 1.0e-03 * kg, "gram")

    # Time units
    s = ut.create_unit("s", Unit(1.0, TIME), "second")
    
    ut.create_unit("ns", 1.0e-09 * s, "nanosecond")
    ut.create_unit("ps", 1.0e-12 * s, "picosecond")
    ut.create_unit("fs", 1.0e-15 * s, "femtosecond")

    return ut


UT = create_unit_table(__codata_version__)


def convert(unit1: str | Unit, unit2: str | Unit) -> float:
    unit1 = UT[unit1] if isinstance(unit1, str) else unit1
    unit2 = UT[unit2] if isinstance(unit2, str) else unit2

    if unit1.dim != unit2.dim:
        raise UnitDimension

    return float(unit1 / unit2)


# ------- LAMMPS CONVERSIONT TABLES


LAMMPS_UNIT_URL = "https://docs.lammps.org/units.html"

LAMMPS_UNIT_SET = {
    "metal": {
        TIME: UT["ps"]
    },
    "real": {
        TIME: UT["fs"]
    }
}

class InvalidLammpsUnitError(Exception):
    def __init__(self, unit_system: str) -> None:
        super().__init__(
            f"Invalid LAMMPS unit system '{unit_system}'\n"
            f"Supported unit system currently are {tuple(LAMMPS_UNIT_SET.keys())}\n"
            f"For more informations see: {LAMMPS_UNIT_URL}"
        )


class InvalidLammpsQuantityError(Exception):
    def __init__(self, quantity: str) -> None:
        super().__init__(
            f"Invalid LAMMPS quantity '{quantity}'\n"
            f"Supported quantity are {tuple(_QUANTITY.keys())}\n"
            f"For more informations see: {LAMMPS_UNIT_URL}"
        )


def lammps_unit(quantity: str, usys: str = "metal"):
    try:
        unit_set = LAMMPS_UNIT_SET[usys]
    except KeyError as err:
        raise InvalidLammpsUnitError(str(usys)) from err

    try:
        return unit_set[_QUANTITY[quantity]]
    except KeyError as err:
        raise InvalidLammpsQuantityError(str(quantity)) from err
