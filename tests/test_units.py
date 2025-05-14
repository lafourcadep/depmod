import pytest

from depmod.units import UT as u
from depmod.units import (
        UnitDimension,
        UnitDimensionError,
        InvalidUnitSymbolError,
        convert,
        lammps_unit,
)

def test_units_access():
    u["m"]
    # Access with unregistered symbols
    with pytest.raises(InvalidUnitSymbolError):
        u["O"]

def test_units_operations():
    m = u["m"]

    # Addition
    assert (10 * m + 5 * m).val == 15.0
    assert (10 * m + 5 * m).dim == UnitDimension(L=1)

    with pytest.raises(UnitDimensionError):
        u["m"] + u["kg"]
    
    # Substraction
    assert (10 * m - 5 * m).val == 5.0
    assert (10 * m - 5 * m).dim == UnitDimension(L=1)
    
    with pytest.raises(UnitDimensionError):
        u["m"] - u["kg"]
    
    # Multiplication
    assert (10 * m).val == 10.0
    assert (10 * m).dim == UnitDimension(L=1)
    
    assert (3 * m * m).val == 3.0
    assert (3 * m * m).dim == UnitDimension(L=2)

    # Division
    assert (m / 10).val == 0.1
    assert (m / 10).dim == UnitDimension(L=1)
    
    assert (3 * m / m).val == 3.0
    assert (3 * m / m).dim == UnitDimension()
    
    # Power
    assert (3 * m ** 2).val == 3.0
    assert (3 * m ** 2).dim == UnitDimension(L=2)
    
    # Units
    assert (10 * u["kg"] / (u["m"] ** 3)).val == 10.0
    assert (10 * u["kg"] / (u["m"] ** 3)).dim == UnitDimension(M=1, L=-3)


def test_units_conversion():
    pass

def test_units_lammps():

    assert lammps_unit("time", "metal") == u["ps"]
    assert convert( lammps_unit("time", "metal"), "ns") == 1.0e-3
    assert convert( lammps_unit("time", "real"), "ns") == 1.0e-6
