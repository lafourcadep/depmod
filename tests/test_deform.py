import pytest

import numpy as np

from depmod.deformation import (
    PureShear,
    TractionCompression,
    Traction,
    Compression
)

from depmod._lib.deformation import (
    _lib_Deformation,
    _lib_TractionCompression,
    _lib_TractionCompressionIsoV,
    _lib_PureShear
)

class Test_python_deformation:

    def test_TractionCompression_flags(self):
        a = TractionCompression()
        assert (a.comp, a.isoV) == (False, False)
        
        for flags in [(False, False), (True, False), (False, True), (True, True)]:
            a = TractionCompression(*flags)
            assert (a.comp, a.isoV) == flags

        assert isinstance(TractionCompression(isoV=False).from_axis([1, 0, 0]), _lib_TractionCompression)
        assert isinstance(TractionCompression(isoV=True).from_axis([1, 0, 0]), _lib_TractionCompressionIsoV)

    def test_TractionCompression_from_axis(self):
        deform1 = TractionCompression().from_axis([1, 0, 0])
        assert np.allclose(deform1.S, np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=float))        

        deform2 = TractionCompression().from_axis([2, 0, 0])
        assert np.allclose(deform2.S, np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=float))        

        deform3 = TractionCompression().from_axis([1, 1, 2])
        deform4 = TractionCompression().from_axis([2, 2, 4])

        assert np.allclose(deform3.S, deform4.S)

    def test_TractionCompression_from_angles(self):
        deform1 = TractionCompression().from_angles(0, 0)
        assert np.allclose(deform1.S, np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=float))        

    def test_Compression(self):
        z = [0, 0, 1]

        assert isinstance(Compression().from_axis(z), _lib_TractionCompression)
        assert Compression().from_axis(z).flags() == (True, False)

        assert isinstance(Compression(isoV=True).from_axis(z), _lib_TractionCompressionIsoV)
        assert Compression(isoV=True).from_axis(z).flags() == (True, True)
    
    def test_Traction(self):
        z = [0, 0, 1]

        assert isinstance(Traction().from_axis(z), _lib_TractionCompression)
        assert Traction().from_axis(z).flags() == (False, False)

        assert isinstance(Traction(isoV=True).from_axis(z), _lib_TractionCompressionIsoV)
        assert Traction(isoV=True).from_axis(z).flags() == (False, True)

    def test_PureShear_from_angles(self):
        deform1 = PureShear().from_angles(0, 0, degree=True)
        assert np.allclose(deform1.S, np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0]], dtype=float))        

        deform2 = PureShear().from_angles(45, 45, degree=True)
        deform3 = PureShear().from_angles(np.deg2rad(45), np.deg2rad(45), degree=False)
        assert np.allclose(deform2.S, deform3.S)

        # Should raise S property is readonly
        with pytest.raises(AttributeError):
            deform3.S = np.eye(3)


    def test_PureShear_from_axis(self):
        deform1 = PureShear().from_axis([1, 0, 0], [0, 0, 1])
        assert np.allclose(deform1.S, np.array([[0, 0, 1], [0, 0, 0], [0, 0, 0]], dtype=float))        

        # check if vector are properly normalized
        deform2 = PureShear().from_axis([3, 0, 0], [0, 0, 2])
        assert np.allclose(deform2.S, np.array([[0, 0, 1], [0, 0, 0], [0, 0, 0]], dtype=float))

        # Should raise error with non orthogonal vectors
        with pytest.raises(ValueError):
            deform = PureShear().from_axis([1, 0, 0], [1, 0, 0])
