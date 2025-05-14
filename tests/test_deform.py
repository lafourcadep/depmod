import numpy as np
import pytest

from depmod.deformation import Compression, PureShear, Traction, TractionCompression


class Test_python_deformation:
    def test_TractionCompression(self):
        deform1 = TractionCompression.from_axis([1, 0, 0])
        assert np.allclose(np.array(deform1.S(), dtype=float), np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=float))
        deform2 = TractionCompression.from_axis([2, 0, 0])
        assert np.allclose(np.array(deform2.S(), dtype=float), np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=float))
        assert np.allclose(TractionCompression.from_axis([1, 1, 2]).S(), TractionCompression.from_axis([2, 2, 4]).S())
        assert np.allclose(TractionCompression.from_angles(0, 0).S(), np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=np.float64))

    def test_Compression(self):
        z = [0, 0, 1]
        assert np.allclose(Compression.from_axis(z).S(), np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=np.float64))
        assert np.allclose(Compression.from_angles(0.0, 0.0).S(), np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=np.float64))
        assert np.allclose(Compression.from_axis(z, isoV = True).S(), np.array([[-0.5, 0, 0], [0, -0.5, 0], [0, 0, 1]], dtype=np.float64))
        assert np.allclose(Compression.from_angles(0.0, 0.0, isoV = True).S(), np.array([[-0.5, 0, 0], [0, -0.5, 0], [0, 0, 1]], dtype=np.float64))

    def test_Traction(self):
        z = [0, 0, 1]
        assert np.allclose(Traction.from_axis(z).S(), np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=np.float64))
        assert np.allclose(Traction.from_angles(0.0, 0.0).S(), np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=np.float64))
        assert np.allclose(Traction.from_axis(z, isoV = True).S(), np.array([[-0.5, 0, 0], [0, -0.5, 0], [0, 0, 1]], dtype=np.float64))
        assert np.allclose(Traction.from_angles(0.0, 0.0, isoV = True).S(), np.array([[-0.5, 0, 0], [0, -0.5, 0], [0, 0, 1]], dtype=np.float64))

    def test_PureShear(self):
        assert np.allclose(PureShear.from_angles(0.0, 0.0, degree=True).S(), np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0]], dtype=np.float64))
        assert np.allclose(PureShear.from_angles(0.0, 0.0, symmetric=True, degree=True).S(), np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=np.float64))
        assert np.allclose(PureShear.from_axis([1, 0, 0], [0, 0, 1]).S(), np.array([[0, 0, 1], [0, 0, 0], [0, 0, 0]], dtype=np.float64))
        assert np.allclose(PureShear.from_axis([3, 0, 0], [0, 0, 2]).S(), np.array([[0, 0, 1], [0, 0, 0], [0, 0, 0]], dtype=np.float64))
        # Should raise error with non orthogonal vectors
        with pytest.raises(ValueError):
            PureShear.from_axis([1, 0, 0], [1, 0, 0], strict = True)
