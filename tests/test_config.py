import pytest

from depmod.config import init_config

from depmod._lib.config import _lib_Configuration

def test_params_assert_init():

    config = init_config(
        gammadot=1e9,
        t_min=0,
        t_max=100e-12,
        npts=100,
        kpts=1
    )

    assert isinstance(config, _lib_Configuration)
