import re
import sys
import glob

import setuptools
from setuptools import setup, Extension

class get_pybind_include():
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

depmod_cpp_module = Extension(
    "depmod._lib",
    glob.glob("depmod/srclib/*.cpp"),
    include_dirs=[
        "depmod/lib/",
        "external/eigen/",
        get_pybind_include(),
        get_pybind_include(True),
    ],
    language="c++",
    extra_compile_args=["-std=c++14", "-O3", "-flto"],
    extra_link_args=[],
)

if __name__ == "__main__":
    setup(ext_modules=[depmod_cpp_module])
