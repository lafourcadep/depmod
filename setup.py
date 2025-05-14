import glob

import setuptools


class get_pybind_include:
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11

        return pybind11.get_include(self.user)


depmod_cpp_module = setuptools.Extension(
    "depmod._lib",
    glob.glob("src/cpp/*.cpp"),
    include_dirs=["src/cpp/", "external/eigen/", get_pybind_include(), get_pybind_include(True)],
    language="c++",
    extra_compile_args=["-std=c++2a", "-O3", "-flto", "-fPIC", "-Wall", "-Wextra", "-Wattributes"],
    extra_link_args=[],
)


if __name__ == "__main__":
    setuptools.setup(ext_modules=[depmod_cpp_module], zip_safe=False)
