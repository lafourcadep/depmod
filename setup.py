import glob

from setuptools.command import build_ext
import setuptools


class get_pybind_include:
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11

        return pybind11.get_include(self.user)


# class custom_build_ext(build_ext):
#     def build_extensions(self):
#         # Override the compiler executables. Importantly, this
#         # removes the "default" compiler flags that would
#         # otherwise get passed on to to the compiler, i.e.,
#         # distutils.sysconfig.get_var("CFLAGS").
#         self.compiler.set_executable("compiler_so", "c++")
#         self.compiler.set_executable("compiler_cxx", "c++")
#         self.compiler.set_executable("linker_so", "c++")
#         build_ext.build_extensions(self)


depmod_cpp_module = setuptools.Extension(
    "depmod._lib",
    glob.glob("src/cpp/*.cpp"),
    include_dirs=[
        "src/cpp/",
        "external/eigen/",
        get_pybind_include(),
        get_pybind_include(True),
    ],
    language="c++",
    extra_compile_args=[
        "-std=c++2a",
        "-O3",
        "-flto",
        "-fPIC",
        "-Wall",
        "-Wextra",
        "-Wattributes",
    ],
    extra_link_args=[],
)


if __name__ == "__main__":
    setuptools.setup(
        ext_modules=[depmod_cpp_module],
        zip_safe=False,
    )
