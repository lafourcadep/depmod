import glob
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class get_pybind_include:
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11

        return pybind11.get_include(self.user)


class CustomBuildExt(build_ext):
    def build_extensions(self):
        kw = {"zlib": "-lz", "bzip2": "-lbz2", "lzma": "-llzma"}

        for k, v in kw.items():
            if subprocess.call(["bash", f"./scripts/check_{k}.sh"]) == 0:
                print(f"{k.upper()} found. Enabling {k} support.")
                for ext in self.extensions:
                    ext.define_macros.append((f"USE_{k.upper()}", "1"))
                    ext.extra_link_args.append(v)
            else:
                print(f"{k.upper()} not found. Building without {k} support.")

        super().build_extensions()


depmod_cpp_module = Extension(
    "depmod._lib",
    glob.glob("src/cpp/*.cpp"),
    include_dirs=["src/cpp/", "external/eigen/", get_pybind_include(), get_pybind_include(True)],
    language="c++",
    extra_compile_args=["-std=c++2a", "-O3", "-flto=auto", "-fPIC", "-Wall", "-Wextra", "-Wattributes"],
    extra_link_args=[],
)


if __name__ == "__main__":
    setup(ext_modules=[depmod_cpp_module], zip_safe=False, cmdclass={"build_ext": CustomBuildExt})
