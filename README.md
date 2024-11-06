# DEformation Paths for MOlecular Dynamics (DEPMOD)

## Installation
```sh
# --recursive to ensure third party libraries are properly cloned
git clone --recursive  https://github.com/lafourcadep/depmod.git
```
```sh
# Installing the package using pip
pip install ./depmod
# To be able to run the tests dev optional modules are required
pip install ./depmod[dev]
```
### Special case when using a proxy
```sh
# Pip can not propagate the --proxy variable to subprocess used for build system.
# It is therefore needed to install build dependencies manually.
pip install -U pip setuptools wheel pybind11 --proxy=https://proxy:0000
# Then install the package
pip install ./depmod[dev] --no-build-isolation --proxy=https://proxy:0000
```
