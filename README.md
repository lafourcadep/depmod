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
# pip can not propagate the --proxy variable to subprocess used for build system.
# It is therefore needed to install build dependencies manually.
pip install -U pip setuptools wheel pybind11 --proxy=https://proxy:0000
# Then install the package
pip install ./depmod[dev] --no-build-isolation --proxy=https://proxy:0000
```

## Usage

### Script
The following code can be used to generante LAMMPS input file of a volume conserving traction
```python
# example/traction_isoV/main.py
from depmod import DeformationPath, Traction, read_lattice
from depmod.core import lmp

if __name__ == "__main__":

    lmp.write_deformation_module(
        # Generate the time evolution of the box parameter
        lmp.generate_box_evolution_data(
            # read H(0) from lammps file
            lattice=read_lattice("./dir.input/Cu_fcc_5x5x5.lmp.gz"),
            # Create a deformation path: Uni-axial isochoric traction
            deformation=DeformationPath(
                Traction.from_axis([1, -1, 2], isoV=True),
                strain_rate=1e10,
                tmax=200e-12,
                npts=200
            ),
            filename="dir.depmod/box_evolution_data.csv",
        ),
        units="metal",
        mod_file="dir.depmod/lmp_fix_deform.mod",
        fit_file="dir.depmod/box_evolution_fit.csv",
        poly_file="dir.depmod/poly_coeffs.csv",
    )
```
To apply the deformation path during a molecular dynamics simulation, just add `include dir.depmod/lmp_fix_deform.mod` in your LAMMPS input script.

<!-- ### Command Line Utility -->

<!-- Deformation path can also be generate directly from the terminal by using the `depmod` command line utility. The following code will produce the same output as above. -->
<!-- ```sh -->
<!-- depmod deform --lmp -units metal -traction ax 1-12 -rate 1e10 -tmax 200e-12 -npts 200 --lattice "./dir.input/Cu_fcc_5x5x5.lmp.gz" -o lmp_fix_deform.mod -->
<!-- ``` -->
<!-- Note that only the `lmp_fix_deform.mod` will be generated. -->
