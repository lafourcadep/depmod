import time

from depmod.config import init_config
from depmod.deformation import Traction
from depmod.io import read_atom

from depmod.lammps import (
    lammps_generate_box_evolution_data,
    lammps_write_fix_deform_module
)


if __name__ == "__main__":
    
    # define the deformation
    # unixial compression in the x direction
    deform = Traction(isoV=True).from_axis([1, 0, 0])

    # define the simulation parameters
    config = init_config(
        gammadot=1e9,
        t_max=200e-12,
        npts=250,
        kpts=1000
    )

    # Read lattice from lammps file
    atoms = read_atom("dir.input/Cu_fcc_5x5x5.lmp.gz")

    # Generate box parameters time evolution
    box_evolution_data = lammps_generate_box_evolution_data(
        atoms.cell,
        config,
        deform,
        method="brute",
        filename="box_evolution_data.csv"
    )

    # Write lammps deformation module
    lammps_write_fix_deform_module(
        box_evolution_data,
        filename="lmp_fix_deform.mod",
        dumpfit="box_evolution_fit.csv"
    )
