import warnings
warnings.filterwarnings("ignore")

from depmod.core import *
from depmod.config import init_config
from depmod.deformation import Traction
from depmod.io import read_atom
from depmod.units import convert

if __name__ == "__main__":

    # define the deformation
    # unixial traction in the x direction
    deform = Traction(isoV=True).from_axis([1, 0, 0])

    # define the simulation parameters
    config = init_config(
        gammadot= 1e10,
        t_max=200e-12,
        npts=200,
        kpts=5000,
    )

    # Read lattice from lammps file
    atoms = read_atom("dir.input/Cu_fcc_5x5x5.lmp.gz")

    # Generate box parameters time evolution
    box_evolution_data = lammps_generate_box_evolution_data(
        atoms.cell,
        config,
        deform,
        method="brute",
        filename="dir.depmod/box_evolution_data.csv"
    )

    # Write lammps deformation module
    lammps_write_deformation_module(
        box_evolution_data,
        mod_file="dir.depmod/lmp_fix_deform.mod",
        fit_file="dir.depmod/box_evolution_fit.csv",
        poly_file="dir.depmod/poly_coeffs.csv"
    )

    # Write exastamp deformation module
    # exastamp_write_deformation_module(
    #     exastamp_generate_box_evolution_data(
    #         config = config,
    #         deform = deform
    #     )
    # )
