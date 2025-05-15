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
                Traction.from_axis([1, 0, 0], isoV=True), strain_rate=1e10, tmax=200e-12, npts=200
            ),
            filename="dir.depmod/box_evolution_data.csv",
        ),
        mod_file="dir.depmod/lmp_fix_deform.mod",
        fit_file="dir.depmod/box_evolution_fit.csv",
        poly_file="dir.depmod/poly_coeffs.csv",
    )
