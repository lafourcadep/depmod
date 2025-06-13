import pytest

from depmod import read_lattice

from .conftest import template_dir


@pytest.fixture(scope="module")
def atomreader_sample_dir(template_dir):
    return template_dir / "atom_samples"


class Test_AtomReader_ReadAtom_LammpsData:
    
    fmt = ("lmp", "xyz")

    def assert_orthorombic_file(self, testfile):
        a = 3.58
        lattice = read_lattice(testfile)

        for inds, r in zip([(0, 0), (1, 1), (2, 2)], [2, 3, 4]):
            assert (lattice[inds] - (r * a)) <= 1e-8

        for inds in [(1, 0), (0, 1), (2, 0), (0, 2), (2, 1), (1, 2)]:
            assert (lattice[inds] - 0.0) <= 1e-8

    def assert_triclinic_file(self, testfile):
        a = 3.58
        e = 0.01

        lattice = read_lattice(testfile)

        for ij, r in zip([(0, 0), (1, 1), (2, 2)], [2, 3, 4]):
            assert (lattice[ij] - (r * a)) <= 1e-8

        for ij, v in zip([(0, 1), (0, 2), (1, 2)], [0.0, 0.0, 4 * a * e]):
            assert (lattice[ij] - v) <= 1e-8

        for ij in [(1, 0), (2, 0), (2, 1)]:
            assert (lattice[ij] - 0.0) <= 1e-8

    def test_read(self, atomreader_sample_dir):
        for fmt in self.fmt:
            self.assert_orthorombic_file(atomreader_sample_dir / f"Cu_fcc_2x3x4.{fmt}")
            self.assert_triclinic_file(atomreader_sample_dir / f"Cu_fcc_2x3x4_triclinic.{fmt}")
    
    def test_ortho_read_gzip(self, atomreader_sample_dir):
        for fmt in self.fmt:
            self.assert_orthorombic_file(atomreader_sample_dir / f"Cu_fcc_2x3x4.{fmt}.gz")
            self.assert_triclinic_file(atomreader_sample_dir / f"Cu_fcc_2x3x4_triclinic.{fmt}.gz")
    
    def test_ortho_read_bzip2(self, atomreader_sample_dir):
        for fmt in self.fmt:
            self.assert_orthorombic_file(atomreader_sample_dir / f"Cu_fcc_2x3x4.{fmt}.bz2")
            self.assert_triclinic_file(atomreader_sample_dir / f"Cu_fcc_2x3x4_triclinic.{fmt}.bz2")

    def test_ortho_read_lzma(self, atomreader_sample_dir):
        for fmt in self.fmt:
            self.assert_orthorombic_file(atomreader_sample_dir / f"Cu_fcc_2x3x4.{fmt}.xz")
            self.assert_triclinic_file(atomreader_sample_dir / f"Cu_fcc_2x3x4_triclinic.{fmt}.xz")
