import pytest

from depmod.io import read_atom, guess_filetype, guess_atom_format

from .conftest import template_dir


@pytest.fixture(scope="module")
def atomreader_sample_dir(template_dir):
    return template_dir / "atom_samples"


class Test_AtomReader_GuessFileType:
    def test_filetype(self, atomreader_sample_dir):
        testfile = atomreader_sample_dir / "Cu_fcc_2x3x4.lmp"
        assert guess_filetype(testfile) == (False, False)

    def test_filetype_tar(self, atomreader_sample_dir):
        testfile = atomreader_sample_dir / "Cu_fcc_2x3x4.lmp.tar"
        filetype = guess_filetype(testfile)
        assert guess_filetype(testfile) == (False, True)

    def test_filetype_gz(self, atomreader_sample_dir):
        testfile = atomreader_sample_dir / "Cu_fcc_2x3x4.lmp.gz"
        assert guess_filetype(testfile) == (True, False)

    def test_filetype_tar_gz(self, atomreader_sample_dir):
        testfile = atomreader_sample_dir / "Cu_fcc_2x3x4.lmp.tar.gz"
        assert guess_filetype(testfile) == (True, True)


class Test_AtomReader_GuessAtomFormat:
    def test_filetype_lmp(self, atomreader_sample_dir):
        testfile = atomreader_sample_dir / "Cu_fcc_2x3x4.lmp"
        assert guess_atom_format(testfile) == "lmp"
    
    def test_filetype_lmp_tar(self, atomreader_sample_dir):
        testfile = atomreader_sample_dir / "Cu_fcc_2x3x4.lmp.tar"
        assert guess_atom_format(testfile) == "lmp"
    
    def test_filetype_lmp_gz(self, atomreader_sample_dir):
        testfile = atomreader_sample_dir / "Cu_fcc_2x3x4.lmp.gz"
        assert guess_atom_format(testfile) == "lmp"
    
    def test_filetype_lmp_tar_gz(self, atomreader_sample_dir):
        testfile = atomreader_sample_dir / "Cu_fcc_2x3x4.lmp.tar.gz"
        assert guess_atom_format(testfile) == "lmp"


class Test_AtomReader_ReadAtom_LammpsData:

    def assert_orthorombic_file(self, testfile):
        a = 3.58
        infos = read_atom(testfile)

        assert infos.filepath == testfile
        assert infos.natom == 96
        
        # check diagonal terms
        for inds, r in zip([(0, 0), (1, 1), (2, 2)], [2, 3, 4]):
            assert (infos.cell[inds] - (r * a)) <= 1e-8

        # check all diagonal terms that should be zero
        for inds in [(1, 0), (0, 1), (2, 0), (0, 2), (2, 1), (1, 2)]:
            assert (infos.cell[inds] - 0.0) <= 1e-8

    def assert_triclinic_file(self, testfile):
        a = 3.58
        e = 0.01

        infos = read_atom(testfile)
        
        assert infos.filepath == testfile
        assert infos.natom == 96
        
        # check diagonal terms
        for ij, r in zip([(0, 0), (1, 1), (2, 2)], [2, 3, 4]):
            assert (infos.cell[ij] - (r * a)) <= 1e-8

        # check upper diag term
        for ij, v in zip([(0, 1), (0, 2), (1, 2)], [0.0, 0.0, 4*a*e]):
            assert (infos.cell[ij] - v) <= 1e-8

        # check lower diag term (sould be zeros)
        for ij in [(1, 0), (2, 0), (2, 1)]:
            assert (infos.cell[ij] - 0.0) <= 1e-8

    def test_ortho_read(self, atomreader_sample_dir):
        self.assert_orthorombic_file(atomreader_sample_dir / "Cu_fcc_2x3x4.lmp")
    
    def test_ortho_read_tar(self, atomreader_sample_dir):
        self.assert_orthorombic_file(atomreader_sample_dir / "Cu_fcc_2x3x4.lmp.tar")
    
    def test_ortho_read_gz(self, atomreader_sample_dir):
        self.assert_orthorombic_file(atomreader_sample_dir / "Cu_fcc_2x3x4.lmp.gz")
    
    def test_ortho_read_tar_gz(self, atomreader_sample_dir):
        self.assert_orthorombic_file(atomreader_sample_dir / "Cu_fcc_2x3x4.lmp.tar.gz")

    def test_triclinic_read(self, atomreader_sample_dir):
        self.assert_triclinic_file(atomreader_sample_dir / "Cu_fcc_2x3x4_triclinic.lmp")

    def test_triclinic_read_tar(self, atomreader_sample_dir):
        self.assert_triclinic_file(atomreader_sample_dir / "Cu_fcc_2x3x4_triclinic.lmp.tar")
    
    def test_triclinic_read_gz(self, atomreader_sample_dir):
        self.assert_triclinic_file(atomreader_sample_dir / "Cu_fcc_2x3x4_triclinic.lmp.gz")
    
    def test_triclinic_read_tar_gz(self, atomreader_sample_dir):
        self.assert_triclinic_file(atomreader_sample_dir / "Cu_fcc_2x3x4_triclinic.lmp.tar.gz")
