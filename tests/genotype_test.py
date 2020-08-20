import pytest

from varifier.genotype import Genotype, GenotypeError


class TestGenotypeIsNull:
    def test_isNull_returnsTrue(self):
        gt = Genotype(-1, -1)

        assert gt.is_null()

    def test_isHom_returnsFalse(self):
        gt = Genotype(-1, 1)

        assert not gt.is_null()


class TestGenotypeIsHom:
    def test_isHom_returnsTrue(self):
        gt = Genotype(0, 0)

        assert gt.is_hom()

    def test_isHet_returnsFalse(self):
        gt = Genotype(1, 0)

        assert not gt.is_hom()

    def test_isNull_returnsFalse(self):
        gt = Genotype(-1, -1)

        assert not gt.is_hom()

    def test_oneAlleleIsNull_returnsTrue(self):
        gt = Genotype(-1, 1)

        assert gt.is_hom()


class TestGenotypeIsHet:
    def test_isHet_returnsTrue(self):
        gt = Genotype(2, 1)

        assert gt.is_het()

    def test_isHom_returnsFalse(self):
        gt = Genotype(0, 0)

        assert not gt.is_het()

    def test_oneAlleleIsNull_returnsFalse(self):
        gt = Genotype(-1, 1)

        assert not gt.is_het()

    def test_isNull_returnsFalse(self):
        gt = Genotype(-1, -1)

        assert not gt.is_het()


class TestGenotypeIsHomRef:
    def test_isHomRef_returnsTrue(self):
        gt = Genotype(0, 0)

        assert gt.is_hom_ref()

    def test_isHomRefOneNull_returnsTrue(self):
        gt = Genotype(-1, 0)

        assert gt.is_hom_ref()

    def test_isHomAlt_returnsFalse(self):
        gt = Genotype(2, 2)

        assert not gt.is_hom_ref()


class TestGenotypeIsHomAlt:
    def test_isHomAlt_returnsTrue(self):
        gt = Genotype(1, 1)

        assert gt.is_hom_alt()

    def test_isHomAltOneNull_returnsTrue(self):
        gt = Genotype(-1, 3)

        assert gt.is_hom_alt()

    def test_isHomRef_returnsFalse(self):
        gt = Genotype(0, 0)

        assert not gt.is_hom_alt()


class TestGenotypeAltIndex:
    def test_isHomAlt_returnsIndex(self):
        gt = Genotype(1, 1)

        actual = gt.alt_index()
        expected = 0

        assert actual == expected

    def test_isHomAltWithOneNull_returnsIndex(self):
        gt = Genotype(-1, 5)

        actual = gt.alt_index()
        expected = 4

        assert actual == expected

    def test_isHomRef_returnsNone(self):
        gt = Genotype(0, 0)

        assert gt.alt_index() is None

    def test_isHet_returnsNone(self):
        gt = Genotype(0, 3)

        assert gt.alt_index() is None


class TestGenotypeAlleleIndex:
    def test_isHomAlt_returnsIndex(self):
        gt = Genotype(1, 1)

        actual = gt.allele_index()
        expected = 1

        assert actual == expected

    def test_isHomAltWithOneNull_returnsIndex(self):
        gt = Genotype(-1, 5)

        actual = gt.allele_index()
        expected = 5

        assert actual == expected

    def test_isHomRef_returnIndex(self):
        gt = Genotype(0, 0)

        actual = gt.allele_index()
        expected = 0

        assert actual == expected

    def test_isHet_returnsNone(self):
        gt = Genotype(0, 3)

        assert gt.allele_index() is None

    def test_bothNull_returnsNone(self):
        gt = Genotype(-1, -1)

        assert gt.allele_index() is None


class TestGenotypeFromArray:
    def test_emptyArr_returnsNull(self):
        arr = []

        actual = Genotype.from_arr(arr)
        expected = Genotype()

        assert actual == expected

    def test_arrWithTooManyAlleles_raisesError(self):
        arr = [2, 4, 5]

        with pytest.raises(GenotypeError):
            Genotype.from_arr(arr)

    def test_arrWithNonInt_nonIntFilteredOut(self):
        arr = [2, 4, False]

        actual = Genotype.from_arr(arr)
        expected = Genotype(2, 4)

        assert actual == expected

    def test_arrWithNoInts_returnsNull(self):
        arr = [False, "foo"]

        actual = Genotype.from_arr(arr)
        expected = Genotype()

        assert actual == expected

    def test_arrWithOneAllele_returnsNullInSecondAllele(self):
        arr = [0]

        actual = Genotype.from_arr(arr)
        expected = Genotype(0, -1)

        assert actual == expected

    def test_arrIsSetInstead_returnsGenotype(self):
        arr = {0, 1}

        actual = Genotype.from_arr(arr)
        expected = Genotype(0, 1)

        assert actual == expected

    def test_arrIsString_raisesError(self):
        arr = "0/3"

        with pytest.raises(GenotypeError):
            Genotype.from_arr(arr)
