from typing import NamedTuple, Optional, Iterable


class GenotypeError(Exception):
    pass


class Genotype(NamedTuple):
    """Class to handle genotypes from cyvcf2"""

    allele1: int = -1
    allele2: int = -1

    def is_null(self) -> bool:
        """Is the genotype null. i.e. ./."""
        return self.allele1 == -1 and self.allele2 == -1

    def is_hom(self) -> bool:
        """Is the genotype homozygous"""
        if self.is_null():
            return False
        if self.allele1 == -1 or self.allele2 == -1:
            return True
        return self.allele1 == self.allele2

    def is_het(self) -> bool:
        """Is the genotype heterozyhous"""
        return not self.is_null() and not self.is_hom()

    def is_hom_ref(self) -> bool:
        """Is genotype homozygous reference?"""
        return self.is_hom() and (self.allele1 == 0 or self.allele2 == 0)

    def is_hom_alt(self) -> bool:
        """Is genotype homozygous alternate?"""
        return self.is_hom() and (self.allele1 > 0 or self.allele2 > 0)

    def alt_index(self) -> Optional[int]:
        """If the genotype is homozygous alternate, returns the 0-based index of the
        alt allele in the alternate allele array.
        """
        if not self.is_hom_alt():
            return None
        return self.allele_index() - 1

    def allele_index(self) -> Optional[int]:
        """If the genotype is homozygous, returns the index of the genotype, otherwise
        it returns None"""
        if not self.is_hom():
            return None
        return max(self.allele1, self.allele2)

    @staticmethod
    def from_arr(arr: Iterable) -> "Genotype":
        """Create a Genotype from an array of genotypes from cyvcf2"""
        if type(arr) is str:
            raise GenotypeError("String not supported in Genotype.from_arr")

        alleles = [a for a in arr if type(a) is int]
        if len(alleles) > 2:
            raise GenotypeError(
                f"Got too many alleles. Only two alleles are currently supported, "
                f"got: {alleles}"
            )
        return Genotype(*alleles)
