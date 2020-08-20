import math
import os
from unittest.mock import patch

import pytest

from varifier import vcf_stats

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "vcf_stats")


class TestFrsFromVcfRecord:
    @patch("cyvcf2.Variant", create=True, autospec=True)
    def test_frsPresentAndFloat(self, mocked_record):
        frs = 0.42
        frs_key = "FRS"
        mocked_record.FORMAT = [frs_key]
        mocked_record.format.return_value = [[frs]]

        actual = vcf_stats._frs_from_vcf_record(mocked_record)
        mocked_record.format.assert_called_once_with(frs_key)
        expected = frs

        assert actual == expected

    @patch("cyvcf2.Variant", create=True, autospec=True)
    def test_frsPresentAndIsNull(self, mocked_record):
        frs = math.nan
        frs_key = "FRS"
        mocked_record.FORMAT = [frs_key]
        mocked_record.format.return_value = [[frs]]

        actual = vcf_stats._frs_from_vcf_record(mocked_record)
        mocked_record.format.assert_called_once_with(frs_key)
        expected = "NA"

        assert actual == expected

    @patch("cyvcf2.Variant", create=True, autospec=True)
    def test_frsNotPresentAndCovNotPresent(self, mocked_record):
        mocked_record.FORMAT = []

        actual = vcf_stats._frs_from_vcf_record(mocked_record)
        expected = "NA"

        assert actual == expected

    @patch("cyvcf2.Variant", create=True, autospec=True)
    def test_frsNotPresentButCovPresent(self, mocked_record):
        cov_key = "COV"
        cov = [1, 3]
        mocked_record.FORMAT = [cov_key]
        mocked_record.genotypes = [[0, 0, False]]
        mocked_record.format.return_value = [cov]

        actual = vcf_stats._frs_from_vcf_record(mocked_record)
        mocked_record.format.assert_called_once_with(cov_key)
        expected = 0.25

        assert actual == expected

    @patch("cyvcf2.Variant", create=True, autospec=True)
    def test_frsNotPresentCovPresentButGenotypeIsNull(self, mocked_record):
        cov_key = "COV"
        cov = [1, 3]
        mocked_record.FORMAT = [cov_key]
        mocked_record.genotypes = [[-1, -1]]
        mocked_record.format.return_value = [cov]

        actual = vcf_stats._frs_from_vcf_record(mocked_record)
        expected = "NA"

        assert actual == expected

    @patch("cyvcf2.Variant", create=True, autospec=True)
    def test_frsNotPresentCovPresentButGenotypeIsHet(self, mocked_record):
        cov_key = "COV"
        cov = [1, 3]
        mocked_record.FORMAT = [cov_key]
        mocked_record.genotypes = [[0, 1]]
        mocked_record.format.return_value = [cov]

        actual = vcf_stats._frs_from_vcf_record(mocked_record)
        expected = "NA"

        assert actual == expected

    @patch("cyvcf2.Variant", create=True, autospec=True)
    def test_frsNotPresentCovPresentButCovAllZero(self, mocked_record):
        cov_key = "COV"
        cov = [0, 0]
        mocked_record.FORMAT = [cov_key]
        mocked_record.genotypes = [[0, 0]]
        mocked_record.format.return_value = [cov]

        actual = vcf_stats._frs_from_vcf_record(mocked_record)
        expected = 0

        assert actual == expected


def test_format_dict_to_edit_dist_scores():
    format_dict = {
        "VFR_RESULT": "CANNOT_USE_GT",
    }
    assert (None, None) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "TP",
        "VFR_ED_RA": 1,
        "VFR_ED_TR": 1,
        "VFR_ED_TA": 0,
    }
    assert (1, 1) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "TP",
        "VFR_ED_RA": 1,
        "VFR_ED_TR": 0,
        "VFR_ED_TA": 0,
    }
    assert (1, 1) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "TP",
        "VFR_ED_RA": 2,
        "VFR_ED_TR": 6,
        "VFR_ED_TA": 5,
    }
    assert (1, 6) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "FP",
        "VFR_ED_RA": 2,
        "VFR_ED_TR": 0,
        "VFR_ED_TA": 0,
    }
    assert (0, 2) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "FP",
        "VFR_ED_RA": 2,
        "VFR_ED_TR": 6,
        "VFR_ED_TA": 5,
    }
    assert (1, 6) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "TP",
        "VFR_ED_RA": 2,
        "VFR_ED_TR": 0,
        "VFR_ED_TA": 1,
    }
    assert (0, 2) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)

    format_dict = {
        "VFR_RESULT": "TP",
        "VFR_ED_RA": 1,
        "VFR_ED_TR": "NA",
        "VFR_ED_TA": 5,
        "VFR_ALLELE_MATCH_FRAC": 0.75,
    }
    assert (0.75, 1) == vcf_stats.format_dict_to_edit_dist_scores(format_dict)


def test_per_record_stats_from_vcf_file():
    infile = os.path.join(data_dir, "per_record_stats_from_vcf_file.vcf")
    expect = [
        {
            "CHROM": "ref1",
            "DP": 10,
            "DPF": 1.1,
            "FRS": 0.6,
            "GT_CONF": 100.0,
            "GT_CONF_PERCENTILE": 0.12,
            "POS": 1,
            "VFR_IN_MASK": 1,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 0,
            "VFR_ED_TA": 1,
            "VFR_ALLELE_LEN": 1,
            "VFR_ALLELE_MATCH_COUNT": 0,
            "VFR_ALLELE_MATCH_FRAC": 0.0,
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "FP",
        },
        {
            "CHROM": "ref2",
            "DP": 70,
            "DPF": 0.97,
            "FRS": 0.99,
            "GT_CONF": 200.0,
            "GT_CONF_PERCENTILE": 0.95,
            "POS": 2,
            "VFR_IN_MASK": 0,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 0,
            "VFR_ED_TA": 1,
            "VFR_ALLELE_LEN": 1,
            "VFR_ALLELE_MATCH_COUNT": 0,
            "VFR_ALLELE_MATCH_FRAC": 0.0,
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "FP",
        },
    ]
    got = vcf_stats.per_record_stats_from_vcf_file(infile)
    for got_dict, expect_dict in zip(got, expect):
        for k, v in got_dict.items():
            if isinstance(v, str):
                assert v == expect_dict[k]
            else:
                assert pytest.approx(v) == expect_dict[k]


def test_summary_stats_from_per_record_stats():
    record_stats = [
        {
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "TP",
            "VFR_ALLELE_MATCH_FRAC": 0.1,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 1,
            "VFR_ED_TA": 0,
        },
        {
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "TP",
            "VFR_ALLELE_MATCH_FRAC": 0.12,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 1,
            "VFR_ED_TA": 0,
        },
        {
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "FP",
            "VFR_ALLELE_MATCH_FRAC": 0.2,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 0,
            "VFR_ED_TA": 1,
        },
        {
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "Partial_TP",
            "VFR_ALLELE_MATCH_FRAC": 0.25,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 2,
            "VFR_ED_TA": 1,
        },
        {
            "VFR_FILTER": "PASS",
            "VFR_RESULT": "FP",
            "VFR_ALLELE_MATCH_FRAC": 0.3,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 0,
            "VFR_ED_TA": 1,
        },
        {
            "VFR_FILTER": "FAIL",
            "VFR_RESULT": "FP",
            "VFR_ALLELE_MATCH_FRAC": 0.4,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 3,
            "VFR_ED_TA": 1,
        },
        {
            "VFR_FILTER": "FAIL_BUT_TEST",
            "VFR_RESULT": "TP",
            "VFR_ALLELE_MATCH_FRAC": 0.5,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 1,
            "VFR_ED_TA": 0,
        },
        {
            "VFR_FILTER": "PASS",
            "VFR_IN_MASK": 1,
            "VFR_RESULT": "TP",
            "VFR_ALLELE_MATCH_FRAC": 0.1,
            "VFR_ED_RA": 1,
            "VFR_ED_TR": 1,
            "VFR_ED_TA": 0,
        },
    ]
    expect = {
        "UNUSED": {"CONFLICT": 0, "OTHER": 1, "MASKED": 1},
        "ALL": {
            "TP": {"Count": 3, "SUM_ALLELE_MATCH_FRAC": 0.72, "SUM_EDIT_DIST": 3},
            "FP": {"Count": 3, "SUM_ALLELE_MATCH_FRAC": 0.75, "SUM_EDIT_DIST": 3},
            "EDIT_DIST_COUNTS": {"numerator": 4, "denominator": 7},
        },
        "FILT": {
            "TP": {"Count": 2, "SUM_ALLELE_MATCH_FRAC": 0.22, "SUM_EDIT_DIST": 2},
            "FP": {"Count": 3, "SUM_ALLELE_MATCH_FRAC": 0.75, "SUM_EDIT_DIST": 3},
            "EDIT_DIST_COUNTS": {"numerator": 3, "denominator": 6},
        },
    }
    got = vcf_stats.summary_stats_from_per_record_stats(record_stats)
    assert got == expect

    for all_or_filt in "ALL", "FILT":
        expect[all_or_filt]["FN"] = expect[all_or_filt]["FP"]
        del expect[all_or_filt]["FP"]
    got = vcf_stats.summary_stats_from_per_record_stats(record_stats, for_recall=True)
    assert got == expect
