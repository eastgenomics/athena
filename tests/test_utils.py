import os
import sys
from unittest import TestCase


sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from bin import utils


class TestCleanIndication(TestCase):
    """
    Tests for utils.clean_indication

    Function takes in a string and improves formatting for displaying
    clinical indication in the report summary text through removing
    suffixes and underscore after the R code
    """

    def test_r_code_underscore_suffix_removed(self):
        """
        Test that if the indication string has an R code that the
        following underscore is removed
        """
        self.assertEqual(
            utils.clean_indication("R1.1_some_indication"),
            "R1.1 some_indication"
        )

    def test_suffix_removed(self):
        """
        Test that _G and _P suffixes get correctly removed
        """
        with self.subTest():
            self.assertEqual(
                utils.clean_indication("R1.1_some_indication_G"),
                "R1.1 some_indication"
            )

        with self.subTest():
            self.assertEqual(
                utils.clean_indication("R1.1_some_indication_P"),
                "R1.1 some_indication"
            )

    def test_hgnc_only_handled_correctly(self):
        """
        Test that HGNC IDs have underscore prefixes removed
        """
        self.assertEqual(
            utils.clean_indication("_HGNC:1234;_HGNC:5678"),
            "HGNC:1234; HGNC:5678"
        )

    def tes_indication_and_gene_symbol(self):
        """
        Test that indication and gene symbols together handled correctly
        """
        with self.subTest():
            self.assertEqual(
                utils.clean_indication("R1.1_some_indication_G;_HGNC:1234"),
                "R1.1 some_indication; HGNC:1234"
            )


    def test_indication_not_changed(self):
        """
        Test that other indication string without what is being removed
        remain unchanged
        """
        self.assertEqual(
            utils.clean_indication("some_other_indication"),
            "some_other_indication"
        )
