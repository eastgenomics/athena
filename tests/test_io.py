import os
from pathlib import Path
from uuid import uuid4
from unittest.mock import patch

import polars as pl
import polars.testing as pl_testing
import pytest

from athena.utils import io
from tests import TEST_DATA_DIR
from tests.test_data import calculate_test_data as test_data


class TestReadFile:
    def test_file_contents_correctly_read_from_file(self, simple_test_file):
        read_contents = io.read_file(simple_test_file)

        assert read_contents == "foo\nbar\nbaz\n"

    def test_file_not_found_error_raised_when_file_does_not_exists(self):
        with pytest.raises(FileNotFoundError):
            io.read_file("not_a_file.txt")


class TestReadImage:
    def test_image_contents_correctly_read_from_file(self):
        read_contents = io.read_image(
            Path(TEST_DATA_DIR).joinpath("single_pixel.png")
        )

        expected_contents = (
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABAQMAAAAl21bKAAAAA1BMVEUAAACnej3a"
            "AAAAAXRSTlMAQObYZgAAAApJREFUCNdjYAAAAAIAAeIhvDMAAAAASUVORK5CYII="
        )

        assert read_contents == expected_contents

    def test_file_not_found_error_raised_when_file_does_not_exists(self):
        with pytest.raises(FileNotFoundError):
            io.read_image("not_a_file.png")


class TestReadAnnotatedBed:
    pass
