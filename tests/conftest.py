import os

import pytest


# additional fixtures defined in each test case scope
pytest_plugins = [
    "tests.test_data.io_read_annotated_bed_data",
]


@pytest.fixture
def simple_test_file(tmp_path):
    target_output = os.path.join(tmp_path, "simple_test_file.txt")

    with open(target_output, "w+") as fh:
        fh.write("foo\nbar\nbaz\n")

    return target_output
