from glob import glob
import gzip
import os

import pytest


# additional fixtures defined in each test case scope, hoover all files
# up that might contain fixtures and register them in required format
pytest_plugins = [
    x.replace("/", ".").replace(".py", "")
    for x in glob("tests/test_data/" + "/**/*.py", recursive=True)
]


@pytest.fixture
def simple_test_file(tmp_path):
    target_output = os.path.join(tmp_path, "simple_test_file.txt")

    with open(target_output, "w+") as fh:
        fh.write("foo\nbar\nbaz\n")

    return target_output


@pytest.fixture
def simple_compressed_test_file(tmp_path):
    target_output = os.path.join(
        tmp_path, "simple_compressed_test_file.txt.gz"
    )

    with gzip.open(target_output, "wb") as fh:
        fh.write("foo\nbar\nbaz\n".encode())

    return target_output
