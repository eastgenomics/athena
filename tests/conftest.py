import os

import pytest


@pytest.fixture
def simple_test_file(tmp_path):
    target_output = os.path.join(tmp_path, "simple_test_file.txt")

    with open(target_output, "w+") as fh:
        fh.write("foo\nbar\nbaz\n")

    return target_output
