import sys
from pathlib import Path

TEST_DATA_DIR = Path(__file__).absolute().parent.joinpath("test_data")

sys.path.append(Path(__file__).absolute().parent.parent)
