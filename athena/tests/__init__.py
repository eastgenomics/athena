import sys
from pathlib import Path

TEST_DATA_DIR = Path(__file__).absolute().parent.joinpath("test_data")
# sys.path.append(TEST_DATA_DIR)

# add in parent athena dir to path for importing
# sys.path.append(
#     os.path.abspath(os.path.join(os.path.realpath(__file__), "../../"))
# )

sys.path.append(Path(__file__).absolute().parent.parent)
