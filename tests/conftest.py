from dataclasses import dataclass
from pathlib import Path
import shutil
import pytest


def pytest_addoption(parser):
    parser.addoption("--all", action="store_true", default=False, help="run slow tests")
    parser.addoption(
        "--update-golden",
        action="store_true",
        help="Update golden test master files to match current output. "
        "ONLY DO THIS IF YOU'VE INTENTIONALLY CHANGED OUTPUTS.",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--all"):
        # --all given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --all option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture
def tests_dir():
    return Path(__file__).parent

@dataclass
class GoldenTest:
    output_file: Path
    updating: bool

    def check_file(self, new_file: Path):
        if self.updating:
            self.output_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(str(new_file), str(self.output_file))
        else:
            expected_output = self.output_file.read_bytes()
            actual_output = new_file.read_bytes()
            assert actual_output == expected_output


    def check_bytes(self, actual_output: bytes):
        if self.updating:
            self.output_file.parent.mkdir(parents=True, exist_ok=True)
            self.output_file.write_bytes(actual_output)
        else:
            expected_output = self.output_file.read_bytes()
            assert actual_output == expected_output

@dataclass
class GoldenTestFixture:
    updating: bool

    def test(self, outfile):
        return GoldenTest(output_file=Path(outfile), updating=self.updating)

@pytest.fixture
def golden(request):
    return GoldenTestFixture(updating=request.config.getoption("update_golden"))
