import os
import pytest
from pathlib import Path

TEST_DIR = Path(__file__).parent


def pytest_configure(config):
    """Change to the test directory so all relative test_files/ paths resolve."""
    os.chdir(TEST_DIR)


@pytest.fixture(scope="session")
def structures_dir():
    return TEST_DIR / "test_files" / "structures"
