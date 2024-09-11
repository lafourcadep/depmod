import pytest
import pathlib


@pytest.fixture(scope="session")
def template_dir():
    return pathlib.Path(__file__).parent.resolve() / "dir.template"
