"""Copyright (c) Crawford Mace.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

from typing import Any

import pytest
from click.testing import Result
from typer.testing import CliRunner

from aeg_scrambler.cli import app


class CLI:
    """CLI fixcture wrapper."""

    def __init__(self) -> None:
        """Constructor."""
        self.runner = CliRunner()

    def invoke(self, *args: Any, **kwargs: Any) -> Result:
        """Invoke a CLI command.

        Args:
            *args: Positional arguments to pass to the command
            **kwargs: Keyword arguments to pass to the command

        Returns:
            Result: Exit code from the command.
        """
        return self.runner.invoke(app, *args, **kwargs)


@pytest.fixture
def cli() -> CLI:
    """Return a CLI object as a fixture.

    Returns:
        CLI: A newly instantiated CLI object.
    """
    return CLI()
