"""Copyright (c) Crawford Mace.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import typer

app = typer.Typer()


@app.command()
def hello(name: str) -> None:
    """Say hello to someone.

    Args:
        name: Name to display in the greeting.

    Returns:
        None
    """
    typer.echo(f"Hello {name}")


@app.command()
def goodbye(name: str, formal: bool = False) -> None:
    """Say goodbye.

    Args:
        name: Name to display in the greeting.
        formal: Whether to use formal or informal greeting.

    Returns:
        None
    """
    if formal:
        typer.echo(f"Goodbye {name}. Have a good day.")
    else:
        typer.echo(f"Bye {name}!")


def main() -> None:
    """Application's main entrypoint.

    Returns:
        None
    """
    app()


if __name__ == "__main__":
    main()
