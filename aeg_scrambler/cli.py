"""Copyright (c) Ronnie Crawford.

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
def rank() -> None:
    """Prioritises genes based on weights of various factors

    Args:
        name: Name to display in the greeting.

    Returns:
        None
    """


@app.command()
def design() -> None:
    """Taking a list of ranked genes, will design pegRNA sequences for the 
    regions between enhancers and call PRIDICT to find the most efficient
    in each inter-enhancer sequence

    Args:
        name: Name to display in the greeting.
        formal: Whether to use formal or informal greeting.

    Returns:
        None
    """
        
@app.command()
def explore() -> None:
    """Taking a list of genes will display their distribution and 
    generate several sets of graphs for several features

    Args:
        name: Name to display in the greeting.

    Returns:
        None
    """

@app.command()
def exportX() -> None:
    """Export either a wig or bed file of plateau regions,
    convolution signal etc as decided by options

    Args:
        name: Name to display in the greeting.

    Returns:
        None
    """


def main() -> None:
    """Application's main entrypoint.

    Returns:
        None
    """
    app()


if __name__ == "__main__":
    main()
