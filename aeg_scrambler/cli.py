"""
Copyright (c) Ronnie Crawford.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import typer
from typing import Optional

from config import Config
from gene_annotations import GeneAnnotations
from gene_expression import GeneExpression
from regulatory_annotations import RegulatoryAnnotations
from metrics import Metrics
from coordinates import Coordinates
from sequences import Sequences

app = typer.Typer()

@app.command()
def hello(name: Optional[str] = None):
    if name:
        typer.echo(f"Hello {name}")
    else:
        typer.echo("Hello World!")

@app.command()
def configure_settings(path: Optional[str] = None):
    
    """
    Looks for and reads config file, then updates user settings.
    """
    configuration = Config()

    if path:
        configuration.set_config_from_file(path)
    else:
        configuration.set_config_from_file("../../config.yaml")
    view_settings()

@app.command()
def view_settings():
    
    """
    Shows current user settings.
    """

    configuration = Config()
    typer.echo("Current settings:")
    configuration.print_config()

@app.command()
def rank(configuration):
    
    """
    Prioritises genes based on weights of various factors.
    """
    
    instance_gene_annotations = GeneAnnotations(configuration)
    instance_gene_expressions = GeneExpression(configuration)
    instance_regulatory_annotations = \
        RegulatoryAnnotations(configuration)

    instance_metrics = Metrics(
        configuration,
        instance_gene_annotations,
        instance_regulatory_annotations,
        instance_gene_expressions
    )
    instance_metrics.print_metrics()

@app.command()
def tune(gene, metrics):
    """
    Tune the weights to find the local maxima rank of a specific gene
    """
    pass

@app.command()
def design(configuration, metrics):
    
    """
    Taking a list of ranked genes, will design pegRNA sequences for the 
    regions between enhancers and call PRIDICT to find the most efficient
    in each inter-enhancer sequence
    """
    
    coordinates = Coordinates(configuration, metrics)
    sequences = Sequences(configuration, coordinates)

@app.command()
def explore():
    
    """
    Taking a list of genes will display their distribution and 
    generate several sets of graphs for several features
    """
    
    print("Exploring")
    pass

@app.command()
def export():
    
    """
    Export either a wig or bed file of plateau regions,
    convolution signal etc as decided by options
    """

    print("Exporting")

def main():
    
    app()
    configuration = Config()
    

if __name__ == "__main__":
    main()
