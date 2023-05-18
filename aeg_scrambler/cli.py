"""Copyright (c) Ronnie Crawford.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import typer

from .config import Config
from .gene_annotations import GeneAnnotations
from .gene_expression import GeneExpression
from .regulatory_element_annotations import RegulatoryElementAnnotations
from .metrics import Metrics
from .coordinates import Coordinates
from .sequences import Sequences

app = typer.Typer()

@app.command()
def configure(configuration, path):
    # Look for and read config file, then update config
    
    configuration.set_config_from_file(configuration, path)

@app.command()
def rank(configuration):
    
    """
    Prioritises genes based on weights of various factors
    """
    
    instance_gene_annotations = GeneAnnotations(configuration)
    instance_gene_expressions = GeneExpression(configuration)
    instance_regulatory_element_annotations = \
        RegulatoryElementAnnotations(configuration)

    instance_metrics = Metrics(
        configuration,
        instance_gene_annotations,
        instance_regulatory_element_annotations,
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
