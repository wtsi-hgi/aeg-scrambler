"""
Copyright (c) Ronnie Crawford.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import typer
from typing import Optional
import pickle

from .config import Config
from .input_data import InputData
from .gene_expression import GeneExpression
from .regulatory_annotations import RegulatoryAnnotations
from .metrics import Metrics
from .coordinates import Coordinates
from .sequences import Sequences

app = typer.Typer()
working_directory = "working/"
tracker = []

@app.command()
def view_tracker():
    
    pass

@app.command()
def set_config(path = None):
    
    """
    Initialises config for program,
    can read from file otherwise will be set to defaults.
    """
    
    config = Config(path)
    print("Created new config: " + config.unique_id)
    pickle_object(config)
   
@app.command() 
def view_config(config_id = None):
    
    """
    Used to view settings of given config.
    """
    
    config = unpickle_object(config_id)
    print(config)

@app.command()
def test(path = None):
    
    # config --set
    config = Config("config.json")
    
    # rank
    
    #gene_expressions = GeneExpression(config)
    #regulatory_annotations = RegulatoryAnnotations(config)
    #metrics = Metrics(config,
    #                  gene_annotations,
    #                  regulatory_annotations,
    #                  gene_expressions)
    
    # rank --tune
    
    # rank --view
    
    # rank --export
    #metrics.export_gene_scores_report(config)
    
    # cluster
    #coordinates = Coordinates(config, metrics)
    
    # cluster --export
    #coordinates.export_convolutions(config)
    
    # design
    #sequences = Sequences(config, coordinates)
    
    #design --export

@app.command()
def rank(configuration):
    
    """
    Prioritises genes based on weights of various factors.
    """
    
    instance_gene_annotations = GeneAnnotations(configuration)
    instance_gene_expressions = GeneExpression(configuration)
    instance_regulatory_annotations = RegulatoryAnnotations(configuration)

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

    pass
    
def pickle_object(object):
    
    with open(working_directory + object.unique_id, "wb") as file:
        
        pickle.dump(object, file)
        
    tracker = tracker.append(object.unique_id)
    
def unpickle_object(object_id):
    
    with open(working_directory + object_id, "rb") as file:
    
        return pickle.load(file)

def main():
    
    app()
    
if __name__ == "__main__":
    main()
