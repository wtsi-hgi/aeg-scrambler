"""
Copyright (c) Ronnie Crawford.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import typer
from typing import Optional
import pickle

from .gradientDescent import GradientDescent
from .config import Config
from .input_data import (
    CCLEExpression,
    ExperimentalExpression,
    GeneAnnotations,
    RegulatoryAnnotations
)
from .metrics import Metrics
from .coordinates import Coordinates
from .sequences import Sequences

app = typer.Typer()
working_directory = "working/"

def load_data_from_config(config):
    """Loads the data found at the locations specified within the config,
    merges these into a dataframe, and converts this into a Metrics object.

    Args:
        config - path which points towards the config file location.

    Returns:
        Metrics - a Metrics object which contains all the merged data found 
        at the paths specified by the config.
    """
    config = Config(config)
    
    print('Finding CCLE Expression...')
    ccle_expression = CCLEExpression(config)

    print('Finding Experimental Expression...')
    experimental_expression = ExperimentalExpression(config)

    print('Finding Gene Annotations...')
    gene_annotations = GeneAnnotations(config)

    print('Finding Regulatory Annotations...')
    regulatory_annotations = RegulatoryAnnotations(config)
    
    print('Merging Data...')
    metrics = Metrics(
        config,
        gene_annotations,
        regulatory_annotations,
        ccle_expression,
        experimental_expression
    )

    return metrics, config

@app.command()
def prioritise(config:str, genes:list[str], agnostic=True):
    """Given a set of genes of interest, will attempt to reorganise the 
    dataframe so that these genes will be prioritised and will appear 
    more towards the top of the dataframe. This means that other genes not 
    explicitly mentioned, with similar features to those mentioned, will
    naturally drift towards the top of the dataframe as well

    Args:
        config - path which points towards the config file location.
        
        genes - genes of interest.

        agnostic - if False, the 0th gene of interest will receive more 
        priority than the 1st, the 1st more so than the 2nd, and so on. 
        If True, all genes of interest are assigned equal priority.

    Returns:
        DataFrame - a dataframe with a priority towards those genes of interest.
    """

    metrics, config = load_data_from_config(config)
    df = metrics.data

    model = GradientDescent(df, genes, config)
    model.assign_gene_priority(genes)
    model.optimise_weights()

    print(model.df.head(15))
    return

@app.command()
def rank(config = None):    
    """
    Ranks the genes.
    """
    
    config = Config(config)
    
    ccle_expression = CCLEExpression(config)
    experimental_expression = ExperimentalExpression(config)
    gene_annotations = GeneAnnotations(config)
    regulatory_annotations = RegulatoryAnnotations(config)
    
    metrics = Metrics(
        config,
        gene_annotations,
        regulatory_annotations,
        ccle_expression,
        experimental_expression
    )
    print(metrics.printable_ranks())
    metrics.export_gene_scores_report(config)

@app.command()
def design(config = None):
    
    """
    
    """
    
    config = Config(config)
    
    ccle_expression = CCLEExpression(config)
    experimental_expression = ExperimentalExpression(config)
    gene_annotations = GeneAnnotations(config)
    regulatory_annotations = RegulatoryAnnotations(config)
    
    metrics = Metrics(
        config,
        gene_annotations,
        regulatory_annotations,
        ccle_expression,
        experimental_expression
    )
    
    coordinates = Coordinates(config, metrics)
    sequences = Sequences(config, coordinates)

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
def load_user_data(config_id = None):
    
    """
    Loads in user supplied data as specified in the config file.
    """
    
    config = unpickle_object(config_id)
    
    ccle_expression = CCLEExpression(config)
    experimental_expression = ExperimentalExpression(config)
    gene_annotations = GeneAnnotations(config)
    regulatory_annotations = RegulatoryAnnotations(config)
    
    pickle_object(ccle_expression)
    pickle_object(experimental_expression)
    pickle_object(gene_annotations)
    pickle_object(regulatory_annotations)

@app.command()
def generate_metrics(
    config_id,
    ccle_expression_id,
    experimental_expression_id,
    gene_annotations_id,
    regulatory_annotations_id
):
    
    """
    Using entered data, finds various metrics for each
    gene and ranks them according to weights in table.
    """
    config = unpickle_object(config_id)
    ccle_expression = unpickle_object(ccle_expression_id)
    experimental_expression = unpickle_object(experimental_expression_id)
    gene_annotations = unpickle_object(gene_annotations_id)
    regulatory_annotations = unpickle_object(regulatory_annotations_id)
    
    metrics = Metrics(
        config,
        gene_annotations,
        regulatory_annotations,
        ccle_expression,
        experimental_expression
    )
   
@app.command()
def tune(gene, metrics):
    """
    Tune the weights to find the local maxima rank of a specific gene
    """
    pass

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
