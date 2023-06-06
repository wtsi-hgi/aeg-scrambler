"""
Copyright (c) Ronnie Crawford.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import typer

from .config import Config
from .coordinates import Coordinates
from .gradient_descent import GradientDescent
from .input_data import (
    CCLEExpression,
    ExperimentalExpression,
    GeneAnnotations,
    RegulatoryAnnotations,
)
from .metrics import Metrics
from .sequences import Sequences

app = typer.Typer()
working_directory = "working/"


def main():
    app()


@app.command()
def rank(config=None):
    """Ranks the genes based on the weights provided within the config file.

    Args:
        config - path which points towards the config file location.
    """

    metrics, config = load_data_from_config(config)

    print(metrics.printable_ranks())
    metrics.export_gene_scores_report(config)


@app.command()
def tune(config: str, genes: list[str], agnostic: bool):
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
    model.assign_gene_priority(genes, agnostic)
    model.optimise_weights()

    print(model.df.head(20))
    return


@app.command()
def design(config=None):
    """

    Args:
        config - path which points towards the config file location.
    """

    metrics, config = load_data_from_config(config)

    print('Finding coordinates...')
    coordinates = Coordinates(config, metrics)

    print('Finding sequences...')
    Sequences(config, coordinates)


def load_data_from_config(config: str):
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
        experimental_expression,
    )

    return metrics, config


if __name__ == "__main__":
    main()
