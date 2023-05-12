"""Copyright (c) Ronnie Crawford.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import typer

import config
import gene_annotations
import gene_expression
import regulatory_element_annotations
import metrics
import coordinates
import sequences


app = typer.Typer()

@app.command()
def configure(configuration, path):
    # Look for and read config file, then update config
    
    configuration.set_config_from_file(configuration, path)

@app.command()
def rank(configuration):
    
    # Prioritises genes based on weights of various factors
    
    instance_gene_annotations = gene_annotations.GeneAnnotations(configuration)
    instance_gene_expressions = gene_expression.GeneExpression(configuration)
    instance_regulatory_element_annotations = \
        regulatory_element_annotations \
            .RegulatoryElementAnnotations(configuration)

    instance_metrics = metrics.Metrics(configuration,
                                       instance_gene_annotations,
                                       instance_regulatory_element_annotations,
                                       instance_gene_expressions)
    instance_metrics.print_metrics()

@app.command()
def tune(gene, metrics):
    
    # Tune the weights to find the local maxima rank of a specific gene
    
    

@app.command()
def design(configuration, metrics):
    
    # Taking a list of ranked genes, will design pegRNA sequences for the 
    # regions between enhancers and call PRIDICT to find the most efficient
    # in each inter-enhancer sequence
    
    coordinates = coordinates.Coordinates(configuration, metrics)
    sequences = sequences.Sequences(configuration, coordinates)
    
        
@app.command()
def explore():
    
    # Taking a list of genes will display their distribution and 
    # generate several sets of graphs for several features
    
    print("Exploring")

@app.command()
def export():
    
    # Export either a wig or bed file of plateau regions,
    # convolution signal etc as decided by options

    print("Exporting")

def main():
    
    app()
    configuration = config.Config()

if __name__ == "__main__":
    main()
