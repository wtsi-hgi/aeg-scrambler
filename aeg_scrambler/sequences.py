import pandas as pd
import pyranges as pr
import re
import subprocess

class Sequences:
    
    def __init__(self, config, coordinates):
        
        self.export_plateaus_for_each_gene(config)

    def export_plateaus_for_each_gene(self, config):
        
        for index, gene in self.data.head(
            config.convolution_limit
        ).iterrows():

            plateaus = pd.DataFrame(
                {"Start" : self.data.loc[index, "Plateau_starts"], 
                "End" : self.data.loc[index, "Plateau_ends"]
                }
            )
            plateaus["Gene_name"] = gene["Gene_name"]
            plateaus["Chromosome"] = "chr" + gene["Chromosome"]
            plateaus["Strand"] = gene["Strand"]
        
            print(plateaus)