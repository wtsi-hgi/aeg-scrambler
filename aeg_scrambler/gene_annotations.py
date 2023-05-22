import pandas as pd
import re

class GeneAnnotations:
    
    def __init__(self, config):
        
        self.read_gene_annotations(config)
        self.clean_gene_annotations(config)
            
    def read_gene_annotations(self, config):
        
        """
        Reads gene annotations gtf file into pandas dataframe
        """
        
        try:
            
            self.data = \
                pd.read_csv(config.gene_annotation_reference,
                            sep = "\t",
                            names = ["Chromosome",
                                     "Source",
                                     "Type",
                                    "Start",
                                    "End",
                                    "Score",
                                    "Strand",
                                    "Phase",
                                    "Attributes"],
                            skiprows = 5,
                            dtype = {"Chromosome" : str,
                                     "Source" : str,
                                    "Type" : str,
                                    "Start" : int,
                                    "End" : int,
                                    "Score" : str,
                                    "Strand" : str,
                                    "Phase" : str,
                                    "Attributes" : str}
                                        )
        
        except:
            
            print("ERROR: Gene annotations file could not be read.")
            
    def clean_gene_annotations(self, config):
        
        """
        Puts annotation data into correct format and removes unecessary data
        """
        
        self.data = self.data.loc[self.data["Chromosome"]
                                  .isin(config.chromosomes_of_interest)]
        self.data = self.data \
            .drop(self.data[self.data["Type"] != "gene"].index)
        self.data["Gene_biotype"] = self.data["Attributes"]\
            .apply(lambda x : re.findall("gene_biotype \"(.*?)\"", x)[0] if \
                re.search("gene_name \"(.*?)\"", x) != None else "None")
        self.data = self.data \
            .drop(self.data[self.data["Gene_biotype"] != \
                "protein_coding"].index)
        self.data["Gene_name"] = self.data["Attributes"]\
            .apply(lambda x : re.findall("gene_name \"(.*?)\"", x)[0] if\
                re.search("gene_name \"(.*?)\"", x) != None else "None")
        self.data = self.data.drop(["Source", "Type", "Score",
                                                "Phase", "Attributes",
                                                "Gene_biotype"], axis = 1)
        self.data = self.data \
            .drop_duplicates(keep = False, subset = ["Gene_name"])
        self.data = self.data \
            .rename(columns = {"Start" : "Gene_start", "End" : "Gene_end"})