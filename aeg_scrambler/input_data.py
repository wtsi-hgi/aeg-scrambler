import re
import pandas as pd
import numpy as np
from abc import abstractmethod, ABC

from .config import Config

class InputData:
    
    def __init__(self, config) -> None:
        
        """Initialises genetic data object.

        Args: 
            filename - path which points to the location of the data.
        
        Returns:
            None.
        """
        
        self.load(config)
        self.clean(config)
        self.assign_unique_id()
        
    def __str__(self) -> str:
        
        return f"""Dataframe of type {self.__class__.__name__}, 
        {self.data}"""
    
    @abstractmethod
    def load(self, config: Config) -> pd.DataFrame:

        return pd.read_csv(
            self.filename, 
            names = self.columns, 
            skiprows = self.skiprows, 
            sep = self.separator
        )
    
    @abstractmethod
    def clean(self, data: pd.DataFrame) -> pd.DataFrame:
        
        """
        Every derived class must implement this
        """
        
        pass
    
    def assign_unique_id(self):
        
        """
        Generates a unique ID for each dataframe.
        """
        
        self.unique_id = self.__class__.__name__ + str(hash(self.data))
    
    def pickle_dataframe(self, config):
        
        """
        Serialises dataframe and saves as file.
        """
        
        self.data.to_pickle(
            config.working_directory +
            self.unique_id
        )
    
    def unpickle_dataframe(self, config):
        
        """
        Unserialises dataframe and loads from file.
        """
        
        pd.read_pickle(
            config.working_directory +
            self.unique_id
        )
    
class CCLEExpression(InputData):
    
    separator = ","

    def clean(self, config: Config) -> pd.DataFrame:
        
        """
        Cleans the general expression data by dropping irrelevant non-numeric
        columns, converting the first row to headers, setting gene name as
        index, renaming columns, finding the means and standard deviation of
        expression for each gene and the z-score of the cell line of
        interest's expression for each gene, dropping irrelevant columns,
        and dropping duplicate genes
        """

        self.data.drop(
            ["depmapID",
                "primary_disease"],
            axis = 0,
            inplace = True
        )
        self.data.columns = self.data.iloc[0]
        self.data = self.data.iloc[1:]
        self.data = self.data.reset_index() \
            .rename(columns = {"index" : "Gene_name"})
        self.data = self.data. \
            rename(columns = 
                {config.cell_line_of_interest : "General_gene_expression"})
        self.find_mean()
        self.find_std()
        self.find_anomalous_score_of_gene_expression()
        self.data = self.data[["Gene_name",
                                            "General_gene_expression",
                                            "Mean",
                                            "Std",
                                            "Anomalous_score"]]
        self.data.drop_duplicates(
            keep = False,
            subset = ["Gene_name"],
            inplace = True
        )

class ExperimentalExpression(InputData):
        
    def clean(self):

        """
        Cleans the expression data specific to the cell line of interest by
        turning minus infinite strings into a floating point representation
        of negative infinity, duplicate genes are dropped
        """
        
        self.data["Specific_gene_expression"] = \
            self.data["Specific_gene_expression"].apply(lambda expression :
                0 if expression == "-Inf" else pow(2, expression))
                
        self.data.drop_duplicates(
            keep = False,
            subset = ["Gene_name"],
            inplace = True
        )
        
class GeneAnnotations(InputData):
    
    def clean(self, config):
    
        """
        Puts annotation data into correct format
        and removes unecessary data
        """
        
        self.data = self.data.loc[self.data["Chromosome"]
                                .isin(config.chromosomes_of_interest)]
        self.data.drop(
            self.data[self.data["Type"] != "gene"].index,
            inplace = True
        )
        self.data["Gene_biotype"] = self.data["Attributes"]\
            .apply(lambda x : re.findall("gene_biotype \"(.*?)\"", x)[0] if \
                re.search("gene_name \"(.*?)\"", x) != None else "None")
        self.data.drop(
            self.data[self.data["Gene_biotype"] != "protein_coding"].index,
            inplace = True
        )
        self.data["Gene_name"] = self.data["Attributes"]\
            .apply(lambda x : re.findall("gene_name \"(.*?)\"", x)[0] if\
                re.search("gene_name \"(.*?)\"", x) != None else "None")
        self.data.drop([
            "Source",
            "Type",
            "Score",
            "Phase",
            "Attributes",
            "Gene_biotype"],
                        axis = 1,
                        inplace = True
                        )
        self.data.drop_duplicates(
            keep = False,
            subset = ["Gene_name"],
            inplace = True
        )
        self.data.rename(
            columns = {"Start" : "Gene_start", "End" : "Gene_end"},
            inplace = True
        )
        
class RegulatoryAnnotations(InputData):
    
    def clean(self, config):
        
        """
        Put regulatory data into correct format and remove unecessary data
        """
            
        self.data["Chromosome"] = self.data["Chromosome"].apply(
            lambda x : x[3:]
        )
        self.data = self.data[
            self.data["Flag"].isin(
                config.enhancer_epigenetic_flags_of_interest
        )]
        self.data.drop(["Flag"], axis = 1, inplace = True)