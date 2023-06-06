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
        """
        
        self.assign_unique_id()
        self.load(config)
        self.clean(config)
        
    def __str__(self) -> str:
        return f"""
        Dataframe of type {self.__class__.__name__}, 
        {self.data}"""
    
    def assign_unique_id(self):
        """
        Generates a unique ID for each dataframe.
        """
        
        self.unique_id = self.__class__.__name__ + str(hash(self))
    
    @abstractmethod
    def load(self, config: Config) -> pd.DataFrame:
        """
        Every derived class must implement this
        """

        pass
    
    @abstractmethod
    def clean(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Every derived class must implement this
        """
        
        pass
    
    @abstractmethod
    def export(self, config):
        """Every derived class must implement this
        """
    

class CCLEExpression(InputData):
    
    def load(self, config: Config) -> pd.DataFrame:
        """Loads the CCLEExpression data into a dataframe.

        Args:
            config - user's configuration file.
        
        Returns:
            pd.DataFrame - loaded dataframe.
        """

        self.data = pd.read_csv(config.ccle_expression_path).transpose()

    def clean(self, config: Config) -> pd.DataFrame:
        """Cleans general expression data. 
        
        Args:
            config - user's configuration file.
        
        Returns:
            pd.DataFrame - cleaned general expression dataframe.
        """

        self.data.drop(
            ["depmapID", "primary_disease"],
            axis = 0,
            inplace = True
        )
        self.data.columns = self.data.iloc[0]
        self.data = self.data.iloc[1:]
        self.data = self.data.reset_index().rename(
            columns = {"index" : "Gene_name"}
        )
        self.data = self.data.rename(
            columns = {config.cell_line_of_interest : "General_gene_expression"}
        )
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
        
    def find_mean(self) -> None:
        """Adds an element-wise mean to the given expression dataframe.
        """

        self.data["Mean"] = self.data.loc[
            :, self.data.columns != "Gene_name"
        ].mean(axis = 1)
        
    def find_std(self) -> None:
        """Adds the standard deviation to the given expression dataframe, 
        giving the standard deviation across expression of each gene in all
        provided cell types.
        """

        self.data["Std"] = self.data.loc[
            :, self.data.columns != "Gene_name"
        ].std(axis = 1)
    
    def find_anomalous_score_of_gene_expression(self) -> None:
        """Adds the z-score of each gene based on its expression
        in the cell line of interest compared to all others.
        """
        
        self.data["Anomalous_score"] = self.data.apply(
            lambda gene : (gene["General_gene_expression"] - 
                           gene["Mean"]) / gene["Std"],
            axis = 1
        )

class ExperimentalExpression(InputData):
    
    def load(self, config) -> None:
        self.data = pd.read_csv(
            config.experimental_expression_path,
            sep = "\t",
            names = ["Gene_name", "Specific_gene_expression"],
            skiprows = 1
        )
       
    def clean(self, _) -> None:
        """
        Cleans the expression data specific to the cell line of interest by
        turning minus infinite strings into a floating point representation
        of negative infinity, duplicate genes are dropped.
        """
        
        self.data["Specific_gene_expression"] = np.where(
            self.data["Specific_gene_expression"] == "-Inf", 
            0, 
            2 ** self.data["Specific_gene_expression"]
        )

        self.data.drop_duplicates(subset = ["Gene_name"], inplace = True)
        
class GeneAnnotations(InputData):
    
    def load(self, config) -> None:
        """
        Load GeneAnnotation data into a dataframe.
        """
        
        column_names = {
            "Chromosome": str,
            "Source": str,
            "Type": str,
            "Start": int,
            "End": int,
            "Score": str,
            "Strand": str,
            "Phase": str,
            "Attributes": str
        }
        
        self.data = pd.read_csv(
            config.gene_annotation_path,
            sep = "\t",
            names = column_names.keys(),
            skiprows = 5,
            dtype = column_names
        )
    
    def clean(self, config) -> None:
        """Puts annotation data into correct format and removes unecessary 
        data.

        Args:
            config - user's configuration file.
        """
    
        self.data = self.data[
            self.data["Chromosome"].isin(config.chromosomes_of_interest)
        ]

        self.data = self.data.query('Type == "gene"')

        # parse Attributes column for content after 'gene_biotype'.
        self.data["Gene_biotype"] = self.data["Attributes"].str.extract(
            'gene_biotype "(.*?)"', expand=False).fillna("None")

        self.data = self.data[self.data["Gene_biotype"] == "protein_coding"]

        # parse Attributes column for content after 'gene_name'.
        self.data["Gene_name"] = self.data["Attributes"].str.extract(
            'gene_name "(.*?)"', expand=False).fillna("None")
        
        unwanted_cols = [
            "Source", "Type", "Score", "Phase", "Attributes", "Gene_biotype"
        ]
        self.data.drop(
            unwanted_cols,
            axis = 1,
            inplace = True
        )

        self.data.drop_duplicates(subset=["Gene_name"], inplace = True)

        self.data.rename(
            columns={"Start": "Gene_start", "End": "Gene_end"}, 
            inplace = True
        )
        
    def export(self, config):
        
        self.data.to_csv(
            config.results_directory + "gene_annotations.bed",
            sep = "\t",
            columns = ["Chromosome", "Gene_start", "Gene_end", "Gene_name"]
        )
        
        return config.results_directory + "gene_annotations.bed"

class RegulatoryAnnotations(InputData):
    
    def load(self, config) -> None:
        """Loads RegulatoryAnnotation data into a dataframe.

        Args:
            config - user's configuration file.
        """
        
        self.data = pd.read_csv(
            config.regulatory_elements_path,
            sep = "\t",
            names = ["Chromosome", "Start", "End", "Flag"]
        )
    
    def clean(self, config) -> None:
        """Puts regulatory data into correct format and remove unecessary 
        data.

        Args:
            config - user's configuration file.
        """
            
        self.data["Chromosome"] = self.data["Chromosome"].str[3:]

        self.data = self.data[
            self.data["Flag"].isin(config.flags_of_interest)
        ]

        self.data.drop("Flag", axis = 1)
        
    def export(self, config):
        
        self.data.to_csv(
            config.results_directory + "regulatory_annotations.bed",
            sep = "\t",
            columns = ["Chromosome", "Start", "End"]
        )
        
        return config.results_directory + "regulatory_annotations.bed"