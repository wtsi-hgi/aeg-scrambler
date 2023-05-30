"""Copyright (c) Crawford Mace.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import enum
import re
import pandas as pd
import numpy as np 
from pathlib import Path
from abc import abstractmethod, ABC

class GeneData(ABC):
    skiprows = 0
    columns = None
    separator = "\t"
    
    def __init__(self, filename: Path) -> None:
        """Initialises genetic data object.

        Args: 
            filename - path which points to the location of the data.
        
        Returns:
            None.
        """

        self.filename = filename
        self.interesting_chromosomes = [str(c) for c in range(23)] + ['X','Y']
        self.data = self.clean(self.load())

    def load(self) -> pd.DataFrame:
        print("Loading data...")

        return pd.read_csv(
            self.filename, 
            names=self.columns, 
            skiprows=self.skiprows, 
            sep=self.separator
        )
    
    def __repr__(self) -> str:
        return f"""
        Dataframe of type {self.__class__.__name__}, 
        {self.data}"""
    
    @abstractmethod
    def clean(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Every derived class must implement this
        """
        pass

class GeneralData(GeneData):
    separator = ","

    def clean(self, data: pd.DataFrame) -> pd.DataFrame:
        return data
    

class SpecificData(GeneData):
    skiprows = 1
    columns = ['Gene_name', 'Specific_gene_expression']

    def clean(self, data: pd.DataFrame) -> pd.DataFrame:
        """Convert -inf to 0; convert log2 to normal base."""

        data = data.drop_duplicates(keep=False, subset='Gene_name')
        data['Specific_gene_expression'] = np.where(
            data['Specific_gene_expression'] == '-Inf',
            0,
            2 ** data['Specific_gene_expression']
        )
        return data

        
class AnnotationData(GeneData):
    skip_rows = 5
    columns = [
        'Chromosome', 
        'Source', 
        'Type', 
        'Start', 
        'End', 
        'Score',
        'Strand', 
        'Phase', 
        'Attributes'
    ]

    def clean(self, data: pd.DataFrame) -> pd.DataFrame:
        geneNameRegex = "gene_name \"(.*?)\""
        geneBiotypeRegex = "gene_biotype \"(.*?)\""

        data = data.query('Chromosome in @self.interesting_chromosomes')
        data = data[data["Type"] == "gene"]
        data["Gene_biotype"] = data["Attributes"].apply(
            lambda x : re.findall(geneBiotypeRegex, x)[0]
            if re.search(geneNameRegex, x) else "None"
        )
        data = data[data["Gene_biotype"] == "protein_coding"]
        
        data["Gene_name"] = data["Attributes"].apply(
            lambda x : re.findall(geneNameRegex, x)[0] 
            if re.search(geneNameRegex, x) else "None"
        )
        data = data.drop(["Source", "Type", "Score", "Phase", "Attributes", "Gene_biotype"], axis = 1)
        data = data.drop_duplicates(keep = False, subset = ["Gene_name"])
        data = data.rename(columns = {"Start" : "Gene_start", "End" : "Gene_end"})
    
        return data