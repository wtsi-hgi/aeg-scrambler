"""Copyright (c) Crawford Mace.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from enum import Enum

class Datatype(Enum):
    GENERAL = 'general'
    SPECIFIC = 'specific'
    ANNOTATION = 'annotation'

    def is_general(self):
        return self == Datatype.GENERAL
    
    def is_specific(self):
        return self == Datatype.SPECIFIC
    
    def is_annotation(self):
        return self == Datatype.ANNOTATION
    
class GeneData:
    interesting_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']

    def __init__(self, name, reference, datatype=Datatype.GENERAL) -> None:
        """Constructor."""
        assert isinstance(datatype, Datatype), "Invalid datatype"
        assert isinstance(reference, str), "Invalid reference path"

        self.name = name        
        self.reference = Path(reference)
        self.datatype = datatype
        self.data = self.read_data()

    def read_data(self) -> pd.DataFrame:
        """Load data at path reference into a pandas dataframe."""

        names = skiprows = dtype = None
        sep = ',' if self.datatype.is_general() else '\t'

        if self.datatype.is_annotation():
            names = ['Chromosome', 'Source', 'Type', 'Start', 'End', 'Score',
                    'Strand', 'Phase', 'Attributes']
            dtype = {name: str for name in names}
            skiprows = 5
        
        elif self.datatype.is_specific():
            names = ['Gene_name', 'Specific_gene_expression']
            dtype = {'Gene_name': str, 'Specific_gene_expression': float}
            skiprows = 1
        
        
        dtype = {name: str for name in names} if names else None
        
                
        dtype = {name: str for name in names} if names else None
        
        data = pd.read_csv(
            self.reference,
            sep=sep,
            names=names, 
            skiprows=skiprows, 
            dtype=dtype
            )
        data = data.transpose() if self.datatype.is_general() else data

        return data
    
    def clean_data(self, data) -> pd.DataFrame:
        if self.datatype.is_specific():
            return self.clean_specific_data(data)
        
        elif self.datatype.is_general():
            return self.clean_general_data(data)
        
        elif self.datatype.is_annotation():
            return self.clean_annotation_data(data)
        
        raise Exception('Trying to clean invalid datatype')
    
    def clean_specific_data(self, data) -> pd.DataFrame:
        """Convert -inf to 0; convert log2 to normal base."""

        data = data.drop_duplicates(keep=False, subset='Gene_name')

        data['Specific_gene_expression'] = np.where(
            data['Specific_gene_expression'] == '-Inf',
            0,
            2 ** data['Specific_gene_expression']
        )

        return data
    
    def clean_general_data(self, data) -> pd.DataFrame:
        pass
    
    def clean_annotation_data(self, data) -> pd.DataFrame:
        pass

    def __str__(self) -> str:
        return f"""Genedata called {self.name}, head \n{self.data.head()}!"""
    
    def __repr__(self) -> str:
        return {self.name, self.datatype, self.reference}
    
    def mean(self):
        pass
    
specificPath = "specifc/path/here"
generalPath = "general/path/here"
annotationPath = "annotation/path/here"

hap1Cells = GeneData('hap1', annotationPath, Datatype('annotation'))
print(hap1Cells)