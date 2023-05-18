"""Copyright (c) Crawford Mace.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import enum
import pandas as pd
import numpy as np 
from pathlib import Path
from abc import abstractmethod, ABC

class GeneData(ABC):
    skiprows = 0
    columns = None
    separator = ","
    
    def __init__(self, filename: Path) -> None:
        print(f"{self.__class__.__name__}")
        self.filename = filename
        self.data = self.clean(self.load())

    def load(self) -> pd.DataFrame:
        print(f"loading data: sep={self.separator} columns={self.columns}, skiprows={self.skiprows}") 
        return pd.read_csv(self.filename, names=self.columns, skiprows=self.skiprows, sep=self.separator)
    
    @abstractmethod
    def clean(self, data: pd.DataFrame) -> pd.DataFrame:
        """Every derived class must implement this"""
        pass

class GeneralData(GeneData):
    separator = "\t"

    def clean(self, data: pd.DataFrame) -> pd.DataFrame:
        print("cleaning GeneralData")
        pass
    

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
        print("cleaning AnnotationData")
        pass

    
class GeneDataLoader:
    class Type(enum.Enum):
        general = enum.auto()
        specific = enum.auto()
        annotation = enum.auto()
        
    loaders = {
        Type.general: GeneralData,
        Type.specific: SpecificData,
        Type.annotation: AnnotationData
    }

    def load(self, filename: Path, data_type: Type) -> pd.DataFrame:
        loader = self.loaders[data_type](filename)
        return loader.load()
    
    
data_loader = GeneDataLoader()

#a = data_loader.load("a.csv", GeneDataLoader.Type.general)
#b = data_loader.load("b.csv", GeneDataLoader.Type.specific)
#c = data_loader.load("c.csv", GeneDataLoader.Type.annotation)    

specificPath = "/specific/path/here"
generalPath = "/general/path/here"
annotationPath = "/annotation/path/here"

a = GeneralData(generalPath)
b = SpecificData(specificPath)
c = AnnotationData(annotationPath)
print(a)