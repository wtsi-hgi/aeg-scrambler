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
    def __init__(self, filename: Path) -> None:
        self.filename = filename
    
    def load(self) -> pd.DataFrame:
        data = pd.read_csv(self.filename)
        return self.clean(data)

    @abstractmethod
    def clean(self, data: pd.DataFrame) -> pd.DataFrame:
        """Every derived class must implement this"""
        pass

class GeneralData(GeneData):
    def clean(self, data: pd.DataFrame) -> pd.DataFrame:
        print("cleaning GeneralData")
    

class SpecificData(GeneData):
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
    def clean(self, data: pd.DataFrame) -> pd.DataFrame:
        print("cleaning AnnotationData")

    
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

a = data_loader.load("a.csv", GeneDataLoader.Type.general)
b = data_loader.load("b.csv", GeneDataLoader.Type.specific)
c = data_loader.load("c.csv", GeneDataLoader.Type.annotation)    