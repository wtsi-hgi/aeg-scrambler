"""Copyright (c) Crawford Mace.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import pandas as pd
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
        names, skiprows = None, None
        sep = ',' if self.datatype.is_general() else '\t'

        if self.datatype.is_annotation():
            names = ['Chromosome','Source','Type','Start','End','Score', 
                    'Strand','Phase','Attributes']
            skiprows = 5
        
        elif self.datatype.is_specific():
            names = ['Gene_name','Specific_gene_expression']
            skiprows = 1
        
        dtype = {name: str for name in names} if names else None
        
        data = pd.read_csv(
            self.reference,
            sep=sep,
            names=names, 
            skiprows=skiprows, 
            dtype=dtype
            )
        data = data.transpose() if self.datatype.is_general() else data
        
        return data.head()

    def __str__(self) -> str:
        return f"Genedata called {self.name}!"
    
    def __repr__(self) -> str:
        return {self.name, self.datatype, self.reference}
    
    def mean(self):
        pass
    
specificPath = "specifc/path/here"
generalPath = "general/path/here"
annotationPath = "annotation/path/here"

hap1Cells = GeneData('hap1', generalPath, Datatype('general'))
print(hap1Cells.data)