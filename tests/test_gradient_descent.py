import pytest
import pandas as pd

from aeg_scrambler.config import Config
from aeg_scrambler.gradient_descent import GradientDescent

def test_object_initialisation():
    data = {
        'Gene_name': ['GAPDH', 'PTMA', 'RPL13A', 'SMAD3', 'TLN2', 'THSD4'],
        'Scaled_std': [0.702, 0.786, 0.695, 1.511, 1.211, 1.657],
        'Scaled_anomalous_score': [-0.976, 0.274, 1.598, -0.154, 0.484, 0.032],
        'Scaled_enhancer_count': [1, 25, 2, 94, 130, 122],
        'Scaled_enhancer_proportion': [0.114, 0.099, 0.118, 0.168, 0.073, 0.083],
        'Scaled_specific_gene_expression': [2318.383, 1791.166, 1718.482, 58.241, 56.301, 6.605],
    }

    df = pd.DataFrame(data)
    genes = ['TLN2']
    config = Config()

    model = GradientDescent(df, genes, config)
    
    assert model.df.equals(df)
