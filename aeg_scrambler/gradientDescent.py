import torch
import numpy as np
import pandas as pd

class GradientDescent:
    
    def __init__(
            self, 
            df,
            weight_names,
            learning_rate=1e-4, 
            epochs=1e5
        ):
        """
            Finds the optimal weights for calculating the interest score,
            given a list of genes of interest within the dataframe.

            Args:
                df - a dataframe of genes.

                weight_names - a list of weight names.

                genes_of_interest - a list of gene names of interest.
                
                learning_rate=1e-4 - by how much the model nudges the weights
                in the direction of the loss.

                epochs=1e5 - how many iterations the model will run.
            
            Returns:
                weights - a list of optimal weights.
        """

        self.columns = [
            'Gene_name',
            'Interest_score',
            'Desired_score',
            'Scaled_std',
            'Scaled_anomalous_score',
            'Scaled_enhancer_count',
            'Scaled_enhancer_proportion',
            'Scaled_specific_gene_expression',
            'Scaled_gene_size'
        ]

        self.df = df
        self.weight_names = weight_names
        self.learning_rate = learning_rate
        self.epochs = epochs
        self.weights = torch.randn(5, requires_grad=True)

    def assign_gene_priority(self, genes_of_interest, agnostic=True):
        """
            Assigns a score distribution to a dataframe containing genes. 
            The distribution is the 'desired' interest scores of the genes,
            and will thus be used to calculate the 'loss' of the gradient 
            descent model.

            Args:
                genes_of_interest - list of genes the user is interested in,
                the model will assign more priority to these and aim to 
                optimise weights such that they end up with a relatively high 
                interest score.

                agnostic = True - if set to False, the genes_of_interest will 
                not receieve equal priority amongst themselves. It will be 
                assumed instead that although all 10 genes are important, the
                first is more important than the second, the second than the 
                third, and so on.
            
            Returns:
                score_distribution - a distribution of scores which represents
                the 'desired' gene interest scores
        """
        assert len(genes_of_interest) < 11

        genes_of_interest_set = set(genes_of_interest)
        self.df['Desired_score'] = 0 
        self.df.loc[self.df['Gene_name'].isin(genes_of_interest_set), 'Desired_score'] = 10

    def calculate_loss(self) -> float:
        return sum((self.df['Interest_score'] - self.df['Desired_score'])**2)

    def optimise_weights(self) -> pd.DataFrame:
        """
            Find optimal weights to satisfy a desired distribution of genes 
            within a dataframe.

            Args:
                None.

            Returns:
                Dataframe with rearranged genes.
        """

        # model output
        def forward(x):
            return torch.sum(self.weights * x, axis=1)

        # loss = MSE
        def loss(y, y_pred):
            return torch.mean((y_pred - y) ** 2)

        # Extract relevant columns
        X = torch.tensor(df[
            ['Scaled_std', 
             'Scaled_anomalous_score', 
             'Scaled_enhancer_count',
             'Scaled_enhancer_proportion', 
             'Scaled_specific_gene_expression']].values,
             requires_grad=True
        )
        Y = torch.tensor(df['Desired_score'].values)

        # Training
        for epoch in range(int(self.epochs)):
            # predict = forward pass
            y_pred = forward(X)

            # loss
            l = loss(Y, y_pred)

            # calculate gradients = backward pass
            l.backward()

            # update weights
            with torch.no_grad():
                self.weights -= self.learning_rate * self.weights.grad

            # zero the gradients after updating
            self.weights.grad.zero_()

            if epoch % 10_000 == 0:
                print(f'{epoch}...')
                print(f'Loss: {l.data}')

        print(f'New weights: {self.weights}')

        self.update_interest_score()
        self.loss = self.calculate_loss()

        return self.df

    def update_interest_score(self) -> None:
        """
            Calculates the the interest score for each row in the model, 
            using the current weights of the model.
        """

        def interest_score(row):
            score = sum([self.weights[val] * row[weight_col] 
                        for val, weight_col in enumerate(self.weight_names)])
            
            return score.item()

        self.df['Interest_score'] = self.df.apply(
            interest_score,
            axis=1
        )
        self.df = self.df.sort_values(
            ['Interest_score'], 
            ascending = False
        )
        self.df = self.df.loc[:, self.columns]
    
    def __str__(self):
        return f"""Gradient Descent model with
        weights={self.weights},
        lr={self.learning_rate},
        epochs={self.epochs}"""

# Read csv file
df = pd.read_csv(
    '/',
    skiprows=56,
    sep='\t')

# Drop first column
df = df.drop(df.columns[0], axis=1)

weight_names = [
    'Scaled_std', 
    'Scaled_anomalous_score', 
    'Scaled_enhancer_count', 
    'Scaled_enhancer_proportion', 
    'Scaled_specific_gene_expression'
]

model = GradientDescent(df, weight_names)
model.assign_gene_priority(['TPM1', 'SOX2'])
print(model.df)

model.optimise_weights()
print(model.df)