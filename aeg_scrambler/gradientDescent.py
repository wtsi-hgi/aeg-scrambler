import torch
import numpy as np
import pandas as pd

desired_score_col = 'Desired_score'
interest_score_col = 'Interest_score'
gene_name_col = 'Gene_name'

class GradientDescent:
    
    def __init__(
            self, 
            df,
            columns,
            learning_rate=1e-4, 
            epochs=1e5
        ):
        """Finds the optimal weights for calculating the interest score,
            given a list of genes of interest within the dataframe.

        Args:
            df - a dataframe of genes.

            columns - a list of columns to be included as part of the 
            optimisation.
            
            learning_rate=1e-4 - by how much the model nudges the weights
            in the direction of the loss.

            epochs=1e5 - how many iterations the model will run.
        
        Returns:
            weights - a list of optimal weights.
        """

        assert len(columns) != 0
        assert df is not None

        self.df = df
        self.columns = columns
        self.learning_rate = learning_rate
        self.epochs = epochs
        self.weights = torch.randn(5, requires_grad=True)

    def assign_gene_priority(self, genes, agnostic=True) -> pd.DataFrame:
        """Assigns a score distribution to a dataframe containing genes. 
            The distribution is the 'desired' interest scores of the genes,
            and will thus be used to calculate the 'loss' of the gradient 
            descent model.

        Args:
            genes - list of genes the user is interested in,
            the model will assign more priority to these and aim to 
            optimise weights such that they end up with a relatively high 
            interest score.

            agnostic = True - if set to False, the genes will 
            not receieve equal priority amongst themselves. It will be 
            assumed instead that although they are all important, the 
            first is more important than the second, the second than the
            third, and so on.
            
        Returns:
            self.df - the objects dataframe, which contains the new column.
        """
        assert len(genes) < 11

        genes = set(genes)
        
        # define column
        self.df[desired_score_col] = 0 
        self.df.loc[
            self.df[gene_name_col].isin(genes), desired_score_col
        ] = 10

        return self.df

    def calculate_loss(self) -> float:
        return sum((self.df[interest_score_col] - self.df[desired_score_col])**2)

    def optimise_weights(self) -> pd.DataFrame:
        """Find optimal weights to satisfy a desired distribution of genes 
            within a dataframe.

        Args:
            None.

        Returns:
            self.df - Dataframe with rearranged genes.
        """

        print('Optimising weights...')

        # model output
        def forward(x):
            return torch.sum(self.weights * x, axis=1)

        # loss = MSE
        def loss(y, y_pred):
            return torch.mean((y_pred - y) ** 2)

        # Extract relevant columns
        X = torch.tensor(self.df[self.columns].values, requires_grad=True)
        Y = torch.tensor(self.df[desired_score_col].values)

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

        self.update_interest_score()
        return self.df

    def update_interest_score(self) -> None:
        """Updates the the interest score for each row of the model, using
            its current weights.

        Returns:
            None.
        """

        def interest_score(row):
            score = sum([self.weights[val] * row[col] 
                        for val, col in enumerate(self.columns)])
            
            return score.item()

        self.df[interest_score_col] = self.df.apply(
            interest_score,
            axis=1
        )

        self.df = self.df.sort_values(
            [interest_score_col], 
            ascending = False
        )

        self.df = self.df.drop(columns=[desired_score_col])
    
    def __str__(self):
        return f"""Gradient Descent model with
        weights={self.weights},
        lr={self.learning_rate},
        epochs={self.epochs}"""
    