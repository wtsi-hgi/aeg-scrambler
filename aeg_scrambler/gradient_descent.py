import numpy as np
import pandas as pd
import torch

desired_score_col = 'Desired_score'
interest_score_col = 'Interest_score'
gene_name_col = 'Gene_name'


class GradientDescent:
    def __init__(self, df, genes, config, learning_rate=1e-4, epochs=1e5):
        """Finds the optimal weights for calculating the interest score,
            given a list of genes of interest within the dataframe.

        Args:
            df - a dataframe of genes.

            genes - a list of genes of interest.

            learning_rate=1e-4 - by how much the model nudges the weights
            in the direction of the loss.

            epochs=1e5 - how many iterations the model will run.

        Returns:
            weights - a list of optimal weights.
        """

        assert len(genes) != 0
        assert df is not None

        self.columns = config.columns

        # get the names of those columns where the dtype is not str
        self.numerical_cols = [
            col_name
            for col_name, col_type in self.columns.items()
            if col_type != str
        ]

        column_names = [i for i, _ in self.columns.items()]
        self.df = df[column_names]

        self.learning_rate = learning_rate
        self.epochs = epochs

        # generate randomised weights
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

        # define column
        self.df[desired_score_col] = 0

        if agnostic:
            self.df.loc[
                self.df[gene_name_col].isin(genes), desired_score_col
            ] = 100
        else:
            distribution = np.linspace(100, 10, len(genes))
            new_val = {gene: score for gene, score in zip(genes, distribution)}
            print(new_val)

            self.df[desired_score_col] = self.df[gene_name_col].map(new_val)

        self.df = self.df.fillna(0)

        return self.df

    def calculate_loss(self) -> float:
        return sum(
            (self.df[interest_score_col] - self.df[desired_score_col]) ** 2
        )

    def optimise_weights(self) -> pd.DataFrame:
        """Find optimal weights to satisfy a desired distribution of genes
            within a dataframe.

        Returns:
            self.df - Dataframe with rearranged genes.
        """

        print('Optimising weights...')

        # model output
        def forward(x):
            return torch.sum(self.weights * x, dim=1)

        # loss = Mean Squared Error (MSE)
        def loss(y, y_pred):
            return torch.mean((y_pred - y) ** 2)

        # Extract relevant columns
        numerical_vals = self.df[self.numerical_cols].astype(float)

        X = torch.tensor(numerical_vals.values, requires_grad=True)
        Y = torch.tensor(self.df[desired_score_col].values)

        # Training
        for _ in range(int(self.epochs)):
            # predict = forward pass
            y_pred = forward(X)

            # calculate gradients = backward pass
            l = loss(Y, y_pred)
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
            score = sum(
                [
                    self.weights[val] * row[col]
                    for val, col in enumerate(self.numerical_cols)
                ]
            )

            return score.item()

        self.df[interest_score_col] = self.df.apply(interest_score, axis=1)

        self.df = self.df.sort_values([interest_score_col], ascending=False)

        self.df = self.df.drop(columns=[desired_score_col])

    def __str__(self):
        return f"""Gradient Descent model with
        weights={self.weights},
        lr={self.learning_rate},
        epochs={self.epochs}"""
