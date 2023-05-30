import torch
import pandas as pd

class GradientDescent:
    
    def __init__(
            self, 
            df,
            weight_names,
            #genes_of_interest, 
            *column_weights, 
            learning_rate=1e-4, 
            epochs=1e4
        ):
        """
            Finds the optimal weights for calculating the interest score,
            given a list of genes of interest within the dataframe.

            Args:
                df - a dataframe of genes.

                weight_names - a list of weight names.

                genes_of_interest - a list of gene names of interest.

                *column_weights - 
                
                learning_rate=1e-4 - by how much the model nudges the weights
                in the direction of the loss.

                epochs=1e4 - how many iterations the model will run.
            
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

        self.df = df[self.columns]
        self.weight_names = weight_names
        self.learning_rate = learning_rate
        self.epochs = epochs
        self.weights = torch.randn(5, requires_grad=True)
        self.loss = self.calculate_loss()


    def calculate_loss(self) -> float:
        return (self.df['Interest_score'] - self.df['Desired_score'])**2

    def optimise_weights(self) -> pd.DataFrame:
        """
        """

        # model output
        def forward(x):
            return torch.sum(self.weights * x, axis=1)

        # loss = MSE
        def loss(y, y_pred):
            return torch.mean((y_pred - y) ** 2)  # Compute mean squared error instead of element-wise squared differences

        # Extract relevant columns as NumPy array
        X = torch.tensor(df[['Scaled_std', 'Scaled_anomalous_score', 'Scaled_enhancer_count',
                            'Scaled_enhancer_proportion', 'Scaled_specific_gene_expression']].values,
                        requires_grad=True)
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

            # Quite arbitrary loss decay that slightly improves the model
            if epoch > 5000:
                self.learning_rate *= 0.1

            if epoch % 1000 == 0:
                print(f'{epoch}...')

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
            calc = (
                self.weights[0] * row['Scaled_std'] + 
                self.weights[1] * row['Scaled_anomalous_score'] + 
                self.weights[2] * row['Scaled_enhancer_count'] + 
                self.weights[3] * row['Scaled_enhancer_proportion'] +
                self.weights[4] * row['Scaled_specific_gene_expression']
            )
            
            return calc.item()

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

# Get first 10 rows
df = df[:10]

# Add desired score column
df['Desired_score'] = [0, 0, 9, 1, 0, 8, 0, 0, 7, 0]

weight_names = [
    'Scaled_std', 
    'Scaled_anomalous_score', 
    'Scaled_enhancer_count', 
    'Scaled_enhancer_proportion', 
    'Scaled_specific_gene_expression'
]

model = GradientDescent(df, weight_names)
model.optimise_weights()
weights = model.weights
print(weights)
print(model.df)

loss = sum(df['Desired_score'] - df['Interest_score'])**2
print(f'Loss beforehand: {loss}')


loss = sum(df['Desired_score'] - df['New_Score'])**2
df = df.sort_values(['New_Score'], ascending = False)
print(f'Overall loss: {loss}')
print()

print(df)