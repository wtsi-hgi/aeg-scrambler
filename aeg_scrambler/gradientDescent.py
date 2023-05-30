import torch
import pandas as pd

# Read csv file
df = pd.read_csv(
    'dataframe/file/here',
    skiprows=56,
    sep='\t')

# Drop first column
df = df.drop(df.columns[0], axis=1)

# Get first 10 rows
df = df[:20]

# Add desired score column
df['Desired_score'] = [8, 10, 11, 4, 3, 6, 3, 8, 7, 10, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4]

def gradientDescent(weights, learning_rate=1e-4, epochs=1e4):
    """
        Performs gradient descent on a dataframe, which must contain 
        an interest score and a desired interest score column.
        
        Args:
            weights - a tensor of the initialised weights.

            learning_rate=1e-4 - by how much the model nudges the 
            weights in the direction of the gradient of the loss.
            
            epochs=1e5 - how many iterations the model will run.
            
        Returns:
            weights - a new set of optimal weights.
    """
    
    # model output
    def forward(x):
        return sum(weights * x)

    # loss = MSE
    def loss(y, y_pred):
        return (y_pred - y)**2

    print('Weights before: ', weights.data)

    # Training
    for epoch in range(int(epochs)):
        for _, row in df.iterrows():
            X = torch.tensor(
                [row['Scaled_std'], 
                row['Scaled_anomalous_score'], 
                row['Scaled_enhancer_count'], 
                row['Scaled_enhancer_proportion'],
                row['Scaled_specific_gene_expression']
                ], 
                requires_grad=True)
            
            Y = torch.tensor(row['Desired_score'])
            
            # predict = forward pass
            y_pred = forward(X)

            # loss
            l = loss(Y, y_pred)

            # calculate gradients = backward pass
            l.backward()

            # update weights
            with torch.no_grad():
                weights -= learning_rate * weights.grad

            # zero the gradients after updating
            weights.grad.zero_()
        
        # Quite arbitrary loss decay that slightly improves the model
        if epoch > 5000:
            learning_rate *= 0.1
            
        if epoch % 1000 == 0:
            print(f'{epoch}...')

    print()
    print('New weights: ', weights.data)
    return weights

loss = sum(df['Desired_score'] - df['Interest_score'])**2
print(f'Loss beforehand: {loss}')

weights = torch.tensor([1, 0.5, 0.25, 0.75, 0.8], dtype=torch.float32, requires_grad=True)
weights = gradientDescent(weights)

def calc_interest_score(row):
    calc = (
        row['Scaled_std']*weights[0] + 
        row['Scaled_anomalous_score']*weights[1] + 
        row['Scaled_enhancer_count']*weights[2] + 
        row['Scaled_enhancer_proportion']*weights[3] +
        row['Scaled_specific_gene_expression']*weights[4]
    )
    
    return calc.item()

df['New_Score'] = df.apply(calc_interest_score, axis=1)
df = df[['Gene_name', 
    'Interest_score',
    'Desired_score',
    'New_Score',
    'Scaled_std',
    'Scaled_anomalous_score',
    'Scaled_enhancer_count',
    'Scaled_enhancer_proportion',
    'Scaled_specific_gene_expression',
    'Scaled_gene_size']
]

loss = sum(df['Desired_score'] - df['New_Score'])**2
print(f'Overall loss: {loss}')
print()

print(df)