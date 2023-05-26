# Read csv file
df = pd.read_csv(
    'datafram/path/here',
    skiprows=56,
    sep='\t')

# Drop first column
df = df.drop(df.columns[0], axis=1)

# Get first 10 rows
df = df[:10]

# Add desired score column
df['Desired_score'] = [8, 10, 11, 4, 3, 6, 3, 8, 7, 10]

# Create new df
df = df[
    ['Gene_name', 
     'Interest_score',
     'Desired_score',
     'Scaled_std',
     'Scaled_anomalous_score',
     'Scaled_enhancer_count',
     'Scaled_enhancer_proportion',
     'Scaled_specific_gene_expression',
     'Scaled_gene_size'
    ]]

weights = {
    'Scaled_std': 1,
    'Scaled_anomalous_score': 0.5,
    'Scaled_enhancer_count': 0.25,
    'Scaled_enhancer_proportion': 0.75,
    'Scaled_specific_gene_expression': 0.8
}

total_loss = []

def gradientDescent(learning_rate=0.001, epochs=1000, **weights):
    """
        Performs gradient descent on a dataframe, which must contain 
        an interest score and a desired interest score column.
        
        Args:
            learning_rate=0.001 - by how much the model nudges the 
            weights in the direction of the gradient of the loss.
            
            epochs=1000 - how many iterations the model will run.
            
            **weights - a dictionary of the form {weightname : weight}.
        
        Returns:
            weights - a new set of optimal weights.
    """
    
    weight_names, weights = zip(*weights['weights'].items())
    weight_names, weights = list(weight_names), list(weights)
    
    for e in range(epochs):                
        for _, row in df.iterrows():
            
            # Calculate interest score with current weights
            interest_score_calc = sum([weights[i] * row[weight_names[i]] for i in range(len(weights))])  
            
            # Calculate loss
            loss = (row['Desired_score'] - interest_score_calc)

            for i in range(len(weights)):
                # Calculate partial derivative
                partial_derivative = -2 * weights[i] * loss

                # Update weight
                weights[i] -= (learning_rate * partial_derivative)
            
        loss = sum((df['Desired_score'] - interest_score_calc) ** 2)
        total_loss.append(loss)
        
        if e % 100 == 0: 
            print(e)
            loss = sum((df['Desired_score'] - interest_score_calc) ** 2)
            print('loss: ', loss)
            print('Weights: ', weights)
            print()

    return weights
        
gradientDescent(weights=weights)