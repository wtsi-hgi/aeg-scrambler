
import numpy as np
import pandas as pd
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

# Load the gene dataframe
df = pd.read_csv('gene_data.txt', sep='\t')

# drop last row
df = df[:-1]

# drop first column
df = df.drop(df.columns[0], axis=1)
print(df)

# Define the preferred order of genes
preferred_order = ['TMSB4X', 'SOX2', 'ID1', 'HES1', 'TMSB10', 'KRT18', 'DDIT4', 'GADD45G', 
                    'CEBPB', 'PEG10', 'ARL4D', 'SOX4', 'FOS', 'GPR27', 'UCHL1', 'BAMBI', 
                    'RGS2', 'BTG2', 'TPM1', 'MYC', 'PIM1']

# Select the relevant columns for training
feature_columns = ['Std', 'Anomalous_score', 'Specific_gene_expression', 'Enhancer_count', 'Enhancer_proportion', 'Gene_size']
target_column = 'Interest_score'

# Extract feature and target data
X = df[feature_columns].values
y = df[target_column].values

# Scale the feature data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Create a mapping of gene names to indices
gene_indices = {gene: index for index, gene in enumerate(df['Gene_name'])}

# Create pairwise training data
pairwise_data = []
for i in range(len(df)):
    if i // 100:
        print(i)
    gene_i = df.loc[i, 'Gene_name']
    for j in range(i + 1, len(df)):
        gene_j = df.loc[j, 'Gene_name']
        pairwise_data.append((X_scaled[i], X_scaled[j], gene_indices[gene_i] > gene_indices[gene_j]))

print(pairwise_data)

# Split the pairwise data into train and test sets
X_train, X_test, y_train, y_test = train_test_split(pairwise_data, y, test_size=0.2, random_state=42)

# Convert the pairwise data to numpy arrays
X_train1, X_train2, y_train_labels = zip(*X_train)
X_train1 = np.array(X_train1)
X_train2 = np.array(X_train2)
y_train_labels = np.array(y_train_labels, dtype=int)

# Create and train the RankNet model
model = MLPRegressor(hidden_layer_sizes=(64, 32), activation='relu', random_state=42)
model.fit([X_train1, X_train2], y_train_labels)

# Predict the pairwise rankings on the test set
X_test1, X_test2, y_test_labels = zip(*X_test)
X_test1 = np.array(X_test1)
X_test2 = np.array(X_test2)
y_test_labels = np.array(y_test_labels, dtype=int)
y_test_pred_labels = model.predict([X_test1, X_test2])
y_test_pred_labels = np.round(y_test_pred_labels).astype(int)

# Calculate the mean squared error
mse = mean_squared_error(y_test_labels, y_test_pred_labels)
print(f"Mean Squared Error: {mse}")

# Get the learned weights from the RankNet model
learned_weights = model.coefs_[0].flatten()

# Sort the feature columns based on the learned weights
sorted_columns = [x for _, x in sorted(zip(learned_weights, feature_columns), reverse=True)]
sorted_df = df[['Gene_name'] + sorted_columns + [target_column]]

# Display the sorted dataframe
print(sorted_df)
