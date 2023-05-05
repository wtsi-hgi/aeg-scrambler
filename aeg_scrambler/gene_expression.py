import pandas as pd

class GeneExpression:
    
    def __init__(self, config):
        
        self.read_general_expression_data(self, config)
        self.clean_general_expression_data(self, config)
        self.read_specific_expression_data(self, config)
        self.clean_specific_expression_data(self)
    
    def read_general_expression_data(self, config):
        
        # Assigns the expression-for-many-cell-types
        # csv file to a pandas dataframe.

        try:
            
            self.data = \
                pd.read_csv(
                    config.general_expression_by_cell_line_reference_path) \
                        .transpose()  
        
        except:
            
            print("ERROR: General expression data could not be read.")
            
    def clean_general_expression_data(self, config):
    
        # Cleans the general expression data by dropping irrelevant non-numeric
        # columns, converting the first row to headers, setting gene name as
        # index, renaming columns, finding the means and standard deviation of
        # expression for each gene and the z-score of the cell line of
        # interest's expression for each gene, dropping irrelevant columns,
        # and dropping duplicate genes

        self.general_data = self.general_data \
            .drop(["depmapID", "primary_disease"], axis = 0)
        self.general_data.columns = self.data.iloc[0]
        self.general_data = self.general_data.iloc[1:]
        self.general_data = self.general_data.reset_index() \
            .rename(columns = {"index" : "Gene_name"})
        self.general_data = self.data. \
            rename(columns = 
                   {config.cell_line_of_interest : "General_gene_expression"})
        self.general_data = self.find_mean(self.general_data)
        self.general_data = self.find_std(self.general_data)
        self.general_data = self. \
            find_anomalous_score_of_gene_expression(self.general_data)
        self.general_data = self.general_data[["Gene_name",
                                               "General_gene_expression",
                                               "Mean",
                                               "Std",
                                               "Anomalous_score"]]
        self.general_data = self.data \
            .drop_duplicates(keep = False, subset = ["Gene_name"])
            
    def find_mean(data):
        
        # Adds the mean of gene expression to the given expression data frame,
        # giving a mean for each gene
        
        data["Mean"] = data.loc[:, data.columns != "Gene_name"].mean(axis = 1)
        
        return data
        
    def find_std(data):
        
        #Adds the standard deviation to the given expression dataframe, giving the
        #standard deviation across expression of each gene in all provided cell
        #types
        
        data["Std"] = data.loc[:, data.columns != "Gene_name"].std(axis = 1)
        
        return data
    
    def find_anomalous_score_of_gene_expression(data):
    
        # Adds the z-score of each gene based on its expression
        # in the cell line of interest compared to all others
        
        data["Anomalous_score"] = data.apply(lambda gene :
            (gene["General_gene_expression"] -
             gene["Mean"]) / gene["Std"], axis = 1)

        return data
    
    def read_specific_expression_data(self, config):
    
        # Assigns the expression-for-cell-line-of-interest csv file to a pandas
        # dataframe.

        try:
            self.specific_data = \
                pd.read_csv(config. \
                    specific_expression_by_cell_line_reference_path,
                    sep = "\t", 
                    names = ["Gene_name", "Specific_gene_expression"],
                    skiprows = 1)   
        
        except:
            print("ERROR: Specific expression data could not be read.")
        
    def clean_specific_expression_data(self):
        
        #Cleans the expression data specific to the cell line of interest by
        # turning minus infinite strings into a floating point representation
        # of negative infinity, duplicate genes are dropped
        
        self.specific_data["Specific_gene_expression"] = \
            self.specific_data["Specific_gene_expression"]. \
                apply(lambda expression : \
                    0 if expression == "-Inf" else pow(2, expression))
        self.specific_data = \
            self.specific_data. \
                drop_duplicates(keep = False, subset = ["Gene_name"])