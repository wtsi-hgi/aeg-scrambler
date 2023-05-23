import pandas as pd

class RegulatoryAnnotations:
    
    def __init__(self, config):
    
        self.read_regulatory_annotations(config)
        self.clean_regulatory_elements(config)
    
    def read_regulatory_annotations(self, config):
        
        """
        Assigns the regulatory elements bed file to a pandas dataframe
        """

        try:
            
            self.data = \
                pd.read_csv(
                    config.regulatory_elements_reference, sep = "\t",
                    names = ["Chromosome", "Start", "End", "Flag"]
                    )
        
        except:
            
            print("ERROR: Could not read regulatory elements file.")
            
    def clean_regulatory_elements(self, config):
    
        """
        Put regulatory data into correct format and remove unecessary data
        """
            
        self.data["Chromosome"] = \
            self.data["Chromosome"].apply(lambda x : x[3:])
        self.data = self.data[self.data["Flag"]
                              .isin(config
                                    .enhancer_epigenetic_flags_of_interest)]
        self.data = self.data.drop(["Flag"], axis = 1)
        
    def pickle_regulatory_annotations(self, config):
        
        """
        Serialises dataframe and saves as file
        """
        
        self.data.to_pickle(
            config.working_directory +
            "regulatory_annotations"
            )
    
    def unpickle_regulatory_annotations(self, config):
        
        """
        Unserialises dataframe and loads from file
        """
        
        pd.read_pickle(
            config.working_directory +
            "regulatory_annotations"
        )