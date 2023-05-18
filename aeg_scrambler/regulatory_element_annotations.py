import pandas as pd

class RegulatoryElementAnnotations:
    
    def __init__(self, config):
    
        self.read_regulatory_element_annotations(
            config.regulatory_elements_reference)
        self.clean_regulatory_elements(
            self.data, config.enhancer_epigenetic_flags_of_interest)
    
    def read_regulatory_element_annotations(self, file_path):
        
        """
        Assigns the regulatory elements bed file to a pandas dataframe
        """

        try:
            
            self.data = \
                pd.read_csv(
                    file_path, sep = "\t",
                    names = ["Chromosome", "Start", "End", "Flag"]
                    )
        
        except:
            
            print("ERROR: Could not read regulatory elements file.")
            
    def clean_regulatory_elements(self, flags_of_interest):
    
        """
        Put regulatory data into correct format and remove unecessary data
        """
            
        self.data["Chromosome"] = \
            self.data["Chromosome"].apply(lambda x : x[3:])
        self.data = self.data[self.data["Flag"].isin(flags_of_interest)]
        self.data = self.data.drop(["Flag"], axis = 1)