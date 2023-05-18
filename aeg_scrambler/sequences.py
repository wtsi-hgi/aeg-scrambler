import pandas as pd
import pyranges as pr
import re
import subprocess

class Sequences:
    
    def __init__(self, config, coordinates):
        
        self.data = coordinates.data
        self.initialise_possible_plateau_insertions(self)
        
    def initialise_possible_plateau_insertions(self):
        
        """
        Creates empty dataframe ready to accept further inputs
        """
        
        self.possible_plateau_insertions = \
            pd.DataFrame(columns = ["Sequence_name", "Insertion_sequence", 
                                    "Insertion_location", "Plateau_sequence"])
        
    def iterate_gene_plateaus(self, config):
    
        """
        export_plateaus saves plateaus associated with each gene as a bed file
        """
        
        for index, gene in self.data.head(config.convolution_limit).iterrows():
            
            plateaus = pd.DataFrame(
                {"Start" : self.data.loc[index, "Plateau_starts"], 
                "End" : self.data.loc[index, "Plateau_ends"]
                }
            )
            plateaus["Gene_name"] = gene["Gene_name"]
            plateaus["Chromosome"] = "chr" + gene["Chromosome"]
            plateaus["Strand"] = gene["Strand"]
            
            self.find_fasta(plateaus, config)
            self.generate_pridict_input(plateaus)
        
    def find_fasta(plateaus, config):
        
        """
        Takes coordinates of plateaus and returns FASTA
        sequence from reference genome
        """
        
        plateaus_pr = pr.PyRanges(plateaus)
        
        seq = pr.get_sequence(plateaus_pr, config.reference_genome)
        plateaus_pr.seq = seq
        plateaus = plateaus_pr.df
        plateaus = plateaus.rename(columns = {"seq" : "Sequence"})

    def generate_pridict_input(self, config):
        
        """
        For each plateau, finds partial insertion sites, 
        exports sequences with insertions to csv
        """
        
        self.plateaus.apply(self.generate_insertion_prefixes_and_suffixes, axis = 1)
        
        self.possible_plateau_insertions["Sequence"] = self.possible_plateau_insertions\
            .apply(self.insert_insertion_sequence, axis = 1)
        
        self.possible_plateau_insertions\
            .to_csv((config.results_directory + "sequences_for_pridict.csv"),
                    index = False, columns = ["Sequence_name", "Sequence"],
                    mode = "w", header = False)

    def generate_insertion_prefixes_and_suffixes(self, plateau, config):
        
        """
        Iteratively finds all possible partial insertions in order of length
        """
        
        plateau_specific_suggested_insertion_sites = \
            pd.DataFrame(columns = ["Sequence_name", "Insertion_sequence",
                                    "Insertion_location", "Plateau_sequence"])
        
        for number_of_bases_absent in range(0, len(config.inserted_sequence)):
        
            for insertion_sequence_position_being_checked in \
                range(0,(len(config.inserted_sequence) - number_of_bases_absent)):
                
                absent_sequence = \
                    config.inserted_sequence\
                        [insertion_sequence_position_being_checked:\
                            (insertion_sequence_position_being_checked + 
                            number_of_bases_absent)]
                present_sequence = \
                    config.inserted_sequence\
                        [:insertion_sequence_position_being_checked] + \
                            config.inserted_sequence\
                                [(insertion_sequence_position_being_checked +
                                number_of_bases_absent):]
        
                insertion_positions = \
                    self.find_prefix_suffix_in_plateau(plateau, present_sequence)
                
                for position in insertion_positions:
                
                    new_row = pd.Series({"Sequence_name" : \
                        (plateau["Gene_name"] + " " + plateau["Chromosome"] +
                        " " + plateau["Strand"] + " "  + str(plateau["Start"]) +
                        "-" + str(plateau["End"])),
                        "Insertion_sequence" : absent_sequence,
                        "Insertion_location" : (position + \
                            insertion_sequence_position_being_checked),
                        "Plateau_sequence" : plateau["Sequence"]})
                    new_df = pd.DataFrame([new_row])
                    plateau_specific_suggested_insertion_sites = \
                        pd.concat(\
                            [plateau_specific_suggested_insertion_sites, new_df],
                            axis = 0, ignore_index = True)
                    
                    if len(plateau_specific_suggested_insertion_sites.index) > \
                        config.partial_insertions_per_region:
                        
                        possible_plateau_insertions = \
                            pd.concat([possible_plateau_insertions, \
                                plateau_specific_suggested_insertion_sites], \
                                    axis = 0, ignore_index = True)
                        
                        break
        
    def find_prefix_suffix_in_plateau(plateau, present_sequence):
        
        """
        Searches plateau sequences for partial insertion sequences
        """
        
        insertions = \
            re.finditer(pattern = present_sequence, string = plateau["Sequence"])
        insertion_positions = [index.start() for index in insertions]
        
        return insertion_positions

    def insert_insertion_sequence(row):
        
        """
        Adds missing insertion sequence into partial insertions found in plateaus
        """
        
        return (row["Plateau_sequence"][:row["Insertion_location"]] + 
                    "(+" +
                    row["Insertion_sequence"] +
                    ")" + 
                    row["Plateau_sequence"][row["Insertion_location"]:])
        
    def run_pridict(config, pridict_input_path):
        
        subprocess.run(
            [
                "python",
                "pridict_pegRNA_design.py",
                "batch",
                "--input-fname",
                pridict_input_path,
                "--output-fname batchseqs",
                config.pridict_output_path
            ]
            )
        
    def read_pridict_output(self, configuration):
        
        self.output = pd.read_csv(configuration.pridict_output_path)
        
    def clean_pridict_output(self):
        
        self.output = self.output.drop([["Original_sequences",
                                         "Edited_sequences"]],
                                       axis = 1)
        self.output = self.output[self.output["Editing_Position"] >= 10]
        self.output = self.output.sort_values(
            "PRIDICT_editing_Score_deep", ascending = False).reset_index()