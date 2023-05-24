import pandas as pd
import pyranges as pr
import re
import subprocess

class Sequences:
    
    def __init__(self, config, coordinates):
        """"""
        self.data = coordinates.data
        self.initialise_plateau_insertions()
        self.iterate_gene_plateaus(config)
        
    def initialise_plateau_insertions(self):
        """
        Creates empty dataframe ready to accept further inputs.
        """
        
        self.plateau_insertions = pd.DataFrame(
            columns = [
                "Sequence_name", 
                "Insertion_sequence", 
                "Insertion_location", 
                "Plateau_sequence"
            ])
        
    def iterate_gene_plateaus(self, config):
        """
        Saves plateaus associated with each gene as a bed file.
        """
        
        for index, gene in self.data.head(
            config.convolution_limit
        ).iterrows():

            self.plateaus = pd.DataFrame(
                {"Start" : self.data.loc[index, "Plateau_starts"], 
                "End" : self.data.loc[index, "Plateau_ends"]
                }
            )
            self.plateaus["Gene_name"] = gene["Gene_name"]
            self.plateaus["Chromosome"] = "chr" + gene["Chromosome"]
            self.plateaus["Strand"] = gene["Strand"]
            
            self.find_fasta(config)
            self.generate_pridict_input(config)
        
    def find_fasta(self, config):
        """
        Takes coordinates of plateaus and returns FASTA
        sequence from reference genome
        """
        
        plateaus_pr = pr.PyRanges(self.plateaus)
        
        seq = pr.get_sequence(plateaus_pr, config.reference_genome)
        plateaus_pr.seq = seq
        self.plateaus = plateaus_pr.df
        self.plateaus.rename(columns = {"seq" : "Sequence"}, inplace = True)

    def generate_pridict_input(self, config):
        """
        For each plateau, finds partial insertion sites, 
        exports sequences with insertions to csv
        """
        
        self.plateaus.apply(
            self.generate_insertion_prefixes_and_suffixes,
            args = (config),
            axis = 1,
        )
        
        self.plateau_insertions["Sequence"] = self.plateau_insertions.apply(
            self.insert_insertion_sequence,
            axis = 1
        )
        
        self.plateau_insertions.to_csv(
            (config.results_directory + "sequences_for_pridict.csv"),
            index = False, columns = ["Sequence_name", "Sequence"],
            mode = "w", header = False
        )

    def generate_insertion_prefixes_and_suffixes(self, plateau, config):
        """
        Iteratively finds all possible partial insertions in order of length.
        """
        
        plateau_insertion_sites = pd.DataFrame(
            columns = [
                "Sequence_name", 
                "Insertion_sequence",
                "Insertion_location", 
                "Plateau_sequence"
            ])
        
        for num_bases_absent in range(0, len(config.inserted_sequence)):
            for checked_position in range(
                0, (len(config.inserted_sequence) - num_bases_absent)):

                absent_sequence = config.inserted_sequence[
                    checked_position : (checked_position + num_bases_absent)
                ]
                present_sequence = (
                    config.inserted_sequence[:checked_position] 
                    + config.inserted_sequence[
                        (checked_position + num_bases_absent):
                    ]
                )
                insertion_positions = self.find_prefix_suffix_in_plateau(
                    plateau,
                    present_sequence
                )
                
                for position in insertion_positions:
                    new_row = pd.Series({
                        "Sequence_name" : (
                            plateau["Gene_name"] + 
                            " " + plateau["Chromosome"] + 
                            " " + plateau["Strand"] + 
                            " "  + str(plateau["Start"]) + 
                            "-" + str(plateau["End"])
                        ),
                        "Insertion_sequence" : absent_sequence,
                        "Insertion_location" : (position + checked_position),
                        "Plateau_sequence" : plateau["Sequence"]
                    })
                    new_df = pd.DataFrame([new_row])
                    plateau_insertion_sites = pd.concat(
                        [plateau_insertion_sites, new_df],
                        axis = 0, ignore_index = True
                    )
                    
                    if len(plateau_insertion_sites.index) > config.partial_insertions_per_region:
                        plateau_insertions = pd.concat(
                            [plateau_insertions, plateau_insertion_sites], 
                            axis = 0, ignore_index = True
                        )

                        break
        
    def find_prefix_suffix_in_plateau(plateau, present_sequence):
        """
        Searches plateau sequences for partial insertion sequences
        """
        
        insertions = re.finditer(
            pattern = present_sequence, 
            string = plateau["Sequence"]
        )
        insertion_positions = [index.start() for index in insertions]
        
        return insertion_positions
    
    def insert_insertion_sequence(row):
        """
        Adds missing insertion sequence into partial insertions
        found in plateaus
        """
        
        return (
            row["Plateau_sequence"][:row["Insertion_location"]] + 
            "(+" + row["Insertion_sequence"] + ")" + 
            row["Plateau_sequence"][row["Insertion_location"]:]
        )
        
    def run_pridict(config, pridict_input_path):
        subprocess.run([
            "python",
            "pridict_pegRNA_design.py",
            "batch",
            "--input-fname",
            pridict_input_path,
            "--output-fname batchseqs",
            config.pridict_output_path
        ])
        
    def read_pridict_output(self, configuration):
        self.output = pd.read_csv(configuration.pridict_output_path)
        
    def clean_pridict_output(self):
        self.output = self.output.drop([["Original_sequences",
                                         "Edited_sequences"]],
                                       axis = 1)
        self.output = self.output[self.output["Editing_Position"] >= 10]
        self.output = self.output.sort_values(
            "PRIDICT_editing_Score_deep", 
            ascending = False
        ).reset_index()