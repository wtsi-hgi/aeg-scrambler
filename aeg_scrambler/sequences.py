import pandas as pd
import pyranges as pr
import re
import subprocess

class Sequences:
    
    def __init__(self, config, coordinates):
        
        self.pridict_image_path = config.pridict_image_path
        self.pridict_path = config.pridict_path
        
        self.coordinates = coordinates.data
        self.results_directory = config.results_directory
        self.inserted_sequence = config.inserted_sequence
        self.partial_insertions_per_region = \
            config.partial_insertions_per_region
        self.iterate_plateaus_for_each_gene(config)

    def iterate_plateaus_for_each_gene(self, config):
        
        for index, gene in self.coordinates.head(
            config.convolution_limit
        ).iterrows():

            plateaus = pd.DataFrame({
                "Start" : self.coordinates.loc[index, "Plateau_starts"], 
                "End" : self.coordinates.loc[index, "Plateau_ends"]
            })
            plateaus["Start"] = plateaus["Start"].astype(int)
            plateaus["End"] = plateaus["End"].astype(int)
            plateaus = plateaus[plateaus["Start"] < plateaus["End"]]
            plateaus["Gene_name"] = gene["Gene_name"]
            plateaus["Chromosome"] = "chr" + gene["Chromosome"]
            plateaus["Strand"] = gene["Strand"]
            plateaus["Plateau_name"] = plateaus.apply(
                self.name_plateau, axis = 1
            )
            plateaus = self.find_fasta(plateaus, config)
            plateaus.apply(self.generate_pridict_input, axis = 1)
                
    def name_plateau(self, plateau):
        
        return (
            str(plateau["Gene_name"]) +
            ":" +
            str(int(plateau["Start"])) +
            "-" +
            str(int(plateau["End"]))
        )
    
    def find_fasta(self, plateaus, config):
        
        plateaus_pr = pr.PyRanges(plateaus)
        plateaus["Sequence"] = pr.get_sequence(plateaus_pr, config.reference_genome_path)
        
        return(plateaus)
    
    def generate_pridict_input(self, plateau):
        
        plateau_specific_insertions = pd.DataFrame(
            columns = [
                "sequence_name",
                "editseq"
            ]
        )
        
        for insertion_name, insertion_sequence in self.generate_insertions(plateau):
            
            insertion = pd.DataFrame({
                "sequence_name" : [insertion_name],
                "editseq" : [insertion_sequence]
            })
            
            plateau_specific_insertions = pd.concat(
                [plateau_specific_insertions, insertion],
                axis = 0,
                ignore_index = True
            )
            
            if len(
                plateau_specific_insertions
            ) > self.partial_insertions_per_region:
                
                    plateau_specific_insertions.to_csv(
                        (self.results_directory +
                         "sequences_for_pridict.csv"),
                        index = False,
                        header = True
                    )

                    self.run_pridict()

                    break

    def generate_insertions(self, plateau):
        
        for found_sequence, absent_sequence in \
            self.generate_partial_sequences():
            
            insertion_indicies = self.find_partial_sequences(
                plateau,
                found_sequence
            )
            
            for insertion_index in insertion_indicies:
                
                insertion_name = self.name_insertion(plateau, insertion_index, absent_sequence)
                insertion_sequence = self.insert_insertion_sequence(
                    plateau,
                    insertion_index,
                    absent_sequence
                )
                
                yield insertion_name, insertion_sequence

    def generate_partial_sequences(self):
        
        for absent_seq_size in range (0, len(self.inserted_sequence)):
            
            for insertion_seq_index in range(
                0,
                (len(self.inserted_sequence) - absent_seq_size)
            ):
                
                absent_sequence = self.inserted_sequence[
                    insertion_seq_index : (
                        insertion_seq_index + absent_seq_size
                )]
                found_sequence = (
                    self.inserted_sequence[:insertion_seq_index] +
                    self.inserted_sequence[(
                        insertion_seq_index +
                        absent_seq_size
                        ):
                ])
                
                yield found_sequence, absent_sequence
    
    def find_partial_sequences(self, plateau, present_sequence):
        
        """Finds already existing partial sequences of insertion sequence within plateaus
        """
        
        found_instances = re.finditer(
            pattern = present_sequence, 
            string = plateau["Sequence"]
        )
        insertion_indicies = [index.start() for index in found_instances]
        
        return insertion_indicies        
             
    def name_insertion(self, plateau, insertion_index, absent_sequence):
        
        
        return (
            str(plateau["Gene_name"]) +
            ":" +
            str(plateau["Start"]) +
            "<" +
            str((plateau["Start"] + insertion_index)) +
            "+" +
            str(len(absent_sequence)) +
            ">" +
            str(plateau["End"])
        )
             
    def insert_insertion_sequence(
        self,
        plateau,
        insertion_index,
        absent_sequence
    ):
        
        return str(
            plateau["Sequence"][:insertion_index] + 
            "(+" +
            absent_sequence +
            ")" + 
            plateau["Sequence"][insertion_index:]
        )

    def run_pridict(self):
        
        subprocess.run([
            
            "singularity",
            "exec",
            "--bind",
            "/lustre",
            self.pridict_image_path,
            self.pridict_path,
            "batch",
            "--input-fname",
            self.results_directory + "sequences_for_pridict.csv",
            "--output-fname",
            "sequences_from_pridict",
            "--output-dir",
            self.results_directory
            
        ])