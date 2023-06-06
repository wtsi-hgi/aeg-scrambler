import os
import re
import subprocess

import pandas as pd
import pyranges as pr


class Sequences:
    def __init__(self, config, coordinates):
        """Sets up working directory, imports data from coordinate object,
        calls iterate_plateaus_for_each_gene to work on each plateau.

        args:
            config, coordinates
        """

        self.pridict_image_path = config.pridict_image_path
        self.pridict_path = config.pridict_path

        self.hashkey = str(hash(self))
        self.working_directory = "./working/" + self.hashkey + "/"
        if not os.path.exists("./working/"):
            os.mkdir("./working/")
        if not os.path.exists(self.working_directory):
            os.mkdir(self.working_directory)
        self.coordinates = coordinates.data
        self.inserted_sequence = config.inserted_sequence
        self.partial_insertions_per_region = (
            config.partial_insertions_per_region
        )
        self.iterate_plateaus_for_each_gene(config)
        os.rmdir(self.working_directory)

    def iterate_plateaus_for_each_gene(self, config):
        """Iterates through each gene, and applies methods to the gene-associated plateaus as necessary.

        args:
            config
        """

        for index, gene in self.coordinates.head(
            config.convolution_limit
        ).iterrows():
            plateaus = pd.DataFrame(
                {
                    "Start": self.coordinates.loc[index, "Plateau_starts"],
                    "End": self.coordinates.loc[index, "Plateau_ends"],
                },
                dtype=int,
            )
            plateaus = plateaus[plateaus["Start"] < plateaus["End"]]
            plateaus["Gene_name"] = gene["Gene_name"]
            plateaus["Chromosome"] = "chr" + gene["Chromosome"]
            plateaus["Strand"] = gene["Strand"]
            plateaus["Plateau_name"] = plateaus.apply(
                self.name_plateau, axis=1
            )
            plateaus = self.find_fasta(plateaus, config)
            plateaus.apply(self.pridict_each_plateau, axis=1)

    def name_plateau(self, plateau):
        """Generates a unique name for each plateau to help identify it.

        args:
            plateau

        returns:
            plateau_name
        """

        return (
            str(plateau["Gene_name"])
            + ":"
            + str(int(plateau["Start"]))
            + "-"
            + str(int(plateau["End"]))
        )

    def find_fasta(self, plateaus, config):
        """For a given set of plateaus, finds the FASTA sequence of each plateau.

        args:
            plateaus, config

        returns:
            plateaus
        """

        plateaus_pr = pr.PyRanges(plateaus)
        plateaus["Sequence"] = pr.get_sequence(
            plateaus_pr, config.reference_genome_path
        )

        return plateaus

    def pridict_each_plateau(self, plateau):
        """Collates a list of suggested insertions within each plateau based on
        partial sequences found within the plateau, and then prepares each for PRIDICT,
        runs PRIDICT and then parses the output of each.

        args:
            plateau
        """

        plateau_specific_insertions = pd.DataFrame(
            columns=["sequence_name", "editseq"]
        )

        previous_insertions = []

        for insertion_name, insertion_sequence in self.generate_insertions(
            plateau
        ):
            if insertion_name not in previous_insertions:
                previous_insertions.append(insertion_name)

                insertion = pd.DataFrame(
                    {
                        "sequence_name": [insertion_name],
                        "editseq": [insertion_sequence],
                    }
                )

                plateau_specific_insertions = pd.concat(
                    [plateau_specific_insertions, insertion],
                    axis=0,
                    ignore_index=True,
                )

                if (
                    len(plateau_specific_insertions)
                    > self.partial_insertions_per_region
                ):
                    plateau_specific_insertions.to_csv(
                        (self.working_directory + "sequences_for_pridict.csv"),
                        index=False,
                        header=True,
                    )

                    self.run_pridict()

                    plateau_specific_output = pd.DataFrame()

                    for insertion in self.plateau_specific_insertions[
                        "sequence_name"
                    ]:
                        pridict_output = self.read_pridict_output(insertion)
                        pridict_output = self.clean_pridict_output(
                            pridict_output
                        )
                        plateau_specific_output = pd.concat(
                            [plateau_specific_output, pridict_output.head(1)],
                            axis=0,
                            ignore_index=True,
                        )

                        plateau_specific_output.to_csv(
                            self.results_directory
                            + plateau["Plateau_name"]
                            + "_suggested_insertions.tsv",
                            sep="\t",
                        )

                    break

    def generate_insertions(self, plateau):
        """Collates all the insertions found within a given plateau.

        args:
            plateau

        yields:
            insertion_name, insertion_sequence
        """

        for (
            found_sequence,
            absent_sequence,
        ) in self.generate_partial_sequences():
            insertion_indicies = self.find_partial_sequences(
                plateau, found_sequence
            )

            for insertion_index in insertion_indicies:
                insertion_name = self.name_insertion(
                    plateau, insertion_index, absent_sequence
                )
                insertion_sequence = self.insert_insertion_sequence(
                    plateau, insertion_index, absent_sequence
                )

                yield insertion_name, insertion_sequence

    def generate_partial_sequences(self):
        """Produces a size-sorted list of partial sequences of a given sequence."""

        for absent_seq_size in range(0, len(self.inserted_sequence)):
            for insertion_seq_index in range(
                0, (len(self.inserted_sequence) - absent_seq_size)
            ):
                absent_sequence = self.inserted_sequence[
                    insertion_seq_index : (
                        insertion_seq_index + absent_seq_size
                    )
                ]
                found_sequence = (
                    self.inserted_sequence[:insertion_seq_index]
                    + self.inserted_sequence[
                        (insertion_seq_index + absent_seq_size) :
                    ]
                )

                yield found_sequence, absent_sequence

    def find_partial_sequences(self, plateau, present_sequence):
        """Finds already existing partial sequences of insertion sequence within plateaus.

        args:
            plateau, present_sequence

        returns:
            insertion_indicies
        """

        found_instances = re.finditer(
            pattern=present_sequence, string=plateau["Sequence"]
        )
        insertion_indicies = [index.start() for index in found_instances]

        return insertion_indicies

    def name_insertion(self, plateau, insertion_index, absent_sequence):
        """Produces a unique name for each insertion.

        args:
            plateau, insertion_index, absent_sequence

        returns:
            insertion_name
        """

        return (
            str(plateau["Gene_name"])
            + ":"
            + str(plateau["Start"])
            + "<"
            + str((plateau["Start"] + insertion_index))
            + "+"
            + str(len(absent_sequence))
            + ">"
            + str(plateau["End"])
        )

    def insert_insertion_sequence(
        self, plateau, insertion_index, absent_sequence
    ):
        """Inserts a sequence into another in a PRIDICT-parsable manner.

        args:
            plateau, insertion_index, absent_sequence

        returns:
            new_sequence
        """

        return str(
            plateau["Sequence"][:insertion_index]
            + "(+"
            + absent_sequence
            + ")"
            + plateau["Sequence"][insertion_index:]
        )

    def run_pridict(self):
        """Calls the terminal to run PRIDICT."""

        subprocess.run(
            [
                "singularity",
                "exec",
                "--bind",
                "/lustre",
                self.pridict_image_path,
                self.pridict_path,
                "batch",
                "--input-fname",
                self.working_directory + "sequences_for_pridict.csv",
                "--output-fname",
                "sequences_from_pridict",
                "--output-dir",
                self.working_directory,
            ]
        )

    def read_pridict_output(self, insertion_name):
        """Generates a dataframe from the csv produced by PRIDICT.

        args:
            insertion_name

        returns:
            pridict_output
        """

        pridict_output_path = (
            self.working_directory
            + insertion_name
            + "_pegRNA_Pridict_full.csv"
        )
        pridict_output = pd.read_csv(pridict_output_path)

        os.remove(pridict_output_path)

        return pridict_output

    def clean_pridict_output(self, pridict_output):
        """Makes the dataframe of PRIDICT output human-readable.

        args:
            pridict_output

        returns:
            pridict_output
        """

        pridict_output.drop(
            [["Original_sequences", "Edited_sequences"]], axis=1, inplace=True
        )
        pridict_output = pridict_output[self.output["Editing_Position"] >= 10]
        pridict_output = pridict_output.sort_values(
            "PRIDICT_editing_Score_deep", ascending=False
        ).reset_index()
        return pridict_output
