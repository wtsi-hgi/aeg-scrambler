import json
import os


class Config:
    def __init__(self, path=None) -> None:
        """
        Assign default values to configuration variables.
        """

        # Assign ID
        self.unique_id = self.assign_unique_id()

        # Important dataframe columns
        self.columns = {
            'Gene_name': str,
            'Scaled_std': float,
            'Scaled_anomalous_score': float,
            'Scaled_enhancer_count': float,
            'Scaled_enhancer_proportion': float,
            'Scaled_specific_gene_expression': float,
        }

        # File paths
        reference_dir = "./reference_data/"
        results_dir = "./results/"

        self.results_directory = results_dir
        if not os.path.exists(self.results_directory):
            os.mkdir(self.results_directory)
        self.gene_report_directory = results_dir

        self.gene_annotation_path = (
            reference_dir + "Homo_sapiens.GRCh38.108.gtf"
        )
        self.regulatory_elements_path = reference_dir + "HAP1_15_segments.bed"
        self.ccle_expression_path = reference_dir + "CCLE_expressions_22Q4.csv"
        self.experimental_expression_path = (
            reference_dir + "Expression_HAP1_c12.tsv"
        )
        self.reference_genome_path = reference_dir + "genome.fa"
        self.hic_path = reference_dir + "4DNFI1E6NJQJ.hic"

        # Experimental specific settings
        self.cell_line_of_interest = "HAP1"
        self.chromosomes_of_interest = [str(i) for i in range(23)] + ['X', 'Y']
        self.flags_of_interest = ["E11"]

        # Search settings
        self.search_type = "whole_gene"
        self.upstream_search = 500000
        self.downstream_search = 500000

        # Filters
        self.std_max = False
        self.std_min = False
        self.anomalous_expression_max = False
        self.anomalous_expression_min = False
        self.enhancer_count_max = False
        self.enhancer_count_min = False
        self.enhancer_proportion_max = False
        self.enhancer_proportion_min = False
        self.cell_line_expression_max = False
        self.cell_line_expression_min = False
        self.gene_size_max = False
        self.gene_size_min = False
        self.symmetry_max = False
        self.symmetry_min = False

        # Weights
        self.std_weight = 1
        self.anomalous_expression_weight = 1
        self.enhancer_count_weight = 1
        self.enhancer_proportion_weight = 1
        self.cell_line_expression_weight = 1
        self.gene_size_weight = 1
        self.symmetry_weight = 1

        # Enhancer kernel
        self.enhancer_kernel_shape = "guassian"
        self.enhancer_kernel_size_type = "relative"
        self.absolute_enhancer_kernel_size = 500
        self.absolute_enhancer_kernel_sigma = 3
        self.relative_enhancer_kernel_size = 0.15
        self.relative_enhancer_kernel_sigma = 0.005

        # Interferring gene settings
        self.specific_expression_threshold = 0.01
        self.interferring_gene_overlaps = False

        # Convolution settings
        self.convolution_limit = 3
        self.plateau_threshold = 0.1

        # Sequence inserting settings
        self.inserted_sequence = "ATAACTTCGTATAATGTACATTATACGAAGTTAT"
        self.partial_insertions_per_region = 10

        # Temporary Pridict paths
        self.pridict_image_path = ""  # temp
        self.pridict_path = ""  # temp
        self.pridict_output_path = ""  # temp

        # Load config from file
        self.set_config_from_file(path)

        print("Current Config is " + self.unique_id[:14])

    def set_config_from_file(self, path) -> None:
        """
        If user supplies a path, reads and changes variables as necessary.
        """

        try:
            with open(path, "r") as config_file:
                settings = json.load(config_file)
                for key, value in settings.items():
                    if hasattr(self, key):
                        setattr(self, key, value)

        except FileNotFoundError:
            print(
                """WARN: Config file not found.
            Using default config values instead"""
            )
            print()
        except Exception as e:
            print(
                f"""WARN: Failed to read config file: {e}
            Using default config values instead."""
            )
            print()

    def __str__(self) -> str:
        """
        Returns a string representation of the config object.
        """

        current_config = vars(self)
        config_str = ""
        for key, value in current_config.items():
            config_str += f"{key}: {value}\n"

        return config_str

    def assign_unique_id(self):
        """
        Generates a unique ID for each dataframe.
        """

        return self.__class__.__name__ + str(hash(self))

    def get_weights(self) -> list:
        """
        Returns the weights associated with the config object.
        """

        return [
            self.std_weight,
            self.anomalous_expression_weight,
            self.enhancer_count_weight,
            self.enhancer_proportion_weight,
            self.cell_line_expression_weight,
            self.gene_size_weight,
            self.symmetry_weight,
        ]

    def get_filters(self) -> list:
        """
        Returns the filters associated with the config object.
        """

        return [
            self.std_max,
            self.std_min,
            self.anomalous_expression_max,
            self.anomalous_expression_min,
            self.enhancer_count_max,
            self.enhancer_count_min,
            self.enhancer_proportion_max,
            self.enhancer_proportion_min,
            self.cell_line_expression_max,
            self.cell_line_expression_min,
            self.gene_size_max,
            self.gene_size_min,
            self.symmetry_max,
            self.symmetry_min,
        ]
