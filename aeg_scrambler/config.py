import json

class Config:
    
    def __init__(self):
    
        self.results_directory = "../results/"
        self.gene_prioritisation_report_directory = "../results/"
        self.gene_annotation_reference = ""
        self.regulatory_elements_reference = ""
        self.general_expression_by_cell_line_reference_path = ""
        self.specific_expression_by_cell_line_reference_path = ""
        self.hic_path = ""
        self.reference_genome = ""
        
        self.cell_line_of_interest = "HAP1"
        self.chromosomes_of_interest = \
            ["1", "2", "3", "4", "5", "6", "7", "8",
             "9", "10", "11", "12", "13", "14", "15",
             "16", "17", "18", "19", "20", "21", "22",
             "X", "Y"]
        self.enhancer_epigenetic_flags_of_interest = ["E11"]
        self.quiescent_epigenetic_flags_of_interest = ["E8"]

        self.search_type = "whole_gene"
        self.upstream_search = 500000
        self.downstream_search = 500000

        self.std_hard_filter_max = False
        self.std_hard_filter_min = False
        self.anomalous_expression_hard_filter_max = False
        self.anomalous_expression_hard_filter_min = False
        self.enhancer_count_hard_filter_max = False
        self.enhancer_count_hard_filter_min = False
        self.enhancer_proportion_hard_filter_max = False
        self.enhancer_proportion_hard_filter_min = False
        self.cell_line_expression_hard_filter_max = False
        self.cell_line_expression_hard_filter_min = False
        self.gene_size_hard_filter_max = False
        self.gene_size_hard_filter_min = False
        self.symmetry_hard_filter_max = False
        self.symmetry_hard_filter_min = False

        self.relative_std_weight = 1
        self.relative_anomalous_expression_weight = 1
        self.relative_enhancer_count_weight = 1
        self.relative_enhancer_proportion_weight = 1
        self.relative_cell_line_expression_weight = 1
        self.relative_gene_size_weight = 1
        self.relative_symmetry_weight = 1

        self.enhancer_kernel_shape = "guassian"
        self.enhancer_kernel_size_type = "relative"
        self.absolute_enhancer_kernel_size = 500
        self.absolute_enhancer_kernel_sigma = 3
        self.relative_enhancer_kernel_size = 0.15
        self.relative_enhancer_kernel_sigma = 0.005

        self.quiescent_kernel_shape = "guassian"
        self.quiescent_kernel_size_type = "relative"
        self.absolute_quiescent_kernel_size = 500
        self.absolute_quiescent_kernel_sigma = 3
        self.relative_quiescent_kernel_size = 0.15
        self.relative_quiescent_kernel_sigma = 0.015

        self.cell_line_specific_expression_threshold = 0.01
        self.interferring_gene_overlaps = False

        self.convolution_limit = 2
        self.enhancer_convolution_weight = 1
        self.quiescent_convolution_weight = 1
        self.plateau_threshold = 0.1

        self.inserted_sequence = "ATAACTTCGTATAATGTACATTATACGAAGTTAT"
        self.partial_insertion_per_region = 100
        
        self = self.set_config_from_file(self)
        
    def set_config_from_file(self, path):
        """
        Reads user generated config file and changes variables as necessary
        """    
            
            
        try:
        
            with open(path, "r") as config_file:
                
                self.settings = json.load(config_file)
                
                self.results_directory = \
                    self.settings["results_directory"]
                self.gene_prioritisation_report_directory = \
                    self.settings["gene_prioritisation_report_directory"]
                self.gene_annotation_reference = \
                    self.settings["gene_annotation_reference"]
                self.regulatory_elements_reference = \
                    self.settings["regulatory_elements_reference"]
                self.general_expression_by_cell_line_reference_path = \
                    self.settings["general_expression_by_cell_line_reference_path"]
                self.specific_expression_by_cell_line_reference_path = \
                    self.settings["specific_expression_by_cell_line_reference_path"]
                self.hic_path = \
                    self.settings["hic_path"]
                self.reference_genome = \
                    self.settings["reference_genome"]

                self.cell_line_of_interest = \
                    self.settings["cell_line_of_interest"]
                self.chromosomes_of_interest = \
                    self.settings["chromosomes_of_interest"]
                self.enhancer_epigenetic_flags_of_interest = \
                    self.settings["enhancer_epigenetic_flags_of_interest"]
                self.quiescent_epigenetic_flags_of_interest = \
                    self.settings["quiescent_epigenetic_flags_of_interest"]

                self.search_type = \
                    self.settings["search_type"]
                self.upstream_search = \
                    self.settings["upstream_search"]
                self.downstream_search = \
                    self.settings["downstream_search"]

                self.std_hard_filter_max = \
                    self.settings["std_hard_filter_max"]
                self.std_hard_filter_min = \
                    self.settings["std_hard_filter_min"]
                self.anomalous_expression_hard_filter_max = \
                    self.settings["anomalous_expression_hard_filter_max"]
                self.anomalous_expression_hard_filter_min = \
                    self.settings["anomalous_expression_hard_filter_min"]
                self.enhancer_count_hard_filter_max = \
                    self.settings["enhancer_count_hard_filter_max"]
                self.enhancer_count_hard_filter_min = \
                    self.settings["enhancer_count_hard_filter_min"]
                self.enhancer_proportion_hard_filter_max = \
                    self.settings["enhancer_proportion_hard_filter_max"]
                self.enhancer_proportion_hard_filter_min = \
                    self.settings["enhancer_proportion_hard_filter_min"]
                self.cell_line_expression_hard_filter_max = \
                    self.settings["cell_line_expression_hard_filter_max"]
                self.cell_line_expression_hard_filter_min = \
                    self.settings["cell_line_expression_hard_filter_min"]
                self.gene_size_hard_filter_max = \
                    self.settings["gene_size_hard_filter_max"]
                self.gene_size_hard_filter_min = \
                    self.settings["gene_size_hard_filter_min"]
                self.symmetry_hard_filter_max = \
                    self.settings["symmetry_hard_filter_max"]
                self.symmetry_hard_filter_min = \
                    self.settings["symmetry_hard_filter_min"]

                self.relative_std_weight = \
                    self.settings["relative_std_weight"]
                self.relative_anomalous_expression_weight = \
                    self.settings["relative_anomalous_expression_weight"]
                self.relative_enhancer_count_weight = \
                    self.settings["relative_enhancer_count_weight"]
                self.relative_enhancer_proportion_weight = \
                    self.settings["relative_enhancer_proportion_weight"]
                self.relative_cell_line_expression_weight = \
                    self.settings["relative_cell_line_expression_weight"]
                self.relative_gene_size_weight= \
                    self.settings["relative_gene_size_weight"]
                self.relative_symmetry_weight = \
                    self.settings["relative_symmetry_weight"]

                self.enhancer_kernel_shape = \
                    self.settings["enhancer_kernel_shape"]
                self.enhancer_kernel_size_type = \
                    self.settings["enhancer_kernel_size_type"]
                self.absolute_enhancer_kernel_size = \
                    self.settings["absolute_enhancer_kernel_size"]
                self.absolute_enhancer_kernel_sigma = \
                    self.settings["absolute_enhancer_kernel_sigma"]
                self.relative_enhancer_kernel_size = \
                    self.settings["relative_enhancer_kernel_size"]
                self.relative_enhancer_kernel_sigma = \
                    self.settings["relative_enhancer_kernel_sigma"]

                self.quiescent_kernel_shape = \
                    self.settings["quiescent_kernel_shape"]
                self.quiescent_kernel_size_type = \
                    self.settings["quiescent_kernel_size_type"]
                self.absolute_quiescent_kernel_size = \
                    self.settings["absolute_quiescent_kernel_size"]
                self.relative_quiescent_kernel_size = \
                    self.settings["relative_quiescent_kernel_size"]
                self.relative_quiescent_kernel_sigma = \
                    self.settings["relative_quiescent_kernel_sigma"]

                self.cell_line_specific_expression_threshold = \
                    self.settings["cell_line_specific_expression_threshold"]
                self.interferring_gene_overlaps = \
                    self.settings["interferring_gene_overlaps"]

                self.convolution_limit = \
                    self.settings["convolution_limit"]
                self.enhancer_convolution_weight = \
                    self.settings["enhancer_convolution_weight"]
                self.quiescent_convolution_weight = \
                    self.settings["quiescent_convolution_weight"]
                self.plateau_threshold = \
                    self.settings["plateau_threshold"]
                
                self.inserted_sequence = \
                    self.settings["inserted_sequence"]
                self.partial_insertion_per_region = \
                    self.settings["partial_insertion_per_region"]
                    
        except:
            
            print("ERROR: Config file could not be read.")
            
    def print_config(self):
        
        """
        Prints the current config that the program will use
        """
        
        current_config = vars(self)
        print("\n".join("%s: %s" % setting for
                        setting in
                        current_config.items()))
        
    def generate_config_checksum(self):

        return str(hash(str(self.settings)))