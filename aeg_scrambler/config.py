import json

class Config:
    
    def __init__(self):
    
        self.results_directory = "../results/"
        self.gene_prioritisation_report_directory = "../results/"
        self.gene_annotation_reference = ""
        self.regulatory_elements_reference = ""
        self.general_expression_by_cell_line_reference_path = ""
        self.specific_expression_by_cell_line_reference_path = ""
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
        
    def set_config_from_file():
            
        try:
        
            with open("../reference/config.json", "r") as config_file:
                
                settings = json.load(config_file)
                
                results_directory = \
                    settings["results_directory"]
                gene_prioritisation_report_directory = \
                    settings["gene_prioritisation_report_directory"]
                gene_annotation_reference = \
                    settings["gene_annotation_reference"]
                regulatory_elements_reference = \
                    settings["regulatory_elements_reference"]
                general_expression_by_cell_line_reference_path = \
                    settings["general_expression_by_cell_line_reference_path"]
                specific_expression_by_cell_line_reference_path = \
                    settings["specific_expression_by_cell_line_reference_path"]
                reference_genome = \
                    settings["reference_genome"]

                cell_line_of_interest = \
                    settings["cell_line_of_interest"]
                chromosomes_of_interest = \
                    settings["chromosomes_of_interest"]
                enhancer_epigenetic_flags_of_interest = \
                    settings["enhancer_epigenetic_flags_of_interest"]
                quiescent_epigenetic_flags_of_interest = \
                    settings["quiescent_epigenetic_flags_of_interest"]

                search_type = \
                    settings["search_type"]
                upstream_search = \
                    settings["upstream_search"]
                downstream_search = \
                    settings["downstream_search"]

                std_hard_filter_max = \
                    settings["std_hard_filter_max"]
                std_hard_filter_min = \
                    settings["std_hard_filter_min"]
                anomalous_expression_hard_filter_max = \
                    settings["anomalous_expression_hard_filter_max"]
                anomalous_expression_hard_filter_min = \
                    settings["anomalous_expression_hard_filter_min"]
                enhancer_count_hard_filter_max = \
                    settings["enhancer_count_hard_filter_max"]
                enhancer_count_hard_filter_min = \
                    settings["enhancer_count_hard_filter_min"]
                enhancer_proportion_hard_filter_max = \
                    settings["enhancer_proportion_hard_filter_max"]
                enhancer_proportion_hard_filter_min = \
                    settings["enhancer_proportion_hard_filter_min"]
                cell_line_expression_hard_filter_max = \
                    settings["cell_line_expression_hard_filter_max"]
                cell_line_expression_hard_filter_min = \
                    settings["cell_line_expression_hard_filter_min"]
                gene_size_hard_filter_max = \
                    settings["gene_size_hard_filter_max"]
                gene_size_hard_filter_min = \
                    settings["gene_size_hard_filter_min"]
                symmetry_hard_filter_max = \
                    settings["symmetry_hard_filter_max"]
                symmetry_hard_filter_min = \
                    settings["symmetry_hard_filter_min"]

                relative_std_weight = \
                    settings["relative_std_weight"]
                relative_anomalous_expression_weight = \
                    settings["relative_anomalous_expression_weight"]
                relative_enhancer_count_weight = \
                    settings["relative_enhancer_count_weight"]
                relative_enhancer_proportion_weight = \
                    settings["relative_enhancer_proportion_weight"]
                relative_cell_line_expression_weight = \
                    settings["relative_cell_line_expression_weight"]
                relative_gene_size_weight= \
                    settings["relative_gene_size_weight"]
                relative_symmetry_weight = \
                    settings["relative_symmetry_weight"]

                enhancer_kernel_shape = \
                    settings["enhancer_kernel_shape"]
                enhancer_kernel_size_type = \
                    settings["enhancer_kernel_size_type"]
                absolute_enhancer_kernel_size = \
                    settings["absolute_enhancer_kernel_size"]
                absolute_enhancer_kernel_sigma = \
                    settings["absolute_enhancer_kernel_sigma"]
                relative_enhancer_kernel_size = \
                    settings["relative_enhancer_kernel_size"]
                relative_enhancer_kernel_sigma = \
                    settings["relative_enhancer_kernel_sigma"]

                quiescent_kernel_shap = \
                    settings["quiescent_kernel_shape"]
                quiescent_kernel_size_type = \
                    settings["quiescent_kernel_size_type"]
                absolute_quiescent_kernel_size = \
                    settings["absolute_quiescent_kernel_size"]
                relative_quiescent_kernel_size = \
                    settings["relative_quiescent_kernel_size"]
                relative_quiescent_kernel_sigma = \
                    settings["relative_quiescent_kernel_sigma"]

                cell_line_specific_expression_threshold = \
                    settings["cell_line_specific_expression_threshold"]
                interferring_gene_overlaps = \
                    settings["interferring_gene_overlaps"]

                convolution_limit = \
                    settings["convolution_limit"]
                enhancer_convolution_weight = \
                    settings["enhancer_convolution_weight"]
                quiescent_convolution_weight = \
                    settings["quiescent_convolution_weight"]
                plateau_threshold = \
                    settings["plateau_threshold"]
                
                inserted_sequence = \
                    settings["inserted_sequence"]
                partial_insertion_per_region = \
                    settings["partial_insertion_per_region"]
                    
        except:
            
            print("ERROR: Config file could not be read.")
            
    def print_config():
        
        