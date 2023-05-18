import numpy as np
import pandas as pd
import pyranges as pr
from scipy import stats
from sklearn.preprocessing import StandardScaler
import hashlib

class Metrics:
    
    def __init__(self,
                 config,
                 gene_annotations,
                 regulatory_element_annotations,
                 gene_expression):
        
        self.interesting_features = [
            "Std", "Anomalous_score", "Enhancer_count",
            "Enhancer_proportion", "Specific_gene_expression",
            "Gene_size", "Symmetry_ratio"
            ]
        self.regulatory_data = regulatory_element_annotations
        self.merge_genetic_data(self, gene_annotations, gene_expression)
        self.find_gene_sizes(self)
        self.find_interferring_genes(self, config)
        self.find_search_windows(self, config)
        self.find_element_overlaps_within_search_window(self)
        self.count_overlaps_per_gene(self, element_type)
        self.find_nearby_enhancer_densities(self)
        self.find_symmetry_of_elements(self)
        self.calculate_interest_score(self, config)
        self.interesting_features = ["Std",
                                     "Anomalous_score",
                                     "Enhancer_count",
                                     "Enhancer_proportion",
                                     "Specific_gene_expression",
                                     "Gene_size",
                                     "Symmetry_ratio"]
        
    def merge_genetic_data(self, gene_annotations, gene_expression):
        
        self.data = pd.merge(gene_annotations.data,
                         gene_expression.general_data,
                         on = "Gene_name", how = "inner")
        self.data = pd.merge(self.data,
                         gene_expression.specific_data,
                         on = "Gene_name", how = "inner") 
        
    def find_gene_sizes(self):
        
        """
        Adds the size of each gene based on its start and end point
        """
        
        self.data["Gene_size"] = self.data["Gene_end"] - self.data["Gene_start"]
    
    def find_interferring_genes(self, config):
    
        """
        For each gene, finds the nearest gene upstream and downstream
        
        "Start" and "End" are generated for use with PyRanges
        module which requires temporary columns "Start" and "End",
        they are removed at the end of function
        """
        
        self.data["Start"] = self.data["Gene_start"]
        self.data["End"] = self.data["Gene_end"]
        
        interferring_genes_search = pr.PyRanges(
            self.data.loc[self.data["Specific_gene_expression"] > \
                config.cell_line_specific_expression_threshold])
        
        gene_search = pr.PyRanges(self.data)
        genes_nearest_upstream_pr = gene_search.nearest(
            interferring_genes_search, how = "upstream", 
            suffix = "_upstream_interferring_gene", 
            overlap = config.interferring_gene_overlaps)
        
        genes_nearest_downstream_pr = gene_search.nearest(
            interferring_genes_search, how = "downstream", 
            suffix = "_downstream_interferring_gene", 
            overlap = config.interferring_gene_overlaps)

        genes_nearest_upstream = genes_nearest_upstream_pr.df
        genes_nearest_downstream = genes_nearest_downstream_pr.df
        
        self.data = pd.merge(
            self.data, 
            genes_nearest_upstream.loc[
                :, 
                ["Gene_name",
                "Start_upstream_interferring_gene",
                "End_upstream_interferring_gene",
                "Gene_name_upstream_interferring_gene"
                ]
            ],
            on = "Gene_name",
            how = "inner"
        )
        
        self.data = pd.merge(
            self.data, 
            genes_nearest_downstream.loc[
                :, 
                ["Gene_name",
                "Start_downstream_interferring_gene",
                "End_downstream_interferring_gene",
                "Gene_name_downstream_interferring_gene"
                ]
            ],
            on = "Gene_name",
            how = "inner"
        )
        
        if not config.interferring_gene_overlaps:
            self.data = self.data.loc[
                self.data["End_upstream_interferring_gene"] < \
                    self.data["Gene_start"]]
            self.data = self.data.loc[
                self.data["Start_downstream_interferring_gene"] > \
                    self.data["Gene_end"]]
            
        self.data.drop(["Start", "End"], axis = 1) 
        
    def find_search_windows(self, config):

        """
        Defines a search window for each gene, based on the number of bases
        upstream and downstream specified, the starting point of the window,
        and whether the gene itself is included. Search window is
        foreshorterned if an interferring gene is within the window, or if
        the window would include negative bases.
        """
        
        if (config.search_type == "whole_gene"):
            downstream_search_start = "Gene_end"
            
        elif (config.search_type == "start_site"):
            downstream_search_start = 'Gene_start'
            
        else: 
            print("ERROR : Invalid search type.")
            
        self.data["Search_window_start"] = self.data.apply(
            lambda gene : gene["Gene_start"] - config.upstream_search 
            if gene["Strand"] == "+" 
            else gene["Gene_start"] - config.downstream_search, axis = 1)
        
        self.data["Search_window_end"] = self.data.apply(
            lambda gene : gene[downstream_search_start] + config.downstream_search 
            if gene["Strand"] == "+" 
            else gene[downstream_search_start] + config.upstream_search, axis = 1)
        
        self.data["Search_window_start"] = self.data.apply(
            lambda gene : 0 
            if gene["Gene_start"] < 0 
            else gene["Search_window_start"], axis = 1)
        
        self.data["Search_window_start"] = self.data.apply(
            lambda gene : gene["End_upstream_interferring_gene"] 
            if gene["Search_window_start"] < \
                gene["End_upstream_interferring_gene"] 
            else gene["Search_window_start"], axis = 1)
        
        self.data["Search_window_end"] = self.data.apply(
            lambda gene : gene["Start_downstream_interferring_gene"] 
            if gene["Search_window_end"] > \
                gene["Start_downstream_interferring_gene"] 
            else gene["Search_window_end"], axis = 1)
            
        self.data["Search_window_size"] = \
            (self.data["Search_window_end"] - self.data["Search_window_start"])
    
    def find_element_overlaps_within_search_window(self):
        
        """
        PyRanges is used to find specified element type overlaps within search
        window given for each gene. "Start" and "End" are generated for
        PyRanges and removed subsequently. Overlaps are stored in a new
        dataframe
        """
        
        self.data["Start"] = self.data["Search_window_start"]
        self.data["End"] = self.data["Search_window_end"]
        
        gene_search = pr.PyRanges(self.data)
        elements_search = pr.PyRanges(self.regulatory_data)
        self.overlaps = gene_search.intersect(elements_search, strandedness = False)
        self.overlaps = self.overlaps.df
        
        self.data.drop(["Start", "End"], axis = 1)
    
    def count_overlaps_per_gene(self, element_type):
        
        """
        Number of specified element overlaps are counted for each gene.
        """

        #overlaps.drop(["Start", "End"], axis = 1)
        self.data = pd.merge(
            self.data, 
            self.overlaps.groupby("Gene_name").size().reset_index(
            name = (element_type + "_count")), 
                on = "Gene_name", 
                how = "inner"
            )
        
    def find_nearby_enhancer_densities(self):

        """
        Density of specifed element overlaps within the given search window is 
        calculated for each gene.
        """
        
        self.overlaps["Enhancer_proportion"] = (
            self.overlaps.loc[:, "End"] - self.overlaps.loc[:, "Start"]) /\
                self.overlaps.loc[:, "Search_window_size"]
        
        self.overlaps = self.overlaps.loc[
            :, 
            ["Gene_name", "Enhancer_proportion"]].groupby(
            ["Gene_name"], 
            as_index = False)["Enhancer_proportion"].sum().reset_index()
        
        self.data = pd.merge(self.data, self.overlaps, on = "Gene_name")
        
    def find_symmetry_of_elements(self):
        
        """
        For each gene calculates the number of elements on one side compared 
        to the other and produces a score of symmetry
        """
        
        self.overlaps["Overlaps_upstream"] = self.overlaps\
            .apply(lambda overlap : True if overlap["End"] <= 
                overlap["Gene_start"] else False, axis = 1)
            
        self.overlaps_upstream = self.overlaps. \
            loc[:, ["Gene_name", "Overlaps_upstream"]] \
                .groupby(["Gene_name"], as_index = False)["Overlaps_upstream"]\
                .sum().reset_index()
                
        self.data = pd.merge(self.data, self.overlaps_upstream,
                             on = "Gene_name")
        
        self.data["Symmetry_ratio"] = \
            2 * (np.absolute((self.data["Overlaps_upstream"] /
                            self.data["Enhancer_count"]) - 0.5))

        self.data.drop(["Overlaps_upstream"], axis = 1)
        
    def calculate_interest_score(self, config):
        
        """
        Various attributes of each gene are scaled and normallised,
        before being weighted and combined into an interest score
        """
        
        scaler = StandardScaler()
        scaled_genes = self.data.loc[:, (["Gene_name"] +
                                         self.interesting_features)]
        scaled_genes.loc[:, self.interesting_features] = \
            scaler.fit_transform(scaled_genes. \
                loc[:, self.interesting_features])
        
        for feature in self.interesting_features:
            self.data["Z-" + feature] = stats.zscore(self.data[feature])
        
        scaled_genes = scaled_genes.assign(
            Interest_score = (
                scaled_genes["Std"] * \
                    config.relative_std_weight +
                scaled_genes["Anomalous_score"] * \
                    config.relative_anomalous_expression_weight +
                scaled_genes["Enhancer_count"] * \
                    config.relative_enhancer_count_weight +
                scaled_genes["Enhancer_proportion"] * \
                    config.relative_enhancer_proportion_weight +
                scaled_genes["Specific_gene_expression"] * \
                    config.relative_cell_line_expression_weight +
                (config.relative_gene_size_weight * \
                    pow((2), (-scaled_genes["Gene_size"] * \
                        config.relative_gene_size_weight * \
                            config.relative_gene_size_weight))) +
                scaled_genes["Symmetry_ratio"] * \
                    config.relative_symmetry_weight
            )
        ).sort_values("Interest_score", ascenconfigng=False)
        
        scaled_genes = scaled_genes.rename(
            columns = {
                "Std" : "Scaled_std",
                "Anomalous_score" : "Scaled_anomalous_score",
                "Enhancer_count" : "Scaled_enhancer_count",
                "Enhancer_proportion" : "Scaled_enhancer_proportion",
                "Specific_gene_expression" : "Scaled_specific_gene_expression",
                "Gene_size" : "Scaled_gene_size",
                "Symmetry_ratio" : "Scaled_symmetry_ratio"
            }
        )
        
        self.data = pd.merge(
            self.data, 
            scaled_genes.loc[:, [
                "Gene_name", 
                "Interest_score", 
                "Scaled_std",
                "Scaled_anomalous_score", 
                "Scaled_enhancer_count",
                "Scaled_enhancer_proportion",
                "Scaled_specific_gene_expression", 
                "Scaled_gene_size",
                "Scaled_symmetry_ratio"
                ]
            ],
            on = "Gene_name"
        )
        self.data = self.iterate_through_hard_filters()
        self.data = self.data.sort_values(
            "Interest_score", ascending = False).reset_index()
        
        for feature in self.interesting_features:
            self.data["Z-" + feature] = stats.zscore(self.data[feature])

    def iterate_through_hard_filters(self, config):

        """
        Calls apply_hard_filter for each feature's min and max filter
        """
        
        max_filters = [
            config.std_hard_filter_max, 
            config.anomalous_expression_hard_filter_max, 
            config.enhancer_count_hard_filter_max, 
            config.enhancer_proportion_hard_filter_max, 
            config.cell_line_expression_hard_filter_max, 
            config.gene_size_hard_filter_max,
            config.symmetry_hard_filter_max
        ]
        
        min_filters = [
            config.std_hard_filter_min, 
            config.anomalous_expression_hard_filter_min, 
            config.enhancer_count_hard_filter_min, 
            config.enhancer_proportion_hard_filter_min, 
            config.cell_line_expression_hard_filter_min, 
            config.gene_size_hard_filter_min,
            config.symmetry_hard_filter_min
        ]
        
        for feature in self.interesting_features:
            self.data = self.apply_hard_filter(
                self.data,
                max_filters[self.interesting_features.index(feature)], 
                    feature, 
                    "max"
            )
            
            self.data = self.apply_hard_filter(
                self.data,
                min_filters[self.interesting_features.index(feature)],
                feature, "min"
            )

    def apply_hard_filter(self, filter, feature, minmax):

        """
        Drops data above max filter or below min filter
        """

        if minmax == "max":
            if filter is not False: self.data = self.data.drop(
                    self.data[self.data[feature] > filter].index
                )
            
        elif minmax == "min":
            if filter is not False: self.data = self.data.drop(
                    self.data[self.data[feature] < filter].index
                )
        
        else: 
            print("ERROR : Could not identify minmax.")

    def export_gene_scores_report(self, configuration):
        
        """
        Not ready for use
        
        Md5 checksum of config file is generated. Gene prioritisation report
        file is created and checksum is included in name to differentiate
        different configs. Report saved in given location.
        """
        
        checksum = self.generate_config_checksum()
        
        with open("""something here?""", "r") as config:
            
            report_name = \
                "gene_prioritisation_report_" + checksum.hexdigest() + ".txt"
            report = \
                open((configuration.gene_prioritisation_report_directory + report_name), "w")
            report.write(config.read() + "\n")
            report.close()
            report = \
                open((configuration.gene_prioritisation_report_directory + report_name), "a")
            self.data.loc[:, (["Gene_name"] +
                            ["Interest_score"] + 
                            self.interesting_features +
                            ["Scaled_std",
                            "Scaled_anomalous_score",
                            "Scaled_enhancer_count",
                            "Scaled_enhancer_proportion",
                            "Scaled_specific_gene_expression",
                            "Scaled_gene_size",
                            "Scaled_symmetry_ratio",
                            "Z-Std",
                            "Z-Anomalous_score",
                            "Z-Enhancer_count",
                            "Z-Enhancer_proportion",
                            "Z-Specific_gene_expression",
                            "Z-Gene_size",
                            "Z-Symmetry_ratio"])].to_csv(
                (configuration.gene_prioritisation_report_directory + report_name),
                sep = "\t", index = True, mode = "a")            
            report.close()
            
    def generate_config_checksum():

        checksum = hashlib.md5()
        
        with open(sys.argv[1], "rb") as config:
            
            for chunk in iter(lambda: config.read(4096), b""):
                
                checksum.update(chunk)
            
        return checksum
            
        
        