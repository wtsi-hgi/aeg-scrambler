import os
import numpy as np
import pandas as pd
import pyranges as pr
from scipy import stats
from sklearn.preprocessing import StandardScaler

class Metrics:
    
    def __init__(
        self,
        config,
        gene_annotations,
        regulatory_element_annotations,
        ccle_expression,
        experimental_expression
    ):
        
        self.assign_unique_id()
        self.interesting_features = ["Std",
                                     "Anomalous_score",
                                     "Enhancer_count",
                                     "Enhancer_proportion",
                                     "Specific_gene_expression",
                                     "Gene_size",
                                     "Symmetry_ratio"]
        self.regulatory_data = regulatory_element_annotations.data
        self.merge_genetic_data(
            gene_annotations,
            ccle_expression,
            experimental_expression
        )
        self.find_gene_sizes()
        self.find_interferring_genes(config)
        self.find_search_windows(config)
        self.find_element_overlaps_within_search_window()
        self.count_overlaps_per_gene()
        self.find_nearby_enhancer_densities()
        self.find_symmetry_of_elements()
        self.calculate_interest_score(config)
        
    def merge_genetic_data(
        self, gene_annotations,
        ccle_expression,
        experimental_expression
    ):
        
        self.data = pd.merge(
            gene_annotations.data,
            ccle_expression.data,
            on = "Gene_name",
            how = "inner"
        )
        self.data = pd.merge(
            self.data,
            experimental_expression.data,
            on = "Gene_name",
            how = "inner"
        ) 
        
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
                config.specific_expression_threshold])
        
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
        self.overlaps = gene_search.intersect(elements_search,
                                              strandedness = False)
        self.overlaps = self.overlaps.df
        
        self.data.drop(["Start", "End"], axis = 1)
    
    def count_overlaps_per_gene(self):
        
        """
        Number of specified element overlaps are counted for each gene.
        """

        #overlaps.drop(["Start", "End"], axis = 1)
        self.data = pd.merge(
            self.data, 
            self.overlaps.groupby("Gene_name").size().reset_index(
            name = ("Enhancer_count")), 
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
        
        summed_overlaps = self.overlaps.loc[
            :, 
            ["Gene_name", "Enhancer_proportion"]].groupby(
            ["Gene_name"], 
            as_index = False)["Enhancer_proportion"].sum().reset_index()
        
        self.data = pd.merge(self.data, summed_overlaps, on = "Gene_name")
        
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
        scaled_genes = self.data.loc[:, (
            ["Gene_name"] + self.interesting_features
        )]
        scaled_genes.loc[:, self.interesting_features] = scaler.fit_transform(
            scaled_genes.loc[:, self.interesting_features]
        )
        
        #for feature in self.interesting_features:
        #    
        #    print(feature)
        #    print(type(self.data[feature]))
        #    print(self.data[feature])
        #    print(self.data[feature].to_numpy())
        #    self.data["Z-" + feature] = stats.zscore(self.data[feature]
        #                                             .to_numpy())
        #    print(self.data["Z-" + feature])
        
        scaled_genes = scaled_genes.assign(
            Interest_score = (
                scaled_genes["Std"] * \
                    config.std_weight +
                scaled_genes["Anomalous_score"] * \
                    config.anomalous_expression_weight +
                scaled_genes["Enhancer_count"] * \
                    config.enhancer_count_weight +
                scaled_genes["Enhancer_proportion"] * \
                    config.enhancer_proportion_weight +
                scaled_genes["Specific_gene_expression"] * \
                    config.cell_line_expression_weight +
                scaled_genes["Gene_size"] * \
                    config.gene_size_weight +
                scaled_genes["Symmetry_ratio"] * \
                    config.symmetry_weight
            )
        ).sort_values("Interest_score", ascending = False)
        
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
        
        self.iterate_through_hard_filters(config)
        self.data.sort_values("Interest_score",
                              ascending = False,
                              inplace = True)
        self.data.reset_index(inplace = True)
        
        #for feature in self.interesting_features:
        #    self.data["Z-" + feature] = stats.zscore(self.data[feature])
        
    def assign_unique_id(self):
        
        """
        Generates a unique ID for each dataframe.
        """
        
        self.unique_id = self.__class__.__name__ + str(hash(self))

    def iterate_through_hard_filters(self, config):

        """
        Calls apply_hard_filter for each feature's min and max filter
        """
        
        max_filters = [
            config.std_max, 
            config.anomalous_expression_max, 
            config.enhancer_count_max, 
            config.enhancer_proportion_max, 
            config.cell_line_expression_max, 
            config.gene_size_max,
            config.symmetry_max
        ]
        
        min_filters = [
            config.std_min, 
            config.anomalous_expression_min, 
            config.enhancer_count_min, 
            config.enhancer_proportion_min, 
            config.cell_line_expression_min, 
            config.gene_size_min,
            config.symmetry_min
        ]
        
        for feature in self.interesting_features:
            
            self.apply_hard_filter(
                max_filters[self.interesting_features.index(feature)],
                feature,
                "max"
            )
            
            self.apply_hard_filter(
                min_filters[self.interesting_features.index(feature)],
                feature,
                "min"
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

    def printable_ranks(self):
        
        """
        Returns a dataframe with only the important columns
        """
        
        return self.data.loc[
            :, (["Gene_name"] + ["Interest_score"] + self.interesting_features)
        ].head(50)
        

    def export_gene_scores_report(self, config):
        
        """
        Md5 checksum of config file is generated. Gene prioritisation report
        file is created and checksum is included in name to differentiate
        different configs. Report saved in given location.

        Please note that if you wish to read this file, with pd.read_csv(), 
        then you will need to pass skiprows=56.
        """
        
        id = config.unique_id[:14]
        report_path = config.gene_report_directory + \
            "gene_rankings:<" + id + ">.txt"
        
        open(report_path, "w")
        with open(report_path, "w") as report:

            report.write(config.__str__())
            
        self.data.loc[:, (
            ["Gene_name"] +
            ["Interest_score"] + 
            self.interesting_features +
            ["Scaled_std",
            "Scaled_anomalous_score",
            "Scaled_enhancer_count",
            "Scaled_enhancer_proportion",
            "Scaled_specific_gene_expression",
            "Scaled_gene_size",
            "Scaled_symmetry_ratio",
            #"Z-Std",
            #"Z-Anomalous_score",
            #"Z-Enhancer_count",
            #"Z-Enhancer_proportion",
            #"Z-Specific_gene_expression",
            #"Z-Gene_size",
            #"Z-Symmetry_ratio"
            ])].to_csv(
                report_path,
                sep = "\t",
                index = True,
                mode = "a"
            )