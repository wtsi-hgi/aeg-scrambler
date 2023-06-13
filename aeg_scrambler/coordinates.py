from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

convolutions_ext = "_convolutions.wig"


class Coordinates:
    def __init__(self, config, metrics):
        """ """

        self.data = metrics.data
        self.overlaps = metrics.overlaps
        self.generate_step_function_of_overlaps()
        self.convolve_step_function(config)
        self.find_plateaus(config)

    def generate_step_function_of_overlaps(self):
        """
        Defines two new columns on the passed dataframe, to be used in an
        enhancer step function.
        """

        self.data['Enhancer_step_function_x'] = self.data.apply(
            self.step_function_x, axis=1
        )
        self.data['Enhancer_step_function_y'] = self.data.apply(
            self.step_function_y, axis=1
        )

        self.data = self.data.sort_values(
            "Interest_score", ascending=False
        ).reset_index(drop=True)

    def step_function_x(self, row):
        """
        step_function_mask returns an array of genome coordinates within
        the search window
        """

        return np.arange(row['Search_window_start'], row['Search_window_end'])

    def step_function_y(self, row):
        """
        Generates a mask for each overlap, and combines them into a step
        function for the search window
        """

        step_function = np.zeros(
            len(row['Enhancer_step_function_x']),
            dtype=int
        )
        overlapping_elements = self.overlaps.loc[
            self.overlaps["Gene_name"] == row["Gene_name"]
        ]

        for i in range(len(overlapping_elements)):
            start = overlapping_elements.at[
                overlapping_elements.index[i], 'Start'
            ]
            stop = overlapping_elements.at[
                overlapping_elements.index[i], 'End'
            ]

            in_range = np.logical_and(
                row['Enhancer_step_function_x'] >= start,
                row['Enhancer_step_function_x'] <= stop,
            )

            step_function = np.logical_or(step_function, in_range).astype(int)

        return step_function

    def convolve_step_function(self, config):
        """
        X and Y coordinates are generated for convolution
        of element step function with chosen kernel
        """

        self.data[("Enhancer_convolution_x")] = [
            np.empty(0, dtype=float)
        ] * len(self.data)
        self.data[("Enhancer_convolution_y")] = [
            np.empty(0, dtype=float)
        ] * len(self.data)

        self.data = self.data.sort_values(
            "Interest_score", ascending=False
        ).reset_index(drop=True)

        for index, gene in self.data.head(config.convolution_limit).iterrows():
            kernel = self.get_kernel(gene, config)

            convolution_y = np.convolve(
                kernel, gene[("Enhancer_step_function_y")]
            )

            convolution_x = np.arange(
                (gene["Search_window_start"] - (len(kernel) // 2)),
                (
                    gene["Search_window_start"]
                    - (len(kernel) // 2)
                    + len(convolution_y)
                ),
            )

            gene, convolution_x, convolution_y = self.trim_convolution_ends(
                gene, config, convolution_x, convolution_y
            )

            self.data.at[index, ("Enhancer_convolution_x")] = convolution_x
            self.data.at[index, ("Enhancer_convolution_y")] = convolution_y

    def get_kernel(self, gene, config):
        """
        Kernel is generated as numpy array depending on desired shape and size
        """

        size = int(
            (
                config.relative_enhancer_kernel_size
                * (gene["Search_window_end"] - gene["Search_window_start"])
            )
        )
        sigma = int(
            config.relative_enhancer_kernel_sigma
            * (gene["Search_window_end"] - gene["Search_window_start"])
        )

        shape = config.enhancer_kernel_shape

        if shape == "flat":
            kernel = np.ones(size)

        elif shape == "guassian":
            kernel = np.zeros(size)
            np.put(kernel, (size // 2), 10)
            kernel = gaussian_filter1d(kernel, sigma)

        else:
            print("ERROR: Kernel shape is neither Flat nor Guassian")

        return kernel

    def trim_convolution_ends(
        self, gene, config, convolution_x, convolution_y
    ):
        """
        Trims the ends of the convolutions which overlap the ends of the step
        functions.
        """

        print(gene["Gene_name"], "is being trimmed")

        upstream_cut_off = gene["Enhancer_step_function_x"][0]
        downstream_cut_off = gene["Enhancer_step_function_x"][-1]

        upstream_cut_off_index = np.where(convolution_x == upstream_cut_off)
        downstream_cut_off_index = np.where(
            convolution_x == downstream_cut_off
        )

        convolution_x = convolution_x[
            upstream_cut_off_index[0][0] : downstream_cut_off_index[0][0]
        ]
        convolution_y = convolution_y[
            upstream_cut_off_index[0][0] : downstream_cut_off_index[0][0]
        ]

        if convolution_y[0] >= config.plateau_threshold:
            convolution_y[0] = 0
        if convolution_y[-1] >= config.plateau_threshold:
            convolution_y[-1] = 0

        return gene, convolution_x, convolution_y

    def find_plateaus(self, config):
        """
        find_plateaus takes convolved coordinates, and applies a
        threshold to separate the search window into regions based on
        the y-value of each convolved base.
        """

        self.data["Plateau_coordinates"] = ""
        self.data["Plateau_starts"] = ""
        self.data["Plateau_ends"] = ""

        self.data = self.data.sort_values("Interest_score", ascending=False)

        for index, gene in self.data.head(config.convolution_limit).iterrows():
            print(
                f"""Finding plateaus for gene {gene['Gene_name']} 
                ({index + 1} of {config.convolution_limit})..."""
            )

            convolved_x = gene["Enhancer_convolution_x"]
            convolved_y = gene["Enhancer_convolution_y"]
            convolved_y = np.append(convolved_y, 0)
            boolean_below_threshold = convolved_y < config.plateau_threshold
            boolean_below_threshold = np.concatenate(
                (
                    boolean_below_threshold[
                        : (int((gene["Gene_start"] - convolved_x[0])))
                    ],
                    np.full((gene["Gene_end"] - gene["Gene_start"]), False),
                    boolean_below_threshold[
                        ((int(gene["Gene_end"] - convolved_x[0]))) :
                    ],
                )
            )

            boolean_below_threshold = (
                boolean_below_threshold[:-1] != boolean_below_threshold[1:]
            )
            plateau_coordinates = convolved_x[boolean_below_threshold]
            plateau_coordinates = np.concatenate(
                [[convolved_x[0]], plateau_coordinates, [convolved_x[-1]]]
            )

            self.data.at[index, "Plateau_coordinates"] = plateau_coordinates
            self.data.at[index, "Plateau_starts"] = plateau_coordinates[::2]
            self.data.at[index, "Plateau_ends"] = plateau_coordinates[1::2]

            self.data.drop(["Plateau_coordinates"], axis=1)

    def export_convolutions(self, config):
        """
        Coordinates of convolutions are exported to wig file, for each gene
        """

        for index, gene in self.data.head(
            config.convolution_limit
        ).iterrows():
            
            id = config.unique_id[:14]
            convolution_path = f"{config.results_directory}Convolution.{id}.{gene.Gene_name}>.wig"
            
            with open(convolution_path, "w") as convolution_file:
                convolution_file.write(
                    "fixedStep chrom=chr"
                    + gene["Chromosome"]
                    + " start="
                    + str(gene["Enhancer_convolution_x"][0])
                    + " step=1\n"
                )
        
            convolution_signal = pd.DataFrame(
                {
                    "Convolution_signal": self.data.loc[
                        index, "Enhancer_convolution_y"
                    ]
                }
            )

            convolution_signal.to_csv(
                convolution_path,
                sep="\t",
                index=False,
                mode="a",
                header=False,
            )