import subprocess
#import pygenometracks


class Tracks:
    def __init__(
        self,
        config,
        gene_annotations,
        regulatory_annotations,
        metrics,
        coordinates,
        sequences,
    ):
        self.metrics = metrics.data
        
        self.hic = config.hic_path
        self.all_gene_annotations = gene_annotations.export(config)
        self.all_regulatory_annotations = regulatory_annotations.export(config)
        
        self.iterate_genes(config)
        
    def iterate_genes(self, config):
        
        """Iterates through the methods which are specific to each gene,
        generating a gene-specific region, convolution file, plateaus file,
        and generates a file of the best insertion in each plateau
        """
        
        for gene in config.genes_of_interest:
            
            id = config.unique_id[:14]
            region = self.generate_region(gene)
            density_convolution = f"{config.results_directory}convolution.{id.lower()}.{gene.Gene_name}.wig"
            density_convolution = self.wig_to_bigwig(density_convolution)
            insertions = self.get_suggested_insertions(config)
        
    def generate_region(self, gene):
        """Takes info from metrics dataframe to find region to look at one gene."""

        chromosome = self.genes.loc[
            self.metrics.loc[self.metrics["Gene_name"] == gene]
        ]["Chromosome"]
        start = self.genes.loc[
            self.metrics.loc[self.metrics["Gene_name"] == gene]
        ]["Search_window_start"]
        end = self.metrics.loc[
            self.metrics.loc[self.metrics["Gene_name"] == gene]
        ]["Search_window_end"]

        return str(chromosome) + ":" + str(start) + "-" + str(end)
    
    def wig_to_bigwig(self):
        
        pass

    def get_suggested_insertions(self, config):
        
        id = config.unique_id[:14]
        insertions_wildcard = f"{config.results_directory}all_insertions.{id.lower()}.*"
        all_insertions = subprocess.run([
            "ls",
            insertions_wildcard
        ])
        for file in all_insertions.stdout:
            suggested_insertion = subprocess.run([
                "sed",
                "-n",
                "2p",
                file
            ])
            

    def generate_track_config(self):
        """Writes file which PyGenomeTracks can then read."""

        tracks_config = open("tracks.ini", "w")
        tracks_config.write("[HiC data]\n")
        tracks_config.write("file = " + self.hic + "\n")
        tracks_config.write("title = HiC data\n")
        tracks_config.write("\n")
        tracks_config.write("[spacer]\n")
        tracks_config.write("\n")
        tracks_config.write("[Genes]\n")
        tracks_config.write("file = " + self.gene_anotations + "\n")
        tracks_config.write("title = Gene annotations\n")
        tracks_config.write("\n")
        tracks_config.write("[spacer]\n")
        tracks_config.write("\n")
        tracks_config.write("[Enhancers]\n")
        tracks_config.write("file = " + self.regulatory_anotations + "\n")
        tracks_config.write("title = Enhancers\n")

    def generate_track_image(self):
        """Runs PyGenomeTracks."""

        subprocess.run(
            [
                "pyGenomeTracks",
                "--tracks",
                "tracks.ini",
                "--input-fname",
                "--region",
                self.region,
            ]
        )
