import subprocess
#import pygenometracks

class Tracks:
    
    def __init__(
        self,
        genes,
        config,
        gene_annotations,
        regulatory_annotations,
        metrics,
        coordinates,
        sequences
        ):
        
        self.genes = genes
        self.region = self.generate_region(self.genes)
        self.genes = metrics.data
        self.hic = config.hic_path
        self.gene_annotations = gene_annotations.export(config)
        self.regulatory_annotations = regulatory_annotations.export(config)
        self.density_convolution = coordinates.export_convolutions(config)
        self.plateaus = coordinates.export_plateaus(config) # used to exist?
        self.preferred_insertion_sites = sequences.export_insertion_sites(config) # does not exist yet
        
    def generate_region(self):
        
        """Takes info from metrics dataframe to find region to look at one gene.
        """
        
        chromosome = self.genes.loc[self.genes.loc[self.genes["Gene_name"] == self.gene]]["Chromosome"]
        start = self.genes.loc[self.genes.loc[self.genes["Gene_name"] == self.gene]]["Search_window_start"]
        end = self.genes.loc[self.genes.loc[self.genes["Gene_name"] == self.gene]]["Search_window_end"]
        
        return str(chromosome) + ":" + str(start) + "-" + str(end)
    
    def generate_track_config(self):
        
        """Writes file which PyGenomeTracks can then read.
        """
        
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
        
        """Runs PyGenomeTracks.
        """
        
        subprocess.run([
                "pyGenomeTracks",
                "--tracks",
                "tracks.ini",
                "--input-fname",
                "--region",
                self.region,
            ])