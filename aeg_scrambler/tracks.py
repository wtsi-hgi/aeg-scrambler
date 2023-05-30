import subprocess

class Tracks:
    
    def __init__(self, config):
        
        self.chromosome
        self.start
        self.end
    
    def generate_region(chromosome, start, end):
        
        return str(chromosome) + ":" + str(start) + "-" + str(end)
    
    def assign_track_paths(self):
    
        self.hic = # Supplied by user
        self.gene_annotations = # Needs to be generated after metrics
        self.regulatory_annotations = # Needs to be generated from overlaps
        self.preferred_insertion_sites = # Exported from coordinates
        self.density_convolution = # Exported from coordinates
        self.sequence = # Exported from sequences
    
    def generate_track_config(self):
        
        config = open("tracks.ini", "w")
        config.write("[HiC data]")
        config.write("file = " + self.hic)
        config.write("depth = 5000")
        config.write("title = HiC data")
        config.write("\n")
        config.write("[spacer]")
        config.write("\n")
        config.write("[Genes]")
        config.write("file = " + self.gene_anotations)
        config.write("height = ")
        
    def generate_track_image(self):
        
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