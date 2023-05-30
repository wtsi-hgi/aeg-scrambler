import metrics
import config

class Tuner:
    
    def __init__(self, config, chosen_genes, learning_rate):
        self.chosen_genes = chosen_genes
        self.learning_rate = learning_rate        
        self.current_weights = config.get_weights()
        
    def gradient_descent(self):
        """
        Performs iterations to improve weights based on chosen genes,
        finds current residuals of chosen genes, calculates new
        probes from current set of weights, finds the gradients of
        each probe and replaces current weights with the weights of the best
        """
        
        for loop in loops:
            self.current_residuals = self.find_residuals(self, metrics)
            self.calculate_probe_weights(self)
            self.probe_scramble_space(self)
            self.current_weights = self.new_weights
    
    def find_residuals(self, metrics):
        """
        Takes a list of ranked genes and an array of chosen genes, finds the
        mean squared 'error' of their indicies
        """
        
        mse = (metrics.data[metrics.data["Gene_name"] \
            .isin(self.chosen_genes)]["index"] ** 2).mean()
        
        return mse
    
    def calculate_probe_weights(self):
        """
        Applies the learning rate to each weight,
        creating a probe in the positive and negative direction
        of each dimension
        """

        self.probe_weights = []

        for weight in self.weights:
            self.probe_weights.append(
                self.weights[:weight] + 
                (self.weights[weight] + self.learning_rate) +
                self.weights[weight:]
            )
            
            self.probe_weights.append(
                self.weights[:weight] +
                (self.weights[weight] - self.learning_rate) +
                self.weights[weight:]
            )
                
    def probe_scramble_space(self):
        """
        Iterates over list of probes, creates ranked list for each probe,
        finds the residuals and gradients of each and
        selects the steepest gradient
        """
        
        steepest_probe_gradients = [0] * len(self.current_residuals)
        
        for probe in self.probe_weights:
            probe_config = config.Config(probe)
            probe_metrics = metrics.Metrics(probe_config)
            probe_residuals = self.find_residuals(probe_metrics)
            probe_gradients = self.current_residuals - probe_residuals

            if probe_gradients.mean() > steepest_probe_gradients.mean():
                steepest_probe_gradients = probe_gradients
                self.new_weights = probe