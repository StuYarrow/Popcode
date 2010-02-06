function display(obj)
% DISPLAY(obj) Displays summary for Neurons object

if ~isscalar(obj)
	retStr = 'Array of Neurons objects';
	return
end

disp(sprintf(['Population size: %.0f\nStimulus dimensionality: %.0f'], obj.Neurons.popSize, double(obj.Neurons.dimensionality)))
disp('Preferred stimuli:')
disp(obj.Neurons.preferredStimulus)
disp('Width:')
disp(obj.width)
disp('Maximum firing rate (Hz):')
disp(obj.Neurons.maxRate)
disp('Background firing rate (Hz):')
disp(obj.Neurons.backgroundRate)
disp(sprintf(['Integration time: %.3f s\n'...
			 'Variability distribution: %s\n'...
			 'Variability coefficient a: %.3f\n'...
			 'Variability exponent alpha: %.3f'], obj.Neurons.integrationTime, obj.Neurons.distribution, obj.Neurons.a, obj.Neurons.alpha))
disp('Correlation matrix R:')
disp(obj.Neurons.R)