classdef CosNeurons < Neurons
	
	properties
		
	end
	
	methods
		
		function obj = CosNeurons(varargin)
		% COSNEURONS/COSNEURONS Constructor for CosNeurons object - pppulation of neurons with raised cosine tuning curves
		% obj = CosNeurons(preferredStimulus, backgroundRate, integrationTime, variabilityScheme, variabilityOpts)
		%
		% preferredStimulus - preferred stimulus value
		% backgroundRate - background (spontaneous) firing rate (Hz)
		% integrationTime - spike counting time per trial
		% variability - a Variability object
		% variabilityOpts - vector of options
		%
		% preferredStimulus and backgroundFiringRate can be scalars or vectors of length popSize.
		% CosNeurons accept only 1-D stimuli at present.
		
			if nargin ~= 5
				error('Wrong number of arguments')
			end
			
			obj = obj@Neurons(1, varargin{1}, 0, varargin{2}, varargin{3}, varargin{4}, varargin{5});
		end
		
		function r = meanR(obj, stim)
		% MEANR calculates mean responses to a set of stimuli
		% r = meanR(obj, stimulusEnsemble)
		%
		% r = 1/0.86 .* max(0, cos(degToRad(stimulus - preferredStimulus)) - 0.14) + backgroundRate;

			if ~isa(stim, 'StimulusEnsemble') 
				error([inputname(2) ' is not a valid StimulusEnsemble object'])
			end

			if stim.dimensionality ~= 1
				error('CosNeurons only supports 1-D stimuli at present')
			end

			stims = repmat(stim.ensemble, obj.popSize, 1);
			backgroundRate = repmat(obj.backgroundRate, 1, stim.n);
			centre = repmat(obj.preferredStimulus, 1, stim.n);
			
			r = 1/0.86 .* max(0, cos(degToRad(stims - centre)) - 0.14) + backgroundRate;
		end
		
		function dr = dMeanR(obj, stim)
			if ~isa(stim, 'StimulusEnsemble') 
				error([inputname(2) ' is not a valid StimulusEnsemble object'])
			end

			if stim.dimensionality ~= 1
				error('CosNeurons only supports 1-D stimuli at present')
			end

			stims = repmat(stim.ensemble, obj.popSize, 1);
			centre = repmat(obj.preferredStimulus, 1, stim.n);
			
			dr = pi/180 .* (1.16279 .* sin(degToRad(centre - stims)) .* double((cos(degToRad(stims - centre)) - 0.14 > 0.0)));
		end
	end
end