classdef GaussNeurons < Neurons
	
	properties
		width = [];
	end
	
	methods
		
		function obj = GaussNeurons(varargin)
		% GAUSSNEURONS/GAUSSNEURONS Constructor for GaussNeurons object - ppulation of neurons with Gaussian tuning curves
		% n = GaussNeurons(preferredStimulus, width, maxRate, backgroundRate, integrationTime, variabilityScheme, variabilityOpts)
		%
		% popSize - population size (number of neurons)
		% preferredStimulus - preferred stimulus value
		% width - tuning curve width (this would be the variance if the curve was a probability distribution)
		% maxRate - maximum firing rate (Hz)
		% backgroundRate - background (spontaneous) firing rate (Hz)
		% integrationTime - spike counting time per trial
		% variability - a Variability object
		% variabilityOpts - vector of options
		%
		% preferredStimulus, width, maxFiringRate and backgroundFiringRate can be scalars or vectors of length popSize.
		% GaussNeurons accept only 1-D stimuli at present.
			
			if nargin ~= 7
				error('Wrong number of arguments')
			end
			
			obj = obj@Neurons(1, varargin{1}, varargin{3}, varargin{4}, varargin{5}, varargin{6}, varargin{7});

			if length(varargin{3}) == 1 && isnumeric(varargin{3})
				obj.width = double(varargin{3}(ones(obj.popSize, 1)));
			elseif length(varargin{2}) == obj.popSize && isvector(varargin{2}) && isnumeric(varargin{2})
				gn.width = reshape(double(varargin{2}), n.popSize, 1);
			else
				error([inputname(3) ' is not a valid maximum firing rate value or vector for population size ' n.popSize])
			end
		end
		
		function r = meanR(obj, stim)
		% MEANR calculates mean responses to a set of stimuli
		% r = meanR(obj, stimulusEnsemble)
		%
		% r = maxRate * exp(-((stimulus - preferredStimulus)^2 / (2 * width^2))) + backgroundRate

			if ~isa(stim, 'StimulusEnsemble') 
				error([inputname(2) ' is not a valid StimulusEnsemble object'])
			end

			if stim.dimensionality ~= 1
				error('GaussNeurons only supports 1-D stimuli at present')
			end

			stims = repmat(stim.ensemble, obj.popSize, 1);

			maxRate = repmat(obj.maxRate, 1, stim.n);
			backgroundRate = repmat(obj.backgroundRate, 1, stim.n);
			centre = repmat(obj.preferredStimulus, 1, stim.n);
			width = repmat(obj.width, 1, stim.n);

			r = maxRate .* exp(-((double(stims - centre) .^ 2) ./ 2.0 .* width .^ 2)) + backgroundRate;
		end
		
		function dr = dMeanR(obj, stim)
			if ~isa(stim, 'StimulusEnsemble') 
				error([inputname(2) ' is not a valid StimulusEnsemble object'])
			end

			if stim.dimensionality ~= 1
				error('CircGaussNeurons only supports 1-D stimuli at present')
			end

			stims = repmat(stim.ensemble, [obj.Neurons.popSize 1]);
			
			% need to update this stuff for linear Gaussian TC \/ \/ \/
			
			%maxRate = repmat(obj.Neurons.maxRate, 1, stim.n);
			%backgroundRate = repmat(obj.Neurons.backgroundRate, 1, stim.n);
			%centre = repmat(obj.Neurons.preferredStimulus, 1, stim.n);
			%width = repmat(obj.width, 1, stim.n);
			%beta = 1.0 ./ deg2rad(width).^2;

			%dr = -maxRate .* exp(-beta) .* beta .* sind(double(stims - centre)) .* exp(beta .* cosd(double(stims - centre)));
		end
		
	end
end	