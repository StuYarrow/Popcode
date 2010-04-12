classdef CosNeurons < Neurons
	
	properties
		
	end
	
	methods
		
		function obj = GaussNeurons(varargin)
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
			
			if nargin ~= 7
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
		
		function obj = remove(obj, nMarg)			
			% Call superclass method
			[obj margMask] = remove@Neurons(obj, nMarg);
			
			% Deal with subclass properties
			if length(obj.width) > 1
				obj.width = obj.width(margMask);
			end
		end
		
	end
end













classdef CircGaussNeurons < Neurons
	
	properties
		width = [];
	end
	
	methods

		function obj = CircGaussNeurons(varargin)
		% CIRCGAUSSNEURONS/CIRCGAUSSNEURONS Constructor for CircGaussNeurons object - population of neurons with circular Gaussian tuning curves
		% obj = CircGaussNeurons(popSize, preferredStimulus, width, maxRate, backgroundRate, integrationTime, variabilityScheme, variabilityOpts)
		%
		% popSize - population size (number of neurons)
		% dimensionality - stimulus dimensionality (only 1-D stimuli currently supported)
		% preferredStimulus - preferred stimulus value
		% width - tuning curve width (this would be the variance if the curve was a probability distribution)
		% maxRate - maximum firing rate (Hz)
		% backgroundRate - background (spontaneous) firing rate (Hz)
		% integrationTime - spike counting time per trial
		% variability - a Variability object
		% variabilityOpts - vector of options
		%
		% preferredStimulus, width, maxFiringRate and backgroundFiringRate can be scalars or vectors of length popSize.
		% CircGaussNeurons accept only 1-D stimuli at present.

			% Default constructor
			obj = obj@Neurons(1, varargin{1}, varargin{3}, varargin{4}, varargin{5}, varargin{6}, varargin{7});

			if length(varargin{2}) == 1 && isnumeric(varargin{2})
				obj.width = double(varargin{2}(ones(obj.popSize, 1)));
			elseif length(varargin{2}) == obj.popSize & isvector(varargin{2}) & isnumeric(varargin{2})
				obj.width = reshape(double(), obj.popSize, 1);
			else
				error([inputname(3) ' is not a valid maximum firing rate value or vector for population size ' obj.popSize])
			end
		end
		
		function r = meanR(obj, stim)
		% CIRCGAUSSNEURONS/MEANR calculates mean responses to a set of stimuli
		% r = meanR(obj, stimulusEnsemble)
		%
		% r = maxRate * exp(-((stimulus - preferredStimulus)^2 / (2 * width^2))) + backgroundRate

			if ~isa(stim, 'StimulusEnsemble') 
				error([inputname(2) ' is not a valid StimulusEnsemble object'])
			end

			if stim.dimensionality ~= 1
				error('CircGaussNeurons only supports 1-D stimuli at present')
			end

			stims = repmat(stim.ensemble, obj.popSize, 1);

			maxRate = repmat(obj.maxRate, 1, stim.n);
			backgroundRate = repmat(obj.backgroundRate, 1, stim.n);
			centre = repmat(obj.preferredStimulus, 1, stim.n);
			width = repmat(obj.width, 1, stim.n);
			
			r = maxRate .* exp(-(1.0 ./ degToRad(width).^2) .* (1 - cosd(double(stims - centre)))) + backgroundRate;
		end		
		
		function dr = dMeanR(obj, stim)
		%
		%
		%

			if ~isa(stim, 'StimulusEnsemble') 
				error([inputname(2) ' is not a valid StimulusEnsemble object'])
			end

			if stim.dimensionality ~= 1
				error('CircGaussNeurons only supports 1-D stimuli at present')
			end

			stims = repmat(stim.ensemble, [obj.popSize 1]);

			maxRate = repmat(obj.maxRate, 1, stim.n);
			centre = repmat(obj.preferredStimulus, 1, stim.n);
			width = repmat(obj.width, 1, stim.n);
			beta = 1.0 ./ degToRad(width).^2;

			dr = -maxRate .* exp(-beta) .* beta .* sind(double(stims - centre)) .* exp(beta .* cosd(double(stims - centre)));
			% Would be better to modify this to return the derivative in
			% spikes/s/deg and to remove the rad->deg conversion in the
			% Neurons.fisher() function.
		end
		
		function obj = widthadapt(obj, width, amnt, centre)
		%
		%
		%

			obj = set(obj, 'width', obj.width .* (1 - amnt .* exp(-(1.0 ./ degToRad(width).^2) .* (1 - cosd(double(obj.preferredStimulus - centre))))));
		end
		
		function obj = remove(obj, nMarg)			
			% Call superclass method
			[obj margMask] = remove@Neurons(obj, nMarg);
			
			% Deal with subclass properties
			if length(obj.width) > 1
				obj.width = obj.width(margMask);
			end
		end
		
	end
end