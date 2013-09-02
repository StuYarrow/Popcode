classdef CosNeurons < Neurons
	
	properties
        backgroundRate = [];
    end
	
	methods
		
		function obj = CosNeurons(varargin)
		% COSNEURONS/COSNEURONS Constructor for CosNeurons object - population of neurons with raised cosine tuning curves
		% obj = CosNeurons(preferredStimuli, maxRate, backgroundRate, integrationTime, variabilityScheme, variabilityOpts)
		%
		% preferredStimuli - preferred stimulus value
		% width - tuning curve width (this would be the variance if the curve was a probability distribution)
		% maxRate - maximum firing rate (Hz)
		% backgroundRate - background (spontaneous) firing rate (Hz)
		% integrationTime - spike counting time per trial
		% variabilityScheme - type of variability model
		% variabilityOpts - vector of options
		%
		% preferredStimulus, maxFiringRate and backgroundFiringRate can be scalars or vectors of length popSize.
		% CosNeurons accepts only 1-D stimuli at present.
            
            switch nargin
                case 5
                    preferredStimuli = varargin{1};
                    backgroundRateIn = varargin{2};
                    integrationTime = varargin{3};
                    variabilityScheme = varargin{4};
                    variabilityOpts = varargin{5};
                otherwise
                    error('Wrong number of inputs')
            end
            
			% Superclass constructor
            obj = obj@Neurons(1, preferredStimuli, integrationTime, variabilityScheme, variabilityOpts);
            
            if isscalar(backgroundRateIn) && isnumeric(backgroundRateIn)
				obj.backgroundRate = double(backgroundRateIn(ones(obj.popSize, 1)));
			elseif length(backgroundRateIn) == obj.popSize && isvector(backgroundRateIn) && isnumeric(backgroundRateIn)
				obj.backgroundRate = reshape(double(backgroundRateIn), obj.popSize, 1);
			else
				error('Invalid background firing rate value or vector for population size %n', obj.popSize)
            end
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
        % CIRCGAUSSNEURONS/DMEANR calculates the derivative of the tuning curve
		% dr dr = dMeanR(obj, stim)

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
            
            if length(obj.backgroundRate) > 1
				obj.backgroundRate = obj.backgroundRate(margMask);
            end
        end
        
    end
end