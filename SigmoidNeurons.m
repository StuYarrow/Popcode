classdef SigmoidNeurons < Neurons
	
	properties
        width = [];
        maxRate = [];
        backgroundRate = [];
    end
	
	methods
		
		function obj = SigmoidNeurons(varargin)
		% SIGMOIDNEURONS/SIGMOIDNEURONS Constructor for SigmoidNeurons object - population of neurons with sigmoidal tuning curves
		% obj = SigmoidNeurons(preferredStimuli, width, maxRate, backgroundRate, integrationTime, variabilityScheme, variabilityOpts)
		%
		% preferredStimuli - transition stimulus value; point of max gradient
		% width - steepness parameter
		% maxRate - maximum firing rate (Hz)
		% backgroundRate - background (spontaneous) firing rate (Hz)
		% integrationTime - spike counting time per trial
		% variabilityScheme - type of variability model
		% variabilityOpts - vector of options
		%
		% preferredStimulus, width, maxFiringRate and backgroundFiringRate can be scalars or vectors of length popSize.
		% SigmoidNeurons accepts only 1-D stimuli at present.
            
            switch nargin
                case 0
                    % Do nothing
                case 7
                    preferredStimuli = varargin{1};
                    widthIn = varargin{2};
                    maxRateIn = varargin{3};
                    backgroundRateIn = varargin{4};
                    integrationTime = varargin{5};
                    variabilityScheme = varargin{6};
                    variabilityOpts = varargin{7};
                otherwise
                    error('Wrong number of inputs')
            end
            
			% Superclass constructor
            obj = obj@Neurons(1, preferredStimuli, integrationTime, variabilityScheme, variabilityOpts);
            
            if isscalar(widthIn) && isnumeric(widthIn)
				obj.width = double(widthIn(ones(obj.popSize, 1)));
			elseif length(widthIn) == obj.popSize && isvector(widthIn) && isnumeric(widthIn)
				obj.width = reshape(double(widthIn), obj.popSize, 1);
			else
				error('Invalid width value or vector for population size %n', obj.popSize)
            end
            
            if isscalar(maxRateIn) && isnumeric(maxRateIn)
				obj.maxRate = double(maxRateIn(ones(obj.popSize, 1)));
			elseif length(maxRateIn) == obj.popSize && isvector(maxRateIn) && isnumeric(maxRateIn)
				obj.maxRate = reshape(double(maxRateIn), obj.popSize, 1);
			else
				error('Invalid max firing rate value or vector for population size %n', obj.popSize)
            end
            
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
		% r = maxRate ./ (1 + exp(-(stims - centre) ./ width)) + backgroundRate;
            
            if isa(stim, 'StimulusEnsemble')
                assert(stim.dimensionality == 1, 'SigmoidNeurons only supports 1-D stimuli at present')
                stims = repmat(stim.ensemble, [obj.popSize 1]);
                maxRates = repmat(obj.maxRate, [1 stim.n]);
                backgroundRates = repmat(obj.backgroundRate, [1 stim.n]);
                centres = repmat(obj.preferredStimulus, [1 stim.n]);
                widths = repmat(obj.width, [1 stim.n]);
            elseif isa(stim, 'double')
                stims = repmat(stim(:)', [obj.popSize 1]);
                maxRates = repmat(obj.maxRate, [1 length(stim)]);
                backgroundRates = repmat(obj.backgroundRate, [1 length(stim)]);
                centres = repmat(obj.preferredStimulus, [1 length(stim)]);
                widths = repmat(obj.width, [1 length(stim)]);
            else
                error('Invalid stimulus: stim may be a StimulusEnsemble object or vector of stimulus values only')
            end
			
            r = maxRates ./ (1 + exp(-(stims - centres) ./ widths)) + backgroundRates;
		end
		
		function dr = dMeanR(obj, stim)
        % CIRCGAUSSNEURONS/DMEANR calculates the derivative of the tuning curve
		% dr dr = dMeanR(obj, stim)
        
            if isa(stim, 'StimulusEnsemble')
                assert(stim.dimensionality == 1, 'SigmoidNeurons only supports 1-D stimuli at present')
                stims = repmat(stim.ensemble, [obj.popSize 1]);
                maxRates = repmat(obj.maxRate, [1 stim.n]);
                centres = repmat(obj.preferredStimulus, [1 stim.n]);
                widths = repmat(obj.width, [1 stim.n]);
            elseif isa(stim, 'double')
                stims = repmat(stim(:)', [obj.popSize 1]);
                maxRates = repmat(obj.maxRate, [1 length(stim)]);
                centres = repmat(obj.preferredStimulus, [1 length(stim)]);
                widths = repmat(obj.width, [1 length(stim)]);
            else
                error('Invalid stimulus: stim may be a StimulusEnsemble object or vector of stimulus values only')
            end
			
			dr = maxRates ./ (2 .* widths .* cosh((centres - stims) ./ widths) + 2 .* widths);
        end
        
		function obj = remove(obj, nMarg)			
			% Call superclass method
			[obj margMask] = remove@Neurons(obj, nMarg);
            
            if ~isscalar(obj.width)
                obj.width = obj.width(margMask);
            end
            
            if ~isscalar(obj.maxRate)
                obj.maxRate = obj.maxRate(margMask);
            end
            
            if ~isscalar(obj.backgroundRate)
				obj.backgroundRate = obj.backgroundRate(margMask);
            end
        end
        
    end
end