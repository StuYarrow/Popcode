classdef GaussNeurons < Neurons
	
	properties
		width = [];
        maxRate = [];
        backgroundRate = [];
	end
	
	methods
		
        function obj = GaussNeurons(varargin)
		% GAUSSNEURONS/GAUSSNEURONS Constructor for GaussNeurons object - population of neurons with Gaussian tuning curves
		% obj = GaussNeurons(preferredStimuli, width, maxRate, backgroundRate, integrationTime, variabilityScheme, variabilityOpts)
		%
		% preferredStimuli - preferred stimulus value
		% width - tuning curve width (this would be the variance if the curve was a probability distribution)
		% maxRate - maximum firing rate (Hz)
		% backgroundRate - background (spontaneous) firing rate (Hz)
		% integrationTime - spike counting time per trial
		% variabilityScheme - type of variability model
		% variabilityOpts - vector of options
		%
		% preferredStimulus, width, maxFiringRate and backgroundFiringRate can be scalars or vectors of length popSize.
		% GaussNeurons accept only 1-D stimuli at present.
            
            switch nargin
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
		% r = maxRate * exp(-((stimulus - preferredStimulus)^2 / (2 * width^2))) + backgroundRate

			if ~isa(stim, 'StimulusEnsemble') 
				error([inputname(2) ' is not a valid StimulusEnsemble object'])
			end

			if stim.dimensionality ~= 1
				error('GaussNeurons only supports 1-D stimuli at present')
			end

			stims = repmat(stim.ensemble, obj.popSize, 1);

			maxRateArr = repmat(obj.maxRate, 1, stim.n);
			backgroundRateArr = repmat(obj.backgroundRate, 1, stim.n);
			centreArr = repmat(obj.preferredStimulus, 1, stim.n);
			widthArr = repmat(obj.width, 1, stim.n);

			r = maxRateArr .* exp(-((stims - centreArr).^2) ./ (2.0 .* widthArr.^2)) + backgroundRateArr;
		end
		
		function dr = dMeanR(obj, stim)
			if ~isa(stim, 'StimulusEnsemble') 
				error([inputname(2) ' is not a valid StimulusEnsemble object'])
			end

			if stim.dimensionality ~= 1
				error('GaussNeurons only supports 1-D stimuli at present')
			end

			stims = repmat(stim.ensemble, [obj.Neurons.popSize 1]);
						
			maxRateArr = repmat(obj.maxRate, 1, stim.n);
			centreArr = repmat(obj.preferredStimulus, 1, stim.n);
			widthArr = repmat(obj.width, 1, stim.n);
            
            dr = maxRateArr .* (centreArr - stims) ./ widthArr.^2 .* exp(-(centreArr - stims).^2 ./ (2 .* widthArr.^2));
        end
		
        function obj = remove(obj, nMarg)
        % CIRCGAUSSNEURONS/DMEANR calculates the derivative of the tuning curve
		% dr dr = dMeanR(obj, stim)

			% Call superclass method
			[obj margMask] = remove@Neurons(obj, nMarg);
			
			% Deal with subclass properties
            if length(obj.width) > 1
				obj.width = obj.width(margMask);
            end
            
            if length(obj.maxRate) > 1
				obj.maxRate = obj.maxRate(margMask);
            end
            
            if length(obj.backgroundRate) > 1
				obj.backgroundRate = obj.backgroundRate(margMask);
            end
		end
	end
end	