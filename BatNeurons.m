classdef BatNeurons < Neurons
	
	properties
        gau = [];
        sat = [];
        
        A = [];
        mu = [];
        sg = [];
        
        D = [];
        lambda = [];
        ss = [];
	end
	
	methods

		function obj = BatNeurons(varargin)
		% BATNEURONS/BATNEURONS Constructor for BatNeurons object - population of neurons with fitted Gaussian/sigmoid tuning curves
		% obj = CircGaussNeurons(dataStruct, integrationTime)
		%
        % dataStruct - structure containing fitted tuning curves
		% integrationTime - spike counting time per trial
            
            switch nargin
                case 2
                    fits = varargin{1};
                    integrationTime = varargin{2};
                otherwise
                    error('Wrong number of inputs')
            end
            
            % There is no consistent system of labelling with this mixed
            % TC model, so use dummy data for the stimulus labels.  This variable
            % is not used anyway.
            preferredStimuli = zeros(1,length(fits));
            
			% Superclass constructor
			obj = obj@Neurons(1, preferredStimuli, integrationTime, 'poisson', []);
                        
            obj.gau = strcmp({fits.type}, 'gau');
            obj.sat = strcmp({fits.type}, 'sat');
            
            parms = vertcat(fits.parms);
            consts = vertcat(fits.const);
            
            obj.A = parms(obj.gau,1);
            obj.mu = parms(obj.gau,2);
            obj.sg = abs(parms(obj.gau,3)) + consts(obj.gau);
            
            obj.D = parms(obj.sat,1);
            obj.lambda = parms(obj.sat,2);
            obj.ss = abs(parms(obj.sat,3)) + consts(obj.sat);            
		end
		
		function r = meanR(obj, stim)
		% BATNEURONS/MEANR calculates mean responses to a set of stimuli
		% r = meanR(obj, stimulusEnsemble)
        
            if isa(stim, 'StimulusEnsemble') 
				stims = repmat(stim.ensemble, [obj.popSize 1]);
                dim2 = stim.n;
            elseif isa(stim, 'double')
                stims = repmat(stim(:)', [obj.popSize 1]);
                dim2 = length(stim);
            else
                error('Invalid stimulus: stim may be a StimulusEnsemble object or vector of stimulus values only')
            end
            
            Aar = repmat(obj.A, [1 dim2]);
            muAr = repmat(obj.mu, [1 dim2]);
            sgAr = repmat(obj.sg, [1 dim2]);
            r(obj.gau,:) = Aar .* exp( -(stims(obj.gau,:) - muAr).^2 ./ (2 .* sgAr.^2) );
            
            Dar = repmat(obj.D, [1 dim2]);
            lambdaAr = repmat(obj.lambda, [1 dim2]);
            ssAr = repmat(obj.ss, [1 dim2]);
			r(obj.sat,:) = Dar ./ (1 + exp( -(stims(obj.sat,:) - lambdaAr) ./ ssAr ));
		end		
		
		function dr = dMeanR(obj, stim)
		% BATNEURONS/DMEANR calculates derivative of tuning curve
		% r = dMeanR(obj, stimulusEnsemble)

            if isa(stim, 'StimulusEnsemble') 
				stims = repmat(stim.ensemble, [obj.popSize 1]);
                dim2 = stim.n;
            elseif isa(stim, 'double')
                stims = repmat(stim(:)', [obj.popSize 1]);
                dim2 = length(stim);
            else
                error('Invalid stimulus: stim may be a StimulusEnsemble object or vector of stimulus values only')
            end
            
            Aar = repmat(obj.A, [1 dim2]);
            muAr = repmat(obj.mu, [1 dim2]);
            sgAr = repmat(obj.sg, [1 dim2]);            
            dr(obj.gau,:) = Aar .* (muAr - stims(obj.gau,:)) ./ sgAr.^2 .* exp(-(muAr - stims(obj.gau,:)).^2 ./ (2 .* sgAr.^2));
            
            Dar = repmat(obj.D, [1 dim2]);
            lambdaAr = repmat(obj.lambda, [1 dim2]);
            ssAr = repmat(obj.ss, [1 dim2]);            
            dr(obj.sat,:) = Dar ./ (2 .* ssAr .* cosh((lambdaAr - stims(obj.sat,:)) ./ ssAr) + 2 .* ssAr);
        end
		
		function obj = remove(obj, nMarg)			
			% Call superclass method
			[obj margMask] = remove@Neurons(obj, nMarg);
			
			% Deal with subclass properties
            gauMask = margMask(obj.gau);
            satMask = margMask(obj.sat);
            
            obj.gau = obj.gau(margMask);
            obj.sat = obj.sat(margMask);

            obj.A = obj.A(gauMask);
            obj.mu = obj.mu(gauMask);
            obj.sg = obj.sg(gauMask);

            obj.D = obj.D(satMask);
            obj.lambda = obj.lambda(satMask);
            obj.ss = obj.ss(satMask);
		end
		
	end
end