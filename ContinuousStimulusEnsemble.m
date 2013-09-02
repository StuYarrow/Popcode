classdef ContinuousStimulusEnsemble < StimulusEnsemble
	
	properties
		continuous = true;
        width = [];
        lowerLimit = 0;
        upperLimit = 0;
	end
	
	methods

		function obj = ContinuousStimulusEnsemble(varargin)
		% ContinuousStimulusEnsemble/ContinuousStimulusEnsemble Constructor for ContinuousStimulusEnsemble object
		% obj = ContinuousStimulusEnsemble()
        %
            
            assert(nargin == 0, 'Wrong number of arguments')
            
			% Superclass constructor
            obj = obj@StimulusEnsemble();
        end
        
    end
end